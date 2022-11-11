/**
 * @file   BeamGateInfoFromTracks_module.cc
 * @brief  Writes a collection of sim::BeamGateInfo from track times.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 22, 2021
 */

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time
#include "lardataalg/Utilities/intervals_fhicl.h" // nanoseconds from FHiCL
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds, ...
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <memory> // std::make_unique()
#include <utility> // std::move()


//------------------------------------------------------------------------------
namespace icarus::trigger { class BeamGateInfoFromTracks; }

/**
 * @brief Writes a set collection of beam gates into each event.
 * 
 * This module writes a list of `sim::BeamGateInfo` based on the time associated
 * to a selection of reconstructed tracks.
 * 
 * It may be used as input to modules which require to operate on beam gates,
 * to select time(s) around the reconstructed (and selected) tracks.
 * 
 * 
 * Input data products
 * ====================
 * 
 * This module acts on a _selection_ of tracks, which implies that the input is
 * a set of _pointers_ to tracks rather than to an actual track collection.
 * For each track, an associated time is required.
 * 
 * * `std::vector<art::Ptr<recob::PFParticle>>` (tag from `T0selProducer`):
 *   list of the selected particles to be considered
 * * `art::Assns<anab::T0, recob::PFParticle>` (tag from `T0Producer`):
 *   association between particles and their time, which must include all
 *   pointers to the tracks listed in `T0selProducer`: each of those particles
 *   must be associated with exactly one `anab::T0`. The time convention is
 *   described in @ref icarus_BeamGateInfoFromTracks_times "its own section".
 * 
 * 
 * Output data products
 * =====================
 *
 * * `std::vector<sim::BeamGateInfo>`: a collection of as many
 *   `sim::BeamGateInfo` as the `BeamGates` configuration parameter entries,
 *   with the content from them.
 * * `art::Assns<sim::BeamGateInfo, recob::PFParticle>`: courtesy association
 *   between the selected particle and its gate, in the same order as the
 *   gates in their data product.
 * * `art::Assns<sim::BeamGateInfo, anab::T0>`: courtesy association between
 *   the particle time and its gate, in the same order as the gates in their
 *   data product.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse online description of the parameters is printed by running
 * `lar --print-description BeamGateInfoFromTracks`.
 * 
 * * `T0selProducer` (input tag, mandatory): the list of pointers to the
 *   particles to be considered.
 * * `T0Producer` (input tag, mandatory): the association of particles with
 *   their times.
 * * `GateStartOffset` (time string, mandatory): the time of the opening of
 *   the gate with respect to the time of the particle; a positive offset means
 *   that the gate starts _after_ the time of the particle.
 * * `GateEndOffset` (time string, mandatory): the time of the closing of
 *   the gate with respect to the time of the particle.
 * * `GateType` (string, default: "unknown"): type of the gates being written;
 *   see online description for the configuration keys representing the values
 *   of `sim::GateType`.
 * * `LogCategory` (string, default: `BeamGateInfoFromTracks`): name of the
 *     output stream category for console messages (managed by MessageFacility
 *     library).
 *
 *
 * @note Time strings are strings with a value and its mandatory unit. For
 *       example, 5.5 microseconds are expressed as `"5.5 us"` or `"5500 ns"`.
 *
 * Time scales
 * ============
 *
 * @anchor icarus_BeamGateInfoFromTracks_times
 * 
 * This module was introduced as a mean to determine the proper time intervals
 * when to apply a trigger simulation.
 * It is evident that for this use case a precise result requires the timing
 * of the gates carefully aligned with the timing of the PMT.
 * (the convention is that the PMT waveforms have a timestamp in the
 * @ref DetectorClocksOpticalElectronicsTime "electronics time reference").
 * 
 * The time of the gates created by this module is offered as
 * `sim::BeamGateInfo`; LArSoft prescribes it to be be specified in nanoseconds
 * and in @ref DetectorClocksSimulationTime "simulation time reference").
 * Designed for simulation, this time scale refers to the opening of the beam
 * gate (or its surrogate definition in samples with no beam to be gated).
 * Or so is my understanding of it.
 * 
 * The input `anab::T0`, on the other end, as produced by the Pandora pattern
 * recognition algorithm, and ultimately relative to the trigger time. In a
 * perfectly aligned setup, this is equivalent to the trigger time in real data.
 * 
 * The output of this module attempts to adhere to that convention, and it
 * defines the gate boundaries with respect to the beam gate time.
 * 
 */
class icarus::trigger::BeamGateInfoFromTracks: public art::SharedProducer {
  
    public:
  
  using nanoseconds = util::quantities::intervals::nanoseconds;
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    enum class GateType_t { // we need to translate enum into a strong type
        kUnknown = sim::kUnknown
      , kBNB     = sim::kBNB
      , kNuMI    = sim::kNuMI
    };
    
    /// Selector for `Type` parameter.
    static util::MultipleChoiceSelection<GateType_t> const GateTypeSelector;
    
    
    fhicl::Atom<art::InputTag> T0selProducer {
      Name("T0selProducer"),
      Comment
        ("tag of the selected particles (as a collection of art::Ptr)")
      // mandatory
      };
    
    fhicl::Atom<art::InputTag> T0Producer {
      Name("T0Producer"),
      Comment("tag of the input track time (t0) information")
      // mandatory
      };
    
    
    fhicl::Atom<nanoseconds> GateStartOffset {
      Name("GateStartOffset"),
      Comment("offset from time track to gate start")
      };

    fhicl::Atom<nanoseconds> GateEndOffset {
      Name("GateEndOffset"),
      Comment("offset from time track to gate end")
      };
    
    fhicl::Atom<std::string> GateType {
      Name("GateType"),
      Comment("beam gate type: " + GateTypeSelector.optionListString()),
      GateTypeSelector.get(GateType_t::kUnknown).name()
      };
    
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "BeamGateInfoFromTracks" // default
      };
    
    
    sim::BeamType_t getGateType() const
      {
        try {
          return static_cast<sim::BeamType_t>
            (GateTypeSelector.parse(GateType()).value());
        }
        catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e)
        {
          throw art::Exception(art::errors::Configuration)
            << "Invalid value for '" << GateType.name()
            << "' parameter: '" << e.label() << "'; valid options: "
            << GateTypeSelector.optionListString() << ".\n";
        }
      } // getGateType()
  
  
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit BeamGateInfoFromTracks
    (Parameters const& config, art::ProcessingFrame const&);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Fills the plots. Also extracts the information to fill them with.
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fT0selProducer; ///< Input particles.
  art::InputTag const fT0Producer; ///< Input particle/time associations.
  
  /// Offset of gate start from particle time.
  nanoseconds const fGateStartOffset;
  nanoseconds const fGateDuration; ///< Width of the gate being created.
  
  sim::BeamType_t const fBeamGateType; ///< Type of gate saved.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
}; // class icarus::trigger::BeamGateInfoFromTracks


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& data)
    { return std::make_unique<T>(std::move(data)); }

} // local namespace


//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  util::MultipleChoiceSelection<BeamGateInfoFromTracks::Config::GateType_t> const
  BeamGateInfoFromTracks::Config::GateTypeSelector
    {
        { GateType_t::kUnknown, "unknown" }
      , { GateType_t::kBNB,     "BNB" }
      , { GateType_t::kNuMI,    "NuMI" }
    };
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
icarus::trigger::BeamGateInfoFromTracks::BeamGateInfoFromTracks
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedProducer(config)
  // configuration
  , fT0selProducer(config().T0selProducer())
  , fT0Producer(config().T0Producer())
  , fGateStartOffset(config().GateStartOffset())
  , fGateDuration(config().GateEndOffset() - fGateStartOffset)
  , fBeamGateType(config().getGateType())
  , fLogCategory(config().LogCategory())
{
  
  async<art::InEvent>();
  
  //
  // output data declaration
  //
  produces<std::vector<sim::BeamGateInfo>>();
  //produces<art::Assns<sim::BeamGateInfo, recob::PFParticle>>();
  produces<art::Assns<sim::BeamGateInfo, recob::Track>>();
  produces<art::Assns<sim::BeamGateInfo, anab::T0>>();
  
  //
  // configuration report (short)
  //
  auto const& beamType
    = Config::GateTypeSelector.get(Config::GateType_t{ fBeamGateType });
  
  mf::LogInfo{ fLogCategory }
    << "Configuration:"
    << "\n - particle selection: '" << fT0selProducer.encode() << '\''
    << "\n - associated times: '" << fT0Producer.encode() << '\''
    << "\n - gate around particle time: " << fGateStartOffset
    << " -- " << (fGateStartOffset+fGateDuration)
    << " (" << fGateDuration << "; gate type: \"" << beamType.name()
    << "\" [#" << fBeamGateType << "])"
    ;
  
} // icarus::trigger::BeamGateInfoFromTracks::BeamGateInfoFromTracks()


//------------------------------------------------------------------------------
void icarus::trigger::BeamGateInfoFromTracks::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // fetch input
  //
  //auto const& particles
  //= event.getProduct<std::vector<art::Ptr<recob::PFParticle>>>(fT0selProducer);

  auto const& tracks                                                                                             
    = event.getProduct<std::vector<art::Ptr<recob::Track>>>(fT0selProducer);  

  mf::LogDebug(fLogCategory)
    << "Writing gates for " << tracks.size() << " tracks.";

  art::FindOneP<anab::T0> const trackT0(tracks, event, fT0Producer);
  
  detinfo::DetectorTimings const detTimings {
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  
  // add this offset to a time vs. trigger to make it vs. beam gate
  nanoseconds const triggerToBeamGate
    = detTimings.TriggerTime() - detTimings.BeamGateTime();
  
  //
  // create the content
  //
  std::vector<sim::BeamGateInfo> gates;
  art::Assns<sim::BeamGateInfo, recob::Track> gateToParticle;
  art::Assns<sim::BeamGateInfo, anab::T0> gateToTime;
  
  art::PtrMaker<sim::BeamGateInfo> const makeGatePtr { event };
  
  for (auto const& [ iTracks, trackPtr ]: util::enumerate(tracks)) {
    
    art::Ptr<anab::T0> const t0Ptr = trackT0.at(iTracks);
    if (t0Ptr.isNull()) {
      art::Exception e { art::errors::NotFound };
      e << "Selected track #" << iTracks << " (";
      if (trackPtr) e << "ID=" << trackPtr->ID() << ", ";
      e << "#" << trackPtr.key()
        << " in its collection) has no associated time."
        ;
      throw e << '\n';
    } // if no T0
    
    // t0 is stored in trigger time
    detinfo::timescales::trigger_time const t0
      { util::quantities::points::nanosecond(t0Ptr->Time()) };
    
    // 1DetectorTimings1 does not handle the beam gate time when converting to
    // simulation time, so I need to explicitly add the difference
    detinfo::timescales::simulation_time const gateStart
      = detTimings.toSimulationTime(t0 + triggerToBeamGate + fGateStartOffset);
    
    // time conversion is currently redundant here, left as documentation
    sim::BeamGateInfo gate {
        gateStart.value()                                 // start
      , fGateDuration.convertInto<nanoseconds>().value()  // width
      , fBeamGateType                                     // type
      };
    
    {
      mf::LogTrace log{ fLogCategory };
      log << "Gate for selected track #" << iTracks << " (";
      if (trackPtr) log << "ID=" << trackPtr->ID() << ", ";
      log << "#" << trackPtr.key()
        << " in its collection): time = " << t0 << " => "
        << gate.Start() << " -- " << (gate.Start() + gate.Width())
        << " ns (" << gate.Width() << " ns)"
        ;
    }
    
    gates.push_back(std::move(gate));
    
    art::Ptr<sim::BeamGateInfo> gatePtr { makeGatePtr(gates.size() - 1) };
    gateToParticle.addSingle(gatePtr, trackPtr);
    gateToTime.addSingle(gatePtr, t0Ptr);
    
  } // for
  
  //
  // store output
  //
  event.put(moveToUniquePtr(gates));
  event.put(moveToUniquePtr(gateToParticle));
  event.put(moveToUniquePtr(gateToTime));
  
} // icarus::trigger::BeamGateInfoFromTracks::produce()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::BeamGateInfoFromTracks)


//------------------------------------------------------------------------------
