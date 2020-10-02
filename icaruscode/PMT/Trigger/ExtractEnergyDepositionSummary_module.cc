/**
 * @file   ExtractEnergyDepositionSummary_module.cc
 * @brief  Module producing discriminated waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 4, 2019
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateStruct.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h"
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()...
#include "icaruscode/Utilities/ChangeMonitor.h" // ThreadSafeChangeMonitor
#include "icaruscode/IcarusObj/SimEnergyDepositSummary.h"
// #include "icaruscode/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
/*
#include "art/Framework/Principal/Handle.h"
*/
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <cassert>
/*
#include <map>
*/

//------------------------------------------------------------------------------
namespace icarus::trigger { class ExtractEnergyDepositionSummary; }

/**
 * @brief Produces energy deposition summary data products.
 * 
 * This data product is an alternatice to storing the whole energy deposition,
 * which can be huge.
 * 
 * Thresholds are ultimately chosen by the tool in charge of actually run the
 * discrimination algorithm. Out of the thresholds that this algorithm produces,
 * it is possible to choose only a subset of them (`SelectThresholds`).
 * 
 * 
 * Output data products
 * =====================
 * 
 * * a single `icarus::SimEnergyDepositSummary`
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<sim::SimEnergyDeposits>`: sums of energies in a spill time,
 *   in a pre-spill time, and total, everywhere and in active volume only.
 * 
 * 
 * Service requirements
 * ---------------------
 * 
 * The following services are _required_:
 * 
 * * `Geometry` for the determination of the active volume(s)
 * * `DetectorClocksService` for the determination of beam gate
 * * _art_ message facility
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description ExtractEnergyDepositionSummary`.
 * 
 * * `EnergyDepositTags`
 *     (list of input tags, default: `[ "largeant:TPCActive" ]`): a list of
 *     data products with energy depositions;
 * * `BeamGateDuration` (time, _mandatory_): the duration of the beam
 *     gate; _the time requires the unit to be explicitly specified_: use
 *     `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *     `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *     `trigger_icarus.fcl`);
 * * `BeamGateStart` (time, default: 0 &micro;s): open the beam gate this long
 *      after the nominal beam gate time;
 * * `PreSpillWindow` (time, default: 10 &micro;s): duration of the pre-spill
 *      window;
 * * `PreSpillWindowGap` (time, default: 0 &micro;s): gap from the end of
 *      pre-spill window to the start of beam gate;
 * * `OutputCategory` (string, default: `"ExtractEnergyDepositionSummary"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 */
class icarus::trigger::ExtractEnergyDepositionSummary: public art::EDProducer {
  
    public:
  
  using microseconds = util::quantities::intervals::microseconds;
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> EnergyDepositTags {
      Name("EnergyDepositTags"),
      Comment("label of energy deposition data product(s) in the detector"),
      std::vector<art::InputTag>{ "largeant:TPCActive" }
      };

    fhicl::Atom<microseconds> BeamGateDuration {
      Name("BeamGateDuration"),
      Comment("length of time interval when optical triggers are accepted")
      };

    fhicl::Atom<microseconds> BeamGateStart {
      Name("BeamGateStart"),
      Comment("open the beam gate this long after the nominal beam gate time"),
      microseconds{ 0.0 }
      };

    fhicl::Atom<microseconds> PreSpillWindow {
      Name("PreSpillWindow"),
      Comment("duration of the pre-spill window"),
      microseconds{ 10.0 }
      };

    fhicl::Atom<microseconds> PreSpillWindowGap {
      Name("PreSpillWindowGap"),
      Comment("gap from the end of pre-spill window to the start of beam gate"),
      microseconds{ 0.0 }
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("tag of the module output to console via message facility"),
      "ExtractEnergyDepositionSummary"
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit ExtractEnergyDepositionSummary(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  ExtractEnergyDepositionSummary(ExtractEnergyDepositionSummary const&) = delete;
  ExtractEnergyDepositionSummary(ExtractEnergyDepositionSummary&&) = delete;
  ExtractEnergyDepositionSummary& operator=(ExtractEnergyDepositionSummary const&) = delete;
  ExtractEnergyDepositionSummary& operator=(ExtractEnergyDepositionSummary&&) = delete;
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Duration of the gate during with global optical triggers are accepted.
  microseconds fBeamGateDuration;
  
  /// Start of the beam gate with respect to `BeamGate()`.
  microseconds fBeamGateStart;
  
  microseconds fPreSpillWindow; ///< Duration of the pre-spill gate.
  
  microseconds fPreSpillStart; ///< Start of the pre-spill gate.
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  geo::GeometryCore const& fGeom; ///< Access to detector geometry information.
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  /// Functor returning whether a gate has changed.
  icarus::ns::util::ThreadSafeChangeMonitor<icarus::trigger::BeamGateStruct>
    fBeamGateChangeCheck;

  details::EventInfoExtractorMaker const fEventInfoExtractorMaker;
  
  // --- END Algorithms --------------------------------------------------------
  
  
  
}; // icarus::trigger::ExtractEnergyDepositionSummary



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::ExtractEnergyDepositionSummary
//------------------------------------------------------------------------------
icarus::trigger::ExtractEnergyDepositionSummary::ExtractEnergyDepositionSummary
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fBeamGateDuration     (config().BeamGateDuration())
  , fBeamGateStart        (config().BeamGateStart())
  , fPreSpillWindow       (config().PreSpillWindow())
  , fPreSpillStart
      (fBeamGateStart - config().PreSpillWindowGap() - fPreSpillWindow)
  , fLogCategory(config().OutputCategory())
  // services
  , fGeom(*(lar::providerFrom<geo::Geometry>()))
  // algorithms
  , fEventInfoExtractorMaker(
    {},                                    // truthTags: none
    { config().EnergyDepositTags() },      // edepTags
    fGeom,                                 // geom
    fLogCategory,                          // logCategory
    consumesCollector()                    // consumesCollector
    )
{
  produces<icarus::SimEnergyDepositSummary>();
} // icarus::trigger::ExtractEnergyDepositionSummary::ExtractEnergyDepositionSummary()


//------------------------------------------------------------------------------
void icarus::trigger::ExtractEnergyDepositionSummary::produce(art::Event& event)
{
  
  //
  // prepare the information extractor and extract the information
  //
  
  // we need to convert the two relevant gates with the proper parameters
  auto const detTimings = icarus::ns::util::makeDetTimings(event);
  
  auto const beamGate = icarus::trigger::makeBeamGateStruct
    (detTimings, fBeamGateDuration, fBeamGateStart);
  
  if (auto oldGate = fBeamGateChangeCheck(beamGate); oldGate) {
    mf::LogWarning(fLogCategory)
      << "Beam gate has changed from " << oldGate->asOptTickRange()
      << " to " << beamGate.asOptTickRange() << " (optical tick)!";
  }
  
  details::EventInfo_t const info = fEventInfoExtractorMaker(
    beamGate.asSimulationRange(),
    icarus::trigger::makeBeamGateStruct
      (detTimings, fPreSpillStart, fPreSpillStart + fPreSpillWindow)
      .asSimulationRange()
    )(event);
  assert(info.hasDepEnergy());
  
  //
  // move the information into the data product
  //
  icarus::SimEnergyDepositSummary summary;
  summary.Total          = info.DepositedEnergy().value();
  summary.Spill          = info.DepositedEnergyInSpill().value();
  summary.PreSpill       = info.DepositedEnergyInPreSpill().value();
  summary.Active         = info.DepositedEnergyInActiveVolume().value();
  summary.SpillActive    = info.DepositedEnergyInSpillInActiveVolume().value();
  summary.PreSpillActive = info.DepositedEnergyInPreSpillInActiveVolume().value();

  //
  // finish
  //
  event.put
    (std::make_unique<icarus::SimEnergyDepositSummary>(std::move(summary)));
  
} // icarus::trigger::ExtractEnergyDepositionSummary::produce()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::ExtractEnergyDepositionSummary)


//------------------------------------------------------------------------------
