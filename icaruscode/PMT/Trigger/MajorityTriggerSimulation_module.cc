/**
 * @file   MajorityTriggerSimulation_module.cc
 * @brief  Plots of efficiency for triggers based on PMT channel global count.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 9, 2020
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/ApplyBeamGate.h"
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // FillTriggerGates()
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h"
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icarusalg/Utilities/ROOTutils.h" // util::ROOT
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()...
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "icarusalg/Utilities/ChangeMonitor.h" // ThreadSafeChangeMonitor
#include "icarusalg/Utilities/rounding.h" // icarus::ns::util::roundup()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_tick...
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds, ...
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/get_elements.h" // util::get_elements()
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger
#include "lardataobj/RawData/OpDetWaveform.h" // raw::ADC_Count_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::CryostatID

// framework libraries
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TEfficiency.h"
#include "TH1F.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <ostream>
#include <algorithm> // std::fill()
#include <map>
#include <vector>
#include <memory> // std::make_unique()
#include <string>
#include <atomic>
#include <optional>
#ifdef __cpp_lib_source_location
#  include <source_location>
#endif //  __cpp_lib_source_location
#include <utility> // std::pair<>, std::move()
#include <cmath> // std::ceil()
#include <cstddef> // std::size_t
#include <cassert>


//------------------------------------------------------------------------------
using namespace util::quantities::time_literals;


// TODO Sort this mess

//------------------------------------------------------------------------------
namespace icarus::trigger { class MajorityTriggerCombiner; }

/// Combines a group of trigger gates for majority trigger. Glorified `Sum()`.
class icarus::trigger::MajorityTriggerCombiner
  : protected icarus::ns::util::mfLoggingClass
{
  
    public:
  
  MajorityTriggerCombiner
    (std::string const& logCategory = "MajorityTriggerCombiner")
    : icarus::ns::util::mfLoggingClass(logCategory) {}
  
  /// Combines all the gates (by cryostat) in a single majority gate.
  template <typename GateObj>
  GateObj combine(std::vector<GateObj> const& gates) const
    { return icarus::trigger::sumGates(gates); }

  
    private:
  
}; // class icarus::trigger::MajorityTriggerCombiner


//------------------------------------------------------------------------------
namespace icarus::trigger { class CryostatTriggerCombiner; }

/// Combines cryostat triggers via OR. Glorified `Max()`.
class icarus::trigger::CryostatTriggerCombiner
  : protected icarus::ns::util::mfLoggingClass
{
  
    public:
  
  CryostatTriggerCombiner
    (std::string const& logCategory = "CryostatTriggerCombiner")
    : icarus::ns::util::mfLoggingClass(logCategory) {}
  
  /// Combines all the gates (by cryostat) in a single majority gate.
  template <typename GateObj>
  GateObj combine(std::vector<GateObj> const& cryoGates) const
    { return icarus::trigger::maxGates(cryoGates); }
  
    private:
  
}; // class icarus::trigger::CryostatTriggerCombiner


//------------------------------------------------------------------------------
namespace icarus::trigger { class GeometryChannelSplitter; }


/// Combines cryostat triggers via OR. Glorified `Max()`.
class icarus::trigger::GeometryChannelSplitter
  : protected icarus::ns::util::mfLoggingClass
{
  
    public:
  
  GeometryChannelSplitter(
    geo::GeometryCore const& geom,
    std::string const& logCategory = "GeometryChannelSplitter"
    );
  
  /// Splits the gates by cryostat.
  template <typename GateObj>
  std::vector<std::vector<GateObj>> byCryostat
    (std::vector<GateObj>&& gates) const;
  
    private:
  
  unsigned int const fNCryostats; ///< Number of cryostats in the detector.
  
  /// Map: optical channel ID -> number of the cryostat with that channel.
  std::vector<geo::CryostatID> const fChannelCryostat;
  
  
  /// Creates a map like `fChannelCryostat` from the geometry information.
  static std::vector<geo::CryostatID> makeChannelCryostatMap
    (geo::GeometryCore const& geom);
  
}; // class icarus::trigger::GeometryChannelSplitter


//------------------------------------------------------------------------------
icarus::trigger::GeometryChannelSplitter::GeometryChannelSplitter(
  geo::GeometryCore const& geom,
  std::string const& logCategory /* = "GeometryChannelSplitter" */
  )
  : icarus::ns::util::mfLoggingClass(logCategory)
  , fNCryostats(geom.Ncryostats())
  , fChannelCryostat(makeChannelCryostatMap(geom))
{}


//------------------------------------------------------------------------------
template <typename GateObj>
std::vector<std::vector<GateObj>>
icarus::trigger::GeometryChannelSplitter::byCryostat
  (std::vector<GateObj>&& gates) const
{
  std::vector<std::vector<GateObj>> gatesPerCryostat{ fNCryostats };
  
  for (auto& gate: gates) {
    assert(gate.hasChannels());
    gatesPerCryostat[fChannelCryostat.at(gate.channels().front()).Cryostat]
      .push_back(std::move(gate));
  } // for gates
  
  return gatesPerCryostat;
} // icarus::trigger::GeometryChannelSplitter::byCryostat()


//------------------------------------------------------------------------------
auto icarus::trigger::GeometryChannelSplitter::makeChannelCryostatMap
  (geo::GeometryCore const& geom) -> std::vector<geo::CryostatID>
{
  
  auto const nOpChannels = geom.NOpChannels();
  
  std::vector<geo::CryostatID> channelCryostatMap(nOpChannels);
  
  for (auto const opChannel: util::counter(nOpChannels)) {
    if (!geom.IsValidOpChannel(opChannel)) continue;
    channelCryostatMap.at(opChannel)
      = geom.OpDetGeoFromOpChannel(opChannel).ID();
  } // for all channels
  
  return channelCryostatMap;
  
} // icarus::trigger::GeometryChannelSplitter::makeChannelCryostatMap()


//------------------------------------------------------------------------------
namespace icarus::trigger { class MajorityTriggerSimulation; }
/**
 * @brief Simulates a "majority" trigger.
 * 
 * A trigger primitive is a two-level function of time which describes when
 * that primitive is on and when it is off. Trigger primitives are given as
 * input to this module and their origin may vary, but the standard source in
 * ICARUS is @ref ICARUSPMTTriggerGlossary "single trigger request".
 * 
 * This module simulates a trigger requesting a minimum number of single trigger
 * requests in the event, and saves the result as `raw::Trigger` data products.
 * While only one trigger definition is used, inputs with different thresholds
 * may be specified to have the different responses.
 * 
 * 
 * Configuration
 * ==============
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the trigger primitives to be used as input;
 *     it must not include any instance name, as the instance names will be
 *     automatically added from `Thresholds` parameter.
 *     The typical trigger primitives used as input may be LVDS discriminated
 *     output (e.g. from `icarus::trigger::LVDSgates` module) or combinations
 *     of them (e.g. from `icarus::trigger::SlidingWindowTrigger` module).
 * * `Thresholds` (list of integers, mandatory): list of the discrimination
 *     thresholds to consider, in ADC counts. A data product containing a
 *     digital signal is read for each one of the thresholds, and the tag of the
 *     data product is expected to be the module label `TriggerGatesTag` with as
 *     instance name the value of the threshold (e.g. for a threshold of 60 ADC
 *     counts the data product tag might be `LVDSgates:60`).
 * * `MinimumPrimitives` (integer, _mandatory_): the required number of
 *   single trigger requests in order for the trigger to fire;
 * * `BeamGateDuration` (time, _mandatory_): the duration of the beam
 *     gate; _the time requires the unit to be explicitly specified_: use
 *     `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *     `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *     `trigger_icarus.fcl`);
 * * `TriggerTimeResolution` (time, default: `8 ns`): time resolution for the
 *     trigger primitives;
 * * `LogCategory` (string, default `TriggerEfficiencyPlots`): name of category
 *     used to stream messages from this module into message facility.
 * 
 * An example job configuration is provided as
 * `simulatemajoritytriggers_icarus.fcl`.
 * 
 * 
 * Output data products
 * =====================
 * 
 * * `std::vector<raw::Trigger>` (one instance per ADC threshold):
 *   list of triggers fired according to the configured trigger definition;
 *   there is one collection (and data product) per ADC threshold, and the
 *   data product has the same instance name as the input data one
 *   (see `TriggerGatesTag` and `Thresholds` configuration parameters);
 *   currently only at most one trigger is emitted, with time stamp matching
 *   the first time the trigger criteria are satisfied.
 * 
 * 
 * 
 * Trigger logic algorithm
 * ========================
 * 
 * @anchor MajorityTriggerSimulation_Algorithm
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::MajorityTriggerSimulation` and its assumptions.
 * 
 * The algorithm keeps the trigger primitives from the different cryostats
 * separate for the most time.
 * Within each cryostat, all trigger primitives are treated equally, whether
 * they originate from one or from two channels (or 10 or 30), and wherever
 * their channels are in the cryostat.
 * The trigger primitives in each cryostat are combined in a multi-level gate by
 * adding them, so that the level of the resulting gate matches at any time how
 * many trigger primitives are on at that time.
 * Finally, the maximum number of trigger primitives open in any of the
 * cryostats at each time is the level to be compared to the trigger
 * requirements.
 * 
 * This multi-level gate is set in coincidence with the beam gate by multiplying
 * the multi-level and the beam gates.
 * The beam gate opens at a time configured in `DetectorClocks` service provider
 * (`detinfo::DetectorClocks::BeamGateTime()`) and has a duration configured
 * in this module (`BeamGateDuration`).
 * 
 * At this point, the trigger gate is a multi-level gate suppressed everywhere
 * except than during the beam gate.
 * The trigger requirement is simply how many trigger primitives must be open at
 * the same time in a single cryostat for the trigger to fire. The requirement
 * is set in the configuration (`MinimumPrimitives`).
 * To determine whether a trigger with the given minimum number of primitives
 * open at the same time has fired, the gate combined as described above is
 * scanned to find _the first tick_ where the level of the gate reaches
 * or passes this minimum required number. If such tick exists, the trigger is
 * considered to have fired, and at that time.
 * 
 * While there _is_ a parameter describing the time resolution of the trigger
 * (`TriggerTimeResolution`), this is currently only used for aesthetic purposes
 * to choose the binning of some plots: the resolution is _not_ superimposed
 * to the gates (yet).
 * 
 * 
 * Output plots
 * -------------
 * 
 * Plots are directly stored in the producer folder of the `TFileService` ROOT
 * output file.
 * 
 * Summary plots are generated:
 * * `NTriggers`: total number of triggers, per threshold
 * * `Eff`: fraction of events with at least one trigger, per threshold
 * * `TriggerTick`: trigger time distribution (optical tick), per threshold
 * 
 * The plots marked "per threshold" have a "threshold" axis where each bin
 * corresponds to the threshold specified in the `Thresholds` configuration
 * parameter. Note that the numerical value of the axis on that bin does not
 * match the threshold value.
 * 
 */
class icarus::trigger::MajorityTriggerSimulation
  : public art::EDProducer
  , private lar::UncopiableAndUnmovableClass
{

    public:
  
  using microseconds = util::quantities::intervals::microseconds;
  using nanoseconds = util::quantities::intervals::nanoseconds;
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Sequence<raw::ADC_Count_t> Thresholds {
      Name("Thresholds"),
      Comment("thresholds to consider [ADC counts]")
      };

    fhicl::Atom<unsigned int> MinimumPrimitives {
      Name("MinimumPrimitives"),
      Comment("minimum required number of trigger primitives for the trigger")
      };
    
    fhicl::Atom<microseconds> BeamGateDuration {
      Name("BeamGateDuration"),
      Comment("length of time interval when optical triggers are accepted")
      };

    fhicl::Atom<std::uint32_t> BeamBits {
      Name("BeamBits"),
      Comment("bits to be set in the trigger object as beam identified")
      };

    fhicl::Atom<nanoseconds> TriggerTimeResolution {
      Name("TriggerTimeResolution"),
      Comment("resolution of trigger in time"),
      8_ns
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "SlidingWindowTrigger" // default
      };
    
  }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit MajorityTriggerSimulation(Parameters const& config);
/*
  // Plugins should not be copied or assigned.
  MajorityTriggerSimulation(MajorityTriggerSimulation const&) = delete;
  MajorityTriggerSimulation(MajorityTriggerSimulation&&) = delete;
  MajorityTriggerSimulation& operator=(MajorityTriggerSimulation const&) = delete;
  MajorityTriggerSimulation& operator=(MajorityTriggerSimulation&&) = delete;
*/
  // --- END Constructors ------------------------------------------------------


  // --- BEGIN Framework hooks -------------------------------------------------

  /// Initializes the plots.
  virtual void beginJob() override;
  
  /// Runs the simulation and saves the results into the _art_ event.
  virtual void produce(art::Event& event) override;
  
  /// Prints end-of-job summaries.
  virtual void endJob() override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  using TriggerInfo_t = details::TriggerInfo_t; ///< Type alias.
  
  /// Type of trigger gate extracted from the input event.
  using InputTriggerGate_t = icarus::trigger::MultiChannelOpticalTriggerGate;
  
  /// A list of trigger gates from input.
  using TriggerGates_t = std::vector<InputTriggerGate_t>;

  /// Type of gate data without channel information.
  using TriggerGateData_t = InputTriggerGate_t::GateData_t;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<icarus::trigger::ADCCounts_t, art::InputTag> fADCthresholds;
  
  /// Minimum number of trigger primitives for a trigger to happen.
  unsigned int const fMinimumPrimitives;
  
  /// Duration of the gate during with global optical triggers are accepted.
  microseconds fBeamGateDuration;
  
  nanoseconds fTriggerTimeResolution; ///< Trigger resolution in time.
  
  std::uint32_t fBeamBits; ///< Bits for the beam gate being simulated.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  // --- BEGIN Service variables -----------------------------------------------

  geo::GeometryCore const& fGeom;

  /// ROOT directory where all the plots are written.
  art::TFileDirectory fOutputDir;

  // --- END Service variables -------------------------------------------------

  
  // --- BEGIN Internal variables ----------------------------------------------
  
  MajorityTriggerCombiner const fCombiner; ///< Algorithm to combine primitives.
  
  /// Algorithm to sort trigger gates by cryostat or TPC.
  GeometryChannelSplitter fChannelSplitter;
  
  /// All plots in one practical sandbox.
  icarus::trigger::PlotSandbox fPlots;

  ///< Count of fired triggers, per threshold.
  std::vector<std::atomic<unsigned int>> fTriggerCount;
  std::atomic<unsigned int> fTotalEvents { 0U }; ///< Count of fired triggers.
  
  
  // TODO this is not multithread-safe, needs a mutex
  /// Functor returning whether a gate has changed.
  icarus::ns::util::ThreadSafeChangeMonitor<icarus::trigger::ApplyBeamGateClass>
    fGateChangeCheck;

  // --- END Internal variables ------------------------------------------------
  

  
  // --- BEGIN Derived class methods -------------------------------------------
  
  /// @brief Initializes the full set of plots (all ADC thresholds).
  void initializePlots();
  
  /**
   * @brief Performs the simulation for the specified ADC threshold.
   * @param event _art_ event to read data from and put results into
   * @param iThr index of the threshold in the configuration
   * @param thr value of the threshold (ADC counts)
   * @return a simple copy of the trigger response information
   * 
   * For the given threshold, the simulation of the configured trigger is
   * performed.
   * The input data is read from the event (the source tag is from the module
   * configuration), simulation is performed, auxiliary plots are drawn and
   * a `raw::Trigger` collection is stored into the event.
   * 
   * The stored collection contains either one or zero `raw::Trigger` elements.
   * 
   * The simulation itself is performed by the `simulate()` method.
   */
  TriggerInfo_t produceForThreshold(
    art::Event& event,
    detinfo::DetectorTimings const& detTimings,
    ApplyBeamGateClass const& beamGate,
    std::size_t const iThr, icarus::trigger::ADCCounts_t const thr
    );
  
  /**
   * @brief Performs the simulation of the configured trigger on `gates` input.
   * @param gates the input to the trigger simulation
   * @return the outcome and details of the trigger simulation
   * 
   * The simulation is performed using the input single trigger requests
   * (or trigger primitives) from the `gates` collection.
   * 
   * The gates are split by cryostat, and the simulation is performed
   * independently on each cryostat (`simulateCryostat()`).
   * Finally, the cryostat triggers are combined (OR) into the final trigger
   * decision, bearing as time the earliest one.
   */
  TriggerInfo_t simulate(ApplyBeamGateClass const& clockData,
                         TriggerGates_t const& gates) const;
  
  /**
   * @brief Simulates the trigger response within a single cryostat.
   * @param gates the trigger primitives to be considered
   * @return the outcome and details of the trigger simulation
   * 
   * The simulation computes the count of trigger `gates` open at any time,
   * sets it in coincidence with the beam gate, and fires a trigger if within
   * that gate the count of open `gates` is equal or larger than the threshold
   * configured (`MinimumPrimitives`).
   * The time is the earliest one when that requirement is met.
   */
  TriggerInfo_t simulateCryostat(ApplyBeamGateClass const& clockData,
                                 TriggerGates_t const& gates) const;
  
  
  /**
   * @brief Converts the trigger information into a `raw::Trigger` object.
   * @param triggerNumber the unique number to assign to this trigger
   * @param info the information about the fired trigger
   * @return a `raw::Trigger` object with all the information encoded
   * 
   * The trigger described by `info` is encoded into a `raw::Trigger` object.
   * The trigger _must_ have fired.
   */
  raw::Trigger triggerInfoToTriggerData
    (detinfo::DetectorTimings const& detTimings,
     unsigned int triggerNumber, TriggerInfo_t const& info) const;
  
  /// Fills the plots for threshold index `iThr` with trigger information.
  void plotTriggerResponse
    (std::size_t iThr, TriggerInfo_t const& triggerInfo) const;
  
  
  /// Prints the summary of fired triggers on screen.
  void printSummary() const;
  
  
  //@{ 
  /// Shortcut to create an `ApplyBeamGate` with the current configuration.
  icarus::trigger::ApplyBeamGateClass makeMyBeamGate
    (art::Event const* event = nullptr) const
    {
      return makeApplyBeamGate(
        fBeamGateDuration,
        icarus::ns::util::makeDetClockData(event),
        fLogCategory
        );
    }
  icarus::trigger::ApplyBeamGateClass makeMyBeamGate
    (art::Event const& event) const { return makeMyBeamGate(&event); }
  //@}
  
  
  /// Reads a set of input gates from the `event`
  /// @return trigger gates, converted into `InputTriggerGate_t`
  static TriggerGates_t readTriggerGates
    (art::Event const& event, art::InputTag const& dataTag);
  

}; // icarus::trigger::MajorityTriggerSimulation



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
icarus::trigger::MajorityTriggerSimulation::MajorityTriggerSimulation
  (Parameters const& config)
  : art::EDProducer       (config)
  // configuration
  , fMinimumPrimitives    (config().MinimumPrimitives())
  , fBeamGateDuration     (config().BeamGateDuration())
  , fTriggerTimeResolution(config().TriggerTimeResolution())
  , fBeamBits             (config().BeamBits())
  , fLogCategory          (config().LogCategory())
  // services
  , fGeom      (*lar::providerFrom<geo::Geometry>())
  , fOutputDir (*art::ServiceHandle<art::TFileService>())
  // internal and cached
  , fCombiner       (fLogCategory)
  , fChannelSplitter(fGeom, fLogCategory)
  , fPlots(
     fOutputDir, "", "minimum primitives: " + std::to_string(fMinimumPrimitives)
    )
{
  
  //
  // more complex parameter parsing
  //
  std::string const discrModuleLabel = config().TriggerGatesTag();
  for (raw::ADC_Count_t threshold: config().Thresholds()) {
    fADCthresholds[icarus::trigger::ADCCounts_t{threshold}]
      = art::InputTag{ discrModuleLabel, util::to_string(threshold) };
  }
  
  // initialization of a vector of atomic is not as trivial as it sounds...
  fTriggerCount = std::vector<std::atomic<unsigned int>>(fADCthresholds.size());
  std::fill(fTriggerCount.begin(), fTriggerCount.end(), 0U);


  //
  // input data declaration
  //
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // trigger primitives
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    consumes<std::vector<OpticalTriggerGateData_t>>(inputDataTag);
    consumes<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag);
  } // for
  
  //
  // output data declaration
  //
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds))
    produces<std::vector<raw::Trigger>>(inputDataTag.instance());
  
  
  mf::LogInfo(fLogCategory)
    << "Requirement of minimum " << fMinimumPrimitives << " primitives.";
  
  
} // icarus::trigger::MajorityTriggerSimulation::MajorityTriggerSimulation()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerSimulation::beginJob() {
  
  initializePlots();
  
} // icarus::trigger::MajorityTriggerSimulation::beginJob()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerSimulation::produce(art::Event& event) {
  
  mf::LogDebug log(fLogCategory);
  log << "Event " << event.id() << ":";
  
  auto const clockData
    = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
  detinfo::DetectorTimings const detTimings{clockData};
  auto const beamGate = makeMyBeamGate(event);

  if (auto oldGate = fGateChangeCheck(beamGate); oldGate) {
    mf::LogWarning(fLogCategory)
      << "Beam gate has changed from " << *oldGate << " to " << beamGate << "!";
  }


  for (auto const [ iThr, thr ]
    : util::enumerate(util::get_elements<0U>(fADCthresholds))
  ) {
    
    TriggerInfo_t const triggerInfo = produceForThreshold(event, detTimings, beamGate, iThr, thr);
    
    log << "\n * threshold " << thr << ": ";
    if (triggerInfo) log << "trigger at " << triggerInfo.atTick();
    else             log << "not triggered";
    
  } // for
  
  ++fTotalEvents;
  
} // icarus::trigger::MajorityTriggerSimulation::produce()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerSimulation::endJob() {
  
  printSummary();
  
} // icarus::trigger::MajorityTriggerSimulation::endJob()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerSimulation::initializePlots() {
  
  //
  // overview plots with different settings
  //
  
  std::vector<std::string> thresholdLabels;
  thresholdLabels.reserve(size(fADCthresholds));
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds))
    thresholdLabels.push_back(inputDataTag.instance());
  
  auto const beamGate = makeMyBeamGate();
  fGateChangeCheck(beamGate);
  mf::LogInfo(fLogCategory)
    << "Beam gate for plots: " << beamGate.asSimulationTime()
    << " (simulation time), " << beamGate.tickRange()
    << " (optical ticks)"
    ;

  //
  // Triggering efficiency vs. ADC threshold.
  //
  auto* NTriggers = fPlots.make<TH1F>(
    "NTriggers",
    "Number of triggering events"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";events",
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );
  util::ROOT::applyAxisLabels(NTriggers->GetXaxis(), thresholdLabels);
  
  auto* Eff = fPlots.make<TEfficiency>(
    "Eff",
    "Efficiency of triggering"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";trigger efficiency",
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );
  // people are said to have earned hell for things like this;
  // but TEfficiency really does not expose the interface to assign labels to
  // its axes, which supposedly could be done had we chosen to create it by
  // histograms instead of directly as recommended.
  util::ROOT::applyAxisLabels
    (const_cast<TH1*>(Eff->GetTotalHistogram())->GetXaxis(), thresholdLabels);
  
  detinfo::timescales::optical_time_ticks const triggerResolutionTicks{
    icarus::ns::util::makeDetTimings().toOpticalTicks(fTriggerTimeResolution)
    };
  
  auto const& beamGateTicks = beamGate.tickRange();
  auto* TrigTime = fPlots.make<TH2F>(
    "TriggerTick",
    "Trigger time tick"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";events",
    static_cast<int>(std::ceil(beamGate.lengthTicks()/triggerResolutionTicks)),
    beamGateTicks.start().value(),
    icarus::ns::util::roundup
     (beamGateTicks.start() + beamGate.lengthTicks(), triggerResolutionTicks)
     .value(),
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );
  util::ROOT::applyAxisLabels(TrigTime->GetYaxis(), thresholdLabels);
  
  
} // icarus::trigger::MajorityTriggerSimulation::initializePlots()


//------------------------------------------------------------------------------
auto icarus::trigger::MajorityTriggerSimulation::produceForThreshold(
  art::Event& event,
  detinfo::DetectorTimings const& detTimings,
  ApplyBeamGateClass const& beamGate,
  std::size_t const iThr, icarus::trigger::ADCCounts_t const thr
) -> TriggerInfo_t {
  
  //
  // get the input
  //
  art::InputTag const dataTag = fADCthresholds.at(thr);
  auto const& gates = readTriggerGates(event, dataTag);
  
  //
  // simulate the trigger response
  //
  TriggerInfo_t const triggerInfo = simulate(beamGate, gates);
  if (triggerInfo) ++fTriggerCount[iThr]; // keep the unique count
  
  //
  // fill the plots
  //
  plotTriggerResponse(iThr, triggerInfo);

  //
  // create and store the data product
  //
  auto triggers = std::make_unique<std::vector<raw::Trigger>>();
  if (triggerInfo.fired()) {
    triggers->push_back
      (triggerInfoToTriggerData(detTimings, fTriggerCount[iThr], triggerInfo));
  } // if
  event.put(std::move(triggers), dataTag.instance());
  
  return triggerInfo;
  
} // icarus::trigger::MajorityTriggerSimulation::produceForThreshold()


//------------------------------------------------------------------------------
auto icarus::trigger::MajorityTriggerSimulation::simulate
  (ApplyBeamGateClass const& beamGate,
   TriggerGates_t const& gates) const -> TriggerInfo_t
{

  /* 
   * 1. split the input by cryostat
   * 2. simulate the cryostat trigger
   * 3. combine the responses (earliest wins)
   */
  
  // to use the splitter we need a *copy* of the gates
  auto const& cryoGates = fChannelSplitter.byCryostat(TriggerGates_t{ gates });
  
  // NOTE to allow for distinction between cryostats, the logic needs to be reworked
  TriggerInfo_t triggerInfo; // not fired by default
  for (auto const& gatesInCryo: cryoGates) {
    
    triggerInfo.replaceIfEarlier(simulateCryostat(beamGate, gatesInCryo));
    
  } // for gates in cryostat
  
  return triggerInfo;
  
} // icarus::trigger::MajorityTriggerSimulation::simulate()


//------------------------------------------------------------------------------
auto icarus::trigger::MajorityTriggerSimulation::simulateCryostat
  (ApplyBeamGateClass const& beamGate, TriggerGates_t const& gates) const
  -> TriggerInfo_t
{

  /* 
   * 1. combine the trigger primitives
   * 2. apply the beam gate on the combination
   * 3. compute the trigger response
   */
  
  auto const combinedCount = beamGate.apply(fCombiner.combine(gates));
  
  // the first tick with enough opened gates:
  TriggerGateData_t::ClockTick_t const tick
    = combinedCount.findOpen(fMinimumPrimitives);
  bool const fired = (tick != TriggerGateData_t::MaxTick);
  
  TriggerInfo_t triggerInfo;
  if (fired) triggerInfo.emplace(optical_tick{ tick });
  return triggerInfo;
  
} // icarus::trigger::MajorityTriggerSimulation::simulateCryostat()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerSimulation::plotTriggerResponse
  (std::size_t iThr, TriggerInfo_t const& triggerInfo) const
{
  
  fPlots.demand<TEfficiency>("Eff").Fill(triggerInfo.fired(), iThr);
  
  if (triggerInfo.fired()) {
    fPlots.demand<TH1>("NTriggers").Fill(iThr);
    fPlots.demand<TH2>("TriggerTick").Fill(triggerInfo.atTick().value(), iThr);
  }
  
} // icarus::trigger::MajorityTriggerSimulation::plotTriggerResponse()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerSimulation::printSummary() const {
  
  //
  // summary from our internal counters
  //
  mf::LogInfo log(fLogCategory);
  log
    << "Summary of triggers requiring " << fMinimumPrimitives
    << "+ primitives for " << fTriggerCount.size() << " ADC thresholds:"
    ;
  for (auto const& [ count, thr ]
    : util::zip(fTriggerCount, util::get_elements<0U>(fADCthresholds)))
  {
    log << "\n  ADC threshold " << thr
      << ": " << count << " events triggered";
    if (fTotalEvents > 0U)
      log << " (" << (double(count) / fTotalEvents * 100.0) << "%)";
  } // for
  
} // icarus::trigger::MajorityTriggerSimulation::printSummary()


//------------------------------------------------------------------------------
raw::Trigger
icarus::trigger::MajorityTriggerSimulation::triggerInfoToTriggerData
  (detinfo::DetectorTimings const& detTimings,
   unsigned int triggerNumber, TriggerInfo_t const& info) const
{
  assert(info.fired());
  
  return {
    triggerNumber,                                        // counter
    double(detTimings.toElectronicsTime(info.atTick())), // trigger time
    double(detTimings.BeamGateTime()), // beam gate in electronics time scale
    fBeamBits                                             // bits 
    };
  
} // icarus::trigger::MajorityTriggerSimulation::triggerInfoToTriggerData()


//------------------------------------------------------------------------------
auto icarus::trigger::MajorityTriggerSimulation::readTriggerGates
  (art::Event const& event, art::InputTag const& dataTag)
  -> TriggerGates_t
{

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // currently the associations are a waste of time memory...
  auto const& gates
    = *(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>(dataTag));
  auto const& gateToWaveforms = *(
    event.getValidHandle
      <art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>(dataTag)
    );
  
  try {
    return icarus::trigger::FillTriggerGates<InputTriggerGate_t>
      (gates, gateToWaveforms);
  }
  catch (cet::exception const& e) {
    throw cet::exception("MajorityTriggerSimulation", "", e)
      << "Error encountered while reading data products from '"
      << dataTag.encode() << "'\n";
  }

} // icarus::trigger::MajorityTriggerSimulation::readTriggerGates()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::MajorityTriggerSimulation)


//------------------------------------------------------------------------------
