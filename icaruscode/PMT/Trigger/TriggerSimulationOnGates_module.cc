/**
 * @file   TriggerSimulationOnGates_module.cc
 * @brief  Plots of efficiency for triggers based on PMT sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 27, 2021
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.h"
#include "icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.h" // WindowTopologyManager
#include "icaruscode/PMT/Trigger/Algorithms/WindowPatternConfig.h"
#include "icaruscode/PMT/Trigger/Algorithms/WindowPattern.h"
#include "icaruscode/PMT/Trigger/Algorithms/ApplyBeamGate.h"
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // FillTriggerGates()
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icarusalg/Utilities/ROOTutils.h" // util::ROOT
#include "icarusalg/Utilities/BinningSpecs.h"
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()...
#include "icarusalg/Utilities/FixedBins.h"
#include "icarusalg/Utilities/PassCounter.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "icarusalg/Utilities/ChangeMonitor.h" // ThreadSafeChangeMonitor
#include "icarusalg/Utilities/rounding.h" // icarus::ns::util::roundup()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
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
#include "lardataobj/Simulation/BeamGateInfo.h"
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
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TEfficiency.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

// C/C++ standard libraries
#include <ostream>
#include <algorithm> // std::fill()
#include <map>
#include <vector>
#include <memory> // std::make_unique()
#include <string>
#include <atomic>
#include <optional>
#include <utility> // std::pair<>, std::move()
#include <cmath> // std::ceil()
#include <cstddef> // std::size_t
#include <cassert>


//------------------------------------------------------------------------------
using namespace util::quantities::time_literals;


//------------------------------------------------------------------------------
namespace icarus::trigger { class TriggerSimulationOnGates; }
/**
 * @brief Simulates a sliding window trigger at specified gate times.
 * 
 * This module produces `raw::Trigger` objects each representing the outcome of
 * some trigger logic applied to discriminated optical detector input
 * ("trigger primitives").
 * The logic is applied to each event at multiple times, according to a list
 * of time intervals read from each event.
 * 
 * The main purpose of this module is to simulate the trigger logic at times of
 * special interest, typically the times some track is believed to have crossed
 * the detector.
 * 
 * A trigger primitive is a two-level function of time which describes when
 * that primitive is on and when it is off. Trigger primitives are given as
 * input to this module and their origin may vary, but the standard source in
 * ICARUS is @ref ICARUSPMTTriggerGlossary "single trigger request (LVDS)".
 * 
 * This module applies a sliding window pattern to the input: the pattern
 * consists of a requirement on the main window and optional additional
 * requirements on the neighbouring windows. This module rebases the configured
 * pattern on each of the available windows, evaluates the requirement of the
 * pattern in that configuration, and decides whether those requirements are
 * met. The general trigger is considered passed if _any_ of the rebased
 * patterns satisfies the requirement at any time, and no special treatment is
 * performed in case multiple windows fulfil them, except that the trigger time
 * is driven by the earliest of the satisfied patterns.
 * 
 * A single trigger pattern is configured for each instance of the module,
 * while multiple input sets (e.g. with different discrimination thresholds)
 * can be processed on the same pattern by the same module instance.
 * Conversely, testing a different pattern requires the instantiation of a new
 * module.
 * 
 * 
 * Configuration
 * ==============
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module instance which
 *     produced the trigger primitives to be used as input; it must not include
 *     any instance name, as the instance names will be automatically added from
 *     `Thresholds` parameter.
 *     The typical trigger primitives used as input are LVDS discriminated
 *     output combined into trigger windows (e.g. from
 *     `icarus::trigger::SlidingWindowTrigger` module).
 * * `Thresholds` (list of names, mandatory): list of the discrimination
 *     thresholds to consider. A data product containing a digital signal is
 *     read for each one of the thresholds, and the tag of the data product is
 *     expected to be the instance name in this configuration parameter for the
 *     module label set in `TriggerGatesTag` (e.g. for a threshold of
 *     `"60"`, supposedly 60 ADC counts, and with `TriggerGatesTag` set to
 *     `"TrigSlidingWindows"`, the data product tag would be
 *     `TrigSlidingWindows:60`).
 * * `KeepThresholdName` (flag, optional): by default, output data products have
 *     each an instance name according to their threshold (from the `Threshold`
 *     parameter), unless there is only one threshold specified. If this
 *     parameter is specified as `true`, the output data product always
 *     includes the threshold instance name, even when there is only one
 *     threshold specified. If this parameter is specified as `false`, if there
 *     is only one threshold the default behaviour (of not adding an instance
 *     name) is confirmed; otherwise, it is a configuration error to have this
 *     parameter set to `false`.
 * * `Pattern` (configuration table, mandatory): describes the sliding window
 *     pattern; the configuration format for a pattern is described under
 *     `icarus::trigger::ns::fhicl::WindowPatternConfig`.
 * * `BeamGates` (input tag, _mandatory_): the data product with the beam gates
 * *   to run the simulation in;
 * * `BeamBits` (bitmask as 32-bit integral number): bits to be set in the
 *     produced `raw::Trigger` objects (see also `daq::TriggerDecoder` tool).
 * * `LogCategory` (string, default `TriggerSimulationOnGates`): name of
 *     category used to stream messages from this module into message facility.
 * 
 * An example job configuration is provided as
 * `simulate_sliding_window_trigger_icarus.fcl`.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `TriggerGatesTag` + `Thresholds`: input gate collections.
 * * `BeamGates` (`std::vector<sim::BeamGateInfo>`): the beam gate intervals
 *     to run the simulation on; one trigger result is produced and saved for
 *     each of the gates in this data product. The gates are interpreted
 *     following LArSoft convention for the simulation, with the times in
 *     nanoseconds and in
 *     @ref DetectorClocksSimulationTime "simulation time reference".
 * 
 * 
 * Output data products
 * =====================
 * 
 * * `std::vector<raw::Trigger>` (one instance per ADC threshold):
 *   list of triggers fired according to the configured trigger definition;
 *   there is one collection (and data product) per ADC threshold, and the
 *   data product has the same instance name as the input data one, unless
 *   there is only one threshold (see `TriggerGatesTag`, `Thresholds` and
 *   `KeepThresholdName` configuration parameters);
 *   one trigger object is produced for each of the beam gates found in the
 *   input data product specified by the `BeamGates` parameter.
 *   Each trigger object has the time stamp matching the first time the trigger
 *   criteria are satisfied. All triggers feature the bits specified in
 *   `BeamBits` configuration parameter.
 * 
 * 
 * 
 * Trigger logic algorithm
 * ========================
 * 
 * @anchor TriggerSimulationOnGates_Algorithm
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::TriggerSimulationOnGates` and its assumptions.
 * Nevertheless, more up-to-date information can be found in
 * `SlidingWindowTrigger` module (for the combination of the LVDS signals into
 * window-wide gates) and in `icarus::trigger::SlidingWindowPatternAlg`,
 * which applies the configured pattern logic to the input.
 * 
 * The module receives as input a multi-level trigger gate for each of the
 * windows to be considered.
 * On the first input (i.e. the first event), that input is parsed to learn
 * the windows and their relative position from the input trigger gates.
 * This topology will be used to apply the configured patterns. On the following
 * events, their input is checked to confirm the compatibility of the
 * composition of its windows with the one from that first event (both aspects
 * are handled by an `icarus::trigger::WindowTopologyManager` object).
 * 
 * All multi-level gates are set in coincidence with the beam gate by
 * multiplying the multi-level and the beam gates. Because of this, trigger
 * gates are suppressed everywhere except than during the beam gate (see below).
 * The reference time for the beam gates is the time configured in
 * `DetectorClocks` service provider
 * (`detinfo::DetectorClocks::BeamGateTime()`).
 * 
 * The algorithm handles independently multiple trigger patterns.
 * On each input, each configured pattern is applied based on the window
 * topology. Each pattern describes a minimum level of the trigger
 * gate in the window, that usually means the number of LVDS signals in
 * coincidence at any given time ("majority"). A pattern may have requirements
 * on the neighbouring windows in addition to the main one. The pattern is
 * satisfied if all involved windows pass their specific requirements at the
 * same time (coincidence between windows).
 * Each pattern is applied in turn to each of the windows (which is the "main"
 * window). The neighborhood described in the pattern is applied with respect to
 * that main window. The trigger fires if one or more of the windows satisfy
 * the pattern, and the trigger time is the one of the earliest satisfied
 * pattern (more precisely, the earliest tick when the coincidence required
 * by that pattern is satisfied).* 
 * All windows in the detector are considered independently, but the supported
 * patterns may only include components in the same cryostat. Therefore,
 * triggers are effectively on a single cryostat.
 * An object of class `icarus::trigger::SlidingWindowPatternAlg` applies this
 * logic: see its documentation for the most up-to-date details.
 * 
 * Eventually, for each event there are as many different trigger responses as
 * how many different patterns are configured (`Patterns` configuration
 * parameter), _times_ how many ADC thresholds are provided in input,
 * configured in `Thresholds`.
 * 
 * 
 * Beam gates
 * -----------
 * 
 * A single instance of this module can perform the simulation on several beam
 * gates. The values of these beam gates are picked from the data product
 * specified in `BeamGates`, event by event. The specified beam gate times are
 * on beam gate time scale, i.e. their reference time `0` is the time of the
 * beam gate as known by `detinfo::DetectorClocks::BeamGateTime()`.
 * In case the same beam gate is desired for all events, such data product can
 * be produced by `icarus::trigger::WriteBeamGateInfo` module.
 * The trigger data product collection produced by this module has the same
 * number of entries as the beam gates in the data product, and they match
 * one-to-one.
 * 
 * 
 * 
 * Technical aspects of the module
 * --------------------------------
 * 
 * @anchor TriggerSimulationOnGates_Tech
 * 
 * This module does not build the trigger gates of the sliding windows, but
 * rather it takes them as input (see e.g. `SlidingWindowTrigger` module).
 * Window topology (size of the windows and their relations) is stored in
 * `icarus::trigger::WindowChannelMap` objects, and its construction is
 * delegated to `icarus::trigger::WindowTopologyAlg` (under the hood of the
 * `WindowTopologyManager` class) which learns it from the actual trigger gate
 * input rather than on explicit configuration. Pattern definitions and
 * configuration are defined in `icarus::trigger::WindowPattern` and
 * `icarus::trigger::ns::fhicl::WindowPatternConfig` respectively. Trigger
 * simulation is delegated to `icarus::trigger::SlidingWindowPatternAlg`.
 * 
 * 
 * @todo Plots need to be thought and implemented.
 * 
 */
class icarus::trigger::TriggerSimulationOnGates
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

    fhicl::Sequence<std::string> Thresholds {
      Name("Thresholds"),
      Comment("tags of the thresholds to consider")
      };

    fhicl::OptionalAtom<bool> KeepThresholdName {
      Name("KeepThresholdName"),
      Comment
        ("add threshold to output product tag even with only one threshold")
      };

    icarus::trigger::ns::fhicl::WindowPatternTable Pattern {
      Name("Pattern"),
      Comment("trigger requirements as a trigger window pattern")
      };
 
    fhicl::Atom<art::InputTag> BeamGates {
      Name("BeamGates"),
      Comment("data product with all beam gates to run simulation into")
      };
    
    fhicl::Atom<std::uint32_t> BeamBits {
      Name("BeamBits"),
      Comment("bits to be set in the trigger object as beam identified")
      };

    fhicl::Atom<nanoseconds> TriggerTimeResolution {
      Name("TriggerTimeResolution"),
      Comment("resolution of trigger in time"),
      25_ns
      };
    
    fhicl::Atom<double> EventTimeBinning {
      Name("EventTimeBinning"),
      Comment("binning for the trigger time plot [second]"),
      300 // 5 minutes
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "TriggerSimulationOnGates" // default
      };
    
  }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit TriggerSimulationOnGates(Parameters const& config);

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
  using InputTriggerGate_t
    = icarus::trigger::SlidingWindowPatternAlg::InputTriggerGate_t;
  
  /// List of trigger gates.
  using TriggerGates_t
    = icarus::trigger::SlidingWindowPatternAlg::TriggerGates_t;
  
  /// Data structure to communicate internally a trigger response.
  using WindowTriggerInfo_t
    = icarus::trigger::SlidingWindowPatternAlg::AllTriggerInfo_t;
  
  /// Content for future histograms, binned.
  using BinnedContent_t = icarus::ns::util::FixedBins<double>;
  
  /// All information needed to generate plots for a specific threshold.
  struct ThresholdPlotInfo_t {
    BinnedContent_t eventTimes;
    BinnedContent_t HWtrigTimeVsBeam;
    BinnedContent_t triggerTimesVsHWtrig;
    BinnedContent_t triggerTimesVsBeam;
  };
  
  /// Type of list of gates to simulate trigger into.
  using BeamGates_t = std::vector<sim::BeamGateInfo>;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Name of ADC thresholds to read, and the input tag connected to their data.
  std::map<std::string, art::InputTag> fADCthresholds;
  
  /// Configured sliding window requirement pattern.
  WindowPattern const fPattern;
  
  art::InputTag const fBeamGateTag; ///< Data product of beam gates to simulate.
  
  std::uint32_t fBeamBits; ///< Bits for the beam gate being simulated.
  
  nanoseconds fTriggerTimeResolution; ///< Trigger resolution in time.
  
  double fEventTimeBinning; ///< Trigger time plot binning [s]
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------

  /// ROOT directory where all the plots are written.
  art::TFileDirectory fOutputDir;

  // --- END Service variables -------------------------------------------------

  
  // --- BEGIN Internal variables ----------------------------------------------
  
  /// Output data product instance names (same order as `fADCthresholds`).
  std::vector<std::string> fOutputInstances;
  
  /// Mapping of each sliding window with location and topological information.
  // mutable = not thread-safe
  mutable icarus::trigger::WindowTopologyManager fWindowMapMan;
  
  /// Pattern algorithm.
  std::optional<icarus::trigger::SlidingWindowPatternAlg> fPatternAlg;
  
  /// All plots in one practical sandbox.
  icarus::trigger::PlotSandbox fPlots;
  
  /// Proto-histogram information in a convenient packet; event-wide.
  ThresholdPlotInfo_t fEventPlotInfo;
  
  /// Proto-histogram information in a not-so-practical array; per threshold.
  std::vector<ThresholdPlotInfo_t> fThresholdPlots;

  ///< Count of fired triggers, per threshold.
  std::vector<std::atomic<unsigned int>> fTriggerCount;
  std::atomic<unsigned int> fTotalGates { 0U }; ///< Count of opened gates.
  
  
  // --- END Internal variables ------------------------------------------------
  

  // --- BEGIN --- Plot infrastructure -----------------------------------------
  
  /// @brief Initializes the full set of plots (all ADC thresholds).
  void initializePlots();
  
  /// Creates summary plots from proto-histogram data.
  void finalizePlots();
  
  /// Creates in `plots` sandbox all plots for threshold `threshold` from data
  /// in `plotInfo`.
  void makeThresholdPlots(
    std::string const& threshold,
    icarus::trigger::PlotSandbox& plots,
    ThresholdPlotInfo_t const& plotInfo
    );
  
  /// Creates in the main sandbox all event-wide plots.
  void makeEventPlots();
  
  /// Fills event-wide plots.
  void plotEvent(
    art::Event const& event, detinfo::DetectorTimings const& detTimings,
    std::vector<icarus::trigger::ApplyBeamGateClass> const& gates
    );
  
  /// Fills the plots for threshold index `iThr` with trigger information.
  void plotTriggerResponse(
    std::size_t iThr, std::string const& thrTag,
    WindowTriggerInfo_t const& triggerInfo,
    detinfo::DetectorTimings const& detTimings
    );
  
  // --- END ----- Plot infrastructure -----------------------------------------
  
  /**
   * @brief Performs the simulation for the specified ADC threshold.
   * @param event _art_ event to read data from and put results into
   * @param detTimings detector clocks service provider proxy
   * @param beamGates list of all beam gates to evaluate
   * @param iThr index of the threshold in the configuration
   * @param thr value of the threshold (ADC counts)
   * @param firstTriggerNumber the next unassigned trigger number
   * @return the trigger response information
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
  std::vector<WindowTriggerInfo_t> produceForThreshold(
    art::Event& event,
    detinfo::DetectorTimings const& detTimings,
    std::vector<ApplyBeamGateClass> const& beamGates,
    std::size_t const iThr, std::string const& thrTag,
    unsigned int firstTriggerNumber
    );
  
  /**
   * @brief Converts the trigger information into a `raw::Trigger` object.
   * @param triggerNumber the unique number to assign to this trigger
   * @param info the information about the fired trigger
   * @return a `raw::Trigger` object with all the information encoded
   * 
   * The trigger described by `info` is encoded into a `raw::Trigger` object.
   * The trigger _must_ have fired.
   */
  raw::Trigger triggerInfoToTriggerData(
    detinfo::DetectorTimings const& detTimings,
    ApplyBeamGateClass const& beamGate,
    unsigned int triggerNumber, WindowTriggerInfo_t const& info
    ) const;
  
  
  /// Prints the summary of fired triggers on screen.
  void printSummary() const;
  
  /// Creates and returns a 1D histogram filled with `binnedContent`.
  TH1* makeHistogramFromBinnedContent(
    icarus::trigger::PlotSandbox& plots,
    std::string const& name, std::string const& title,
    BinnedContent_t const& binnedContent
    ) const;

  /// Shortcut to create an `ApplyBeamGate` with the specified `gate`.
  icarus::trigger::ApplyBeamGateClass makeMyBeamGate
    (detinfo::DetectorTimings const& detTimings, sim::BeamGateInfo const& gate) const
    {
      // the input gate is assumed to be relative to the global beam gate
      // opening (which is the implicit convention of simulation time and of
      // sim::BeamGateInfo in my understanding - [petrillo@slac.stanford.edu])
      // so it does not need further processing here
      return makeApplyBeamGate(
        nanoseconds{ gate.Width() }, nanoseconds{ gate.Start() },
        detTimings.clockData(), fLogCategory
        );
    }
  
  //@{ 
  /// Shortcut to create `ApplyBeamGate` from a list of gates.
  std::vector<icarus::trigger::ApplyBeamGateClass> makeMyBeamGates
    (detinfo::DetectorTimings const& detTimings, BeamGates_t const& gates) const
    {
      std::vector<icarus::trigger::ApplyBeamGateClass> applyGates;
      for (sim::BeamGateInfo const& gate: gates)
        applyGates.push_back(makeMyBeamGate(detTimings, gate));
      return applyGates;
    }
  std::vector<icarus::trigger::ApplyBeamGateClass> makeMyBeamGates
    (art::Event const* event, BeamGates_t const& gates) const
    { return makeMyBeamGates(icarus::ns::util::makeDetTimings(event), gates); }
  std::vector<icarus::trigger::ApplyBeamGateClass> makeMyBeamGates
    (art::Event const& event, BeamGates_t const& gates) const
    { return makeMyBeamGates(&event, gates); }
  //@}
  
  
  //@{
  /// Returns the time of the event in seconds from The Epoch.
  static double eventTimestampInSeconds(art::Timestamp const& time);
  static double eventTimestampInSeconds(art::Event const& event);
  //@}

}; // icarus::trigger::TriggerSimulationOnGates



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
icarus::trigger::TriggerSimulationOnGates::TriggerSimulationOnGates
  (Parameters const& config)
  : art::EDProducer       (config)
  // configuration
  , fPattern              (config().Pattern())
  , fBeamGateTag          (config().BeamGates())
  , fBeamBits             (config().BeamBits())
  , fTriggerTimeResolution(config().TriggerTimeResolution())
  , fEventTimeBinning     (config().EventTimeBinning())
  , fLogCategory          (config().LogCategory())
  // services
  , fOutputDir (*art::ServiceHandle<art::TFileService>())
  // internal and cached
  , fWindowMapMan
    { *lar::providerFrom<geo::Geometry>(), fLogCategory + "_WindowMapManager" }
  , fPlots(
     fOutputDir, "", "requirement: " + fPattern.description()
    )
  , fEventPlotInfo{
        BinnedContent_t{ fEventTimeBinning }  // eventTimes
      , BinnedContent_t{                      // HWtrigTimeVsBeam
          fTriggerTimeResolution.convertInto
            <detinfo::timescales::trigger_time::interval_t>().value()
        }
      , BinnedContent_t{                      // triggerTimesVsHWtrig
          fTriggerTimeResolution.convertInto
            <detinfo::timescales::trigger_time::interval_t>().value()
        }
      , BinnedContent_t{                      // triggerTimesVsBeam
          fTriggerTimeResolution.convertInto
            <detinfo::timescales::trigger_time::interval_t>().value()
        }
    }
{
  
  //
  // more complex parameter parsing
  //
  std::string const& discrModuleLabel = config().TriggerGatesTag();
  for (std::string const& threshold: config().Thresholds())
    fADCthresholds[threshold] = art::InputTag{ discrModuleLabel, threshold };
  
  // initialization of a vector of atomic is not as trivial as it sounds...
  fTriggerCount = std::vector<std::atomic<unsigned int>>(fADCthresholds.size());
  std::fill(fTriggerCount.begin(), fTriggerCount.end(), 0U);

  //
  // input data declaration
  //
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // trigger primitives
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    icarus::trigger::TriggerGateReader<>{ inputDataTag }
      .declareConsumes(consumesCollector());
  } // for
  
  //
  // output data declaration
  //
  // keepThresholdName is true if we write instance name in output data products
  bool const keepThresholdName
    = config().KeepThresholdName().value_or(config().Thresholds().size() > 1);
  if (!keepThresholdName && (config().Thresholds().size() > 1)) {
    throw art::Exception(art::errors::Configuration)
      << config().KeepThresholdName.name()
      << " can be set to `true` only when a single threshold is specified ("
      << config().Thresholds.name() << " has " << config().Thresholds().size()
      << ")";
  }
  
  for (auto const& inputDataTag: util::const_values(fADCthresholds)) {
    std::string const outputInstance
      = keepThresholdName? inputDataTag.instance(): "";
    produces<std::vector<raw::Trigger>>(outputInstance);
    fOutputInstances.push_back(outputInstance);
  }
  
  {
    mf::LogInfo log(fLogCategory);
    log << "\nConfigured " << fADCthresholds.size() << " thresholds (ADC):";
    for (auto const& [ thresholdTag, dataTag ]: fADCthresholds)
      log << "\n * " << thresholdTag << " (from '" << dataTag.encode() << "')";
    
  } // local block
  
  
} // icarus::trigger::TriggerSimulationOnGates::TriggerSimulationOnGates()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::beginJob() {
  
  initializePlots();
  
} // icarus::trigger::TriggerSimulationOnGates::beginJob()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::produce(art::Event& event)
{
  
  //
  // prepare all the gates to run the simulation on
  //
  detinfo::DetectorTimings const detTimings {
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  
  std::vector<icarus::trigger::ApplyBeamGateClass> const beamGates {
    makeMyBeamGates(
      detTimings,
      event.getProduct<std::vector<sim::BeamGateInfo>>(fBeamGateTag)
      )
    };

  
  { // BEGIN local block
    mf::LogDebug log { fLogCategory };
    log << "Trigger simulation for " << beamGates.size() << " gates";
    if (!beamGates.empty()) {
      log << " ('" << fBeamGateTag.encode() << "'):";
      for (auto const& [iGate, gate]: util::enumerate(beamGates)) {
        log << "\n [" << iGate << "]  " << gate;
      } // for
    } // if
  } // END local block
  
  
  //
  // run the simulation on each threshold in turn
  //
  mf::LogDebug log(fLogCategory); // this will print at the end of produce()
  log << "Event " << event.id() << ":";
  
  // FIXME these two operations should be atomic
  unsigned int const firstTriggerNumber
    = fTotalGates.fetch_add(beamGates.size());
  
  for (auto const& [ iThr, thrTag ]
    : util::enumerate(util::get_elements<0U>(fADCthresholds))
  ) {
    
    std::vector<WindowTriggerInfo_t> const triggers = produceForThreshold
      (event, detTimings, beamGates, iThr, thrTag, firstTriggerNumber);
    
    log << "\n * threshold " << thrTag << ": ";
    icarus::ns::util::PassCounter gateResults;
    for (WindowTriggerInfo_t const& triggerInfo: triggers)
      gateResults.add(triggerInfo.info.fired());
    log << gateResults.passed() << "/" << gateResults.total()
      << " gates triggered";
    
  } // for
  
  //
  // event-level plots
  //
  plotEvent(event, detTimings, beamGates);
  
} // icarus::trigger::TriggerSimulationOnGates::produce()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::endJob() {
  
  finalizePlots();
  
  printSummary();
  
} // icarus::trigger::TriggerSimulationOnGates::endJob()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::initializePlots() {
  
#if 0
  //
  // overview plots with different settings
  //
  
  std::vector<std::string> thresholdLabels;
  thresholdLabels.reserve(size(fADCthresholds));
  for (std::string thr: util::get_elements<0U>(fADCthresholds))
    thresholdLabels.push_back(std::move(thr));
  
  auto const beamGate = makeMyBeamGate();
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
    "Triggering pass fraction"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";trigger pass fraction",
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
  
  
  // we allow some fixed margin in the plot, just in case:
  constexpr microseconds beamPlotPadding { 4_us };
  
  using detinfo::timescales::trigger_time;
  
  // hardware trigger may happen at any place within the beam gate,
  // and in this plot range I want to include the full beam gate;
  // so I take a beam gate worth of time before the trigger time,
  // and as much after it; since this plot is relative to the hardware trigger,
  // the hardware trigger time itself is 0
  icarus::ns::util::BinningSpecs const HWtrigBinning = alignBinningTo(
    icarus::ns::util::BinningSpecs{
      (- beamGate.length() - beamPlotPadding).value(),
      (beamGate.length() + beamPlotPadding).value(),
      fTriggerTimeResolution.convertInto<trigger_time::interval_t>().value()
      },
      0.0
    );
  fPlots.make<TH2F>(
    "TriggerTimeVsHWTrig",
    "Time of the trigger"
      ";trigger time (relative to hardware trigger)  [ #mus ]"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";triggers",
    HWtrigBinning.nBins(), HWtrigBinning.lower(), HWtrigBinning.upper(),
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );

  icarus::ns::util::BinningSpecs const beamGateBinning = alignBinningTo(
    icarus::ns::util::BinningSpecs{
      (-beamPlotPadding).value(),
      (beamGate.length() + beamPlotPadding).value(),
      fTriggerTimeResolution.convertInto<trigger_time::interval_t>().value()
      },
      0.0
    );
  fPlots.make<TH2F>(
    "TriggerTimeVsBeamGate",
    "Time of the trigger"
      ";trigger time (relative to beam gate opening)  [ #mus ]"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";triggers",
    beamGateBinning.nBins(), beamGateBinning.lower(), beamGateBinning.upper(),
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );

  // 
  // per-threshold plots; should this initialization be set into its own method?
  // 
  for (auto const& [ thr, info ]
    : util::zip(util::get_elements<0U>(fADCthresholds), fThresholdPlots))
  {
    icarus::trigger::PlotSandbox& plots
      = fPlots.addSubSandbox("Thr" + thr, "Threshold: " + thr);
    
    plots.make<TGraph>(
      "TriggerTimeVsHWTrigVsBeam",
      "Time of the trigger: emulated vs. hardware"
        ";hardware trigger time (relative to beam gate opening)  [ #mus ]"
        ";emulated trigger time (relative to beam gate opening)  [ #mus ]"
      );
    
  } // for thresholds
  
  fThresholdPlots.resize(
    size(fADCthresholds),
    {
      BinnedContent_t{ fEventTimeBinning },         // eventTimes
      BinnedContent_t{ HWtrigBinning.binWidth() },  // HWtrigTimeVsBeam
      BinnedContent_t{ HWtrigBinning.binWidth() },  // triggerTimesVsHWtrig
      BinnedContent_t{ beamGateBinning.binWidth() } // triggerTimesVsBeam
    }
    );
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::initializePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::finalizePlots() {
  
#if 0
  
  for (auto const& [ thr, info ]
    : util::zip(util::get_elements<0U>(fADCthresholds), fThresholdPlots))
  {
    icarus::trigger::PlotSandbox& plots = fPlots.demandSandbox("Thr" + thr);
    makeThresholdPlots(thr, plots, info);
    if (plots.empty()) fPlots.deleteSubSandbox(plots.name());
  } // for thresholds
  
  makeEventPlots();
  
#endif // 0

} // icarus::trigger::TriggerSimulationOnGates::finalizePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::makeThresholdPlots(
  std::string const& threshold,
  icarus::trigger::PlotSandbox& plots,
  ThresholdPlotInfo_t const& plotInfo
) {
  
#if 0
  
  BinnedContent_t const* content;
  
  content = &(plotInfo.eventTimes);
  makeHistogramFromBinnedContent(plots,
    "TriggerTime",
    "Time of the triggered events"
      ";time"
      ";triggered events  [ / " + std::to_string(content->binWidth())
        + "\" ]",
    *content
    );
  
  content = &(plotInfo.HWtrigTimeVsBeam);
  makeHistogramFromBinnedContent(plots,
    "HWTrigVsBeam",
    "Time of the hardware trigger"
      ";trigger time (relative to beam gate)  [ #mus ]"
      ";events  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
  content = &(plotInfo.triggerTimesVsHWtrig);
  makeHistogramFromBinnedContent(plots,
    "TriggerTimeVsHWTrig",
    "Time of the trigger"
      ";trigger time (relative to hardware trigger)  [ #mus ]"
      ";triggers  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
  content = &(plotInfo.triggerTimesVsBeam);
  makeHistogramFromBinnedContent(plots,
    "TriggerTimeVsBeamGate",
    "Time of the trigger"
      ";trigger time (relative to beam gate opening)  [ #mus ]"
      ";triggers  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::makeThresholdPlots()



//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::makeEventPlots() {
  
#if 0
  
  BinnedContent_t const* content;
  
  content = &(fEventPlotInfo.eventTimes);
  makeHistogramFromBinnedContent(fPlots,
    "EventTime",
    "Time of the events"
      ";time"
      ";events  [ / " + std::to_string(content->binWidth())
        + "\" ]",
    *content
    );
  
  content = &(fEventPlotInfo.HWtrigTimeVsBeam);
  makeHistogramFromBinnedContent(fPlots,
    "HWTrigVsBeam",
    "Time of the hardware trigger"
      ";trigger time (relative to beam gate)  [ #mus ]"
      ";events  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::makeEventPlots()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::produceForThreshold(
  art::Event& event,
  detinfo::DetectorTimings const& detTimings,
  std::vector<ApplyBeamGateClass> const& beamGates,
  std::size_t const iThr, std::string const& thrTag,
  unsigned int firstTriggerNumber
) -> std::vector<WindowTriggerInfo_t> {
  
//   auto& plotInfo = fThresholdPlots[iThr];
  
  //
  // get the input
  //
  art::InputTag const& dataTag = fADCthresholds.at(thrTag);
  auto const& gates = icarus::trigger::ReadTriggerGates(event, dataTag);
  
  
  // extract or verify the topology of the trigger windows
  if (fWindowMapMan(gates))
    fPatternAlg.emplace(*fWindowMapMan, fPattern, fLogCategory);
  assert(fPatternAlg);
  
  //
  // simulate the trigger response on all beam gates
  //
  std::vector<WindowTriggerInfo_t> allTriggerInfo; // one per gate
  auto triggers = std::make_unique<std::vector<raw::Trigger>>();
  unsigned int triggerNumber = firstTriggerNumber;
  for (auto const& beamGate: beamGates) {
    WindowTriggerInfo_t const triggerInfo
      = fPatternAlg->simulateResponse(beamGate.applyToAll(gates));
    
    // FIXME what do we do with statistics and plots?
    if (triggerInfo) {
      ++fTriggerCount[iThr]; // keep the unique count
//       plotInfo.eventTimes.add(eventTimestampInSeconds(event));
    }
    
    //
    // fill the plots
    //
    plotTriggerResponse(iThr, thrTag, triggerInfo, detTimings);

    //
    // create and store the data product
    //
    triggers->push_back(
      triggerInfoToTriggerData
        (detTimings, beamGate, triggerNumber++, triggerInfo)
      );
    allTriggerInfo.push_back(std::move(triggerInfo));
    
  } // for beam gates
  
  event.put(std::move(triggers), fOutputInstances[iThr]);
  
  return allTriggerInfo;
  
} // icarus::trigger::TriggerSimulationOnGates::produceForThreshold()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::plotEvent(
  art::Event const& event, detinfo::DetectorTimings const& detTimings,
  std::vector<icarus::trigger::ApplyBeamGateClass> const& gates
) {
  
#if 0
  
  detinfo::timescales::trigger_time const beamGateTime
    { detTimings.toTriggerTime(detTimings.BeamGateTime()) };
  
  fEventPlotInfo.eventTimes.add(eventTimestampInSeconds(event));
  fEventPlotInfo.HWtrigTimeVsBeam.add(-beamGateTime.value());
  
  // `gates` is currently unused; it may be used e.g. to show how many gates
  // were tested in each event
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::plotEvent()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::plotTriggerResponse(
  std::size_t iThr, std::string const& thrTag,
  WindowTriggerInfo_t const& triggerInfo,
  detinfo::DetectorTimings const& detTimings
) {
  
#if 0
  
  bool const fired = triggerInfo.info.fired();
  
  fPlots.demand<TEfficiency>("Eff").Fill(fired, iThr);
  
  if (fired) {
    using namespace detinfo::timescales;
    
    // time of the beam gate in hardware trigger time scale
    trigger_time const beamGateTime
      { detTimings.toTriggerTime(detTimings.BeamGateTime()) };
    
    optical_tick const thisTriggerTick { triggerInfo.info.atTick() };
    trigger_time const thisTriggerTimeVsHWtrig
      { detTimings.toTriggerTime(thisTriggerTick) };
    time_interval const thisTriggerTimeVsBeamGate
      { thisTriggerTimeVsHWtrig - beamGateTime };
    
    mf::LogTrace(fLogCategory)
      << "Trigger " << fPattern.tag() << " at tick " << thisTriggerTick
      << " (" << thisTriggerTimeVsHWtrig << " vs. HW trigger, "
      << thisTriggerTimeVsBeamGate << " vs. beam gate)"
      ;
    
    fPlots.demand<TH1>("NTriggers").Fill(iThr);
    fPlots.demand<TH2>("TriggerTick").Fill(thisTriggerTick.value(), iThr);
    fPlots.demand<TH2>("TriggerTimeVsHWTrig").Fill
      (thisTriggerTimeVsHWtrig.value(), iThr);
    fPlots.demand<TH2>("TriggerTimeVsBeamGate").Fill
      (thisTriggerTimeVsBeamGate.value(), iThr);
    
    icarus::trigger::PlotSandbox& plots{ fPlots.demandSandbox("Thr" + thrTag) };
//     plots.demand<TGraph>("TriggerTimeVsHWTrigVsBeam").AddPoint( // ROOT 6.24?
    TGraph& graph = plots.demand<TGraph>("TriggerTimeVsHWTrigVsBeam");
    graph.SetPoint(graph.GetN(),
      -beamGateTime.value(), thisTriggerTimeVsBeamGate.value()
      );
    
    ThresholdPlotInfo_t& plotInfo { fThresholdPlots[iThr] };
    plotInfo.HWtrigTimeVsBeam.add(-beamGateTime.value());
    plotInfo.triggerTimesVsHWtrig.add(thisTriggerTimeVsHWtrig.value());
    plotInfo.triggerTimesVsBeam.add(thisTriggerTimeVsBeamGate.value());
    
  }
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::plotTriggerResponse()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::printSummary() const {
  
  //
  // summary from our internal counters
  //
  mf::LogInfo log(fLogCategory);
  log
    << "Summary of triggers for " << fTriggerCount.size()
    << " thresholds (ADC) with pattern: " << fPattern.description()
    ;
  for (auto const& [ count, thr ]
    : util::zip(fTriggerCount, util::get_elements<0U>(fADCthresholds)))
  {
    log << "\n  threshold " << thr
      << ": " << count;
    if (fTotalGates > 0U) {
      log << "/" << fTotalGates
        << " (" << (double(count) / fTotalGates * 100.0) << "%)";
    }
    else log << " gates triggered";
  } // for
  
} // icarus::trigger::TriggerSimulationOnGates::printSummary()


//------------------------------------------------------------------------------
raw::Trigger
icarus::trigger::TriggerSimulationOnGates::triggerInfoToTriggerData(
  detinfo::DetectorTimings const& detTimings,
  ApplyBeamGateClass const& beamGate,
  unsigned int triggerNumber, WindowTriggerInfo_t const& info
) const {
  
  return {
    triggerNumber,                      // counter
    info.info.fired()                   // trigger time
      ? double(detTimings.toElectronicsTime(info.info.atTick()))
      : std::numeric_limits<double>::lowest()
      ,
    double(detTimings.toElectronicsTime(beamGate.tickRange().start())),
                                        // beam gate in electronics time scale
    (info.info.fired()? fBeamBits: 0)   // bits
    };
  
} // icarus::trigger::TriggerSimulationOnGates::triggerInfoToTriggerData()


//------------------------------------------------------------------------------
TH1*
icarus::trigger::TriggerSimulationOnGates::makeHistogramFromBinnedContent(
  icarus::trigger::PlotSandbox& plots,
  std::string const& name, std::string const& title,
  BinnedContent_t const& binnedContent
) const {
  
  if (binnedContent.empty()) return nullptr;
  
  TH1* hist = plots.make<TH1F>(
    name, title,
    binnedContent.nBins(), binnedContent.min(), binnedContent.max()
    );
  
  // directly transfer the content bin by bin
  unsigned int total = 0U;
  for (auto [ iBin, count ]: util::enumerate(binnedContent)) {
    hist->SetBinContent(iBin + 1, count);
    total += count;
  }
  hist->SetEntries(static_cast<double>(total));
  return hist;
} // icarus::trigger::TriggerSimulationOnGates::makeHistogramFromBinnedContent


//------------------------------------------------------------------------------
double icarus::trigger::TriggerSimulationOnGates::eventTimestampInSeconds
  (art::Timestamp const& time)
{
  // high value: seconds from the Epoch (Jan 1, 1970 UTC?);
  // low value: nanoseconds after that the start of that second
  return static_cast<double>(time.timeHigh())
    + static_cast<double>(time.timeHigh()) * 1e-9;
} // icarus::trigger::TriggerSimulationOnGates::eventTimestampInSeconds()


//------------------------------------------------------------------------------
double icarus::trigger::TriggerSimulationOnGates::eventTimestampInSeconds
  (art::Event const& event)
  { return eventTimestampInSeconds(event.time()); }


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::TriggerSimulationOnGates)


//------------------------------------------------------------------------------
