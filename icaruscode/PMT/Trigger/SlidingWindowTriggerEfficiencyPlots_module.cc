/**
 * @file   SlidingWindowTriggerEfficiencyPlots_module.cc
 * @brief  Plots of efficiency for triggers based on PMT channel sliding window.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 9, 2020
 * @see    icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h"
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.h"
#include "icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.h" // WindowTopologyManager
#include "icaruscode/PMT/Trigger/Algorithms/WindowPatternConfig.h"
#include "icaruscode/PMT/Trigger/Algorithms/WindowPattern.h"
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h"
#include "icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h" // gatesIn()
#include "icarusalg/Utilities/ROOTutils.h" // util::ROOT
#include "icarusalg/Utilities/sortBy.h" // also icarus::util::sortCollBy()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/TensorIndices.h" // util::MatrixIndices
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_time_ticks..
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

// ROOT libraries
#include "TTree.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <algorithm> // std::sort(), std::transform(), ...
#include <vector>
#include <array>
#include <memory> // std::unique_ptr
#include <utility> // std::pair<>, std::move()
#include <limits> // std::numeric_limits<>
#include <type_traits> // std::is_pointer_v, ...
#include <cstddef> // std::size_t
#include <cassert>


//------------------------------------------------------------------------------
namespace {
  
  /// Moves all elements of `src` at the end of `dest`.
  template <typename DestColl, typename SrcColl>
  DestColl& appendCollection(DestColl& dest, SrcColl&& src);
  
} // local namespace


//------------------------------------------------------------------------------

// --- BEGIN -- ROOT tree helpers ----------------------------------------------
/**
 * @brief Class managing the serialization of trigger responses in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignResponse()`, the proper branch address is assigned the specified
 * trigger response (`true` or `false`).
 *
 * The branch structure is: a `RespTxxSxx/O` branch for each threshold and
 * settings label, with a single branch per element.
 *
 */
struct ResponseTree: public icarus::trigger::details::TreeHolder {

  // `std::vector<bool>` is too special for us. Let's pack it to go.
  util::MatrixIndices indices;
  std::unique_ptr<bool[]> RespTxxSxx;

  
  /**
   * @brief Constructor: accommodates that many thresholds and trigger settings.
   * @tparam Thresholds an iterable type yielding objects convertible to numbers
   * @tparam Settings an iterable type yielding objects convertible to string
   * @param tree the ROOT tree to add branches to (managed elsewhere)
   * @param thresholds collection of thresholds to be included
   * @param settings collection of trigger settings to be included
   * 
   * The `thresholds` must be convertible to ADC counts (i.e. numbers), while
   * the `settings` elements must support conversion to string via `to_string()`
   * function call (more precisely, `util::to_string()` from
   * `larcorealg/CoreUtils/StdUtils.h`).
   */
  template <typename Thresholds, typename Settings>
  ResponseTree
    (TTree& tree, Thresholds const& thresholds, Settings const& settings);

  /// Assigns the response for the specified trigger.
  void assignResponse(std::size_t iThr, std::size_t iSettings, bool resp);

}; // struct ResponseTree


// --- END -- ROOT tree helpers ------------------------------------------------


//------------------------------------------------------------------------------
namespace icarus::trigger { class SlidingWindowTriggerEfficiencyPlots; }
/**
 * @brief Produces plots about trigger simulation and trigger efficiency.
 * 
 * This module is an implementation of `TriggerEfficiencyPlotsBase`
 * for a trigger defined as a pattern of sliding windows.
 * 
 * Note that the multi-level logical waveforms from the sliding windows are
 * expected to be provided as input.
 * 
 * The single sliding window with the highest activity is picked as a reference.
 * A requirement on the number of trigger primitives "on" in that window is
 * imposed. Additional requirements may be imposed on three other sliding
 * windows: the upstream one (if any), the downstream one (if any) and the
 * opposite one in the same cryostat.
 * 
 * As usual for `TriggerEfficiencyPlotsBase` based modules, this happens for
 * every configured PMT discrimination threshold.
 * 
 * 
 * Trigger logic algorithm
 * ========================
 * 
 * @anchor SlidingWindowTriggerEfficiencyPlots_Algorithm
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::SlidingWindowTriggerEfficiencyPlots` and its assumptions.
 * 
 * The module receives as input a multi-level trigger gate for each of the
 * windows to be considered.
 * On the first input (i.e. the first event), that input is parsed to learn
 * the windows and their relative position from the input trigger gates
 * (`initializeTopologicalMaps()`). This topology will be used to apply the
 * configured patterns. On the following events, their input is checked to
 * confirm the compatibility of the composition of its windows with the one from
 * that first event (`verifyTopologicalMap()`).
 * 
 * All multi-level gates are set in coincidence with the beam gate by
 * multiplying the multi-level and the beam gates. Beacuse of this, trigger
 * gates are suppressed everywhere except than during the beam gate.
 * The beam gate opens at a time configured in `DetectorClocks` service provider
 * (`detinfo::DetectorClocks::BeamGateTime()`) and has a duration configured
 * in this module (`BeamGateDuration`).
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
 * that main window. The method `applyWindowPattern()` performs this
 * combination.
 * The trigger fires if one or more of the windows satisfy the pattern, and the
 * trigger time is the one of the earliest satisfied pattern (more precisely,
 * the earliest tick when the coincidence required by that pattern is
 * satisfied).
 * 
 * All windows in the detector are considered independently, but the supported
 * patterns may only include components in the same cryostat. Therefore,
 * triggers are effectively on a single cryostat.
 * 
 * Eventually, for each event there are as many different trigger responses as
 * how many different patterns are configured (`Patterns` configuration
 * parameter), _times_ how many ADC thresholds are provided in input,
 * configured in `Thresholds`.
 * 
 * While there _is_ a parameter describing the time resolution of the trigger
 * (`TriggerTimeResolution`), this is currently only used for aesthetic purposes
 * to choose the binning of some plots: the resolution is _not_ superimposed
 * to the gates (yet).
 * 
 * The set of plots and their organization are described in the documentation of
 * `icarus::trigger::TriggerEfficiencyPlotsBase`.
 * In the following documentation only the additions are described.
 * 
 * 
 * Output plots
 * -------------
 * 
 * A generic "setting" of `icarus::trigger::TriggerEfficiencyPlotsBase` is
 * in this module represented by a tag encoding the characteristics of the
 * pattern (see `WindowPattern::tag()`). The folders and plots will
 * identify each requirement with tags like `M8` or `M5O2`.
 * 
 * There are different "types" of plots. Some
 * @ref SlidingWindowTriggerEfficiencyPlots_SelectionPlots "do not depend on triggering at all",
 * like the deposited energy distribution. Others
 * @ref SlidingWindowTriggerEfficiencyPlots_MultiTriggerPlots "cross different trigger definitions",
 * like the trigger efficiency as function of trigger requirement. Others still
 * @ref SlidingWindowTriggerEfficiencyPlots_SingleTriggerPlots "assume a single trigger definition":
 * this is the case of trigger efficiency plots versus energy. Finally, there are
 * @ref SlidingWindowTriggerEfficiencyPlots_SingleTriggerResponsePlots "plots that depend on a specific trigger definition and outcome":
 * this is the case of all the plots including only triggering or non-triggering
 * events.
 * 
 * There are a few plots that are produced by this module in addition to the
 * ones in `TriggerEfficiencyPlotsBase`. They are described below.
 * 
 * All the plots are always relative to a specific optical detector channel
 * threshold (ADC) and a broad event category.
 * 
 * 
 * ### Plots independent of the triggers (selection plots)
 * 
 * @anchor SlidingWindowTriggerEfficiencyPlots_SelectionPlots
 * 
 * Only @ref TriggerEfficiencyPlotsBase_SelectionPlots "the standard plots"
 * from `TriggerEfficiencyPlotsBase` are produced in this category.
 * 
 * 
 * ### Plots including different trigger requirements
 * 
 * @anchor SlidingWindowTriggerEfficiencyPlots_MultiTriggerPlots
 * 
 * In addition to @ref TriggerEfficiencyPlotsBase_MultiTriggerPlots "the plots"
 * from `TriggerEfficiencyPlotsBase`, the following plots are also produced:
 * 
 * * `Eff`: trigger efficiency defined as number of triggered events over the
 *   total number of events, as function of the pattern (as encoded above);
 *   uncertainties are managed by `TEfficiency`.
 * * `Triggers`: trigger count as function of the pattern (as encoded above).
 * * `TriggerTick`: distribution of the time of the earliest trigger for the
 *   event, as function of the pattern (as in `Eff`). Each event appears at most
 *   once for each trigger pattern.
 * 
 * 
 * ### Plots depending on a specific trigger definition
 * 
 * @anchor SlidingWindowTriggerEfficiencyPlots_SingleTriggerPlots
 * 
 * Only @ref TriggerEfficiencyPlotsBase_SingleTriggerPlots "the standard plots"
 * from `TriggerEfficiencyPlotsBase` are produced in this category.
 * 
 *
 * ### Plots depending on a specific trigger definition and response
 *
 * @anchor SlidingWindowTriggerEfficiencyPlots_SingleTriggerResponsePlots
 * 
 * Only @ref TriggerEfficiencyPlotsBase_SingleTriggerResponsePlots "the standard plots"
 * from `TriggerEfficiencyPlotsBase` are produced in this category.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * @anchor SlidingWindowTriggerEfficiencyPlots_Configuration
 * 
 * In addition to
 * @ref TriggerEfficiencyPlotsBase_Configuration "all the configuration parameters"
 * from `TriggerEfficiencyPlotsBase`, the following one is also present:
 * 
 * * `Patterns` (list of pattern specifications, _mandatory_): a list of
 *     alternative pattern requirements for the definition of a trigger; each
 *     value is a table on its own, with the following elements:
 *     * `inMainWindow` (integer, mandatory): the minimum number of primitives
 *         to be fired in the central ("main") window of the pattern
 *     * `inDownstreamWindow` (integer, default: `0`): the minimum number of
 *         primitives to be fired in the window downstream of the main one
 *         (downstream is farther from the face beam neutrinos enter the
 *         detector through, i.e. larger _z_)
 *     * `inUpstreamWindow` (integer, default: `0`): the minimum number of
 *         primitives to be fired in the window upstream of the main one
 *         (upstream is closer to the face beam neutrinos enter the
 *         detector through, i.e. smaller _z_)
 *     * `inOppositeWindow` (integer, default: `0`): the minimum number of
 *         primitives to be fired in the window opposite to the main one
 *     * `requireDownstreamWindow` (flag, default: `false`): if set, this
 *         pattern is applied only on main windows that have a window downstream
 *         of them, i.e. farther from the face of the detector the beam enters
 *         through (larger _z_ coordinate); if set to `false`, if a main window
 *         has no downstream window, the downstream window requirement is always
 *         considered passed
 *     * `requireUpstreamWindow` (flag, default: `false`): if set, this pattern
 *         is applied only on main windows that have a window upstream of them,
 *         in the same way as for the `requireDownstreamWindow` setting
 *         described above
 * 
 * An example job configuration is provided as
 * `makeslidingwindowtriggerplots_icarus.fcl`.
 * 
 * 
 * Technical description of the module
 * ====================================
 * 
 * This module class is derived from
 * `icarus::trigger::TriggerEfficiencyPlotsBase`, which provides a backbone to
 * perform the simulation of triggers and plotting of their efficiency.
 * 
 * There is no superior design involved in this separation, but rather the goal
 * to share most code possible between different modules which simulate
 * different trigger patterns and as a consequence might have specific plots to
 * fill.
 * 
 * This module redefines:
 * 
 * * `initializePlotSet()` to define the
 *   @ref SlidingWindowTriggerEfficiencyPlots_MultiTriggerPlots "few additional plots"
 *   needed;
 * * `simulateAndPlot()`, which must always be defined, and which connects
 *   the simulation pieces and the plotting.
 * 
 * It does not redefine `initializeEfficiencyPerTriggerPlots()` nor
 * `initializeEventPlots()` because there are no additional plots of the types
 * these functions deal with.
 * 
 * The event categories are from the default list (`DefaultPlotCategories`) too.
 * 
 */
class icarus::trigger::SlidingWindowTriggerEfficiencyPlots
  : public art::EDAnalyzer
  , private icarus::trigger::TriggerEfficiencyPlotsBase
{
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config: public TriggerEfficiencyPlotsBase::Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    icarus::trigger::ns::fhicl::WindowPatternSequence Patterns {
      Name("Patterns"),
      Comment("sliding window pattern requirements")
      };
    
  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit SlidingWindowTriggerEfficiencyPlots(Parameters const& config);

  // --- END Constructors ------------------------------------------------------


  // --- BEGIN Framework hooks -------------------------------------------------

  /// Initializes the plots.
  virtual void beginJob() override;
  
  /// Fills the plots. Also extracts the information to fill them with.
  virtual void analyze(art::Event const& event) override;
  
  /// Prints end-of-job summaries.
  virtual void endJob() override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  /// List of configured patterns.
  using WindowPatterns_t = icarus::trigger::WindowPatterns_t;
  
  using TriggerInfo_t = details::TriggerInfo_t; // type alias

  /// Data structure to communicate internally a trigger response.
  using WindowTriggerInfo_t
    = icarus::trigger::SlidingWindowPatternAlg::AllTriggerInfo_t;

  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Configured sliding window requirement patterns.
  WindowPatterns_t const fPatterns;
  
  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Internal variables ----------------------------------------------
  
  /// Mapping of each sliding window with location and topological information.
  // mutable = not thread-safe; optional to allow delayed construction
  mutable icarus::trigger::WindowTopologyManager fWindowMapMan;
  
  /// All algorithm instances, one per pattern.
  std::vector<icarus::trigger::SlidingWindowPatternAlg> fPatternAlgs;
  
  std::unique_ptr<ResponseTree> fResponseTree; ///< Handler of ROOT tree output.
  
  // --- END Internal variables ------------------------------------------------

  
  // @{
  /// Access to the helper.
  SlidingWindowTriggerEfficiencyPlots const& helper() const { return *this; }
  SlidingWindowTriggerEfficiencyPlots& helper() { return *this; }
  // @}
  

  // --- BEGIN Derived class methods -------------------------------------------
  /**
   * @brief Initializes full set of plots for (ADC threshold + category).
   * 
   * This customization of `TriggerEfficiencyPlotsBase::initializePlotSet()`
   * adds some trigger-definition specific plots and some overview plots
   * across different trigger definitions.
   */
  virtual void initializePlotSet
    (PlotSandbox& plots, std::vector<SettingsInfo_t> const& settings) const
    override;
  
  /**
   * @brief Simulates all trigger minimum requirements plots the results.
   * @param thresholdIndex the index of the PMT threshold of input primitives
   * @param gates the trigger primitives used to simulate the trigger response
   * @param eventInfo general information about the event being simulated
   * @param selectedPlots list of boxes containing plots to be filled
   * 
   * This method is expected to perform the following steps for each trigger
   * primitive requirement in `MinimumPrimitives`:
   * 
   * 1. combine the trigger primitives: `combineTriggerPrimitives()`;
   * 2. apply the beam gate: `applyBeamGateToAll()` on the combined primitives;
   * 3. generate the trigger response: in `plotResponse()`;
   * 4. fill all plots: also in in `plotResponse()`.
   * 
   * Details are in the documentation of the relevant methods.
   * 
   * This method is invoked once per PMT threshold.
   */
  virtual void simulateAndPlot(
    std::size_t const thresholdIndex,
    TriggerGatesPerCryostat_t const& gates,
    EventInfo_t const& eventInfo,
    detinfo::DetectorClocksData const& clockData,
    PlotSandboxRefs_t const& selectedPlots
    ) override;
    
  // --- END Derived class methods ---------------------------------------------

  /**
   * @brief Fills plots with the specified trigger response.
   * @param iThr index of PMT threshold (used in tree output)
   * @param threshold PMT threshold tag (for printing)
   * @param iPattern index of the pattern being plotted
   * @param pattern the pattern being plotted
   * @param plotSets set of plot boxes to fill (from `initializePlotSet()`)
   * @param eventInfo event information for plotting
   * @param triggerInfo the information about the response of this trigger
   * 
   * This method fills all the relevant plots for the specified trigger pattern
   * and threshold. The trigger response is passed as a parameter.
   */
  void plotResponse(
    std::size_t iThr, std::string const& threshold,
    std::size_t iPattern, WindowPattern const& pattern,
    PlotSandboxRefs_t const& plotSets,
    EventInfo_t const& eventInfo,
    PMTInfo_t const& PMTinfo,
    WindowTriggerInfo_t const& triggerInfo
    ) const;

  /// Constructs all the pattern algorithms.
  /// Must be called after setting the window topology.
  void initializePatternAlgorithms();

  /// Fills all event plots with data from `eventInfo` as in `fillEventPlots()`.
  void fillAllEventPlots
    (PlotSandboxRefs_t const& plotSets, EventInfo_t const& eventInfo) const;

  /// Fills all PMY plots with data from `PMTinfo` as in `fillPMTplots()`.
  void fillAllPMTplots
    (PlotSandboxRefs_t const& plotSets, PMTInfo_t const& PMTinfo) const;

}; // icarus::trigger::SlidingWindowTriggerEfficiencyPlots



//------------------------------------------------------------------------------
//---  Implementation
//------------------------------------------------------------------------------
//--- Local namespace
//------------------------------------------------------------------------------
namespace {
  template <typename DestColl, typename SrcColl>
  DestColl& appendCollection(DestColl& dest, SrcColl&& src) {
    dest.insert(
      dest.end(), std::move_iterator(src.begin()), std::move_iterator(src.end())
      );
    return dest;
  } // appendCollection()
} // local namespace


//------------------------------------------------------------------------------
//--- ResponseTree
//------------------------------------------------------------------------------
template <typename Thresholds, typename Settings>
ResponseTree::ResponseTree
  (TTree& tree, Thresholds const& thresholds, Settings const& settings)
  : TreeHolder(tree)
  , indices(std::size(thresholds), std::size(settings))
  , RespTxxSxx{ std::make_unique<bool[]>(indices.size()) }
{

  for (auto [ iThr, thresholdTag]: util::enumerate(thresholds)) {

    for (auto [ iSetting, setting ]: util::enumerate(settings)) {

      std::string const branchName
        = "RespT" + thresholdTag + "S" + util::to_string(setting);

      this->tree().Branch
        (branchName.c_str(), &(RespTxxSxx[indices(iThr, iSetting)]));

    } // for all requirements

  } // for all thresholds

} // ResponseTree::ResponseTree()


//------------------------------------------------------------------------------
void ResponseTree::assignResponse
  (std::size_t iThr, std::size_t iSettings, bool resp)
{
  RespTxxSxx[indices(iThr, iSettings)] = resp;
} // ResponseTree::assignResponse()


//------------------------------------------------------------------------------
//--- icarus::trigger::SlidingWindowTriggerEfficiencyPlots
//------------------------------------------------------------------------------
icarus::trigger::SlidingWindowTriggerEfficiencyPlots::SlidingWindowTriggerEfficiencyPlots
  (Parameters const& config)
  : art::EDAnalyzer           (config)
  , TriggerEfficiencyPlotsBase(config(), consumesCollector())
  // configuration
  , fPatterns(config().Patterns())
  // internal variables
  , fWindowMapMan{
      helper().geometry(),
      helper().logCategory() + "_WindowMapManager"
      }
{
  
  if (fPatterns.empty()) {
    throw art::Exception(art::errors::Configuration)
      << "At least one 'MinimumPrimitives' requirement... required.";
  }
  
  std::size_t iPattern [[maybe_unused]] = 0U; // NOTE: incremented only in DEBUG
  for (auto const& pattern: fPatterns) {
    std::size_t const index [[maybe_unused]]
      = createCountersForPattern(pattern.tag());
    assert(index == iPattern++);
    
  } // for patterns
  
  //
  // more complex parameter parsing
  //
  if (helper().eventTree()) {

    fResponseTree = std::make_unique<ResponseTree>
      (*(helper().eventTree()), helper().ADCthresholds(), fPatterns);

  } // if make tree

  {
    mf::LogInfo log(helper().logCategory());
    
    log
      << "Requirement of sliding window patterns ("
      << fPatterns.size() << "):";
    for (auto const& pattern: fPatterns)
      log << "\n [" << pattern.tag() << "] " << pattern.description();
    
  }
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::SlidingWindowTriggerEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::beginJob() {

  // hook helper and framework;
  // actual action is in the (overridden) virtual functions `initializeXxx()`.
  
  // NOTE this action can't happen in constructor because relies on polymorphism
  std::vector<SettingsInfo_t> settings;
  for (auto const& [ iPattern, pattern ]: util::enumerate(fPatterns)) {
    settings.emplace_back(
      iPattern,             // index
      pattern.tag(),        // tag
      pattern.description() // description
      );
  } // for
  
  // we use the default plot categories defined in
  // `TriggerEfficiencyPlotsBase::DefaultPlotCategories`
  helper().initializePlots(settings);
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::beginJob()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::analyze
  (art::Event const& event)
{
  
  // hook helper and framework;
  // actual action is in the (overridden) virtual function `plotResponse()`
  // and `fillXxxPlots()`.
  helper().process(event);
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::analyze()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::endJob() {
  
  helper().deleteEmptyPlots(); // don't keep plots with no entries
  
  // hook helper and framework
  helper().printSummary();
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::endJob()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::initializePlotSet
  (PlotSandbox& plots, std::vector<SettingsInfo_t> const& settings) const
{

  //
  // Selection-related plots
  //
  
  // (inherited)
  helper().TriggerEfficiencyPlotsBase::initializePlotSet(plots, settings);
  
  //
  // overview plots with different settings
  //
  
  std::vector<std::string> patternLabels; // each pattern has tag as label
  std::transform(
    fPatterns.cbegin(), fPatterns.cend(), back_inserter(patternLabels),
    std::mem_fn(&WindowPattern::tag)
    );
  
  //
  // Triggering efficiency vs. requirements.
  //
  auto const [ detTimings, beamGate, preSpillWindow ] = makeGatePack();
  detinfo::timescales::optical_time_ticks const triggerResolutionTicks
    { detTimings.toOpticalTicks(helper().triggerTimeResolution()) };
  
  auto const& beamGateOpt = beamGate.asOptTickRange();
  
  auto* TrigTime = plots.make<TH2F>(
    "TriggerTick",
    "Trigger time tick"
      ";window pattern"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]",
    fPatterns.size(), 0.0, double(fPatterns.size()),
    beamGateOpt.duration() / triggerResolutionTicks,
    beamGateOpt.first.value(), beamGateOpt.second.value()
    );
  
  util::ROOT::applyAxisLabels(TrigTime->GetXaxis(), patternLabels);
  
  auto* Triggers = plots.make<TH1F>(
    "Triggers",
    "Triggered events"
      ";window pattern"
      ";triggered events",
    fPatterns.size(), 0.0, double(fPatterns.size())
    );
  
  util::ROOT::applyAxisLabels(Triggers->GetXaxis(), patternLabels);
  
  auto* Eff = plots.make<TEfficiency>(
    "Eff",
    "Efficiency of triggering"
      ";window pattern"
      ";trigger efficiency",
    fPatterns.size(), 0.0, double(fPatterns.size())
//    fMinimumPrimitives.back(), 0, fMinimumPrimitives.back() + 1
    );
  
  // people are said to have earned hell for things like this;
  // but TEfficiency really does not expose the interface to assign labels to
  // its axes, which supposedly could be done had we chosen to create it by
  // histograms instead of directly as recommended.
  // Also need to guess which is the relevant histogram.
  util::ROOT::applyAxisLabels
    (const_cast<TH1*>(Eff->GetTotalHistogram())->GetXaxis(), patternLabels);
  util::ROOT::applyAxisLabels
    (const_cast<TH1*>(Eff->GetPassedHistogram())->GetXaxis(), patternLabels);
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::initializePlotSet()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::simulateAndPlot(
  std::size_t const thresholdIndex,
  TriggerGatesPerCryostat_t const& gates,
  EventInfo_t const& eventInfo,
  detinfo::DetectorClocksData const& clockData,
  PlotSandboxRefs_t const& selectedPlots
) {
  
  auto const threshold = helper().ADCthresholdTag(thresholdIndex);
  
  /*
   * 0. initialize or verify the topology of the input
   * 1. apply the beam gate to each input gate
   * 2. for each pattern: apply the pattern, plot the trigger outcome
   * 3. fill all trigger-independent plots
   * 4. fill all PMT plots (threshold-dependent)
   */
  
  //
  // 0. initialize or verify the topology of the input
  //
  
  // throws exception on verification failure;
  // pattern algorithms are constructed here because they require window mapping
  if (fWindowMapMan(gates)) initializePatternAlgorithms();

  
  auto const& beamGate = helper().makeMyBeamGate(clockData);
  
  //
  // 1. apply the beam gate to each input gate
  //    (it's ok to lose provenance information since we have the map)
  //
  TriggerGates_t inBeamGates;
  for (auto const& cryoGates: gates)
    appendCollection(inBeamGates, beamGate.applyToAll(cryoGates));
  
  // --- BEGIN DEBUG -----------------------------------------------------------
  {
    mf::LogTrace log(helper().logCategory());
    log << "Input for threshold " << threshold << ": " << inBeamGates.size()
      << " primitives. After beam gate:";
    unsigned int nOpen = 0U;
    using icarus::trigger::gatesIn;
    for (auto const& [ iWindow, gate ]: util::enumerate(gatesIn(inBeamGates))) {
      auto const maxTick = gate.findMaxOpen();
      if (maxTick == gate.MinTick) continue;
      ++nOpen;
      log << "\n  window #" << iWindow << ": maximum "
        << gate.openingCount(maxTick) << " at tick " << maxTick;
    } // for
    if (!nOpen) log << "  nothing.";
  }
  // --- END DEBUG -------------------------------------------------------------
  
  // get which gates are active during the beam gate
  PMTInfo_t const PMTinfo
    { threshold, helper().extractActiveChannels(gates) };
  
  //
  // 2. for each pattern:
  //
  for (auto const& [ iPattern, pattern ]: util::enumerate(fPatterns)) {

    auto& patternAlg = fPatternAlgs[iPattern];
    
    WindowTriggerInfo_t const triggerInfo
      = patternAlg.simulateResponse(inBeamGates);
    
    registerTriggerResult(thresholdIndex, iPattern, triggerInfo.info);

    plotResponse(
      thresholdIndex, threshold,
      iPattern, pattern,
      selectedPlots,
      eventInfo, PMTinfo, triggerInfo
      );
    
  } // for window patterns
  
  //
  // 3. fill all trigger-independent plots (one copy per threshold... meh)
  //
  fillAllEventPlots(selectedPlots, eventInfo);
  
  //
  // 4. fill all PMT plots (threshold-dependent)
  //
  fillAllPMTplots(selectedPlots, PMTinfo);
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::simulateAndPlot()


//------------------------------------------------------------------------------
void
icarus::trigger::SlidingWindowTriggerEfficiencyPlots::initializePatternAlgorithms
  ()
{
  fPatternAlgs.clear();
  for (auto const& pattern: fPatterns)
    fPatternAlgs.emplace_back(*fWindowMapMan, pattern, helper().logCategory());
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::initializePatternAlgorithms()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::fillAllEventPlots
  (PlotSandboxRefs_t const& plotSets, EventInfo_t const& eventInfo) const
{
  /*
   * Now fill the plots independent of the trigger response:
   * the same value is plotted in all plot sets.
   * (again for all pertinent event categories, e.g. charged currents, etc.)
   */
  for (PlotSandbox const& plotSet: plotSets) {
    
    //
    // general plots, independent of trigger definition details
    //
    fillEventPlots(eventInfo, plotSet);
    
  } // for
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::fillAllEventPlots()



//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::fillAllPMTplots
  (PlotSandboxRefs_t const& plotSets, PMTInfo_t const& PMTinfo) const
{
  /*
   * Now fill the plots independent of the trigger response:
   * the same value is plotted in all plot sets.
   * (again for all pertinent event categories, e.g. charged currents, etc.)
   */
  for (PlotSandbox const& plotSet: plotSets) {
    
    //
    // general plots, independent of trigger definition but dependent on
    // threshold
    //
    fillPMTplots(PMTinfo, plotSet);
    
  } // for
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::fillAllPMTplots()


void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::plotResponse(
  std::size_t iThr, std::string const& threshold,
  std::size_t iPattern, WindowPattern const& pattern,
  PlotSandboxRefs_t const& plotSets,
  EventInfo_t const& eventInfo,
  PMTInfo_t const& PMTinfo,
  WindowTriggerInfo_t const& triggerInfo
) const {
  
  using namespace std::string_literals;
  
  bool const fired = triggerInfo.info.fired();
  
  if (fResponseTree) fResponseTree->assignResponse(iThr, iPattern, fired);
  
  
  std::string const patternTag { pattern.tag() };
  
  // go through all the plot categories this event qualifies for
  // (for example: charged currents, muon neutrinos, ...)
  for (icarus::trigger::PlotSandbox const& plotSet: plotSets) {
    
    //
    // overview plots from different thresholds
    //
    
    HistGetter const get { plotSet };
    
    // simple efficiency
    get.Eff("Eff"s).Fill(fired, iPattern);
    
    // simple count
    if (fired) get.Hist("Triggers"s).Fill(iPattern);
    
    // trigger time (if any)
    if (fired) {
      get.Hist2D("TriggerTick"s).Fill
        (iPattern, triggerInfo.info.atTick().value());
    }
    
    //
    // plots depending on the trigger response
    // (but not caring of the trigger definition details)
    //
    
    // efficiency plots
    // (including event plots in the triggered or non-triggered category)
    helper().fillAllEfficiencyPlots
      (eventInfo, PMTinfo, triggerInfo.info, plotSet.demandSandbox(patternTag));
    
    //
    // add here further trigger-specific plots
    //
    
  } // for all qualifying plot categories
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::plotResponse()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::SlidingWindowTriggerEfficiencyPlots)


//------------------------------------------------------------------------------
