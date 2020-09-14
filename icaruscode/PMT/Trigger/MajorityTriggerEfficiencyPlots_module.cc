/**
 * @file   MajorityTriggerEfficiencyPlots_module.cc
 * @brief  Plots of efficiency for triggers based on PMT channel global count.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 9, 2020
 * @see    icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h"
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h" // sumGates()
#include "icaruscode/PMT/Trigger/Utilities/ROOTutils.h" // util::ROOT

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/TensorIndices.h" // util::MatrixIndices
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_time_ticks..
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"

// ROOT libraries
#include "TTree.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <algorithm> // std::sort()
#include <vector>
#include <memory> // std::unique_ptr
#include <utility> // std::pair<>, std::move()
#include <cstddef> // std::size_t


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
 * The branch structure is: a `RespTxxRxx/O` branch for each threshold (numeric)
 * and requirement (also numeric).
 * with a single branch per element.
 *
 */
struct ResponseTree: public icarus::trigger::details::TreeHolder {

  // `std::vector<bool>` is too special for us. Let's pack it to go.
  util::MatrixIndices indices;
  std::unique_ptr<bool[]> RespTxxRxx;

  
  /// Constructor: accommodates that many thresholds and requirements.
  template <typename Thresholds, typename Requirements>
  ResponseTree
    (TTree& tree, Thresholds const& thresholds, Requirements const& minReqs);

  /// Assigns the response for the specified trigger.
  void assignResponse(std::size_t iThr, std::size_t iReq, bool resp);

}; // struct ResponseTree


// --- END -- ROOT tree helpers ------------------------------------------------


//------------------------------------------------------------------------------
namespace icarus::trigger { class MajorityTriggerEfficiencyPlots; }
/**
 * @brief Produces plots about trigger simulation and trigger efficiency.
 * 
 * This module is an implementation of `TriggerEfficiencyPlotsBase`
 * for a trigger defined as a minimum number of trigger primitives beyond
 * threshold.
 * 
 * A trigger primitive is a two-level function of time which describes when
 * that primitive is on and when it is off. Trigger primitives are given as
 * input to this module and their origin may vary, but the standard source in
 * ICARUS is the pairing with AND or OR of two optical detector channels
 * discriminated against a certain ADC count threshold.
 * 
 * 
 * Trigger logic algorithm
 * ========================
 * 
 * @anchor MajorityTriggerEfficiencyPlots_Algorithm
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::MajorityTriggerEfficiencyPlots` and its assumptions.
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
 * The algorithm handles multiple trigger primitive requirements.
 * Each requirement is simply how many trigger primitives must be open at the
 * same time in a single cryostat for the trigger to fire. The values of these
 * requirements are set in the configuration (`MinimumPrimitives`).
 * To determine whether a trigger with a given requirement, i.e. with a required
 * minimum number of trigger primitives open at the same time, has fired, the
 * gate combined as described above is scanned to find _the first tick_ where
 * the level of the gate reaches or passes this minimum required number. If such
 * tick exists, the trigger is considered to have fired, and at that time.
 * 
 * As a consequence, there are for each event as many different trigger
 * responses as how many different requirements are configured in
 * `MinimumPrimitives`, _times_ how many ADC thresholds are provided in input,
 * configured in `Thresholds`.
 * 
 * While there _is_ a parameter describing the time resolution of the trigger
 * (`TriggerTimeResolution`), this is currently only used for aesthetic purposes
 * to choose the binning of some plots: the resolution is _not_ superimposed
 * to the gates (yet).
 * 
 * The combination algorithm is implemented in `combineTriggerPrimitives()`
 * while the requirement evaluation and plotting are implemented in
 * `plotResponses()`.
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
 * in this module represented by the single trigger primitive requirement _N_
 * (configuration parameter: `MinimumPrimitives`). The folders and plots will
 * identify each requirement with the tag `ReqN` (e.g. `Req5` when requesting
 * at least 5 trigger primitives for an event trigger).
 * 
 * There are a few plots that are produced by this module in addition to the
 * ones in `TriggerEfficiencyPlotsBase`.
 * 
 * There are different "types" of plots. Some
 * @ref MajorityTriggerEfficiencyPlots_SelectionPlots "do not depend on triggering at all",
 * like the deposited energy distribution. Others
 * @ref MajorityTriggerEfficiencyPlots_MultiTriggerPlots "cross different trigger definitions",
 * like the trigger efficiency as function of trigger requirement. Others still
 * @ref MajorityTriggerEfficiencyPlots_SingleTriggerPlots "assume a single trigger definition":
 * this is the case of trigger efficiency plots versus energy. Finally, there are
 * @ref MajorityTriggerEfficiencyPlots_SingleTriggerResponsePlots "plots that depend on a specific trigger definition and outcome":
 * this is the case of all the plots including only triggering or non-triggering
 * events.
 * 
 * A list of additional plots follows for each plot type.
 * All the plots are always relative to a specific optical detector channel
 * threshold (ADC) and a broad event category.
 * 
 * 
 * ### Plots independent of the triggers (selection plots)
 * 
 * @anchor MajorityTriggerEfficiencyPlots_SelectionPlots
 * 
 * Only @ref TriggerEfficiencyPlotsBase_SelectionPlots "the standard plots"
 * from `TriggerEfficiencyPlotsBase` are produced in this category.
 * 
 * 
 * ### Plots including different trigger requirements
 * 
 * @anchor MajorityTriggerEfficiencyPlots_MultiTriggerPlots
 * 
 * In addition to @ref TriggerEfficiencyPlotsBase_MultiTriggerPlots "the plots"
 * from `TriggerEfficiencyPlotsBase`, the following plots are also produced:
 * 
 * * `Eff`: trigger efficiency defined as number of triggered events over the
 *   total number of events, as function of the minimum number of trigger
 *   primitives (`MinimumPrimitives`) to define a firing trigger; uncertainties
 *   are managed by `TEfficiency`.
 * * `TriggerTick`: distribution of the time of the earliest trigger for the
 *   event, as function of the minimum number of trigger primitives (as in
 *   `Eff`). It may happen that the event is such that there is e.g. a
 *   20-primitive flash, then subsiding, and then another 30-primitive flash.
 *   In such a case, in the trigger requirement "&geq; 15 primitives" such event
 *   will show at the time of the 20-primitive flash, while in the trigger
 *   requirement "&geq; 25 primitives" it will show at the time of the 
 *   30-primitive flash. Each event appears at most once for each trigger
 *   requirement, and it may not appear at all if does not fire a trigger.
 * * `NPrimitives`: the maximum number of primitives "on" at any time.
 * 
 * 
 * ### Plots depending on a specific trigger definition
 * 
 * @anchor MajorityTriggerEfficiencyPlots_SingleTriggerPlots
 * 
 * Only @ref TriggerEfficiencyPlotsBase_SingleTriggerPlots "the standard plots"
 * from `TriggerEfficiencyPlotsBase` are produced in this category.
 * 
 *
 * ### Plots depending on a specific trigger definition and response
 *
 * @anchor MajorityTriggerEfficiencyPlots_SingleTriggerResponsePlots
 * 
 * Only @ref TriggerEfficiencyPlotsBase_SingleTriggerResponsePlots "the standard plots"
 * from `TriggerEfficiencyPlotsBase` are produced in this category.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * @anchor MajorityTriggerEfficiencyPlots_Configuration
 * 
 * In addition to
 * @ref TriggerEfficiencyPlotsBase_Configuration "all the configuration parameters"
 * from `TriggerEfficiencyPlotsBase`, the following one is also present:
 * 
 * * `MinimumPrimitives` (list of integers, _mandatory_): a list of alternative
 *     requirements for the definition of a trigger; each value is the number
 *     of trigger primitives needed to be "on" at the same time for the trigger
 *     to fire;
 * 
 * An example job configuration is provided as `maketriggerplots_icarus.fcl`.
 * 
 * 
 * Technical description of the module
 * ====================================
 * 
 * This module class is derived from
 * `icarus::trigger::TriggerEfficiencyPlotsBase`, which provides a backbone to
 * perform the simulation of triggers and plotting of their efficiency.
 * 
 * There is not any superior design involved in this separation, but just the
 * desire to share most code possible between different modules which simulate
 * different trigger patterns and as a consequence might have specific plots to
 * fill.
 * 
 * This module redefines:
 * 
 * * `initializePlotSet()` to define the
 *   @ref MajorityTriggerEfficiencyPlots_MultiTriggerPlots "few additional plots"
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
class icarus::trigger::MajorityTriggerEfficiencyPlots
  : public art::EDAnalyzer
  , private icarus::trigger::TriggerEfficiencyPlotsBase
{

    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config: public TriggerEfficiencyPlotsBase::Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<unsigned int> MinimumPrimitives {
      Name("MinimumPrimitives"),
      Comment("minimum required number of trigger primitives for a trigger")
      };
    
  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit MajorityTriggerEfficiencyPlots(Parameters const& config);

  // Plugins should not be copied or assigned.
  MajorityTriggerEfficiencyPlots(MajorityTriggerEfficiencyPlots const&) = delete;
  MajorityTriggerEfficiencyPlots(MajorityTriggerEfficiencyPlots&&) = delete;
  MajorityTriggerEfficiencyPlots& operator=(MajorityTriggerEfficiencyPlots const&) = delete;
  MajorityTriggerEfficiencyPlots& operator=(MajorityTriggerEfficiencyPlots&&) = delete;

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
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Minimum number of trigger primitives for a trigger to happen.
  std::vector<unsigned int> fMinimumPrimitives;
  
  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Internal variables ----------------------------------------------
  
  std::unique_ptr<ResponseTree> fResponseTree; ///< Handler of ROOT tree output.
  
  // --- END Internal variables ------------------------------------------------

  
  // @{
  /// Access to the helper.
  MajorityTriggerEfficiencyPlots const& helper() const { return *this; }
  MajorityTriggerEfficiencyPlots& helper() { return *this; }
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
   * 2. apply the beam gate: `applyBeamGate()` on the combined primitives;
   * 3. generate the trigger response: in `plotResponses()`;
   * 4. fill all plots: also in in `plotResponses()`.
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
    ) const override;
    
  // --- END Derived class methods ---------------------------------------------

  
  /**
   * @brief Completes the event trigger simulation and fills the plots.
   * @param iThr index of PMT threshold (used in tree output)
   * @param threshold PMT threshold in ADC counts (for printing)
   * @param plotSets set of plot boxes to fill (from `initializePlotSet()`)
   * @param eventInfo event information for plotting
   * @param combinedTrigger combined trigger primitive
   * 
   * For each of the trigger requirements (`MinimumPrimitives`), this method:
   * 
   * 1. applies the requirement to the `combinedTrigger` trigger primitive
   * 2. computes the event trigger
   * 3. fills all plots in all the plot sets for this requirement accordingly
   * 
   * The input combined trigger primitive contains the maximum number of
   * trigger primitives active at each optical clock tick.
   * It is assumed that the beam gate has already been "applied" so that outside
   * it no trigger primitive is considered open.
   * 
   * A trigger with requirement of minimum trigger primitives _N_ is fired if
   * this combined primitive reaches or passes _N_, i.e. if there are at least
   * _N_ open trigger primitives at any time.
   * The time of the trigger is taken as the first tick at which the requirement
   * is met.
   * 
   * Extra plots are an overview of the trigger efficiency for different
   * requirements, per threshold, the distribution of the trigger time for
   * different requirements, and maximum number of open primitive in a cryostat,
   * per event.
   * 
   * Note that there is no information about which cryostat is meeting the
   * trigger requirement.
   */
  void plotResponses(
    std::size_t iThr, ADCCounts_t const threshold,
    PlotSandboxRefs_t const& plotSets, EventInfo_t const& eventInfo,
    detinfo::DetectorClocksData const& clockData,
    TriggerGateData_t const& combinedTrigger
    ) const;
  
  /**
   * @brief Computes the trigger response from primitives with the given
   *        `threshold`.
   * @param cryoGates collections of trigger primitives, one coll. per cryostat
   * @param threshold PMT threshold of the primitives (for printing purposes)
   * @return a list of combined trigger primitives, one combination per cryostat
   * 
   * The input trigger primitives are already grouped by cryostat.
   * For each cryostat, the primitives are combined and one combination is
   * returned. The combination is just the "total" of the primitives opened at
   * each tick.
   * 
   * The event trigger is not finalized here, and the cryostat trigger
   * primitives are all returned.
   */
  TriggerGateData_t combineTriggerPrimitives(
    TriggerGatesPerCryostat_t const& cryoGates,
    ADCCounts_t const threshold
    ) const;

  
}; // icarus::trigger::MajorityTriggerEfficiencyPlots



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- ResponseTree
//------------------------------------------------------------------------------
template <typename Thresholds, typename Requirements>
ResponseTree::ResponseTree
  (TTree& tree, Thresholds const& thresholds, Requirements const& minReqs)
  : TreeHolder(tree)
  , indices(std::size(thresholds), std::size(minReqs))
  , RespTxxRxx{ std::make_unique<bool[]>(indices.size()) }
{

  for (auto [ iThr, threshold]: util::enumerate(thresholds)) {
    std::string const thrStr = util::to_string(raw::ADC_Count_t(threshold));

    for (auto [ iReq, req ]: util::enumerate(minReqs)) {

      std::string const branchName
        = "RespT" + thrStr + "R" + util::to_string(req);

      this->tree().Branch
        (branchName.c_str(), &(RespTxxRxx[indices(iThr, iReq)]));

    } // for all requirements

  } // for all thresholds

} // ResponseTree::ResponseTree()


//------------------------------------------------------------------------------
void ResponseTree::assignResponse
  (std::size_t iThr, std::size_t iReq, bool resp)
{
  RespTxxRxx[indices(iThr, iReq)] = resp;
} // ResponseTree::assignResponse()


//------------------------------------------------------------------------------
//--- icarus::trigger::MajorityTriggerEfficiencyPlots
//------------------------------------------------------------------------------
icarus::trigger::MajorityTriggerEfficiencyPlots::MajorityTriggerEfficiencyPlots
  (Parameters const& config)
  : art::EDAnalyzer           (config)
  , TriggerEfficiencyPlotsBase(config(), consumesCollector())
  // configuration
  , fMinimumPrimitives(config().MinimumPrimitives())
{
  std::sort(fMinimumPrimitives.begin(), fMinimumPrimitives.end());
  
  //
  // more complex parameter parsing
  //
  if (helper().eventTree()) {

    fResponseTree = std::make_unique<ResponseTree>
      (*(helper().eventTree()), helper().ADCthresholds(), fMinimumPrimitives);

  } // if make tree

  if (fMinimumPrimitives.empty()) {
    throw art::Exception(art::errors::Configuration)
      << "At least one 'MinimumPrimitives' requirement... required.";
  }
  
  {
    mf::LogInfo log(helper().logCategory());
    log
      << "Requirement of minimum trigger primitives ("
      << fMinimumPrimitives.size() << "):";
    for (auto const& req: fMinimumPrimitives) log << " " << req;
    log << ".";
  }
  
} // icarus::trigger::MajorityTriggerEfficiencyPlots::MajorityTriggerEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerEfficiencyPlots::beginJob() {

  // hook helper and framework;
  // actual action is in the (overridden) virtual functions `initializeXxx()`.
  
  // NOTE this action can't happen in constructor because relies on polymorphism
  std::vector<SettingsInfo_t> settings;
  for (auto [ iReq, minCount ]: util::enumerate(fMinimumPrimitives)) {
    std::string const minCountStr { std::to_string(minCount) };
    settings.emplace_back(
      iReq,                              // index
      "Req" + minCountStr,               // tag
      minCountStr + " channels required" // description
      );
  } // for
  
  // we use the default plot categories defined in
  // `TriggerEfficiencyPlotsBase::DefaultPlotCategories`
  helper().initializePlots(settings);
  
} // icarus::trigger::MajorityTriggerEfficiencyPlots::beginJob()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerEfficiencyPlots::analyze
  (art::Event const& event)
{
  
  // hook helper and framework;
  // actual action is in the (overridden) virtual function `plotResponses()`
  // and `fillXxxPlots()`.
  helper().process(event);
  
} // icarus::trigger::MajorityTriggerEfficiencyPlots::analyze()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerEfficiencyPlots::endJob() {
  
  // hook helper and framework
  helper().printSummary();
  
} // icarus::trigger::MajorityTriggerEfficiencyPlots::endJob()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerEfficiencyPlots::initializePlotSet
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
  
  // a variable binning for the required number of trigger primitives
  auto [ minimumPrimBinning, minimumPrimBinningLabels ]
    = util::ROOT::makeVariableBinningAndLabels(fMinimumPrimitives);
  assert(minimumPrimBinning.size() == minimumPrimBinningLabels.size() + 1U);

  {
    mf::LogTrace log(helper().logCategory());
    log << "MajorityTriggerEfficiencyPlots (plots '"
      << plots.name() << "') variable binning including the "
      << fMinimumPrimitives.size() << " points {";
    for (auto value: fMinimumPrimitives) log << " " << value;
    log << " } => " << minimumPrimBinningLabels.size() << " bins: ";
    for (auto const& [ value, label ]
      : util::zip<1U>(minimumPrimBinning, minimumPrimBinningLabels))
    {
      log << " " << value << " (\"" << label << "\") =>";
    } // for
    log << " " << minimumPrimBinning.back();
  } // debug output block

  //
  // Triggering efficiency vs. requirements.
  //

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  detinfo::timescales::optical_time_ticks const triggerResolutionTicks
    { helper().detTimings(clockData).toOpticalTicks(helper().triggerTimeResolution()) };
  auto const& beamGateOpt = helper().beamGateTickRange();
  auto* TrigTime = plots.make<TH2F>(
    "TriggerTick",
    "Trigger time tick"
      ";minimum requested number of trigger primitives"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]",
    minimumPrimBinning.size() - 1U, minimumPrimBinning.data(),
//    fMinimumPrimitives.back(), 0, fMinimumPrimitives.back() + 1
    (beamGateOpt.second - beamGateOpt.first) / triggerResolutionTicks,
    beamGateOpt.first.value(), beamGateOpt.second.value()
    );
  
  util::ROOT::applyAxisLabels(TrigTime->GetXaxis(), minimumPrimBinningLabels);
  
  auto* Eff = plots.make<TEfficiency>(
    "Eff",
    "Efficiency of triggering"
      ";minimum requested number of trigger primitives"
      ";trigger efficiency",
    minimumPrimBinning.size() - 1U, minimumPrimBinning.data()
//    fMinimumPrimitives.back(), 0, fMinimumPrimitives.back() + 1
    );
  
  // people are said to have earned hell for things like this;
  // but TEfficiency really does not expose the interface to assign labels to
  // its axes, which supposedly could be done had we chosen to create it by
  // histograms instead of directly as recommended.
  util::ROOT::applyAxisLabels(
    const_cast<TH1*>(Eff->GetTotalHistogram())->GetXaxis(),
    minimumPrimBinningLabels
    );
  
  //
  // plots independent of the trigger primitive requirements
  //
  plots.make<TH1F>(
    "NPrimitives",
    "Number of trigger primitives (\"channels firing at once\")"
    ";maximum trigger primitives at the same time on a single cryostat"
    ";events",
    192, 0.0, 192.0 // large number, zoom in presentations!
    );
  
  
} // icarus::trigger::MajorityTriggerEfficiencyPlots::initializePlotSet()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerEfficiencyPlots::simulateAndPlot(
  std::size_t const thresholdIndex,
  TriggerGatesPerCryostat_t const& gates,
  EventInfo_t const& eventInfo,
  detinfo::DetectorClocksData const& clockData,
  PlotSandboxRefs_t const& selectedPlots
) const {
  
  auto const threshold = helper().ADCthreshold(thresholdIndex);

  /* 
   * 1. combine the trigger primitives (`combineTriggerPrimitives()`)
   * 2. apply the beam gate on the combination (`applyBeamGate()`)
   * 3. and compute the trigger response (`plotResponses()`)
   * 4. fill plots with the result (also `plotResponses()`)
   */
  plotResponses(
    thresholdIndex, threshold, selectedPlots, eventInfo,
    clockData,
    helper().applyBeamGate(combineTriggerPrimitives(gates, threshold))
    );
  
} // icarus::trigger::MajorityTriggerEfficiencyPlots::simulateAndPlot()


//------------------------------------------------------------------------------
void icarus::trigger::MajorityTriggerEfficiencyPlots::plotResponses(
  std::size_t iThr,
  icarus::trigger::ADCCounts_t const threshold,
  PlotSandboxRefs_t const& plotSets,
  EventInfo_t const& eventInfo,
  detinfo::DetectorClocksData const& clockData,
  TriggerGateData_t const& combinedCount
) const {
  
  /*
   * This function plots according to the configured minimum number of trigger
   * primitives: for each requirement of minimum number of primitives, the
   * earliest time where that requirement is met is found, and that is
   * considered as the trigger time.
   * 
   * The following quantities are drawn per ADC threshold and per plot category:
   * 
   * * per minimum number of primitives:
   *    * trigger efficiency (i.e. whether there was _any time_ a number of
   *      primitives fulfilling the requirement)
   *    * trigger time in ticks (distribution as 2D histogram)
   * * maximum number of trigger primitives present at any time
   * * deposited energy during beam spill
   * 
   */
  using namespace std::string_literals;
  
  using ClockTick_t = TriggerGateData_t::ClockTick_t;
  using OpeningCount_t = TriggerGateData_t::OpeningCount_t;
  
  using PrimitiveCount_t = std::pair<ClockTick_t, OpeningCount_t>;
  
  auto const maxPrimitiveTime { combinedCount.findMaxOpen() };
  PrimitiveCount_t const maxPrimitives
    { maxPrimitiveTime, combinedCount.openingCount(maxPrimitiveTime) };

  mf::LogTrace(helper().logCategory())
    << "Max primitive count in " << threshold << ": "
    << maxPrimitives.second << " at tick " << maxPrimitives.first << " ("
    << helper().detTimings(clockData).toElectronicsTime
      (detinfo::DetectorTimings::optical_tick{ maxPrimitives.first })
    << ")"
    ;
  
  /*
   * Fill all the histograms for all the minimum primitive requirements
   * (filling the information whether or not the trigger fired),
   * for all the qualifying categories.
   * Note that in this type of plots each event appears in all bins
   * (may be with "fired" or "not fired" on each bin)
   */
  PrimitiveCount_t lastMinCount { TriggerGateData_t::MinTick, 0 };
  bool fired = true; // the final trigger response (changes with requirement)
  
  
  for (auto [ iReq, minCount ]: util::enumerate(fMinimumPrimitives)) {
    
    // in this check, `fired` remembers the outcome from the previous threshold
    if (fired && (lastMinCount.second < minCount)) {
      // if we haven't passed this minimum yet
      ClockTick_t const time = combinedCount.findOpen(minCount);
      if (time == TriggerGateData_t::MaxTick) {
        mf::LogTrace(helper().logCategory())
          << "Never got at " << minCount << " primitives or above.";
        fired = false;
      }
      else {
        lastMinCount = { time, combinedCount.openingCount(time) };
        mf::LogTrace(helper().logCategory())
          << "Reached " << minCount << " primitives or above ("
          << lastMinCount.second << ") at " << lastMinCount.first << ".";
      }
    } // if
    
    TriggerInfo_t triggerInfo;
    if (fired) triggerInfo.emplace(optical_tick{ lastMinCount.first });
    
    // at this point we know we have minCount or more trigger primitives,
    // and the time of this one is in lastMinCount.first (just in case)
    
    if (fResponseTree) fResponseTree->assignResponse(iThr, iReq, fired);
    
    std::string const minCountStr { "Req" + std::to_string(minCount) };
    
    // go through all the plot categories this event qualifies for
    // (for example: charged currents, muon neutrinos, ...)
    for (icarus::trigger::PlotSandbox const& plotSet: plotSets) {
      
      //
      // overview plots from different thresholds
      //
      
      HistGetter const get { plotSet };
      
      // simple efficiency
      get.Eff("Eff"s).Fill(fired, minCount);
      
      // trigger time (if any)
      if (fired) {
        get.Hist2D("TriggerTick"s).Fill(minCount, triggerInfo.atTick().value());
      }
      
      //
      // plots depending on the trigger response
      // (but not caring of the trigger definition details)
      //
      
      // efficiency plots
      // (including event plots in the triggered or non-triggered category)
      helper().fillAllEfficiencyPlots
        (eventInfo, triggerInfo, plotSet.demandSandbox(minCountStr));
      
      //
      // add here further trigger-specific plots
      //
      
    } // for all qualifying plot categories
    
  } // for all requirements
  
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
    
    //
    // trigger-definition specific plots
    //
    HistGetter const get(plotSet);
    
    // number of primitives
    get.Hist("NPrimitives"s).Fill(maxPrimitives.second);
    
  } // for

} // icarus::trigger::MajorityTriggerEfficiencyPlots::plotResponses()


//------------------------------------------------------------------------------
auto icarus::trigger::MajorityTriggerEfficiencyPlots::combineTriggerPrimitives(
  TriggerGatesPerCryostat_t const& cryoGates,
  icarus::trigger::ADCCounts_t const threshold
) const -> TriggerGateData_t {

  //
  // simple count
  //
  std::vector<TriggerGateData_t> cryoCombinedGate;
  cryoCombinedGate.reserve(cryoGates.size());

  for (auto const& [ iCryo, gates ]: util::enumerate(cryoGates)) {
    geo::CryostatID const cryoID(iCryo);
    
    mf::LogTrace(helper().logCategory())
      << "Simulating trigger response with ADC threshold " << threshold
      << " for " << cryoID << " (" << gates.size() << " primitives)";

    if (gates.empty()) { // this is unexpected...
      mf::LogWarning(helper().logCategory())
        << "No trigger primitive found for threshold " << threshold
        << " in " << cryoID;
      return {};
    } // if no gates

    cryoCombinedGate.push_back(icarus::trigger::sumGates(gates));
  } // for

  //
  // largest number of trigger primitives at any time for any cryostat
  //
  return icarus::trigger::maxGates(cryoCombinedGate);
  
} // icarus::trigger::MajorityTriggerEfficiencyPlots::combineTriggerPrimitives()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::MajorityTriggerEfficiencyPlots)


//------------------------------------------------------------------------------
