/**
 * @file   SlidingWindowTriggerEfficiencyPlots_module.cc
 * @brief  Plots of efficiency for triggers based on PMT channel sliding window.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 9, 2020
 * @see    icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h"
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h"
#include "icaruscode/PMT/Trigger/Utilities/ROOTutils.h" // util::ROOT
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
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// ROOT libraries
#include "TTree.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <algorithm> // std::sort(), std::transform(), ...
#include <vector>
#include <array>
#include <optional>
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
namespace icarus::trigger::details {
  
  class WindowChannelMap;
  
  // [cryostat][window index in cryostat] => list of channels in window
  using WindowChannelsPerCryostat_t
    = std::vector<std::vector<std::vector<raw::Channel_t>>>;
  
  WindowChannelMap makeWindowChannelMap(
    WindowChannelsPerCryostat_t const& allWindowChannels,
    geo::GeometryCore const& geom
    );
  
} // namespace icarus::trigger::details
class icarus::trigger::details::WindowChannelMap {
  
    public:
  
  using WindowIndex_t = std::size_t; ///< Type of window index
  
  /// Special index denoting an invalid window.
  static const WindowIndex_t InvalidWindowIndex
    = std::numeric_limits<WindowIndex_t>::max();
  
  struct WindowInfo {
    
    WindowIndex_t index; ///< Index of the window this information is about.
    
    geo::Point_t center; ///< Center of the window.
    
    geo::CryostatID cryoid; ///< Which cryostat the channels are in.
    
    /// Optical detector channels covered by this window.
    std::vector<raw::Channel_t> channels;
    
    /// Index of the window opposite to this one.
    WindowIndex_t opposite { InvalidWindowIndex };
    
    /// Index of the window upstream of this one.
    WindowIndex_t upstream { InvalidWindowIndex };
    
    /// Index of the window downstream of this one.
    WindowIndex_t downstream { InvalidWindowIndex };
    
    /// Returns whether the window is in a single, known cryostat.
    bool hasCryostat() const { return cryoid.isValid; }
    
    /// Returns whether the main window has another upstream of it.
    bool hasUpstreamWindow() const { return isValidWindow(upstream); }
    
    /// Returns whether the main window has another downstream of it.
    bool hasDownstreamWindow() const { return isValidWindow(downstream); }
    
    /// Returns whether the main window has another downstream of it.
    bool hasOppositeWindow() const { return isValidWindow(opposite); }
    
    /// Prints the information content (single line).
    template <typename Stream>
    void dump(Stream&& out, std::string const& indent = "") const;
    
  }; // struct WindowInfo
  
  
  /// Construction: moves the `windows` information into the map.
  explicit WindowChannelMap(std::vector<WindowInfo>&& windows)
    : fWindows(std::move(windows)) {}
  
  // --- BEGIN Access to window information ------------------------------------
  /// @name Access to window information.
  /// @{
  
  /// Number of sliding windows.
  std::size_t nWindows() const { return fWindows.size(); }
  
  /// Returns whether a window with the specified index is present.
  bool hasWindow(WindowIndex_t index) const { return index < nWindows(); }
  
  /// Returns the information for the window with specified `index` (unchecked).
  WindowInfo const& info(WindowIndex_t index) const { return fWindows[index]; }
  
  /// @}
  // --- END Access to window information --------------------------------------

  // @{
  /// Prints the content of the full mapping.
  template <typename Stream>
  void dump
    (Stream&& out, std::string const& indent, std::string const& firstIndent)
    const;
  template <typename Stream>
  void dump(Stream&& out, std::string const& indent = "") const
    { dump(std::forward<Stream>(out), indent, indent); }
  // @}
  

  
  /// @name Iterable standard interface
  /// @{
  using value_type = WindowInfo;
  using reference_type = WindowInfo const&;
  using pointer_type = WindowInfo const*;
  
  std::size_t size() const { return nWindows(); }
  bool empty() const { return fWindows.empty(); }
  auto begin() const { return fWindows.begin(); }
  auto end() const { return fWindows.end(); }
  auto cbegin() const { return begin(); }
  auto cend() const { return end(); }
  /// @}
  
  
  /// Returns whether the specified index is not the invalid window index.
  static bool isValidWindow(WindowIndex_t index)
    { return index != InvalidWindowIndex; }
  
  
    private:
  std::vector<WindowInfo> fWindows; /// Information for each window.
  
  friend WindowChannelMap makeWindowChannelMap
    (WindowChannelsPerCryostat_t const&, geo::GeometryCore const&);
  
}; // icarus::trigger::details::WindowChannelMap


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
 * Note that the logical waveforms from the sliding windows are expected to
 * be provided as input.
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
 * TODO
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::SlidingWindowTriggerEfficiencyPlots` and its assumptions.
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
 * `plotResponse()`.
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
 * A list of additional plots follows for each plot type.
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
    
    struct WindowPattern {
      
      fhicl::Atom<unsigned int> inMainWindow {
        Name("inMainWindow"),
        Comment("minimum fired primitives in the main sliding window")
        };
      fhicl::Atom<unsigned int> inUpstreamWindow {
        Name("inUpstreamWindow"),
        Comment(
          "minimum fired primitives in the sliding window upstream of main one"
          ),
        0U // default
        };
      fhicl::Atom<unsigned int> inDownstreamWindow {
        Name("inDownstreamWindow"),
        Comment(
         "minimum fired primitives in the sliding window downstream of main one"
         ),
        0U // default
        };
      fhicl::Atom<unsigned int> inOppositeWindow {
        Name("inOppositeWindow"),
        Comment(
          "minimum fired primitives in the sliding window opposite of main one"
          ),
        0U // default
        };
      
      fhicl::Atom<bool> requireUpstreamWindow {
        Name("requireUpstreamWindow"),
        Comment("an upstream window must be present (no border main window)"),
        false
        };
      fhicl::Atom<bool> requireDownstreamWindow {
        Name("requireDownstreamWindow"),
        Comment("a downstream window must be present (no border main window)"),
        false
        };
      
    }; // WindowPattern
    
    fhicl::Sequence<fhicl::Table<WindowPattern>> Patterns {
      Name("Patterns"),
      Comment("sliding window pattern requirements")
      };
    
  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit SlidingWindowTriggerEfficiencyPlots(Parameters const& config);

  // Plugins should not be copied or assigned.
  SlidingWindowTriggerEfficiencyPlots(SlidingWindowTriggerEfficiencyPlots const&) = delete;
  SlidingWindowTriggerEfficiencyPlots(SlidingWindowTriggerEfficiencyPlots&&) = delete;
  SlidingWindowTriggerEfficiencyPlots& operator=(SlidingWindowTriggerEfficiencyPlots const&) = delete;
  SlidingWindowTriggerEfficiencyPlots& operator=(SlidingWindowTriggerEfficiencyPlots&&) = delete;

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
  
  /// Specification of the requirement of sliding window firing pattern.
  struct WindowPattern {
    
    /// @name Minimum required number of open trigger primitives per window.
    /// @{
    unsigned int minInMainWindow;
    unsigned int minInUpstreamWindow;
    unsigned int minInDownstreamWindow;
    unsigned int minInOppositeWindow;
    /// @}
    
    /// Whether a window location with no upstream window should be discarded.
    bool requireUpstreamWindow = false;
    
    /// Whether a window location with no downstream window should be discarded.
    bool requireDownstreamWindow = false;
    
    /// Returns a tag summarizing the pattern.
    std::string tag() const;
    
    /// Returns a description of the pattern.
    std::string description() const;
    
  }; // WindowPattern
  
  
  friend std::string to_string
    (SlidingWindowTriggerEfficiencyPlots::WindowPattern const& pattern);
  
  /// List of configured patterns.
  using WindowPatterns_t = std::vector<WindowPattern>;
  
  using TriggerInfo_t = details::TriggerInfo_t; // type alias
  
  
  /// Data structure to communicate internally a trigger response.
  struct WindowTriggerInfo {
    
    std::size_t windowIndex = std::numeric_limits<std::size_t>::max();
    TriggerInfo_t info;
    
    bool fired() const { return info.fired(); }
    
    operator bool() const { return bool(info); }
    bool operator! () const { return !info; }
    
    void emplace(std::size_t index, TriggerInfo_t info)
      { windowIndex = index; this->info = std::move(info); }
    
  }; // WindowTriggerInfo
  

  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Configured sliding window requirement patterns.
  WindowPatterns_t const fPatterns;
  
  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Internal variables ----------------------------------------------
  
  /// Mapping of each sliding window with location and topological information.
  // mutable = not thread-safe; optional to allow delayed construction
  mutable std::optional<details::WindowChannelMap> fWindowMap;
  
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
    ) const override;
    
  // --- END Derived class methods ---------------------------------------------

  /**
   * @brief Returns the trigger time for the specified sliding window, if any.
   * @param pattern the pattern used to decide whether the trigger fires
   * @param iWindow the index of the main sliding window in the pattern
   * @param gates all the sliding window gates
   * @return the trigger information for the specified sliding window
   * 
   * This method applies the specified pattern using `iWindow` as the main
   * window of the pattern, and assigning the upstream, downstream and opposite
   * windows with the current topological map.
   * 
   * The return contains the earliest optical time tick at which the pattern
   * requirements are all satisfied. If that never happens, the returned value
   * has instead `fired()` returning `false`.
   */
  TriggerInfo_t applyWindowPattern(
    WindowPattern const& pattern,
    std::size_t iWindow,
    TriggerGates_t const& gates
    ) const;

  
  /**
   * @brief Fills plots with the specified trigger response.
   * @param iThr index of PMT threshold (used in tree output)
   * @param threshold PMT threshold in ADC counts (for printing)
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
    std::size_t iThr, ADCCounts_t const threshold,
    std::size_t iPattern, WindowPattern const& pattern,
    PlotSandboxRefs_t const& plotSets,
    EventInfo_t const& eventInfo,
    WindowTriggerInfo const& triggerInfo
    ) const;

  /// Fills all event plots with data from `eventInfo` as in `fillEventPlots()`.
  void fillAllEventPlots
    (PlotSandboxRefs_t const& plotSets, EventInfo_t const& eventInfo) const;

  /**
   * @brief Builds the channel maps from the specified gates.
   * @param gates the combined sliding window trigger gates, per cryostat
   * @param force (default: `false`) if `true`, any existing map is discarded
   * @return `true` if initialization happened, `false` if it was already there
   * 
   * The `gates` are parsed to determine the position and composition of each
   * sliding window and to order them and determine where their neighbours are.
   * 
   * The internal map is initialized. If that map was already initialized,
   * no action happens if `force` is `false`, while the map is recomputed
   * from scratch if `force` is `true`.
   * 
   * NOTE Multithreading: this needs locking.
   */
  bool initializeTopologicalMaps
    (TriggerGatesPerCryostat_t const& gates, bool force = false) const;
  
  /**
   * @brief Verifies that the current channel maps are compatible with `gates`.
   * @param gates the combined sliding window trigger gates, per cryostat
   * @throw IncompatibleMap (category: `SlidingWindowTriggerEfficiencyPlots`)
   *        or derived, if an incompatibility is found
   * 
   * The method verifies that the current channel mapping is compatible with the
   * gates.
   * 
   * This currently means that the `gates` are in the expected order and have
   * the expected channel content.
   */
  void verifyTopologicalMap(TriggerGatesPerCryostat_t const& gates) const;
  
  
  /// Builds a `WindowPattern` object from its FHiCL configuration.
  static WindowPattern makeWindowPattern(Config::WindowPattern const& config);
  
  /// Builds a sequence of `WindowPattern` objects from FHiCL configuration.
  static WindowPatterns_t makeWindowPatterns
    (std::vector<Config::WindowPattern> const& config);
  
}; // icarus::trigger::SlidingWindowTriggerEfficiencyPlots


namespace icarus::trigger {
  
  /// Converts a `pattern` into a compact string.
  std::string to_string
    (SlidingWindowTriggerEfficiencyPlots::WindowPattern const& pattern);
  
} // namespace icarus::trigger



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
//--- icarus::trigger::details::WindowChannelMap
//------------------------------------------------------------------------------
template <typename Stream>
void icarus::trigger::details::WindowChannelMap::WindowInfo::dump
  (Stream&& out, std::string const& indent /* = "" */) const
{
  out << indent << "window #" << index << " at " << center << " cm";
  if (hasCryostat()) out << " in " << cryoid;
  else               out << " (cryostat undefined)";
  out << " includes " << channels.size() << " channels";
  if (!channels.empty()) {
    auto iChannel = channels.begin();
    auto const cend = channels.end();
    out << " (" << *iChannel;
    while (++iChannel != cend) out << ", " << *iChannel;
    out << ")";
  } // if channels
  if (hasOppositeWindow()) out << " opposite to [#" << opposite << "]";
  if (hasUpstreamWindow()) out << " downstream of [#" << upstream << "]";
  if (hasDownstreamWindow()) out << " upstream of [#" << downstream << "]";
} // icarus::trigger::details::WindowChannelMap::WindowInfo::dump()


//------------------------------------------------------------------------------
template <typename Stream>
void icarus::trigger::details::WindowChannelMap::dump
  (Stream&& out, std::string const& indent, std::string const& firstIndent)
  const
{
  out << firstIndent << "Map has " << nWindows() << " windows:";
  for (WindowInfo const& info: fWindows) {
    out << "\n  "; // additional indentation
    info.dump(std::forward<Stream>(out), indent);
  } // for
  out << '\n';
} // icarus::trigger::details::WindowChannelMap::dump()


//------------------------------------------------------------------------------
auto icarus::trigger::details::makeWindowChannelMap(
  WindowChannelsPerCryostat_t const& allWindowChannels,
  geo::GeometryCore const& geom
  ) -> WindowChannelMap
{
  /*
   * 1.     For each cryostat:
   * 1.1.     fill the window information with local information
   * 1.2.     sort the windows in drift plane (first cryostat TPC as reference)
   * 1.3.     split the windows per plane
   * 1.4.     for each plane:
   * 1.4.1.     sort windows on beam direction (TPC width direction)
   * 1.4.2.     fill the neighbour information on each window
   * 1.5.     add the information to the complete map
   * 2.     return a map with the collected information
   */
  std::vector<WindowChannelMap::WindowInfo> windows;
  
  //
  // 1.     For each cryostat:
  //
  std::size_t iWindow = 0U;
  for (auto const& [ iCryo, cryoWindows ]: util::enumerate(allWindowChannels)) {
    
    // window number required to be even (each window has an unique opposite)
    assert(cryoWindows.size() % 2 == 0);
    
    geo::CryostatID const cryoid(iCryo);
    assert(geom.HasCryostat(cryoid));
    
    // use the first TPC of the cryostat as reference for directions
    geo::CryostatGeo const& cryo = geom.Cryostat(cryoid);
    assert(cryo.NTPC() > 0U);
    geo::TPCGeo const& refTPC = cryo.TPC(0U);
    
    
    //
    // 1.1. fill the window information with local information
    //
    using WindowInfoPtrs_t = std::vector<WindowChannelMap::WindowInfo*>;
    
    windows.reserve(windows.size() + cryoWindows.size());
    
    WindowInfoPtrs_t cryoWindowInfo;
    cryoWindowInfo.reserve(cryoWindows.size());
    
    for (auto const& channels: cryoWindows) {
      
      WindowChannelMap::WindowInfo wInfo;
      
      wInfo.index = iWindow++;
      wInfo.channels = channels;
      std::sort(wInfo.channels.begin(), wInfo.channels.end());
      wInfo.cryoid = channels.empty()
        ? geo::CryostatID{}
        : geom.OpDetGeoFromOpChannel(channels.front()).ID().asCryostatID()
        ;
      
      geo::vect::MiddlePointAccumulator middlePoint;
      for (raw::Channel_t channel: channels) {
        // documentation of OpDetGeoFromOpChannel() does not say what on error..
        geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpChannel(channel);
        middlePoint.add(opDet.GetCenter());
        if (opDet.ID() != wInfo.cryoid) wInfo.cryoid = geo::CryostatID{};
      } // for channel
      wInfo.center = middlePoint.middlePoint();
      
      windows.push_back(std::move(wInfo));
      cryoWindowInfo.push_back(&windows.back());
      
    } // for windows
    
    //
    // 1.2. sort the windows in drift plane (first cryostat TPC as reference)
    //
    auto const normalProjection = [&refTPC](auto const* info)
     { return refTPC.DistanceFromReferencePlane(info->center); };
    WindowInfoPtrs_t const windowsByNormal
      = util::sortCollBy(cryoWindowInfo, normalProjection);
    
    //
    // 1.3. split the windows per plane
    //
    // split the list in two; there is a good deal of faith here
    auto const beamCoordinate
      = [&refPlane=refTPC.ReferencePlane()](auto const* info)
        { return refPlane.PointWidthDepthProjection(info->center).X(); }
      ;

    //
    // 1.4. for each plane:
    //
    // 1.4.1. sort windows on beam direction (TPC width direction)
    //
    auto const iMiddleWindow
      = std::next(windowsByNormal.cbegin(), windowsByNormal.size() / 2U);
    std::array<WindowInfoPtrs_t, 2U> const windowsByPlane = {
      util::sortBy(windowsByNormal.cbegin(), iMiddleWindow, beamCoordinate),
      util::sortBy(iMiddleWindow, windowsByNormal.cend(), beamCoordinate)
    };
    
    for (auto const& [ iPlane, planeWindows ]: util::enumerate(windowsByPlane))
    {
      //
      // 1.4.2.     fill the neighbour information on each window
      //
      auto const& otherPlaneWindows
        = windowsByPlane.at(windowsByPlane.size() - 1U - iPlane);
      std::size_t const iLastPlaneWindow = planeWindows.size() - 1U;
      for (auto [ iPlaneWindow, windowInfo ]: util::enumerate(planeWindows)) {
        
        // assumes all topology information is InvalidWindowIndex by default
        windowInfo->opposite = otherPlaneWindows[iPlaneWindow]->index;
        
        if (iPlaneWindow > 0U)
          windowInfo->upstream = planeWindows[iPlaneWindow - 1U]->index;
        
        if (iPlaneWindow < iLastPlaneWindow)
          windowInfo->downstream = planeWindows[iPlaneWindow + 1U]->index;
        
      } // for window in plane
    } // for planes
    
    //
    // 1.5.     add the information to the complete map
    //
    // having acted on pointers, all information is already updated;
    // and the window information was created in the right position from start
    
  } // for cryostats
  
  //
  // 2. return a map with the collected information
  //
  return WindowChannelMap{ std::move(windows) };
  
} // icarus::trigger::details::makeWindowChannelMap()


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

  for (auto [ iThr, threshold]: util::enumerate(thresholds)) {
    std::string const thrStr = util::to_string(raw::ADC_Count_t(threshold));

    for (auto [ iSetting, setting ]: util::enumerate(settings)) {

      std::string const branchName
        = "RespT" + thrStr + "S" + util::to_string(setting);

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
std::string
icarus::trigger::SlidingWindowTriggerEfficiencyPlots::WindowPattern::tag
  () const
{
  using namespace std::string_literals;
  
  std::string s;
  
  s += "M"s + std::to_string(minInMainWindow);
  
  if (minInOppositeWindow > 0U)
    s += "O"s + std::to_string(minInOppositeWindow);
  
  if ((minInDownstreamWindow > 0U) || requireDownstreamWindow) {
    s += "D"s + std::to_string(minInDownstreamWindow);
    if (requireDownstreamWindow) s+= "req"s;
  } // if downstream
  
  if ((minInUpstreamWindow > 0U) || requireUpstreamWindow) {
    s += "U"s + std::to_string(minInUpstreamWindow);
    if (requireUpstreamWindow) s+= "req"s;
  } // if upstream
  
  return s;
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::WindowPattern::description()


//------------------------------------------------------------------------------
std::string
icarus::trigger::SlidingWindowTriggerEfficiencyPlots::WindowPattern::description
  () const
{
  using namespace std::string_literals;
  
  std::string s = "required:";
  
  s += " "s + std::to_string(minInMainWindow);
  
  if (minInOppositeWindow > 0U)
    s += " + " + std::to_string(minInOppositeWindow) + " (opposite)"s;
  
  if ((minInDownstreamWindow > 0U) || requireDownstreamWindow) {
    s += " + " + std::to_string(minInDownstreamWindow) + " (downstream";
    if (requireDownstreamWindow) s+= ", mandatory)"s;
    s += ")"s;
  } // if downstream
  
  if ((minInUpstreamWindow > 0U) || requireUpstreamWindow) {
    s += " + " + std::to_string(minInUpstreamWindow) + " (upstream";
    if (requireUpstreamWindow) s+= ", mandatory)"s;
    s += ")"s;
  } // if upstream
  
  return s;
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::WindowPattern::description()


//------------------------------------------------------------------------------
std::string icarus::trigger::to_string
  (icarus::trigger::SlidingWindowTriggerEfficiencyPlots::WindowPattern const& pattern)
  { return pattern.tag(); }


//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  std::string to_string
    (SlidingWindowTriggerEfficiencyPlots::WindowPattern const& pattern);
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
icarus::trigger::SlidingWindowTriggerEfficiencyPlots::SlidingWindowTriggerEfficiencyPlots
  (Parameters const& config)
  : art::EDAnalyzer           (config)
  , TriggerEfficiencyPlotsBase(config(), consumesCollector())
  // configuration
  , fPatterns(makeWindowPatterns(config().Patterns()))
{
  
  if (fPatterns.empty()) {
    throw art::Exception(art::errors::Configuration)
      << "At least one 'MinimumPrimitives' requirement... required.";
  }
  
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
  util::ROOT::applyAxisLabels
    (const_cast<TH1*>(Eff->GetTotalHistogram())->GetXaxis(), patternLabels);
  
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::initializePlotSet()


// -----------------------------------------------------------------------------
bool icarus::trigger::SlidingWindowTriggerEfficiencyPlots::initializeTopologicalMaps
  (TriggerGatesPerCryostat_t const& gates, bool force /* = false */) const
{
  if (!force && fWindowMap) return false;
  
  // extract the channel numbers from the windows
  // [cryostat][window index in cryostat] => list of channels in window
  details::WindowChannelsPerCryostat_t channels;
  channels.reserve(gates.size());
  
  for (TriggerGates_t const& cryoGates: gates) {
    // all... in this cryostat:
    details::WindowChannelsPerCryostat_t::value_type allWindowChannels;
    allWindowChannels.reserve(cryoGates.size());
    
    for (InputTriggerGate_t const& gate: cryoGates) {
      auto const& channels = gate.channels();
      allWindowChannels.emplace_back(channels.begin(), channels.end());
    } // for windows
    
    channels.push_back(std::move(allWindowChannels));
  } // for cryostats
  
  fWindowMap.emplace
    (details::makeWindowChannelMap(channels, helper().geometry()));
  
  {
    mf::LogTrace log { helper().logCategory() };
    log << "Window map: ";
    fWindowMap->dump(log);
  }
  
  return true;
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::initializeTopologicalMaps()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::verifyTopologicalMap
  (TriggerGatesPerCryostat_t const& gates) const
{
  /*
   *    * @brief Verifies that the current channel maps are compatible with `gates`.
   * @param gates the combined sliding window trigger gates, per cryostat
   * @throw IncompatibleMap (category: `SlidingWindowTriggerEfficiencyPlots`)
   *        or derived, if an incompatibility is found
   * 
   * The method verifies that the current channel mapping is compatible with the
   * gates.
   * 
   * This currently means that the `gates` are in the expected order and have
   * the expected channel content.
   */

  std::size_t iWindow = 0U;
  std::string errorMsg; // if this stays `empty()` there is no error
  for (auto const& cryoGates: gates) {
    for (auto const& gate: cryoGates) {
      
      std::string windowError; // if this stays `empty()` there is no error
      
      details::WindowChannelMap::WindowInfo const& windowInfo
        = fWindowMap->info(iWindow++);
      
      auto const channelInWindow
        = [begin=windowInfo.channels.cbegin(),end=windowInfo.channels.cend()]
        (raw::Channel_t channel)
        { return std::binary_search(begin, end, channel); }
        ;
      
      for (raw::Channel_t const channel: gate.channels()) {
        if (channelInWindow(channel)) continue;
        if (windowError.empty()) {
          windowError =
            "channels not in window #" + std::to_string(windowInfo.index)
            + ":";
        } // if first error
        windowError += " " + std::to_string(channel);
      } // for all channels in gate
      
      if (!windowError.empty()) errorMsg += windowError + '\n';
    } // for gates in cryostat
  } // for cryostats
  
  if (errorMsg.empty()) return;
  
  // put together the exception message and throw it.
  cet::exception e("SlidingWindowTriggerEfficiencyPlots");
  e << "Some channels from trigger gates do not match the previous window allocation:\n"
    << errorMsg
    << "\n"
    << "Window allocation: ";
  fWindowMap->dump(e, "  ");
  throw e;
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::verifyTopologicalMap()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::simulateAndPlot(
  std::size_t const thresholdIndex,
  TriggerGatesPerCryostat_t const& gates,
  EventInfo_t const& eventInfo,
  detinfo::DetectorClocksData const& clockData,
  PlotSandboxRefs_t const& selectedPlots
) const {
  
  auto const threshold = helper().ADCthreshold(thresholdIndex);
  
  /*
   * 0.   initialize or verify the topology of the input
   * 1.   apply the beam gate to each input gate
   * 2.   for each pattern:
   * 2.1.   for each main window, apply the pattern
   * 2.2.   pick the main window with the earliest successful response, if any;
   *        that defines location and time of the trigger
   * 2.3.   plot the trigger outcome
   * 3.   fill all trigger-independent plots
   * 
   */
  
  //
  // 0. initialize or verify the topology of the input
  //
  if (!initializeTopologicalMaps(gates)) verifyTopologicalMap(gates);
  
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
    
    for (auto const& [ iWindow, gate ]: util::enumerate(inBeamGates)) {
      auto const maxTick = gate.findMaxOpen();
      if (maxTick == gate.MinTick) continue;
      log << "\n  window #" << iWindow << ": maximum "
        << gate.openingCount(maxTick) << " at tick " << maxTick;
    } // for
  }
  // --- END DEBUG -------------------------------------------------------------
  
  //
  // 2.   for each pattern:
  //
  std::size_t const nWindows = fWindowMap->nWindows();
  for (auto const& [ iPattern, pattern ]: util::enumerate(fPatterns)) {
    
    //
    // 2.1.   for each main window, apply the pattern
    //
    WindowTriggerInfo triggerInfo; // start empty
    for (std::size_t const iWindow: util::counter(nWindows)) {
      
      TriggerInfo_t const windowResponse
        = applyWindowPattern(pattern, iWindow, inBeamGates);
      
      if (!windowResponse) continue;
      
      mf::LogTrace(helper().logCategory())
        << "Pattern '" << pattern.tag() << "' on window #" << iWindow
        << " (threshold: " << threshold
        << ") fired at tick " << windowResponse.atTick() << " ("
        << detinfo::DetectorTimings(clockData).toElectronicsTime
          (detinfo::DetectorTimings::optical_tick{ windowResponse.atTick() })
        << ")"
        ;
      
      //
      // 2.2.   pick the main window with the earliest successful response, if any;
      //        that defines location and time of the trigger
      //
      if (!triggerInfo || triggerInfo.info.atTick() > windowResponse.atTick())
        triggerInfo.emplace(iWindow, windowResponse);
      
    } // main window choice
    
    //
    // 2.3.   plot the trigger outcome
    //
    plotResponse(
      thresholdIndex, threshold,
      iPattern, pattern,
      selectedPlots,
      eventInfo, triggerInfo
      );
    
  } // for window patterns
  
  //
  // 3. fill all trigger-independent plots (one copy per threshold... meh)
  //
  fillAllEventPlots(selectedPlots, eventInfo);
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::simulateAndPlot()


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
/**
  * @brief Returns the trigger time for the specified sliding window, if any.
  * @param pattern the pattern used to decide whether the trigger fires
  * @param iWindow the index of the main sliding window in the pattern
  * @param gates all the sliding window gates
  * @return the trigger information for the specified sliding window
  * 
  * This method applies the specified pattern using `iWindow` as the main
  * window of the pattern, and assigning the upstream, downstream and opposite
  * windows with the current topological map.
  * 
  * The return contains the earliest optical time tick at which the pattern
  * requirements are all satisfied. If that never happens, the returned value
  * has instead `fired()` returning `false`.
  */
auto icarus::trigger::SlidingWindowTriggerEfficiencyPlots::applyWindowPattern(
  WindowPattern const& pattern, std::size_t iWindow, TriggerGates_t const& gates
  ) const -> TriggerInfo_t
{
  
  /*
   * 1. check that the pattern can be applied; if not, return no trigger
   * 2. discriminate all the relevant gates against their required minimum count
   * 3. combine them in AND
   * 4. find the trigger time, fill the trigger information accordingly
   */
  assert(fWindowMap);
  details::WindowChannelMap::WindowInfo const& windowInfo
    = fWindowMap->info(iWindow);
  assert(windowInfo.index == iWindow);
  assert(windowInfo.hasOppositeWindow());
  
  TriggerInfo_t res; // no trigger by default
  assert(!res);

  //
  // 1. check that the pattern can be applied; if not, return no trigger
  //
  
  // check that the pattern centered into iWindow has all it needs:
  if (pattern.requireUpstreamWindow && !windowInfo.hasUpstreamWindow())
    return res;
  if (pattern.requireDownstreamWindow && !windowInfo.hasDownstreamWindow())
    return res;
  
  
  //
  // 2. discriminate all the relevant gates against their required minimum count
  // 3. combine them in AND
  //
  
  // main window
  TriggerGateData_t trigPrimitive
    = discriminate(gates[windowInfo.index], pattern.minInMainWindow);
  
  // add opposite window requirement (if any)
  if ((pattern.minInOppositeWindow > 0U) && windowInfo.hasOppositeWindow()) {
    trigPrimitive.Mul
      (discriminate(gates[windowInfo.opposite], pattern.minInOppositeWindow));
  } // if
  
  // add upstream window requirement (if any)
  if ((pattern.minInUpstreamWindow > 0U) && windowInfo.hasUpstreamWindow()) {
    trigPrimitive.Mul
      (discriminate(gates[windowInfo.upstream], pattern.minInUpstreamWindow));
  } // if
  
  // add downstream window requirement (if any)
  if ((pattern.minInDownstreamWindow > 0U) && windowInfo.hasDownstreamWindow())
  {
    trigPrimitive.Mul(
      discriminate(gates[windowInfo.downstream], pattern.minInDownstreamWindow)
      );
  } // if
  
  //
  // 4. find the trigger time, fill the trigger information accordingly
  //
  auto const trigTick = trigPrimitive.findOpen(); // first trigger
  if (trigTick != trigPrimitive.MaxTick) {
    res.emplace(detinfo::timescales::optical_tick{ trigTick });
    assert(res);
  }
  
  return res;
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::applyWindowPattern()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTriggerEfficiencyPlots::plotResponse(
  std::size_t iThr, icarus::trigger::ADCCounts_t const threshold,
  std::size_t iPattern, WindowPattern const& pattern,
  PlotSandboxRefs_t const& plotSets,
  EventInfo_t const& eventInfo,
  WindowTriggerInfo const& triggerInfo
) const {
  
  using namespace std::string_literals;
  
  bool const fired = triggerInfo.fired();
  
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
      (eventInfo, triggerInfo.info, plotSet.demandSandbox(patternTag));
    
    //
    // add here further trigger-specific plots
    //
    
  } // for all qualifying plot categories
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::plotResponse()


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowTriggerEfficiencyPlots::makeWindowPatterns
  (std::vector<Config::WindowPattern> const& config) -> WindowPatterns_t
{
  
  WindowPatterns_t patterns;
  patterns.reserve(config.size());
  for (Config::WindowPattern const& patternConfig: config) {
    
    WindowPattern const pattern {
      patternConfig.inMainWindow(),           // minInMainWindow
      patternConfig.inUpstreamWindow(),       // minInUpstreamWindow
      patternConfig.inDownstreamWindow(),     // minInDownstreamWindow
      patternConfig.inOppositeWindow(),       // minInOppositeWindow
      patternConfig.requireUpstreamWindow(),  // requireUpstreamWindow
      patternConfig.requireDownstreamWindow() // requireDownstreamWindow
    };
    patterns.push_back(std::move(pattern));
    
  } // for
  
  return patterns;
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::makeWindowPatterns()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::SlidingWindowTriggerEfficiencyPlots)


//------------------------------------------------------------------------------
