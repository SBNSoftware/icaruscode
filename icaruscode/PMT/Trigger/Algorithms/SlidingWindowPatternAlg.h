/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.h
 * @brief  Applies sliding window trigger patterns.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWPATTERNALG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWPATTERNALG_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/WindowChannelMap.h"
#include "icaruscode/PMT/Trigger/Algorithms/WindowPattern.h"
#include "icaruscode/PMT/Trigger/Algorithms/ApplyBeamGate.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h"
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h" // sbn::OpDetWaveformMeta

// framework libraries
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <vector>
#include <optional>
#include <string>
#include <limits> // std::numeric_limits<>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace icarus::trigger { class SlidingWindowPatternAlg; }
/**
 * @brief Applies sliding window patterns to discriminated waveforms.
 * 
 * This algorithm takes as input one trigger gate for each window.
 * Each window is identified by an index, and the input gates are assigned one
 * per window, in order.
 * 
 * The window mapping that defines the topology of the windows is passed to the
 * algorithm.
 * 
 * If a beam gate is present, it is applied to all input before simulating the
 * pattern. Otherwise, the full input is used.
 * 
 * For the definition of the windows, see `icarus::trigger::WindowChannelMap`.
 * 
 */
class icarus::trigger::SlidingWindowPatternAlg
  : public icarus::ns::util::mfLoggingClass
{
  
    public:
  
  /// Record of the trigger response.
  using TriggerInfo_t = icarus::trigger::details::TriggerInfo_t;
  
  /// Type of trigger gate provided as input.
  using InputTriggerGate_t
    = icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>;
  
  /// A list of trigger gates from input.
  using TriggerGates_t = std::vector<InputTriggerGate_t>;

  /// Type of gate data without channel information (gate levels only).
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGateData_t;
  
  /// Type holding information about composition and topology of all windows.
  using WindowTopology_t = icarus::trigger::WindowChannelMap;
  
  /// Type representing the requirement pattern for a window.
  using WindowPattern_t = icarus::trigger::WindowPattern;
  
  /// Additional information on the trigger.
  struct MoreInfo_t {
    
    /// Index of the window which led the trigger.
    std::size_t windowIndex = std::numeric_limits<std::size_t>::max();
    
  }; // struct MoreInfo_t
  
  /// Complete information from this algorithm, standard + non-standard (extra).
  struct AllTriggerInfo_t {
    TriggerInfo_t info; ///< Standard trigger information.
    MoreInfo_t extra; ///< Extra trigger information.
    /// Returns whether this trigger fired.
    operator bool() const { return info.fired(); }
    /// Returns whether this trigger did not fire.
    bool operator!() const { return !info.fired(); }
  }; // AllTriggerInfo_t
  
  /**
   * @brief Constructor: configures window topology and times.
   * @param windowTopology full composition and topology description of windows
   * @param windowPattern the pattern that this algorithm applies
   * @param beamGate object applying the beam gate to trigger gates
   * @param logCategory category tag for algorithm messages on screen
   * 
   * Window topology can be computed with `icarus::trigger::WindowTopologyAlg`
   * starting from a (complete) set of trigger gates.
   * 
   * A beam gate application class can be constructed via
   * `icarus::trigger::makeApplyBeamGate()` helper function.
   * 
   */
  SlidingWindowPatternAlg(
    WindowTopology_t windowTopology,
    WindowPattern_t windowPattern,
    icarus::trigger::ApplyBeamGateClass beamGate,
    std::string const& logCategory = "SlidingWindowPatternAlg"
    );
  
  /**
   * @brief Constructor: configures window topology and times.
   * @param windowTopology full composition and topology description of windows
   * @param windowPattern the pattern that this algorithm applies
   * @param logCategory category tag for algorithm messages on screen
   * 
   * Window topology can be computed with `icarus::trigger::WindowTopologyAlg`
   * starting from a (complete) set of trigger gates.
   * 
   * The beam gate is set to be empty; it must be set before simulating the
   * response (`simulateResponse()`) via `setBeamGate()`.
   * 
   */
  SlidingWindowPatternAlg(
    WindowTopology_t windowTopology,
    WindowPattern_t windowPattern,
    std::string const& logCategory = "SlidingWindowPatternAlg"
    );
  
  /**
   * @brief Returns the trigger response from the specified `gates`.
   * @param gates the trigger gates to be used as input, one per window
   * @return the response to the configured pattern
   * @see `applyWindowPattern(WindowChannelMap::WindowInfo const&, WindowPattern_t const&, TriggerGates_t const&)`
   * 
   * The return value comprises a "standard" `TriggerInfo_t` object, which
   * stores whether the pattern fired and when, and some `extra` information,
   * including which window made the trigger fire first (in case of ties,
   * the window with the smallest index is reported).
   * 
   * The list `gates` must contain a single trigger gate for each window, with
   * the gate index in the list matching the window index in the configured
   * window topology.
   * 
   * See the static version of `applyWindowPattern()` for more details.
   * 
   */
  AllTriggerInfo_t simulateResponse(TriggerGates_t const& gates) const;
  

  /// Returns a new collection of gates, set each in coincidence with beam gate.
  TriggerGates_t applyBeamGate(TriggerGates_t const& gates) const;
  
  /// Returns whether a beam gate is being applied.
  bool hasBeamGate() const;
  
  /// Changes the beam gate to the specified value.
  void setBeamGate(icarus::trigger::ApplyBeamGateClass beamGate);
  
  /// Do not apply any beam gate.
  void clearBeamGate();
  
  
  /**
   * @brief Returns the trigger response for the specified window pattern.
   * @param windowInfo the topology of the windows
   * @param pattern the trigger requirement pattern
   * @param gates trigger gates, one per window
   * @return a `TriggerInfo_t` record with the response of the pattern
   * 
   * The input `gates` represent all the trigger gates relevant to the trigger
   * determination, in a list by window index: the first gate belongs to the
   * window `0`, the second one to the window `1` and so on.
   * The topology `windowInfo` identifies a window and its topology, i.e. its
   * neighbours, by window index. It is required and assumed that the values of
   * the window indices match the ones in the `gates` list.
   * This function applies the requirements in the specified `pattern` to the
   * window identified by `windowInfo`, using the actual trigger gates from
   * the `gates` list as needed.
   * 
   * The return value includes the first tick at which the requirements have
   * been all met.
   * Note that if a beam gate needs to be applied, it should be applies to
   * all `gates` before calling this function (see e.g. `applyBeamGate()`).
   * 
   */
  TriggerInfo_t applyWindowPattern(
    WindowTopology_t::WindowInfo_t const& windowInfo,
    WindowPattern_t const& pattern,
    TriggerGates_t const& gates
    ) const;
  
  
    private:
  
  /// Data structure to communicate internally a trigger response.
  struct WindowTriggerInfo_t {
    
    std::size_t windowIndex = std::numeric_limits<std::size_t>::max();
    TriggerInfo_t info;
    
    bool fired() const { return info.fired(); }
    
    operator bool() const { return bool(info); }
    bool operator! () const { return !info; }
    
    void emplace(std::size_t index, TriggerInfo_t info)
      { windowIndex = index; this->info = std::move(info); }
    
  }; // WindowTriggerInfo_t
  
  
  /// Definition of the neighborhood of each window in terms of window indices.
  WindowTopology_t const fWindowTopology;
  
  /// Requirement pattern to be applied to each window.
  WindowPattern_t const fWindowPattern;
  
  /// Time interval when to evaluate the trigger.
  std::optional<icarus::trigger::ApplyBeamGateClass const> fBeamGate;
  
  
  /**
   * @brief Returns the trigger response for the specified window pattern.
   * @param pattern the trigger requirement pattern
   * @param iWindow index of the window to be tested
   * @param gates trigger gates, one per window
   * @return a `TriggerInfo_t` record with the response of the pattern
   * @see `applyWindowPattern(WindowChannelMap::WindowInfo const&, WindowPattern_t const&, TriggerGates_t const&)`
   * 
   * Applies the specified `pattern` to the window with index `iWindow` of the
   * configured window topology.
   * 
   * See the static version of `applyWindowPattern()` for details.
   */
  TriggerInfo_t applyWindowPattern(
    WindowPattern_t const& pattern, std::size_t iWindow,
    TriggerGates_t const& gates
    ) const;
  
  /**
   * @brief Checks `gates` are compatible with the current window configuration.
   * @param gates the combined sliding window trigger gates, per cryostat
   * @throw cet::exception (category: `SlidingWindowTriggerEfficiencyPlots`)
   *        or derived, if an incompatibility is found
   * 
   * The method verifies that the current channel mapping is compatible with the
   * gates.
   * 
   * This currently means that the `gates` are in the expected order and have
   * the expected channel content.
   */
  void verifyInputTopology(TriggerGates_t const& gates) const;
  
}; // class icarus::trigger::SlidingWindowPatternAlg


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWPATTERNALG_H
