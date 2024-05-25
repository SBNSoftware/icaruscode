/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h
 * @brief  Helper class to store transient trigger result.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 15, 2020
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_TRIGGERINFO_T_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_TRIGGERINFO_T_H

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h" // gateIn()
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_tick
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h" // 

// C++ standard libraries
#include <vector>
#include <algorithm>
#include <optional>
#include <utility> // std::pair, std::tie()
#include <limits> // std::numeric_limits<>
#include <utility> // std::forward()
#include <type_traits> // std::decay_t
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus::trigger::details {
  struct TriggerInfo_t;
  template <typename Gate> class GateOpeningInfoExtractor;
} // namespace icarus::trigger::details

// -----------------------------------------------------------------------------
/**
 * @brief Helper data structure to store transient trigger result.
 * 
 * This record maintains a list of all openings (`add()`), and has a special one
 * (`main()`) which represents the "global" trigger.
 * 
 * Note that unless no opening is `add()`'ed, there is always a `main()`
 * trigger.
 * 
 * Each trigger is described by the information in a `OpeningInfo_t` record.
 * The records can be registered with `add()`; `replace()` changes the `main()`
 * trigger information unconditionally, while `addAndReplaceIfEarlier()`
 * replaces `main()` only if the argument proposes a trigger earlier than the
 * current `main()`.
 * 
 */
struct icarus::trigger::details::TriggerInfo_t {
  
  using optical_tick = detinfo::timescales::optical_tick; ///< Type alias.
  
  /// Type of gate opening level.
  using Opening_t
    = icarus::trigger::OpticalTriggerGateData_t::GateData_t::OpeningCount_t;
  
  using LocationID_t = std::size_t; ///< Type for ID of trigger location.
  
  struct OpeningInfo_t {
    
    /// ID for a trigger in unknown location.
    static constexpr LocationID_t UnknownLocation
      = std::numeric_limits<std::size_t>::max();
    
    /// Tick at which the trigger fired.
    optical_tick tick = std::numeric_limits<optical_tick>::min();
    
    /// Tick at which trigger requirements are not met any more.
    optical_tick endTick = std::numeric_limits<optical_tick>::min();
    
    Opening_t level = 0U; ///< Maximum level on the main trigger opening.
    
    /// Identified of the trigger location.
    LocationID_t locationID = UnknownLocation;
    
    
    /// Default constructor.
    OpeningInfo_t() = default;
    
    /// Constructor: specify `tick`, the rest is optional.
    OpeningInfo_t(
      optical_tick tick, optical_tick endTick,
      Opening_t level = 0U, std::size_t locationID = UnknownLocation
      )
      : tick(tick), endTick(endTick), level(level), locationID(locationID) {}
    
    /// Returns whether the location is set.
    bool hasLocation() const { return locationID != UnknownLocation; }
    
    /// Comparison: order from time.
    bool operator< (OpeningInfo_t const& other) const
      { return tick < other.tick; }
    
  }; // OpeningInfo_t
  
  
  // --- BEGIN -- Construction -------------------------------------------------
  
  TriggerInfo_t() = default; // no trigger
  TriggerInfo_t(OpeningInfo_t info) { replace(std::move(info)); }
  
  // --- END -- Construction ---------------------------------------------------
  
  
  // --- BEGIN -- Query whether the trigger fired ------------------------------
  /// @name Query whether the trigger fired.
  /// @{
  
  /// Returns whether the trigger fired.
  bool fired() const { return !fAll.empty(); }
  
  /// Returns whether there is trigger information.
  operator bool() const { return fired(); }
  
  /// Returns whether there is no trigger information.
  bool operator! () const { return !fired(); }
  
  /// @}
  // --- END -- Query whether the trigger fired --------------------------------
  
  
  // --- BEGIN -- Modify trigger information -----------------------------------
  /// @name Modify trigger information
  /// @{
  
  /// Sets `info` as the new `main()` trigger
  void replace(OpeningInfo_t info)
    { fMain = info; add(std::move(info)); }
  
  /**
   * @brief If `other` has fired, and at an earlier tick, set a new `main()`.
   * @param other another set of triggers
   * @return whether the `main()` information from `other` was copied
   * 
   * If `other` has `fired()` earlier than this one (or if this one hasn't fired
   * at all), the main trigger of `other` is adopted.
   * In any case all openings from `other` are `add()`'ed to this one.
   */
  bool addAndReplaceIfEarlier(TriggerInfo_t const& other);
  
  /**
   * @brief If `info` is earlier than `main()`, it is set as new `main()`.
   * @param info the opening information candidate as new `main()`
   * @return whether `main()` was updated.
   * 
   * If `info` is earlier tick than `main()` (or if there is no `main()` yet,
   * i.e. it has not `fired()`), `info` becomes the new `main()`.
   * In all cases `info` is added (`add()`).
   */
  bool addAndReplaceIfEarlier(OpeningInfo_t const& info);
    
  
  /// Adds an opening to `all` list (`main` is not affected). Not sorted.
  /// If no trigger is marked as `main()`, this becomes it.
  void add(OpeningInfo_t info)
    { if (fAll.empty()) fMain = info; fAll.push_back(std::move(info)); }
  
  /// Sorts `all` openings by time.
  void sortOpenings() { std::sort(fAll.begin(), fAll.end()); }
  
  /// @}
  // --- END ---- Modify trigger information -----------------------------------
  
  
  // --- BEGIN -- Access to trigger information --------------------------------
  /**
   * @name Access to trigger information
   * 
   * If the trigger did not fire, the result and behaviour of these methods are
   * undefined.
   */
  /// @{
  
  /// Returns the information of the main trigger (undefined if `!fired()`).
  OpeningInfo_t const& info() const { return main(); }
  
  /// Returns the time of the trigger (undefined if `!fired()`).
  optical_tick atTick() const { return main().tick; }
  
  /// Returns the time of the trigger (undefined if `!fired()`).
  optical_tick endTick() const { return main().endTick; }
  
  /// Returns the opening level of the trigger (undefined if `!fired()`).
  Opening_t level() const { return main().level; }
  
  /// Returns the ID of the location of the trigger (undefined if `!fired()`).
  LocationID_t location() const { return main().locationID; }
  
  /// Returns if the location of the trigger is set (undefined if `!fired()`).
  bool hasLocation() const { return main().hasLocation(); }
  
  /// Returns the full data (undefined if `!fired()`).
  OpeningInfo_t const& main() const { return fMain; }
  
  /// Returns all the registered opening in the current order (not resorted).
  /// @see `sortOpenings()`
  std::vector<OpeningInfo_t> const& all() const { return fAll; }
  
  /// Returns the number of registered triggers.
  std::size_t nTriggers() const { return all().size(); }
  
  /// @}
  // --- END -- Access to trigger information (if fired) -----------------------
  
    private:
  
  /// Main trigger (also found in `fAll`), if any.
  OpeningInfo_t fMain;
  
  ///< Information about all global trigger candidates.
  std::vector<OpeningInfo_t> fAll;
  
}; // icarus::trigger::details::TriggerInfo_t


//------------------------------------------------------------------------------
/**
 * @brief Helper to extract `OpeningInfo_t` from a trigger gate.
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarus::trigger::details::TriggerInfo_t triggerInfo;
 * icarus::trigger::details::GateOpeningInfoExtractor extractOpeningInfo
 *   { gate, { 6U } };
 * while (extractOpeningInfo) {
 *   auto info = extractOpeningInfo();
 *   if (info) triggerInfo.add(info.value());
 * } // while
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * fills `triggerInfo` with all the openings equal or above level `6U`.
 * Each opening is defined as from when `gate` reaches a specified threshold
 * ("opening threshold") to when it reaches or goes below another one ("closing
 * threshold"), with no dead time afterward.
 * The time of the opening is the time when threshold is passed, but
 * _the reported level is the maximum in the opening range_.
 * By default, the closing threshold is one less than the opening one (i.e. as
 * soon as the level goes below the opening threshold, the gate closes).
 * 
 * This class does not support multi-threading (it does have a mutable state).
 */
template <typename Gate>
class icarus::trigger::details::GateOpeningInfoExtractor {
  using Gate_t = Gate;
  using GateData_t = std::decay_t<decltype(gateDataIn(std::declval<Gate_t>()))>;
  
    public:
  using ClockTick_t = typename GateData_t::ClockTick_t;
  using OpeningCount_t = icarus::trigger::details::TriggerInfo_t::Opening_t;
  using OpeningInfo_t = icarus::trigger::details::TriggerInfo_t::OpeningInfo_t;
  using LocationID_t = icarus::trigger::details::TriggerInfo_t::LocationID_t;
  
  /// Configuration of the algorithm.
  struct Config_t {
    OpeningCount_t openThreshold { 1U };
    OpeningCount_t closeThreshold { 0U };
    unsigned int minWidth { 1U };
    unsigned int minGap { 0U };
    LocationID_t location { OpeningInfo_t::UnknownLocation };
    
    Config_t() = default;
    
    Config_t(
      OpeningCount_t openThreshold,
      OpeningCount_t closeThreshold,
      unsigned int minWidth = 1U,
      unsigned int minGap = 0U,
      LocationID_t location = OpeningInfo_t::UnknownLocation
      )
      : openThreshold(openThreshold), closeThreshold(closeThreshold)
      , minWidth(minWidth), minGap(minGap)
      , location(location)
      {}
    
    Config_t(OpeningCount_t threshold): Config_t(threshold, threshold - 1) {}
    
  }; // struct Config_t
  
  
  /**
   * @brief Constructor: uses default configuration.
   * @param gate the gate to operate on
   * @param config configuration of the algorithm
   * @see `configure()`
   */
  GateOpeningInfoExtractor(Gate_t const& gate)
    : GateOpeningInfoExtractor(gate, Config_t{})
    {}
  
  /**
   * @brief Constructor: sets all configuration parameters.
   * @param gate the gate to operate on
   * @param config configuration of the algorithm (see `configure()`)
   */
  GateOpeningInfoExtractor(Gate_t const& gate, Config_t config)
    : gateSrc(gate), config(std::move(config)), gate(gateDataIn(gateSrc))
    { restart(); }
  
  /**
   * @brief Constructor: sets all configuration parameters.
   * @param gate the gate to operate on
   * @param threshold both opening and closing threshold
   * @see `configure()`
   */
  GateOpeningInfoExtractor(Gate_t const& gate, OpeningCount_t threshold)
    : GateOpeningInfoExtractor(gate, Config_t{ threshold })
    {}
  
  /**
   * @brief Sets the configuration.
   * @brief Constructor: sets all configuration parameters.
   * @param gate the gate to operate on
   * @param config configuration of the algorithm
   * 
   * The configuration includes:
   *  * `openThreshold`, `closeThreshold`: the opening and closing thresholds
   *  * `minWidth`: minimum width of each opening range
   *  * `minGap`: minimum gap between successive opening ranges
   *  * `location`:  location ID to add to the produced `OpeningInfo_t` records
   */
  void configure(Config_t config) { config = std::move(config); }
  
  //@{
  /// The returned gate is at least `minWidth` ticks wide.
  std::optional<OpeningInfo_t> operator() () { return findNextOpening(); }
  
  std::optional<OpeningInfo_t> findNextOpening();
  //@}
  
  bool atEnd() const { return nextStart == MaxTick; }
  operator bool() const { return !atEnd(); }
  bool operator! () const { return atEnd(); }
  
  /// Resets the search from the specified time tick (beginning by default).
  void restart(ClockTick_t fromTick = MinTick)
    { nextStart = findOpen(fromTick); }
  
  
  /// @name Configuration access
  /// @{
  
  Config_t const& configuration() const { return config; }
  
  OpeningCount_t openThreshold() const { return config.openThreshold; }
  OpeningCount_t closeThreshold() const { return config.closeThreshold; }
  unsigned int minGap() const { return config.minGap; }
  unsigned int minWidth() const { return config.minWidth; }
  LocationID_t location() const { return config.location; }
  
  void setLocation(LocationID_t location) { config.location = location; }
  
  /// @}
  
  
    private:
  
  static constexpr ClockTick_t MinTick = GateData_t::MinTick;
  static constexpr ClockTick_t MaxTick = GateData_t::MaxTick;
  
  
  // --- BEGIN -- Configuration ------------------------------------------------
  Gate_t const& gateSrc;
  Config_t config;
  // --- END ---- Configuration ------------------------------------------------
  
//   GateData_t const& gate() { return gateDataIn(gateSrc); }
  GateData_t const& gate;
  
  ClockTick_t nextStart = MinTick;
  
  
  ClockTick_t findOpen(ClockTick_t start) const
    { return gate.findOpen(openThreshold(), start); }
  ClockTick_t findClose(ClockTick_t start) const
    { return gate.findClose(closeThreshold() + 1, start); }
  
  /// Returns the first closing and reopening above threshold from `start` on.
  std::pair<ClockTick_t, ClockTick_t> findNextCloseAndOpen
    (ClockTick_t start) const;
  
}; // class icarus::trigger::details::GateOpeningInfoExtractor<>


// -----------------------------------------------------------------------------
// ---  Inline implementation
// -----------------------------------------------------------------------------
// --- icarus::trigger::details::TriggerInfo_t
// -----------------------------------------------------------------------------
inline bool icarus::trigger::details::TriggerInfo_t::addAndReplaceIfEarlier
  (TriggerInfo_t const& other)
{
  if (!other.fired()) return false;
  
  bool const hadFired = fired();
  
  // register all triggers anyway; do not sort, do not resolve duplicates
  fAll.reserve(nTriggers() + other.nTriggers());
  for (OpeningInfo_t const& info: other.all()) add(info);
  
  if (hadFired && (other.atTick() >= atTick())) return false;
  
  fMain = other.main();
  return true;
} // icarus::trigger::details::TriggerInfo_t::addAndReplaceIfEarlier()


// -----------------------------------------------------------------------------
inline bool icarus::trigger::details::TriggerInfo_t::addAndReplaceIfEarlier
  (OpeningInfo_t const& info)
{
  add(info);
  if (info.tick >= atTick()) return false;
  
  fMain = info;
  return true;
} // icarus::trigger::details::TriggerInfo_t::addAndReplaceIfEarlier()


// -----------------------------------------------------------------------------
// ---  icarus::trigger::details::GateOpeningInfoExtractor<>
// -----------------------------------------------------------------------------
template <typename Gate>
auto icarus::trigger::details::GateOpeningInfoExtractor<Gate>::findNextOpening()
  -> std::optional<OpeningInfo_t>
{
  if (atEnd()) return {};
  
  using ClockDiff_t = decltype( ClockTick_t{} - ClockTick_t{} );
  
  ClockTick_t const start = nextStart;
  ClockTick_t closing;
  do {
    std::tie(closing, nextStart) = findNextCloseAndOpen(nextStart);
    if (nextStart == MaxTick) break;
  } while(
    (closing - start < static_cast<ClockDiff_t>(minWidth()))
    || (nextStart - closing < static_cast<ClockDiff_t>(minGap()))
  );
  
  return std::optional<OpeningInfo_t>{ std::in_place,
    detinfo::timescales::optical_tick{ start },
    detinfo::timescales::optical_tick{ closing },
    gate.openingCount(gate.findMaxOpen(start, closing)),
    location()
    };
} // icarus::trigger::details::GateOpeningInfoExtractor<>::findNextOpening()


// -----------------------------------------------------------------------------
template <typename Gate>
auto icarus::trigger::details::GateOpeningInfoExtractor<Gate>::findNextCloseAndOpen
  (ClockTick_t start) const ->  std::pair<ClockTick_t, ClockTick_t>
{
  ClockTick_t const closing = (gate.openingCount(start) > closeThreshold())
    ? findClose(start): start;
  return { closing, findOpen(closing) };
} // icarus::trigger::details::GateOpeningInfoExtractor<>::findNextCloseAndOpen()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_TRIGGERINFO_T_H
