/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h
 * @brief  Helper class to store transient trigger result.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 15, 2020
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_TRIGGERINFO_T_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_TRIGGERINFO_T_H

// ICARUS libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_tick
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h" // 

// C++ standard libraries
#include <vector>
#include <algorithm>
#include <limits> // std::numeric_limits<>
#include <utility> // std::forward()


// -----------------------------------------------------------------------------
namespace icarus::trigger::details { struct TriggerInfo_t; }

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
    
    Opening_t level = 0U; ///< Maximum level on the main trigger opening.
    
    /// Identified of the trigger location.
    LocationID_t locationID = UnknownLocation;
    
    
    /// Default constructor.
    OpeningInfo_t() = default;
    
    /// Constructor: specify `tick`, the rest is optional.
    OpeningInfo_t(
      optical_tick tick,
      Opening_t level = 0U, std::size_t locationID = UnknownLocation
      )
      : tick(tick), level(level), locationID(locationID) {}
    
    
    /// Comparison: order from time.
    bool operator< (OpeningInfo_t const& other) const
      { return tick < other.tick; }
    
  }; // OpeningInfo_t
  
  
  // --- BEGIN -- Construction -------------------------------------------------
  
  TriggerInfo_t() = default; // no trigger
  TriggerInfo_t(OpeningInfo_t info) { replace(std::move(info)); }
  
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
  
  /// Returns the opening level of the trigger (undefined if `!fired()`).
  Opening_t level() const { return main().level; }
  
  /// Returns the ID of the location of the trigger (undefined if `!fired()`).
  LocationID_t location() const { return main().locationID; }
  
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


// -----------------------------------------------------------------------------
// ---  Inline implementation
// -----------------------------------------------------------------------------
bool icarus::trigger::details::TriggerInfo_t::addAndReplaceIfEarlier
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
bool icarus::trigger::details::TriggerInfo_t::addAndReplaceIfEarlier
  (OpeningInfo_t const& info)
{
  add(info);
  if (info.tick >= atTick()) return false;
  
  fMain = info;
  return true;
} // icarus::trigger::details::TriggerInfo_t::addAndReplaceIfEarlier()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_TRIGGERINFO_T_H
