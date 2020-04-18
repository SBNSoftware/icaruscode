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

// C++ standard libraries
#include <optional>
#include <utility> // std::forward()


// -----------------------------------------------------------------------------
namespace icarus::trigger::details { struct TriggerInfo_t; }

// -----------------------------------------------------------------------------
/// Helper data structure to store transient trigger result.
struct icarus::trigger::details::TriggerInfo_t {
  
  using optical_tick = detinfo::timescales::optical_tick; ///< Type alias.
  
  /// Optional trigger information; present if the trigger fired.
  struct Info_t {
    optical_tick tick; ///< Tick at which the trigger fired.
  }; // Info_t
  
  
  // --- BEGIN -- Construction -------------------------------------------------
  
  TriggerInfo_t() = default; // no trigger
  TriggerInfo_t(Info_t const& info): fInfo(info) {}
  TriggerInfo_t(Info_t&& info): TriggerInfo_t(info) {} // copy
  
  /// Reinitializes the object by constructing a `Info_t` with `args`.
  template <typename... Args>
  void emplace(Args&&... args)
    { fInfo.emplace(Info_t{ std::forward<Args>(args)... }); }
  
  
  /// If `other` has fired, and at an earlier tick, copy `other`'s information.
  /// @return Whether the information from `other` was copied.
  bool replaceIfEarlier(TriggerInfo_t const& other);
    
  
  // --- END -- Construction ---------------------------------------------------
  
  
  // --- BEGIN -- Query whether the trigger fired ------------------------------
  /// @name Query whether the trigger fired.
  /// @{
  
  /// Returns whether the trigger fired.
  bool fired() const { return fInfo.has_value(); }
  
  /// Returns whether there is trigger information.
  operator bool() const { return fired(); }
  
  /// Returns whether there is no trigger information.
  bool operator! () const { return !fired(); }
  
  /// @}
  // --- END -- Query whether the trigger fired --------------------------------
  
  
  // --- BEGIN -- Access to trigger information --------------------------------
  /**
   * @name Access to trigger information
   * 
   * If the trigger did not fire, the result and behaviour of these methods are
   * undefined.
   */
  /// @{
  
  /// Returns the full data (undefined behaviour if `!fired()`).
  Info_t const& info() const { return fInfo.value(); }
  
  /// Returns the time of the trigger (undefined behaviour if `!fired()`).
  optical_tick atTick() const { return fInfo->tick; }
  
  /// @}
  // --- END -- Access to trigger information (if fired) -----------------------
  
    private:
  
  std::optional<Info_t> fInfo; ///< Fired trigger information, present if fired.
  
}; // icarus::trigger::details::TriggerInfo_t


// -----------------------------------------------------------------------------
// ---  Inline implementation
// -----------------------------------------------------------------------------
bool icarus::trigger::details::TriggerInfo_t::replaceIfEarlier
  (TriggerInfo_t const& other)
{
  if (!other.fired()) return false;
  if (fired() && (other.atTick() >= atTick())) return false;
  
  fInfo = other.info();
  return true;
} // icarus::trigger::details::TriggerInfo_t::replaceIfEarlier()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_TRIGGERINFO_T_H
