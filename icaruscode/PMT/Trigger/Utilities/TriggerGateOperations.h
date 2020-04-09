/**
 * @file   icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h
 * @brief  Utilities for the conversion of trigger gate data formats.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 3, 2020
 * 
 * This library is header-only.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEOPERATIONS_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEOPERATIONS_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
// #include "larcorealg/CoreUtils/enumerate.h"

// C/C++ standard libraries
// #include <vector>
// #include <tuple>
// #include <utility> // std::move()


namespace icarus::trigger {

  // ---------------------------------------------------------------------------

  /**
   * @brief Returns a discriminated version of `gate`.
   * @tparam GateObj type of gate being discriminated (and returned)
   * @param gate the gate to be discriminated
   * @param threshold (_default: `1`_) the discrimination threshold
   * @param pass (_default: `1`_) discriminated gate value on test pass
   * @param fail (_default: `0`_) discriminated gate value on test fail
   * @return a gate resulting from the discrimination of `gate`.
   * @see `gateAbove()`, `gateBelow()`
   * 
   * A new gate of the same type as the input `gate` is returned.
   * This gate has two opening count values: either `pass`, in the time
   * intervals where the `gate` is at or above `threshold`, or `fail`,
   * in the time intervals where the `gate` is below `threshold`.
   * 
   * Examples:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const discrGate = discriminate(gate, 5U, 0U, 5U);
   * auto const discrGateNeg = discriminate(gate, 5U, 5U, 0U);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will have `discrGate` as a gate with value `0U` wherever `gate` has opening
   * count 4 or less, and `5U` where `gate` has opening count 5 or more.
   * The gate` discrGateNeg` has the two opening values swapped and therefore
   * results the complement of `discrGate`.
   */
  template <typename GateObj>
  GateObj discriminate(
    GateObj const& gate,
    typename GateObj::OpeningCount_t threshold = 1U,
    typename GateObj::OpeningCount_t pass = 1U,
    typename GateObj::OpeningCount_t fail = 0U
    );
  
  // ---------------------------------------------------------------------------

} // namespace icarus::trigger


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename GateObj>
GateObj icarus::trigger::discriminate(
  GateObj const& gate,
  typename GateObj::OpeningCount_t threshold /* = 1U */,
  typename GateObj::OpeningCount_t pass /* = 1U */,
  typename GateObj::OpeningCount_t fail /* = 0U */
  )
{
  // we copy the gate hoping that the copy constructor takes care of everything
  // else that we don't want to touch...
  auto discrGate { gate };
  
  auto const closeToOpen
    = static_cast<typename GateObj::OpeningDiff_t>(pass - fail);
  
  // ... except for the data, which we do want to touch
  
  // set the starting level according to the discrimination
  auto lastTick = discrGate.MinTick;
  if (gate.openingCount(lastTick) < threshold) {
    // gate starts from fail (most "normal" case)
    discrGate.setOpeningAt(lastTick, fail);
  }
  else {
    // gate starts from pass!
    discrGate.setOpeningAt(lastTick, pass);
    
    // bring it to fail
    lastTick = gate.findClose(threshold, ++lastTick);
    if (lastTick != gate.MaxTick) discrGate.closeAt(lastTick, closeToOpen);
  } // if started from above threshold
  
  // we are at a fail state now
  while (lastTick < gate.MaxTick) {
    
    // bring it to pass...
    lastTick = gate.findOpen(threshold, ++lastTick);
    if (lastTick == gate.MaxTick) break;
    discrGate.openAt(lastTick, closeToOpen);
    
    // ... back to fail...
    lastTick = gate.findClose(threshold, ++lastTick);
    if (lastTick == gate.MaxTick) break;
    discrGate.closeAt(lastTick, closeToOpen);
    
  } // while last tick
  
  return discrGate;
} // icarus::trigger::discriminate()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEOPERATIONS_H
