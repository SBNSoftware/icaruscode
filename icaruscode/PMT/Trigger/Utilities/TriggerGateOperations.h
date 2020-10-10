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

  /// --- BEGIN -- Gate operations ---------------------------------------------
  /**
   * @name Gate operations
   *
   * Currently the following operations are supported:
   * 
   * * discrimination against a threshold: `discriminate()`;
   * * sum of trigger gates: `sumGates()`, `sumGatesSequence()`;
   * * maximum of trigger gates: `maxGates()`, `maxGatesSequence()`.
   *
   * Operations on more than one gate can take a sequence (begin and end
   * iterators), a collection or an arbitrary number of gates.
   */
  /// @{

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
  
  
  // --- BEGIN -- Gate operations: sum  ----------------------------------------
  
  /**
   * @brief Sums a sequence of gates.
   * @tparam BIter type of iterator to the gates to add
   * @tparam EIter type of iterator past-the-end of the sequence of gates to add
   * @param begin iterator to the first gate to add
   * @param end iterator past-the-last gate to add
   * @return a new gate sum of all the specified ones
   * @see `sumGates()`
   * 
   * All gates from the first on are summed.
   * The returned gate has the same type as the first gate.
   */
  template <typename BIter, typename EIter>
  auto sumGatesSequence(BIter const begin, EIter const end);


  /**
   * @brief Sums all the gates in a collection.
   * @tparam GateColl type of collection of gates
   * @param gates the collection of gates to sum
   * @return a new gate sum of all the `gates`
   * @see `sumGatesSequence()`
   * 
   * All gates from the first on are summed.
   * The returned gate has the same type as the gates in the collection.
   */
  template <typename GateColl>
  auto sumGates(GateColl const& gates);


  /**
   * @brief Sums all the specified gates.
   * @tparam AGate type of the first gate to sum; it determines the return type
   * @tparam BGate type of the second gate to sum
   * @tparam OGates type of other gates to sum
   * @param A first gate to be summed
   * @param B second gate to be summed
   * @param others other gates to be summed
   * @return a new gate sum of all the argument gates, of the same type as `A`
   * @see `sumGatesSequence()`
   * 
   * All arguments are combined in the sum of the gates.
   * The returned gate has the same type as the first gate.
   * 
   * @note There must be _at least two_ gates in the arguments.
   */
  template <typename AGate, typename BGate, typename... OTrigGates>
  auto sumGates(AGate A, BGate const& B, OTrigGates const&... others);
  
  
  /// --- END -- Gate operations: sum ------------------------------------------
  
  
  // --- BEGIN -- Gate operations: max  ----------------------------------------
  
  /**
   * @brief Computes the gate with the maximum opening of a sequence of gates.
   * @tparam BIter type of iterator to the gates to add
   * @tparam EIter type of iterator past-the-end of the sequence of gates to add
   * @param begin iterator to the first gate to add
   * @param end iterator past-the-last gate to add
   * @return a new gate maximum of all the specified ones
   * @see `maxGates()`
   * 
   * For each tick, the maximum opening among all the gates in the sequence is
   * picked.
   * The returned gate has the same type as the first gate.
   */
  template <typename BIter, typename EIter>
  auto maxGatesSequence(BIter const begin, EIter const end);


  /**
   * @brief Computes the gate with the maximum opening of gates from collection.
   * @tparam GateColl type of collection of gates
   * @param gates the collection of gates to sum
   * @return a new gate maximum of all the `gates`
   * @see `maxGatesSequence()`
   * 
   * For each tick, the maximum opening among all the gates in the collection is
   * picked.
   * The returned gate has the same type as the gates in the collection.
   */
  template <typename GateColl>
  auto maxGates(GateColl const& gates);


  /**
   * @brief Computes the gate with the maximum opening of the specified gates.
   * @tparam AGate type of the first gate; it determines the return type
   * @tparam BGate type of the second gate
   * @tparam OGates type of other gates
   * @param A first gate
   * @param B second gate
   * @param others other gates
   * @return a new gate maximum of all the argument gates, of same type as `A`
   * @see `maxGatesSequence()`
   * 
   * For each tick, the maximum opening among all the specified gates is picked.
   * The returned gate has the same type as the first gate.
   * 
   * @note There must be _at least two_ gates in the arguments.
   */
  template <typename AGate, typename BGate, typename... OTrigGates>
  auto maxGates(AGate A, BGate const& B, OTrigGates const&... others);
  
  
  /// --- END -- Gate operations: sum ------------------------------------------
  
  
  /// @}
  /// --- END -- Gate operations -----------------------------------------------
  
  
} // namespace icarus::trigger


// =============================================================================
// ===  template implementation
// =============================================================================
// ---  Discriminate
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
// ---  Sum
// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::sumGatesSequence(BIter const begin, EIter const end) {
  
  // if `gates` is empty return a default-constructed gate of the contained type
  if (begin == end) return decltype(*begin){};
  
  auto iGate = begin;
  auto resGate = *iGate;
  while (++iGate != end) resGate.Sum(*iGate);
  
  return resGate;
  
} // icarus::trigger::sumGatesSequence()


// -----------------------------------------------------------------------------
template <typename GateColl>
auto icarus::trigger::sumGates(GateColl const& gates)
  { return sumGatesSequence(begin(gates), end(gates)); }


// -----------------------------------------------------------------------------
template <typename AGate, typename BGate, typename... OGates>
auto icarus::trigger::sumGates(AGate A, BGate const& B, OGates const&... others)
{
  if constexpr(sizeof...(others) == 0U) { // only two operands: end of recursion
    A.Sum(B); // A is already a copy...
    return A;
  }
  else {
    return sumGates(sumGates(A, B), others...);
  }
} // icarus::trigger::sumGates()


// -----------------------------------------------------------------------------
// ---  Max
// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::maxGatesSequence(BIter const begin, EIter const end) {
  
  // if `gates` is empty return a default-constructed gate of the contained type
  if (begin == end) return decltype(*begin){};
  
  auto iGate = begin;
  auto resGate = *iGate;
  while (++iGate != end) resGate.Max(*iGate);
  
  return resGate;
  
} // icarus::trigger::maxGatesSequence()


// -----------------------------------------------------------------------------
template <typename GateColl>
auto icarus::trigger::maxGates(GateColl const& gates)
  { return maxGatesSequence(begin(gates), end(gates)); }


// -----------------------------------------------------------------------------
template <typename AGate, typename BGate, typename... OGates>
auto icarus::trigger::maxGates(AGate A, BGate const& B, OGates const&... others)
{
  if constexpr(sizeof...(others) == 0U) { // only two operands: end of recursion
    A.Max(B); // A is already a copy...
    return A;
  }
  else {
    return maxGates(maxGates(A, B), others...);
  }
} // icarus::trigger::maxGates()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEOPERATIONS_H
