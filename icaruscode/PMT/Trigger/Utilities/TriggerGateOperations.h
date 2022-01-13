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
#include "icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"

// C/C++ standard libraries
#include <iterator> // std::iterator_traits


namespace icarus::trigger {

  // --- BEGIN -- Gate operations ----------------------------------------------
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
  [[nodiscard]] GateObj discriminate(
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
  [[nodiscard]] auto sumGatesSequence(BIter const begin, EIter const end);


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
  [[nodiscard]] auto sumGates(GateColl const& gates);


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
  [[nodiscard]] auto sumGates
    (AGate A, BGate const& B, OTrigGates const&... others);
  
  
  // --- END -- Gate operations: sum -------------------------------------------
  
  
  // --- BEGIN -- Gate operations: multiplication ------------------------------
  
  /**
   * @brief Multiplies a sequence of gates.
   * @tparam BIter type of iterator to the gates to add
   * @tparam EIter type of iterator past-the-end of the sequence of gates to add
   * @param begin iterator to the first gate to add
   * @param end iterator past-the-last gate to add
   * @return a new gate product of all the specified ones
   * @see `mulGates()`
   * 
   * All gates from the first on are multiplied.
   * The returned gate has the same type as the first gate.
   */
  template <typename BIter, typename EIter>
  [[nodiscard]] auto mulGatesSequence(BIter const begin, EIter const end);


  /**
   * @brief Multiplies all the gates in a collection.
   * @tparam GateColl type of collection of gates
   * @param gates the collection of gates to multiply
   * @return a new gate product of all the `gates`
   * @see `mulGatesSequence()`
   * 
   * All gates from the first on are summed.
   * The returned gate has the same type as the gates in the collection.
   */
  template <typename GateColl>
  [[nodiscard]] auto mulGates(GateColl const& gates);


  /**
   * @brief Multiplies all the specified gates.
   * @tparam AGate type of the first gate to sum; it determines the return type
   * @tparam BGate type of the second gate to sum
   * @tparam OGates type of other gates to sum
   * @param A first gate to be multiplied
   * @param B second gate to be multiplied
   * @param others other gates to be multiplied
   * @return a new gate product of all the argument gates, of the same type as
   *         `A`
   * @see `mulGatesSequence()`
   * 
   * All arguments are combined in the product of the gates.
   * The product of two gates is at each tick the product of the level of
   * opening of the two gates at that tick.
   * The returned gate has the same type as the first gate.
   * 
   * @note There must be _at least two_ gates in the arguments.
   */
  template <typename AGate, typename BGate, typename... OTrigGates>
  [[nodiscard]] auto mulGates
    (AGate A, BGate const& B, OTrigGates const&... others);
  
  
  // --- END -- Gate operations: multiplication --------------------------------
  
  
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
  [[nodiscard]] auto maxGatesSequence(BIter const begin, EIter const end);


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
  [[nodiscard]] auto maxGates(GateColl const& gates);


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
  [[nodiscard]] auto maxGates
    (AGate A, BGate const& B, OTrigGates const&... others);
  
  
  // --- END -- Gate operations: max -------------------------------------------
  
  
  // --- BEGIN -- Gate operations: min  ----------------------------------------
  
  /**
   * @brief Computes the gate with the minimum opening of a sequence of gates.
   * @tparam BIter type of iterator to the gates to add
   * @tparam EIter type of iterator past-the-end of the sequence of gates to add
   * @param begin iterator to the first gate to add
   * @param end iterator past-the-last gate to add
   * @return a new gate minimum of all the specified ones
   * @see `minGates()`
   * 
   * For each tick, the minimum opening among all the gates in the sequence is
   * picked.
   * The returned gate has the same type as the first gate.
   */
  template <typename BIter, typename EIter>
  [[nodiscard]] auto minGatesSequence(BIter const begin, EIter const end);


  /**
   * @brief Computes the gate with the minimum opening of gates from collection.
   * @tparam GateColl type of collection of gates
   * @param gates the collection of gates to sum
   * @return a new gate minimum of all the `gates`
   * @see `minGatesSequence()`
   * 
   * For each tick, the minimum opening among all the gates in the collection is
   * picked.
   * The returned gate has the same type as the gates in the collection.
   */
  template <typename GateColl>
  [[nodiscard]] auto minGates(GateColl const& gates);


  /**
   * @brief Computes the gate with the minimum opening of the specified gates.
   * @tparam AGate type of the first gate; it determines the return type
   * @tparam BGate type of the second gate
   * @tparam OGates type of other gates
   * @param A first gate
   * @param B second gate
   * @param others other gates
   * @return a new gate minimum of all the argument gates, of same type as `A`
   * @see `minGatesSequence()`
   * 
   * For each tick, the minimum opening among all the specified gates is picked.
   * The returned gate has the same type as the first gate.
   * 
   * @note There must be _at least two_ gates in the arguments.
   */
  template <typename AGate, typename BGate, typename... OTrigGates>
  [[nodiscard]] auto minGates
    (AGate A, BGate const& B, OTrigGates const&... others);
  
  
  // --- END -- Gate operations: min -------------------------------------------


  // --- BEGIN -- Gate operations: generic -------------------------------------
  
  /**
   * @brief Applies an operation to two gates.
   * @tparam Op type of binary operation `AGate (Op)(AGate, BGate)`
   * @tparam AGate type of the first gate; it determines the return type
   * @tparam BGate type of the second gate
   * @tparam OGates type of other gates
   * @param op the binary operation to apply (copied)
   * @param A first gate
   * @param B second gate
   * @param others other gates
   * @return a new gate of same type as `A`, result of the operation.
   * @see `OpGatesSequence()`, `OpGateColl()`
   * 
   * For each tick, the operation `op` is performed between `A`, `B` and the
   * `others`, and the result is stored into the return value.
   * The order of operation is the same as in the arguments:
   * `((A op B) op others...)`.
   * The returned gate has the same type as the first gate (`A`).
   * 
   * Standard operations are provided in namespace `GateOps`.
   * 
   * A specific action is performed if the gates are "tracking" (i.e. instances
   * of `TrackedTriggerGate`), in which case tracking will also be propagated
   * so that the result includes all the tracked objects of all the operands.
   * 
   * @note There must be _at least two_ gates in the arguments.
   */
  template <typename Op, typename AGate, typename BGate, typename... OGates>
  [[nodiscard]] AGate OpGates
    (Op op, AGate A, BGate const& B, OGates const&... others);
  
  /**
   * @brief Computes the result of an operation on all gates from collection.
   * @tparam Op type of binary operation `AGate (Op)(AGate, BGate)`
   * @tparam GateColl type of collection of gates
   * @param op the binary operation to apply (copied)
   * @param gates the collection of gates to sum
   * @return a new gate result of `op` on all the `gates`
   * @see `OpGates()`, `OpGatesSequence()`
   * 
   * This is the sequential application of `op` to all `gates` via `OpGates()`,
   * according to their order in the collection.
   * The returned gate has the same type as the gates in the collection.
   */
  template <typename Op, typename GateColl>
  [[nodiscard]] auto OpGateColl(Op op, GateColl const& gates);
  
  /**
   * @brief Computes the result of an operation on all gates in the sequence.
   * @tparam Op type of binary operation `AGate (Op)(AGate, BGate)`
   * @tparam BIter type of iterator to the gates to add
   * @tparam EIter type of iterator past-the-end of the sequence of gates to add
   * @param op the binary operation to apply (copied)
   * @param begin iterator to the first gate to add
   * @param end iterator past-the-last gate to add
   * @return a new gate result of `op` on all the `gates`
   * @see `OpGates()`, `OpGateColl()`
   * 
   * This is the sequential application of `op` to all gates from `begin` to
   * `end` via `OpGates()`, following their order in the collection.
   * The returned gate has the same type as the one pointed by the begin
   * iterator (`BIter::value_type` for standard iterators).
   */
  template <typename Op, typename BIter, typename EIter>
  [[nodiscard]] auto OpGatesSequence(Op op, BIter const begin, EIter const end);
  
  // --- END ---- Gate operations: generic -------------------------------------
  
  
  /// Gate operations expressed as generic functions.
  namespace GateOps {
    
    /**
     * @name Binary gate operations.
     * 
     * The return gate is the same type as the first operand.
     */
    /// @{ 
    
    inline constexpr auto Min = [](auto A, auto const& B){ A.Min(B); return A; };
    
    inline constexpr auto Max = [](auto A, auto const& B){ A.Max(B); return A; };
    
    inline constexpr auto Sum = [](auto A, auto const& B){ A.Sum(B); return A; };
    
    inline constexpr auto Mul = [](auto A, auto const& B){ A.Mul(B); return A; };
    
    /// @}
    
  } // namespace GateOps
  
  
  /// @}
  // --- END -- Gate operations ------------------------------------------------
  
  
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
  // that includes tracking information, if any, which is preserved untouched
  auto discrGate { gate };
  
  auto const closeToOpen
    = static_cast<typename GateObj::OpeningDiff_t>(pass - fail);
  
  // ... except for the data, which we do want to touch:
  discrGate.clear();
  
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
// ---  Generic binary operations A -> A op B
// -----------------------------------------------------------------------------
template <typename Op, typename AGate, typename BGate, typename... OGates>
AGate icarus::trigger::OpGates
  (Op op, AGate A, BGate const& B, OGates const&... others)
{
  
  if constexpr(sizeof...(others) == 0U) { // two operands: end of recursion
    
    if constexpr(isTrackedTriggerGate_v<AGate>) {
      if constexpr(isTrackedTriggerGate_v<BGate>) {
        A.tracking().add(B.tracking());
      }
      A.gate() = op(std::move(A.gate()), gateIn(B));
      return A;
    } // if tracking
    else return op(std::move(A), B);
  }
  else {
    return OpGates(op, OpGates(op, std::move(A), B), others...);
  }
  
} // icarus::trigger::OpGates()


// -----------------------------------------------------------------------------
template <typename Op, typename GateColl>
auto icarus::trigger::OpGateColl(Op op, GateColl const& gates)
  { return OpGatesSequence(std::move(op), begin(gates), end(gates)); }


// -----------------------------------------------------------------------------
template <typename Op, typename BIter, typename EIter>
auto icarus::trigger::OpGatesSequence(Op op, BIter const begin, EIter const end)
{
  
  using Gate_t = typename std::iterator_traits<BIter>::value_type;
  
  // if `gates` is empty return a default-constructed gate
  if (begin == end) return Gate_t{};
  
  BIter iGate = begin;
  Gate_t resGate = *iGate;
  while (++iGate != end)
    resGate = OpGates(std::move(op), std::move(resGate), *iGate);
  
  return resGate;
  
} // icarus::trigger::OpGatesSequence()


// -----------------------------------------------------------------------------
// ---  Sum
// -----------------------------------------------------------------------------
template <typename AGate, typename BGate, typename... OGates>
auto icarus::trigger::sumGates(AGate A, BGate const& B, OGates const&... others)
  { return OpGates(GateOps::Sum, A, B, others...); }


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::sumGatesSequence(BIter const begin, EIter const end)
  { return OpGatesSequence(GateOps::Sum, begin, end); }


// -----------------------------------------------------------------------------
template <typename GateColl>
auto icarus::trigger::sumGates(GateColl const& gates)
  { return sumGatesSequence(begin(gates), end(gates)); }


// -----------------------------------------------------------------------------
// ---  Mul
// -----------------------------------------------------------------------------
template <typename AGate, typename BGate, typename... OGates>
auto icarus::trigger::mulGates(AGate A, BGate const& B, OGates const&... others)
  { return OpGates(GateOps::Mul, A, B, others...); }


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::mulGatesSequence(BIter const begin, EIter const end)
  { return OpGatesSequence(GateOps::Mul, begin, end); }


// -----------------------------------------------------------------------------
template <typename GateColl>
auto icarus::trigger::mulGates(GateColl const& gates)
  { return mulGatesSequence(begin(gates), end(gates)); }


// -----------------------------------------------------------------------------
// ---  Max
// -----------------------------------------------------------------------------
template <typename AGate, typename BGate, typename... OGates>
auto icarus::trigger::maxGates(AGate A, BGate const& B, OGates const&... others)
  { return OpGates(GateOps::Max, A, B, others...); }


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::maxGatesSequence(BIter const begin, EIter const end)
  { return OpGatesSequence(GateOps::Max, begin, end); }


// -----------------------------------------------------------------------------
template <typename GateColl>
auto icarus::trigger::maxGates(GateColl const& gates)
  { return maxGatesSequence(begin(gates), end(gates)); }


// -----------------------------------------------------------------------------
// ---  Min
// -----------------------------------------------------------------------------
template <typename AGate, typename BGate, typename... OGates>
auto icarus::trigger::minGates(AGate A, BGate const& B, OGates const&... others)
  { return OpGates(GateOps::Min, A, B, others...); }


// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::minGatesSequence(BIter const begin, EIter const end)
  { return OpGatesSequence(GateOps::Min, begin, end); }


// -----------------------------------------------------------------------------
template <typename GateColl>
auto icarus::trigger::minGates(GateColl const& gates)
  { return minGatesSequence(begin(gates), end(gates)); }


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEOPERATIONS_H
