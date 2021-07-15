/**
 * @file   icaruscode/PMT/Trigger/Data/TriggerGateData.tcc
 * @brief  A logical multilevel gate for triggering (template implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Data/TriggerGateData.h`
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_TCC
#define ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_TCC

// LArSoft libraries
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()

// C/C++ standard libraries
#include <ostream>
#include <stdexcept> // std::runtime_error
#include <algorithm> // std::min(), std::max(), std::upper_bound()
#include <utility> // std::move(), std::swap()
#include <functional> // std::plus<>, std::multiplies<>
#include <iterator> // std::prev(), std::next()
#include <cassert>


// make "sure" this header is not included directly
#ifndef ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_H
# error "TriggerGateData.tcc must not be directly included!"\
        " #include \"icaruscode/PMT/Trigger/Data/TriggerGateData.h\" instead."
#endif // !ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_H


//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::TriggerGateData<>
//------------------------------------------------------------------------------
template <typename TK, typename TI>
typename icarus::trigger::TriggerGateData<TK, TI>::Status const
icarus::trigger::TriggerGateData<TK, TI>::NewGateStatus
  { EventType::Set, MinTick, 0U };


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findOpen(
  OpeningCount_t minOpening /* = 1U */,
  ClockTick_t start /* = MinTick */, ClockTick_t end /* = MaxTick */
) const -> ClockTick_t
{
  auto iStatus = findOpenStatus(minOpening, start, end);
  return (iStatus == fGateLevel.end())? end: iStatus->tick;
} // icarus::trigger::TriggerGateData<>::findOpen()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findClose(
  OpeningCount_t minOpening /* = 1U */,
  ClockTick_t start /* = MinTick */, ClockTick_t end /* = MaxTick */
) const -> ClockTick_t
{
  auto iStatus = findCloseStatus(minOpening, start, end);
  return (iStatus == fGateLevel.end())? end: iStatus->tick;
} // icarus::trigger::TriggerGateData<>::findClose()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findMaxOpen
  (ClockTick_t start /* = MinTick */, ClockTick_t end /* = MaxTick */) const
  -> ClockTick_t
{
  auto iStatus = findMaxOpenStatus(start, end);
  // the iterator could point to a status change before the target interval...
  return (iStatus == fGateLevel.end())? end: std::max(iStatus->tick, start);
} // icarus::trigger::TriggerGateData<>::findMaxOpen()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
void icarus::trigger::TriggerGateData<TK, TI>::compact() {
  
  /*
   * Removes:
   * * multiple sets at the same tick
   * * shifts that don't change opening level
   * 
   * Keeps:
   * * unknown stati (because they are unknown and they might have reason to be)
   * 
   */
  
  auto const send = fGateLevel.end();
  auto iLast = fGateLevel.begin(); // last good status
  auto iDest = std::next(iLast); // the status next to be assigned
  auto iTest = iDest; // the next status candidate
  
  // declare the current candidate status good:
  auto accept = [&iDest, &iTest, &iLast]()
    { if (iDest != iTest) *iDest = *iTest; iLast = iDest; ++iDest; };
  
  while (iTest != send) {
    
    switch (iTest->event) {
      case EventType::Shift:
        if (iLast->opening != iTest->opening) accept();
        break;
      case EventType::Set:
        if (iTest->tick > iLast->tick) accept();
        break;
      case EventType::Unknown:
        accept();
        break;
    } // switch
    
    ++iTest;
    
  } // while
  
  fGateLevel.erase(iDest, send);
  
} // icarus::trigger::TriggerGateData<>::compact()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::openingRange
  (ClockTick_t start, ClockTick_t end) const
  -> std::pair<OpeningCount_t, OpeningCount_t> 
{
  /*
   * parse from the first to the last relevant status of the gate and collect
   * the extrema of their opening counts
   */
  assert(!fGateLevel.empty());
  
  if (start >= end) return { {}, {} }; // this should rank pretty high in the
                                       // competition of weird C++ constructs
  
  // look where to start from; if the start tick is before the start of the
  // gate, fast forward to the actual start of the gate
  auto const send = fGateLevel.end();
  auto const maybeStatusIter = findLastStatusFor(start);
  auto iStatus = maybeStatusIter? maybeStatusIter.value(): fGateLevel.begin();
  assert(iStatus != fGateLevel.end());
  
  // so the requested range is fully before the start of the gate: bail out!
  if (iStatus->tick >= end) return { {}, {} };
  
  OpeningCount_t const startOpening = iStatus->opening;
  std::pair<OpeningCount_t, OpeningCount_t> limits
    { startOpening, startOpening + 1 };
  while (++iStatus != send) {
    if (iStatus->tick >= end) break; // we reached the end of the tick interval
    
    auto const opening = iStatus->opening;
    if (opening < limits.first) limits.first = opening;
    if (opening >= limits.second) limits.second = opening + 1;
    
  } // while
  
  return limits;
} // icarus::trigger::TriggerGateData<>::openingRange()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::openingCount
  (ClockTick_t tick) const -> OpeningCount_t
{
  return findLastStatusForTickOrThrow(tick)->opening;
} // icarus::trigger::TriggerGateData<>::openingCount()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
void icarus::trigger::TriggerGateData<TK, TI>::setOpeningAt
  (ClockTick_t tick, OpeningCount_t openingCount)
{
  //
  // first find where to start acting
  //
  auto iStatus = findLastStatusForTickOrThrow(tick); // may be before the tick
  
  //
  // set the current status (possibly a new one)
  //
  if (iStatus->tick == tick) { // overwrite the previous action
    iStatus->event = EventType::Set;
    iStatus->opening = openingCount;
  }
  else { // no status exactly at this tick, iStatus is just before
    // insert the new status after iStatus
    iStatus = fGateLevel.insert
      (++iStatus, { EventType::Set, tick, openingCount });
  }
  
  //
  // update all the following stati;
  // the "Set" action stops when another set happens, but overrides shifts
  //
  auto const send = fGateLevel.end();
  while (++iStatus != send) {
    switch (iStatus->event) {
      case EventType::Shift: // just remove the shift event
        iStatus = fGateLevel.erase(iStatus);
        break;
      case EventType::Set: // the later Set event takes over, we are done
        return;
      case EventType::Unknown: // not sure about this... let's keep going
        break;
    } // switch event type
  }
  
} // icarus::trigger::TriggerGateData<>::setOpeningAt()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
void icarus::trigger::TriggerGateData<TK, TI>::openBetween
  (ClockTick_t start, ClockTick_t end, OpeningDiff_t count /* = 1 */)
{
  /*
   * The plan is:
   * 1) find the even just before (or at) the tick where we want the gate opened
   * 2) open the gate at that location
   * 3) find the last existing event in the gate that we need to change, by tick
   * 4) update all the events to reflect the level change
   *    1) if the event is a (re)set event, it takes priority and we are done
   * 5) insert the gate close event at the right tick
   */
  if (start >= end) return; // weird, yet valid
  
  //
  // (1) first find where to start acting
  //
  auto iStatus = findLastStatusForTickOrThrow(start); // may be before the start
  
  //
  // (2) set the current status (possibly a new one)
  //
  if (iStatus->tick == start) { // there is already something here
    switch (iStatus->event) {
      case EventType::Unknown: // the shift overrides... we don't even know what
        iStatus->event = EventType::Shift;
        [[fallthrough]];
      case EventType::Shift:
      case EventType::Set: // in case of Set, we change the set value
        iStatus->opening += count;
        break;
    } // switch
  }
  else { // no status exactly at start tick, iStatus is just before
    // insert the new status after iStatus
    auto const opening = iStatus->opening + count;
    iStatus = fGateLevel.insert
      (++iStatus, { EventType::Shift, start, opening });
  }
  
  //
  // (3) also find where to stop acting
  //
  auto send = findLastStatusForTickOrThrow(end); // may be before `end`
  // if `send` is stricly earlier than `end`, we want it affected too
  if (send->tick < end) ++send;
  
  //
  // (4) update all the following stati;
  // the "Shift" action stops when a set happens, and stacks with other shifts
  //
  while (++iStatus != send) {
    switch (iStatus->event) {
      case EventType::Shift: // change the resulting opening
        if ((count < 0) && (iStatus->opening < OpeningCount_t(-count))) {
          throw std::runtime_error(
            "icarus::trigger::TriggerGateData::openBetween(): "
            "asked to close " + util::to_string(-count)
            + " gate counts starting at "
            + util::to_string(start) + " but at time "
            + util::to_string(iStatus->tick) + " only "
            + util::to_string(iStatus->opening) + " are still open"
            );
        }
        iStatus->opening += count;
        break;
      case EventType::Set: // (4.1) Set event takes over, we are done
        return;
      case EventType::Unknown: // not sure about this... let's keep going
        break;
    } // switch event type
  } // while
  
  //
  // (5) finally close the gate: mark the end, and leave the following alone
  // Note that gate is unlimited in time.
  //
  
  if ((iStatus != fGateLevel.end()) && (iStatus->tick == end)) {
    // there is already something here, and it must be at the right level
    switch (iStatus->event) {
      case EventType::Unknown: // the shift overrides... we don't even know what
        iStatus->event = EventType::Shift;
        [[fallthrough]];
      case EventType::Shift:
      case EventType::Set:
        break;
    } // switch
  }
  else { // no status exactly at this tick, just need to insert one
    // insert the new status before iStatus;
    // the correct opening now is the last one, minus what we added;
    // it is also the status in iStatus, if it exists
    assert(iStatus != fGateLevel.begin()); // we must have added one status!
    
    auto const opening = (iStatus == fGateLevel.end())
      ? (std::prev(iStatus)->opening - count): (iStatus->opening);
    
    // checking that if the two branches above are both available, they match;
    // the branch above just takes the easier way when available
    assert((iStatus == fGateLevel.end())
      || (opening == (std::prev(iStatus)->opening - count))
      );
    
    iStatus = fGateLevel.insert(iStatus, { EventType::Shift, end, opening });
  }
  // if we get here, now iStatus contains the gate closing (for what we care)
  
} // icarus::trigger::TriggerGateData<>::openFor()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Min
  (TriggerGateData const& other) -> triggergatedata_t&
{
  fGateLevel = std::move(TriggerGateData::Min(*this, other).fGateLevel);
  return *this;
} // icarus::trigger::TriggerGateData<>::Min()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Max
  (TriggerGateData const& other) -> triggergatedata_t&
{
  fGateLevel = std::move(TriggerGateData::Max(*this, other).fGateLevel);
  return *this;
} // icarus::trigger::TriggerGateData<>::Max()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Sum
  (TriggerGateData const& other) -> triggergatedata_t&
{
  fGateLevel = std::move(TriggerGateData::Sum(*this, other).fGateLevel);
  return *this;
} // icarus::trigger::TriggerGateData<>::Sum()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Mul
  (TriggerGateData const& other) -> triggergatedata_t&
{
  fGateLevel = std::move(TriggerGateData::Mul(*this, other).fGateLevel);
  return *this;
} // icarus::trigger::TriggerGateData<>::Mul()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Min
  (TriggerGateData const& a, TriggerGateData const& b) -> triggergatedata_t
{
  return SymmetricCombination
    ([](OpeningCount_t a, OpeningCount_t b){ return std::min(a, b); }, a, b);
} // icarus::trigger::TriggerGateData<>::Min()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Max
  (TriggerGateData const& a, TriggerGateData const& b) -> triggergatedata_t
{
  return SymmetricCombination
    ([](OpeningCount_t a, OpeningCount_t b){ return std::max(a, b); }, a, b);
} // icarus::trigger::TriggerGateData<>::Max()

  
//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Sum
  (TriggerGateData const& a, TriggerGateData const& b) -> triggergatedata_t
{
  return SymmetricCombination(std::plus<OpeningCount_t>(), a, b);
} // icarus::trigger::TriggerGateData<>::Sum()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::Mul
  (TriggerGateData const& a, TriggerGateData const& b) -> triggergatedata_t
{
  return SymmetricCombination(std::multiplies<OpeningCount_t>(), a, b);
} // icarus::trigger::TriggerGateData<>::Mul()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
template <typename Op>
auto icarus::trigger::TriggerGateData<TK, TI>::SymmetricCombination(
  Op&& op, triggergatedata_t const& a, triggergatedata_t const& b,
  ClockTicks_t aDelay /* = { 0 } */, ClockTicks_t bDelay /* { = 0 } */
) -> triggergatedata_t {
  /*
   * Combines two gates into a new one.
   * The new gate is created directly event by event.
   * All events following tick #0 are changed to `Shift` if they shift level,
   * or otherwise discarded.
   * 
   */
  
  struct State {
    
    /// Keeps track of information travelling with a gate.
    struct GateStatus_t {
      // protecting against Clang bug 33298 (at least until Clang 8)
      // (https://bugs.llvm.org/show_bug.cgi?id=33298)
#if defined(__clang__) && (__clang_major__ < 9)
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-local-typedef"
#endif // __clang_major__ < 9
      using iterator_t = triggergatedata_t::status_const_iterator;
#if defined(__clang__) && (__clang_major__ < 9)
# pragma clang diagnostic pop
#endif // __clang_major__ < 9
      
      triggergatedata_t const& triggerData;
      iterator_t iStatus;
      ClockTicks_t delay { 0 };
      OpeningCount_t prevLevel = 0U;
      
      iterator_t operator->() const { return iStatus; }
      decltype(auto) operator*() const { return *iStatus; }
      
      GateStatus_t(
        triggergatedata_t const& triggerData,
        ClockTicks_t delay = ClockTicks_t{ 0 },
        OpeningCount_t prevLevel = 0U
        )
        : triggerData(triggerData)
        , iStatus(begin())
        , delay(delay)
        , prevLevel(prevLevel)
        {
          assert(iStatus != end());
          assert(iStatus->event == EventType::Set);
        }
      
      /// Returns the tick of the current status (delay is properly applied).
      ClockTick_t tick() const { return delayedTick(iStatus->tick); }
      
      /// Returns the current opening level.
      OpeningCount_t opening() const { return iStatus->opening; }
      
      /// Returns whether current opening level is different from the previous.
      bool changed() const { return opening() != prevLevel; }
      
      GateStatus_t& operator++()
        { prevLevel = opening(); ++iStatus; return *this; }
      
      auto begin() const { return triggerData.fGateLevel.begin(); }
      auto end() const { return triggerData.fGateLevel.end(); }
      
      bool atEnd() const { return iStatus == end(); }
      
      ClockTick_t delayedTick(ClockTick_t tick) const { return tick + delay; }
      
    }; // GateStatus_t
    
    
    State(
      triggergatedata_t const& left,
      triggergatedata_t const& right,
      ClockTicks_t leftDelay = ClockTicks_t{ 0 },
      ClockTicks_t rightDelay = ClockTicks_t{ 0 }
      )
      : leftStatus(left, leftDelay)
      , rightStatus(right, rightDelay)
      , pCurrent(&leftStatus)
      , pOther(&rightStatus)
      {
        // choose which gate to start from
        if (pCurrent->tick() > pOther->tick()) swap();
      } // State()
    
    /// Returns the gate status to be considered.
    GateStatus_t& current() { return *pCurrent; }
    
    /// Returns the gate status that is not current.
    GateStatus_t& other() { return *pOther; }
    
    /// Applies the specified operation on the current status.
    decltype(auto) apply(Op& op) const
      { return op(pOther->prevLevel, pCurrent->opening()); }
    
    /// Returns whether both sides are over.
    bool atEnd() const { return pCurrent->atEnd(); }
    
    /// Prepares for the next entry, and returns false if we are done.
    bool next()
      {
        ++(*pCurrent);
        if (pCurrent->atEnd()) {
          swap();
          return !pCurrent->atEnd();
        }
        if (!pOther->atEnd() && (pOther->tick() < pCurrent->tick())) swap();
        return true;
      }
    
      private:
    
    GateStatus_t leftStatus; ///< Helper for the left gate.
    GateStatus_t rightStatus; ///< Helper for the right gate.
    GateStatus_t* pCurrent; ///< To the gate with the next time tick.
    GateStatus_t* pOther; ///< To the gate not with the next time tick.
    
    /// Swap the current and the other gate status.
    void swap() { std::swap(pCurrent, pOther); }
    
  }; // struct State
  
  
  // prepare the container of the combination:
  triggergatedata_t result;
  auto& resultLevels = result.fGateLevel;
  resultLevels.reserve(a.fGateLevel.size() + b.fGateLevel.size() - 1U);
  
  // the state automatically selects the earliest as start
  State state(a, b, aDelay, bDelay);
  
  // we set the starting level in a "special" way, not combining anything
  OpeningCount_t const startOpening = state.current().opening();
  resultLevels.back().opening = startOpening;
  // equalize the previous gate opening;
  // this is relevant if the first gate has more than one status before the
  // other starts, in which case the next action will be an attempt to combine
  // the second status of the early gate with the "previous" opening value of
  // the late gate; we choose that value to be the same as the other gate.
  state.other().prevLevel = startOpening;
  
  while (state.next()) {
    
    auto const& newStatus = state.current();
    
    // quick way out: if the event does not change its gate level, it's out
    if (!newStatus.changed()) continue;
    
    // add a new status if more recent that the current one
    OpeningCount_t const newLevel
      = op(state.other().prevLevel, newStatus.opening());
    
    ClockTick_t const newTick = newStatus.tick();
    if (newTick == resultLevels.back().tick) { // update the status
      auto& currentStatus = resultLevels.back();
      currentStatus.event = EventType::Shift;
      currentStatus.opening = newLevel;
    }
    else { // add a new status
      resultLevels.emplace_back(EventType::Shift, newTick, newLevel);
    }
    
  } // while
  
  result.compact();
  resultLevels.shrink_to_fit();
  
  return result;
} // icarus::trigger::TriggerGateData<>::SymmetricCombination()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
bool icarus::trigger::TriggerGateData<TK, TI>::operator ==
  (TriggerGateData const& other) const
{
  return fGateLevel == other.fGateLevel;
} // bool icarus::trigger::TriggerGateData<>::operator==


//------------------------------------------------------------------------------
template <typename TK, typename TI>
bool icarus::trigger::TriggerGateData<TK, TI>::operator !=
  (TriggerGateData const& other) const
{
  return fGateLevel != other.fGateLevel;
} // bool icarus::trigger::TriggerGateData<>::operator==


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findLastStatusFor
  (ClockTick_t tick) -> std::optional<status_iterator>
{
  assert(!fGateLevel.empty());
  // this comes *after* tick
  auto const sbegin = fGateLevel.begin();
  auto iStatus
    = std::upper_bound(sbegin, fGateLevel.end(), tick, CompareTick());
  if (iStatus == sbegin) return {}; // no value
  else return { --iStatus };
} // icarus::trigger::TriggerGateData<>::findLastStatusFor()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findLastStatusFor
  (ClockTick_t tick) const -> std::optional<status_const_iterator>
{
  assert(!fGateLevel.empty());
  // this comes *after* tick
  auto const sbegin = fGateLevel.begin();
  auto iStatus
    = std::upper_bound(sbegin, fGateLevel.end(), tick, CompareTick());
  if (iStatus == sbegin) return {}; // no value
  else return { --iStatus };
} // icarus::trigger::TriggerGateData<>::findLastStatusFor() const


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findLastStatusForTickOrThrow
  (ClockTick_t tick) -> status_iterator
{
  auto const iStatus = findLastStatusFor(tick); // status may be before the tick
  if (iStatus) return iStatus.value();
  // this should be not even possible
  throw std::runtime_error(
    "icarus::trigger::TriggerGateData: requested time " + util::to_string(tick)
    + " is before the gate channel was created (at "
    + util::to_string(fGateLevel.front().tick)
    + ")"
    );
} // icarus::trigger::TriggerGateData<>::findLastStatusForTickOrThrow()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findLastStatusForTickOrThrow
  (ClockTick_t tick) const -> status_const_iterator
{
  auto const iStatus = findLastStatusFor(tick); // status may be before the tick
  if (iStatus) return iStatus.value();
  // this should be not even possible
  throw std::runtime_error(
    "icarus::trigger::TriggerGateData: requested time " + util::to_string(tick)
    + " is before the gate channel was created (at "
    + util::to_string(fGateLevel.front().tick)
    + ")"
    );
} // icarus::trigger::TriggerGateData<>::findLastStatusForTickOrThrow() const


//------------------------------------------------------------------------------
template <typename TK, typename TI>
template <typename Op>
auto icarus::trigger::TriggerGateData<TK, TI>::findStatus
  (Op op, ClockTick_t start /* = MinTick */, ClockTick_t end /* = MaxTick */)
  const -> status_const_iterator
{
  assert(!fGateLevel.empty()); // by construction we have a opening status at 0#
  
  // this is the status before (or at) `start` tick:
  auto const ppStartStatus = findLastStatusFor(start);
  
  // if there is no status at or before `start`, no worry: we just look later;
  // if the status is before `start`, by agreement we ignore it and start from
  // next
  auto iStatus = ppStartStatus
    ? ((ppStartStatus.value()->tick >= start)
      ? ppStartStatus.value()
      : std::next(ppStartStatus.value())
      )
    : fGateLevel.begin()
    ;
  
  auto const send = fGateLevel.end();
  while (iStatus != send) {
    if (iStatus->tick >= end) break;
    switch (iStatus->event) {
      case EventType::Shift:
      case EventType::Set:
        if (op(*iStatus)) return iStatus;
        [[fallthrough]];
      case EventType::Unknown:
        break;
    } // switch
    ++iStatus;
  } // while
  return send;
  
} // icarus::trigger::TriggerGateData<>::findOpenStatus()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findOpenStatus(
  OpeningCount_t minOpening /* = 1U */,
  ClockTick_t start /* = MinTick */, ClockTick_t end /* = MaxTick */
) const -> status_const_iterator
{
  auto const isOpen
    = [minOpening](Status const& status){ return status.opening >= minOpening; }
    ;
  return findStatus(isOpen, start, end);
  
#if 0
  assert(!fGateLevel.empty()); // by construction we have a opening status at 0#
  
  // this is the status before (or at) `start` tick:
  auto const ppStartStatus = findLastStatusFor(start);
  
  // if there is no status at or before `start`, no worry: we just look later;
  // if the status is before `start`, by agreement we ignore it and start from
  // next
  auto iStatus = ppStartStatus
    ? ((ppStartStatus.value()->tick >= start)
      ? ppStartStatus.value()
      : std::next(ppStartStatus.value())
      )
    : fGateLevel.begin()
    ;
  
  auto const send = fGateLevel.end();
  while (iStatus != send) {
    if (iStatus->tick >= end) break;
    switch (iStatus->event) {
      case EventType::Shift:
      case EventType::Set:
        if (iStatus->opening >= minOpening) return iStatus;
        [[fallthrough]];
      case EventType::Unknown:
        break;
    } // switch
    ++iStatus;
  } // while
  return send;
#endif // 0
} // icarus::trigger::TriggerGateData<>::findOpenStatus()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findCloseStatus(
  OpeningCount_t minOpening /* = 1U */,
  ClockTick_t start /* = MinTick */, ClockTick_t end /* = MaxTick */
) const -> status_const_iterator
{
  auto const isClose
    = [minOpening](Status const& status){ return status.opening < minOpening; };
  return findStatus(isClose, start, end);
} // icarus::trigger::TriggerGateData<>::findCloseStatus()


//------------------------------------------------------------------------------
template <typename TK, typename TI>
auto icarus::trigger::TriggerGateData<TK, TI>::findMaxOpenStatus
  (ClockTick_t start /* = MinTick */, ClockTick_t end /* = MaxTick */) const
  -> status_const_iterator
{
  assert(!fGateLevel.empty()); // by construction we have a opening status at 0#
  
  auto const send = fGateLevel.end();
  
  // this is the status before (or at) `start` tick:
  auto const ppStartStatus = findLastStatusFor(start);
  
  // if there is no status at or before `start`, no worry
  auto iStatus = ppStartStatus? ppStartStatus.value(): fGateLevel.begin();
  assert(iStatus != send);
  
  auto iMaxOpening = iStatus;
  
  while (iStatus != send) {
    if (iStatus->tick >= end) break;
    switch (iStatus->event) {
      case EventType::Shift:
      case EventType::Set:
        if (iStatus->opening > iMaxOpening->opening) iMaxOpening = iStatus;
        [[fallthrough]];
      case EventType::Unknown:
        break;
    } // switch
    ++iStatus;
  } // while
  return iMaxOpening;
  
} // icarus::trigger::TriggerGateData<>::findMaxOpenStatus()


//------------------------------------------------------------------------------
//--- output functions
//------------------------------------------------------------------------------
template <typename TK, typename TI>
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, icarus::trigger::TriggerGateData<TK, TI> const& gate)
{
  assert(!gate.fGateLevel.empty());
  auto const send = gate.fGateLevel.end();
  auto iStatus = gate.findOpenStatus(); // first opening, any level
  out << "gate";
  if (iStatus == send) {
    out << " not open at all";
  }
  else {
    auto current = iStatus->opening;
    out << " opened at " << iStatus->tick;
    while (++iStatus != send) {
      if (iStatus->opening == current) continue;
      if ((current == 0) && (iStatus->opening > current)) {
        out << ", reopened";
        if (iStatus->opening > 1) out << " to " << iStatus->opening;
      }
      else if ((current > 0) && (iStatus->opening == 0)) {
        out << ", closed";
      }
      else {
        out << ", changed to " << iStatus->opening;
      }
      out << " at " << iStatus->tick;
      current = iStatus->opening;
    } // while
  } // if ... else
  return out;
} // icarus::trigger::operator<< (icarus::trigger::TriggerGateData)


//------------------------------------------------------------------------------
template <typename TK, typename TI>
std::ostream& icarus::trigger::operator<< (
  std::ostream& out,
  typename icarus::trigger::TriggerGateData<TK, TI>::Status const& status
) {
  using EventType
    = typename icarus::trigger::TriggerGateData<TK, TI>::EventType;
  out << "at " << status.tick;
  switch (status.event) {
    case EventType::Shift:   out << " shift to"; break;
    case EventType::Set:     out << " set to"; break;
    case EventType::Unknown: out << " (unknown)"; break;
  } // switch
  out << " opening: " << status.opening;
  return out;
} // icarus::trigger::operator<< (icarus::trigger::TriggerGateData<>::Status)


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_TCC

