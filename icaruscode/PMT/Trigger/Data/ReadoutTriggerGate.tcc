/**
 * @file   icaruscode/PMT/Trigger/Data/ReadoutTriggerGate.tcc
 * @brief  A trigger gate data object associated to one or more channels.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 16, 2019
 * @see    icaruscode/PMT/Trigger/Data/ReadoutTriggerGate.h
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_TCC
#define ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_TCC

// C/C++ standard libraries
#include <ostream>
#include <string> // std::to_string()
#include <algorithm> // std::unique(), std::lower_bound(), std::sort()...
#include <iterator> // std::next()
#include <type_traits> // std::is_base_of_v, std::enable_if_t, std::decay_t


// make "sure" this header is not included directly
#ifndef ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_H
# error "ReadoutTriggerGate.tcc must not be directly included!"\
        " #include \"icaruscode/PMT/Trigger/Data/ReadoutTriggerGate.h\" instead."
#endif // !ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_H


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::ReadoutTriggerGate
  (std::initializer_list<ChannelID_t> channels)
  : ReadoutTriggerGate()
{
  addChannels(channels);
} // icarus::trigger::ReadoutTriggerGate<>::ReadoutTriggerGate(ChannelID_t)


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
decltype(auto)
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::channels
  () const
{
  return fChannels;
} // icarus::trigger::ReadoutTriggerGate<>::channels()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::channel
  () const -> ChannelID_t 
{
  if (!hasChannels())
    { throw NoChannelError("No channel associated to trigger data."); }
  if (nChannels() > 1U) {
    throw MoreThanOneChannelError(
      std::to_string(nChannels()) + " channels associated to trigger data.",
      nChannels()
      );
  }
  return fChannels.front();
} // icarus::trigger::ReadoutTriggerGate<>::channel()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::addChannel
  (ChannelID_t const channel) -> This_t&
{
  auto const iNearest
    = std::lower_bound(fChannels.begin(), fChannels.end(), channel);
  if ((iNearest == fChannels.end()) || (*iNearest != channel))
    fChannels.insert(iNearest, channel);
  return *this;
} // icarus::trigger::ReadoutTriggerGate<>::addChannel()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
void
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::associateChannels
  (std::initializer_list<ChannelID_t> const& moreChannels)
{
  mergeChannelsInto(fChannels, moreChannels.begin(), moreChannels.end());
} // icarus::trigger::ReadoutTriggerGate<>::associateChannels()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
void
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::associateChannels
  (ChannelList_t const& moreChannels)
{
  mergeChannelsInto(fChannels, moreChannels.begin(), moreChannels.end());
} // icarus::trigger::ReadoutTriggerGate<>::associateChannels()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::normalizeChannels
  (ChannelList_t& channels) -> ChannelList_t&
{
  std::sort(channels.begin(), channels.end());
  return normalizeSortedChannels(channels);
} // icarus::trigger::ReadoutTriggerGate<>::normalizeChannels()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::normalizeChannels
  (ChannelList_t&& channels) -> ChannelList_t
{
  std::sort(channels.begin(), channels.end());
  normalizeSortedChannels(channels);
  return std::move(channels);
} // icarus::trigger::ReadoutTriggerGate<>::normalizeChannels()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::normalizeSortedChannels
  (ChannelList_t& channels) -> ChannelList_t&
{
  // move non-unique elements to the end, and delete them afterwards
  channels.erase
    (std::unique(channels.begin(), channels.end()), channels.end());
  channels.shrink_to_fit();
  return channels;
} // icarus::trigger::ReadoutTriggerGate<>::normalizeChannels()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
template <typename BIter, typename EIter>
auto
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::mergeSortedChannelsInto
  (ChannelList_t& channels, BIter b, EIter e) -> ChannelList_t&
{
  // non-optimized implementation:
  std::size_t nOrigChannels = channels.size();
  channels.insert(channels.end(), b, e);
  std::inplace_merge(
    channels.begin(), std::next(channels.begin(), nOrigChannels), channels.end()
    );
  return normalizeSortedChannels(channels);
} // icarus::trigger::ReadoutTriggerGate<>::mergeSortedChannelsInto()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
template <typename BIter, typename EIter>
auto
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::mergeChannelsInto
  (ChannelList_t& channels, BIter b, EIter e) -> ChannelList_t&
{
  // non-optimized implementation:
  channels.insert(channels.end(), b, e);
  return normalizeChannels(channels);
} // icarus::trigger::ReadoutTriggerGate<>::mergeChannelsInto()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::mergeChannels
  (ChannelList_t const& a, ChannelList_t const& b) -> ChannelList_t
{
  // non-optimized implementation:
  ChannelID_t merged { a };
  return mergeChannelsInto(a, b.begin(), b.end());
} // icarus::trigger::ReadoutTriggerGate<>::mergeChannels()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Min
  (This_t const& other) -> This_t& 
{
  GateData_t::Min(other);
  associateChannels(other.channels());
  return *this;
} // icarus::trigger::ReadoutTriggerGate<>::Min()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Max
  (This_t const& other) -> This_t& 
{
  GateData_t::Max(other);
  associateChannels(other.channels());
  return *this;
} // icarus::trigger::ReadoutTriggerGate<>::Max()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Sum
  (This_t const& other) -> This_t& 
{
  GateData_t::Sum(other);
  associateChannels(other.channels());
  return *this;
} // icarus::trigger::ReadoutTriggerGate<>::Sum()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Mul
  (This_t const& other) -> This_t& 
{
  GateData_t::Mul(other);
  associateChannels(other.channels());
  return *this;
} // icarus::trigger::ReadoutTriggerGate<>::Mul()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Min
  (This_t const& a, This_t const& b) -> This_t
  { auto combination { a }; combination.Min(b); return combination; }


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Max
  (This_t const& a, This_t const& b) -> This_t
  { auto combination { a }; combination.Max(b); return combination; }


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Sum
  (This_t const& a, This_t const& b) -> This_t
  { auto combination { a }; combination.Sum(b); return combination; }


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
auto icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::Mul
  (This_t const& a, This_t const& b) -> This_t
  { auto combination { a }; combination.Mul(b); return combination; }


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
template <typename Op>
auto
icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::SymmetricCombination(
  Op&& op, This_t const& a, This_t const& b,
  ClockTicks_t aDelay /* = ClockTicks_t{ 0 } */,
  ClockTicks_t bDelay /* = ClockTicks_t{ 0 } */
  ) -> This_t
{
  return {
    GateData_t::SymmetricCombination
      (std::forward<Op>(op), a, b, aDelay, bDelay),
    mergeChannels(a.channels(), b.channels())
  };
  
} // icarus::trigger::ReadoutTriggerGate<>::SymmetricCombination()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
bool icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::operator==
  (ReadoutTriggerGate const& other) const
{
  return
    (gateLevels() == other.gateLevels()) && (channels() == other.channels());
} // icarus::trigger::ReadoutTriggerGate<>::operator==()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
bool icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>::operator!=
  (ReadoutTriggerGate const& other) const
{
  return
    (gateLevels() != other.gateLevels()) || (channels() != other.channels());
} // icarus::trigger::ReadoutTriggerGate<>::operator!=()


//------------------------------------------------------------------------------
template <typename Tick, typename TickInterval, typename ChannelIDType>
std::ostream& icarus::trigger::operator<< (
  std::ostream& out,
  icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType> const& gate
) {
  out << "[";
  if (gate.hasChannel()) {
    out << "channel: " << gate.channel();
  }
  else if (gate.hasChannels()) {
    out << gate.nChannels() << " channels:";
    for (auto channel: gate.channels()) out << " " << channel;
  }
  else out << "no channel";
  out << "] " << gate.gateLevels();
  return out;
} // icarus::trigger::operator<< (icarus::trigger::ReadoutTriggerGate)


//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  namespace details {
    
    template <typename Gate, typename = void>
    struct isReadoutTriggerGateImpl: std::false_type {};
    
    
    template <typename Gate>
    struct isReadoutTriggerGateImpl
      <Gate, std::enable_if_t<std::is_base_of_v<ReadoutTriggerGateTag, Gate>>>
      : std::true_type {};
    
  } // namespace details
  
  template <typename Gate>
  struct isReadoutTriggerGate
    : details::isReadoutTriggerGateImpl<std::decay_t<Gate>>
  {};
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_H
