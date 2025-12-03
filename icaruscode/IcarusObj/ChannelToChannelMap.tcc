/**
 * @file   icaruscode/IcarusObj/ChannelToChannelMap.tcc
 * @brief  Object storing a channel map (template implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/IcarusObj/ChannelToChannelMap.h
 * 
 */

#ifndef ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_TCC
#define ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_TCC

#ifndef ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_H
# error "ChannelToChannelMap.tcc must not be included directly."\
        " #include \"icaruscode/IcarusObj/ChannelToChannelMap.h\" instead."
#endif // ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_H


// C/C++ standard libraries
#include <algorithm> // std::sort(), std::lower_bound()


// -----------------------------------------------------------------------------
template <typename Key, typename Value>
typename icarus::ChannelToChannelMap<Key, Value>::MappedChannels_t const
icarus::ChannelToChannelMap<Key, Value>::NoChannels;


// -----------------------------------------------------------------------------
template <typename Key, typename Value>
auto icarus::ChannelToChannelMap<Key, Value>::operator[]
  (KeyChannel_t channel) const noexcept -> MappedChannels_t const&
{
  pointer item = find(channel);
  return item? item->second: NoChannels;
}


// -----------------------------------------------------------------------------
template <typename Key, typename Value>
auto icarus::ChannelToChannelMap<Key, Value>::at(KeyChannel_t channel) const
  -> MappedChannels_t const&
{
  pointer item = find(channel);
  if (!item)
    throw std::out_of_range{ "Unknown channel: " + std::to_string(channel) };
  return item->second;
}


// -----------------------------------------------------------------------------
template <typename Key, typename Value>
bool icarus::ChannelToChannelMap<Key, Value>::contains
  (KeyChannel_t channel) const noexcept
{
  return find(channel) != nullptr;
}


// -----------------------------------------------------------------------------
template <typename Key, typename Value>
auto icarus::ChannelToChannelMap<Key, Value>::lowestChannel() const -> KeyChannel_t
  { return fMap.front().first; }


  // -----------------------------------------------------------------------------
template <typename Key, typename Value>
auto icarus::ChannelToChannelMap<Key, Value>::highestChannel() const -> KeyChannel_t
  { return fMap.back().first; }


// -----------------------------------------------------------------------------
template <typename Key, typename Value>
auto icarus::ChannelToChannelMap<Key, Value>::addChannel
  (KeyChannel_t channel, MappedChannels_t toChannels) -> ChannelToChannelMap&
{
  std::sort(toChannels.begin(), toChannels.end());
  auto const insPoint
    = std::lower_bound(fMap.begin(), fMap.end(), channel, CompareByChannel{});
  if ((insPoint == fMap.cend()) || (insPoint->first != channel)) // not found
    fMap.emplace(insPoint, channel, std::move(toChannels));
  else insPoint->second = std::move(toChannels); // exists: replace
  return *this;
} // icarus::ChannelToChannelMap<Key, Value>::addChannel()


// -----------------------------------------------------------------------------
template <typename Key, typename Value>
auto icarus::ChannelToChannelMap<Key, Value>::find
  (KeyChannel_t channel) const noexcept -> const_pointer
{
  auto const it
    = std::lower_bound(fMap.cbegin(), fMap.cend(), channel, CompareByChannel{});
  return (it == fMap.cend() || it->first != channel)? nullptr: &*it;
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_TCC
