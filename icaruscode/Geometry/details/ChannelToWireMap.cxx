/**
 * @file   icaruscode/Geometry/details/ChannelToWireMap.cxx
 * @brief  Channel-to-wire mapping data structure (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 16, 2019
 * @see    `icaruscode/Geometry/details/ChannelToWireMap.h`
 */

// library header
#include "icaruscode/Geometry/details/ChannelToWireMap.h"

// C/C++ libraries
#include <iterator> // std::prev()
#include <algorithm> // std::upper_bound(), std::find_if()


//------------------------------------------------------------------------------
auto icarus::details::ChannelToWireMap::find
  (raw::ChannelID_t channel) const -> ChannelsInROPStruct const*
{
  // find the first element whose channel is _not_ included
  auto const dbegin = fROPfirstChannel.begin(), dend = fROPfirstChannel.end();
  auto const iNextData = std::upper_bound
    (dbegin, dend, channel, ChannelsInROPStruct::Compare{});
  assert(iNextData != dbegin);
  
  return ((iNextData == dend) && (channel >= endChannel()))
    ? nullptr: &*std::prev(iNextData);
  
} // icarus::details::ChannelToWireMap::find(raw::ChannelID_t)

//------------------------------------------------------------------------------
auto icarus::details::ChannelToWireMap::find
  (readout::ROPID const& ropid) const -> ChannelsInROPStruct const*
{
  // find the first element whose channel is _not_ included
  auto const dend = fROPfirstChannel.end();
  auto const iData = std::find_if(
    fROPfirstChannel.begin(), dend,
    [&ropid](ChannelsInROPStruct const& data){ return data.ropid == ropid; }
    );
  return (iData == dend)? nullptr: &*iData;
} // icarus::details::ChannelToWireMap::find(readout::ROPID)


// -----------------------------------------------------------------------------
void icarus::details::ChannelToWireMap::clear() {
  fROPfirstChannel.clear();
  fEndChannel = raw::ChannelID_t{ 0 };
} // icarus::details::ChannelToWireMap::clear()


// -----------------------------------------------------------------------------
