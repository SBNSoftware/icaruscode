/**
 * @file   icaruscode/Geometry/details/ChannelToWireMap.h
 * @brief  Channel-to-wire mapping data structure.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 16, 2019
 * @see    `icaruscode/Geometry/details/ChannelToWireMap.cxx`
 */

#ifndef ICARUSCODE_GEOMETRY_DETAILS_CHANNELTOWIREMAP_H
#define ICARUSCODE_GEOMETRY_DETAILS_CHANNELTOWIREMAP_H

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// C/C++ standard libraries
#include <vector>
#include <cassert>
#include <utility> // std::pair<>


// -----------------------------------------------------------------------------
// --- icarus::details::ChannelToWireMap
// -----------------------------------------------------------------------------
namespace icarus::details { class ChannelToWireMap; }
/**
 * @brief Mapping of ROP channels into wires.
 * 
 * For each readout plane, we store its first channel (the last one can be
 * found from the next ROP), and for all wire planes after the first one,
 * the first channel that is shared (with the previous plane, that is) and
 * the first channel on that plane that is not shared.
 */
class icarus::details::ChannelToWireMap {
  
    public:
  struct ChannelsInROPStruct {
    
    /// First channel in ROP.
    raw::ChannelID_t firstChannel = raw::InvalidChannelID;
    unsigned int nChannels = 0U;
    readout::ROPID ropid; ///< ID of the ROP we cover.
    
    /// Compares with a channel ID (std::less<>).
    struct Compare {
      template <typename A, typename B>
      constexpr bool operator() (A const& a, B const& b) const
        { return cmp(a, b); }
      
      static constexpr raw::ChannelID_t channelOf(raw::ChannelID_t channel)
        { return channel; }
      static constexpr raw::ChannelID_t channelOf
        (ChannelsInROPStruct const& data)
        { return data.firstChannel; }
      
      template <typename A, typename B>
      static constexpr bool cmp(A const& a, B const& b)
        { return channelOf(a) < channelOf(b); }
      
    }; // struct Compare
    
    ChannelsInROPStruct() = default;
    ChannelsInROPStruct(
      raw::ChannelID_t firstChannel, unsigned int nChannels,
      readout::ROPID const& ropid
      )
      : firstChannel(firstChannel)
      , nChannels(nChannels)
      , ropid(ropid)
      {}
    
    /// Returns if `channel` is not in a lower ROP than this one.
    constexpr bool isChannelAbove(raw::ChannelID_t channel) const
      { return (channel >= firstChannel); }
    
    /// Strict ordering according to the first channel in ROP.
    constexpr bool operator< (ChannelsInROPStruct const& other) const
      { return Compare{}(*this, other); }
    
  }; // struct ChannelsInROPStruct
  
  /// Returns data of the ROP including `channel`, `nullptr` if none.
  ChannelsInROPStruct const* find(raw::ChannelID_t channel) const;
  
  /// Returns data of the ROP `ropid`, `nullptr` if none.
  ChannelsInROPStruct const* find(readout::ROPID const& ropid) const;
  
  /// Returns the ID of the first invalid channel (the last channel, plus 1).
  raw::ChannelID_t endChannel() const { return fEndChannel; }
  
  /// Returns the number of mapped channels.
  unsigned int nChannels() const
    { return endChannel() - raw::ChannelID_t{0}; }
  
  /// Adds the next ROP with its first channel ID
  /// (must be larger than previous).
  void addROP(
    readout::ROPID const& rid,
    raw::ChannelID_t firstROPchannel, unsigned int nChannels
    )
    {
      assert(
        fROPfirstChannel.empty()
        || (firstROPchannel > fROPfirstChannel.back().firstChannel)
        );
      fROPfirstChannel.emplace_back(firstROPchannel, nChannels, rid);
    }
  
  /// Sets the ID of the channels after the last valid one.
  void setEndChannel(raw::ChannelID_t channel) { fEndChannel = channel; }
  
  /// Resets the data of the map to like just constructed.
  void clear();
  
    private:
  
  /// Collection of channel ROP channel information, sorted by first channel.
  std::vector<ChannelsInROPStruct> fROPfirstChannel;
  
  raw::ChannelID_t fEndChannel = 0; ///< ID of the first invalid channel.
  
}; // class icarus::details::ChannelToWireMap


#endif // ICARUSCODE_GEOMETRY_DETAILS_CHANNELTOWIREMAP_H
