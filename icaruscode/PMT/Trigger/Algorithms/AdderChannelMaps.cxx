/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.cxx
 * @brief  Utilities for assigning IDs to adder channels (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 6, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.h
 * 
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.h"

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapDataTypes.h"
#include "icarusalg/Utilities/IntegerRanges.h"

// LArSoft/framework libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/values.h"

// C/C++ standard libraries
#include <algorithm> // std::sort()
#include <cassert>
#include <ios> // std::dec, std::hex
#include <ostream>
#include <string>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::trigger::details {
  
  
  struct AdderInfo_t {
    
    std::vector<icarusDB::PMTChannelInfo_t> PMTs;
    
    bool empty() const { return PMTs.empty(); }
    icarusDB::PMTChannelInfo_t const& front() const
      { assert(!empty()); return PMTs.at(0); }
    
    std::string const& digitizerLabel() const { return front().digitizerLabel; }
    
    void sortPMTs()
      {
        std::sort(
          PMTs.begin(), PMTs.end(),
          [](icarusDB::PMTChannelInfo_t const& A, icarusDB::PMTChannelInfo_t const& B)
            { return A.channelID < B.channelID; }
          );
      }
    
    bool operator< (AdderInfo_t const& other) const
      {
        assert(!empty() && !other.empty());
        return front().channelID < other.front().channelID;
      }
    
  }; // AdderInfo_t
  
  
  /// Adder channel ID information. Handles both new and old standard.
  struct AdderChannelIDinfo {
    
    AdderChannelID channel{}; ///< The standard adder channel ID.
    
    short int inRack = -1; ///< Whether its crate is at top (0) or bottom (1) of the rack.
    short int slot = 0; ///< Which position in the rack the board is (`1`, `2`, `3`).
    
    /// Returns the cryostat this channel ID belongs to (`0` for east, `1` for west).
    constexpr int cryostat() const
      { return PMTwall() >> 1; }
    
    /// Returns the side of cryostat this channel ID belongs to (`0` for east, `1` for west).
    constexpr int TPC() const
      { return PMTwall() & 1; }
    
    /// Returns the PMT wall this channel ID belongs to: 0 ("EE") to 3 ("WW").
    constexpr int PMTwall() const
      { return (static_cast<raw::Channel_t>(channel) & 0x0FFF) / 6; }
    
    /// Returns the position in its rack of the crate serving this channel
    /// (`0` for bottom, `1` for top).
    constexpr int rackCrate() const { return inRack; }
    
    /// Position of the slot of readout board serving this channel (1 to 3).
    constexpr int crateBoard() const { return slot; }
    
    /// Sequence number of this adder in its PMT wall (`0` to `5`).
    constexpr int sequenceInWall() const
      { return (static_cast<raw::Channel_t>(channel) & 0x0FFF) % 6; }
    
    /// Default constructor: invalid channel.
    constexpr AdderChannelIDinfo() = default;
    
    /// Constructor: specifies ID and optionally location in crate.
    constexpr AdderChannelIDinfo
      (AdderChannelID channel, short int endSide = -1, short int slot = 0)
      : channel{ channel }, inRack{ endSide }, slot{ slot }
      {}
    
    
    /// Returns this channel as the "old standard ID".
    constexpr AdderChannelID oldStandardID() const
      { return 0x8000 + (((PMTwall() * 2 + rackCrate()) & 0xF) << 4) + crateBoard(); }
    
  }; // AdderChannelIDinfo


} // local namespace


// -----------------------------------------------------------------------------
// ---  icarus::trigger::AdderChannelMapBuilder
// -----------------------------------------------------------------------------
icarus::trigger::AdderChannelMapBuilder::AdderChannelMapBuilder
  (std::string logCategory /* = "AdderChannelMapBuilder" */)
  : icarus::ns::util::mfLoggingClass{ std::move(logCategory) }
  {}


// -----------------------------------------------------------------------------
icarus::trigger::AdderChannelMap icarus::trigger::AdderChannelMapBuilder::build(
  icarusDB::IICARUSChannelMapProvider const& PMTchannelMap,
  [[maybe_unused]] int run /* = 0 */
) const {
  
  using details::AdderInfo_t;
  
  //
  // collect the information, grouped by adder bit
  //
  
  // { cryostat, connector, bit }
  struct AdderBit: std::tuple<short int, unsigned short int, unsigned short int> {
    short int cryostat() const { return std::get<0>(*this); }
    unsigned short int connector() const { return std::get<1>(*this); }
    unsigned short int bit() const { return std::get<2>(*this); }
  };
  
  std::map<AdderBit, AdderInfo_t> AdderChannelsByBit;
  
  // collect information by the PMT fragment ID
  for (unsigned int const fragmentID: PMTchannelMap.PMTfragmentIDs()) {

    for (icarusDB::PMTChannelInfo_t const& channelInfo
      : PMTchannelMap.getPMTchannelInfo(fragmentID)
    ) {
      
      AdderBit const adderBit {{
        channelInfo.digitizerInfo.cryostat,
        channelInfo.adderConnector,
        channelInfo.adderBit
      }};
      
      AdderInfo_t& adderInfo = AdderChannelsByBit[adderBit];
      
      if (!adderInfo.empty()
        && (adderInfo.digitizerLabel() != channelInfo.digitizerLabel)
      ) {
        throw cet::exception{ "buildAdderChannelMap" }
          << "Adder on C:" << adderBit.cryostat()
          << " connector " << adderBit.connector() << " bit #" << adderBit.bit()
          << " mixes PMTs on digitizer " << adderInfo.digitizerLabel()
          << " (with PMT channels like " << adderInfo.front().channelID
          << ") to ones  on digitizer " << channelInfo.digitizerLabel
          << " (with PMT channels like " << channelInfo.channelID << ")!\n";
      } // if inconsistent digitizer
      
      adderInfo.PMTs.push_back(channelInfo);
      mfLogTrace() << "Adder connector " << channelInfo.adderConnector
        << " bit#" << channelInfo.adderBit << " gained PMT channel #"
        << channelInfo.channelID;
      
    } // for channels in fragment
  } // for PMT fragment ID

  
  //
  // sort the groups by PMT channel
  //
  
  // compact the information (loosing the adder bit key)
  std::vector<AdderInfo_t> AdderChannels;
  AdderChannels.reserve(AdderChannelsByBit.size());
  for (AdderInfo_t& adderInfo: util::values(AdderChannelsByBit)) {
    adderInfo.sortPMTs();
    AdderChannels.push_back(std::move(adderInfo));
  }
  
  // sort by the value of the first (lowest) PMT channel in the adder channel:
  std::sort(AdderChannels.begin(), AdderChannels.end());
  
  //
  // assign adder channels an ID
  //
  AdderChannelMap channelMap;
  
  // number of channels already assigned, per wall index
  
  for (auto const& [ adderIndex, adderInfo ]: util::enumerate(AdderChannels)) {
    
    // assume that the adder board lies beside the PMT readout boards it serves:
    // learn the digitizer information from the first PMT channel
    icarusDB::PMTChannelInfo_t::DigitizerInfo_t const& digitizerInfo
      = adderInfo.front().digitizerInfo;
    
    AdderChannelID const channel
      = AdderChannelID::BaseAdderChannelID + adderIndex;
    
    details::AdderChannelIDinfo const channelInfo{
      /* channel = */ static_cast<raw::Channel_t>(adderIndex),
      /* endSize = */ digitizerInfo.endSide,
      /* slot =    */ static_cast<short int>(digitizerInfo.slot + 1)
      };
    
    assert(channelInfo.channel == channel); // they should be consistent
      
    AdderChannelInfo_t::PMTchannelList_t PMTs;
    PMTs.reserve(adderInfo.PMTs.size());
    for (icarusDB::PMTChannelInfo_t const& channelInfo: adderInfo.PMTs)
      PMTs.push_back(channelInfo.channelID);
    
    AdderChannelInfo_t adderChannelInfo {
      /* .channel = */ channel,
      /* .PMTs    = */ std::move(PMTs)
      };
    
    mfLogTrace() << "Adder channel " << adderChannelInfo.channel
      << " (" << channelInfo.channel << "/" << channelInfo.oldStandardID() << ")"
      << " assigned " << adderChannelInfo.PMTs.size()
      << " PMT channels: "
      << icarus::IntegerRanges{ adderChannelInfo.PMTs.cbegin(), adderChannelInfo.PMTs.cend() };

    channelMap.emplace(channel, std::move(adderChannelInfo));
  } // for
  
  return channelMap;
  
} // icarus::trigger::buildAdderChannelMap()


// -----------------------------------------------------------------------------
