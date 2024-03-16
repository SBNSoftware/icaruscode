/**
 * @file   icaruscode/PMT/Trigger/Algorithms/LVDSbitMaps.cxx
 * @brief  Utilities for mapping CAEN V1730B LVDS channels to PMT channels.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 10, 2024
 * @see    icaruscode/PMT/Trigger/Algorithms/LVDSbitMaps.h
 * 
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/LVDSbitMaps.h"

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapDataTypes.h"

// framework libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <ostream>
#include <algorithm> // std::sort()
#include <iomanip> // std::dec, std::hex
#include <cstdlib> // std::div()
#include <cassert>


// -----------------------------------------------------------------------------
// ---  icarus::trigger::LVDSbitID implementation
// -----------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, LVDSbitID const& bit)
{
  static constexpr char Side[] = "EW";
  auto const fillch = out.fill(); // width is reset to 0, fill character is not
  out << (bit.hasCryostat()? Side[bit.cryostat]: '?')
    << '/' << (bit.hasPMTwall()? Side[bit.PMTwall]: '?')
    << '/' << std::setw(2) << std::setfill('0') << bit.statusBit
    << std::setfill(fillch);
  return out;
}


// -----------------------------------------------------------------------------
auto icarus::trigger::LVDSbitID::fromIndex(Index_t index) -> LVDSbitID
{
  auto const [ cryostat, wb ]
    = std::div(index, NPMTwallsPerCryostat * NStatusWordBits);
  auto const [ wall, bit ] = std::div(wb, NStatusWordBits);
  return { (std::size_t) cryostat, (std::size_t) wall, (StatusBit_t) bit };
} // icarus::trigger::LVDSbitID::fromIndex()


// -----------------------------------------------------------------------------
// ---  icarus::trigger::LVDSHWbitID implementation
// -----------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, LVDSHWbitID const& bit)
{
  if (!bit) return out << "<invl>";
  auto const fillch = out.fill(); // width is reset to 0, fill character is not
  out << bit.cryostat << '/' << bit.connector << '/'
    << std::setw(2) << std::setfill('0') << bit.bit << std::setfill(fillch);
  return out;
}


// -----------------------------------------------------------------------------
auto icarus::trigger::LVDSHWbitID::fromIndex(Index_t index) -> LVDSHWbitID {
  auto const [ cryostat, cb ] = std::div(index, NConnectorsPerCryostat * 32);
  auto const [ connector, bit ] = std::div(cb, 32);
  return {
    (unsigned short int) cryostat,
    (Connector_t) connector,
    (ConnectorBit_t) bit
    };
} // icarus::trigger::LVDSHWbitID::fromIndex()


// -----------------------------------------------------------------------------
// ---  icarus::trigger::PMTpairBitMap implementation
// -----------------------------------------------------------------------------
auto icarus::trigger::PMTpairBitMap::bitSource(LVDSbitID bit) const
  -> LVDSinfo_t
{
  return fLVDSinfoMap.at(bit.index());
}


// -----------------------------------------------------------------------------
std::vector<std::vector<raw::Channel_t>>
icarus::trigger::PMTpairBitMap::channelPairs() const {
  
  std::vector<std::vector<raw::Channel_t>> pairs;
  for (LVDSinfo_t const& info: fLVDSinfoMap) {
    if (!info.channels.empty()) pairs.push_back(info.channels);
  }
  
  return pairs;
}

// -----------------------------------------------------------------------------
auto icarus::trigger::PMTpairBitMap::buildLVDSsourceMap
  (icarusDB::IICARUSChannelMapProvider const& channelMap) const
  -> std::vector<LVDSinfo_t>
{
  //
  // get the LVDS connector/bit information channel by channel
  //
  std::vector<LVDSchannel_t> PMTchannelDB;
  
  mf::LogTrace debugLog{ "PMTpairBitMap" };
  debugLog << "PMTpairBitMap::buildLVDSsourceMap() collecting channel info:";
  
  for (unsigned int const fragmentID: util::counter(0x2000, 0x2018)) {
    if (!channelMap.hasPMTDigitizerID(fragmentID)) {
      throw Exception{ Errors::NoPMTfragment }
        << "PMT fragment ID 0x" << std::hex << fragmentID << std::dec
        << " not available in the channel mapping database.\n";
    }
    
    for (icarusDB::PMTChannelInfo_t const& info
      : channelMap.getPMTchannelInfo(fragmentID)
    ) {
      unsigned short int const cryostat = info.channelID / NChannelsPerCryostat;
      
      if (!info.hasLVDSinfo()) {
        throw Exception{ Errors::NoPairBitInformation }
          << "PMT fragment ID 0x" << std::hex << fragmentID << std::dec
          << " has no LVDS information in the channel mapping database.\n";
      }
      
      debugLog << "\n[#" << PMTchannelDB.size() << "] CH=" << info.channelID
        << " " << cryostat << "/" << info.LVDSconnector << "/" << info.LVDSbit;
      
      PMTchannelDB.push_back(LVDSchannel_t{
          info.channelID        // .channel
        , LVDSHWbitID{          // .LVDSsource
            cryostat              // .cryostat
          , info.LVDSconnector    // .connector
          , info.LVDSbit          // .bit
          }
        });
    } // for channel in fragment
  } // for fragment ID
  
  assert(PMTchannelDB.size() == NChannels); // insider knowledge...
  
  std::sort(PMTchannelDB.begin(), PMTchannelDB.end());
  
  //
  // fill the logic LVDS bit positions with their source information
  //
  
  /*
   * The logic of the filling is a bit complicated. Class documentation explains
   * it and shows it to be:
   * 
   * east cryostat, east wall: `0x00aAbBcC'00dDeEfF`, channel `0` the most
   * significant bit in `a` (bit #55) and channel `45` the one in `d` (bit 23);
   * east cryostat, west wall: `0x00aAbBcC'00dDeEfF`, channel `90` the most
   * significant bit in `a` (bit #55) and channel `135` the one in `d` (bit 23);
   * west cryostat, east wall: `0x00aAbBcC'00dDeEfF`, channel `180` the most
   * significant bit in `a` (bit #55) and channel `225` the one in `d` (bit 23);
   * west cryostat, west wall: `0x00aAbBcC'00dDeEfF`, channel `270` the most
   * significant bit in `a` (bit #55) and channel `315` the one in `d` (bit 23).
   * 
   * How we get there:
   *  * we sort the filling with an ordered progression of channels, lowest
   *    first, and because most LVDS bit are associated to two channels, we
   *    track the channels we already covered when processing the other PMT
   *    channel of the pair, and skip them when met those LVDS bit again
   *  * the destination bit is also an ordered progression within the word,
   *    but inverted: note in fact that if we inverted the position of the bits
   *    in each 32-bit word
   */
  std::vector<LVDSinfo_t> sourceMap{ LVDSbitID::NIndices() };
  std::vector<LVDSbitID> HWtoLogicID{ LVDSHWbitID::NIndices() };
  
  unsigned int nBits = 0;
  for (LVDSchannel_t const& info: PMTchannelDB) { // in PMT channel order
    
    LVDSHWbitID::Index_t const sourceIndex = info.LVDSsource.index();
    
    LVDSbitID logicBit = HWtoLogicID[sourceIndex];
    
    LVDSinfo_t* mapInfo = logicBit? &(sourceMap[logicBit.index()]): nullptr;
    if (!mapInfo) { // new bit: set the source
      // which logic bit to set:
      logicBit.cryostat = info.LVDSsource.cryostat;
      assert(logicBit.cryostat == info.channel / NChannelsPerCryostat);
      logicBit.PMTwall
        = info.channel / (NChannelsPerCryostat/LVDSbitID::NPMTwallsPerCryostat)
        % LVDSbitID::NCryostats;
      logicBit.statusBit = 47 - (nBits++ % 48); // 0 => 47, 1 => 46,..., 47 => 0
      if (logicBit.statusBit >= 24)
        logicBit.statusBit += 8; // 0 => 0,..., 23 => 23, 24 => 32, 25 => 33,...
      
      assert(logicBit);
      HWtoLogicID[sourceIndex] = logicBit;
      mapInfo = &(sourceMap[logicBit.index()]);
      mapInfo->source = info.LVDSsource;
    }
    assert(mapInfo);
    assert(mapInfo->source == info.LVDSsource);
    
    mapInfo->channels.push_back(info.channel);
    // check that they are sorted (the input list should be sorted already..._
    assert((mapInfo->channels.size() == 1)
      || (
        mapInfo->channels[mapInfo->channels.size()-2]
        < mapInfo->channels[mapInfo->channels.size()-1]
      ));
    
  } // for PMTchannelDB
  
  
  // --- BEGIN DEBUG -----
  debugLog << "\nFinal map:";
  for (auto const& [ bit, info ]: util::enumerate(sourceMap)) {
    debugLog << "\n[ " << std::setw(3) << bit << "] logic="
      << LVDSbitID::fromIndex(bit) << " => HW=";
    if (info.source) {
      debugLog << info.source << " (ix="
        << std::setw(3) << info.source.index() << ") "
        << info.channels.size() << " channels {";
      for (raw::Channel_t const channel: info.channels)
        debugLog << " " << std::setw(3) << channel;
      debugLog << " }";
    }
    else {
      debugLog << "(invalid)";
    }
  }
  // --- END DEBUG -------
  
  assert(nBits == 192);
  return sourceMap;
  
} // icarus::trigger::PMTpairBitMap::buildLVDSsourceMap()


// -----------------------------------------------------------------------------
std::string icarus::trigger::PMTpairBitMap::errorMessage(Errors code) {
  using namespace std::string_literals;
  switch (code) {
    case Errors::NoError:
      return "no error"s;
    case Errors::Error:
      return "(unspecified) error"s;
    case Errors::NoDatabase:
      return "channel mapping database not configured"s;
    case Errors::NoPMTfragment:
      return "channel mapping database does not have PMT fragment information"s;
    case Errors::NoPairBitInformation:
      return "channel mapping database does not have PMT pair bit information"s;
    case Errors::Unknown:
      return "unknown error"s;
    default:
      return "unexpected error code "s + std::to_string(static_cast<int>(code));
  } // switch
} // icarus::trigger::PMTpairBitMap::errorMessage()


// -----------------------------------------------------------------------------
