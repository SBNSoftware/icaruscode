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
#include <algorithm> // std::sort()
#include <iomanip> // std::dec, std::hex
#include <cstdlib> // std::div()
#include <cassert>


// -----------------------------------------------------------------------------
// ---  icarus::trigger::PMTpairBitID implementation
// -----------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, PMTpairBitID const& bit)
{
  return out << static_cast<PMTpairBitID::BaseID_t const&>(bit);
}


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
// ---  icarus::trigger::LVDSbitMaps implementation
// -----------------------------------------------------------------------------
bool icarus::trigger::LVDSbitMaps::hasMap(Map map) const {
  switch (map) {
    case Map::PMTpairs: return !fLVDSinfoMap.empty();
    case Map::Adders:   return !fAdderInfoMap.empty();
    default:
      throw Exception{ Errors::Error }
        << "Logic error: hasMap() does not support map type #"
        << static_cast<int>(map) << "!\n";
  } // switch
} // icarus::trigger::LVDSbitMaps::hasMap()


// -----------------------------------------------------------------------------
auto icarus::trigger::LVDSbitMaps::bitSource(PMTpairBitID const& bit) const
  -> LVDSinfo_t
{
  return fLVDSinfoMap.at(bit.index());
}


// -----------------------------------------------------------------------------
auto icarus::trigger::LVDSbitMaps::bitSource(AdderBitID const& bit) const
  -> AdderInfo_t
{
  return fAdderInfoMap.at(bit.index());
}


// -----------------------------------------------------------------------------
std::vector<std::vector<raw::Channel_t>>
icarus::trigger::LVDSbitMaps::channelPairs() const {
  
  std::vector<std::vector<raw::Channel_t>> pairs;
  for (LVDSinfo_t const& info: fLVDSinfoMap) {
    if (!info.channels.empty()) pairs.push_back(info.channels);
  }
  
  return pairs;
}


// -----------------------------------------------------------------------------
void icarus::trigger::LVDSbitMaps::buildAllMaps
  (icarusDB::IICARUSChannelMapProvider const& channelMap)
{
  buildMaps(channelMap, { Map::PMTpairs, Map::Adders });
}


// -----------------------------------------------------------------------------
void icarus::trigger::LVDSbitMaps::buildMaps(
  icarusDB::IICARUSChannelMapProvider const& channelMap,
  std::initializer_list<Map> maps
) {

  fPMTchannelInfo = buildChannelInfo(channelMap);

  for (Map map: maps) {
    switch (map) {
      case Map::PMTpairs:
        fLVDSinfoMap = buildPMTpairSourceMap(channelMap);
        break;
      case Map::Adders:
        fAdderInfoMap = buildAdderSourceMap(channelMap);
        break;
    } // switch
  } // for
  
} // icarus::trigger::LVDSbitMaps::buildMaps()


// -----------------------------------------------------------------------------
auto icarus::trigger::LVDSbitMaps::buildChannelInfo
  (icarusDB::IICARUSChannelMapProvider const& channelMap) const
  -> std::vector<LVDSchannel_t>
{
  std::vector<LVDSchannel_t> PMTchannelDB;
  
  mf::LogTrace debugLog{ "LVDSbitMaps" };
  debugLog << "LVDSbitMaps::buildChannelInfo() collecting channel info:";
  
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
      
      debugLog << "\n[#" << PMTchannelDB.size() << "] CH=" << info.channelID
        << " " << cryostat << "/" << info.LVDSconnector << "/" << info.LVDSbit;
      
      LVDSHWbitID const PMTpairBit = info.hasLVDSinfo()
        ? LVDSHWbitID{ cryostat, info.LVDSconnector, info.LVDSbit }
        : LVDSHWbitID{}
        ;
      LVDSHWbitID const adderBit = info.hasAdderInfo()
        ? LVDSHWbitID{ cryostat, info.adderConnector, info.adderBit }
        : LVDSHWbitID{}
        ;
      
      PMTchannelDB.push_back(LVDSchannel_t{
          info.channelID        // .channel
        , PMTpairBit            // .PMTpairSource
        , adderBit              // .adderSource
        });
    } // for channel in fragment
  } // for fragment ID
  
  assert(PMTchannelDB.size() == NChannels); // insider knowledge...
  
  std::sort(PMTchannelDB.begin(), PMTchannelDB.end());
  
  return PMTchannelDB;
  
} // icarus::trigger::LVDSbitMaps::buildChannelInfo()


// -----------------------------------------------------------------------------
auto icarus::trigger::LVDSbitMaps::buildPMTpairSourceMap
  (icarusDB::IICARUSChannelMapProvider const& channelMap) const
  -> std::vector<LVDSinfo_t>
{
  // note that this function must support also the case where there is no PMT
  // pair information in the channel mapping database, in which case all
  // connector information will be "invalid" and protocol demands an empty map
  
  assert(fPMTchannelInfo.size() == NChannels); // insider knowledge...
  
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
   *    but inverted because the _most_ significant bit is assigned to the
   *    _lowest_ channels
   */
  std::vector<LVDSinfo_t> sourceMap{ PMTpairBitID::NIndices() };
  std::vector<PMTpairBitID> HWtoLogicID{ LVDSHWbitID::NIndices() };
  
  unsigned int nBits = 0;
  for (LVDSchannel_t const& info: fPMTchannelInfo) { // in PMT channel order
    
    if (!info.PMTpairSource) continue; // skip invalid
    
    LVDSHWbitID::Index_t const sourceIndex = info.PMTpairSource.index();
    
    PMTpairBitID logicBit = HWtoLogicID[sourceIndex];
    
    LVDSinfo_t* mapInfo = logicBit? &(sourceMap[logicBit.index()]): nullptr;
    if (!mapInfo) { // new bit: set the source
      // which logic bit to set:
      logicBit.cryostat = info.PMTpairSource.cryostat;
      assert(logicBit.cryostat == info.channel / NChannelsPerCryostat);
      logicBit.PMTwall
        = info.channel/(NChannelsPerCryostat/PMTpairBitID::NPMTwallsPerCryostat)
        % PMTpairBitID::NCryostats;
      logicBit.statusBit = 47 - (nBits++ % 48); // 0 => 47, 1 => 46,..., 47 => 0
      if (logicBit.statusBit >= 24)
        logicBit.statusBit += 8; // 0 => 0,..., 23 => 23, 24 => 32, 25 => 33,...
      
      assert(logicBit);
      HWtoLogicID[sourceIndex] = logicBit;
      mapInfo = &(sourceMap[logicBit.index()]);
      mapInfo->ID = logicBit;
      mapInfo->source = info.PMTpairSource;
    }
    assert(mapInfo);
    assert(mapInfo->source == info.PMTpairSource);
    
    mapInfo->channels.push_back(info.channel);
    // check that they are sorted (the input list should be sorted already...)
    assert((mapInfo->channels.size() == 1)
      || (
        mapInfo->channels[mapInfo->channels.size()-2]
        < mapInfo->channels[mapInfo->channels.size()-1]
      ));
    
  } // for fPMTchannelInfo
  
  
  // --- BEGIN DEBUG -----
  mf::LogTrace debugLog{ "LVDSbitMaps" };
  debugLog << "Final PMT pair bit map:";
  for (auto const& [ bit, info ]: util::enumerate(sourceMap)) {
    debugLog << "\n[" << std::setw(3) << bit << "] logic="
      << PMTpairBitID::fromIndex(bit) << " => HW=";
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
  
  assert((nBits == 0) || (nBits == 192)); // 0 if no info in mapping database
  return (nBits == 0)? std::vector<LVDSinfo_t>{}: sourceMap;
  
} // icarus::trigger::LVDSbitMaps::buildPMTpairSourceMap()


// -----------------------------------------------------------------------------
auto icarus::trigger::LVDSbitMaps::buildAdderSourceMap
  (icarusDB::IICARUSChannelMapProvider const& channelMap) const
  -> std::vector<AdderInfo_t>
{
  // note that this function must support also the case where there is no adder
  // information in the channel mapping database, in which case all connector
  // information will be "invalid" and the protocol demands an empty map
  
  assert(fPMTchannelInfo.size() == NChannels); // insider knowledge...
  
  /*
   * The logic of the filling is a bit complicated. Class documentation explains
   * it and shows it to be:
   * 
   * east cryostat, east wall: `0x0a000000'0A000000`, channel `0` the most
   * significant bit in `a` (bit #58) and channel `45` the one in `A` (bit 26);
   * east cryostat, west wall: `0x0a000000'0A000000`, channel `90` the most
   * significant bit in `a` (bit #58) and channel `135` the one in `d` (bit 26);
   * west cryostat, east wall: `0x0a000000'0A000000`, channel `180` the most
   * significant bit in `a` (bit #58) and channel `225` the one in `d` (bit 26);
   * west cryostat, west wall: `0x0a000000'0A000000`, channel `270` the most
   * significant bit in `a` (bit #58) and channel `315` the one in `d` (bit 26).
   * 
   * How we get there:
   *  * we sort the filling with an ordered progression of channels, lowest
   *    first, and because adder bit are associated to 15 channels, we
   *    track the channels we already covered when processing the other PMT
   *    channels of the adder, and skip them when met those LVDS bit again
   *  * the destination bit is also an ordered progression within the word,
   *    but inverted because the _most_ significant bit is assigned to the
   *    _lowest_ channels
   */
  std::vector<AdderInfo_t> sourceMap{ AdderBitID::NIndices() };
  std::vector<AdderBitID> HWtoLogicID{ LVDSHWbitID::NIndices() };
  
  unsigned int nBits = 0;
  for (LVDSchannel_t const& info: fPMTchannelInfo) { // in PMT channel order
    
    if (!info.adderSource) continue; // skip invalid
    
    LVDSHWbitID::Index_t const sourceIndex = info.adderSource.index();
    
    AdderBitID logicBit = HWtoLogicID[sourceIndex];
    
    AdderInfo_t* mapInfo = logicBit? &(sourceMap[logicBit.index()]): nullptr;
    if (!mapInfo) { // new bit: set the source
      // which logic bit to set:
      logicBit.cryostat = info.adderSource.cryostat;
      assert(logicBit.cryostat == info.channel / NChannelsPerCryostat);
      logicBit.PMTwall
        = info.channel / (NChannelsPerCryostat/AdderBitID::NPMTwallsPerCryostat)
        % AdderBitID::NCryostats;
      logicBit.statusBit = 5 - (nBits++ % 6); // 0/6/... => 5, 1/7/... => 4, ...
      
      assert(logicBit);
      HWtoLogicID[sourceIndex] = logicBit;
      mapInfo = &(sourceMap[logicBit.index()]);
      mapInfo->ID = logicBit;
      mapInfo->source = info.adderSource;
    }
    assert(mapInfo);
    assert(mapInfo->source == info.adderSource);
    
    mapInfo->channels.push_back(info.channel);
    // check that they are sorted (the input list should be sorted already...)
    assert((mapInfo->channels.size() == 1)
      || (
        mapInfo->channels[mapInfo->channels.size()-2]
        < mapInfo->channels[mapInfo->channels.size()-1]
      ));
    
  } // for PMTchannelDB
  
  
  // --- BEGIN DEBUG -----
  mf::LogTrace debugLog{ "LVDSbitMaps" };
  debugLog << "Final adder bit map:";
  for (auto const& [ bit, info ]: util::enumerate(sourceMap)) {
    debugLog << "\n[" << std::setw(3) << bit << "] logic="
      << AdderBitID::fromIndex(bit) << " => HW=";
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
  
  assert((nBits == 0) || (nBits == 24)); // 0 if no info in mapping database
  return (nBits == 0)? std::vector<AdderInfo_t>{}: sourceMap;
  
} // icarus::trigger::LVDSbitMaps::buildAdderSourceMap()


// -----------------------------------------------------------------------------
std::string icarus::trigger::LVDSbitMaps::errorMessage(Errors code) {
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
    case Errors::NoAdderBitInformation:
      return "channel mapping database does not have adder bit information"s;
    case Errors::Unknown:
      return "unknown error"s;
    default:
      return "unexpected error code "s + std::to_string(static_cast<int>(code));
  } // switch
} // icarus::trigger::LVDSbitMaps::errorMessage()


// -----------------------------------------------------------------------------
