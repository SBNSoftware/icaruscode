/**
 * @file   icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.cxx
 * @brief  Interface class for hardware/software channel mapping for ICARUS.
 * @see    icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h
 */

// library header
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"

// C/C++ standard libraries
#include <map>
#include <regex>


// -----------------------------------------------------------------------------
auto icarusDB::IICARUSChannelMapProvider::unpackDigitizerInfo
  (std::string const& label, int format /* = 0 */)
  -> icarusDB::PMTChannelInfo_t::DigitizerInfo_t
{
  static std::regex const LabelPattern{ R"(([EW])([EW])-(TOP|BOT)-([A-C]))" };
  
  static std::map<std::string, short int> const EastWestMap
    = { { "E", 0 }, { "W", 1 } };
  
  static std::map<std::string, short int> const EndSideMap
    = { { "BOT", 0 }, { "TOP", 1 } };
  
  static std::map<std::string, short int> const SlotMap
    = { { "A", 0 }, { "B", 1 }, { "C", 2 } };
  
  
  std::smatch match;
  if (!std::regex_match(label, match, LabelPattern)) return {};
  
  return {
    /* .cryostat = */ EastWestMap.at(match.str(1)),
    /* .PMTwall  = */ EastWestMap.at(match.str(2)),
    /* .endSide  = */ EndSideMap.at(match.str(3)),
    /* .slot     = */ SlotMap.at(match.str(4))
  };
  
} // icarusDB::IICARUSChannelMapProvider::unpackDigitizerInfo()


// -----------------------------------------------------------------------------
