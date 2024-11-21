/**
 * @file   icaruscode/Decode/ChannelMapping/ChannelMapDumper.cxx
 * @brief  Utility dumping the content of PMT channel mapping on screen.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 * This utility can be run with any configuration file including a configuration
 * for `IICARUSChannelMap` service featuring a few hard-coded backends.
 * 
 * It may be using _art_ facilities for tool loading, but it does not run in
 * _art_ environment. So it may break without warning and without solution.
 * 
 */


// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/Legacy/ICARUSChannelMapProvider.h"
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapSQLiteProvider.h"
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapPostGresProvider.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

// LArSoft and framework libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/TestUtils/unit_test_base.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <iomanip> // std::setw()
#include <iostream>
#include <algorithm>
#include <memory>
#include <numeric> // std::iota()
#include <array>


// -----------------------------------------------------------------------------
/// Helper to split text in rows.
class Pager {
  unsigned int fLinesPerPage;
  
  unsigned int fLinesLeft = 1;
  
    public:
  Pager(unsigned int linesPerPage): fLinesPerPage{ linesPerPage } {}
  
  /// Increments the number of line and says if it starts on a new page.
  bool nextHasNewLine()
    {
      if (--fLinesLeft) return false;
      fLinesLeft = fLinesPerPage;
      return true;
    }
  bool operator() () { return nextHasNewLine(); }
  
}; // Pager


// -----------------------------------------------------------------------------
void dumpTPCmapping(icarusDB::IICARUSChannelMapProvider const& mapping) {
  
  constexpr unsigned int FragmentIDoffset = 0x1000;
  
  mf::LogVerbatim log("ChannelMapDumper");
  
  // by fragment
  unsigned int nFragmentIDs = mapping.nTPCfragmentIDs();
  log << "Reported " << nFragmentIDs << " TPC fragment IDs.";
  unsigned int nFragmentsFound = 0;
  for (
    unsigned int fragmentID = FragmentIDoffset;
    nFragmentsFound < nFragmentIDs;
    ++fragmentID
  ) {
    
    if (!mapping.hasFragmentID(fragmentID)) continue;
    ++nFragmentsFound;
    
    icarusDB::ReadoutIDVec const& boardIDs
      = mapping.getReadoutBoardVec(fragmentID);
      
    log << "\nFragment ID 0x" << std::hex << fragmentID << std::dec
      << " on crate " << mapping.getCrateName(fragmentID) << " has "
      << boardIDs.size() << " boards:";
    for (unsigned int boardID: boardIDs) log << " " << boardID;
    
  } // for fragment
  
  // by board
  unsigned int const nBoards = mapping.nTPCboardIDs();
  log << "\n\nReported " << nBoards << " TPC readout board IDs.";
  unsigned int nFoundBoards = 0;
  for (unsigned int boardID = 0; nFoundBoards < nBoards; ++boardID) {
    
    if (!mapping.hasBoardID(boardID)) continue;
    ++nFoundBoards;

    /// Returns the number of TPC board IDs known to the service.
    unsigned int const slot = mapping.getBoardSlot(boardID);
    
    icarusDB::ChannelPlanePairVec const& channelPlanes
      = mapping.getChannelPlanePair(boardID);
    log << "\nBoard " << boardID << " on slot #" << slot << " has "
      << channelPlanes.size() << " channels:";
    Pager pager{ 4 };
    for (auto [ channel, plane ]: channelPlanes) {
      if (pager.nextHasNewLine()) log << "\n     ";
      log << " CH=" << channel << " (P:" << plane << ")";
    }
    
  } // for board

  // full mapping
  icarusDB::TPCReadoutBoardToChannelMap const& boardChannelMap
    = mapping.getReadoutBoardToChannelMap(); // map<board, SlotChannelVecPair>
  log << "\n\nReported " << nBoards << " TPC boards in the mapping:";
  
  // board (unsigned int), SlotChannelVecPair = pair<slot, channels>:
  for (auto const& [ board, slotInfo ]: boardChannelMap) {
    
    auto const& [ slot, channelPlanes ] = slotInfo;
    
    log << "\nBoard " << board << " on slot #" << slot
      << " has " << channelPlanes.size() << " channels:";
    
    Pager pager{ 4 };
    for (auto [ channel, plane ]: channelPlanes) {
      if (pager.nextHasNewLine()) log << "\n     ";
      log << " CH=" << channel << " (P:" << plane << ")";
    }
    
  } // for board channel map
  
} // dumpTPCmapping()


// -----------------------------------------------------------------------------
void dumpPMTmapping(icarusDB::IICARUSChannelMapProvider const& mapping) {
  
  mf::LogVerbatim log("ChannelMapDumper");
  unsigned int nFragmentIDs = mapping.nPMTfragmentIDs();
  log << "Reported " << nFragmentIDs << " PMT fragment IDs.";
  
  // hard-coded list of fragment ID; don't like it?
  // ask for an extension of the channel mapping features.
  std::vector<unsigned int> FragmentIDs(nFragmentIDs);
  std::iota(FragmentIDs.begin(), FragmentIDs.end(), 0x2000);
  
  log << "\nPMT fragment IDs:";
  for (auto const [ iFragment, fragmentID]: util::enumerate(FragmentIDs)) {
    unsigned int const effFragmentID = fragmentID & 0xFF;
    
    if (!mapping.hasPMTDigitizerID(effFragmentID)) {
      log << "\n[" << iFragment << "] " << std::hex << fragmentID << std::dec
        << " not found in the database (as "
        << std::hex << effFragmentID << std::dec << ")";
      continue;
    }
    
    icarusDB::PMTdigitizerInfoVec digitizerChannels
      = mapping.getPMTchannelInfo(effFragmentID);
    
    std::sort(digitizerChannels.begin(), digitizerChannels.end());
    
    
    log
      << "\n[" << iFragment << "] 0x" << std::hex << fragmentID << std::dec
      << " includes " << digitizerChannels.size()
      << " LArSoft channels between " << digitizerChannels.front().channelID
      << " and " << digitizerChannels.back().channelID
      << " [board channel info in brackets]:";
    Pager pager{ 3 };
    for(auto const & chInfo: digitizerChannels) {
      if (pager.nextHasNewLine()) log << "\n     ";
      log << " " << std::setw(3) << chInfo.channelID
        << " [" << chInfo.digitizerLabel << "@"
        << std::setfill('0') << std::setw(2) << chInfo.digitizerChannelNo;
      if (chInfo.hasLVDSinfo()) {
        log << ", LVDS:" << chInfo.LVDSconnector << "-"
          << std::setw(2) << chInfo.LVDSbit;
      }
      if (chInfo.hasAdderInfo()) {
        log << ", adder:" << chInfo.adderConnector << "-"
          << std::setw(2) << chInfo.adderBit;
      }
      log << "]";
    } // for channel
    
  } // for fragment
  
} // dumpPMTmapping()


// -----------------------------------------------------------------------------
void dumpCRTmapping(icarusDB::IICARUSChannelMapProvider const& mapping) {
  
  mf::LogVerbatim log("ChannelMapDumper");
  
  constexpr unsigned int MaxHWmac = 1000;
  
  log << "CRT address mapping:";
  for (unsigned int const HWaddress: util::counter(MaxHWmac)) {
    unsigned int const simMac = mapping.getSimMacAddress(HWaddress);
    if (simMac == 0) continue; // not present
    log << "\nSide CRT HW=" << HWaddress << " -> sim=" << simMac;
  }
  for (unsigned int const HWaddress: util::counter(MaxHWmac)) {
    unsigned int const simMac = mapping.gettopSimMacAddress(HWaddress);
    if (simMac == 0) continue; // not present
    log << "\nTop CRT  HW=" << HWaddress << " -> sim=" << simMac;
  }
  
  // the interface does not provide any clue on which mac5 and chan values are
  // available, so we have to skip this
  // std::pair<double, double> const = mapping.getSideCRTCalibrationMap(mac5, chan);
   
} // dumpCRTmapping()


// -----------------------------------------------------------------------------
void dumpMapping(icarusDB::IICARUSChannelMapProvider const& channelMapping) {
  
  mf::LogVerbatim("ChannelMapDumper") << std::string(80, '-');
  dumpTPCmapping(channelMapping);
  
  mf::LogVerbatim("ChannelMapDumper") << std::string(80, '-');
  dumpPMTmapping(channelMapping);
  
  mf::LogVerbatim("ChannelMapDumper") << std::string(80, '-');
  dumpCRTmapping(channelMapping);
  
  mf::LogVerbatim("ChannelMapDumper") << std::string(80, '-');
}


// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
  
  using Environment
    = testing::TesterEnvironment<testing::BasicEnvironmentConfiguration>;
  
  testing::BasicEnvironmentConfiguration config("ChannelMapDumper");

  //
  // parameter parsing
  //
  int iParam = 0;

  // first argument: configuration file (mandatory)
  if (++iParam < argc)
    config.SetConfigurationPath(argv[iParam]);
  else {
    std::cerr << "FHiCL configuration file path required as first argument!"
      << std::endl;
    return 1;
  }
  
  Environment const Env { config };
  
  //
  // create the service provider object
  //
  // limited selection of available service providers, hard-coded
  // 
  std::unique_ptr<icarusDB::IICARUSChannelMapProvider> channelMapping;
  fhicl::ParameterSet channelMapConfig
    = Env.ServiceParameters("IICARUSChannelMap");
  std::string const serviceType = channelMapConfig.get<std::string>
    ("service_provider", "ICARUSChannelMapSQLite");
  channelMapConfig.erase("service_provider"); // if not present, still friends
  mf::LogVerbatim("ChannelMapDumper")
    << "Initializing a provider of type '" << serviceType << "'";
  mf::LogTrace("ChannelMapDumper")
    << "with configuration:\n{" << channelMapConfig.to_indented_string(1)
    << "\n}";
  if (serviceType == "ICARUSChannelMapSQLite") {
    channelMapping = std::make_unique<icarusDB::ICARUSChannelMapSQLiteProvider>
      (channelMapConfig);
  }
  else if (serviceType == "ICARUSChannelMapPostGres") {
    channelMapping
      = std::make_unique<icarusDB::ICARUSChannelMapPostGresProvider>
      (channelMapConfig);
  }
  else if (serviceType == "ICARUSChannelMap") { // legacy
    channelMapping = std::make_unique<icarusDB::ICARUSChannelMapProvider>
      (channelMapConfig);
  }
  else {
    mf::LogError("ChannelMapDumper")
      << "Fatal: don't know what to do with IICARUSChannelMap.service_type: '"
      << serviceType << "'.";
    return 1;
  }
  
  //
  // dump
  //
  for (icarusDB::RunPeriod const period: icarusDB::RunPeriods::All) {
    mf::LogVerbatim("ChannelMapDumper")
      << std::string(80, '=') << "\nRun period #"
      << static_cast<unsigned int>(period) << ":\n";
    channelMapping->forPeriod(period);
    dumpMapping(*channelMapping);
  }
  
  mf::LogVerbatim("ChannelMapDumper")
    << std::string(80, '=')
    << "\nDumped " << icarusDB::RunPeriods::All.size() << " run periods."
    ;
  
  return 0;
} // main()


