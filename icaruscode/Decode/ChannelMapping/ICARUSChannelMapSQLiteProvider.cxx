/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapSQLiteProvider.cxx
 * @author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapSQLiteProvider.h
 */

// library header
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapSQLiteProvider.h"

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"

// C++ standard libraries
#include <string>
#include <iostream>
#include <cassert>



// -----------------------------------------------------------------------------
icarusDB::ICARUSChannelMapSQLiteProvider::ICARUSChannelMapSQLiteProvider
  (Config const& config)
  : fDiagnosticOutput{ config.DiagnosticOutput() }
  , fLogCategory{ config.LogCategory() }
  , fChannelMappingAlg
    { std::make_unique<ChannelMappingAlg_t>(config.ChannelMappingTool()) }
{
}


// -----------------------------------------------------------------------------
bool icarusDB::ICARUSChannelMapSQLiteProvider::forRun(int run) {
  
  RunPeriod period = RunPeriod::NPeriods;
  try {
    period = RunPeriods::withRun(run);
  } catch (...) {
    mf::LogError{ fLogCategory }
      << "Failed to associate a run period number to run " << run << ".";
    throw;
  }
  
  if (period == RunPeriod::NPeriods) {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "forRun(): can't determine the period of run " << run << "!\n";
  }
  
  return forPeriod(period);
  
} // icarusDB::ICARUSChannelMapSQLiteProvider::forRun()


// -----------------------------------------------------------------------------
bool icarusDB::ICARUSChannelMapSQLiteProvider::forPeriod
  (icarusDB::RunPeriod period)
{
  
  // if cache is not invalidated, we don't refresh it
  if (!fChannelMappingAlg->SelectPeriod(period)) return false;
  
  readFromDatabase();
  return true;
}


// -----------------------------------------------------------------------------
bool icarusDB::ICARUSChannelMapSQLiteProvider::hasFragmentID
  (const unsigned int fragmentID) const 
{
  return fFragmentToReadoutMap.find(fragmentID) != fFragmentToReadoutMap.end();
}


// -----------------------------------------------------------------------------
unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::nTPCfragmentIDs() const {
  return fFragmentToReadoutMap.size();
}


// -----------------------------------------------------------------------------
std::string const& icarusDB::ICARUSChannelMapSQLiteProvider::getCrateName
  (const unsigned int fragmentID) const
{
  auto const fragToReadoutItr = fFragmentToReadoutMap.find(fragmentID);

  if (fragToReadoutItr == fFragmentToReadoutMap.end()) {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Fragment ID " << fragmentID
      << " not found in lookup map when looking up crate name \n";
  }

  return fragToReadoutItr->second.first;
}


// -----------------------------------------------------------------------------
auto icarusDB::ICARUSChannelMapSQLiteProvider::getReadoutBoardVec
  (const unsigned int fragmentID) const -> icarusDB::ReadoutIDVec const&
{
  auto const fragToReadoutItr = fFragmentToReadoutMap.find(fragmentID);

  if (fragToReadoutItr == fFragmentToReadoutMap.end()) {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Fragment ID " << fragmentID
      << " not found in lookup map when looking up board vector.\n";
  }

  return fragToReadoutItr->second.second;
}


// -----------------------------------------------------------------------------
auto icarusDB::ICARUSChannelMapSQLiteProvider::getReadoutBoardToChannelMap()
  const -> const TPCReadoutBoardToChannelMap&
{
  return fReadoutBoardToChannelMap;
}


// -----------------------------------------------------------------------------
bool icarusDB::ICARUSChannelMapSQLiteProvider::hasBoardID
  (const unsigned int boardID)  const
{
  return
    fReadoutBoardToChannelMap.find(boardID) != fReadoutBoardToChannelMap.end();
}


// -----------------------------------------------------------------------------
unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::nTPCboardIDs() const {
  return fReadoutBoardToChannelMap.size();
}


// -----------------------------------------------------------------------------
unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::getBoardSlot
  (const unsigned int boardID)  const
{
  auto const readoutBoardItr = fReadoutBoardToChannelMap.find(boardID);

  if (readoutBoardItr == fReadoutBoardToChannelMap.end()) {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Board ID " << boardID
      << " not found in lookup map when looking up board slot.\n";
  }

  return readoutBoardItr->second.first;
}


// -----------------------------------------------------------------------------
auto icarusDB::ICARUSChannelMapSQLiteProvider::getChannelPlanePair
  (const unsigned int boardID) const -> ChannelPlanePairVec const&
{
  auto const readoutBoardItr = fReadoutBoardToChannelMap.find(boardID);

  if (readoutBoardItr == fReadoutBoardToChannelMap.end()) {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Board ID " << boardID
      << " not found in lookup map when looking up channel/plane pair.\n";
  }

  return readoutBoardItr->second.second;
}


// -----------------------------------------------------------------------------
bool icarusDB::ICARUSChannelMapSQLiteProvider::hasPMTDigitizerID
  (const unsigned int fragmentID) const
{
  return findPMTfragmentEntry(fragmentID) != nullptr;
}


// -----------------------------------------------------------------------------
unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::nPMTfragmentIDs() const {
  return fFragmentToDigitizerMap.size();
}


// -----------------------------------------------------------------------------
auto icarusDB::ICARUSChannelMapSQLiteProvider::getChannelIDPairVec
  (const unsigned int fragmentID) const
  -> DigitizerChannelChannelIDPairVec const&
{
  mf::LogTrace{ "ICARUSChannelMapSQLiteProvider" }
    << "Call to: ICARUSChannelMapSQLiteProvider::getChannelIDPairVec("
    << fragmentID << ")";
  
  DigitizerChannelChannelIDPairVec const* digitizerPair
    = findPMTfragmentEntry(fragmentID);
  
  if (digitizerPair) return *digitizerPair;
  throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
    << "Fragment ID " << fragmentID
    << " not found in lookup map when looking for PMT channel info.\n";
  
}


// -----------------------------------------------------------------------------
unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::getSimMacAddress
  (const unsigned int hwmacaddress)  const
{
  for(auto const& [ hw, sim ]
    : util::const_values(fCRTChannelIDToHWtoSimMacAddressPairMap)
  ) {
    if (hw == hwmacaddress) return sim;
  }

  return 0;
}


// -----------------------------------------------------------------------------
unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::gettopSimMacAddress
  (const unsigned int hwmacaddress) const
{
#if 0 // FIXME
  unsigned int   simmacaddress = 0;

  for(const auto& pair : fTopCRTHWtoSimMacAddressPairMap){
    if (pair.first == hwmacaddress)
simmacaddress = pair.second;
  }

  return simmacaddress;
#else
  // untested:
  auto const it = fTopCRTHWtoSimMacAddressPairMap.find(hwmacaddress);
  return (it == fTopCRTHWtoSimMacAddressPairMap.end())? 0: it->second;
  
#endif
}


// -----------------------------------------------------------------------------
auto icarusDB::ICARUSChannelMapSQLiteProvider::getSideCRTCalibrationMap
  (int mac5, int chan) const -> std::pair<double, double>
{
  auto const itGainAndPedestal
    = fSideCRTChannelToCalibrationMap.find({ mac5, chan });
  return (itGainAndPedestal == fSideCRTChannelToCalibrationMap.cend())
    ? std::pair{ -99., -99. }: itGainAndPedestal->second;
}


// -----------------------------------------------------------------------------
auto icarusDB::ICARUSChannelMapSQLiteProvider::findPMTfragmentEntry
  (unsigned int fragmentID) const -> DigitizerChannelChannelIDPairVec const*
{
  auto it = fFragmentToDigitizerMap.find(PMTfragmentIDtoDBkey(fragmentID));
  return (it == fFragmentToDigitizerMap.end())? nullptr: &(it->second);
}


// -----------------------------------------------------------------------------
void icarusDB::ICARUSChannelMapSQLiteProvider::readFromDatabase() {

  mf::LogInfo{ fLogCategory } << "Building the channel mapping";

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // TPC fragment-based mapping
  cet::cpu_timer theClockFragmentIDs;
  theClockFragmentIDs.start();
  fFragmentToReadoutMap.clear();
  if (fChannelMappingAlg->BuildTPCFragmentIDToReadoutIDMap(fFragmentToReadoutMap))
  {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Cannot recover the TPC fragment ID channel map from the database.\n";
  }
  else if (fDiagnosticOutput) {
    
    mf::LogVerbatim log{ "ICARUSChannelMapSQLiteProvider" };
    log << "FragmentID to Readout ID map has " << fFragmentToReadoutMap.size()
      << " elements";
    for(auto const& [ fragmentID, crateAndBoards ]: fFragmentToReadoutMap) {
      log << "\n   Frag: " << std::hex << fragmentID << std::dec << ", Crate: "
        << crateAndBoards.first << ", # boards: "
        << crateAndBoards.second.size();
    }
  }
  theClockFragmentIDs.stop();
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // TPC readout-board-based mapping
  double fragmentIDsTime = theClockFragmentIDs.accumulated_real_time();

  cet::cpu_timer theClockReadoutIDs;
  theClockReadoutIDs.start();

  fReadoutBoardToChannelMap.clear();
  if (fChannelMappingAlg->BuildTPCReadoutBoardToChannelMap(fReadoutBoardToChannelMap))
  {
    mf::LogError{ "ICARUSChannelMapSQLiteProvider" }
      << "******* FAILED TO CONFIGURE CHANNEL MAP ********";
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Failed to read the database.\n";
  }

  theClockReadoutIDs.stop();
  double readoutIDsTime = theClockReadoutIDs.accumulated_real_time();

  mf::LogInfo{ "ICARUSChannelMapSQLiteProvider" }
    << "==> FragmentID map time: " << fragmentIDsTime << ", Readout IDs time: "
    << readoutIDsTime;
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // PMT channel mapping
  fFragmentToDigitizerMap.clear();
  if (fChannelMappingAlg->BuildFragmentToDigitizerChannelMap(fFragmentToDigitizerMap))
  {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Cannot recover the PMT fragment ID channel map from the database.\n";
  }
  else if (fDiagnosticOutput) {
    mf::LogVerbatim log{ "ICARUSChannelMapSQLiteProvider" };
    log << "FragmentID to Readout ID map has " << fFragmentToDigitizerMap.size()
      << " Fragment IDs";
    
    for(auto const& [ fragmentID, digitizers ]: fFragmentToDigitizerMap) {
      log << "\n   Frag: " << std::hex << fragmentID << std::dec
        << ", # pairs: " << digitizers.size();
    }
  }
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Side CRT channel mapping
  fCRTChannelIDToHWtoSimMacAddressPairMap.clear();
  if (fChannelMappingAlg->BuildCRTChannelIDToHWtoSimMacAddressPairMap(fCRTChannelIDToHWtoSimMacAddressPairMap))
  {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Cannot recover the HW MAC Address  from the database.\n";
  }
  else if (fDiagnosticOutput) {
    mf::LogVerbatim log{ "ICARUSChannelMapSQLiteProvider" };
    log << "ChannelID to MacAddress map has "
      << fCRTChannelIDToHWtoSimMacAddressPairMap.size() << " Channel IDs";
    for(auto const& [ channel, addresses ]: fCRTChannelIDToHWtoSimMacAddressPairMap) {
      log <<"\n ChannelID: "<< channel
        << ", hw mac address: " << addresses.first
        << ", sim mac address: " << addresses.second;
    }
  }
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Top CRT channel mapping
  fTopCRTHWtoSimMacAddressPairMap.clear();
  if (fChannelMappingAlg->BuildTopCRTHWtoSimMacAddressPairMap(fTopCRTHWtoSimMacAddressPairMap))
  {
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Cannot recover the Top CRT HW MAC Address  from the database.\n";
  }
  else if (fDiagnosticOutput) {
    mf::LogVerbatim log{ "ICARUSChannelMapSQLiteProvider" };
    log << "Top CRT MacAddress map has " << fTopCRTHWtoSimMacAddressPairMap.size() << " rows";
    for(auto const [ hwaddress, simaddress ]: fTopCRTHWtoSimMacAddressPairMap) {
      log << "\n hw mac address: " << hwaddress
        << ", sim mac address: " << simaddress;
    }
  }
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // CRT Charge Calibration initialization
  fSideCRTChannelToCalibrationMap.clear();
  if (fChannelMappingAlg->BuildSideCRTCalibrationMap(fSideCRTChannelToCalibrationMap))
  {
    mf::LogError{ "ICARUSChannelMapSQLiteProvider" }
      << "******* FAILED TO CONFIGURE CRT Calibration  ********";
    throw cet::exception{ "ICARUSChannelMapSQLiteProvider" }
      << "Cannot recover the CRT charge calibration information from the database.\n";
  }
  else if (fDiagnosticOutput) {
    mf::LogVerbatim log{ "ICARUSChannelMapSQLiteProvider" };
    log << "side crt calibration map has "
      << fSideCRTChannelToCalibrationMap.size() << " list of rows";
    
    for(auto const& [ key, calib ]: fSideCRTChannelToCalibrationMap) {
      log << "\n mac5: "<< key.first << ", chan: " << key.second
        << ", Gain: " << calib.first << ", Pedestal: " << calib.second;
    }
  }
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
} // icarusDB::ICARUSChannelMapSQLiteProvider::readFromDatabase()


// -----------------------------------------------------------------------------
constexpr unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::PMTfragmentIDtoDBkey
  (unsigned int fragmentID)
{
  /*
   * PMT channel mapping database stores the board number (0-23) as key.
   * Fragment ID are currently in the pattern 0x20xx, with xx the board number.
   */
  
  // protest if this is a fragment not from the PMT;
  // but make an exception for old PMT fragment IDs (legacy)
  assert(((fragmentID & ~0xFF) == 0x00) || ((fragmentID & ~0xFF) == 0x20));
  
  return fragmentID & 0xFF;
  
} // ICARUSChannelMapSQLiteProvider::PMTfragmentIDtoDBkey()


// -----------------------------------------------------------------------------
constexpr unsigned int icarusDB::ICARUSChannelMapSQLiteProvider::DBkeyToPMTfragmentID
  (unsigned int DBkey)
{
  /*
   * PMT channel mapping database stores the board number (0-23) as key.
   * Fragment ID are currently in the pattern 0x20xx, with xx the board number.
   */
  
  // protest if this is a fragment not from the PMT;
  // but make an exception for old PMT fragment IDs (legacy)
  assert((DBkey & 0xFF) < 24);
  
  return (DBkey & 0xFF) | 0x2000;
  
} // icarusDB::ICARUSChannelMapSQLiteProvider::PMTfragmentIDtoDBkey()


// -----------------------------------------------------------------------------
