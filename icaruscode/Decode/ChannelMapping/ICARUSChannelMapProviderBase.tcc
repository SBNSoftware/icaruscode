/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.txx
 * @author T. Usher (factorised by G. Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.h
 */

// library header
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.h"

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()

// framework libraries
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"

// C++ standard libraries
#include <string>
#include <cassert>


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::ICARUSChannelMapProviderBase
  (Config const& config)
  : icarus::ns::util::mfLoggingClass
    { config.LogCategory().value_or(config.ChannelMappingTool().LogCategory()) }
  , fDiagnosticOutput  { config.DiagnosticOutput() }
  , fChannelMappingAlg { config.ChannelMappingTool() }
{
  addCacheTags({ "TPC", "PMT", "CRT" }); // all caches updated at the same time
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
bool icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::forRun(int run) {
  
  RunPeriod period = RunPeriod::NPeriods;
  try {
    period = RunPeriods::withRun(run);
  } catch (...) {
    mfLogError()
      << "Failed to associate a run period number to run " << run << ".";
    throw;
  }
  
  if (period == RunPeriod::NPeriods) {
    throw myException()
      << "forRun(): can't determine the period of run " << run << "!\n";
  }
  
  return forPeriod(period);
  
} // icarusDB::ICARUSChannelMapProviderBase<>::forRun()


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
bool icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::forPeriod
  (icarusDB::RunPeriod period)
{
  
  // if cache is not invalidated, we don't refresh it
  if (!fChannelMappingAlg.SelectPeriod(period)) return false;
  
  readFromDatabase();
  return true;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
bool icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::hasFragmentID
  (const unsigned int fragmentID) const
{
  return
    fTPCFragmentToReadoutMap.find(fragmentID) != fTPCFragmentToReadoutMap.end();
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
unsigned int icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::nTPCfragmentIDs() const {
  return fTPCFragmentToReadoutMap.size();
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
std::string const& icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getCrateName
  (const unsigned int fragmentID) const
{
  auto const fragToReadoutItr = fTPCFragmentToReadoutMap.find(fragmentID);

  if (fragToReadoutItr == fTPCFragmentToReadoutMap.end()) {
    throw myException() << "Fragment ID " << fragmentID
      << " not found in lookup map when looking up crate name \n";
  }

  return fragToReadoutItr->second.first;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
auto icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getReadoutBoardVec
  (const unsigned int fragmentID) const -> icarusDB::ReadoutIDVec const&
{
  auto const fragToReadoutItr = fTPCFragmentToReadoutMap.find(fragmentID);

  if (fragToReadoutItr == fTPCFragmentToReadoutMap.end()) {
    throw myException() << "Fragment ID " << fragmentID
      << " not found in lookup map when looking up board vector.\n";
  }

  return fragToReadoutItr->second.second;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
auto icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getReadoutBoardToChannelMap()
  const -> const TPCReadoutBoardToChannelMap&
{
  return fTPCReadoutBoardToChannelMap;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
bool icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::hasBoardID
  (const unsigned int boardID)  const
{
  return
    fTPCReadoutBoardToChannelMap.find(boardID) != fTPCReadoutBoardToChannelMap.end();
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
unsigned int icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::nTPCboardIDs() const {
  return fTPCReadoutBoardToChannelMap.size();
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
unsigned int icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getBoardSlot
  (const unsigned int boardID)  const
{
  auto const readoutBoardItr = fTPCReadoutBoardToChannelMap.find(boardID);

  if (readoutBoardItr == fTPCReadoutBoardToChannelMap.end()) {
    throw myException() << "Board ID " << boardID
      << " not found in lookup map when looking up board slot.\n";
  }

  return readoutBoardItr->second.first;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
auto icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getChannelPlanePair
  (const unsigned int boardID) const -> ChannelPlanePairVec const&
{
  auto const readoutBoardItr = fTPCReadoutBoardToChannelMap.find(boardID);

  if (readoutBoardItr == fTPCReadoutBoardToChannelMap.end()) {
    throw myException() << "Board ID " << boardID
      << " not found in lookup map when looking up channel/plane pair.\n";
  }

  return readoutBoardItr->second.second;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
bool icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::hasPMTDigitizerID
  (const unsigned int fragmentID) const
{
  return findPMTfragmentEntry(fragmentID) != nullptr;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
unsigned int icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::nPMTfragmentIDs()
  const
{
  return fPMTFragmentToDigitizerMap.size();
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
auto icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getPMTchannelInfo
  (unsigned int fragmentID) const -> PMTdigitizerInfoVec const&
{
  mfLogTrace() << "Call to: ICARUSChannelMapProviderBase::getPMTchannelInfo("
    << fragmentID << ")";
  
  PMTdigitizerInfoVec const* digitizerInfo = findPMTfragmentEntry(fragmentID);
  
  if (digitizerInfo) return *digitizerInfo;
  throw myException() << "Fragment ID " << fragmentID
    << " not found in lookup map when looking for PMT channel info.\n";
  
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
unsigned int icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getSimMacAddress
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
template <typename ChMapAlg>
unsigned int
icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::gettopSimMacAddress
  (const unsigned int hwmacaddress) const
{
  auto const it = fTopCRTHWtoSimMacAddressPairMap.find(hwmacaddress);
  return (it == fTopCRTHWtoSimMacAddressPairMap.end())? 0: it->second;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
auto icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::getSideCRTCalibrationMap
  (int mac5, int chan) const -> std::pair<double, double>
{
  auto const itGainAndPedestal
    = fSideCRTChannelToCalibrationMap.find({ mac5, chan });
  return (itGainAndPedestal == fSideCRTChannelToCalibrationMap.cend())
    ? std::pair{ -99., -99. }: itGainAndPedestal->second;
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
auto icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::findPMTfragmentEntry
  (unsigned int fragmentID) const -> PMTdigitizerInfoVec const*
{
  auto it = fPMTFragmentToDigitizerMap.find(PMTfragmentIDtoDBkey(fragmentID));
  return (it == fPMTFragmentToDigitizerMap.end())? nullptr: &(it->second);
}


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
void icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::readFromDatabase() {

  mfLogInfo() << "Building the channel mapping";

  updateCacheID("TPC");
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // TPC fragment-based mapping
  cet::cpu_timer theClockFragmentIDs;
  theClockFragmentIDs.start();
  fTPCFragmentToReadoutMap.clear();
  if (
   fChannelMappingAlg.BuildTPCFragmentIDToReadoutIDMap(fTPCFragmentToReadoutMap)
  ) {
    throw myException()
      << "Cannot recover the TPC fragment ID channel map from the database.\n";
  }
  else if (fDiagnosticOutput) {
    
    auto log = mfLogVerbatim();
    log << "FragmentID to Readout ID map has " << fTPCFragmentToReadoutMap.size()
      << " elements";
    for(auto const& [ fragmentID, crateAndBoards ]: fTPCFragmentToReadoutMap) {
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

  fTPCReadoutBoardToChannelMap.clear();
  if (fChannelMappingAlg.BuildTPCReadoutBoardToChannelMap
    (fTPCReadoutBoardToChannelMap)
  ) {
    mfLogError() << "******* FAILED TO CONFIGURE CHANNEL MAP ********";
    throw myException() << "Failed to read the database.\n";
  }

  theClockReadoutIDs.stop();
  double readoutIDsTime = theClockReadoutIDs.accumulated_real_time();

  mfLogInfo() << "==> FragmentID map time: " << fragmentIDsTime
    << ", Readout IDs time: " << readoutIDsTime;
  
  
  updateCacheID("PMT");
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // PMT channel mapping
  fPMTFragmentToDigitizerMap.clear();
  if (fChannelMappingAlg.BuildPMTFragmentToDigitizerChannelMap(fPMTFragmentToDigitizerMap))
  {
    throw myException() 
      << "Cannot recover the PMT fragment ID channel map from the database.\n";
  }
  else if (fDiagnosticOutput) {
    auto log = mfLogVerbatim();
    log << "FragmentID to Readout ID map has " << fPMTFragmentToDigitizerMap.size()
      << " Fragment IDs";
    
    for(auto const& [ fragmentID, digitizers ]: fPMTFragmentToDigitizerMap) {
      log << "\n   Frag: " << std::hex << fragmentID << std::dec
        << ", # pairs: " << digitizers.size();
    }
  }
  
  
  updateCacheID("CRT");
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Side CRT channel mapping
  fCRTChannelIDToHWtoSimMacAddressPairMap.clear();
  if (fChannelMappingAlg.BuildCRTChannelIDToHWtoSimMacAddressPairMap(fCRTChannelIDToHWtoSimMacAddressPairMap))
  {
    throw myException()
      << "Cannot recover the HW MAC Address  from the database.\n";
  }
  else if (fDiagnosticOutput) {
    auto log = mfLogVerbatim();
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
  if (fChannelMappingAlg.BuildTopCRTHWtoSimMacAddressPairMap(fTopCRTHWtoSimMacAddressPairMap))
  {
    throw myException()
      << "Cannot recover the Top CRT HW MAC Address  from the database.\n";
  }
  else if (fDiagnosticOutput) {
    auto log = mfLogVerbatim();
    log << "Top CRT MacAddress map has " << fTopCRTHWtoSimMacAddressPairMap.size() << " rows";
    for(auto const [ hwaddress, simaddress ]: fTopCRTHWtoSimMacAddressPairMap) {
      log << "\n hw mac address: " << hwaddress
        << ", sim mac address: " << simaddress;
    }
  }
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // CRT Charge Calibration initialization
  fSideCRTChannelToCalibrationMap.clear();
  if (fChannelMappingAlg.BuildSideCRTCalibrationMap
    (fSideCRTChannelToCalibrationMap)
  ) {
    mfLogError() << "******* FAILED TO CONFIGURE CRT Calibration  ********";
    throw myException()
      << "Cannot recover the CRT charge calibration information from the database.\n";
  }
  else if (fDiagnosticOutput) {
    auto log = mfLogVerbatim();
    log << "side crt calibration map has "
      << fSideCRTChannelToCalibrationMap.size() << " list of rows";
    
    for(auto const& [ key, calib ]: fSideCRTChannelToCalibrationMap) {
      log << "\n mac5: "<< key.first << ", chan: " << key.second
        << ", Gain: " << calib.first << ", Pedestal: " << calib.second;
    }
  }
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
} // icarusDB::ICARUSChannelMapProviderBase<>::readFromDatabase()


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
constexpr unsigned int icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::PMTfragmentIDtoDBkey
  (unsigned int fragmentID)
{
  /*
   * PMT channel mapping database stores the board number (0-23) as key.
   * Fragment ID are currently in the pattern 0x20xx, with xx the board number.
   */
  
  // protest if this is a fragment not from the PMT;
  // but make an exception for old PMT fragment IDs (legacy)
  assert(((fragmentID & ~0xFF) == 0x0000) || ((fragmentID & ~0xFF) == 0x2000));
  
  return fragmentID & 0xFF;
  
} // ICARUSChannelMapProviderBase<>::PMTfragmentIDtoDBkey()


// -----------------------------------------------------------------------------
template <typename ChMapAlg>
constexpr unsigned int
icarusDB::ICARUSChannelMapProviderBase<ChMapAlg>::DBkeyToPMTfragmentID
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
  
} // icarusDB::ICARUSChannelMapProviderBase<>::PMTfragmentIDtoDBkey()


// -----------------------------------------------------------------------------
