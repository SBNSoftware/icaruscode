/**
 * @file   icaruscode/PMT/Status/PMTChannelStatusProvider.cxx
 * @brief  Implementation of PMT channel status provider.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @see    icaruscode/PMT/Status/PMTChannelStatusProvider.h
 */

#include "icaruscode/PMT/Status/PMTChannelStatusProvider.h"

// LArSoft libraries
#include "larevt/CalibrationDBI/IOVData/TimeStampDecoder.h"

// Framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <iomanip> // std::setw()


// -----------------------------------------------------------------------------
icarusDB::PMTChannelStatusProvider::PMTChannelStatusProvider
  (const fhicl::ParameterSet& pset)
  : fVerbose    ( pset.get<bool>       ("Verbose",     false                    ) )
  , fLogCategory( pset.get<std::string>("LogCategory", "PMTChannelStatusProvider") )
  , fStatusTag  ( pset.get<std::string>("StatusTag",   ""                       ) )
  , fDB         ( pset.get<std::string>("DBname", "pmt_voltage_data"),
                  "", "", fStatusTag, true, false )
{
  int const defaultStatusInt = pset.get<int>("DefaultStatus", static_cast<int>(kON));
  fDefault.status  = static_cast<PMTChannelStatusValue>(defaultStatusInt);
  fDefault.voltage = pset.get<double>("DefaultVoltage", 1500.0);

  mf::LogInfo(fLogCategory)
    << "PMTChannelStatusProvider connected to '"
    << pset.get<std::string>("DBname", "pmt_voltage_data")
    << "' DB, tag '" << fStatusTag << "'"
    << ", default status: " << defaultStatusInt
    << ", default voltage: " << fDefault.voltage << " V";
}


// -----------------------------------------------------------------------------
void icarusDB::PMTChannelStatusProvider::readStatusFromDB(unsigned int run)
{
  mf::LogInfo(fLogCategory)
    << "Reading PMT channel statuses from database for run " << run;

  bool const ret = fDB.UpdateData( RunToDatabaseTimestamp(run) );
  mf::LogTrace(fLogCategory)
    << "Status" << (ret ? "" : " not") << " updated for run " << run
    << "\nFetched IoV [ " << fDB.CachedStart().DBStamp()
    << " ; " << fDB.CachedEnd().DBStamp()
    << " ] to cover t=" << RunToDatabaseTimestamp(run)
    << " [=" << lariov::TimeStampDecoder::DecodeTimeStamp(RunToDatabaseTimestamp(run)).DBStamp() << "]";

  std::vector<unsigned int> channelList;
  if (int const res = fDB.GetChannelList(channelList); res != 0) {
    throw cet::exception(fLogCategory)
      << "GetChannelList() returned " << res << " on run " << run << " query\n";
  }

  if (channelList.empty()) {
    throw cet::exception(fLogCategory)
      << "Got an empty channel list for run " << run << "\n";
  }

  fDatabaseStatus.clear();

  mf::LogDebug log(fLogCategory);
  log << "Loading status for " << channelList.size() << " channels in run " << run;

  for (auto const channel : channelList) {

    long statusInt = 0;
    int error = fDB.GetNamedChannelData(channel, "status", statusInt);
    if (error) {
      throw cet::exception(fLogCategory)
        << "Error (code " << error << ") reading 'status' for channel "
        << channel << "\n";
    }

    double voltage = 0.0;
    error = fDB.GetNamedChannelData(channel, "voltage", voltage);
    if (error) {
      throw cet::exception(fLogCategory)
        << "Error (code " << error << ") reading 'voltage' for channel "
        << channel << "\n";
    }

    fDatabaseStatus[channel].status  = static_cast<PMTChannelStatusValue>(statusInt);
    fDatabaseStatus[channel].voltage = voltage;

    if (fVerbose)
      log << "\n  channel " << std::setw(3) << channel
          << "  status " << statusInt
          << "  voltage " << voltage << " V";

  } // for channel

} // readStatusFromDB()


// -----------------------------------------------------------------------------
std::uint64_t
icarusDB::PMTChannelStatusProvider::RunToDatabaseTimestamp(unsigned int run) const
{
  // Same convention as other ICARUS DB services:
  // run XXXXX -> (XXXXX + 1'000'000'000) * 1'000'000'000  [nanoseconds]
  // e.g. run XXXXX -> 10000XXXXX * 1e9
  std::uint64_t timestamp = std::uint64_t(run) + 1'000'000'000ULL;
  timestamp *= 1'000'000'000ULL;

  if (fVerbose)
    mf::LogInfo(fLogCategory)
      << "Run " << run << " status from DB timestamp " << timestamp;

  return timestamp;
}


// -----------------------------------------------------------------------------
icarusDB::PMTChannelStatusProvider::ChannelSet_t
icarusDB::PMTChannelStatusProvider::getChannelsWithStatus
  (PMTChannelStatusValue status) const
{
  ChannelSet_t result;
  for (auto const& [channel, info] : fDatabaseStatus)
    if (info.status == status) result.insert(channel);
  return result;
}


// -----------------------------------------------------------------------------
