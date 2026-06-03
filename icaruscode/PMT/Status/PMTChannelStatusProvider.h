/**
 * @file   icaruscode/PMT/Status/PMTChannelStatusProvider.h
 * @brief  Provider for PMT channel status from the conditions database.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @see    icaruscode/PMT/Status/PMTChannelStatusProvider.cxx
 */

#ifndef ICARUSCODE_PMT_STATUS_PMTCHANNELSTATUSPROVIDER_H
#define ICARUSCODE_PMT_STATUS_PMTCHANNELSTATUSPROVIDER_H

// ICARUS libraries
#include "icaruscode/PMT/Status/PMTChannelStatus.h"

// LArSoft libraries
#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// Framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// C++ standard libraries
#include <cstdint>
#include <map>
#include <string>


namespace icarusDB::details {

  /// Internal structure holding the status and voltage for a single channel.
  struct PMTChannelStatusDB {
    PMTChannelStatusValue status  = kOff;
    double                voltage = 0.0;  ///< Nominal voltage [V].
  };

} // namespace icarusDB::details


// -----------------------------------------------------------------------------
namespace icarusDB { class PMTChannelStatusProvider; }
/**
 * @brief Provider for PMT channel status from the conditions database.
 *
 * Reads the channel status (kOff=0, kGood=1, kBad=2) from the
 * `pmt_voltage_data` database table, indexed by run number.
 *
 * The database interface is accessed only on `readStatusFromDB()` calls,
 * and the relevant information is cached.
 *
 * Configuration parameters
 * -------------------------
 * * `DBname` (default: `"pmt_voltage_data"`): database table name.
 * * `StatusTag` (default: `""`): database tag to select.
 * * `DefaultStatus` (default: `0` = kOff): status for channels absent from DB.
 * * `Verbose` (default: `false`): print channel statuses when loading.
 * * `LogCategory` (default: `"PMTChannelStatusProvider"`).
 *
 */
class icarusDB::PMTChannelStatusProvider : public PMTChannelStatus {

public:

  /// Constructor from FHiCL configuration (no database access yet).
  PMTChannelStatusProvider(const fhicl::ParameterSet& pset);

  /// Loads channel statuses for the given run from the database.
  void readStatusFromDB(unsigned int run);

  /// Converts a run number to the internal 19-digit nanosecond timestamp.
  std::uint64_t RunToDatabaseTimestamp(unsigned int run) const;

  // --- PMTChannelStatus interface ---

  PMTChannelStatusValue getChannelStatus(unsigned int channel) const override
    { return getChannelStatusOrDefault(channel).status; }

  ChannelSet_t getOnChannels()  const override { return getChannelsWithStatus(kGood);  }
  ChannelSet_t getOffChannels() const override { return getChannelsWithStatus(kOff); }
  ChannelSet_t getBadChannels() const override { return getChannelsWithStatus(kBad); }

  double getChannelVoltage(unsigned int channel) const override
    { return getChannelStatusOrDefault(channel).voltage; }

  std::string getDatabaseTag() const override { return fStatusTag; }

private:

  using PMTChannelStatusDB = details::PMTChannelStatusDB;

  bool        fVerbose;
  std::string fLogCategory;
  std::string fStatusTag;
  int fOverrideRunNumber; ///< If non-negative, overrides the run number for DB queries.

  PMTChannelStatusDB fDefault; ///< Status used for channels not present in DB.

  lariov::DBFolder fDB; ///< Cached database interface.

  /// Map of channel ID to status information.
  std::map<unsigned int, PMTChannelStatusDB> fDatabaseStatus;

  /// Returns the channel status record, or the default if channel is absent.
  PMTChannelStatusDB const& getChannelStatusOrDefault(unsigned int channel) const {
    auto const it = fDatabaseStatus.find(channel);
    return (it == fDatabaseStatus.end()) ? fDefault : it->second;
  }

  /// Returns the set of all channels with the given status.
  ChannelSet_t getChannelsWithStatus(PMTChannelStatusValue status) const;

}; // class icarusDB::PMTChannelStatusProvider


#endif // ICARUSCODE_PMT_STATUS_PMTCHANNELSTATUSPROVIDER_H
