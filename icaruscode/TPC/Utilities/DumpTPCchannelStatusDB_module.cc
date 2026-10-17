/**
 * @file   icaruscode/TPC/Utilities/DumpTPCchannelStatusDB_module.cc
 * @brief  Provides the `DumpTPCchannelStatusDB` module.
 * @author Gianluca Petrillo (petrillo@slac.stanfird.edu)
 * @date   March 11, 2025
 * 
 */

// SBN/ICARUS libraries
#include "icaruscode/TPC/Utilities/TPCSIOVchannelStatusDBdumper.h"

// LArSoft libraries
#include "larevt/CalibrationDBI/Providers/SIOVChannelStatusProvider.h"
#include "larevt/CalibrationDBI/IOVData/IOVDataError.h"
#include "larevt/CalibrationDBI/Interface/CalibrationDBIFwd.h" // lariov::DBTimeStamp_t

#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/CoreUtils/counter.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/Exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalTable.h"


// C/C++ standard libraries
#include <iomanip> // std::setw()
#include <string>
#include <vector>


/**
 * @brief Dumps on screen the state of all channels at specified times.
 * 
 * Dump the content of the channel status database interfaced with 
 * `SIOVChannelStatusProvider`, for a series of timestamps specified in the
 * configuration.
 * 
 * The dump format is a header in the form:
 * ```
 * Status of XXXX channels for N timestamps:
 * ```
 * (`XXX` the total number of channels in the geometry, N the number of
 * timestamps in the configuration) and then, for each timestamp:
 * ```
 * === BEGIN TIMESTAMP: xxxxxxxxxxxxxxxxxxx =======================================
 * CH=0: status
 * CH=1: status
 * [...]
 * 
 * Counting B BAD channels: BadCH BadCH [...]
 * Counting N NOISY channels: NoisyCH NoisyCh [...]
 * Counting G good channels.
 * === END   TIMESTAMP: xxxxxxxxxxxxxxxxxxx =======================================
 * ```
 * 
 * 
 * Timestamp format
 * -----------------
 * 
 * The timestamp format is really defined by the database conventions.
 * In ICARUS database, for example, that is in nanoseconds and in UTC.
 * 
 * 
 * Dependencies
 * -------------
 * 
 * The following services must be configured:
 *  * `WireReadout` (to discover the number of TPC channels in the detector).
 * 
 * @note This module does not use the `ChannelStatusService` configured
 *       in the _art_ job. The reason is that it does not offer any control on
 *       the timestamp, but rather it binds it to the one of the input event.
 *       Here we need control on the timestamp, which would complicate a lot
 *       the workflow (we would need to generate empty events with specific
 *       timestamps, which should be delegated to `EmptyEvent` with a custom
 *       timestamp plugin where the timestamps should be inserted).
 * 
 * 
 * Configuration
 * --------------
 * 
 * * `Timestamps` (sequence of integers, default: [ `2000000000000000000` ]):
 *    sequence of timestamps to query about. To dump the latest settings, a safe
 *    timestamp is a large one like `2'000'000'000'000'000'000` (2e18, which
 *    is also the only timestamp in the default value).
 * * `PrintIndividualChannels` (flag, default: `false`): if set, one line will
 *    be printed for each channel, showing its status.
 * * `PrintSummary` (flag, default: `true`): if set, the list of noisy and of
 *    bad channels will be printed, together with a count of good channels.
 * * `SIOVChannelStatusProviderConfig` (parameter table): configuration of the
 *    service provider of type `lariov::SIOVChannelStatusProvider` used to
 *    access the database. It can be copied from the official experiment one
 *    if the experiment normally uses that provider.
 * * `LogCategory` (string, default: `DumpTPCchannelStatusDB`): messagefacility
 *    category stream where the information is written (`mf::LogInfo` level).
 * 
 */
class DumpTPCchannelStatusDB: public art::EDAnalyzer {
  
    public:
  
  using DBTimeStamp_t = lariov::DBTimeStamp_t;
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<DBTimeStamp_t> Timestamps {
      Name{ "Timestamps" },
      Comment{ "timestamps to query about" },
      std::vector{ DBTimeStamp_t{ 2'000'000'000'000'000'000 } }
      };
    
    fhicl::OptionalTable<lariov::TPCSIOVchannelStatusDBdumper::Config> DumperSettings {
      Name{ "DumperSettings" },
      Comment{ "configuration of the dumping algorithm" }
      };
    
    fhicl::DelegatedParameter SIOVChannelStatusProviderConfig {
      Name{ "ProviderConfiguration" },
      Comment{ "full configuration of a service provider" }
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of the messagefacility category to send output to" },
      "DumpTPCchannelStatusDB"
      };
    
  }; // Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  
  DumpTPCchannelStatusDB(Parameters const& params);
  
  
  /// Dumps the content of the channel status database.
  virtual void beginJob() override;
  
  /// Does nothing (but it is required).
  virtual void analyze(art::Event const&) override {}
  
  
    private:
  
  // --- BEGIN ---  Configuration parameters  ----------------------------------
  
  std::vector<DBTimeStamp_t> const fTimestamps; ///< Timestamps to query.
  
  std::string const fLogCategory; ///< Messagefacility category to write into.
  
  // --- END -----  Configuration parameters  ----------------------------------
  
  // --- BEGIN ---  Cached service information  --------------------------------
  
  /// Channel status service provider.
  lariov::SIOVChannelStatusProvider fChannelStatus;
  
  // --- END   ---  Cached service information  --------------------------------
  
  lariov::TPCSIOVchannelStatusDBdumper fDumper; ///< The dumping algorithm.
  
}; // class DumpTPCchannelStatusDB


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
DumpTPCchannelStatusDB::DumpTPCchannelStatusDB(Parameters const& params)
  : art::EDAnalyzer{ params }
  , fTimestamps{ params().Timestamps() }
  , fLogCategory{ params().LogCategory() }
  , fChannelStatus
    { params().SIOVChannelStatusProviderConfig.get<fhicl::ParameterSet>() }
  , fDumper{
      params().DumperSettings().value_or(lariov::TPCSIOVchannelStatusDBdumper::Config{}),
      fChannelStatus,
      art::ServiceHandle<geo::WireReadout>()->Get().Nchannels()
    }
{
}


// -----------------------------------------------------------------------------
void DumpTPCchannelStatusDB::beginJob() {
  
  mf::LogVerbatim{ fLogCategory }
    << "Status of " << fDumper.nChannels() << " channels for "
    << fTimestamps.size() << " timestamps:";
  
  for (DBTimeStamp_t const timestamp: fTimestamps) {
    mf::LogVerbatim{ fLogCategory } << fDumper.timestampToStream(timestamp);
  }
  
} // DumpTPCchannelStatusDB::beginJob()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(DumpTPCchannelStatusDB)


// -----------------------------------------------------------------------------
