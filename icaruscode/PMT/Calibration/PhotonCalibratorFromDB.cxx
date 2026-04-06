/**
 * @file   icaruscode/PMT/Calibration/PhotonCalibratorFromDB.cxx
 * @brief  Implementation of optical hit photoelectron calibration from database.
 * @author Michael Carrigan, Matteo Vicenzi
 * @see    icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h
 */

#include "icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h"

// Database interface helpers
#include "larevt/CalibrationDBI/IOVData/TimeStampDecoder.h"

// Message facility
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard library
#include <iomanip> // std::setw()


// -----------------------------------------------------------------------------
namespace icarusDB {
  
  details::PhotonCalibratorInfo convert
    (PhotonCalibratorFromDB::Config::DefaultCalib const& config)
  {
    details::PhotonCalibratorInfo info;
    info.speArea = config.SPEArea();
    // all the rest is left to C++ default
    return info;
  }
  
}


// -----------------------------------------------------------------------------
icarusDB::PhotonCalibratorFromDB::PhotonCalibratorFromDB(const Config& config)
  : fCalibDefaults     ( config.Defaults() )
  , fVerbose           ( config.Verbose() )
  , fLogCategory       ( config.LogCategory() )
  , fAreaTag           ( config.AreaTag() )
  , fOverrideRunNumber ( config.OverrideRunNumber() )
  , fDB                ( config.DBname(), "", "", config.AreaTag(), true, false)
{
  mf::LogInfo(fLogCategory)
    << "PhotonCalibratorFromDB connected to " << config.DBname() << " DB, tag '" << config.AreaTag() << "'";
}


// -----------------------------------------------------------------------------
void icarusDB::PhotonCalibratorFromDB::readCalibrationFromDB(unsigned int run)
{

  if (fOverrideRunNumber >= 0) {
    mf::LogInfo(fLogCategory)
      << "Overriding run number " << run << " with " << fOverrideRunNumber
      << " for DB queries.";
    run = static_cast<unsigned int>(fOverrideRunNumber);
  }

  mf::LogInfo(fLogCategory) << "Reading SPE area calibrations from database for run " << run;

  bool ret = fDB.UpdateData( RunToDatabaseTimestamp(run) ); // select table based on run number
  mf::LogTrace(fLogCategory)
         << "Calibrations" << (ret? "": " not") << " updated for run " << run
         << "\nFetched IoV [ " << fDB.CachedStart().DBStamp() << " ; " << fDB.CachedEnd().DBStamp()
         << " ] to cover t=" << RunToDatabaseTimestamp(run)
         << " [=" << lariov::TimeStampDecoder::DecodeTimeStamp(RunToDatabaseTimestamp(run)).DBStamp() << "]";

  std::vector<unsigned int> channelList;
  if (int res = fDB.GetChannelList(channelList); res != 0) {
    throw cet::exception
      ( "PhotonCalibratorFromDB" ) << "GetChannelList() returned " << res << " on run " << run << " query\n";
  }
  
  if (channelList.empty()) {
    throw cet::exception("PhotonCalibratorFromDB") << "got an empty channel list for run " << run << "\n";
  }

  mf::LogDebug log(fLogCategory);
  log << "Loading calibration for " << channelList.size() << " channels in run " << run;
  for( auto channel : channelList ) {

    // SPE area [ADC x tick]
    double area = 0;
    int error  = fDB.GetNamedChannelData( channel, "area", area );
    if( error ) throw cet::exception( fLogCategory ) << "Error (code " << error << ") reading 'area' for channel " << channel << "\n";
    fDatabaseSPECalibrations[channel].speArea = area; // ADC x tick

    // SPE area uncertainty from fit
    double areaErr = 0;
    error = fDB.GetNamedChannelData( channel, "area_err", areaErr );
    if( error ) throw cet::exception( fLogCategory ) << "Error (code " << error << ") reading 'area_err' for channel " << channel << "\n";
    fDatabaseSPECalibrations[channel].speAreaErr = areaErr; // ADC x tick

    // SPE fit width (sigma) [ADC x tick]
    double width = 0;
    error = fDB.GetNamedChannelData( channel, "width", width );
    if( error ) throw cet::exception( fLogCategory ) << "Error (code " << error << ") reading 'width' for channel " << channel << "\n";
    fDatabaseSPECalibrations[channel].speFitWidth = width; // ADC x tick

    // SPE fit width uncertainty
    double widthErr = 0;
    error = fDB.GetNamedChannelData( channel, "width_err", widthErr );
    if( error ) throw cet::exception( fLogCategory ) << "Error (code " << error << ") reading 'width_err' for channel " << channel << "\n";
    fDatabaseSPECalibrations[channel].speFitWidthErr = widthErr; // ADC x tick

    // Fit quality
    double fitQuality = 0;
    error = fDB.GetNamedChannelData( channel, "fit_quality", fitQuality );
    if( error ) throw cet::exception( fLogCategory ) << "Error (code " << error << ") reading 'fit_quality' for channel " << channel << "\n";
    fDatabaseSPECalibrations[channel].speFitQuality = fitQuality; 

    if (fVerbose)
      log << "\n  channel " << std::setw(3) << channel
          << "  SPE area " << area << " +/- " << areaErr
          << "  sigma " << width << " +/- " << widthErr
          << "  fit quality " << fitQuality;
  }
}


// -----------------------------------------------------------------------------
double icarusDB::PhotonCalibratorFromDB::PE(double adcs, int channel) const
{
  return adcs / getChannelCalibOrDefault(channel).speArea;
}


// -----------------------------------------------------------------------------
bool icarusDB::PhotonCalibratorFromDB::UseArea() const
{
  return true;
}


// -----------------------------------------------------------------------------
uint64_t icarusDB::PhotonCalibratorFromDB::RunToDatabaseTimestamp( unsigned int run ) const{

  // Run number to timestamp used in the db
  // DBFolder.h only takes 19 digit (= timestamp in nano second),
  // but ICARUS tables are currently using run numbers
  // Step 1) Add 1000000000 to the run number; e.g., run XXXXX -> 10000XXXXX
  // Step 2) Multiply 1000000000
  uint64_t runNum = uint64_t(run);
  uint64_t timestamp = runNum+1'000'000'000;
  timestamp *= 1'000'000'000;

  if( fVerbose ) mf::LogInfo(fLogCategory) << "Run " << runNum << " calibrations from DB timestamp " << timestamp;
  
  return timestamp;
} 

