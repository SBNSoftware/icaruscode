#include "icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h"

namespace icarusDB {

  PhotonCalibratorFromDB::PhotonCalibratorFromDB(const fhicl::ParameterSet& pset)
    : fVerbose   ( pset.get<bool>("Verbose", false) )
    , fLogCategory( pset.get<std::string>("LogCategory", "PhotonCalibratorFromDB") )
  {
    fAreaTag  = pset.get<std::string>("AreaTag", "");
  }

  void PhotonCalibratorFromDB::readCalibrationFromDB(unsigned int run)
  {

    mf::LogInfo(fLogCategory) << "Reading SPE area calibrations from database for run " << run;

    const std::string dbname("pmt_speareas_data");
    lariov::DBFolder db(dbname, "", "", fAreaTag, true, false);

    mf::LogDebug(fLogCategory) << "Connecting to " << dbname << " folder";

    bool ret = db.UpdateData( RunToDatabaseTimestamp(run) ); // select table based on run number   
    mf::LogDebug(fLogCategory) << dbname + " corrections" << (ret? "": " not") << " updated for run " << run;
    mf::LogTrace(fLogCategory)
           << "Fetched IoV [ " << db.CachedStart().DBStamp() << " ; " << db.CachedEnd().DBStamp()
           << " ] to cover t=" << RunToDatabaseTimestamp(run)
           << " [=" << lariov::TimeStampDecoder::DecodeTimeStamp(RunToDatabaseTimestamp(run)).DBStamp() << "]";

    std::vector<unsigned int> channelList;
    if (int res = db.GetChannelList(channelList); res != 0) {
      throw cet::exception
        ( "PhotonCalibratorFromDB" ) << "GetChannelList() returned " << res << " on run " << run << " query in " << dbname << "\n";
    }
    
    if (channelList.empty()) {
      throw cet::exception("PhotonCalibratorFromDB") << "got an empty channel list for run " << run << " in " << dbname << "\n";
    }

    for( auto channel : channelList ) {
        
        // Laser correction
        double area = 0;
        int error  = db.GetNamedChannelData( channel, "area", area );
        if( error ) throw cet::exception( "PhotonCalibratorFromDB" ) << "Encountered error (code " << error << ") while trying to access 'area' on table " << dbname << "\n";

        fDatabaseSPECalibrations[channel].speArea = area; 
    }  
  }

  double PhotonCalibratorFromDB::PE(double adcs, int channel) const
  {
    return adcs / getChannelCorrOrDefault(channel).speArea;
  }

  bool PhotonCalibratorFromDB::UseArea() const
  {
    return false;
  }

  uint64_t PhotonCalibratorFromDB::RunToDatabaseTimestamp( uint32_t run ) const{

    // Run number to timestamp used in the db
    // DBFolder.h only takes 19 digit (= timestamp in nano second),
    // but ICARUS tables are currently using run numbers
    // Step 1) Add 1000000000 to the run number; e.g., run XXXXX -> 10000XXXXX
    // Step 2) Multiply 1000000000
    uint64_t runNum = uint64_t(run);
    uint64_t timestamp = runNum+1000000000;
    timestamp *= 1000000000;

    if( fVerbose ) mf::LogInfo(fLogCategory) << "Run " << runNum << " corrections from DB timestamp " << timestamp;
    
    return timestamp;
  } 

} // namespace icarusDB