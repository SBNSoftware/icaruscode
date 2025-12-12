#include "icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h"

namespace icarusDB {

  PhotonCalibratorFromDB::PhotonCalibratorFromDB(const fhicl::ParameterSet& pset)
    : fVerbose   ( pset.get<bool>("Verbose", false) )
    , fLogCategory( pset.get<std::string>("LogCategory", "PhotonCalibratorFromDB") )
  {}

  void PhotonCalibratorFromDB::readCalibrationFromDB(unsigned int run)
  {

    return;
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