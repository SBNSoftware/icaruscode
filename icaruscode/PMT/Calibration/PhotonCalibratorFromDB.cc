#include "icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h"

namespace calib {

  PhotonCalibratorFromDB::PhotonCalibratorFromDB(const fhicl::ParameterSet& pset)
    : fVerbose   ( pset.get<bool>("Verbose", false) )
    , fLogCategory( pset.get<std::string>("LogCategory", "PhotonCalibratorFromDB") )
  {}

  void PhotonCalibratorFromDB::readCalibrationFromDB(unsigned int run)
  {

    return;
  }

  double PhotonCalibratorFromDB::PE(double adcs, int opchannel) const
  {
    return 0.0;
  }

  bool PhotonCalibratorFromDB::UseArea() const
  {
    return false;
  }

} // namespace calib