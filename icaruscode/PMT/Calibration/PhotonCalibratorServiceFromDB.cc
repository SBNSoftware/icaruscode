#include "icaruscode/PMT/Calibration/PhotonCalibratorServiceFromDB.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace calib {

  PhotonCalibratorServiceFromDB::PhotonCalibratorServiceFromDB(
    Parameters const& pset, 
    art::ActivityRegistry& reg
  )
    : fProvider{ pset.get_PSet() }
  {}

  void PhotonCalibratorServiceFromDB::preBeginRun(art::Run const& run)
  {
    fProvider.readCalibrationFromDB(run.run());
  }

} // namespace calib

DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceFromDB,
                                  calib::IPhotonCalibratorService)