/**
 * @file   icaruscode/PMT/Calibration/ICARUSPhotonCalibratorServiceFromDB_service.cc
 * @brief  Framework service interface to `icarusDB::PhotonCalibratorFromDB`.
 * @author Michael Carrigan (micarrig@fnal.gov)
 * @see    icaruscode/PMT/Calibration/ICARUSPhotonCalibratorServiceFromDB.h
 */

#include "icaruscode/PMT/Calibration/ICARUSPhotonCalibratorServiceFromDB.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"


calib::ICARUSPhotonCalibratorServiceFromDB::ICARUSPhotonCalibratorServiceFromDB(
  Parameters const& params,
  art::ActivityRegistry& reg
)
  : fProvider{ params().providerConfig() }
{
  reg.sPreBeginRun.watch(this, &ICARUSPhotonCalibratorServiceFromDB::preBeginRun);
}

void calib::ICARUSPhotonCalibratorServiceFromDB::preBeginRun(art::Run const& run)
{
  fProvider.readCalibrationFromDB(run.run());
}


DEFINE_ART_SERVICE_INTERFACE_IMPL(calib::ICARUSPhotonCalibratorServiceFromDB,
                                  calib::IPhotonCalibratorService)
