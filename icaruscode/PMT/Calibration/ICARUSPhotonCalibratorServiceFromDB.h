/**
 * @file   icaruscode/PMT/Calibration/ICARUSPhotonCalibratorServiceFromDB.h
 * @brief  Framework service interface to `icarusDB::PhotonCalibratorFromDB`.
 * @author Michael Carrigan (micarrig@fnal.gov)
 * @see    icaruscode/PMT/Calibration/ICARUSPhotonCalibratorServiceFromDB_service.cc
 */

#ifndef ICARUSCODE_PMT_CALIBRATION_ICARUSPHOTONCALIBRATORSERVICEFROMDB_H
#define ICARUSCODE_PMT_CALIBRATION_ICARUSPHOTONCALIBRATORSERVICEFROMDB_H

// LArSoft Includes
#include "larreco/Calibrator/IPhotonCalibratorService.h"
#include "icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "fhiclcpp/types/TableFragment.h"

namespace calib {

  /**
   * @brief Photoelectron calibration service using ICARUS database.
   * @see `icarusDB::PhotonCalibratorFromDB`
   * 
   * The service reads the information on the translation from a reconstructed
   * optical hit to photoelectrons from an external database.
   * 
   * This is a _art_ framework service, which only delivers a
   * framework-independent service provider (`icarusDB::PhotonCalibratorFromDB`)
   * that does all the lifting.
   * The provider can be summoned in an _art_ context for example by
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * `calib::IPhotonCalibrator const& calibrator
   *   = *(lar::providerFrom<calib::IPhotonCalibratorService>());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`#include "larcore/CoreUtils/ServiceUtil.h"` for `lar::providerFrom()`).
   * 
   * The service triggers an update of the provider caches on every change of
   * run. This allows the service to be thread-safe in _art_ terms (because
   * _art_ threads do not mix runs).
   * 
   */
  class ICARUSPhotonCalibratorServiceFromDB : public IPhotonCalibratorService {
  public:
    using provider_type = icarusDB::PhotonCalibratorFromDB;

    struct ServiceConfiguration_t {
      
      // this part of the configuration is passed to the provider:
      fhicl::TableFragment<provider_type::Config> providerConfig;
      
    };

    using Parameters = art::ServiceTable<ServiceConfiguration_t>;

    ICARUSPhotonCalibratorServiceFromDB(Parameters const& params, art::ActivityRegistry& reg);

  private:
    provider_type const* provider() const override { return &fProvider; }

    void preBeginRun(art::Run const& run);

    icarusDB::PhotonCalibratorFromDB fProvider;

  };

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::ICARUSPhotonCalibratorServiceFromDB,
                                   calib::IPhotonCalibratorService,
                                   SHARED)

#endif // ICARUSCODE_PMT_CALIBRATION_ICARUSPHOTONCALIBRATORSERVICEFROMDB_H
