////////////////////////////////////////////////////////////////////////
// \file PhotonCalibratorServiceStandard.h
//
// \brief Framework interface to PhotonCalibratorStandard
//
// \author micarrig@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef PHOTONCALIBRATORSERVICEFROMDB
#define PHOTONCALIBRATORSERVICEFROMDB

// LArSoft Includes
#include "larreco/Calibrator/IPhotonCalibratorService.h"
#include "icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"

namespace calib {

  class PhotonCalibratorServiceFromDB : public IPhotonCalibratorService {
  public:
    using provider_type = icarusDB::PhotonCalibratorFromDB;

    struct ServiceConfiguration_t {
      fhicl::Atom<float> SPESize{fhicl::Name("SPESize")};
      fhicl::Atom<float> SPEShift{fhicl::Name("SPEShift")};
      fhicl::Atom<bool> UseArea{fhicl::Name("UseArea")};
      fhicl::Atom<bool>  Verbose{fhicl::Name("Verbose"), false};
      fhicl::Atom<std::string> LogCategory{fhicl::Name("LogCategory"), "PhotonCalibratorServiceFromDB"};
      fhicl::Atom<std::string> AreaTag{fhicl::Name("AreaTag"), ""};
    };

    using Parameters = art::ServiceTable<ServiceConfiguration_t>;

    PhotonCalibratorServiceFromDB(Parameters const& pset, art::ActivityRegistry& reg);

  private:
    provider_type const* provider() const override { return &fProvider; }

    void preBeginRun(art::Run const& run);

    icarusDB::PhotonCalibratorFromDB fProvider;

    bool fVerbose = false; ///< Whether to print the configuration we read.
    std::string fLogCategory; ///< Category tag for messages.
  };

}

DECLARE_ART_SERVICE_INTERFACE_IMPL(calib::PhotonCalibratorServiceFromDB,
                                   calib::IPhotonCalibratorService,
                                   LEGACY)

#endif // PHOTONCALIBRATORSERVICEFROMDB