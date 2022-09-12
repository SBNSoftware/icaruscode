
////////////////////////////////////////////////////////////////////////
// \file SpaceChargeService.h
//
// \brief pure virtual service interface for timing corrections
//
// \author ascarpell@bnl.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGESERVICE_H
#define SPACECHARGESERVICE_H

#include "larevt/SpaceCharge/SpaceCharge.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace spacecharge {
  class SpaceChargeService {
  public:
    using provider_type = spacecharge::SpaceCharge;

    virtual ~SpaceChargeService() = default;
    virtual const spacecharge::SpaceCharge* provider() const = 0;
  };
}

DECLARE_ART_SERVICE_INTERFACE(spacecharge::SpaceChargeService, SHARED)

#endif // SPACECHARGESERVICE_H