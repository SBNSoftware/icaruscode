/**
 * @file   icaruscode/PMT/Status/PMTChannelStatusService_service.cc
 * @brief  Art service for PMT channel status from the conditions database.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 */

#include "icaruscode/PMT/Status/IPMTChannelStatusService.h"
#include "icaruscode/PMT/Status/PMTChannelStatusProvider.h"

// Framework libraries
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class PMTChannelStatusService; }
class icarusDB::PMTChannelStatusService
  : public IPMTChannelStatusService, private PMTChannelStatusProvider {

    void preBeginRun(const art::Run& run);

    /// Returns a constant pointer to the service provider.
    virtual PMTChannelStatusProvider const* do_provider() const override
      { return this; }

  public:

    PMTChannelStatusService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

}; // class icarusDB::PMTChannelStatusService


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarusDB::PMTChannelStatusService::PMTChannelStatusService
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& reg)
  : PMTChannelStatusProvider(pset)
{
  reg.sPreBeginRun.watch(this, &PMTChannelStatusService::preBeginRun);
}


// -----------------------------------------------------------------------------
void icarusDB::PMTChannelStatusService::preBeginRun(const art::Run& run)
{
  readStatusFromDB(run.run());
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::PMTChannelStatusService, icarusDB::IPMTChannelStatusService, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(icarusDB::PMTChannelStatusService, icarusDB::IPMTChannelStatusService)


// -----------------------------------------------------------------------------
