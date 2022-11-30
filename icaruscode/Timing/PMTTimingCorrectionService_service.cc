/**
 * @file   icaruscode/Timing/PMTTimingCorrectionService_service.cc
 * @brief  Wrapper class for 'PMTTimingCorrectionsProvider.h'
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 */

#include "icaruscode/Timing/IPMTTimingCorrectionService.h"
#include "icaruscode/Timing/PMTTimingCorrectionsProvider.h"

// framework libraries
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class PMTTimingCorrectionService; }
class icarusDB::PMTTimingCorrectionService
  : public IPMTTimingCorrectionService, private PMTTimingCorrectionsProvider {
    
      void preBeginRun(const art::Run& run);
    
      /// Returns a constant pointer to the service provider
      virtual PMTTimingCorrectionsProvider const* do_provider() const override
         { return this; }
      
    public:
  
      PMTTimingCorrectionService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class icarusDB::PMTTimingCorrectionService


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarusDB::PMTTimingCorrectionService::PMTTimingCorrectionService
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& reg)
  : PMTTimingCorrectionsProvider(pset)
{
  reg.sPreBeginRun.watch(this, &PMTTimingCorrectionService::preBeginRun);
}


// -----------------------------------------------------------------------------
void icarusDB::PMTTimingCorrectionService::preBeginRun(const art::Run& run)
{
  readTimeCorrectionDatabase(run);
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::PMTTimingCorrectionService, icarusDB::IPMTTimingCorrectionService, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(icarusDB::PMTTimingCorrectionService, icarusDB::IPMTTimingCorrectionService)


// -----------------------------------------------------------------------------
