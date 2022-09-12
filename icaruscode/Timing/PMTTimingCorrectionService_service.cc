/**
 * @file   icaruscode/Timing/PMTTimingCorrectionService_service.cc
 * @brief  Wrapper class for 'PMTTimingCorrectionsProvider.h'
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 */

#include "icaruscode/Timing/PMTTimingCorrectionsProvider.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class PMTTimingCorrectionService; }
class icarusDB::PMTTimingCorrectionService: public PMTTimingCorrectionsProvider {
    
    public:
  
      PMTTimingCorrectionService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class icarusDB::ICARUSChannelMap


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarusDB::PMTTimingCorrectionService::PMTTimingCorrectionService
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& /* reg */)
  : PMTTimingCorrectionsProvider(pset)
  {}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::PMTTimingCorrectionService, icarusDB::PMTTimingCorrections, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(icarusDB::PMTTimingCorrectionService, icarusDB::PMTTimingCorrections)


// -----------------------------------------------------------------------------