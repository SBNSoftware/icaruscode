/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMap_service.cc
 * @brief  Wrapper service for `icarusDB::ICARUSChannelMapProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMap; }
class icarusDB::ICARUSChannelMap: public ICARUSChannelMapProvider {
    
    public:
  
  ICARUSChannelMap(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class icarusDB::ICARUSChannelMap


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarusDB::ICARUSChannelMap::ICARUSChannelMap
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& /* reg */)
  : ICARUSChannelMapProvider(pset)
  {}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap)


// -----------------------------------------------------------------------------
