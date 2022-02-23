/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMap_service.cc
 * @brief  Wrapper service for `icarusDB::ICARUSChannelMapProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#include "icaruscode/TPC/Simulation/DetSim/tools/CoherentNoiseFactorProvider.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace Noise { class CoherentNoiseFactor; }
class Noise::CoherentNoiseFactor: public CoherentNoiseFactorProvider {
    
    public:
  
  CoherentNoiseFactor(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class icarusDB::ICARUSChannelMap


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
Noise::CoherentNoiseFactor::CoherentNoiseFactor
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& /* reg */)
  : CoherentNoiseFactorProvider(pset)
  {}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(Noise::CoherentNoiseFactor, Noise::ICoherentNoiseFactor, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(Noise::CoherentNoiseFactor, Noise::ICoherentNoiseFactor)


// -----------------------------------------------------------------------------
