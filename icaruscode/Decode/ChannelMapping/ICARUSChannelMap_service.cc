/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMap_service.cc
 * @brief  Wrapper service for `icarusDB::ICARUSChannelMapProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h"
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMap; }
class icarusDB::ICARUSChannelMap: public ICARUSChannelMapProvider {
  
  /// Prepares the mapping for the specified run.
  void preBeginRun(art::Run const& run);
  
    public:
  
  ICARUSChannelMap(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class icarusDB::ICARUSChannelMap


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarusDB::ICARUSChannelMap::ICARUSChannelMap
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& reg)
  : ICARUSChannelMapProvider(pset)
{
  reg.sPreBeginRun.watch(this, &ICARUSChannelMap::preBeginRun);
  forPeriod(RunPeriod::Runs0to2); // prepare for some run, in case anybody asks
}


// -----------------------------------------------------------------------------
void icarusDB::ICARUSChannelMap::preBeginRun(art::Run const& run) {
  if (forRun(run.run())) {
    mf::LogDebug{ "ICARUSChannelMap" }
      << "Loaded mapping for run " << run.run();
  }
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap)


// -----------------------------------------------------------------------------
