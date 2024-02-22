/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMap_service.cc
 * @brief  Wrapper service for `icarusDB::ICARUSChannelMapProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapSQLiteProvider.h"
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMapSQLite; }
class icarusDB::ICARUSChannelMapSQLite: public ICARUSChannelMapSQLiteProvider {
  
  /// Prepares the mapping for the specified run.
  void preBeginRun(art::Run const& run);
  
    public:
  
  using provider_type = icarusDB::ICARUSChannelMapSQLiteProvider;
  
  using Parameters = art::ServiceTable<provider_type::Config>;
  
  /// Constructor: configures the provider and hooks to the framework.
  ICARUSChannelMapSQLite(Parameters const& params, art::ActivityRegistry& reg);
  
  /// Returns the service provider (for use with `lar::providerFrom()`).
  provider_type const* provider() const { return this; }
  
}; // class icarusDB::ICARUSChannelMapSQLite


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarusDB::ICARUSChannelMapSQLite::ICARUSChannelMapSQLite
  (Parameters const& params, art::ActivityRegistry& reg)
  : provider_type(params())
{
  reg.sPreBeginRun.watch(this, &ICARUSChannelMapSQLite::preBeginRun);
  forPeriod(RunPeriod::Runs0to2); // prepare for some run, in case anybody asks
}


// -----------------------------------------------------------------------------
void icarusDB::ICARUSChannelMapSQLite::preBeginRun(art::Run const& run) {
  if (forRun(run.run())) {
    mf::LogDebug{ "ICARUSChannelMapSQLite" }
      << "Loaded mapping for run " << run.run();
  }
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL
  (icarusDB::ICARUSChannelMapSQLite, icarusDB::IICARUSChannelMap, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL
  (icarusDB::ICARUSChannelMapSQLite, icarusDB::IICARUSChannelMap)


// -----------------------------------------------------------------------------
