/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapPostGres_service.cc
 * @brief  Wrapper service for `icarusDB::ICARUSChannelMapPostGresProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapPostGresProvider.h"
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMapPostGres; }
/**
 * @brief LArSoft service for ICARUS channel mapping (PostgreSQL DB backend).
 * 
 * This service provides access to ICARUS channel mapping database, using the
 * PostgreSQL database.
 * 
 * This service implements the generic channel mapping access service provider
 * interface `icarusDB::IICARUSChannelMapProvider`.
 * To use the channel mapping, include in your code the header of the service
 * interface `icarusDB::IICARUSChannelMap`, and access the service provider via:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarusDB::IICARUSChannelMapProvider const& channelMapping
 *   = *lar::providerFrom<icarusDB::IICARUSChannelMap>();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * or similar (`lar::providerFrom()` is in `larcore/CoreUtils/ServiceUtils.h`)
 * or directly the service via
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarusDB::IICARUSChannelMapProvider const& channelMapping
 *   = *art::ServiceHandle<icarusDB::IICARUSChannelMap>();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * For details on the interface, see `icarusDB::IICARUSChannelMapProvider`.
 * For details on the implementation, see
 * `icarusDB::ICARUSChannelMapPostGresProvider`.
 * 
 * @note For production and processes that need robustness, using the SQLite
 *       database and backend (`icarusDB::ICARUSChannelMapSQLite`) is suggested
 *       instead. That database is local to the process and not subject to
 *       network interruptions.
 * 
 */
class icarusDB::ICARUSChannelMapPostGres
  : public IICARUSChannelMap, public ICARUSChannelMapPostGresProvider,
{
  
  /// Prepares the mapping for the specified run.
  void preBeginRun(art::Run const& run);
  
    public:
  
  using provider_type = icarusDB::ICARUSChannelMapPostGresProvider;
  
  using Parameters = art::ServiceTable<provider_type::Config>;
  
  /// Constructor: configures the provider and hooks to the framework.
  ICARUSChannelMapPostGres(Parameters const& params, art::ActivityRegistry& reg);
  
  /// Returns the service provider (for use with `lar::providerFrom()`).
  provider_type const* provider() const { return this; }
  
}; // class icarusDB::ICARUSChannelMapPostGres


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarusDB::ICARUSChannelMapPostGres::ICARUSChannelMapPostGres
  (Parameters const& params, art::ActivityRegistry& reg)
  : provider_type(params())
{
  reg.sPreBeginRun.watch(this, &ICARUSChannelMapPostGres::preBeginRun);
  forPeriod(RunPeriod::Runs0to2); // prepare for some run, in case anybody asks
}


// -----------------------------------------------------------------------------
void icarusDB::ICARUSChannelMapPostGres::preBeginRun(art::Run const& run) {
  if (forRun(run.run())) {
    mf::LogDebug{ "ICARUSChannelMapPostGres" }
      << "Loaded mapping for run " << run.run();
  }
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL
  (icarusDB::ICARUSChannelMapPostGres, icarusDB::IICARUSChannelMap, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL
  (icarusDB::ICARUSChannelMapPostGres, icarusDB::IICARUSChannelMap)


// -----------------------------------------------------------------------------
