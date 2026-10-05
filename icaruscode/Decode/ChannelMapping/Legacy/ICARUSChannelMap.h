/**
 * @file   icaruscode/Decode/ChannelMapping/Legacy/ICARUSChannelMap.h
 * @brief  Wrapper service for `icarusDB::ICARUSChannelMapProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_LEGACY_ICARUSCHANNELMAP_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_LEGACY_ICARUSCHANNELMAP_H

#include "icaruscode/Decode/ChannelMapping/Legacy/ICARUSChannelMapProvider.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
// #include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMap; }
/**
 * @brief LArSoft service for ICARUS channel mapping (legacy implementation).
 * 
 * This service provides access to ICARUS channel mapping database, using a
 * backend in a _art_ tool chosen by the configuration.
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
 * `icarusDB::ICARUSChannelMapProvider`.
 * 
 * @note For production and processes that need robustness, using the SQLite
 *       database and backend (`icarusDB::ICARUSChannelMapSQLite`) is suggested
 *       instead. That database is local to the process and not subject to
 *       network interruptions.
 * 
 */
class icarusDB::ICARUSChannelMap: public ICARUSChannelMapProvider {
  
  /// Prepares the mapping for the specified run.
  void preBeginRun(art::Run const& run);
  
    public:
  
  ICARUSChannelMap(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
  
}; // class icarusDB::ICARUSChannelMap


// -----------------------------------------------------------------------------
// DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap, SHARED)


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_LEGACY_ICARUSCHANNELMAP_H
