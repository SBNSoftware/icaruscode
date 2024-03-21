/**
 * @file   icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h
 * @brief  Interface class for hardware/software channel mapping for ICARUS.
 * @author T. Usher
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAP_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAP_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"


namespace icarusDB { class IICARUSChannelMap; };

/**
 * @brief Interface for ICARUS channel mapping service.
 * @see icarusDB::IICARUSChannelMapProvider
 * 
 * This class represents the interface of the _art_ service to access ICARUS
 * channel mapping. The exposed interface is the same as the service provider
 * interface, `icarusDB::IICARUSChannelMapProvider`.
 * 
 * The only reason for the existence of this object is the overhead to make any
 * class a framework service.
 * 
 * To use the channel mapping, include in your code the header of this service
 * interface and access the service provider via:
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
 * The FHiCL configuration will route the service interface to the correct
 * implementation.
 * 
 * For details on the interface, see `icarusDB::IICARUSChannelMapProvider`.
 * For details on the implementation, see the chosen service implementation.
 * 
 * For use in environments other than _art_, use directly the service provider
 * implementation of your choice (any implementation of
 * `icarusDB::IICARUSChannelMapProvider`, for example
 * `icarusDB::IICARUSChannelMapProviderSQLite`).
 * 
 * 
 * Implementation details
 * -----------------------
 * 
 * This interface directly exposes the interface of the service provider,
 * deriving from it. It is foreseen that the implementations of this service
 * derive directly from their framework-independent service provider,
 * and therefore they expose the provider interface already. For this reason
 * both the service provider and the service interface need to have virtual
 * inheritance from the service provider interface, to make sure that there is
 * only one underlying interface object.
 * 
 */
class icarusDB::IICARUSChannelMap
  : virtual public icarusDB::IICARUSChannelMapProvider
{
    public:
  
  /// Type of the service provider.
  using provider_type = icarusDB::IICARUSChannelMapProvider;
  
  /// Returns the service provider.
  provider_type const* provider() const { return this; }
  
};


DECLARE_ART_SERVICE_INTERFACE(icarusDB::IICARUSChannelMap, SHARED)


#endif // ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAP_H
