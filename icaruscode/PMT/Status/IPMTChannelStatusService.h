/**
 * @file   icaruscode/PMT/Status/IPMTChannelStatusService.h
 * @brief  Art service interface for PMT channel status.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 */

#ifndef ICARUSCODE_PMT_STATUS_IPMTCHANNELSTATUSSERVICE_H
#define ICARUSCODE_PMT_STATUS_IPMTCHANNELSTATUSSERVICE_H

// ICARUS libraries
#include "icaruscode/PMT/Status/PMTChannelStatus.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceProviderWrappers.h"


// -----------------------------------------------------------------------------
namespace icarusDB {
  /// The only thing this service does is to return its service provider of type
  /// `icarusDB::PMTChannelStatus`.
  using IPMTChannelStatusService
    = lar::ServiceProviderInterfaceWrapper<PMTChannelStatus>;
}


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE(icarusDB::IPMTChannelStatusService, SHARED)


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_STATUS_IPMTCHANNELSTATUSSERVICE_H
