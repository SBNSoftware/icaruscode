/**
 * @file   icaruscode/Decode/ChannelMapping/Legacy/ICARUSChannelMap.cc
 * @brief  Service registration for `icarusDB::ICARUSChannelMap` service.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#include "icaruscode/Decode/ChannelMapping/Legacy/ICARUSChannelMap.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap)


// -----------------------------------------------------------------------------
