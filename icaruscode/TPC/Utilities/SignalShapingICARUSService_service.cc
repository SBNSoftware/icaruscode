////////////////////////////////////////////////////////////////////////
/// \file   icaruscode/TPC/Utilities/SignalShapingICARUSService_service.cc
/// \author H. Greenlee
/// Modified by X. Qian 1/6/2015
/// if histogram is used, inialize
/// Response_Offset, Response_Sampling, FieldBins from histogram
////////////////////////////////////////////////////////////////////////

#include "icaruscode/TPC/Utilities/SignalShapingICARUSService.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

DEFINE_ART_SERVICE(icarusutil::SignalShapingICARUSService)
