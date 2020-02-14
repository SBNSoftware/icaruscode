////////////////////////////////////////////////////////////////////////////////
/// \file IcarusGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> Icarus/Geometry
#include "icaruscode/Geometry/IcarusGeometryHelper.h"
#include "icaruscode/Geometry/ChannelMapIcarusAlg.h"

#include "larcorealg/Geometry/GeometryCore.h" // larcore. geo::GeometryData_t

namespace Icarus
{

IcarusGeometryHelper::IcarusGeometryHelper(fhicl::ParameterSet const & pset)
  : fPset(pset)
{}

std::unique_ptr<geo::ChannelMapAlg>
IcarusGeometryHelper::doConfigureChannelMapAlg(fhicl::ParameterSet const & sortingParameters,
                                               std::string const& detectorName) const
{
    if ( detectorName.find("icarus") == std::string::npos ) {
        std::cout << __PRETTY_FUNCTION__ << ": WARNING USING CHANNEL MAP ALG WITH NON-ICARUS GEO!" << std::endl;
    }

    return std::make_unique<geo::ChannelMapIcarusAlg>(fPset);
}

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(Icarus::IcarusGeometryHelper, geo::ExptGeoHelperInterface)
