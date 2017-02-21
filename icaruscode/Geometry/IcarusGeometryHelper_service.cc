////////////////////////////////////////////////////////////////////////////////
/// \file IcarusGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> Icarus/Geometry
#include "Icarus/Geometry/IcarusGeometryHelper.h"
#include "Icarus/Geometry/ChannelMapIcarusAlg.h"

//#include "larcore/Geometry/ChannelMapAlg.h"
#include "larcore/Geometry/GeometryCore.h" // larcore. geo::GeometryData_t

#include "TString.h"


namespace Icarus
{

IcarusGeometryHelper::IcarusGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & reg )
:  fPset( pset )
   //fReg( reg )
{}

IcarusGeometryHelper::~IcarusGeometryHelper() throw()
{}  

void IcarusGeometryHelper::doConfigureChannelMapAlg( fhicl::ParameterSet const & sortingParameters, geo::GeometryCore* geom ) 
{
    fChannelMap.reset();
    std::string const detectorName = geom->DetectorName();

    if ( detectorName.find("microboone") == std::string::npos ) {
        std::cout << __PRETTY_FUNCTION__ << ": WARNING USING CHANNEL MAP ALG WITH NON-MICROBOONE GEO!" << std::endl;
    }

//    fChannelMap = std::make_shared<geo::ChannelMapIcarusAlg>( fPset, sortingParameters );
    fChannelMap = std::make_shared<geo::ChannelMapIcarusAlg>( sortingParameters );

    if ( fChannelMap )
    {
        geom->ApplyChannelMap(fChannelMap); // calls Initialize(fGeoData) for us
    }

}

std::shared_ptr<const geo::ChannelMapAlg> IcarusGeometryHelper::doGetChannelMapAlg() const
{
    return fChannelMap;
}

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(Icarus::IcarusGeometryHelper, geo::ExptGeoHelperInterface)
