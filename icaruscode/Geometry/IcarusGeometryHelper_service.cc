/**
 * @file    icaruscode/Geometry/IcarusGeometryHelper_service.cc
 * @brief   Geometry helper service for ICARUS geometries: implementation file.
 * @see     `icaruscode/Geometry/IcarusGeometryHelper.h`
 */

// library header
#include "icaruscode/Geometry/IcarusGeometryHelper.h"

// ICARUS libraries
// #include "icaruscode/Geometry/ICARUSChannelMapAlg.h"
// #include "icaruscode/Geometry/ChannelMapIcarusAlg.h"

// LArSoft libraries
//#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryCore.h" // geo::GeometryData_t

// framework libraries
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>


//------------------------------------------------------------------------------
namespace {
  
  // name of the tool to create the default channel mapping algorithm
  std::string const DefaultChannelMapSetupTool
    = "ICARUSsingleInductionChannelMapSetupTool";
  
} // local namespace


//------------------------------------------------------------------------------
void icarus::IcarusGeometryHelper::doConfigureChannelMapAlg
  (fhicl::ParameterSet const& /* sortingParameters */, geo::GeometryCore* geom)
{
  fChannelMap.reset();

  //
  // detector type check
  //
  std::string const detectorName = geom->DetectorName();
  if (detectorName.find("icarus") == std::string::npos) {
    MF_LOG_WARNING("IcarusGeometryHelper")
      << "Using a ICARUS channel mapping with an unsupported (non-ICARUS?)"
         " detector geometry";
  } // if not ICARUS detector
  
  //
  // channel mapping creation and setup
  //
  fChannelMap = makeChannelMapping(fPset);
  
  //
  // apply channel mapping to geometry
  //
  if (fChannelMap)
    geom->ApplyChannelMap(fChannelMap); // calls Initialize(fGeoData) for us
  
} // icarus::IcarusGeometryHelper::doConfigureChannelMapAlg()


//------------------------------------------------------------------------------
std::shared_ptr<const geo::ChannelMapAlg>
icarus::IcarusGeometryHelper::doGetChannelMapAlg() const
{
  return fChannelMap;
}


//------------------------------------------------------------------------------
std::unique_ptr<geo::ChannelMapAlg>
icarus::IcarusGeometryHelper::makeChannelMapping
  (fhicl::ParameterSet const& parameters) const
{
  fhicl::ParameterSet mapperDefaultSet;
  mapperDefaultSet.put("tool_type", DefaultChannelMapSetupTool);
  auto channelMapSetupTool = art::make_tool<geo::ChannelMapSetupTool>
    (parameters.get<fhicl::ParameterSet>("Mapper", mapperDefaultSet));
  
  return channelMapSetupTool->setupChannelMap();
} // icarus::IcarusGeometryHelper::makeChannelMapping()


//------------------------------------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL
  (icarus::IcarusGeometryHelper, geo::ExptGeoHelperInterface)


//------------------------------------------------------------------------------
