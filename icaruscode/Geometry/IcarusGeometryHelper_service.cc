/**
 * @file    icaruscode/Geometry/IcarusGeometryHelper_service.cc
 * @brief   Geometry helper service for ICARUS geometries: implementation file.
 * @see     `icaruscode/Geometry/IcarusGeometryHelper.h`
 */

// library header
#include "icaruscode/Geometry/IcarusGeometryHelper.h"

// LArSoft libraries
#include "larcore/Geometry/ChannelMapSetupTool.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"

// framework libraries
#include "art/Utilities/make_tool.h"
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
auto icarus::IcarusGeometryHelper::doConfigureChannelMapAlg(
  fhicl::ParameterSet const& /* sortingParameters */,
  std::string const& detectorName
  ) const -> ChannelMapAlgPtr_t
{
  //
  // detector type check
  //
  if (detectorName.find("icarus") == std::string::npos) {
    MF_LOG_WARNING("IcarusGeometryHelper")
      << "Using a ICARUS channel mapping with an unsupported (non-ICARUS?)"
         " detector geometry";
  } // if not ICARUS detector

  //
  // channel mapping creation and setup
  //
  return makeChannelMapping(fPset);

} // icarus::IcarusGeometryHelper::doConfigureChannelMapAlg()


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
