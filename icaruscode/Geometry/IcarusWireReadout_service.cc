/**
 * @file    icaruscode/Geometry/IcarusWireReadout_service.cc
 * @brief   Geometry helper service for ICARUS geometries: implementation file.
 * @see     `icaruscode/Geometry/IcarusWireReadout.h`
 */

// library header
#include "icaruscode/Geometry/IcarusWireReadout.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadoutSetupTool.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

// framework libraries
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>



//------------------------------------------------------------------------------
namespace {

  // name of the tool to create the default channel mapping algorithm
  std::string const DefaultChannelMapSetupTool
    = "ICARUSsplitInductionChannelMapSetupTool";

  std::unique_ptr<geo::WireReadoutGeom>
  makeChannelMapping(fhicl::ParameterSet const& parameters, geo::GeometryCore const* geom)
  {
    fhicl::ParameterSet mapperDefaultSet;
    mapperDefaultSet.put("tool_type", DefaultChannelMapSetupTool);
    auto wireReadoutSetupTool = art::make_tool<geo::WireReadoutSetupTool>
      (parameters.get<fhicl::ParameterSet>("Mapper", mapperDefaultSet));
    auto wireReadoutSorter = art::make_tool<geo::WireReadoutSorter>
      (parameters.get<fhicl::ParameterSet>("SortingParameters"));

    return wireReadoutSetupTool->setupWireReadout(geom, std::move(wireReadoutSorter));
  }
} // local namespace


//------------------------------------------------------------------------------
icarus::IcarusWireReadout::IcarusWireReadout(fhicl::ParameterSet const& pset)
{
  // detector type check
  art::ServiceHandle<geo::Geometry> geom;
  if (geom->DetectorName().find("icarus") == std::string::npos) {
    MF_LOG_WARNING("IcarusWireReadout")
      << "Using a ICARUS channel mapping with an unsupported (non-ICARUS?)"
         " detector geometry";
  } // if not ICARUS detector

  fWireReadoutGeom = makeChannelMapping(pset, geom.get());
}


//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL
  (icarus::IcarusWireReadout, geo::WireReadout)


//------------------------------------------------------------------------------
