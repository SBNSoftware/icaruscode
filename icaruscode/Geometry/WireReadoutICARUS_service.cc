/**
 * @file    icaruscode/Geometry/WireReadoutICARUS_service.cc
 * @brief   Geometry helper service for ICARUS geometries: implementation file.
 * @see     `icaruscode/Geometry/WireReadoutICARUS.h`
 */

// library header
#include "icaruscode/Geometry/WireReadoutICARUS.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadoutSetupTool.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

// framework libraries
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h" // MF_LOG_WARNING()

// C/C++ standard libraries
#include <string>


//------------------------------------------------------------------------------
icarus::WireReadoutICARUS::WireReadoutICARUS(Parameters const& params)
  : fWireReadoutGeom{ makeWireReadout(params()) }
{
  geometryCheck();
}


//------------------------------------------------------------------------------
auto icarus::WireReadoutICARUS::makeWireReadout
  (Config const& config) -> std::unique_ptr<provider_type>
{
  std::unique_ptr wireReadoutSorter = art::make_tool<geo::WireReadoutSorter>
    (config.SortingParameters.get<fhicl::ParameterSet>());
  
  return std::make_unique<provider_type>(
    config.Mapper(),
    lar::providerFrom<geo::Geometry>(), std::move(wireReadoutSorter)
    );
} // icarus::WireReadoutICARUS::makeWireReadout()


//------------------------------------------------------------------------------
void icarus::WireReadoutICARUS::geometryCheck() {
  
  // detector type check
  geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
  if (geom->DetectorName().find("icarus") == std::string::npos) {
    MF_LOG_WARNING("WireReadoutICARUS")
      << "Using a ICARUS channel mapping with an unsupported (non-ICARUS?)"
         " detector geometry '" << geom->DetectorName() << "'";
  } // if not ICARUS detector
  
} // icarus::WireReadoutICARUS::geometryCheck()


//------------------------------------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL
  (icarus::WireReadoutICARUS, geo::WireReadout)


//------------------------------------------------------------------------------
