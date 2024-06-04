/**
 * @file    icaruscode/Geometry/IcarusGeometryHelper.h
 * @brief   Geometry helper service for ICARUS geometries.
 * @see     `icaruscode/Geometry/IcarusGeometryHelper_service.cc`
 *
 * Handles Icarus-specific information for the generic Geometry service
 * within LArSoft. Derived from the `geo::WireReadout` class.
 */

#ifndef ICARUSCODE_GEOMETRY_ICARUSWIREREADOUT_H
#define ICARUSCODE_GEOMETRY_ICARUSWIREREADOUT_H

// LArSoft libraries
#include "larcorealg/Geometry/fwd.h"
#include "larcore/Geometry/WireReadout.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>, std::shared_ptr<>


// -----------------------------------------------------------------------------
// Forward declarations
namespace icarus { class IcarusWireReadout; }

// -----------------------------------------------------------------------------
/**
 * @brief Implementation of `geo::WireReadout` for ICARUS.
 *
 * This service utilizes a _art_ tool to create the proper channel mapper
 * algorithm instance.
 *
 *
 * Configuration
 * --------------
 *
 * * *Mapper* (tool parameter set):
 *     the configuration includes a `tool_type` parameter to identify the tool
 *     to use, and the rest of the configuration is passed to the tool itself.
 *     The standard tool (_not default!_) is `icarus::ICARUSChannelMapAlgTool`
 *     By default, uses `icarus::ICARUSSingleInductionChannelMapAlgTool`
 *     (legacy).
 *
 *
 */
class icarus::IcarusWireReadout: public geo::WireReadout {

    public:

  /// Constructor: records the configuration.
  IcarusWireReadout(fhicl::ParameterSet const& pset);

    private:

  std::unique_ptr<geo::WireReadoutGeom> fWireReadoutGeom;

  // --- BEGIN -- Virtual interface definitions --------------------------------
  geo::WireReadoutGeom const& wireReadoutGeom() const override { return *fWireReadoutGeom; }
  // --- END -- Virtual interface definitions ----------------------------------

}; // icarus::IcarusWireReadout


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarus::IcarusWireReadout,
                                   geo::WireReadout,
                                   SHARED)

#endif // ICARUSCODE_GEOMETRY_ICARUSWIREREADOUT_H
