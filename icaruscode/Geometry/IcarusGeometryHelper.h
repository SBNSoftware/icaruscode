/**
 * @file    icaruscode/Geometry/IcarusGeometryHelper.h
 * @brief   Geometry helper service for ICARUS geometries.
 * @see     `icaruscode/Geometry/IcarusGeometryHelper_service.cc`
 *
 * Handles Icarus-specific information for the generic Geometry service
 * within LArSoft. Derived from the `geo::ExptGeoHelperInterface` class.
 */

#ifndef ICARUSCODE_GEOMETRY_ICARUSGEOMETRYHELPER_H
#define ICARUSCODE_GEOMETRY_ICARUSGEOMETRYHELPER_H

// LArSoft libraries
// #include "larcore/Geometry/ChannelMapSetupTool.h"
#include "larcore/Geometry/ExptGeoHelperInterface.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>, std::shared_ptr<>


// -----------------------------------------------------------------------------
// Forward declarations
namespace geo { class ChannelMapAlg; }
namespace icarus { class IcarusGeometryHelper; }

// -----------------------------------------------------------------------------
/**
 * @brief Implementation of `geo::ExptGeoHelperInterface` for ICARUS.
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
class icarus::IcarusGeometryHelper: public geo::ExptGeoHelperInterface {

    public:

  /// Constructor: records the configuration.
  IcarusGeometryHelper(fhicl::ParameterSet const& pset): fPset(pset) {}

    private:

  fhicl::ParameterSet fPset; ///< Copy of configuration parameter set.

  // --- BEGIN -- Virtual interface definitions --------------------------------
  virtual ChannelMapAlgPtr_t doConfigureChannelMapAlg(
    fhicl::ParameterSet const& /* sortingParameters */,
    std::string const& detectorName
    ) const override;

  // --- END -- Virtual interface definitions ----------------------------------

  /// Creates and returns the channel mapping instance via a _art_ tool.
  std::unique_ptr<geo::ChannelMapAlg> makeChannelMapping
    (fhicl::ParameterSet const& parameters) const;

}; // icarus::IcarusGeometryHelper


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarus::IcarusGeometryHelper,
                                   geo::ExptGeoHelperInterface,
                                   SHARED)

#endif // ICARUSCODE_GEOMETRY_ICARUSGEOMETRYHELPER_H
