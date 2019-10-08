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
#include "larcore/Geometry/ChannelMapSetupTool.h"
#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcorealg/Geometry/GeometryCore.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>, std::shared_ptr<>


// -----------------------------------------------------------------------------
// Forward declarations
namespace geo {
  class ChannelMapAlg;
}
namespace icarus {
  class IcarusGeometryHelper;
} // namespace icarus

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
  IcarusGeometryHelper
    (fhicl::ParameterSet const& pset /* , art::ActivityRegistry& reg */)
    : fPset(pset)
    {}

  /*
   * Public interface for ExptGeoHelperInterface (for reference purposes)
   * 
   * Configure and initialize the channel map.
   * 
   * void ConfigureChannelMapAlg(
   *   std::string const& detectorName,
   *   fhicl::ParameterSet const& sortingParam,
   *   std::vector<geo::CryostatGeo*>& c,
   *   std::vector<geo::AuxDetGeo*>& ad
   *   );
   * 
   * 
   * Returns null pointer if the initialization failed
   * NOTE:  the sub-class owns the ChannelMapAlg object
   * 
   * std::shared_ptr<const geo::ChannelMapAlg> & GetChannelMapAlg() const;
   * 
   */

    private:

  fhicl::ParameterSet fPset; ///< Copy of configuration parameter set.
  
  // This is shared because of how the channel mapping was originally used in
  // `geo::GeometryCore`. In the future this service will be entirely obsolete
  // and this weirdness will disappear. Maybe.
  /// The channel mapping algorithm.
  std::shared_ptr<geo::ChannelMapAlg> fChannelMap;
  
  
  // --- BEGIN -- Virtual interface definitions --------------------------------
  virtual void doConfigureChannelMapAlg
    (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
    override;
  
  virtual ChannelMapAlgPtr_t doGetChannelMapAlg() const override;
  // --- END -- Virtual interface definitions ----------------------------------
  
  
  /// Creates and returns the channel mapping instance via a _art_ tool.
  std::unique_ptr<geo::ChannelMapAlg> makeChannelMapping
    (fhicl::ParameterSet const& parameters) const;
  
  
}; // icarus::IcarusGeometryHelper


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL
  (icarus::IcarusGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_GEOMETRY_ICARUSGEOMETRYHELPER_H
