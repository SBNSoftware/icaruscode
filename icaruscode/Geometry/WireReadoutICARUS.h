/**
 * @file    icaruscode/Geometry/WireReadoutICARUS.h
 * @brief   Geometry helper service for ICARUS geometries.
 * @see     `icaruscode/Geometry/WireReadoutICARUS_service.cc`
 */

#ifndef ICARUSCODE_GEOMETRY_WIREREADOUTICARUS_H
#define ICARUSCODE_GEOMETRY_WIREREADOUTICARUS_H

// ICARUS libraries
#include "icarusalg/Geometry/WireReadoutGeomICARUS.h"

// LArSoft libraries
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/fwd.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Table.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr


// -----------------------------------------------------------------------------
namespace icarus { class WireReadoutICARUS; }
/**
 * @brief Implementation of `geo::WireReadout` for ICARUS.
 *
 * This _art_ service sets up and returns the ICARUS standard wire readout
 * service provider, `icarus::WireReadoutGeomICARUS`.
 *
 * Configuration
 * --------------
 *
 * * *SortingParameters* (tool parameter set): configuration of the
 *   `geo::WireReadoutSorter` _art_ tool to sort the geometry components under
 *   `WireReadoutGeom` (wire planes and their wires).
 * * *Mapper* (`icarus::WireReadoutGeomICARUS` parameter set): full
 *   configuration of `icarus::WireReadoutGeomICARUS`.
 *
 */
class icarus::WireReadoutICARUS: public geo::WireReadout {

    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::DelegatedParameter SortingParameters {
      Name{ "SortingParameters" },
      Comment{ "configuration of the art tool to sort the wire planes and wires" }
      };
    
    fhicl::Table<icarus::WireReadoutGeomICARUS::Config> Mapper {
      Name{ "Mapper" },
      Comment{ "configuration of the ICARUS channel mapping" }
      };
    
  };
  
  using Parameters = art::ServiceTable<Config>;

  /// Constructor: records the configuration.
  WireReadoutICARUS(Parameters const& params);
  
  
  // --- BEGIN  -- Service provider protocol  ----------------------------------
  
  /// Type of the service provider returned.
  using provider_type = icarus::WireReadoutGeomICARUS;
  
  /// Returns the owned service provider.
  provider_type const* provider() const { return fWireReadoutGeom.get(); }
  
  // --- END    -- Service provider protocol  ----------------------------------
  
  
    private:

  /// The service provider (initialized on construction and always ready).
  std::unique_ptr<provider_type> fWireReadoutGeom;

  // --- BEGIN -- Virtual interface definitions --------------------------------
  geo::WireReadoutGeom const& wireReadoutGeom() const override
    { return *fWireReadoutGeom; }
  // --- END -- Virtual interface definitions ----------------------------------
  
  
  /// Creates and initializes the wire readout object.
  static std::unique_ptr<provider_type> makeWireReadout
    (Config const& config);
  
  /// Checks the detector the geometry belongs to. Prints warning on failure.
  static void geometryCheck();

  
}; // icarus::WireReadoutICARUS


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE_INTERFACE_IMPL(icarus::WireReadoutICARUS,
                                   geo::WireReadout,
                                   SHARED)

#endif // ICARUSCODE_GEOMETRY_WIREREADOUTICARUS_H
