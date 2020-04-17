/**
 * @file icaruscode/Geometry/ICARUSstandaloneGeometrySetup.h
 * @brief  Functions to facilitate ICARUS geometry initialization outside _art_.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 16, 2020
 * 
 * This is a header-only library.
 * It does include some static code (inlined).
 * 
 * If this library becomes too heavy, it can be left with only the template
 * definition `icarus::geo::ICARUSStandaloneGeometrySetup()`,
 * while each specialization of `lar::standalone::SetupGeometry()` can end up
 * into its own header file.
 */

#ifndef ICARUSCODE_GEOMETRY_ICARUSSTANDALONEGEOMETRYSETUP_H
#define ICARUSCODE_GEOMETRY_ICARUSSTANDALONEGEOMETRYSETUP_H


// ICARUS libraries
#include "icaruscode/Geometry/ChannelMapIcarusAlg.h"
#include "icaruscode/Geometry/ICARUSChannelMapAlg.h"

// LArSoft libraries
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"

// framework libraries
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <set>
#include <memory> // std::make_unique()
#include <utility> // std::forward()
#include <type_traits> // std::enable_if_t, std::void_t


namespace icarus::geo {
  
  namespace details {
    
    /// Creates a configuration object for `Class` from a parameter set.
    template <typename Class, typename = void>
    struct ConfigObjectMaker;
    
  } // namespace details
  
  
  
  /**
   * @brief Initialization of geometry with ICARUS specific conventions.
   * @tparam ChannelMapClass type of channel mapping to be used
   * @param pset configuration of the geometry service
   * @return an initialized instance of `geo::GeometryCore` service provider
   * 
   * The conventions that are specific to ICARUS are:
   * 
   * * a copy of the channel mapping algorithm configuration is kept in the
   *   `Geometry` service configuration table, as `ChannelMapping` table;
   *   because of the unusual requirement, this item is _mandatory_ (even if
   *   it may well be empty);
   * * channel mapping objects may take a `ChannelMapClass::Config` object as
   *   configuration, which can be wrapped in a FHiCL table;
   * * channel mapping configuration _may_ contain a spurious `tool_type` entry
   *   which is ignored.
   * 
   */
  template <typename ChannelMapClass, typename... Args>
  std::unique_ptr<::geo::GeometryCore>
  SetupICARUSGeometry(fhicl::ParameterSet const& pset, Args&&... args) {
    
    auto const& config = details::ConfigObjectMaker<ChannelMapClass>::make
      (pset.get<fhicl::ParameterSet>("ChannelMapping"));
    
    auto channelMap = std::make_unique<ChannelMapClass>
      (config, std::forward<Args>(args)...);
    
    return lar::standalone::SetupGeometryWithChannelMapping
      (pset, move(channelMap));
    
  } // SetupICARUSGeometry()
  
} // namespace icarus::geo


// specializations
namespace icarus::geo::details {

  /// General implementation: passes the parameter set through
  template <typename Class, typename /* = void */>
  struct ConfigObjectMaker {
    
    static fhicl::ParameterSet const& make(fhicl::ParameterSet const& pset)
      { return pset; }
    
  }; // struct ConfigObjectMaker


  /// Specialization: class with a `Config` configuration data structure.
  template <typename Class>
  struct ConfigObjectMaker<Class, std::void_t<typename Class::Config>> {
    
    using Config = typename Class::Config;
    using Parameters = fhicl::Table<Config>;
    
    static Config make(fhicl::ParameterSet const& pset)
      { return Parameters{ pset, std::set<std::string>{ "tool_type" } }(); }
    
  }; // struct ConfigObjectMaker
  
} // namespace icarus::geo::details


namespace lar::standalone {
  
  // ---------------------------------------------------------------------------
  /// Specialization of `lar::standalone::SetupGeometry()`
  /// for ICARUS channel mapping `geo::ChannelMapIcarusAlg`.
  template <>
  inline std::unique_ptr<geo::GeometryCore>
  SetupGeometry<geo::ChannelMapIcarusAlg>(fhicl::ParameterSet const& pset)
    { return icarus::geo::SetupICARUSGeometry<geo::ChannelMapIcarusAlg>(pset); }
  
  
  // ---------------------------------------------------------------------------
  /// Specialization of `lar::standalone::SetupGeometry()`
  /// for ICARUS channel mapping `icarus::ICARUSChannelMapAlg`.
  template <>
  inline std::unique_ptr<geo::GeometryCore>
  SetupGeometry<icarus::ICARUSChannelMapAlg>(fhicl::ParameterSet const& pset)
    { return icarus::geo::SetupICARUSGeometry<icarus::ICARUSChannelMapAlg>(pset); }
  
  
  // ---------------------------------------------------------------------------
  
  
} // namespace lar::standalone



#endif // ICARUSCODE_GEOMETRY_ICARUSSTANDALONEGEOMETRYSETUP_H
