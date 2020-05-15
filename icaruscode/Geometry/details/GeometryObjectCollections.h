/**
 * @file   icaruscode/Geometry/details/GeometryObjectCollections.h
 * @brief  A few simple data type definitions.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 10, 2019
 */

#ifndef ICARUSCODE_GEOMETRY_DETAILS_GEOMETRYOBJECTCOLLECTIONS_H
#define ICARUSCODE_GEOMETRY_DETAILS_GEOMETRYOBJECTCOLLECTIONS_H


// LArSoft libraries
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// C/C++ standard library
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus::details {
  
  // forward declarations
  struct ChannelRange_t;
  
  //
  // data type definitions
  //
  
  /// Type of collection of TPCs (pointers to `geo::TPCGeo`).
  using TPCColl_t = std::vector<geo::TPCGeo const*>;
  
  /// Type of collection of planes (pointers to `geo::PlaneGeo`).
  using PlaneColl_t = std::vector<geo::PlaneGeo const*>;
  
} // namespace icarus::details


// -----------------------------------------------------------------------------
// --- icarus::details::ChannelRange_t
// -----------------------------------------------------------------------------
namespace icarus::details { struct ChannelRange_t; }


/// A simple range of channels.
struct icarus::details::ChannelRange_t
  : std::pair<raw::ChannelID_t, raw::ChannelID_t>
{
  
  using std::pair<raw::ChannelID_t, raw::ChannelID_t>::pair;
  
  /// Returns the ID of the first channel in the range.
  constexpr raw::ChannelID_t begin() const { return first; }
  
  /// Returns the ID of the channel after the last one in the range.
  constexpr raw::ChannelID_t end() const { return second; }
  
  /// Returns whether this range includes the specified `channel`.
  constexpr bool contains(raw::ChannelID_t channel) const
    { return (channel >= begin()) && (channel < end()); }
  
}; // struct icarus::details::ChannelRange_t


// ----------------------------------------------------------------------------

#endif // ICARUSCODE_GEOMETRY_DETAILS_GEOMETRYOBJECTCOLLECTIONS_H
