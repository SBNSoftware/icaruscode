/**
 * @file   icaruscode/Analysis/trigger/details/CathodeCrossingUtils.cxx
 * @brief  Algorithms dealing with a trajectory and the cathode (implementation)
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 27, 2022
 * @see    icaruscode/Analysis/trigger/details/CathodeCrossingUtils.h
 */

// library header
#include "icaruscode/Analysis/trigger/details/CathodeCrossingUtils.h"

// LArSoft libraries
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator...


// -----------------------------------------------------------------------------
icarus::CathodeDesc_t icarus::findTPCcathode
(geo::Point_t const& point, geo::GeometryCore const& geom)
{
  geo::CryostatGeo const* pCryo = geom.PositionToCryostatPtr(point);
  if (!pCryo) return {}; // all invalid, goodbye
  
  return {
    findCathodeCenter(*pCryo)                   // center
      , -(pCryo->TPC(0).DriftDir())  // normal
      };
  
} // icarus::findTPCcathode()


// -----------------------------------------------------------------------------
double icarus::distance(geo::Point_t const& point, CathodeDesc_t const& cathode)
{ return (point - cathode.center).Dot(cathode.normal); }


// -----------------------------------------------------------------------------
geo::Point_t icarus::findCathodeCenter(geo::CryostatGeo const& cryo) {
  
  // cathode position: assumes one cathode plane shared by all TPC
  geo::vect::MiddlePointAccumulator cathodePos;
  
  for (geo::TPCGeo const& TPC: cryo.IterateTPCs())
    cathodePos.add(TPC.GetCathodeCenter());
  
  return cathodePos.middlePoint();
  
} // icarus::findCathodeCenter()
