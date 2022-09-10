/**
 * @file   icaruscode/Analysis/trigger/details/CathodeCrossingUtils.h
 * @brief  Algorithms dealing with a trajectory and the cathode.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 27, 2022
 * @see    icaruscode/Analysis/trigger/details/CathodeCrossingUtils.cxx
 */

#ifndef ICARUSCODE_ANALYSIS_TRIGGER_DETAILS_CATHODECROSSINGUTILS_H
#define ICARUSCODE_ANALYSIS_TRIGGER_DETAILS_CATHODECROSSINGUTILS_H


// ICARUS libraries
#include "icaruscode/Utilities/TrajectoryUtils.h" // util::findCrossingSegment()

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t ...

// C/C++ standard libraries
#include <vector>
#include <iterator> // std::distance()
#include <type_traits> // std::is_same_v
#include <cmath> // std::abs()


// -----------------------------------------------------------------------------
// forward declarations
namespace geo {
  class GeometryCore;
  class CryostatGeo;
}


// -----------------------------------------------------------------------------
namespace icarus {
  
  // --- BEGIN -- Data structures  ---------------------------------------------
  /// Simple description for the cathode.
  struct CathodeDesc_t {
    geo::Point_t center;
    geo::Vector_t normal;
  };
  
  
  /// Information about the cathode crossing of a path.
  struct CathodeCrossing_t {
    static std::size_t const NoIndex = std::numeric_limits<std::size_t>::max();
    std::size_t indexBefore = NoIndex; ///< Index of the point "before" cathode.
    std::size_t indexAfter = NoIndex; ///< Index of the point "after" cathode.
    double before = 0.0; ///< Length of the path "before" the cathode.
    double after = 0.0; ///< Length of the path "after" the cathode.
    geo::Point_t crossingPoint; ///< Trajectory crossing point.
    
    /// Returns whether the crossing information is valid.
    operator bool() const { return indexBefore != indexAfter; }
    
    /// Returns whether the crossing information is invalid.
    bool operator!() const { return indexBefore == indexAfter; }
    
  }; // CathodeCrossing_t
  
  
  // --- END ---- Data structures  ---------------------------------------------
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns cathode information for cryostat at the specified `point`.
   * 
   * The normal to the cathode is always the same as the normal of the first TPC
   * in the cryostat containing `point`. This ensures consistency within the
   * cryostat.
   */
  CathodeDesc_t findTPCcathode
    (geo::Point_t const& point, geo::GeometryCore const& geom);
  
  
  // ---------------------------------------------------------------------------
  /// Returns the distance of a `point` from the `cathode`.
  double distance(geo::Point_t const& point, CathodeDesc_t const& cathode);
  
  
  // ---------------------------------------------------------------------------
  /// Returns the center of the cathode in the specified cryostat.
  geo::Point_t findCathodeCenter(geo::CryostatGeo const& cryo);
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns the crossing point of a trajectory on the cathode.
   * @tparam Iter type of iterator to a point in space
   * @param begin iterator to the first point of the trajectory
   * @param end iterator past the last point of the trajectory
   * @param cathode description of the cathode plane
   * @return the crossing information
   * 
   * This algorithm intersects a path on a infinite plane in 3D space, and
   * returns the intersection point.
   * 
   * The path is a sequence of segments. For simple trajectories, there is only
   * one crossing, and then the returned point lies on both the path and the
   * plane (it is, indeed, the path/cathode intersection point).
   * For trajectories which cross the plane multiple times, the details of the
   * part of the trajectory between the first and the last crossing are ignored
   * and the two chunks of end paths that lie entirely on one side of the
   * cathode are joint by a segment, shortcutting the details of the path in
   * between.
   * 
   * The function also returns the length of the paths at the two sides of the
   * cathode. Again, in case of multiple crossings, the details of the crossing
   * are ignored and a single segment is used to connect the two sides, and the
   * two lengths reflect that approximation.
   * If no intersection point is found, the two partial lengths are both `0`.
   */
  template <typename Iter>
  CathodeCrossing_t detectCrossing
    (Iter begin, Iter end, CathodeDesc_t const& cathode);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus


// -----------------------------------------------------------------------------
// --- Template implementation
// -----------------------------------------------------------------------------
template <typename Iter>
icarus::CathodeCrossing_t icarus::detectCrossing
  (Iter begin, Iter end, CathodeDesc_t const& cathode)
{
  using Point_t = typename Iter::value_type;
  
  // this function by itself is already generic;
  // but the utilities and data it relies on (including the return value) aren't
  static_assert(std::is_same_v<Point_t, geo::Point_t>,
    "icarus::detectCrossing() only supports geo::Point_t points;"
    " if needed, ask the author for extension"
    );
  
  auto const isNegative = [&cathode]
    (auto const& p) { return distance(p, cathode) < 0.0; };
  
  auto const [ itBeforeCathode, itAfterCathode ]
    = util::findCrossingSegment(begin, end, isNegative);
  
  if (itBeforeCathode == itAfterCathode) {
    std::size_t const nPoints = std::distance(begin, end);
    return {
        nPoints // indexBefore
      , nPoints // indexAfter
        // the rest is default
      }; // no crossing at all
  }
  
  CathodeCrossing_t crossInfo;
  
  crossInfo.indexBefore = std::distance(begin, itBeforeCathode);
  crossInfo.indexAfter = std::distance(begin, itAfterCathode);
  
  // first estimation:
  crossInfo.before = util::pathLength(begin, itBeforeCathode);
  crossInfo.after = util::pathLength(itAfterCathode, end);
  
  auto const step = (*itAfterCathode - *itBeforeCathode);
  double const stepLength = geo::vect::norm(step);
    
  // the two ends share the step in proportion to their cathode distance
  double fBefore = std::abs(distance(*itBeforeCathode, cathode));
  double fAfter = std::abs(distance(*itAfterCathode, cathode));
  double const fTotal = fBefore + fAfter;
  fBefore /= fTotal;
  fAfter /= fTotal;
  
  crossInfo.before += fBefore * stepLength;
  crossInfo.after += fAfter * stepLength;
  
  crossInfo.crossingPoint = *itBeforeCathode + fBefore * step;
  
  return crossInfo;
  
} // icarus::detectCathodeCrossing()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_ANALYSIS_TRIGGER_DETAILS_CATHODECROSSINGUTILS_H
