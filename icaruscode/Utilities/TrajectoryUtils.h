/**
 * @file   icaruscode/Utilities/TrajectoryUtils.h
 * @brief  Algorithms dealing with a trajectory as a sequence of 3D points.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 26, 2022
 * 
 * This is a header-only library.
 * 
 * @todo Move this file into `icarusalg/Utilities`.
 */

#ifndef ICARUSALG_UTILITIES_TRAJECTORYUTILS_H
#define ICARUSALG_UTILITIES_TRAJECTORYUTILS_H


// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::norm()

// C/C++ standard libraries
#include <vector>
#include <utility> // std::pair
#include <iterator> // std::next()


namespace util {
  
  // ---------------------------------------------------------------------------
  /// Data returned by `findMiddlePointInPath()`.
  template <typename Iter>
  struct PathPointInfo_t {
    double length { -1 }; ///< Path length; negative means not available.
    double step { -1 }; ///< Distance from the detected point to the next one.
    double frac { 0.5 }; ///< The point is `frac*step` away from the point.
    Iter itPoint; ///< Iterator to the path point just before the required one.
  }; // PathPointInfo_t
  
  /**
   * @brief Returns information to identify the middle of the specified path.
   * @tparam FIter type of iterator to a point in the path
   * @tparam LIter type of iterator to a point in the path
   * @param itFirst iterator to the first point in the path
   * @param itLast iterator to the last point in the path (must be valid point)
   * @param relTarget (default: `0.5`) fraction of the total length to pursue
   * @return information to identify the desired point
   * 
   * The sequence between the iterators `itFirst` and `itLast`, both included,
   * defines a path. This function returns the point that cuts that path into
   * two subpaths of the same length. The point usually does not match any of
   * the points in the sequence, but rather is in between two of them.
   * 
   * The path is considered a sequence of straight segments connecting the
   * points.
   * 
   * Example of use with a LArSoft `recob::Track` object (`track`):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const& path = track.Trajectory().Positions();
   * geo::Point_t const middle
   *   = util::PathMiddlePoint(path.begin(), std::prev(path.end()));
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (note that here the function actually returns a
   * `recob::tracking::Point_t`).
   * 
   * 
   * Requirements
   * -------------
   * 
   * The `FIter` type is a forward iterator, the type `LIter` is a backward
   * iterator.
   *  * `FIter` must support prefix `operator++`, `LIter` must support prefix
   *    `operator--()`.
   *  * `FIter` and `LIter` must compare equal (`operator == (FIter, LIter)`)
   *    if they point to the same element of the sequence.
   *  * The two iterators must point to the same type of point (called `Point`
   *    in the following text).
   * 
   * The type `Point` must satisfy the following operations.
   *  * `Vector operator- (Point, Point)` returns an object describing the
   *    displacement to go from the second point to the first one.
   *  * `geo::vect::norm(Vector)` is a function returning the magnitude of the
   *    specified displacement vector.
   * 
   * LArSoft `geo::Point_t` type is known to satisfy the requirements of
   * `Point` (`geo::Vector_t` being its corresponding `Vector` type). The
   * modulus function `geo::vect::norm()` is provided for `geo::Vector_t` in the
   * LArSoft header `larcorealg/CoreUtils/geo_vector_utils.h`.
   * 
   */
  template <typename FIter, typename LIter>
  PathPointInfo_t<FIter> findMiddlePointInPath
    (FIter itFirst, LIter itLast, double relTarget = 0.5);
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns the geometric point in the middle of the specified path.
   * @tparam FIter type of iterator to a point in the path
   * @tparam LIter type of iterator to a point in the path
   * @param itFirst iterator to the first point in the path
   * @param itLast iterator to the last point in the path (must be valid point)
   * @param relTarget (default: `0.5`) fraction of the total length to pursue
   * @return the geometric middle point of the path
   * 
   * The sequence between the iterators `itFirst` and `itLast`, both included,
   * defines a path. This function returns the point that cuts that path into
   * two subpaths, the first of which has `relTarget` fraction of the length of
   * the total path (hence, with the default `relTarget` value of `0.5` the two
   * subpaths have the same length). The point usually does not match any of
   * the points in the sequence, but rather is in between two of them.
   * 
   * The path is considered a sequence of straight segments connecting the
   * points.
   * 
   * Example of use with a LArSoft `recob::Track` object (`track`):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const& path = track.Trajectory().Positions();
   * geo::Point_t const middle
   *   = util::PathMiddlePoint(path.begin(), std::prev(path.end()));
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (note that here the function actually returns a
   * `recob::tracking::Point_t`).
   * 
   * 
   * Requirements
   * -------------
   * 
   * The `FIter` type is a forward iterator, the type `LIter` is a backward
   * iterator.
   *  * `FIter` must support prefix `operator++`, `LIter` must support prefix
   *    `operator--()`.
   *  * `FIter` and `LIter` must compare equal (`operator == (FIter, LIter)`)
   *    if they point to the same element of the sequence.
   *  * The two iterators must point to the same type of point (called `Point`
   *    in the following text).
   * 
   * The type `Point` must satisfy the following operations.
   *  * `Vector operator- (Point, Point)` returns an object describing the
   *    displacement to go from the second point to the first one.
   *  * `Point operator+ (Point, Vector)` returns the point after being
   *    displaced by the specified vector.
   *  * `Vector operator* (double)` returns a displacement vector with its
   *    magnitude rescaled by the specified factor.
   *  * `geo::vect::norm(Vector)` is a function returning the magnitude of the
   *    specified displacement vector.
   * 
   * The returned object is an instance of `Point`.
   * 
   * LArSoft `geo::Point_t` type is known to satisfy the requirements of
   * `Point` (`geo::Vector_t` being its corresponding `Vector` type). The
   * modulus function `geo::vect::norm()` is provided for `geo::Vector_t` in the
   * LArSoft header `larcorealg/CoreUtils/geo_vector_utils.h`.
   * 
   */
  template <typename FIter, typename LIter>
  auto pathMiddlePoint(FIter itFirst, LIter itLast, double relTarget = 0.5);
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns the total length of the specified path.
   * @tparam BIter type of iterator to the start of the trajectory
   * @tparam EIter type of iterator to the end of the trajectory
   * @param begin iterator to the first point of the trajectory
   * @param end iterator past the last point of the trajectory
   * @return the length of the path
   * @see `partialPathLengths()`
   * 
   * The sequence of partial lengths is returned, with one entry per point.
   * The length is the sum of the distances between adjacent points of the
   * sequence.
   * Empty sequences and sequences of a single points have length `0`.
   * 
   * The length are in the same units as the input point coordinates.
   * 
   * 
   * Requirements
   * -------------
   * 
   * The `BIter` type is a forward iterator.
   *  * `BIter` must support prefix `operator++`
   *  * `operator != (BIter, EIter)` must return whether the first iterator
   *    is equivalent to the second (or rather, whether the sequence is over).
   * 
   * The type `Point` must satisfy the following operations.
   *  * `Vector operator- (Point, Point)` returns an object describing the
   *    displacement to go from the second point to the first one.
   *  * `geo::vect::norm(Vector)` is a function returning the magnitude of the
   *    specified displacement vector.
   */
  template <typename BIter, typename EIter>
  double pathLength(BIter begin, EIter end);
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns a sequences of partial lengths for the specified path.
   * @tparam BIter type of iterator to the start of the trajectory
   * @tparam EIter type of iterator to the end of the trajectory
   * @param begin iterator to the first point of the trajectory
   * @param end iterator past the last point of the trajectory
   * @return a sequence of partial lengths, one per point
   * @see `pathLength()`
   * 
   * The sequence of partial lengths is returned, with one entry per point.
   * A partial length up to the point _i_ is the sum of the distances between
   * adjacent points starting from the first point (`*begin`) to the _i_-th,
   * included.
   * The first entry of the sequence is always `0`.
   * 
   * The length are in the same units as the input point coordinates.
   * 
   * 
   * Requirements
   * -------------
   * 
   * The `BIter` type is a forward iterator.
   *  * `BIter` must support prefix `operator++`
   *  * `operator != (BIter, EIter)` must return whether the first iterator
   *    is equivalent to the second (or rather, whether the sequence is over).
   * 
   * The type `Point` must satisfy the following operations.
   *  * `Vector operator- (Point, Point)` returns an object describing the
   *    displacement to go from the second point to the first one.
   *  * `geo::vect::norm(Vector)` is a function returning the magnitude of the
   *    specified displacement vector.
   */
  template <typename BIter, typename EIter>
  std::vector<double> partialPathLengths(BIter begin, EIter end);
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns a path segment with ends on different sides than path ends.
   * @tparam Iter type of iterator to the points in the path
   * @tparam SideFunc type of functor returning the side of a point
   * @param begin iterator to the first point of the path
   * @param end iterator past the last point of the path
   * @param sideOf functor returning the side of the point in argument
   * @return a pair of iterators pointing to the points found, `{ end, end }`
   *         if the path does not change side at all
   * 
   * The path is defined as a sequence of positions of arbitrary types.
   * The functor `sideOf` takes as argument one such position, and returns an
   * also arbitrary value representing the side of that position.
   * This algorithm returns a pair with the first point where the side changes,
   * and the point after the last side change.
   * The two points match only if there is no side change at all (in which case
   * a pair with two `end` iterators is returned). The two returned iterators
   * are adjacent only if there is only one change of side.
   * 
   * Note that this algorithm is effectively not specific to a space path but
   * it's rather a generic classification algorithm: given an ordered sequence,
   * it finds the first and the last change of class.
   * 
   * 
   * Requirements
   * -------------
   * 
   * * `Iter` must be a bidirectional operator (ask if this is too restrictive)
   * * `SideFunc` must support a call like.
   *   `Side SideFunc::operator() (Iter::value_type) const`, i.e. must accept
   *   as argument the position pointed by an iterator of type `Iter` and return
   *   an arbitrary `Side` value.
   * * `Side` type (returned by `sideOf`) must support
   *   `operator!= (Side, Side)`.
   * 
   */
  template <typename Iter, typename SideFunc>
  std::pair<Iter, Iter> findCrossingSegment
    (Iter begin, Iter end, SideFunc sideOf);
    
  // ---------------------------------------------------------------------------
  
} // namespace util


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename FIter, typename LIter>
util::PathPointInfo_t<FIter> util::findMiddlePointInPath
  (FIter itFirst, LIter itLast, double relTarget /* = 0.5 */)
{
  
  // We don't check that `LIter` also points to `Point_t`.
  // The return value will be `Point_t` anyway.
  
  // TODO if this algorithm works as of now, replace SegmentInfo_t with a simple double
  struct SegmentInfo_t {
    double length = 0.0;
    
    bool operator< (SegmentInfo_t const& other) const
      { return length < other.length; }
    
  };
  
  // single point: return that point
  if (itFirst == itLast) {
    return PathPointInfo_t<FIter>{
        0.0      // length
      , 0.0      // step
      , 0.0      // frac
      , itFirst  //itPoint
      };
  }
  
  SegmentInfo_t startChunk, endChunk;
  FIter itStart = itFirst, itPrevStart = itStart;
  LIter itEnd = itLast, itPrevEnd = itEnd;
  
  do {
    // advance the shortest of the chunks
    if (startChunk < endChunk) { // advance start
      itPrevStart = itStart;
      if (++itStart == itEnd) {
        // note that in this case start chunk is already shorter so itPrevStart
        // will not be needed later (at most, the last end step will be undone)
        itStart = itPrevStart; // we went too far: restore to the previous point
        break;
      }
      startChunk.length += geo::vect::norm(*itStart - *itPrevStart);
    }
    else { // advance end
      itPrevEnd = itEnd;
      if (itStart == --itEnd) {
        // note that in this case end chunk is already shorter so itPrevEnd
        // will not be needed later (at most the last start step will be undone)
        itEnd = itPrevEnd; // we went too far: restore to the previous point
        break;
      }
      endChunk.length += geo::vect::norm(*itPrevEnd - *itEnd);
    }
    assert(!(itStart == itEnd));
  } while(true);
  
  // now we have the two points at the sides of a segment,
  // and the amount of path behind each;
  // it may still be that the last step of one of the two ends jumped over the
  // middle point; in that case we move back
  
  assert(!(itStart == itEnd));
  
  auto delta = *itEnd - *itStart;
  double deltaLength = geo::vect::norm(delta);
  
  double const pathLength = startChunk.length + deltaLength + endChunk.length;
  
  double const targetLength = pathLength * relTarget;
  
  if (startChunk.length > targetLength) {
    // give the last step back to end...
    --itEnd;
    endChunk.length += deltaLength;
    assert(itStart == itEnd);
    // ... and take it away from start
    itStart = itPrevStart;
    delta = *itEnd - *itStart;
    deltaLength = geo::vect::norm(delta);
    startChunk.length -= deltaLength;
  }
  else if (endChunk.length > targetLength) {
    // give one step back to start...
    ++itStart;
    startChunk.length += deltaLength;
    assert(itStart == itEnd);
    // ... and take it away from end
    itEnd = itPrevEnd;
    delta = *itEnd - *itStart;
    deltaLength = geo::vect::norm(delta);
    endChunk.length -= deltaLength;
  }
  
  // now we have the two points at the side of the middle point
  assert(!(itStart == itEnd));
  assert(startChunk.length <= targetLength);
  assert(startChunk.length + deltaLength >= targetLength);
  assert(endChunk.length <= pathLength - targetLength);
  assert(endChunk.length + deltaLength >= pathLength - targetLength);
  
  // f is the relative position of the middlepoint between itStart and itEnd:
  // f=0 means itStart, f=1 means itEnd, f=0.5 means exact middle point
  double const f = (deltaLength != 0.0)
    ? ((targetLength - startChunk.length) / deltaLength): 0.0;
  assert(f >= 0.0);
  assert(f <= 1.0);
  
  return PathPointInfo_t<FIter>{
      pathLength   // length
    , deltaLength  // step
    , f            // frac
    , itStart      //itPoint
    };
  
} // util::findMiddlePointInPath()


// -----------------------------------------------------------------------------
template <typename FIter, typename LIter>
auto util::pathMiddlePoint
  (FIter itFirst, LIter itLast, double relTarget /* = 0.5 */)
{
  
  // We don't check that `LIter` also points to `Point_t`.
  // The return value will be `Point_t` anyway.
  
  PathPointInfo_t<FIter> const targetInfo
    = findMiddlePointInPath(itFirst, itLast, relTarget);
  
  // "next" point might be not there (e.g. if length is 0):
  // let's not use it if possible
  if (targetInfo.frac == 0.0) return *(targetInfo.itPoint);

  auto const delta = *std::next(targetInfo.itPoint) - *(targetInfo.itPoint);
  return *(targetInfo.itPoint) + delta * targetInfo.frac;
  
} // util::pathMiddlePoint()


// ---------------------------------------------------------------------------
template <typename BIter, typename EIter>
double util::pathLength(BIter begin, EIter end) {
  
  if (!(begin != end)) return 0.0;
  
  double totalLength = 0.0;
  BIter it = begin, prev = begin;
  while (++it != end) {
    totalLength += geo::vect::norm(*it - *prev);
    prev = it;
  }
  return totalLength;
  
} // util::pathLength()


// ---------------------------------------------------------------------------
template <typename BIter, typename EIter>
std::vector<double> util::partialPathLengths(BIter begin, EIter end) {
  if (!(begin != end)) return {};
  std::vector<double> lengths;
  lengths.reserve(std::distance(begin, end));
  double totalLength = 0.0;
  BIter it = begin, prev = begin;
  do {
    totalLength += geo::vect::norm(*it - *prev);
    lengths.push_back(totalLength);
    prev = it;
  } while (++it != end);
  return lengths;
} // util::partialPathLengths()


// -----------------------------------------------------------------------------
template <typename Iter, typename SideFunc>
std::pair<Iter, Iter> util::findCrossingSegment
  (Iter begin, Iter end, SideFunc sideOf)
{
  
  // find the place from start where the side first changes
  Iter itStart = begin;
  auto const startSide = sideOf(*itStart);
  while (++itStart != end) {
    if (sideOf(*itStart) != startSide) break; // changed side
  }
  
  if (itStart == end) return { end, end }; // path is not crossing side at all
  --itStart; // let have itStart point to the last point on the original side
  
  // find the place from end where the side last changes
  // (if it starts from the same side as the other end, it's weird but we
  //  pretend nothing happened)
  Iter itEnd = std::prev(end);
  auto const endSide = sideOf(*itEnd);
  while (--itEnd != itStart) {
    if (sideOf(*itEnd) != endSide) break; // changed side
  }
  ++itEnd; // ... and itEnd points to the last point on the end original side
  
  return { itStart, itEnd };
  
} // util::findCrossingSegment()


// -----------------------------------------------------------------------------


#endif // ICARUSALG_UTILITIES_TRAJECTORYUTILS_H
