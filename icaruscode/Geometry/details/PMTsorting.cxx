/**
 * @file   icaruscode/Geometry/details/PMTsorting.cxx
 * @brief  PMT sorting functions for ICARUS.
 * @date   April 26, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/Geometry/details/PMTsorting.h
 */


// library header
#include "icaruscode/Geometry/details/PMTsorting.h"

// LArSoft libraries
// #include "larcorealg/CoreUtils/span.h"

// C/C++ standard libraries
#include <vector>
#include <iterator> // std::next(), std::prev()
#include <algorithm> // std::sort(), std::stable_sort()
#include <cassert>


// -----------------------------------------------------------------------------
void icarus::PMTsorterStandard::sort(std::vector<geo::OpDetGeo>& opDets) const {
  assert(opDets.size() % 2 == 0); // must be even!
  
  /*
   * 1. sort all optical detectors by _x_
   * 2. split them by plane
   * 3. sort the detectors within each plane.
   */
  
  //
  // 1. sort all optical detectors by _x_
  //
  std::sort(begin(opDets), end(opDets), fSmallerCenterX);
  
  //
  // 2. split them by plane: we take a horrible shortcut here...
  //
  
  std::vector<OpDetSpan_t> OpDetsPerPlane;
  
  // just split the list in two
  auto const middle = std::next(opDets.begin(), opDets.size() / 2U);
  assert(fSmallerCenterX(*std::prev(middle), *middle));
  
  OpDetsPerPlane.emplace_back(opDets.begin(), middle);
  OpDetsPerPlane.emplace_back(middle, opDets.end());
  assert(OpDetsPerPlane[0].size() == OpDetsPerPlane[1].size());
  
  //
  // 3. sort the detectors within each plane.
  //
  for (auto const& planeOpDets: OpDetsPerPlane)
    sortInPlane(util::make_span(planeOpDets));
  
  // all done in place, no return
  
} // icarus::PMTsorterStandard::sort()


// -----------------------------------------------------------------------------
void icarus::PMTsorterStandard::sortInPlane(OpDetSpan_t const& opDets) const {
  /*
   * 0. assume it's already sorted by _x_
   * 1. sort by vertical coordinate (_y_)
   * 2. stable sort by beam coordinate (_z_)
   */
  
  // 1. sort by vertical coordinate (_y_)
  std::stable_sort(begin(opDets), end(opDets), fSmallerCenterY);
  
  // 2. stable sort by beam coordinate (_z_)
  std::stable_sort(begin(opDets), end(opDets), fSmallerCenterZ);
  
} // icarus::PMTsorterStandard::sortInPlane()


// -----------------------------------------------------------------------------
