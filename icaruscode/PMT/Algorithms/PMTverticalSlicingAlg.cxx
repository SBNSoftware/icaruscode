/**
 * @file   icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.cxx
 * @brief  Algorihtm to group PMTs into piling towers (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 7, 2020
 * @see    icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h
 */

// library header
#include "icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h"
#include "icarusalg/Utilities/SimpleClustering.h" // util::clusterBy()
#include "icarusalg/Utilities/sortLike.h"

// LArSoft libraries
#include "lardataalg/Utilities/StatCollector.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/RealComparisons.h" // lar::util::Vector3DComparisons
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// C/C++ standard library
#include <vector>
#include <functional> // std::less
#include <iterator> // std::make_move_iterator(), std::back_inserter()
#include <limits>
#include <cmath> // std::abs()
#include <cstdint> // std::int_fast8_t
#include <cassert>


//------------------------------------------------------------------------------
struct StandardLArSoftGeometrySorter {
  
  bool operator() (geo::Point_t const& a, geo::Point_t const& b) const
    { return cmp(a, b) < 0; }
  
  static std::int_fast8_t cmp(double a, double b)
    {
      if (a < b) return -1; // std::strong_ordering::less
      if (a > b) return +1; // std::strong_ordering::greater
      return 0;             // std::strong_ordering::equal
    }
  
  static std::int_fast8_t cmp(geo::Point_t const& a, geo::Point_t const& b)
    {
      // C++20: use three-way comparison
      if (int const res = cmp(a.X(), b.X())) return res;
      if (int const res = cmp(b.Z(), a.Z())) return res;
      return cmp(a.Y(), b.Y());
    } // cmp()
  
}; // StandardLArSoftGeometrySorter

//------------------------------------------------------------------------------
namespace {
  
  // Use the projection of the PMT center on the specified direction as
  // clustering key; we need a reference point to project... we pick `origin()`.
  struct PMTprojectorClass {
    geo::Vector_t const dir; ///< Projection direction.
    double operator() (geo::OpDetGeo const* opDet) const
      { return (opDet->GetCenter() - geo::origin()).Dot(dir); };
  }; // PMTprojectorClass
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /**
   * @brief Moves the content of `source` into the end of `dest`.
   * @tparam T type of objects contained in the collections
   * @tparam Src type of source collection
   * @param dest vector to be added elements
   * @param source container with all the elements to be appended to `dest`
   * @return a reference to `dest` itself
   * 
   * The content of `source` is moved away, and `source` itself is moved.
   */
  template <typename T, typename Src>
  std::vector<T>& append(std::vector<T>& dest, Src&& src) {
    
    using std::size;
    using std::begin;
    using std::end;
    
    dest.reserve(dest.size() + size(src));
    dest.insert(dest.end(),
      std::make_move_iterator(begin(src)), std::make_move_iterator(end(src))
      );
    Src{ std::move(src) }; // move src into a temporary (which we don't care of)
    return dest;
  } // append()
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::trigger::PMTverticalSlicingAlg
//------------------------------------------------------------------------------
icarus::trigger::PMTverticalSlicingAlg::PMTverticalSlicingAlg
 (std::string logCategory /* = "PMTverticalSlicingAlg" */)
 : fLogCategory(std::move(logCategory))
 {}


//------------------------------------------------------------------------------
void icarus::trigger::PMTverticalSlicingAlg::appendCryoSlices
  (Slices_t& slices, geo::CryostatGeo const& cryo) const
{
  geo::Vector_t const& driftDir = determineDriftDir(cryo.IterateTPCs());
  geo::Vector_t const& widthDir = determineLengthDir(cryo.IterateTPCs());
  appendSlices(slices, getCryoPMTs(cryo), driftDir, widthDir);
} // icarus::trigger::PMTverticalSlicingAlg::appendCryoSlices()


//------------------------------------------------------------------------------
auto icarus::trigger::PMTverticalSlicingAlg::PMTwalls
  (geo::CryostatGeo const& cryo) const
  -> std::vector<std::pair<double, PMTlist_t>>
{
  return PMTwalls(getCryoPMTs(cryo), determineDriftDir(cryo.IterateTPCs()));
} // icarus::trigger::PMTverticalSlicingAlg::PMTwalls(CryostatGeo)


//------------------------------------------------------------------------------
auto icarus::trigger::PMTverticalSlicingAlg::PMTwalls
  (geo::GeometryCore const& geom) const
  -> std::vector<std::pair<double, PMTlist_t>>
{
  return PMTwalls(getPMTs(geom), determineDriftDir(geom.Iterate<geo::TPCGeo>()));
} // icarus::trigger::PMTverticalSlicingAlg::PMTwalls(GeometryCore)


//------------------------------------------------------------------------------
void icarus::trigger::PMTverticalSlicingAlg::appendSlices(
  Slices_t& slices, PMTlist_t const& PMTs,
  geo::Vector_t const& planeNorm, geo::Vector_t const& clusterDir
  ) const
{

  /*
   * 1. group PMT by plane (cluster on drift direction)
   * 2. group PMT in each plane by width coordinate (cluster by width direction)
   * 3. sort them by central coordinate using LArSoft standards
   */

  //
  // 1. group PMT by plane (cluster on drift direction)
  //
  
  std::vector<PMTlist_t> const PMTplanes = clusterPMTby(PMTs, planeNorm);

  // BEGIN debug
  {
    mf::LogTrace log { fLogCategory };
    log << PMTs.size() << " PMTs grouped into " << PMTplanes.size()
      << " planes in direction " << planeNorm << ":";
    for (auto&& [ iPlane, PMTplane ]: util::enumerate(PMTplanes)) {
      log << "\n [#" << iPlane << "] " << PMTplane.size() << " PMT:";
      for (geo::OpDetGeo const* opDet: PMTplane)
        log << " <" << opDet->ID() << ">";
    } // for
  }
  // END debug

  //
  // 2. group PMT in each plane by width coordinate (cluster by width direction)
  //
  unsigned int NClusteredPMTs [[maybe_unused]] = 0U;
  for (auto const& PMTplane: PMTplanes) {
    slices.push_back(clusterPMTby(PMTplane, clusterDir));

    // BEGIN debug
    auto const& planeSlices = slices.back();
    mf::LogTrace log { fLogCategory };
    log << PMTplane.size() << " PMTs in the plane grouped into "
      << planeSlices.size() << " towers in direction " << clusterDir << ":";
    for (auto&& [ iSlice, slice ]: util::enumerate(planeSlices)) {
      NClusteredPMTs += slice.size();
      log << "\n [#" << iSlice << "] " << slice.size() << " PMT:";
      for (geo::OpDetGeo const* opDet: slice) log << " <" << opDet->ID() << ">";
    } // for planeSlices
    // END debug

  } // for planes
  
  assert(NClusteredPMTs == PMTs.size());
  
  //
  // 3. sort them by central coordinate using LArSoft standards
  //
  // the decision here is to sort globally, including the existing entries
  mf::LogTrace{ fLogCategory } << "Reordering " << slices.size() << " planes";
  std::vector<geo::Point_t> wallCenters;
  for (PMTtowerOnPlane_t const& wall: slices) {
    wallCenters.push_back(PMTwallCenter(wall));
    mf::LogTrace{ fLogCategory }
      << " - wall from " << wall.size() << " towers at " << wallCenters.back();
  }
  util::sortLike(
    slices.begin(), slices.end(), wallCenters.cbegin(), wallCenters.cend(),
    StandardLArSoftGeometrySorter{}
    );

  {
    // BEGIN debug
    mf::LogTrace log { fLogCategory };
    log << PMTs.size() << " PMTs grouped into " << slices.size() << " planes:";
    for (auto const& [ iWall, wall ]: util::enumerate(slices)) {
      log << "\n  [#" << iWall << "] " << wall.size() << " towers:";
      for (auto&& [ iSlice, slice ]: util::enumerate(wall)) {
        log << "\n  [#" << iWall << ":" << iSlice << "] "
          << slice.size() << " PMT:";
        for (geo::OpDetGeo const* opDet: slice)
          log << " <" << opDet->ID() << ">";
      } // for wall
    } // for slices
    // END debug
  }
  
} // icarus::trigger::PMTverticalSlicingAlg::appendSlices()


//------------------------------------------------------------------------------
auto icarus::trigger::PMTverticalSlicingAlg::PMTwalls
  (PMTlist_t const& PMTs, geo::Vector_t const& dir) const
  -> std::vector<std::pair<double, PMTlist_t>>
{
  // fill the walls with PMTs
  std::vector<PMTlist_t> walls = clusterPMTby(PMTs, dir);
  
  // determine their coordinate
  std::vector<std::pair<double, PMTlist_t>> posAndWalls;
  posAndWalls.reserve(walls.size());
  for (PMTlist_t& wall: walls)
    posAndWalls.emplace_back(PMTwallPosition(wall, dir), std::move(wall));
  
  // sort by coordinate
  std::sort(posAndWalls.begin(), posAndWalls.end()); // sort by std::pair::first
  
  return posAndWalls;
  
} // icarus::trigger::PMTverticalSlicingAlg::PMTwalls(CryostatGeo)


//------------------------------------------------------------------------------
double icarus::trigger::PMTverticalSlicingAlg::PMTwallPosition
  (PMTlist_t const& PMTs, geo::Vector_t const& dir)
{
  if (PMTs.empty()) return std::numeric_limits<double>::lowest(); // arbitrary
  
  PMTprojectorClass const PMTcenterProjection { dir };
  lar::util::StatCollector<double> stats;
  for (geo::OpDetGeo const* PMT: PMTs) stats.add(PMTcenterProjection(PMT));
  return stats.Average();
  
} // icarus::trigger::PMTverticalSlicingAlg::PMTwallPosition()


//------------------------------------------------------------------------------
auto icarus::trigger::PMTverticalSlicingAlg::clusterPMTby
  (PMTlist_t const& PMTs, geo::Vector_t const& dir) -> std::vector<PMTlist_t>
{
#if 0
  // use the projection of the PMT center on the specified direction as
  // clustering key; we need a reference point to project... we pick `origin()`.
  auto const PMTcenterProjection = [dir](geo::OpDetGeo const* opDet)
    { return (opDet->GetCenter() - geo::origin()).Dot(dir); };
#endif // 0
  
  PMTprojectorClass const PMTcenterProjection { dir };
  
  // PMTs with a coordinate within 5 cm (kind of) will be clustered together
  constexpr double tol = 5.0; // cm
  auto const closeEnough
    = [](double a, double b){ return std::abs(a - b) < tol; };

  return util::clusterBy(PMTs, PMTcenterProjection, closeEnough, std::less<>());

} // icarus::trigger::PMTverticalSlicingAlg::clusterPMTby()


//------------------------------------------------------------------------------
auto icarus::trigger::PMTverticalSlicingAlg::getCryoPMTs
  (geo::CryostatGeo const& cryo) -> PMTlist_t
{
  PMTlist_t opDets;
  for (auto iOpDet: util::counter(cryo.NOpDet()))
    opDets.push_back(&cryo.OpDet(iOpDet));
  return opDets;
} // icarus::trigger::PMTverticalSlicingAlg::getCryoPMTs()


//------------------------------------------------------------------------------
auto icarus::trigger::PMTverticalSlicingAlg::getPMTs
  (geo::GeometryCore const& geom) -> PMTlist_t
{
  PMTlist_t opDets;
  for (geo::CryostatGeo const& cryo: geom.Iterate<geo::CryostatGeo>())
    append(opDets, getCryoPMTs(cryo));
  return opDets;
} // icarus::trigger::PMTverticalSlicingAlg::getCryoPMTs()


//------------------------------------------------------------------------------
bool icarus::trigger::PMTverticalSlicingAlg::areParallel
  (geo::Vector_t const& a, geo::Vector_t const& b)
{
  lar::util::Vector3DComparison cmp { 1e-4 };
  return cmp.zero(a.Cross(b));
} // icarus::trigger::PMTverticalSlicingAlg::areParallel()


//------------------------------------------------------------------------------
geo::Point_t icarus::trigger::PMTverticalSlicingAlg::PMTwallCenter
  (PMTtowerOnPlane_t const& wall)
{
  geo::vect::MiddlePointAccumulator wallCenter;
  for (PMTtower_t const& tower: wall) {
    for (geo::OpDetGeo const* PMT: tower) {
      wallCenter.add(PMT->GetCenter());
    } // for PMTs
  } // for towers
  assert(!wallCenter.empty()); // don't know what to do with an empty plane
  return wallCenter.middlePoint();
} // icarus::trigger::PMTverticalSlicingAlg::PMTwallCenter()


//------------------------------------------------------------------------------
