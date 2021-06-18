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

// LArSoft libraries
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/RealComparisons.h" // lar::util::Vector3DComparisons

// C/C++ standard library
#include <vector>
#include <functional> // std::less
#include <cmath> // std::abs()
#include <cassert>


//------------------------------------------------------------------------------
icarus::trigger::PMTverticalSlicingAlg::PMTverticalSlicingAlg
 (std::string logCategory /* = "PMTverticalSlicingAlg" */)
 : fLogCategory(std::move(logCategory))
 {}


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
void icarus::trigger::PMTverticalSlicingAlg::appendCryoSlices
  (Slices_t& slices, geo::CryostatGeo const& cryo) const
{
  geo::Vector_t const& driftDir = determineDriftDir(cryo.IterateTPCs());
  geo::Vector_t const& widthDir = determineLengthDir(cryo.IterateTPCs());
  appendSlices(slices, getCryoPMTs(cryo), driftDir, widthDir);
} // icarus::trigger::PMTverticalSlicingAlg::appendCryoSlices()


//------------------------------------------------------------------------------
void icarus::trigger::PMTverticalSlicingAlg::appendSlices(
  Slices_t& slices, PMTlist_t const& PMTs,
  geo::Vector_t const& planeNorm, geo::Vector_t const& clusterDir
  ) const
{

  /*
   * 1. group PMT by plane (cluster on drift direction)
   * 2. group PMT in each plane by width coordinate (cluster by width direction)
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
  unsigned int NClusteredPMTs [[gnu::unused]] = 0U;
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

} // icarus::trigger::PMTverticalSlicingAlg::appendSlices()


//------------------------------------------------------------------------------
auto icarus::trigger::PMTverticalSlicingAlg::clusterPMTby
  (PMTlist_t const& PMTs, geo::Vector_t const& dir) -> std::vector<PMTlist_t>
{
  // use the projection of the PMT center on the specified direction as
  // clustering key; we need a reference point to project... we pick `origin()`.
  auto const PMTcenterProjection = [dir](geo::OpDetGeo const* opDet)
    { return (opDet->GetCenter() - geo::origin()).Dot(dir); };

  // PMTs with a coordinate within 5 cm (kind of) will be clustered together
  constexpr double tol = 5.0; // cm
  auto const closeEnough
    = [](double a, double b){ return std::abs(a - b) < tol; };

  return util::clusterBy(PMTs, PMTcenterProjection, closeEnough, std::less<>());

} // icarus::trigger::PMTverticalSlicingAlg::clusterPMTby()


//------------------------------------------------------------------------------
bool icarus::trigger::PMTverticalSlicingAlg::areParallel
  (geo::Vector_t const& a, geo::Vector_t const& b)
{
  lar::util::Vector3DComparison cmp { 1e-4 };
  return cmp.zero(a.Cross(b));
} // icarus::trigger::PMTverticalSlicingAlg::areParallel()


//------------------------------------------------------------------------------
