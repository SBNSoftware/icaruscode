/**
 * @file   icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h
 * @brief  Algorihtm to group PMTs into piling towers.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 7, 2020
 * @see    icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.cxx
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PMTVERTICALSLICINGALG_H
#define ICARUSCODE_PMT_ALGORITHMS_PMTVERTICALSLICINGALG_H

// LArSoft libraries
#include "larcorealg/CoreUtils/StdUtils.h" // util::begin(), util::end()
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Vector_t

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <vector>
#include <string>


// -----------------------------------------------------------------------------
// --- forward declarations
// ---
namespace geo {
  class GeometryCore;
  class CryostatGeo;
  class OpDetGeo;
} // namespace geo


//------------------------------------------------------------------------------
namespace icarus::trigger { class PMTverticalSlicingAlg; }
/**
 * @brief Algorithm clustering PMT according to their position.
 *
 * The algorithm groups the provided PMT by plane (determined from the drift
 * direction of the TPC) and by the beam-like direction
 * (the one that is not the vertical, i.e. TPC "width" direction).
 *
 */
class icarus::trigger::PMTverticalSlicingAlg {

    public:

  /// Type of optical detector list.
  using PMTlist_t = std::vector<geo::OpDetGeo const*>;

  /// Type of optical detector list in a PMT tower.
  using PMTtower_t = PMTlist_t;

  /// Type of list of PMT towers on a single optical detector plane.
  using PMTtowerOnPlane_t = std::vector<PMTtower_t>;

  /// Type of PMT towers, per plane.
  using Slices_t = std::vector<PMTtowerOnPlane_t>;


  /// Constructor: no configuration parameters so far.
  PMTverticalSlicingAlg(std::string logCategory = "PMTverticalSlicingAlg");


  /// Computes slices from all PMT in `cryo` and appends them to `slices`.
  void appendCryoSlices(Slices_t& slices, geo::CryostatGeo const& cryo) const;

  /**
   * @brief Groups optical detectors under the specified cryostat into walls.
   * @param cryo the cryostat containing the optical detectors
   * @return a list of pairs: each an absolute drift coordinate and PMT list
   * 
   * The algorithm returns a list of walls, with as first element a coordinate
   * representing the wall (drift coordinate) and the second the list of
   * optical detectors in that wall, in no particular order.
   * The walls are sorted by increasing drift coordinate.
   */
  std::vector<std::pair<double, PMTlist_t>> PMTwalls
    (geo::CryostatGeo const& cryo) const;

  /**
   * @brief Groups optical detectors in all the detector into walls.
   * @param geom the geometry description of the detector
   * @return a list of pairs: each an absolute drift coordinate and PMT list
   * @see `PMTwalls(geo::CryostatGeo const&) const`
   * 
   * PMT walls are extracted for each of the cryostats in the detector, and
   * the result is merged into a single collection.
   */
  std::vector<std::pair<double, PMTlist_t>> PMTwalls
    (geo::GeometryCore const& geom) const;


    private:
  
  std::string fLogCategory; ///< Category for message streaming.

  /// Computes slices and appends them to an existing list.
  void appendSlices(
    Slices_t& slices, PMTlist_t const& PMTs,
    geo::Vector_t const& planeNorm, geo::Vector_t const& clusterDir
    ) const;

  /**
   * @brief Groups the specifies optical detectors into walls.
   * @param PMTs the list of PMT to be grouped
   * @param dir the direction normal to the walls
   * @return a list of pairs: each an absolute wall coordinate and PMT list
   * 
   * The algorithm returns a list of walls, with as first element a coordinate
   * representing the wall and the second the list of optical detectors in that
   * wall, in no particular order.
   * The walls are sorted by increasing coordinate.
   */
  std::vector<std::pair<double, PMTlist_t>> PMTwalls
    (PMTlist_t const& PMTs, geo::Vector_t const& dir) const;


  /// Returns an (arbitrary) coordinate along `dir` representing the PMT list.
  static double PMTwallPosition
    (PMTlist_t const& PMTs, geo::Vector_t const& dir);


  template <typename TPCCont>
  static geo::Vector_t determineDriftDir(TPCCont const& TPCcont);

  template <typename TPCCont>
  static geo::Vector_t determineLengthDir(TPCCont const& TPCcont);

  /// Clusters the PMTs along the specified direction.
  static std::vector<PMTlist_t> clusterPMTby
    (PMTlist_t const& PMTs, geo::Vector_t const& dir);

  /// Returns a list of all `geo::OpDetGeo` in cryostat `cryo`.
  static PMTlist_t getCryoPMTs(geo::CryostatGeo const& cryo);

  /// Returns a list of all `geo::OpDetGeo` in the whole geometry.
  static PMTlist_t getPMTs(geo::GeometryCore const& geom);

  /// Returns the geometric center of the PMT `wall`.
  static geo::Point_t PMTwallCenter(PMTtowerOnPlane_t const& wall);

  // pretty sure this is already in some feature branch if not in LArSoft...
  /// Returns whether `a` and `b` are parallel.
  static bool areParallel(geo::Vector_t const& a, geo::Vector_t const& b);

}; // class icarus::trigger::PMTverticalSlicingAlg


//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
template <typename TPCCont>
geo::Vector_t icarus::trigger::PMTverticalSlicingAlg::determineDriftDir
  (TPCCont const& TPCcont)
{

  auto iTPC = util::begin(TPCcont);
  auto const tend = util::end(TPCcont);

  if (iTPC == tend) {
    throw cet::exception("PMTverticalSlicingAlg")
      << "determineDriftDir(): no TPC in the specified object!";
  }

  auto const dir = iTPC->DriftDir();
  while(++iTPC != tend) {
    if (areParallel(iTPC->DriftDir(), dir))
      continue;
    throw cet::exception("PMTverticalSlicingAlg")
      << "determineDriftDir(): TPC " << iTPC->ID()
      << " has a drift direction incompatible with " << dir << ".\n";
  } // while

  return dir;

} // icarus::trigger::PMTverticalSlicingAlg::determineDriftDir()


//------------------------------------------------------------------------------
template <typename TPCCont>
geo::Vector_t icarus::trigger::PMTverticalSlicingAlg::determineLengthDir
  (TPCCont const& TPCcont)
{

  auto iTPC = util::begin(TPCcont);
  auto const tend = util::end(TPCcont);

  if (iTPC == tend) {
    throw cet::exception("PMTverticalSlicingAlg")
      << "determineLengthDir(): no TPC in the specified object!";
  }

  auto const dir = iTPC->LengthDir();
  while(++iTPC != tend) {
    if (areParallel(iTPC->LengthDir(), dir))
      continue;
    throw cet::exception("PMTverticalSlicingAlg")
      << "determineLengthDir(): TPC " << iTPC->ID()
      << " has a width direction incompatible with " << dir << ".\n";
  } // while

  return dir;

} // icarus::trigger::PMTverticalSlicingAlg::determineLengthDir()


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_ALGORITHMS_PMTVERTICALSLICINGALG_H
