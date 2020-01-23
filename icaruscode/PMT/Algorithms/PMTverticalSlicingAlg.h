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


    private:
  std::string fLogCategory; ///< Category for message streaming.

  /// Computes slices and appends them to an existing list.
  void appendSlices(
    Slices_t& slices, PMTlist_t const& PMTs,
    geo::Vector_t const& planeNorm, geo::Vector_t const& clusterDir
    ) const;

  template <typename TPCCont>
  static geo::Vector_t determineDriftDir(TPCCont const& TPCcont);

  template <typename TPCCont>
  static geo::Vector_t determineLengthDir(TPCCont const& TPCcont);

  /// Clusters the PMTs along the specified direction.
  static std::vector<PMTlist_t> clusterPMTby
    (PMTlist_t const& PMTs, geo::Vector_t const& dir);

  /// Returns a list of all `geo::OpDetGeo` in cryostat `cryo`.
  static PMTlist_t getCryoPMTs(geo::CryostatGeo const& cryo);

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

  auto const dir = iTPC->template DriftDir<geo::Vector_t>();
  while(++iTPC != tend) {
    if (areParallel(iTPC->template DriftDir<geo::Vector_t>(), dir))
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

  auto const dir = iTPC->template LengthDir<geo::Vector_t>();
  while(++iTPC != tend) {
    if (areParallel(iTPC->template LengthDir<geo::Vector_t>(), dir))
      continue;
    throw cet::exception("PMTverticalSlicingAlg")
      << "determineLengthDir(): TPC " << iTPC->ID()
      << " has a width direction incompatible with " << dir << ".\n";
  } // while

  return dir;

} // icarus::trigger::PMTverticalSlicingAlg::determineLengthDir()


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_ALGORITHMS_PMTVERTICALSLICINGALG_H
