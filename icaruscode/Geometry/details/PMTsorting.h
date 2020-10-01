/**
 * @file   icaruscode/Geometry/details/PMTsorting.h
 * @brief  Geometry obect sorter with PMT following TPC wire order.
 * @date   April 26, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/Geometry/details/PMTsorting.cxx
 */

#ifndef ICARUSCODE_GEOMETRY_DETAILS_PMTSORTING_H
#define ICARUSCODE_GEOMETRY_DETAILS_PMTSORTING_H


// LArSoft libraries
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/span.h" // util::span
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus { class PMTsorterStandard; }

/**
 * @brief Sorter sorting PMT to follow the same order as TPC (standard).
 * 
 * This class sorts the elements of the LArSoft detector description.
 * The "standard" algorithm (`geo::GeoObjectSorterStandard`) arranges TPC
 * elements with x, then z, then y. The PMT are arranged so that their channels
 * mimic the order of the TPC channels.
 * 
 * The algorithm for assigning channels to the PMT follows the criteria:
 * 
 * * PMT are ordered by increasing _x_ (related to drift direction);
 * * channels are assigned value ranges increasing with _x_ coordinate;
 * * within a wire plane, PMT number increases with its _z_ (beam direction)
 *   coordinate;
 * * in case of same _z_ (in ICARUS all PMT are in 2- or 3-PMT "towers"), an
 *   increasing _y_ order (geographical vertical, toward the sky) is chosen.
 * 
 * PMT channels are assigned by a fixed LArSoft algorithm, cryostat by cryostat
 * with increasing cryostat number (first `C:0`, then `C:1`, ...).
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * This sorter supports the following parameters:
 * 
 * * `ToleranceX`, `ToleranceY`, `ToleranceZ` (double, default: `1.0`):
 *   rounding used for sorting optical detector position, in centimeters;
 *   if the coordinate of two optical detectors are closer than the tolerance
 *   for that coordinate (absolute), they two detectors are considered to be at
 *   the same position in that coordinate
 *   (note that the default value is very generous).
 * 
 */
class icarus::PMTsorterStandard {
  
  /// List of optical detector pointers for sorting.
  using OpDetList_t = std::vector<geo::OpDetGeo>;
  
  /// Part of list of optical detector pointers for sorting.
  using OpDetSpan_t = util::span<OpDetList_t::iterator>;
  
    public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<double> ToleranceX {
      Name("ToleranceX"),
      Comment("tolerance when sorting optical detectors on x coordinate [cm]"),
      1.0 // default
      };
    
    fhicl::Atom<double> ToleranceY {
      Name("ToleranceY"),
      Comment("tolerance when sorting optical detectors on x coordinate [cm]"),
      1.0 // default
      };
    
    fhicl::Atom<double> ToleranceZ {
      Name("ToleranceZ"),
      Comment("tolerance when sorting optical detectors on x coordinate [cm]"),
      1.0 // default
      };
    
  }; // Config
  
  
  /// Constructor: passes the configuration to the base class.
  PMTsorterStandard(Config const& config)
    : fSmallerCenterX{ config.ToleranceX() }
    , fSmallerCenterY{ config.ToleranceY() }
    , fSmallerCenterZ{ config.ToleranceZ() }
    {}
  
  
  // @{
  /**
   * @brief Sorts the specified optical detectors.
   * @param opDets collection of pointers to all optical detectors in a cryostat
   * 
   * The collection `opDets` of optical detectors is sorted in place.
   * Sorting criteria are documented in `icarus::GeoObjectSorterPMTasTPC` class
   * documentation.
   * 
   * This algorithm requires all optical detectors to have their center defined
   * (`geo::OpDetGeo::GetCenter()`). No other information is used.
   * 
   * @note The current implementation is very sensitive to rounding errors!
   * 
   */
  void sort(std::vector<geo::OpDetGeo>& opDets) const;
  
  void operator() (std::vector<geo::OpDetGeo>& opDets) const { sort(opDets); }
  // @}
  
  
    private:
  
  /// `geo::OpDetGeo` comparer according to one coordinate of their center.
  /// Accomodates for some tolerance.
  template <double (geo::Point_t::*Coord)() const>
  struct OpDetGeoCenterCoordComparer {
    
    /// Object used for comparison; includes a tolerance.
    lar::util::RealComparisons<double> const fCmp;
    
    /// Constructor: fixes the tolerance for the comparison.
    OpDetGeoCenterCoordComparer(double tol = 0.0): fCmp(tol) {}
    
    /// Returns whether `A` has a center coordinate `Coord` smaller than `B`.
    bool operator() (geo::OpDetGeo const& A, geo::OpDetGeo const& B) const
      {
        return fCmp.strictlySmaller
            ((A.GetCenter().*Coord)(), (B.GetCenter().*Coord)());
      }
    
  }; // OpDetGeoCenterCoordComparer
  
  
  /// Sorting criterium according to _x_ coordinate of `geo::OpDetGeo` center.
  OpDetGeoCenterCoordComparer<&geo::Point_t::X> const fSmallerCenterX;
  
  /// Sorting criterium according to _y_ coordinate of `geo::OpDetGeo` center.
  OpDetGeoCenterCoordComparer<&geo::Point_t::Y> const fSmallerCenterY;
  
  /// Sorting criterium according to _z_ coordinate of `geo::OpDetGeo` center.
  OpDetGeoCenterCoordComparer<&geo::Point_t::Z> const fSmallerCenterZ;
  
  
  /// Sorts the `geo::OpDetGeo` assuming they belong to the same plane.
  void sortInPlane(OpDetSpan_t const& opDets) const;
  
}; // icarus::PMTsorterStandard


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_GEOMETRY_DETAILS_PMTSORTING_H
