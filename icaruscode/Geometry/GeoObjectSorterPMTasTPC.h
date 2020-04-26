/**
 * @file icaruscode/Geometry/GeoObjectSorterPMTasTPC.h
 * @brief  Geometry obect sorter with PMT following TPC wire order.
 * @date   April 26, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/Geometry/GeoObjectSorterPMTasTPC.cxx
 */

#ifndef ICARUSCODE_GEOMETRY_GEOOBJECTSORTERPMTASTPC_H
#define ICARUSCODE_GEOMETRY_GEOOBJECTSORTERPMTASTPC_H

// LArSoft libraries
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/span.h" // util::span

// framework libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus { class GeoObjectSorterPMTasTPC; }

/**
 * @brief Geometry sorter having PMT channels follow the same order as TPC.
 * 
 * This class sorts the elements of the LArSoft detector description.
 * The TPC elements are sorted according to the "standard" algorithm
 * (`geo::GeoObjectSorterStandard`). The PMT are arranged so that their channels
 * mimic the order of the TPC channels.
 * 
 * The algorithm for assigning channels to the wires follows the criteria:
 * 
 * * TPC are ordered by increasing _x_ (related to drift direction);
 * * channels are assigned value ranges increasing with the TPC number,
 *   i.e. with increasing _x_ coordinate;
 * * within a wire plane, channel number increases with the _z_ (beam direction)
 *   coordinate of the wire(s) behind the channel;
 * * in case of same _z_ (as for ICARUS first induction plane), an increasing
 *   _y_ order (geographical vertical, toward the sky) is chosen.
 * 
 * PMT channels are assigned by a fixed LArSoft algorithm, cryostat by cryostat
 * with increasing cryostat number (first `C:0`, then `C:1`, ...).
 * Each cryostat has its set of optical detectors, sorted by a customizable
 * geometric sorting algorithm, and the channel number assignment follows the
 * sequence of optical detectors as sorted by that algorithm.
 * 
 * This class reimplements the geometric sorting algorithm following criteria
 * similar to the TPC wires:
 * 
 * * optical detectors are split by plane (_x_ direction);
 * * starting with the plane with lower _x_, optical detectors are sorted
 *   by _z_ coordinate, then by _y_ coordinate.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * In addition to the parameters for the standard sorter
 * (`geo::GeoObjectSorterStandard`), this sorter supports the following
 * parameters:
 * 
 * * `ToleranceX`, `ToleranceY`, `ToleranceZ` (double, default: `1.0`):
 *   rounding used for sorting optical detector position, in centimeters;
 *   if the coordinate of two optical detectors are closer than the tolerance
 *   for that coordinate (absolute), they two detectors are considered to be at
 *   the same position in that coordinate
 *   (note that the default value is very generous).
 * 
 */
class icarus::GeoObjectSorterPMTasTPC: public geo::GeoObjectSorterStandard {
  
  /// List of optical detector pointers for sorting.
  using OpDetList_t = std::vector<geo::OpDetGeo>;
  
  /// Part of list of optical detector pointers for sorting.
  using OpDetSpan_t = util::span<OpDetList_t::iterator>;
  
    public:
  
  /*
  // not using configuration validation
  // because we could not pass the configuration to the base class
  
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
  */
  
  
  /// Constructor: passes the configuration to the base class.
  GeoObjectSorterPMTasTPC(fhicl::ParameterSet const& pset)
    : geo::GeoObjectSorterStandard(pset)
    , fSmallerCenterX{ pset.get("ToleranceX", 1.0) }
    , fSmallerCenterY{ pset.get("ToleranceY", 1.0) }
    , fSmallerCenterZ{ pset.get("ToleranceZ", 1.0) }
    {}
  
  
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
  virtual void SortOpDets(std::vector<geo::OpDetGeo>& opDets) const override;
  
  
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
  void sortOpDetsInPlane(OpDetSpan_t const& opDets) const;
  
  
}; // icarus::GeoObjectSorterPMTasTPC


static_assert(std::is_base_of_v<geo::GeoObjectSorter, icarus::GeoObjectSorterPMTasTPC>);


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_GEOMETRY_GEOOBJECTSORTERPMTASTPC_H
