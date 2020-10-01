/**
 * @file icaruscode/Geometry/GeoObjectSorterPMTasTPC.h
 * @brief  Geometry object sorter with PMT following TPC wire order.
 * @date   April 26, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/Geometry/GeoObjectSorterPMTasTPC.cxx
 */

#ifndef ICARUSCODE_GEOMETRY_GEOOBJECTSORTERPMTASTPC_H
#define ICARUSCODE_GEOMETRY_GEOOBJECTSORTERPMTASTPC_H

// ICARUS libraries
#include "icaruscode/Geometry/details/PMTsorting.h" // icarus::PMTsorterStandard

// LArSoft libraries
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/span.h" // util::span

// framework libraries
#include "fhiclcpp/types/Table.h"
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
 * mimic the order of the TPC channels (delegated to `icarus::PMTsorter`
 * algorithm).
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
 * * `OpDetSorter` (configuration table; default: empty): configures the
 *   PMT sorter object (see `icarus::PMTsorter` for details)
 * 
 */
class icarus::GeoObjectSorterPMTasTPC: public geo::GeoObjectSorterStandard {
  
  /// The sorting algorithm we use.
  using PMTsorter_t = icarus::PMTsorterStandard;
  
  using PMTsorterConfigTable = fhicl::Table<PMTsorter_t::Config>;
  
    public:
  
  /// Constructor: passes the configuration to the base class.
  GeoObjectSorterPMTasTPC(fhicl::ParameterSet const& pset)
    : geo::GeoObjectSorterStandard(pset)
    , fPMTsorter
      (PMTsorterConfigTable{ pset.get("OpDetSorter", fhicl::ParameterSet{}) }())
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
  virtual void SortOpDets(std::vector<geo::OpDetGeo>& opDets) const override
    { fPMTsorter.sort(opDets); }
  
  
  /// Custom ICARUS sorting of CRT.
  virtual void SortAuxDets(std::vector<geo::AuxDetGeo>& adgeo) const override;
  
  /// Custom ICARUS sorting of CRT submodules.
  virtual void SortAuxDetSensitive
    (std::vector<geo::AuxDetSensitiveGeo> & adsgeo) const override;
  
    private:
  
  PMTsorter_t fPMTsorter; ///< PMT sorting algorithm.
  
}; // icarus::GeoObjectSorterPMTasTPC


static_assert(std::is_base_of_v<geo::GeoObjectSorter, icarus::GeoObjectSorterPMTasTPC>);


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_GEOMETRY_GEOOBJECTSORTERPMTASTPC_H
