/**
 * @file   icaruscode/Geometry/details/AuxDetSorting.h
 * @brief  Functions for sorting ICARUS CRT modules (auxiliary detectors).
 * @author Chris Hilgenberg, Gianluca Petrillo (refactoring only)
 * @date   August 7, 2018
 * @see    icaruscode/Geometry/details/AuxDetSorting.cxx
 */

#ifndef ICARUSCODE_GEOMETRY_DETAILS_AUXDETSORTING_H
#define ICARUSCODE_GEOMETRY_DETAILS_AUXDETSORTING_H


// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

// C/C++ standard libraries
#include <vector>


namespace icarus {
  
  //----------------------------------------------------------------------------
  /// Sorts ICARUS CRT modules in standard configuration.
  void SortAuxDetsStandard(std::vector<geo::AuxDetGeo>& adgeo);
  
  
  //----------------------------------------------------------------------------
  /// Sorts ICARUS CRT submodules in standard configuration.
  void SortAuxDetSensitiveStandard
    (std::vector<geo::AuxDetSensitiveGeo>& adsgeo);
  
  
  //----------------------------------------------------------------------------
  
} // namespace icarus


#endif // ICARUSCODE_GEOMETRY_DETAILS_AUXDETSORTING_H
