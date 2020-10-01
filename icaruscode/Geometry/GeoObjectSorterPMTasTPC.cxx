/**
 * @file icaruscode/Geometry/GeoObjectSorterPMTasTPC.cxx
 * @brief  Geometry object sorter with PMT following TPC wire order.
 * @date   April 26, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/Geometry/GeoObjectSorterPMTasTPC.h
 * 
 * Nothing, really.
 */


// library header
#include "icaruscode/Geometry/GeoObjectSorterPMTasTPC.h"

// ICARUS libraries
#include "icaruscode/Geometry/details/AuxDetSorting.h"


//------------------------------------------------------------------------------
void icarus::GeoObjectSorterPMTasTPC::SortAuxDets
  (std::vector<geo::AuxDetGeo>& adgeo) const
{
  icarus::SortAuxDetsStandard(adgeo);
}

//------------------------------------------------------------------------------
void icarus::GeoObjectSorterPMTasTPC::SortAuxDetSensitive
  (std::vector<geo::AuxDetSensitiveGeo>& adsgeo) const
{
  icarus::SortAuxDetSensitiveStandard(adsgeo);
}


//------------------------------------------------------------------------------
