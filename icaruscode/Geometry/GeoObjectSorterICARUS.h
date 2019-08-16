////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterICARUS.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  wketchum@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERICARUS_H
#define GEO_GEOOBJECTSORTERICARUS_H

#include <vector>

#include "fhiclcpp/fwd.h"

#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class GeoObjectSorterICARUS : public GeoObjectSorter {

  public:

    GeoObjectSorterICARUS(fhicl::ParameterSet const& p);
    ~GeoObjectSorterICARUS();

    void SortAuxDets        (std::vector<geo::AuxDetGeo*>          & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo)   const;
    void SortCryostats      (std::vector<geo::CryostatGeo*>        & cgeo)     const;
    void SortTPCs     	    (std::vector<geo::TPCGeo*>      	   & tgeo)     const;
    void SortPlanes   	    (std::vector<geo::PlaneGeo*>    	   & pgeo,	      
		      	     geo::DriftDirection_t     	     const & driftDir) const;
    void SortWires    	    (std::vector<geo::WireGeo*>     	   & wgeo)     const;
    
  private:
    
  };

}

#endif // GEO_GEOOBJECTSORTERICARUS_H
