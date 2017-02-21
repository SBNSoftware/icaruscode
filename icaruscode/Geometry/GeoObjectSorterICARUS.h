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

#include "larcore/Geometry/GeoObjectSorterStandard.h"

namespace geo{

  class GeoObjectSorterICARUS : public GeoObjectSorterStandard {

  public:

    GeoObjectSorterICARUS(fhicl::ParameterSet const& p);
    ~GeoObjectSorterICARUS();

    void SortCryostats      (std::vector<geo::CryostatGeo*>        & cgeo)     const;
    void SortTPCs     	    (std::vector<geo::TPCGeo*>      	   & tgeo)     const;
    void SortPlanes   	    (std::vector<geo::PlaneGeo*>    	   & pgeo,	      
		      	     geo::DriftDirection_t     	     const & driftDir) const;
    void SortWires    	    (std::vector<geo::WireGeo*>     	   & wgeo)     const;
    
  private:
    
  };

}

#endif // GEO_GEOOBJECTSORTERICARUS_H
