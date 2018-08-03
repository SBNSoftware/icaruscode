///////////////////////////////////////////////////////////////////////////////
/// \file CRTGeoObjectSorter.h
/// \brief Interface to algorithm class for sorting of AuxDetGeo objects
///
/// Originally ported from AuxDetGeoObjectSorterLArIAT.h (Author: brebel@fnal.gov)
/// and modified for SBND (Author: mastbaum@uchicago.edu)
///
/// \version $Id: 0.0 $
/// \author chilge@rams.colostate.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef ICARUS_CRTGeoObjectSorter_h
#define ICARUS_CRTGeoObjectSorter_h

#include <vector>

#include "larcorealg/Geometry/AuxDetGeoObjectSorter.h"

//namespace icarus {
//namespace crt {
namespace geo{

  class CRTGeoObjectSorter : public AuxDetGeoObjectSorter {
  public:

    CRTGeoObjectSorter(fhicl::ParameterSet const& p);

    ~CRTGeoObjectSorter();

    void SortAuxDets (std::vector<geo::AuxDetGeo*>& adgeo) const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*>& adsgeo) const;

  };

} //namespace crt
//} //namespace icarus

#endif  // ICARUS_CRTGeoObjectSorter_h

