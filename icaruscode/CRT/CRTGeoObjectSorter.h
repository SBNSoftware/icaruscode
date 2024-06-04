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

#include "larcorealg/Geometry/AuxDetGeoObjectSorter.h"
#include "fhiclcpp/fwd.h"

//namespace icarus {
//namespace crt {
namespace geo{

  class CRTGeoObjectSorter : public AuxDetGeoObjectSorter {
  public:
    CRTGeoObjectSorter();
    CRTGeoObjectSorter(fhicl::ParameterSet const&);

  private:
    bool compareAuxDets(AuxDetGeo const& ad1, AuxDetGeo const& ad2) const override;
    bool compareAuxDetSensitives(AuxDetSensitiveGeo const& ads1,
                                 AuxDetSensitiveGeo const& ads2) const override;
  };

} //namespace crt
//} //namespace icarus

#endif  // ICARUS_CRTGeoObjectSorter_h
