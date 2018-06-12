////////////////////////////////////////////////////////////////////////
/// \file CRTGeoObjectSorter.cxx
/// \brief Interface to algorithm class for sorting of AuxDetGeo objects
///
/// Originally ported from AuxDetGeoObjectSorterLArIAT.h (Author: brebel@fnal.gov)
/// and modified for SBND (Author: mastbaum@uchicago.edu) then ICARUS
///
/// \version $Id: 0.0 $
/// \author chilge@rams.colostate.edu
////////////////////////////////////////////////////////////////////////

#include "icaruscode/CRT/CRTGeoObjectSorter.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

namespace geo {

  //----------------------------------------------------------------------------
  // Define sort order for AuxDets in standard configuration
  static bool sortAuxDetICARUS(const AuxDetGeo* ad1, const AuxDetGeo* ad2) {

    // sort based off of GDML name, module number
    std::string ad1name = (ad1->TotalVolume())->GetName();
    std::string ad2name = (ad2->TotalVolume())->GetName();
    // assume volume name is "volAuxDet_Module_###_<region>"
    std::string base = "volAuxDet_Module_";

    int ad1Num = atoi( ad1name.substr( base.size(), 3).c_str() );
    int ad2Num = atoi( ad2name.substr( base.size(), 3).c_str() );
    
    return ad1Num < ad2Num;
  }

  //----------------------------------------------------------------------------
  // Define sort order for AuxDetSensitives in standard configuration
  static bool sortAuxDetSensitiveICARUS(const AuxDetSensitiveGeo* ad1,
                                      const AuxDetSensitiveGeo* ad2)
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string ad1name = (ad1->TotalVolume())->GetName();
    std::string ad2name = (ad2->TotalVolume())->GetName();
    // assume volume name is "volAuxDetSensitive_Module_###_Strip_##"
    std::string baseMod = "volAuxDetSensitive_Module_";
    std::string baseStr = "volAuxDetSensitive_Module_###_Strip_";

    int ad1Num = atoi( ad1name.substr( baseMod.size(), 3).c_str() );
    int ad2Num = atoi( ad2name.substr( baseMod.size(), 3).c_str() );

    if(ad1Num!=ad2Num) return ad1Num < ad2Num; 

    ad1Num = atoi( ad1name.substr( baseStr.size(), 2).c_str() );
    ad2Num = atoi( ad2name.substr( baseStr.size(), 2).c_str() );

    return ad1Num < ad2Num;

  }

  //----------------------------------------------------------------------------
  CRTGeoObjectSorter::CRTGeoObjectSorter(
      fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  CRTGeoObjectSorter::~CRTGeoObjectSorter() {}

  //----------------------------------------------------------------------------
  void CRTGeoObjectSorter::SortAuxDets(
      std::vector<geo::AuxDetGeo*> & adgeo) const {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetICARUS);
  }

  //----------------------------------------------------------------------------
  void CRTGeoObjectSorter::SortAuxDetSensitive(
      std::vector<geo::AuxDetSensitiveGeo*> & adsgeo) const {
    std::sort(adsgeo.begin(), adsgeo.end(), sortAuxDetSensitiveICARUS);
  }

}  // namespace geo

