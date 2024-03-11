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

namespace geo{

  //----------------------------------------------------------------------------
  // Define sort order for AuxDets in standard configuration
  static bool sortAuxDetICARUS(const AuxDetGeo& ad1, const AuxDetGeo& ad2) {

    // sort based off of GDML name, module number
    std::string ad1name = (ad1.TotalVolume())->GetName();
    std::string ad2name = (ad2.TotalVolume())->GetName();
    // assume volume name is "volAuxDet<subsystem>module###<region>"
    std::string modulePrefix = "module";

    //keep compatibility with legacy g4
    ad1name.erase(std::remove(ad1name.begin(), ad1name.end(), '_'), ad1name.end());
    ad2name.erase(std::remove(ad2name.begin(), ad2name.end(), '_'), ad2name.end());

    int ad1Num = atoi( ad1name.substr(ad1name.find(modulePrefix)+modulePrefix.length(),3).c_str() );
    int ad2Num = atoi( ad2name.substr(ad2name.find(modulePrefix)+modulePrefix.length(), 3).c_str() );

    return ad1Num < ad2Num;
  }

  //----------------------------------------------------------------------------
  // Define sort order for AuxDetSensitives in standard configuration
  static bool sortAuxDetSensitiveICARUS(const AuxDetSensitiveGeo& ad1,
                                      const AuxDetSensitiveGeo& ad2)
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string ad1name = (ad1.TotalVolume())->GetName();
    std::string ad2name = (ad2.TotalVolume())->GetName();
    // assume volume name is "volAuxDetSensitive<subsystem>module###(<cut<###>,top or bot>)strip##"
    std::string modulePrefix = "module";
    std::string stripPrefix = "strip";

    //keep compatibility with legacy g4
    ad1name.erase(std::remove(ad1name.begin(), ad1name.end(), '_'), ad1name.end());
    ad2name.erase(std::remove(ad2name.begin(), ad2name.end(), '_'), ad2name.end());

    int ad1Num = atoi( ad1name.substr(ad1name.find(modulePrefix)+modulePrefix.length(), 3).c_str() );
    int ad2Num = atoi( ad2name.substr(ad2name.find(modulePrefix)+modulePrefix.length(), 3).c_str() );

    if(ad1Num!=ad2Num) return ad1Num < ad2Num;

    ad1Num = atoi( ad1name.substr(ad1name.find(stripPrefix)+stripPrefix.length(), 2).c_str() );
    ad2Num = atoi( ad2name.substr(ad2name.find(stripPrefix)+stripPrefix.length(), 2).c_str() );

    return ad1Num < ad2Num;

  }

  //----------------------------------------------------------------------------
  // Define sort order for AuxDets in co-ordinate system

  static bool CRTIncreaseX(const AuxDetGeo& ad1, const AuxDetGeo& ad2) {
    auto const xyz1 = ad1.GetCenter();
    auto const xyz2 = ad2.GetCenter();
    return xyz1.X() < xyz2.X();
  }
  
  static bool CRTDecreaseY(const AuxDetGeo& ad1, const AuxDetGeo& ad2) {
    auto const xyz1 = ad1.GetCenter();
    auto const xyz2 = ad2.GetCenter();
    return xyz1.Y() > xyz2.Y();
  }
  
  static bool CRTIncreaseZ(const AuxDetGeo& ad1, const AuxDetGeo& ad2) {
    return ad1.GetCenter().Z() < ad2.GetCenter().Z();
  }
  

  //----------------------------------------------------------------------------
  // Define sort order for AuxDetSensitive in co-ordinate system

  static bool CRTSensitiveIncreaseX(std::pair<int , geo::AuxDetSensitiveGeo> const& ads1,
                                    std::pair<int , geo::AuxDetSensitiveGeo> const& ads2){
    return ads1.second.GetCenter().X() < ads2.second.GetCenter().X();
  }


  static bool CRTSensitiveDecreaseY(std::pair<int , geo::AuxDetSensitiveGeo> const& ads1,
                                    std::pair<int , geo::AuxDetSensitiveGeo> const& ads2){
    return ads1.second.GetCenter().Y() > ads2.second.GetCenter().Y();
  }

  static bool CRTSensitiveIncreaseZ(std::pair<int , geo::AuxDetSensitiveGeo> const& ads1,
                                    std::pair<int , geo::AuxDetSensitiveGeo> const& ads2){
    return ads1.second.GetCenter().Z() < ads2.second.GetCenter().Z();
  }


  //----------------------------------------------------------------------------
  CRTGeoObjectSorter::CRTGeoObjectSorter(
      fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  void CRTGeoObjectSorter::SortAuxDets(std::vector<geo::AuxDetGeo> & adgeo) const {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetICARUS);
  }

  //----------------------------------------------------------------------------
  void CRTGeoObjectSorter::SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo> & adsgeo) const {
    std::sort(adsgeo.begin(), adsgeo.end(), sortAuxDetSensitiveICARUS);
  }
  
  //----------------------------------------------------------------------------
  // sorting CRT co-ordiantes  decreasing vertical coordinate, 
  // next increasing beam coordinate, and next increasing drift direction

  void CRTGeoObjectSorter::SortCRTs(std::vector<geo::AuxDetGeo> & adgeo) const {

    // 3. stable sort by increasing drift direction (_x_)
    std::stable_sort (adgeo.begin(), adgeo.end(), CRTIncreaseX);


    // 2. stable sort by increasing beam coordinate (_z_)
    std::stable_sort (adgeo.begin(), adgeo.end(), CRTIncreaseZ);

    // 1. sort by decreasing vertical coordinate (_y_)
    std::stable_sort (adgeo.begin(), adgeo.end(), CRTDecreaseY);
  }


  //----------------------------------------------------------------------------
  // sorting CRTs Sensitive with co-ordiantes  decreasing vertical coordinate,
  // next increasing beam coordinate, and next increasing drift direction
  
  void CRTGeoObjectSorter::SortCRTSensitive(std::vector<std::pair <int ,geo::AuxDetSensitiveGeo> > & adsgeo) const {
    
    // 3. stable sort by increasing drift direction (_x_)
    std::stable_sort (adsgeo.begin(), adsgeo.end(), CRTSensitiveIncreaseX);
    

    // 2. stable sort by increasing beam coordinate (_z_)
    std::stable_sort (adsgeo.begin(), adsgeo.end(), CRTSensitiveIncreaseZ);
    
    // 1. sort by decreasing vertical coordinate (_y_)
    std::stable_sort (adsgeo.begin(), adsgeo.end(), CRTSensitiveDecreaseY);
  }
  

} //namespace crt
//} //namespace icarus
