/**
 * @file   icaruscode/Geometry/details/AuxDetSorting.cxx
 * @brief  Functions for sorting ICARUS CRT modules (auxiliary detectors).
 * @author Chris Hilgenberg, Gianluca Petrillo (refactoring only)
 * @date   August 7, 2018
 * @see    icaruscode/Geometry/details/AuxDetSorting.h
 */

// library header
#include "icaruscode/Geometry/details/AuxDetSorting.h"

// LArSoft libraries
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

// C/C++ standard libraries
#include <string>
#include <algorithm> // std::sort()
#include <cstdlib> // std::atoi()


namespace {
  
  //--------------------------------------------------------------------------
  /// Define sort order for CRT modules in standard configuration.
  bool AuxDetStandardSortingRule
    (const geo::AuxDetGeo& ad1, const geo::AuxDetGeo& ad2)
  {
    
    std::string type1 = "", type2 = "";
    switch (ad1.NSensitiveVolume()) {
        case 20 : type1 = "MINOS"; break;
        case 16 : type1 = "CERN"; break;
        case 64 : type1 = "DC"; break;
    }
    switch (ad2.NSensitiveVolume()) {
        case 20 : type2 = "MINOS"; break;
        case 16 : type2 = "CERN"; break;
        case 64 : type2 = "DC"; break;
    }

    // sort based off of GDML name, module number
    std::string ad1name = (ad1.TotalVolume())->GetName();
    std::string ad2name = (ad2.TotalVolume())->GetName();
    // assume volume name is "volAuxDet_<type>_module_###_<region>"
    std::string base1 = "volAuxDet_"+type1+"_module_";
    std::string base2 = "volAuxDet_"+type2+"_module_";

    int ad1Num = std::atoi( ad1name.substr( base1.size(), 3).c_str() );
    int ad2Num = std::atoi( ad2name.substr( base2.size(), 3).c_str() );

    return ad1Num < ad2Num;

  } // AuxDetStandardSortingRule()
  
  
  //----------------------------------------------------------------------------
  /// Define sort order for CRT submodules in standard configuration.
  bool AuxDetSensitiveStandardSortingRule
    (const geo::AuxDetSensitiveGeo& ad1, const geo::AuxDetSensitiveGeo& ad2)
  {
    std::string type1 = "", type2 = "";

    // sort based off of GDML name, assuming ordering is encoded
    std::string ad1name = (ad1.TotalVolume())->GetName();
    std::string ad2name = (ad2.TotalVolume())->GetName();

    if ( ad1name.find("MINOS") != std::string::npos ) type1 = "MINOS";
    if ( ad1name.find("CERN") != std::string::npos ) type1 = "CERN";
    if ( ad1name.find("DC") != std::string::npos ) type1 = "DC";
    if ( ad2name.find("MINOS") != std::string::npos ) type2 = "MINOS";
    if ( ad2name.find("CERN") != std::string::npos ) type2 = "CERN";
    if ( ad2name.find("DC") != std::string::npos ) type2 = "DC";

    // assume volume name is "volAuxDetSensitive_<type>_module_###_strip_##"
    std::string baseMod1 = "volAuxDetSensitive_"+type1+"module_";
    std::string baseStr1 = "volAuxDetSensitive_"+type1+"module_###_strip_";
    std::string baseMod2 = "volAuxDetSensitive_"+type2+"module_";
    std::string baseStr2 = "volAuxDetSensitive_"+type2+"module_###_strip_";

    int ad1Num = atoi( ad1name.substr( baseMod1.size(), 3).c_str() );
    int ad2Num = atoi( ad2name.substr( baseMod2.size(), 3).c_str() );

    if(ad1Num!=ad2Num) return ad1Num < ad2Num;

    ad1Num = std::atoi( ad1name.substr( baseStr1.size(), 2).c_str() );
    ad2Num = std::atoi( ad2name.substr( baseStr2.size(), 2).c_str() );


    return ad1Num < ad2Num;

  } // AuxDetSensitiveStandardSortingRule()
  
  
  //----------------------------------------------------------------------------
  
} // local namespace


//------------------------------------------------------------------------------
void icarus::SortAuxDetsStandard(std::vector<geo::AuxDetGeo> & adgeo) {
  std::sort(adgeo.begin(), adgeo.end(), AuxDetStandardSortingRule);
}


//------------------------------------------------------------------------------
void icarus::SortAuxDetSensitiveStandard
  (std::vector<geo::AuxDetSensitiveGeo>& adsgeo)
{
  std::sort(adsgeo.begin(), adsgeo.end(), AuxDetSensitiveStandardSortingRule);
}


//------------------------------------------------------------------------------
