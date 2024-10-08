#include "icaruscode/CRT/CRTGeoObjectSorter.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "icaruscode/CRT/compareCRTs.h"

#include "art/Utilities/ToolMacros.h"

namespace geo{

  //----------------------------------------------------------------------------
  CRTGeoObjectSorter::CRTGeoObjectSorter() = default;

  CRTGeoObjectSorter::CRTGeoObjectSorter(
      fhicl::ParameterSet const&) {}

  bool CRTGeoObjectSorter::compareAuxDets(AuxDetGeo const& ad1, AuxDetGeo const& ad2) const
  {
    // sort based off of GDML name, module number
    std::string ad1name = (ad1.TotalVolume())->GetName();
    std::string ad2name = (ad2.TotalVolume())->GetName();
    // assume volume name is "volAuxDet_<subsystem>_module_###_<region>"
    std::string modulePrefix = "module_";

    int ad1Num = atoi( ad1name.substr(ad1name.find(modulePrefix)+modulePrefix.length(),3).c_str() );
    int ad2Num = atoi( ad2name.substr(ad2name.find(modulePrefix)+modulePrefix.length(), 3).c_str() );

    return ad1Num < ad2Num;
  }

  bool CRTGeoObjectSorter::compareAuxDetSensitives(AuxDetSensitiveGeo const& ads1,
                                                   AuxDetSensitiveGeo const& ads2) const
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string ad1name = (ads1.TotalVolume())->GetName();
    std::string ad2name = (ads2.TotalVolume())->GetName();
    // assume volume name is "volAuxDetSensitive_<subsystem>_module_###_(<cut<###>,top or bot>_)strip_##"
    std::string modulePrefix = "module_";
    std::string stripPrefix = "strip_";

    int ad1Num = atoi( ad1name.substr(ad1name.find(modulePrefix)+modulePrefix.length(), 3).c_str() );
    int ad2Num = atoi( ad2name.substr(ad2name.find(modulePrefix)+modulePrefix.length(), 3).c_str() );

    if(ad1Num!=ad2Num) return ad1Num < ad2Num;

    ad1Num = atoi( ad1name.substr(ad1name.find(stripPrefix)+stripPrefix.length(), 2).c_str() );
    ad2Num = atoi( ad2name.substr(ad2name.find(stripPrefix)+stripPrefix.length(), 2).c_str() );

    return ad1Num < ad2Num;
  }

} //namespace geo

DEFINE_ART_CLASS_TOOL(geo::CRTGeoObjectSorter)
