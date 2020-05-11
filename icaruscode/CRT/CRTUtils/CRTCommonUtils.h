#ifndef IC_CRTCOMMONUTILS_H
#define IC_CRTCOMMONUTILS_H

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"

#include <map>
#include <vector>
#include <string>

namespace icarus{
 namespace crt {
  namespace CRTCommonUtils {

    char GetAuxDetType(geo::AuxDetGeo const& adgeo);
    std::string GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
    int GetAuxDetRegionNum(std::string reg);
    std::string GetRegionNameFromNum(int num);
    int MacToRegion(int mac);
    char MacToType(int mac);
    std::string MacToRegionName(int mac);
    int ChannelToAuxDetSensitiveID(int mac, int chan);
    int MacToAuxDetID(int mac, int chan);
    std::map<int,std::vector<std::pair<int,int>>> GetFebMap();

  }//CRTCommonUtils
 }//namespace crt
}//namespace icarus







#endif
