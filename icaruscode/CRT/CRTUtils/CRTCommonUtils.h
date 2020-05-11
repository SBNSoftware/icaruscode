#ifndef IC_CRTCOMMONUTILS_H
#define IC_CRTCOMMONUTILS_H

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"

#include "TGeoManager.h"
//#include "Math/GenVector/XYZTVector.h"
//#include "Math/GenVector/LorentzVector.h" 
#include "TLorentzVector.h"

#include <map>
#include <vector>
#include <string>
#include <utility>

using std::string;
using std::map;
using std::vector;
using std::pair;

namespace icarus{
 namespace crt {
  namespace CRTCommonUtils {

    int            ModToTypeCode(geo::AuxDetGeo const& adgeo);
    char           GetAuxDetType(geo::AuxDetGeo const& adgeo);
    string         GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
    int            GetAuxDetRegionNum(std::string reg);
    string         GetRegionNameFromNum(int num);
    pair<uint32_t,uint32_t> ADToMac(const map<int,vector<pair<int,int>>>& febMap, uint32_t adid);
    int            MacToRegion(int mac);
    char           MacToType(int mac);
    string         MacToRegionName(int mac);
    int            ChannelToAuxDetSensitiveID(int mac, int chan);
    int            MacToAuxDetID(int mac, int chan);
    map<int,vector<pair<int,int>>> GetFebMap();
    //ROOT::Math::XYZTVector 
    TLorentzVector AvgIDEPoint(sim::AuxDetIDE ide);
    double         LengthIDE(sim::AuxDetIDE ide);
    int            GetLayerID(geo::GeometryCore const* geoService, sim::AuxDetSimChannel const& adsc);
    int            GetLayerID(geo::GeometryCore const* geoService, const art::Ptr<sim::AuxDetSimChannel> adsc);

  }//CRTCommonUtils
 }//namespace crt
}//namespace icarus







#endif
