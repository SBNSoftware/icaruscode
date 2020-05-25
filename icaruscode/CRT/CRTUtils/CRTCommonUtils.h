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
    class CRTCommonUtils;
 }
}

class icarus::crt::CRTCommonUtils {

 public:
    CRTCommonUtils();

    int            GetAuxDetTypeCode(size_t adid);
    char           GetAuxDetType(size_t adid);
    string         GetAuxDetRegion(size_t adid);
    int            AuxDetRegionNameToNum(string reg);
    string         GetRegionNameFromNum(int num);
    char           GetRegTypeFromRegName(string name);
    int            GetTypeCodeFromRegion(string name);
    pair<uint8_t,uint8_t> ADToMac(size_t adid);
    int            ADToChanGroup(size_t adid);
    int            NFeb(size_t adid);
    string         MacToRegion(uint8_t mac);
    char           MacToType(uint8_t mac);
    int            MacToTypeCode(uint8_t mac);
    int            ChannelToAuxDetSensitiveID(uint8_t mac, int chan);
    size_t         MacToAuxDetID(uint8_t mac, int chan);
    TLorentzVector AvgIDEPoint(sim::AuxDetIDE ide);
    double         LengthIDE(sim::AuxDetIDE ide);
    int            GetLayerID(sim::AuxDetSimChannel const& adsc);
    int            GetLayerID(const art::Ptr<sim::AuxDetSimChannel> adsc);
    int            GetMINOSLayerID(size_t adid);
    TVector3       ChanToLocalCoords(const uint8_t mac, const int chan);
    TVector3       ChanToWorldCoords(const uint8_t mac, const int chan);

 private:

    geo::GeometryCore const* fGeoService;
    map<size_t,vector<pair<uint8_t,int>>> fAuxDetIdToFeb;
    map<uint8_t,vector<size_t>> fFebToAuxDetId;
    map<size_t,char>            fAuxDetIdToType;
    map<size_t,string>          fAuxDetIdToRegion;
    map<string,size_t>          fNameToAuxDetId;
    map<size_t,int>             fAuxDetIdToChanGroup;

    void   FillFebMap();
    void   FillAuxDetMaps();
    string AuxDetNameToRegion(string name);

};//CRTCommonUtils
// }//namespace crt
//}//namespace icarus







#endif
