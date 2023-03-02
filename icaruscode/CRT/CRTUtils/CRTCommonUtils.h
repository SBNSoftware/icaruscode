#ifndef IC_CRTCOMMONUTILS_H
#define IC_CRTCOMMONUTILS_H

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"


// icaruscode includes
#include "sbnobj/Common/CRT/CRTHit.hh"

#include "TGeoManager.h"
//#include "Math/GenVector/XYZTVector.h"
//#include "Math/GenVector/LorentzVector.h" 
#include "TLorentzVector.h"
// ROOT
#include "TVector3.h"


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
    TVector3       WorldToModuleCoords(TVector3 point, size_t adid);
    // Simple distance of closest approach between infinite track and centre of hit
    double SimpleDCA(sbn::crt::CRTHit hit, TVector3 start, TVector3 direction);

    // Minimum distance from infinite track to CRT hit assuming that hit is a 2D square
    double DistToCrtHit(sbn::crt::CRTHit hit, TVector3 start, TVector3 end);

    // Distance between infinite line (2) and segment (1)
    // http://geomalgorithms.com/a07-_distance.html
    double LineSegmentDistance(TVector3 start1, TVector3 end1, TVector3 start2, TVector3 end2);

    // Intersection between axis-aligned cube and infinite line
    // (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
    std::pair<TVector3, TVector3> CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end);


 private:
    geo::AuxDetGeometryCore const* fAuxDetGeom;
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
