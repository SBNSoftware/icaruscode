#include "TPCGeoUtil.h"

#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"

namespace icarus {
  namespace TPCGeoUtil {
    int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> const& hits, const geo::GeometryCore& geometry) {
      // Return tpc of hit collection or -1 if in multiple
	//thanks to Jacob Z for his idea to change this function -tb

   if(hits.empty()) return -1;
//      geo::TPCID tpcID = hits[0]->WireID().asTPCID();
//      const geo::TPCGeo& tpcGeo = geometry.GetElement(tpcID);
   auto const driftDirection_init = geometry.TPC(hits[0]->WireID()).DriftAxisWithSign();
      int tpcid = (int)hits[0]->WireID().TPC;
//     for(size_t i = 0; i < hits.size(); i++){
//       tpcID = hits[i]->WireID().asTPCID();
//       const geo::TPCGeo& tpcGeo_now = geometry.GetElement(tpcID);
   for(art::Ptr<recob::Hit> const& hit: hits){
       auto const driftDirection_now = geometry.TPC(hit->WireID()).DriftAxisWithSign();
       if(driftDirection_now != driftDirection_init) {return -1;}
     }
   return tpcid;  

    }//end definition of int DetectedInTPC

    // Work out the drift limits for a collection of hits
    std::pair<double, double> XLimitsFromHits(geo::GeometryCore const& geometry,
                                              geo::WireReadoutGeom const& wireReadout,
                                              std::vector<art::Ptr<recob::Hit>> const& hits){
      // If there are no hits then return 0
      if(hits.empty()) return std::make_pair(0, 0);
  
      // If the track is stitched (in multiple TPCs) return 0
      if(DetectedInTPC(hits, geometry) == -1) return std::make_pair(0, 0);

      // Work out the drift direction
//      geo::TPCID tpcID = hits[0]->WireID().asTPCID();
//      const geo::TPCGeo& tpcGeo = geometry.GetElement(tpcID);
        const geo::TPCGeo& tpcGeo = geometry.TPC(hits[0]->WireID());
      return std::make_pair(tpcGeo.MinX(), tpcGeo.MaxX());
    }

    int DriftDirectionFromHits(const geo::GeometryCore& geometry, std::vector<art::Ptr<recob::Hit>> const& hits){
      // If there are no hits then return 0
      if(hits.empty()) return 0;
  
      // If the track is stitched (in multiple TPCs) return 0
      if(DetectedInTPC(hits, geometry) == -1) return 0;

      // Work out the drift direction
//      geo::TPCID tpcID = hits[0]->WireID().asTPCID();
//      const geo::TPCGeo& tpcGeo = geometry.GetElement(tpcID);
      auto const [axis, sign] = geometry.TPC(hits[0]->WireID()).DriftAxisWithSign();
      if(axis == geo::Coordinate::X) return 0;
      return to_int(sign);
    }

    // Is point inside given TPC
    bool InsideTPC(geo::Point_t point, const geo::TPCGeo& tpc, double buffer){
      if(point.X() < (tpc.MinX()-buffer) || point.X() > (tpc.MaxX()+buffer)
	 || point.Y() < (tpc.MinY()-buffer) || point.Y() > (tpc.MaxY()+buffer)
	 || point.Z() < (tpc.MinZ()-buffer) || point.Z() > (tpc.MaxZ()+buffer)) return false;
      return true;
    }

  } // namespace TPCGeoUtil
} // namespace icarus
