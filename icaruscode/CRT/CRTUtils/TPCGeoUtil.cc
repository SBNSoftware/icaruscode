#include "TPCGeoUtil.h"

namespace icarus {
  namespace TPCGeoUtil {
    int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits, const geo::GeometryCore *GeometryService){
      // Return tpc of hit collection or -1 if in multiple
	//thanks to Jacob Z for his idea to change this function -tb

   if(hits.size() == 0) return -1;
      geo::TPCID tpcID = hits[0]->WireID().asTPCID();
      const geo::TPCGeo& tpcGeo = GeometryService->GetElement(tpcID);
      double driftDirection_init = tpcGeo.DetectDriftDirection();
      int tpcid = (int)hits[0]->WireID().TPC;
     for(size_t i = 0; i < hits.size(); i++){
       tpcID = hits[i]->WireID().asTPCID();
       const geo::TPCGeo& tpcGeo_now = GeometryService->GetElement(tpcID);
       double driftDirection_now = tpcGeo_now.DetectDriftDirection();
       if(driftDirection_now != driftDirection_init) {return -1;}
     }
   return tpcid;  

    }//end definition of int DetectedInTPC

    // Work out the drift limits for a collection of hits
    std::pair<double, double> XLimitsFromHits(const geo::GeometryCore *GeometryService, std::vector<art::Ptr<recob::Hit>> hits){
      // If there are no hits then return 0
      if(hits.size() == 0) return std::make_pair(0, 0);
  
      // If the track is stitched (in multiple TPCs) return 0
      if(DetectedInTPC(hits, GeometryService) == -1) return std::make_pair(0, 0);

      // Work out the drift direction
      geo::TPCID tpcID = hits[0]->WireID().asTPCID();
      const geo::TPCGeo& tpcGeo = GeometryService->GetElement(tpcID);
      return std::make_pair(tpcGeo.MinX(), tpcGeo.MaxX());
    }

    int DriftDirectionFromHits(const geo::GeometryCore *GeometryService, std::vector<art::Ptr<recob::Hit>> hits){
      // If there are no hits then return 0
      if(hits.size() == 0) return 0;
  
      // If the track is stitched (in multiple TPCs) return 0
      if(DetectedInTPC(hits, GeometryService) == -1) return 0;

      // Work out the drift direction
      geo::TPCID tpcID = hits[0]->WireID().asTPCID();
      const geo::TPCGeo& tpcGeo = GeometryService->GetElement(tpcID);
      double driftDirection = tpcGeo.DetectDriftDirection();
      if(std::abs(driftDirection) != 1) driftDirection = 0;
      return driftDirection;
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

