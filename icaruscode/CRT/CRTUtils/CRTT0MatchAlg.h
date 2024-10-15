#ifndef CRTT0MATCHALG_H_SEEN
#define CRTT0MATCHALG_H_SEEN


///////////////////////////////////////////////
// CRTT0MatchAlg.h
//
// Functions for CRT t0 matching
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
//#include "icaruscode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "icaruscode/CRT/CRTUtils/TPCGeoUtil.h"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

// ROOT
#include <vector>
#include "TVector3.h"
#include "TGeoManager.h"


namespace icarus{


  struct  matchCand {
    sbn::crt::CRTHit thishit;
    double t0 = DBL_MIN;
    double dca = DBL_MIN;
    double extrapLen = DBL_MIN;
    int best_DCA_pos = -1;//-1=no match; 0=startdir is best; 1=enddir is best; 2=both start+end are equally good;
    bool simple_cathodecrosser = false;
    int driftdir = -5;
    double t0min = DBL_MIN;
    double t0max = DBL_MIN;
    double crtTime = DBL_MIN;
    TVector3 startDir{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 endDir{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 tpc_track_start{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 tpc_track_end{DBL_MIN,DBL_MIN,DBL_MIN};

  };

  struct  match_geometry {
    sbn::crt::CRTHit thishit;
    double t0 = DBL_MIN;
    double dca = DBL_MIN;
    double extrapLen = DBL_MIN;
    bool simple_cathodecrosser = false;
    int best_DCA_pos = -1;
    int driftDir = -2;

    //Evtdisp data vars
    double t0min = DBL_MIN;
    double t0max = DBL_MIN;
    double crtTime = DBL_MIN;
	
    int driftdir = -5;	
    int hit_id = INT_MIN;
    long track_id = LONG_MIN;
    long meta_track_id = LONG_MIN;

    TVector3 startDir{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 endDir{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 tpc_track_start{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 tpc_track_end{DBL_MIN,DBL_MIN,DBL_MIN};
    TVector3 crt_hit_pos{DBL_MIN,DBL_MIN,DBL_MIN};
    double simpleDCA_startDir = DBL_MIN;
    double simpleDCA_endDir = DBL_MIN;
  };



  class CRTT0MatchAlg {
  public:

    explicit CRTT0MatchAlg(const fhicl::ParameterSet& pset);
    CRTT0MatchAlg();
    
    void reconfigure(const fhicl::ParameterSet& pset);

    // Utility function that determines the possible x range of a track
    std::pair<double, double> TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
                                           double startX, double endX, int driftDirection, std::pair<double, double> xLimits) const;

    // Calculate the distance of closest approach (DCA) between the end of a track and a crt hit
    double DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
                                 geo::Point_t const& trackPos, TVector3 trackDir, sbn::crt::CRTHit const& crtHit, int driftDirection, double t0) const;

    std::pair<TVector3, TVector3> TrackDirectionAverage(recob::Track const& track, double frac) const;

    std::pair<TVector3, TVector3> TrackDirection(detinfo::DetectorPropertiesData const& detProp,recob::Track const& track, 
						 double frac, double CRTtime, int driftDirection) const;

    std::pair<TVector3, TVector3> TrackDirectionAverageFromPoints(recob::Track const& track, double frac) const;

    // Keeping ClosestCRTHit function for backwards compatibility
    // *** use GetClosestCRTHit instead
    std::pair<sbn::crt::CRTHit, double> ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
						      recob::Track const& tpcTrack, std::pair<double, double> t0MinMax, 
						      std::vector<sbn::crt::CRTHit> const& crtHits, int driftDirection, uint64_t trigger_timestamp) const;

    std::vector<std::pair<sbn::crt::CRTHit, double> >ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
								   recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
								   const art::Event& event, uint64_t trigger_timestamp) const;

    std::pair<sbn::crt::CRTHit, double> ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
						      recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
						      std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp) const;

    // Return the closest CRT hit to a TPC track and the DCA
    matchCand GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
			       std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp, bool IsData) const;

    std::vector<matchCand> GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
					    const art::Event& event, uint64_t trigger_timestamp) const;

    matchCand GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track const& tpcTrack, std::pair<double, double> t0MinMax, 
			       std::vector<sbn::crt::CRTHit> const& crtHits, int driftDirection, uint64_t& trigger_timestamp, bool IsData) const;

    // Match track to T0 from CRT hits
    std::vector<double> T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
				      recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
				      const art::Event& event, uint64_t trigger_timestamp) const;

    double T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
			 recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
			 std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t& trigger_timestamp) const;
    
    // Match track to T0 from CRT hits, also return the DCA
    std::vector<std::pair<double, double> >  T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
								 recob::Track const& tpcTrack, std::vector<sbn::crt::CRTHit> const& crtHits, 
								 const art::Event& event, uint64_t trigger_timestamp) const;

    std::pair<double, double>  T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
						   recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
						   std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t& trigger_timestamp) const;
 
    // Simple distance of closest approach between infinite track and centre of hit
    double SimpleDCA(sbn::crt::CRTHit const& hit, TVector3 start, TVector3 direction) const;

    // Minimum distance from infinite track to CRT hit assuming that hit is a 2D square
    double DistToCrtHit(sbn::crt::CRTHit const& hit, TVector3 start, TVector3 end) const;

    // Distance between infinite line (2) and segment (1)
    // http://geomalgorithms.com/a07-_distance.html
    double LineSegmentDistance(TVector3 start1, TVector3 end1, TVector3 start2, TVector3 end2) const;

    // Intersection between axis-aligned cube and infinite line
    // (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
    std::pair<TVector3, TVector3> CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end) const;

    //Get multiple CRT Hits potentially related to a TPC track based on DCA
    std::vector<icarus::match_geometry> GetClosestCRTHit_geo(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track const& tpcTrack, std::vector<art::Ptr<recob::Hit>> const& hits, 
			       std::vector<sbn::crt::CRTHit> const& crtHits, uint64_t trigger_timestamp, bool IsData) const;

 /*   std::vector<icarus::match_geometry> GetClosestCRTHit_geo(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 
					    const art::Event& event, uint64_t trigger_timestamp, bool IsData);
*/
    std::vector<icarus::match_geometry> GetClosestCRTHit_geo(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track const& tpcTrack, std::pair<double, double> t0MinMax, 
			       std::vector<sbn::crt::CRTHit> const& crtHits, int driftDirection, uint64_t& trigger_timestamp, bool IsData) const;

    double GetCRTTime(sbn::crt::CRTHit const& crthit, uint64_t trigger_timestamp, bool isdata) const;

  private:

    geo::GeometryCore const* fGeometryService;
    spacecharge::SpaceCharge  const* fSCE;

    double fMinTrackLength;
    double fTrackDirectionFrac;
    double fDistanceLimit;
    int    fTSMode;
    double fTimeCorrection;
    int    fDirMethod;
    bool   fSCEposCorr;
    bool   fDCAuseBox;
    bool   fDCAoverLength;
    double fDoverLLimit;
    double fPEcut;
    double fMaxUncert;
    //    double fDistEndpointAVedge;
    std::vector<art::InputTag> fTPCTrackLabel;
//    bool IsData;
  };


}
#endif
