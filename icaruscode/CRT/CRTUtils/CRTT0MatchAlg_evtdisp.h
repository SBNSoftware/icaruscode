#ifndef CRTT0MATCHALG_H_SEEN
#define CRTT0MATCHALG_H_SEEN


///////////////////////////////////////////////
// CRTT0MatchAlg_evtdisp.h
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
#include "lardataobj/RecoBase/PFParticle.h"
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
#include "TVector3.h"
#include "TGeoManager.h"


namespace icarus{


  struct  match_geometry {
    sbn::crt::CRTHit thishit;
    double t0;
    double dca;
    double extrapLen;


	//Evtdisp data vars
	bool simple_cathodecrosser;
	double t0min;
	double t0max;
	double crtTime;

	int hit_id;
	long track_id;
	long meta_track_id;

	TVector3 startDir;
	TVector3 endDir;
	TVector3 tpc_track_start;
	TVector3 tpc_track_end;
	TVector3 crt_hit_pos;
	double simpleDCA_startDir;
	double simpleDCA_endDir;
	bool is_best_DCA_startDir;
	bool is_best_DCA_endDir;

	//Evtdisp MC vars
//	TVector3 MC_tpc_trackstart;
//	TVector3 MC_tpc_trackend;
//	TVector3 MC_CRTHit_pt;
	

  };


  class CRTT0MatchAlg_evtdisp {
  public:

    explicit CRTT0MatchAlg_evtdisp(const fhicl::ParameterSet& pset);
    CRTT0MatchAlg_evtdisp();
    
    //void reconfigure(const Config& config);
     
    void reconfigure(const fhicl::ParameterSet& pset);

    // Utility function that determines the possible x range of a track
    std::pair<double, double> TrackT0Range(detinfo::DetectorPropertiesData const& detProp,
                                           double startX, double endX, int driftDirection, std::pair<double, double> xLimits);

    // Calculate the distance of closest approach (DCA) between the end of a track and a crt hit
    double DistOfClosestApproach(detinfo::DetectorPropertiesData const& detProp,
                                 TVector3 trackPos, TVector3 trackDir, sbn::crt::CRTHit crtHit, int driftDirection, double t0);

    std::pair<TVector3, TVector3> TrackDirectionAverage(recob::Track track, double frac);

    std::pair<TVector3, TVector3> TrackDirection(detinfo::DetectorPropertiesData const& detProp,recob::Track track, 
						 double frac, double CRTtime, int driftDirection);

//    std::pair<TVector3, TVector3> TrackDirectionAverageFromPoints(recob::Track track, double frac);

    // Return the closest CRT hit to a TPC track and the DCA
    std::vector<match_geometry> GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track tpcTrack, std::pair<double, double> t0MinMax, 
			       std::vector<sbn::crt::CRTHit> crtHits, int driftDirection, uint64_t& trigger_timestamp, const art::Event& event, bool IsData);

    std::vector<match_geometry> GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 
					    const art::Event& event, uint64_t trigger_timestamp, bool IsData);

    std::vector<match_geometry> GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, 
			       std::vector<sbn::crt::CRTHit> crtHits, uint64_t trigger_timestamp, const  art::Event& event, bool IsData);

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

    //Determine if the track crosses the cathode
    bool IsCathodeCrosser(const art::Event& event, recob::Track track, double &timestamp);

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
    bool   fIsData;
    //    double fDistEndpointAVedge;
    std::vector<art::InputTag> fTPCTrackLabel;

  };


}
#endif
