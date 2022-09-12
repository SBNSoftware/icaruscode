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
#include "TVector3.h"
#include "TGeoManager.h"


namespace icarus{


  struct  matchCand {
    sbn::crt::CRTHit thishit;
    double t0;
    double dca;
    double extrapLen;
  };


  class CRTT0MatchAlg {
  public:

    explicit CRTT0MatchAlg(const fhicl::ParameterSet& pset);
    CRTT0MatchAlg();
    
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

    std::pair<TVector3, TVector3> TrackDirectionAverageFromPoints(recob::Track track, double frac);

    // Keeping ClosestCRTHit function for backwards compatibility
    // *** use GetClosestCRTHit instead
    std::pair<sbn::crt::CRTHit, double> ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
						      recob::Track tpcTrack, std::pair<double, double> t0MinMax, 
						      std::vector<sbn::crt::CRTHit> crtHits, int driftDirection, uint64_t trigger_timestamp);

    std::vector<std::pair<sbn::crt::CRTHit, double> >ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
								   recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 
								   const art::Event& event, uint64_t trigger_timestamp);

    std::pair<sbn::crt::CRTHit, double> ClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
						      recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, 
						      std::vector<sbn::crt::CRTHit> crtHits, uint64_t trigger_timestamp);

    // Return the closest CRT hit to a TPC track and the DCA
    matchCand GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track tpcTrack, std::pair<double, double> t0MinMax, 
			       std::vector<sbn::crt::CRTHit> crtHits, int driftDirection, uint64_t& trigger_timestamp);

    std::vector<matchCand> GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
					    recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 
					    const art::Event& event, uint64_t trigger_timestamp);

    matchCand GetClosestCRTHit(detinfo::DetectorPropertiesData const& detProp,
			       recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, 
			       std::vector<sbn::crt::CRTHit> crtHits, uint64_t trigger_timestamp);

    // Match track to T0 from CRT hits
    std::vector<double> T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
				      recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 
				      const art::Event& event, uint64_t trigger_timestamp);

    double T0FromCRTHits(detinfo::DetectorPropertiesData const& detProp,
			 recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, 
			 std::vector<sbn::crt::CRTHit> crtHits, uint64_t& trigger_timestamp);
    
    // Match track to T0 from CRT hits, also return the DCA
    std::vector<std::pair<double, double> >  T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
								 recob::Track tpcTrack, std::vector<sbn::crt::CRTHit> crtHits, 
								 const art::Event& event, uint64_t trigger_timestamp);

    std::pair<double, double>  T0AndDCAFromCRTHits(detinfo::DetectorPropertiesData const& detProp,
						   recob::Track tpcTrack, std::vector<art::Ptr<recob::Hit>> hits, 
						   std::vector<sbn::crt::CRTHit> crtHits, uint64_t& trigger_timestamp);
 
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

  };


}
#endif
