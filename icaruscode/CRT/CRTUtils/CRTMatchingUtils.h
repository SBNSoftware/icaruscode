#ifndef CRTMATCHINGUTILS_H
#define CRTMATCHINGUTILS_H

///////////////////////////////////////////////
// CRTMatchingUtils.h
//
// Functions for CRT matching
// Francesco Poppi (poppi@bo.infn.it), October 2024
///////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <sstream>
#include <fstream> 
#include <string>
#include <cetlib/search_path.h>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

#include "canvas/Persistency/Common/Ptr.h" 
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TGraph2D.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Data product includes
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"

// ROOT
#include <vector>
#include "TVector3.h"
#include "TGeoManager.h"

namespace icarus::crt{

double const vdrift=0.157; //cm/us
double const tics=0.4; // 0.4 us, 400 ns
double const minLimitW=61.7; // cm Anode WE position
double const maxLimitW=358.73; // cm Anode WW position
double const minLimitE=-61.7; // cm Anode EW position
double const maxLimitE=-358.73; // cm Anode WW position
double const cathW=210; // cm Cathode W position
double const cathE=-210; // cm Cathode E position
double const exc=2; // cm max displacement out of boundaries

struct CrossPoint
{
    double X;
    double Y;
    double Z;
};

typedef CrossPoint ProjectionPoint;

struct TrackBarycenter
{
    float BarX; // Track Barycenter X coordinate
    float BarY; // Track Barycenter Y coordinate
    float BarZ; // Track Barycenter Z coordinate
    bool isGood; // Track Barycenter quality
};

struct DriftedTrack
{
    std::vector<float> spx; // Drifted Track Hit Points X coordinate
    std::vector<float> spy; // Drifted Track Hit Points Y coordinate
    std::vector<float> spz; // Drifted Track Hit Points Z coordinate
    std::vector<float> spi; // Drifted Track Hit Points integral
    int outbound; // Number of hit points out of the logical volume of the TPC
    //double drifted_startx; // 
    //double drifted_endx;
};

struct Direction
{
    double dirx; // Direction of the Track: X 
    double diry; // Direction of the Track: Y 
    double dirz; // Direction of the Track: Z 
    double meanx; // Mean Point of the Track: X 
    double meany; // Mean Point of the Track: Y 
    double meanz; // Mean Point of the Track: Z
};

struct CandCRT{
    sbn::crt::CRTHit CRThit;
    art::Ptr<sbn::crt::CRTHit> ptrCRThit;
    double distance;
    double deltaX;
    double deltaZ;
};

struct ModuleCenter
{
    double X;
    double Y;
    double Z;
};

struct AffineTrans
{
    double A11;
    double A12;
    double A21;
    double A22;
    double B1;
    double B2;
    double Accuracy;
    double N;
};

typedef int feb_index;

typedef std::map<feb_index, AffineTrans> TopCRTCorrectionMap;

struct TopCrtTransformations
{
    TopCRTCorrectionMap EE;
    TopCRTCorrectionMap EW;
    TopCRTCorrectionMap EastCC;
    TopCRTCorrectionMap WE;
    TopCRTCorrectionMap WW;
    TopCRTCorrectionMap WestCC;
};

typedef std::map<feb_index, ModuleCenter> TopCRTCentersMap;

typedef std::pair<int, double> CrtPlane;

typedef std::pair<double, double> TransformedCrtHit;

TopCRTCentersMap LoadTopCRTCenters();

TransformedCrtHit AffineTransformation(double DX, double DZ,AffineTrans affine);

TopCrtTransformations LoadTopCrtTransformations();


class CRTMatchingAlg {
public:

    explicit CRTMatchingAlg(const fhicl::ParameterSet& pset);
    CRTMatchingAlg();

    void reconfigure(const fhicl::ParameterSet& pset);

    Direction PCAfit (std::vector<float> x, std::vector<float> y, std::vector<float> z);

    CrtPlane DeterminePlane(sbn::crt::CRTHit CRThit);

    ProjectionPoint CalculateProjection(double, double, double, double, double, double, double);

    CrossPoint CalculateForPlane(const Direction& dir, int plane, double position);

    CrossPoint DetermineProjection(const Direction& dir, CrtPlane plane);

    TrackBarycenter GetTrackBarycenter(std::vector<float> hx, std::vector<float> hy, std::vector<float> hz, std::vector<float> hw);

    DriftedTrack DriftTrack(const std::vector<art::Ptr<recob::Hit>>& trkHits, const std::vector<const recob::TrackHitMeta*>& trkHitMetas, const geo::GeometryCore *GeometryService, detinfo::DetectorPropertiesData const& detProp, double time, int cryo, const recob::Track& tpcTrack);

private:

};

}

#endif // CRTMATCHINGUTILS_H