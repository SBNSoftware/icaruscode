#ifndef ICARUSCODE_CRT_CRTUTILS_CRTMATCHINGUTILS_H
#define ICARUSCODE_CRT_CRTUTILS_CRTMATCHINGUTILS_H

/**
 * @file   icaruscode/CRT/CRTUtils/CRTMatchingUtils.h
 * @author Francesco Poppi (poppi@bo.infn.it)
 * @date   January 2025
 */

#include <vector>
#include <utility> // std::pair
#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h" 
// LArSoft headers
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
// Data product headers
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"

namespace icarus::crt{

struct TrackBarycenter
{
    double BarX; // Track Barycenter X coordinate
    double BarY; // Track Barycenter Y coordinate
    double BarZ; // Track Barycenter Z coordinate
    bool isGood; // Track Barycenter quality
};

struct DriftedTrack
{
    std::vector<geo::Point_t> sp; // Drifted space points
    std::vector<double> spi; // Drifted Track Hit Points integral
    int outbound; // Number of hit points out of the logical volume of the TPC
    //double drifted_startx; 
    //double drifted_endx;
};

struct TranslationVector{
    geo::Vector_t dir;
    geo::Point_t mean;
};

struct PCAResults{
    geo::Vector_t eigenVector1; // First EigenVector
    geo::Vector_t eigenVector2; // Second EigenVector
    geo::Vector_t eigenVector3; // Third EigenVector
    double eigenValue1; // First EigenValue
    double eigenValue2; // Second EigenValue
    double eigenValue3; // Third EigenValue
    geo::Point_t mean; // Mean X,Y,Z coordinates
};

using CrossingPoint = geo::Point_t;

struct CandCRT{
    sbn::crt::CRTHit CRThit;
    art::Ptr<sbn::crt::CRTHit> ptrCRThit;
    int plane;
    double distance;
    double deltaX;
    double deltaY;
    double deltaZ;
    CrossingPoint crossPoint;
};

/// CRTPlane corresponds to a pair of:
/// 1st CRT fixed coordinate plane (e.g. 0 is X coordinate)
/// 2nd value of the  CRT fixed coordinate plane (e.g. for Top CRT Horizontal this value is ~618 cm.
using CRTPlane = std::pair<int,double>;

class CRTMatchingAlg {
public:

    explicit CRTMatchingAlg(const fhicl::ParameterSet& pset);
    CRTMatchingAlg();

    void reconfigure(const fhicl::ParameterSet& pset);
    
    /// This function runs the PCA analysis on the spatial coordinates vector and return the PCAResults.
    static PCAResults PCAfit (std::vector<geo::Point_t> const& sp);
    
    /// This function determines the coordinate in which the CRT module is constant.
    /// 0 for fixed Y, 1 for fixed X, 2 for fixed Z,
    /// e.g. in the Top CRT Horizontal Plane, the Y coordinate is fixed, in the Side CRT West Walll, the X Coordinate is fixed, ...
    CRTPlane DeterminePlane(sbn::crt::CRTHit const& CRThit);

    /// This function evaluates the Track Crossing point onto the CRT plane considered.
    /// dir is the track direction in the CRTWall reference system.
    /// mean is the mean value of the track in the CRTWall reference system.
    static CrossingPoint TranslatePointTo(geo::Vector_t dir, geo::Point_t mean, CRTPlane CRTWall);

    /// This function rotate the translation vector (direction and mean position) from the TPC reference system, to the CRTWall one.
    TranslationVector RotateToLocalCRTPlane(const TranslationVector& transl, CRTPlane CRTWall);

    /// This function rotate the CrossingPoint from the CRTWall reference system, to the TPC one.
    CrossingPoint RotateFromLocalCRTPlane(CrossingPoint crossPointCRT, CRTPlane CRTWall);

    /// This function returns the "correct" predicted crossing point given the cosine directors of the track fitted with PCA
    /// and a CRTPlane (fixed coordinate and value of the coordinate).
    CrossingPoint DetermineProjection(const TranslationVector& dir, CRTPlane CRTWall);
    
    /// This function evaluates the TrackBarycenter with weighted mean. posVector is a vector of track spacepoints
    /// and w is the weight (e.g. one could use integral or peak amplitude or more...).
    /// The entries in the i-th index of posVector and hw vectors must correspond to the same point.
    TrackBarycenter GetTrackBarycenter(std::vector<geo::Point_t> posVector, std::vector<double> hw);

    /// Function which drifts the track x coordinates assuming a time (in the trigger reference system).
    /// The function also returns the number of hit points which would be outside of the physical boundaries of the TPCs at the considered time.
    DriftedTrack DriftTrack(const std::vector<art::Ptr<recob::Hit>>& trkHits, const std::vector<const recob::TrackHitMeta*>& trkHitMetas, const geo::GeometryCore *GeometryService, detinfo::DetectorPropertiesData const& detProp, double time, const recob::Track& tpcTrack) const;

private:

    double fTickPeriod;
    double fTickAtAnode;
    double fAllowedOffsetCM;
};

}

#endif // CRTMATCHINGUTILS_H