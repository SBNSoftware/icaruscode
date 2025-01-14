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

struct CrossPoint
{
    double X;
    double Y;
    double Z;
};

using ProjectionPoint = CrossPoint;

struct TrackBarycenter
{
    double BarX; // Track Barycenter X coordinate
    double BarY; // Track Barycenter Y coordinate
    double BarZ; // Track Barycenter Z coordinate
    bool isGood; // Track Barycenter quality
};

struct DriftedTrack
{
    std::vector<double> spx; // Drifted Track Hit Points X coordinate
    std::vector<double> spy; // Drifted Track Hit Points Y coordinate
    std::vector<double> spz; // Drifted Track Hit Points Z coordinate
    std::vector<double> spi; // Drifted Track Hit Points integral
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
    int plane;
    double distance;
    double deltaX;
    double deltaY;
    double deltaZ;
    double crossX;
    double crossY;
    double crossZ;
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

using FebIndex_t = int;

using TopCRTCorrectionMap = std::map<FebIndex_t, AffineTrans>;

struct TopCRTTransformations
{
    TopCRTCorrectionMap EE;
    TopCRTCorrectionMap EW;
    TopCRTCorrectionMap EastCC;
    TopCRTCorrectionMap WE;
    TopCRTCorrectionMap WW;
    TopCRTCorrectionMap WestCC;
};

using TopCRTCentersMap = std::map<FebIndex_t, ModuleCenter>;

/// CRTPlane corresponds to a pair of:
/// 1st CRT fixed coordinate plane (e.g. 0 is X coordinate)
/// 2nd value of the  CRT fixed coordinate plane (e.g. for Top CRT Horizontal this value is ~618 cm.
using CRTPlane = std::pair<int,double>;

/// The transformed CRT Hits are in cm
using TransformedCRTHit = std::pair<double, double>;

/// This function loads the Top CRT modules centers.
TopCRTCentersMap LoadTopCRTCenters();

/// This function performs the affine transformation of the CRT hit points.
/// The AffineTransformation requires input variables in cm.
TransformedCRTHit AffineTransformation(double DX, double DZ, AffineTrans affine);

/// This functions loads the Affine Transformation TXT files.
TopCRTTransformations LoadTopCRTTransformations();


class CRTMatchingAlg {
public:

    explicit CRTMatchingAlg(const fhicl::ParameterSet& pset);
    CRTMatchingAlg();

    void reconfigure(const fhicl::ParameterSet& pset);
    
    /// This function runs the PCA analysis on three vectors of spatial coordinates x,y,z and return the principal eigenvector.
    /// The entries in the i-th index of X, Y and Z vectors must correspond to the same point.
    static Direction PCAfit (std::vector<double> const& x, std::vector<double> const& y, std::vector<double> const& z);
    
    /// This function determines the coordinate in which the CRT module is constant.
    /// e.g. in the Top CRT Horizontal Plane, the Y coordinate is fixed, in the Side CRT West Walll, the X Coordinate is fixed, ...
    CRTPlane DeterminePlane(sbn::crt::CRTHit const& CRThit);

    /// This function evaluates the Track Crossing point onto the CRT plane considered.
    /// dir1, dir2, dir3 are the three cosine directors of a track.
    /// p1, p2 and p3 are the coordinates of the CRT hit position.
    static ProjectionPoint TranslatePointTo(double dir1, double dir2, double dir3, double p1, double p2, double p3, double position);
    //static ProjectionPoint CalculateProjection(double dir1, double dir2, double dir3, double p1, double p2, double p3, double position);

    /// This function runs the CalculateProjection function, after deciding the director cosine order expected from
    /// the CalculateProjection function, assuming the evaluated CRTPlane (fixed coordinate and value of the coordinate).
    CrossPoint CalculateForPlane(const Direction& dir, int plane, double position);

    /// This function returns the "correct" predicted crossing point given the cosine directors of the track fitted with PCA
    /// and a CRTPlane (fixed coordinate and value of the coordinate).
    CrossPoint DetermineProjection(const Direction& dir, CRTPlane plane);

    /// This function evaluates the TrackBarycenter with weighted mean. hx, hy, hz and hw are vectors with the x, y z coordinates
    /// and w is the weight (e.g. one could use integral or peak amplitude or more...).
    /// The entries in the i-th index of hx, hy, hz and hw vectors must correspond to the same point.
    TrackBarycenter GetTrackBarycenter(std::vector<double> hx, std::vector<double> hy, std::vector<double> hz, std::vector<double> hw);

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