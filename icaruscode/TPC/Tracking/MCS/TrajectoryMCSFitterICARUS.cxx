#include "TrajectoryMCSFitterICARUS.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardata/RecoBaseProxy/Track.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"

#include <numeric>

using namespace std;
using namespace trkf;
using namespace recob::tracking;

//main function, return MCS momentum
recob::MCSFitResultGS TrajectoryMCSFitterICARUS::fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst) const {
  //print some fhicl parameters
  cout << "cut mode = " << cutMode_ << endl;
  cout << "dim mode = " << dimMode_ << endl;
  cout << "plane mode = " << planeMode_ << endl;
  cout << "fit mode = " << fitMode_ << endl;

  //geometrical check if track is stopping
  bool stop = GeoStopCheck(traj);
  
  //check if track length is less than a certain value [cm]
  auto lenmin = 100;
  auto lenrec = traj.Length();
  auto lenval = sqrt((traj.LocationAtPoint(traj.LastValidPoint()) - traj.LocationAtPoint(traj.FirstValidPoint())).Mag2());
  if (lenrec < lenmin || lenval < lenmin) {
    cout << "length less than " << lenmin << " cm, stopping fit" << endl;
    return recob::MCSFitResultGS(); }
  
  //for 2d, print 3d and 2d hits of track in chosen plane and last tpc
  cout << "start = " << traj.LocationAtPoint(traj.FirstValidPoint()) << endl;
  cout << "end = " << traj.LocationAtPoint(traj.LastValidPoint()) << endl;
  size_t index = traj.FirstValidPoint();
  while (index < traj.LastValidPoint()) {
    if (isinplane(index, planeMode_) && isintpc(index, lasttpc(traj))) {
      cout << "3d hit = " << traj.LocationAtPoint(index) << endl;
      cout << "2d hit = " << hit2d(traj, index, planeMode_, lasttpc(traj)) << endl; 
      cout << "integral = " << Integral(traj, index) << endl;
      cout << "index = " << index << endl; }
    index = traj.NextValidPoint(index + 1); }

  //define vector breakpoints with segmentation index
  vector<size_t> breakpoints; 
  //define vector seglens and cumseglens with segment lenghts [cm] and cumulative segment lengths [cm]
  vector<float> seglens; vector<float> cumseglens; 
  //define vector seghits and cumseghits with segment number of hits and cumulative segment number of hits
  vector<int> seghits; vector<int> cumseghits;

  vector<int> isDelta; 
  std::cout << " before processing delta rays " << planeMode_ <<  std::endl;
  for (unsigned int jp = 0; jp < traj.NPoints(); jp++) isDelta.push_back(0);
  if (planeMode_ < 3) {
    ProcessDeltaRays(traj, planeMode_, isDelta);
    std::cout << " after processing delta rays " << planeMode_ <<  std::endl;
    }
  for (unsigned int jp = 0; jp < traj.NPoints(); jp++)
   if (isDelta[jp]) std::cout << " delta candidate " << jp << std::endl;


  //break track into segments and populate vectors defined above
  breakTrajInSegments(traj, breakpoints, seglens, cumseglens, seghits, cumseghits);

  //check if number of segments is less than a certain value
  if (seglens.size() < minNSegs_) {
    cout << "number of segments less than " << minNSegs_ << ", stopping fit" << endl;
    return recob::MCSFitResultGS(); }
  //check if total number of hits is less than a certain value
  if (cumseghits.back() < minHits_) {
    cout << "number of hits less than " << minHits_ << ", stopping fit" << endl;
    return recob::MCSFitResultGS(); }

  //define vector dthetaLin with expected scattering angles for linear fit
  vector<float> dthetaLinExp; dthetaLinExp.clear();
  //define vector dthetaLin with measured scattering angles for linear fit
  vector<float> dthetaLin; dthetaLin.clear();
  //define vector pcdir0, pcdir1 with directions of adjacent segments
  Vector_t pcdir0; Vector_t pcdir1; 

  //print range momentum [GeV/c]
  std::cout << "range momentum [GeV/c] " << rangeP << std::endl;
  //compute total energy [GeV] from range momentum [GeV/c]
  const double p2 = pow(rangeP, 2);
  const double m2 = pow(mass(pid), 2);
  const double E2 = sqrt(p2 + m2);

  for (unsigned int p = 0; p < seglens.size(); p++) {
    cout << "linear fit: iteration number " << p << endl;
    //perform linear regression of current segment and memorize its direction in pcdir1
    if (dimMode_ == 2) linearRegression2D(traj, breakpoints[p], breakpoints[p+1], pcdir1);
    else linearRegression(traj, breakpoints[p], breakpoints[p+1], pcdir1);

    if (p == 0) cout << "first iteration, scattering angle impossible to compute" << endl;
    if (p > 0) {
      //recall length [cm] of current segment, return null if less or equal than zero
      double Ls = seglens[p]; 
      cout << "current segment length [cm] = " << Ls << endl;
      if (Ls <= 0) {
        cout << "invalid length of current segment" << endl;
        dthetaLinExp.push_back(0.); }
      else {
        //compute energy [GeV] at current segment, return null if less than energy at rest
        const double Eij = GetE(E2, cumseglens[p], mass(pid));
        cout << "current energy [GeV] = " << Eij << endl;
        const double Eij2 = Eij * Eij;
        if (Eij2 <= m2) {
          cout << "invalid energy at current segment: less than energy at rest" << endl;
          dthetaLinExp.push_back(0.); }
        else {
          //compute momentum [GeV/c] and velocity of particle at current segment
          const double pij = sqrt(Eij2 - m2);
          cout << "current momentum [GeV/c] = " << Eij << endl;
          //compute velocity beta = v/c of particle at current segment
          const double beta = pij / Eij;
          //compute theta MCS [rad] just as claimed by Highland formula
          double thetamcs = 0.001 * ThetaExpected(traj, pij, beta, Ls);
          cout << p << "-th expected scattering angle [mrad] = " << thetamcs << endl;
          if (dimMode_ == 2) thetamcs *= 0.86; 
          dthetaLinExp.push_back(thetamcs); } }

      //check if current direction is backwards with respect to previous direction
      if (pcdir0.Dot(pcdir1) < 0.) pcdir1 *= -1.;

      //check if current or previous segment lengths have anomalous values
      if (seglens[p] <= 0. || seglens[p - 1] <= 0.) {
        cout << "WARNING! invalid segment length found!" << endl;
        dthetaLin.push_back(0); }
      else {
        cout << p << "-th direction = " << pcdir0 << endl;
        cout << p+1 << "-th direction = " << pcdir1 << endl;
        
        //check if current or previous direction have anomalous values
        auto trivial = Vector_t(0, 0, 0);
        if (pcdir0 == trivial || pcdir1 == trivial) dthetaLin.push_back(0);
        else {
          //compute scalar product of two directions, i.e. cosine of scattering angle
          double cosval = pcdir0.X() * pcdir1.X() + pcdir0.Y() * pcdir1.Y() + pcdir0.Z() * pcdir1.Z();
          if (cosval < -1.0) cosval = -1.0; 
          if (cosval > 1.0) cosval = 1.0;

          //compute scattering angle [rad] from its cosine and add it to vector dthetaLin
          double dt = acos(cosval);
          cout << p << "-th scattering angle [rad] = " << dt << endl;
          if (cutMode_ == 0) dthetaLin.push_back(dt);
          else if (cutMode_ == 1) {
            if (p >= seglens.size() - cutAngles_) dthetaLin.push_back(0);
            else dthetaLin.push_back(dt); }
          else if (cutMode_ == 2) {
            if (p > cutAngles_) dthetaLin.push_back(0);
            else dthetaLin.push_back(dt); } } } }

    pcdir0 = pcdir1; 
    cout << " " << endl; }

  //define vector dthetaPoly with expected scattering angles for linear fit
  vector<float> dthetaPolyExp; dthetaPolyExp.clear();
  //define vector dthetaPoly with measured scattering angles for polygonal fit
  vector<float> dthetaPoly; dthetaPoly.clear(); 
  //define vector barycenters, bary with barycenter positions
  vector<Vector_t> barycenters; Vector_t bary; 

  for (unsigned int p = 0; p < seglens.size(); p++) {
    cout << "poligonal fit: iteration number " << p << endl;
    //find barycenter position of current segment and memorize it in bary
    if (dimMode_ == 2) find2DSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary);
    else findSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary);
    barycenters.push_back(bary);

    if (p == 0) cout << "first iteration, scattering angle impossible to compute" << endl;
    if (p == 1) cout << "second iteration, scattering angle impossible to compute" << endl;
    if (p > 1) {
      //recall length [cm] of current segment, return null if less or equal than zero
      double Ls = seglens[p]; 
      cout << "current segment length [cm] = " << Ls << endl;
      if (Ls <= 0) {
        cout << "invalid length of current segment" << endl;
        dthetaPolyExp.push_back(0.); }
      else {
        //compute energy [GeV] at current segment, return null if less than energy at rest
        const double Eij = GetE(E2, cumseglens[p], mass(pid));
        cout << "current energy [GeV] = " << Eij << endl;
        const double Eij2 = Eij * Eij;
        if (Eij2 <= m2) {
          cout << "invalid energy at current segment: less than energy at rest" << endl;
          dthetaPolyExp.push_back(0.); }
        else {
          //compute momentum [GeV/c] and velocity of particle at current segment
          const double pij = sqrt(Eij2 - m2);
          cout << "current momentum [GeV/c] = " << Eij << endl;
          //compute velocity beta = v/c of particle at current segment
          const double beta = pij / Eij;
          //compute theta MCS [rad] just as claimed by Highland formula
          double thetamcs = 0.001 * ThetaExpected(traj, pij, beta, Ls);
          cout << p << "-th expected scattering angle [mrad] = " << thetamcs << endl;
          if (dimMode_ == 2) thetamcs *= 0.74; 
          dthetaPolyExp.push_back(thetamcs); } }

      //check if current or previous segment lengths have anomalous values
      if (seglens[p] <= 0. || seglens[p - 1] <= 0. || seglens[p - 2] <= 0.) {
        cout << "WARNING! invalid segment length found!" << endl;
        dthetaPoly.push_back(0); } 
      else {
        cout << p - 1 << "-th barycenter = " << barycenters[p - 2] << endl;
        cout << p << "-th barycenter = " << barycenters[p - 1] << endl;
        cout << p + 1 << "-th barycenter = " << barycenters[p] << endl;

        //check if current or previous barycenters have anomalous values
        auto trivial = Vector_t(0, 0, 0);
        if (barycenters[p - 2] == trivial || barycenters[p - 1] == trivial || barycenters[p] == trivial) dthetaPoly.push_back(0);
        else {
          //compute direction between first (p-2) and second (p-1) barycenter, then normalize it
          Vector_t dbcm = barycenters[p-1] - barycenters[p-2];
          float normm = sqrt(dbcm.X() * dbcm.X() + dbcm.Y() * dbcm.Y() + dbcm.Z() * dbcm.Z());
          dbcm /= normm;
          cout << "direction between " << p-1 << "-th and " << p << "-th barycenter = " << dbcm << endl;

          //compute direction between second (p-1) and third (p) barycenter, then normalize it
          Vector_t dbcp = barycenters[p] - barycenters[p-1];
          float normp = sqrt(dbcp.X() * dbcp.X() + dbcp.Y() * dbcp.Y() + dbcp.Z() * dbcp.Z());
          dbcp /= normp;
          cout << "direction between " << p << "-th and " << p+1 << "-th barycenter = " << dbcp << endl;

          //compute scalar product of two directions, i.e. cosine of scattering angle
          double cosval = dbcp.X() * dbcm.X() + dbcp.Y() * dbcm.Y() + dbcp.Z() * dbcm.Z();
          if (cosval < -1.0) cosval = -1.0; 
          if (cosval > 1.0) cosval = 1.0;

          //compute scattering angle [rad] from its cosine and add it to vector dthetaPoly
          double dt = acos(cosval);
          cout << p-1 << "-th scattering angle [rad] = " << dt << endl;
          if (cutMode_ == 0) dthetaPoly.push_back(dt);
          else if (cutMode_ == 1) {
            if (p >= seglens.size() - cutAngles_) dthetaPoly.push_back(0);
            else dthetaPoly.push_back(dt); }
          else if (cutMode_ == 2) {
            if (p > cutAngles_) dthetaPoly.push_back(0);
            else dthetaPoly.push_back(dt); } } } }
    cout << " " << endl; }

  //print all vectors defined and populated at this point
  cout << "segment lengths [cm] = "; for (auto i : seglens) cout << i << ' '; cout << endl;
  cout << "cumulative segment lengths [cm] = "; for (auto i : cumseglens) cout << i << ' '; cout << endl;
  cout << "segment number of hits = "; for (auto i : seghits) cout << i << ' '; cout << endl;
  cout << "cumulative segment number of hits = "; for (auto i : cumseghits) cout << i << ' '; cout << endl;
  cout << "expected scattering angles, linear fit [rad] = "; for (auto i : dthetaLinExp) cout << i << ' '; cout << endl;
  cout << "measured scattering angles, linear fit [rad] = "; for (auto i : dthetaLin) cout << i << ' '; cout << endl;
  cout << "expected scattering angles, poligonal fit [rad] = "; for (auto i : dthetaPolyExp) cout << i << ' '; cout << endl;
  cout << "measured scattering angles, poligonal fit [rad] = "; for (auto i : dthetaPoly) cout << i << ' '; cout << endl;
  cout << " " << endl;
  cout << "check if hit is delta = "; for (auto i : isDelta) cout << i << ' '; cout << endl;
  cout << " " << endl;

  /*
  vector<float> dthetaLinClean = CleanThetaLin(dthetaLin);
  vector<float> dthetaPolyClean = CleanThetaLin(dthetaPoly);
  dthetaLin = dthetaLinClean;
  dthetaPoly = dthetaPolyClean;
  */

  //check if there is minimum number of angles > 0
  unsigned int firstseg = 0; unsigned int lastseg = 0; 
  bool checkLin = false; bool checkPoly = false;
  ThetaCheck(dthetaLin, dthetaPoly, checkLin, checkPoly, firstseg, lastseg);
  if (!checkLin) {
    cout << "error: number of scattering angles from linear fit is not greater or equal " << minNAngs_ << endl;
    cout << "fit stops here, return null" << endl;
    return recob::MCSFitResultGS(); }
  if (!checkPoly) {
    cout << "error: number of scattering angles from poligonal fit is not greater or equal " << minNAngs_-1 << endl;
    cout << "fit stops here, return null" << endl;
    return recob::MCSFitResultGS(); }

  //perform a c2fit scan over all the computed angles
  float c2sigma = 0.05; 
  const ScanResult fitResult = C2Fit(traj, breakpoints, seglens, cumseglens, seghits, cumseghits, dthetaLin, dthetaPoly, pid, c2sigma);

  const float c2prange = C2Function(traj, breakpoints, seglens, cumseglens, seghits, cumseghits, dthetaLin, dthetaPoly, pid, rangeP*1000.);
  std::cout << "c2 function at range momentum = " << c2prange << std::endl;

  //recall sigma3p and L1D, L2D, L3D functions
  double sigma3p = PrintD3P();
  double L1D = length1D(traj, planeMode_); 
  double L2D = length2D(traj, planeMode_);
  double L3D = length3D(traj, planeMode_);

  //return pid, bestp and errp, sigma3p, length 1D 2D 3D, segment lens and number of hits, scattering linear-fit and polygonal angles, stopping check
  return recob::MCSFitResultGS(
    pid, 
    fitResult.bestp, fitResult.errp, fitResult.minp, fitResult.maxp,
    fitResult.alpha, fitResult.dalpha, fitResult.beta, fitResult.dbeta,
    fitResult.testp, fitResult.c2function, 
    c2prange, 
    sigma3p, 
    L1D, L2D, L3D, 
    seglens, cumseglens,
    seghits, cumseghits, 
    dthetaLinExp, dthetaLin, 
    dthetaPolyExp, dthetaPoly, 
    stop,
    isDelta); 
}

//break input trajectory into smaller pieces called segments 
void TrajectoryMCSFitterICARUS::breakTrajInSegments(const recob::TrackTrajectory& traj, vector<size_t>& breakpoints, vector<float>& seglens, vector<float>& cumseglens, vector<int>& seghits, vector<int>& cumseghits) const {
  cout << "track length [cm] = " << traj.Length() << endl;
  cout << "segment length [cm] = " << segLen_ << endl;
  cout << " " << endl;

  //initialize vectors cumseglens and cumseghits (to use .back() later)
  cumseglens.push_back(0.); cumseghits.push_back(0);
  double seglen = 0.; int seghit = 0;

  //determine index of first valid point of track, add it to vector breakpoints
  size_t index = traj.FirstValidPoint(); 
  breakpoints.push_back(index);

  //define position of first valid point
  auto pos0 = traj.LocationAtPoint(index);

  //find next valid point index
  index = traj.NextValidPoint(index + 1);

  while (index != recob::TrackTrajectory::InvalidIndex) {
    //define position of current valid point, then update segment length and position of next valid point
    auto pos1 = traj.LocationAtPoint(index);
    seglen += (pos1 - pos0).R();
    pos0 = pos1;

    //update number of hits depending on 2D or 3D
    if ((dimMode_ == 2) && (isinplane(index, planeMode_) && isintpc(index, lasttpc(traj)) && CathodeCheck(traj, index))) seghit++;
    if ((dimMode_ == 3) && CathodeCheck(traj, index)) seghit++;
    
    //break into a new segment if current segment length is greater than required segment length
    if (seglen >= segLen_) {
      //add current valid point index to vector breakpoints
      breakpoints.push_back(index);

      //add current segment length [cm] and cumulative segment length [cm] to vector seglens, cumseglens
      if (seghit >= minHitsPerSegment_) seglens.push_back(seglen);
      else seglens.push_back(0.);
      cumseglens.push_back(cumseglens.back() + seglen);
      seglen = 0.;

      //add current segment number of hits and cumulative number of hits to vector seghits, cumseghits
      seghits.push_back(seghit);
      cumseghits.push_back(cumseghits.back() + seghit);
      seghit = 0; 
    }
    index = traj.NextValidPoint(index + 1);
  }

  //then add last segment
  if (seglen > 0.) {
    //add last valid point index to vector breakpoints
    breakpoints.push_back(traj.LastValidPoint() + 1);

    //add last segment length [cm] and total segment length [cm] to vector seglens, cumseglens
    seglens.push_back(seglen);
    cumseglens.push_back(cumseglens.back() + seglen);

    //add last segment number of hits and total number of hits to vector seghits, cumseghits
    seghits.push_back(seghit);
    cumseghits.push_back(cumseghits.back() + seghit);
  }
  return;
}

//find barycenter position of input segment (useful for poligonal fit) in 3D
void TrajectoryMCSFitterICARUS::findSegmentBarycenter(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& bary) const {
  //initalize number of points in segment
  int npoints = 0;

  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;
  
  //determine index of first valid point of the segment
  size_t index = firstPoint;
  while (index < lastPoint) {
    if (CathodeCheck(traj, index)) {
      //add position of current valid point to vector middlePointCalc
      middlePointCalc.add(traj.LocationAtPoint(index));
      npoints++; }
    index = traj.NextValidPoint(index + 1); }

  //determine position of segment barycenter from position of valid points in the segment
  bary = middlePointCalc.middlePoint();
}

//find average direction of trajectory between firstPoint and lastPoint in 3D
void TrajectoryMCSFitterICARUS::linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir) const {
  //initalize number of points in segment
  int npoints = 0;

  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;

  //determine index of first valid point of the segment
  size_t index = firstPoint;
  while (index < lastPoint) {
    if (CathodeCheck(traj, index)) {
      //add position of current valid point to vector middlePointCalc
      middlePointCalc.add(traj.LocationAtPoint(index));
      npoints++; }
    index = traj.NextValidPoint(index + 1); }

  //check if number of points is greater than zero, otherwise return trivial vector
  if (npoints == 0) pcdir = Vector_t();
  else {
    //determine position of segment barycenter from position of valid points in the segment
    const auto avgpos = middlePointCalc.middlePoint();
    //define normalization factor as inverse of number of valid points
    const double norm = 1./double(npoints);
    //define covariance matrix m as symmetric matrix 3x3
    TMatrixDSym m(3);

    index = firstPoint;
    while (index < lastPoint) {
      if (CathodeCheck(traj, index)) {
        //determine position of current valid point
        auto p = traj.LocationAtPoint(index);

        //compute coordinate differences between current valid point and average point
        const double xxw0 = p.X() - avgpos.X();
        const double yyw0 = p.Y() - avgpos.Y();
        const double zzw0 = p.Z() - avgpos.Z();

        //update covariance matrix values with normalized values of coordinate differences above
        m(0, 0) += xxw0 * xxw0 * norm; m(0, 1) += xxw0 * yyw0 * norm; m(0, 2) += xxw0 * zzw0 * norm;
        m(1, 0) += yyw0 * xxw0 * norm; m(1, 1) += yyw0 * yyw0 * norm; m(1, 2) += yyw0 * zzw0 * norm;
        m(2, 0) += zzw0 * xxw0 * norm; m(2, 1) += zzw0 * yyw0 * norm; m(2, 2) += zzw0 * zzw0 * norm; }

      index = traj.NextValidPoint(index + 1); }

    //compute eigenvalues and eigenvectors of covariance matrix m
    const TMatrixDSymEigen me(m);
    const auto& eigenval = me.GetEigenValues(); const auto& eigenvec = me.GetEigenVectors();
    int maxevalidx = 0; double maxeval = eigenval(0);
    for (int i = 1; i < 3; ++i) {
      if (eigenval(i) > maxeval) {
        maxevalidx = i;
        maxeval = eigenval(i); } } 

    //update segment direction with eigenvector associated to maximum eigenvalue
    pcdir = Vector_t(eigenvec(0, maxevalidx), eigenvec(1, maxevalidx), eigenvec(2, maxevalidx));

    //reverse segment direction if current direction is opposite to trajectory direction in the first valid point
    if (traj.DirectionAtPoint(firstPoint).Dot(pcdir) < 0.) pcdir *= -1.;
  }
}

//find barycenter position of input segment (useful for poligonal fit) in 2D
void TrajectoryMCSFitterICARUS::find2DSegmentBarycenter(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& bary2D) const {
  //initalize number of points in segment
  int npoints = 0;

  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;
  
  //determine index of first valid point of the segment
  size_t index = firstPoint;
  while (index < lastPoint) {
    //check if current valid point is in input plane and last tpc
    if (isinplane(index, planeMode_) && isintpc(index, lasttpc(traj)) && CathodeCheck(traj, index)) {
      //determine position of current valid point
      auto p = hit2d(traj, index, planeMode_, lasttpc(traj));
      //add position of current valid point to vector middlePointCalc
      middlePointCalc.add(p);
      npoints++; }
    index = traj.NextValidPoint(index + 1); }

  //check if number of points is greater than zero, otherwise return trivial vector
  if (npoints == 0) bary2D = Vector_t();
  //determine position of segment barycenter from position of valid points in the segment
  else bary2D = middlePointCalc.middlePoint();
}

//find average direction of trajectory between firstPoint and lastPoint in 2D
void TrajectoryMCSFitterICARUS::linearRegression2D(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir2D) const {
  //initalize number of points in segment
  int npoints = 0;

  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;

  //determine index of first valid point of the segment
  size_t index = firstPoint;
  while (index < lastPoint) {
    //check if current valid point is in input plane and last tpc
    if (isinplane(index, planeMode_) && isintpc(index, lasttpc(traj)) && CathodeCheck(traj, index)) {
      //determine position of current valid point
      auto p = hit2d(traj, index, planeMode_, lasttpc(traj));
      //add position of current valid point to vector middlePointCalc
      middlePointCalc.add(p);
      npoints++; } 
    index = traj.NextValidPoint(index + 1); }

  //check if number of points is greater than zero, otherwise return trivial vector
  if (npoints == 0) pcdir2D = Vector_t();
  else {
    //determine position of segment barycenter from position of valid points in the segment
    const auto avgpos = middlePointCalc.middlePoint();
    //define normalization factor as inverse of number of valid points
    const double norm = 1./double(npoints);
    //define covariance matrix m as symmetric matrix 3x3
    TMatrixDSym m(3);

    index = firstPoint;
    while (index < lastPoint) {
      //check if current valid point is in input plane and last tpc
      if (isinplane(index, planeMode_) && isintpc(index, lasttpc(traj)) && CathodeCheck(traj, index)) {
        //determine position of current valid point
        auto p = hit2d(traj, index, planeMode_, lasttpc(traj));

        //compute coordinate differences between current valid point and average point
        const double xxw0 = p.X() - avgpos.X();
        const double yyw0 = p.Y() - avgpos.Y();
        const double zzw0 = p.Z() - avgpos.Z();

        //update covariance matrix values with normalized values of coordinate differences above
        m(0, 0) += xxw0 * xxw0 * norm; m(0, 1) += xxw0 * yyw0 * norm; m(0, 2) += xxw0 * zzw0 * norm;
        m(1, 0) += yyw0 * xxw0 * norm; m(1, 1) += yyw0 * yyw0 * norm; m(1, 2) += yyw0 * zzw0 * norm;
        m(2, 0) += zzw0 * xxw0 * norm; m(2, 1) += zzw0 * yyw0 * norm; m(2, 2) += zzw0 * zzw0 * norm; }

      index = traj.NextValidPoint(index + 1); }

    //compute eigenvalues and eigenvectors of covariance matrix m
    const TMatrixDSymEigen me(m);
    const auto& eigenval = me.GetEigenValues(); const auto& eigenvec = me.GetEigenVectors();
    int maxevalidx = 0; double maxeval = eigenval(0);
    for (int i = 1; i < 3; ++i) {
      if (eigenval(i) > maxeval) {
        maxevalidx = i;
        maxeval = eigenval(i); } } 
    
    //update segment direction with eigenvector associated to maximum eigenvalue
    pcdir2D = Vector_t(eigenvec(0, maxevalidx), eigenvec(1, maxevalidx), eigenvec(2, maxevalidx));
  }
}

//find most probable value of energy loss according to Landau distribution
double TrajectoryMCSFitterICARUS::energyLossLandau(const double mass2, const double e2, const double x) const {
  //return null energy if travelled distance is negative or null
  if (x <= 0.) return 0.;

  //define inverse of squared ionization potential 1/I2 [MeV-2]
  constexpr double Iinv2 = 1. / (188.E-6 * 188.E-6);
  //define matConst [mol cm-3] as product of density rho, atomic number Z, inverse of mass number 1/A 
  constexpr double matConst = 1.4 * 18. / 40.;
  //define electron mass [MeV/c2]
  constexpr double me = 0.511;
  //define K constant [MeV mol−1 cm2]
  constexpr double kappa = 0.307075;
  //define Landau parameter j
  constexpr double j = 0.200;

  //compute squared velocity beta2 and squared Lorentz factor gamma2
  const double beta2 = (e2 - mass2) / e2;
  const double gamma2 = 1. / (1.0 - beta2);

  //return most probable value of energy [GeV] according to Landau distribution
  const double epsilon = 0.5 * kappa * x * matConst / beta2;
  return 0.001 * epsilon * (log(2. * me * beta2 * gamma2 * epsilon * Iinv2) + j - beta2);
}

//find mean energy loss according to BetheBloch distribution
double TrajectoryMCSFitterICARUS::energyLossBetheBloch(const double mass, const double e2) const {
  //define inverse of ionization potential 1/I [MeV-1]
  constexpr double Iinv = 1. / 188.E-6;
  //define matConst [mol cm-3] as product of density rho, atomic number Z, inverse of mass number 1/A 
  constexpr double matConst = 1.4 * 18. / 40.;
  //define electron mass [MeV/c2]
  constexpr double me = 0.511;
  //define K constant [MeV mol−1 cm2]
  constexpr double kappa = 0.307075;

  //compute squared velocity beta2 and squared Lorentz factor gamma2
  const double beta2 = (e2 - mass * mass) / e2;
  const double gamma2 = 1. / (1.0 - beta2);

  //compute ratio between electron mass and particle mass
  const double massRatio = me / mass;
  //compute argument of logaritmic term 
  const double argument = (2. * me * gamma2 * beta2 * Iinv) / sqrt(1 + 2 * sqrt(gamma2) * massRatio + massRatio * massRatio);
  
  //initialize energy loss [MeV cm-1]
  double dedx = kappa * matConst / beta2;
  //return null energy if particle mass is null or if argument of logaritmic term is not greater than exp(beta2)
  if (mass == 0.0) return (0.0);
  if (argument <= exp(beta2)) return (0.0);
  else {
    dedx *= (log(argument) - beta2);
    //return null energy if this energy loss is negative
    if (dedx < 0.) return (0.0);
  }

  //return energy loss [GeV cm-1] according to Bethe Bloch
  return 0.001 * dedx;
}

//find energy of track at a certain distance travelled
double TrajectoryMCSFitterICARUS::GetE(const double initial_E, const double length_travelled, const double m) const {
  //compute energy if particle is considered MIP with constant energy loss
  if (eLossMode_ == 1) {
    //define energy loss [GeV cm-1] for MIP muon in liquid Argon
    constexpr double kcal = 0.002105;
    //return energy [GeV] as difference between initial energy and kcal * travelled distance [cm]
    return (initial_E - kcal * length_travelled);
  }

  //define step size as ratio between travelled distance [cm] and an arbitrary number
  const double step_size = length_travelled / nElossSteps_;
  //define initial energy [GeV] and squared mass [GeV2 c-4]
  double current_E = initial_E;
  const double m2 = m * m;

  //compute energy [GeV] if particle is considered with variable energy loss
  for (auto i = 0; i < nElossSteps_; ++i) {
    //compute energy [GeV] according to Bethe Bloch, negleting density correction
    if (eLossMode_ == 2) {
      double dedx = energyLossBetheBloch(m, current_E * current_E);
      current_E -= (dedx * step_size);
    }

    //compute energy [GeV] according to Landau distribution (default configuration!)
    else {
      current_E -= energyLossLandau(m2, current_E * current_E, step_size);
    }

    //return null energy if current energy is less or equal than energy at rest
    if (current_E <= m) {return 0.;}
  }

  //return current energy of particle [GeV] in function of travelled distance [cm]
  return current_E;
}

//compute delta3p [mm] used for measurement error along drift direction
void TrajectoryMCSFitterICARUS::ComputeD3P(int plane) {
  //define delta3p for every triplet of consecutive hits
  float d3p;
  //define residual value res for every triplet of consecutive hits
  double res;
  //define vectors h0 that will be filled with residuals, and alf that will be filled with parameter alfa
  vector<double> h0;
  vector<double> alf;

  //define vector hits that contains hits associated to a certain plane (Induction-1, Induction-2, Collection)
  vector<recob::Hit> hits;
  if (plane == 0) hits = hits2dI1;
  if (plane == 1) hits = hits2dI2;
  if (plane == 2) hits = hits2dC;
  //return null if there are no hits in the selected plane
  if (!hits.size()) return;   

  //make a ROOT histogram hd3pv to visualize residuals, 100 bins, range [-5, 5]
  TH1D* hd3pv = new TH1D("hd3pv", "hd3pv", 100, -5., 5.);
  //iterate over every triplet of consecutive hits
  for (unsigned int j = 0; j < hits.size() - 2; j += 3) {
    double a;
    //compute residual and alfa for the current triplet
    res = computeResidual(j, a, hits);
    //save absolute value of residuals if less than 5 mm and parameter alfa
    if (abs(res) < 5) { 
      h0.push_back(res);
      alf.push_back(a);
      hd3pv->Fill(res);
    }
  }

  //make a ROOT file d3phisto.root in UPDATE mode and save residual histogram in this ROOT file
  bool writeHisto = true;
  if (writeHisto) {
    TFile *f = new TFile("d3phisto.root","RECREATE");
    hd3pv->Write();
    f->Close();
    f->Delete();
  }

  //return d3p = 0.4 mm as default value if there are no valid residuals
  if (!h0.size()) d3p = 0.4;
  //otherwise return d3p as rms of residuals saved in histogram 
  else d3p = hd3pv->GetRMS();
  
  //fill d3p with right quantity depending on the plane
  if (plane == 0) d3pI1 = d3p;
  if (plane == 1) d3pI2 = d3p;
  if (plane == 2) d3pC = d3p;
}

//compute delta3p [mm] in 3D used for measurement error along drift direction (not used)
float TrajectoryMCSFitterICARUS::ComputeD3P3D(const recob::TrackTrajectory& tr) const {
  //define delta3p for every triplet of consecutive hits
  float d3p;
  //define residual value res for every triplet of consecutive hits
  double res;
  //define vectors h0 that will be filled with residuals, and alf that will be filled with parameter alfa
  vector<double> h0;
  vector<double> alf;

  //make a ROOT histogram hd3pv to visualize residuals, 100 bins, range [-5, 5]
  TH1D* hd3pv = new TH1D("hd3pv", "hd3pv", 100, -5., 5.);
  //define index of first valid point and iterate over every triplet of consecutive hits
  unsigned int nextValid = tr.NextValidPoint(0);
  while (nextValid != tr.LastValidPoint()) {
    double a;
    //compute residual and alfa for the current triplet of consecutive hits
    res = computeResidual3D(tr, nextValid, a);
    //save absolute value of residuals if less than 5 mm and parameter alfa
    if (abs(res) < 5) {
      h0.push_back(res);
      alf.push_back(a);
      hd3pv->Fill(res);
    }
    nextValid = tr.NextValidPoint(nextValid + 1);
  }

  //make a ROOT file d3phisto.root in UPDATE mode and save residual histogram in this ROOT file
  bool writeHisto = true;
  if(writeHisto) {
    TFile *f = new TFile("d3phisto3D.root","RECREATE");
    hd3pv->Write();
    f->Close();
    f->Delete();
  }

  //return d3p = 0.4 mm as default value if there are no valid residuals
  if (!h0.size()) d3p = 0.4;
  //otherwise return d3p as rms of residuals saved in histogram 
  else d3p = hd3pv->GetRMS();

  //return 3D delta3p
  return d3p;
}

//compute residual used in delta3p computation
double TrajectoryMCSFitterICARUS::computeResidual(int i, double& alfa, vector<recob::Hit> hits2d) const {
  //get three consecutive hits from vector hits2d
  recob::Hit h0 = hits2d.at(i);
  recob::Hit h1 = hits2d.at(i + 1);
  recob::Hit h2 = hits2d.at(i + 2);

  //compute 2D coordinates from wire ID and signal peak time for three consecutive hits
  float x0 = h0.WireID().Wire * 3; auto y0 = h0.PeakTime() * 0.622;
  float x1 = h1.WireID().Wire * 3; auto y1 = h1.PeakTime() * 0.622;
  float x2 = h2.WireID().Wire * 3; auto y2 = h2.PeakTime() * 0.622;

  //interpolate y1 for x = x1 starting from points (x0, y0) and (x2, y2)
  double ym = y0 + (y2 - y0) * (x1 - x0) / (x2 - x0);
  //check if differences between wire coordinates are not too small
  if (abs(x2 - x0) < 0.001) return -999;
  if (abs(x2 - x1) < 0.001) return -999;
  if (abs(x1 - x0) < 0.001) return -999;

  //compute parameter alfa taking into account point geometry
  double K = (x1 - x0) / (x2 - x0);
  alfa = 1 / sqrt(1 + K * K + (1 - K) * (1 - K));

  //compute coordinate differences between measured and expected value
  double dy = y1 - ym;
  return dy * alfa;
}

//compute 3D residual used in delta3p computation (not used)
double TrajectoryMCSFitterICARUS::computeResidual3D(recob::TrackTrajectory tr, int i, double& alfa) const {
  //get three consecutive hits from track 
  auto p0 = tr.LocationAtPoint(i);
  auto p1 = tr.LocationAtPoint(i + 1);
  auto p2 = tr.LocationAtPoint(i + 2);

  //compute 3D coordinates x, y, z for three consecutive hits
  auto x0 = 10 * p0.X(); auto y0 = 10 * p0.Y(); auto z0 = 10 * p0.Z();
  auto x1 = 10 * p1.X(); auto y1 = 10 * p1.Y(); auto z1 = 10 * p1.Z();
  auto x2 = 10 * p2.X(); auto y2 = 10 * p2.Y(); auto z2 = 10 * p2.Z();

  //interpolate y1 for x = x1 starting from points (x0, y0) and (x2, y2)
  double xmy = x0 + (y1 - y0) / (y2 - y0) * (x2 - x0);
  //interpolate z1 for x = x1 starting from points (x0, z0) and (x2, z2)
  double xmz = x0 + (z1 - z0) / (z2 - z0) * (x2 - x0);
  //check if distance between points are not too small
  if ((p0 - p2).Mag2() < 0.001) return -999;
  if ((p1 - p2).Mag2() < 0.001) return -999;
  if ((p0 - p1).Mag2() < 0.001) return -999;

  //compute parameter alfa taking into account point geometry
  double Ky = (y1 - y0) / (y2 - y0); double alfay = 1 / sqrt(1 + Ky * Ky + (1 - Ky) * (1 - Ky));
  double Kz = (z1 - z0) / (z2 - z0); double alfaz = 1 / sqrt(1 + Kz * Kz + (1 - Kz) * (1 - Kz));
  alfa = sqrt(alfay * alfay + alfaz * alfaz);

  //compute coordinate differences between measured and expected value
  double xm = 0.5 * (xmy + xmz);
  double dx = x1 - xm;
	return dx * alfa;
}

//build c2 function used in MCS fit
const double TrajectoryMCSFitterICARUS::C2Function(const recob::TrackTrajectory& traj, std::vector<size_t>& breakpoints, std::vector<float>& seglens, std::vector<float>& cumseglens, std::vector<int>& seghits, std::vector<int>& cumseghits, std::vector<float> dthetaLin, std::vector<float> dthetaPoly, int pid, double p0) const {
  //define thetaexp, thetamcs, thetaerr that represent expected, measured and error angle
  double thetaexp, thetamcs, thetaerr;

  //define ttallLin, ttallPoly 
  vector<float> ttall; vector<float> ttallLin; vector<float> ttallPoly;
  ttall.clear(); ttallLin.clear(); ttallPoly.clear();

  //check if there is minimum number of angles > 0
  unsigned int firstseg = 0; unsigned int lastseg = 0; 
  bool checkLin = false; bool checkPoly = false;
  ThetaCheck(dthetaLin, dthetaPoly, checkLin, checkPoly, firstseg, lastseg);
  bool check = checkLin && checkPoly;
  cout << "first = " << firstseg << " last = " << lastseg << endl;
  if (!check) return 0;

  //compute number of segments as difference (+1) between first and last segment index
  unsigned int nseg = lastseg - firstseg + 1;
  if (nseg < 7) return 0;

  //define covariance matrices for both linear fit and polygonal
  TMatrixDSym mat(2*nseg-3);
  TMatrixDSym matFull(2*nseg-3);

  //define covariance matrices for linear fit
  TMatrixDSym matLin(nseg-1);
  TMatrixDSym matMCSLin(nseg-1);
  TMatrixDSym matErrLin(nseg-1);
  TMatrixDSym matOffDiagLin(nseg-1);

  //define covariance matrices for polygonal fit
  TMatrixDSym matPoly(nseg-2);
  TMatrixDSym matMCSPoly(nseg-2);
  TMatrixDSym matErrPoly(nseg-2);
  TMatrixDSym matOffDiagPoly(nseg-2);

  //define covariance matrix with mix terms
  TMatrixD matMCSMixed(nseg-1, nseg-2);

  //recall measurement error from delta3p [cm]
  double sigma3p = PrintD3P();
  cout << "sigma3p = " << sigma3p << endl;

  //define propagation constant for measurement error associated to delta3p
  double kLin = sqrt(24.); double kPoly = sqrt(6.);

  //convert input momentum p0 [MeV/c] to momentum p1 [GeV/c]
  const double p1 = 0.001 * p0;
  //compute total energy [GeV] from momentum p1 [GeV/c]
  const double p2 = p1 * p1;
  const double m2 = mass(pid) * mass(pid);
  const double Etot = sqrt(p2 + m2);

  //notify computation of scattering angles for linear fit has started
  cout << "starting computation of scattering angles for linear fit ..." << endl;
  cout << " " << endl;

  for (unsigned int jp = firstseg; jp < lastseg; jp++) {
    cout << "linear fit: segment number " << jp << endl;
    
    //compute energy [GeV] at current segment, return null if less than energy at rest
    const double Eij = GetE(Etot, cumseglens[jp], mass(pid));
    cout << "current energy [GeV] = " << Eij << endl;
    const double Eij2 = Eij * Eij;
    if (Eij2 <= m2) {
      cout << "invalid energy at current segment: less than energy at rest" << endl;
      return 0; }
    
    //compute momentum [GeV/c] and velocity of particle at current segment
    const double pij = sqrt(Eij2 - m2);
    cout << "current momentum [GeV/c] = " << pij << endl;

    //compute velocity beta = v/c of particle at current segment
    const double beta = pij / Eij;
    cout << "current beta = v/c = " << beta << endl;

    //recall length [cm] of current segment, return null if less or equal than zero
    double Ls = seglens[jp]; 
    cout << "current segment length [cm] = " << Ls << endl;
    if (Ls <= 0) {
      cout << "invalid length of current segment" << endl;
      return 0; }

    //recall number of hits of current segment, return null if less or equal than zero
    int nhit = seghits[jp];
    cout << "current segment number of hits = " << nhit << endl;
    if (nhit <= 0) {
      cout << "invalid number of hits of current segment" << endl;
      return 0; }

    //update ttall vector with measured scattering angle [mrad]
    thetaexp = 1000. * dthetaLin[jp];
    ttallLin.push_back(thetaexp);
    cout << "current theta MCS measured [mrad] = " << thetaexp << endl;

    //compute theta MCS [mrad] just as claimed by Highland formula
    thetamcs = ThetaExpected(traj, pij, beta, Ls);
    if (dimMode_ == 2) thetamcs *= 0.86;
    cout << "current theta MCS expected from Highland [mrad] = " << thetamcs << endl;

    //compute theta error [mrad] associated to delta3p
    thetaerr = 1000 * (kLin * sigma3p)/(Ls * sqrt(float(nhit))) * ThetaErrFactor(traj, planeMode_);
    cout << "current theta MCS error from delta3p [mrad] = " << thetaerr << endl;

    //update matrix matMCSLin, matErrLin with on-diag MCS and error terms
    matMCSLin(jp-firstseg, jp-firstseg) = thetamcs * thetamcs;
    matErrLin(jp-firstseg, jp-firstseg) = thetaerr * thetaerr;

    //update matrix matOffDiagLin with off-diag MCS and error terms
    FillOffDiagMCSLin(matOffDiagLin, matMCSLin, jp-firstseg);
    FillOffDiagErrLin(matOffDiagLin, matErrLin, jp-firstseg);
    cout << " " << endl; }

  for (unsigned int jp = firstseg; jp < lastseg - 1; jp++) {
    cout << "poligonal fit: segment number " << jp << endl;
    
    //compute energy [GeV] at current segment, return null if less than energy at rest
    const double Eij = GetE(Etot, cumseglens[jp], mass(pid));
    cout << "current energy [GeV] = " << Eij << endl;
    const double Eij2 = Eij * Eij;
    if (Eij2 <= m2) {
      cout << "invalid energy at current segment: less than energy at rest" << endl;
      return 0; }
    
    //compute momentum [GeV/c] and velocity of particle at current segment
    const double pij = sqrt(Eij2 - m2);
    cout << "current momentum [GeV/c] = " << pij << endl;

    //compute velocity beta = v/c of particle at current segment
    const double beta = pij / Eij;
    cout << "current beta = v/c = " << beta << endl;

    //recall length [cm] of current segment, return null if less or equal than zero
    double Ls = seglens[jp]; 
    cout << "current segment length [cm] = " << Ls << endl;
    if (Ls <= 0) {
      cout << "invalid length of current segment" << endl;
      return 0; }

    //recall number of hits of current segment, return null if less or equal than zero
    int nhit = seghits[jp];
    cout << "current segment number of hits = " << nhit << endl;
    if (nhit <= 0) {
      cout << "invalid number of hits of current segment" << endl;
      return 0; }

    //update ttall vector with measured scattering angle [mrad]
    thetaexp = 1000. * dthetaPoly[jp];
    ttallPoly.push_back(thetaexp);
    cout << "current theta MCS measured [mrad] = " << thetaexp << endl;

    //compute theta MCS [mrad] just as claimed by Highland formula
    thetamcs = ThetaExpected(traj, pij, beta, Ls);
    if (dimMode_ == 2) thetamcs *= 0.74;
    cout << "current theta MCS expected from Highland [mrad] = " << thetamcs << endl;

    //compute theta error [mrad] associated to delta3p
    thetaerr = 1000 * (kPoly * sigma3p)/(Ls * sqrt(float(nhit))) * ThetaErrFactor(traj, planeMode_);
    cout << "current theta MCS error from delta3p [mrad] = " << thetaerr << endl;

    //update matrix matMCSPoly, matErrPoly with on-diag MCS and error terms
    matMCSPoly(jp-firstseg, jp-firstseg) = thetamcs * thetamcs;
    matErrPoly(jp-firstseg, jp-firstseg) = thetaerr * thetaerr;

    //update matrix matOffDiagPoly with off-diag MCS and error terms
    FillOffDiagMCSPoly(matOffDiagPoly, matMCSPoly, jp-firstseg); 
    FillOffDiagErrPoly(matOffDiagPoly, matErrPoly, jp-firstseg);
    cout << " " << endl; }

  //fill matrix matLin with on-diag, off-diag MCS and error terms
  for (int jm = 0; jm < matMCSLin.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matMCSLin.GetNrows(); jmm++) {  
      matLin(jm, jmm) = matMCSLin(jm, jmm) + matErrLin(jm, jmm) + matOffDiagLin(jm, jmm); } }

  //fill matrix matPoly with on-diag, off-diag MCS and error terms
  for (int jm = 0; jm < matMCSPoly.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matMCSPoly.GetNrows(); jmm++) {
      matPoly(jm, jmm) = matMCSPoly(jm, jmm) + matErrPoly(jm, jmm) + matOffDiagPoly(jm, jmm); } }

  //fill matrix matMCSMixed with mixed terms
  FillMCSMixedTerms(matMCSMixed, matMCSLin, matMCSPoly);

  //merge matrices matLin, matPoly, matMCSMixed into matFull
  for (int jm = 0; jm < matLin.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matLin.GetNrows(); jmm++) {
      matFull(jm, jmm) = matLin(jm, jmm); } }
  for (int jm = 0; jm < matPoly.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matPoly.GetNrows(); jmm++) {
      matFull(jm + matLin.GetNrows(), jmm + matLin.GetNrows()) = matPoly(jm, jmm); } }
  for (int jm = 0; jm < matPoly.GetNrows(); jm++) {
      for (int jmm = 0; jmm < matLin.GetNrows(); jmm++) {
        matFull(jm + matLin.GetNrows(), jmm) = matMCSMixed(jmm, jm); } }
    for (int jm = 0; jm < matLin.GetNrows(); jm++) {
      for (int jmm = 0; jmm < matPoly.GetNrows(); jmm++) {
        matFull(jm, jmm + matLin.GetNrows()) = matMCSMixed(jm, jmm); } } 

  //in case of both linear and polygonal fit, take matFull
  if (fitMode_ == 0) {
    mat = matFull;
    ttall.insert(ttall.end(), ttallLin.begin(), ttallLin.end());
    ttall.insert(ttall.end(), ttallPoly.begin(), ttallPoly.end()); }
  //in case of only linear fit, take matLin
  if (fitMode_ == 1) {
    mat.ResizeTo(matMCSLin.GetNrows(), matMCSLin.GetNrows());
    mat = matLin; 
    ttall = ttallLin; }
  //in case of only polygonal fit, take matPoly
  if (fitMode_ == 2) {
    mat.ResizeTo(matMCSPoly.GetNrows(), matMCSPoly.GetNrows());
    mat = matPoly;
    ttall = ttallPoly; }
  
  //initialize matrix cov as a copy of matrix mat
  TMatrixD cov = mat;
  //define column matrix vtall with same size of vector ttall, and fill it with elements of vector ttall
  TMatrixD vtall(ttall.size(), 1);
  for (unsigned int jv = 0; jv < ttall.size(); jv++) vtall(jv, 0) = ttall[jv]; 
  //define row matrix tvtall as the transpose of column matrix vtall
  TMatrixD tvtall(TMatrixD::kTransposed, vtall); 
  //define matrix invcov as the inverse of matrix cov
  TMatrixD invcov = cov.Invert();
  //define column matrix vtcov as matrix product of invcov and vtall
  TMatrixD vtcov = invcov * vtall;

  //initialize vtcovmed and vector terms
  double vtcovmed = 0;
  vector<double> terms;
  //iterate over elements of row matrix tvtall and column matrix vtcov
  for (unsigned int jv = 0; jv < ttall.size(); jv++) {
    double term = tvtall(0, jv) * vtcov(jv, 0);
    terms.push_back(term);
    vtcovmed += term; }

  //define column matrix mterms as matrix product of tvtall and vtcov
  TMatrixD mterms(terms.size(), 1);
  for (unsigned int jv = 0; jv < terms.size(); jv++) mterms(jv, 0) = terms[jv]; 
  cout << "column matrix mterms, with c2 function terms = " << endl; mterms.Print();
  //compute c2 function without trunc as average of mterms
  vtcovmed /= (ttall.size());
  cout << "c2 function without trunc = " << vtcovmed << endl;

  //define vector ttrunc that will be filled with accepted terms after cut
  vector<double> ttrunc;
  //define vector tails that will be filles with indexes of rejected terms
  vector<int> tails;
  //define ttrunctot as sum of accepted terms
  double ttrunctot = 0;
  //define vector atrunc that will be filled with ttall values of accepted terms
  vector<double> atrunc;

  //iterate over elements of vector terms
  //double gaus = 3.;
  double threshold = 100; 
  for (unsigned int jt = 0; jt < ttall.size(); jt++) {
    //check if element jt-th of vector terms is under threshold and observed angle is positive
    //if (terms[jt] < vtcovmed * pow(gaus, 2) && ttall[jt] > 0.) {
    if (ttall[jt] < threshold && ttall[jt] > 0.) {
      ttrunctot += terms[jt];
      atrunc.push_back(ttall[jt]);
      ttrunc.push_back(terms[jt]); } 
    else { 
      tails.push_back(jt); 
      cout << "measured angle = " << ttall[jt] << " with index = " << jt << " is greater than " << threshold << " rad" << endl;} }
  
  //define column matrix vtrunc and fill it with ttall values of accepted terms, then print it
  TMatrixD vtrunc(atrunc.size(), 1);
  for (unsigned int jv = 0; jv < atrunc.size(); jv++) vtrunc(jv, 0) = atrunc[jv];
  //define row matrix tvtrunc as the transpose of column matrix vtrunc
  TMatrixD tvtrunc(TMatrixD::kTransposed, vtrunc); 
  //define vtrunctot as sum of elements of vtrunc
  double vtrunctot = 0;
  for (int jt = 0; jt < vtrunc.GetNrows(); jt++) vtrunctot += vtrunc(jt, 0);
  //define vector covs of covariance matrices and fill it with matrix cov
  vector<TMatrix> covs;
  covs.push_back(cov);
  //check if there are rejected terms and remove correspondant row and column
  if (tails.size()) {
    for (int jta = tails.size() - 1; jta >= 0; jta--) {   
      TMatrix covtemp = CleanCovariance(covs[covs.size() - 1], tails[jta]); 
      covs.push_back(covtemp); } }

  //define matrix covcut which is final covariance matrix after cut
  TMatrix covcut = covs[covs.size() - 1];
  if (!covcut.GetNrows()) return -999;
  //define matrix invcovmod as the copy of matrix covcut
  TMatrixD invcovmod = covcut;
  //define matrix vtcovmod as product of matrix invcovmod and column matrix vtrunc, then print it
  TMatrixD vtcovmod = invcovmod * vtrunc;

  //initialize vmediomod and vector tterms
  double vmediomod = 0;
  vector<double> tterms;
  //iterate over elements of row matrix tvtrunc and column matrix vtcovmod
  for (unsigned int jv = 0; jv < ttrunc.size(); jv++) {
    double tterm = tvtrunc(0, jv) * vtcovmod(jv, 0);
    tterms.push_back(tterm);
    vmediomod += tterm; }
    
  //define column matrix mmterms as matrix product of tvtrunc and vtcovmod
  TMatrixD mtterms(tterms.size(), 1);
  for (unsigned int jv = 0; jv < tterms.size(); jv++) mtterms(jv, 0) = tterms[jv]; 
  cout << "column matrix mtterms, with no trunc c2 function terms = " << endl; mtterms.Print();
  //compute c2 function with trunc as average of mtterms
  vmediomod /= (ttrunc.size());
  return vmediomod; 
}

//fill linear matrix with off-diag MCS terms
const void TrajectoryMCSFitterICARUS::FillOffDiagMCSLin(TMatrixDSym& matOD, TMatrixD matdiag, int jp) const {
  double whp = 0.174;
  if (jp >= 1) {
    matOD(jp, jp-1) += sqrt(matdiag(jp, jp) * matdiag(jp-1, jp-1)) * whp; 
    matOD(jp-1, jp) += sqrt(matdiag(jp, jp) * matdiag(jp-1, jp-1)) * whp; } 
}

//fill polygonal matrix with off-diag MCS terms
const void TrajectoryMCSFitterICARUS::FillOffDiagMCSPoly(TMatrixDSym& matOD, TMatrixD matdiag, int jp) const {
  double whp = 0.174;
  double whpp = 0.067;
  if (jp >= 1) {
    matOD(jp, jp-1) += sqrt(matdiag(jp, jp) * matdiag(jp-1, jp-1)) * whp;
    matOD(jp-1, jp) += sqrt(matdiag(jp, jp) * matdiag(jp-1, jp-1)) * whp; }
  if (jp >= 2) { 
    matOD(jp, jp-2) += sqrt(matdiag(jp, jp) * matdiag(jp-2, jp-2)) * whpp;   
    matOD(jp-2, jp) += sqrt(matdiag(jp, jp) * matdiag(jp-2, jp-2)) * whpp; } 
}

//fill linear matrix with off-diag error terms
const void TrajectoryMCSFitterICARUS::FillOffDiagErrLin(TMatrixDSym& matOD, TMatrixD matdiag, int jp) const {
  double whp = -1./2.;
  if (jp >= 1) {
    matOD(jp, jp-1) += matdiag(jp, jp) * whp;
    matOD(jp-1, jp) += matdiag(jp, jp) * whp; } 
}

//fill polygonal matrix with off-diag error terms
const void TrajectoryMCSFitterICARUS::FillOffDiagErrPoly(TMatrixDSym& matOD, TMatrixD matdiag, int jp) const {
  double whp = -2./3.;
  double whpp = 1./6.;
  if (jp >= 1) {
    matOD(jp, jp-1) += matdiag(jp, jp) * whp;
    matOD(jp-1, jp) += matdiag(jp, jp) * whp; }
  if (jp >= 2) {
    matOD(jp, jp-2) += matdiag(jp, jp) * whpp;
    matOD(jp-2, jp) += matdiag(jp, jp) * whpp; } 
}

//fill full matrix with off-diag mixed terms
const void TrajectoryMCSFitterICARUS::FillMCSMixedTerms(TMatrixD& matMix, TMatrixD matLin, TMatrixD matPoly) const {
  double whp = 0.174;
  double whpp = 0.067;
  for (int jp = 0; jp < matPoly.GetNrows(); jp++) {
    matMix(jp, jp) += sqrt(matLin(jp, jp) * matPoly(jp, jp)) * whp;
    if (jp >= 1) matMix(jp-1, jp) += sqrt(matLin(jp, jp) * matPoly(jp-1, jp-1)) * whpp;
    if (jp < matPoly.GetNrows() - 1) matMix(jp+1, jp) += sqrt(matLin(jp, jp) * matPoly(jp+1, jp+1)) * whp; 
    if (jp < matPoly.GetNrows() - 2) matMix(jp, jp+2) += sqrt(matLin(jp, jp) * matPoly(jp+2, jp+2)) * whpp; } 
}

//perform fit of c2 function to determine MCS momentum
const TrajectoryMCSFitterICARUS::ScanResult TrajectoryMCSFitterICARUS::C2Fit(const recob::TrackTrajectory& traj, std::vector<size_t>& breakpoints, std::vector<float>& seglens, std::vector<float>& cumseglens, std::vector<int>& seghits, std::vector<int>& cumseghits, std::vector<float> dthetaLin, std::vector<float> dthetaPoly, int pid, double sigma) const {
  //define number of steps
  int nMom = ceil((pMax_ - pMin_) / pStep_) + 1;
  //define vector momentum, c2 function, error on momentum, error on c2 function 
  vector<float> wmom(nMom);
  vector<float> wc2(nMom);
  vector<float> wsmom(nMom);
  vector<float> wsigma(nMom);

  //iterate over p_test, compute c2 function and populate vectors wmom, wc2, wsmom, wsigma
  int jMom = 0; int firstValid = 0; bool allZero = true;
  for (double p_test = pMin_; p_test <= pMax_; p_test += pStep_) {
    cout << "test momentum [MeV/c] = " << p_test << endl; cout << " " << endl;
    double c2 = C2Function(traj, breakpoints, seglens, cumseglens, seghits, cumseghits, dthetaLin, dthetaPoly, pid, p_test);
    cout << "c2 function with trunc = " << c2 << endl; cout << " " << endl;
    wmom[jMom] = p_test / 1000.; 
    wc2[jMom] = c2; 
    wsmom[jMom] = 1;
    wsigma[jMom] = sigma * c2; 
    if (firstValid < 0.01 && c2 > 0.01) firstValid = jMom;
    if (c2 > 0.) allZero = false;
    jMom++; }
  cout << " " << endl;

  //return null is c2 function is null for all values of p test
  if (allZero) {
    cout << "c2 function is zero everywhere! end fit" << endl;
    return ScanResult(); }

  //define new vectors for momentum, c2 function, error on momentum, error on c2 function only for c2!=0
  vector<float> rmom(nMom - firstValid);
  vector<float> rc2(nMom - firstValid);
  vector<float> rsmom(nMom - firstValid);
  vector<float> rsigma(nMom - firstValid);
  for (int jp = firstValid; jp < jMom; jp++) {
    rmom[jp - firstValid] = wmom[jp];
    rc2[jp - firstValid] = wc2[jp];
    rsmom[jp - firstValid] = wsmom[jp];
    rsigma[jp - firstValid] = wsigma[jp]; }
  cout << "momentum [MeV/c] = "; for (auto i : rmom) cout << i << ' '; cout << endl;
  cout << "c2 function = "; for (auto i : rc2) cout << i << ' '; cout << endl;
  cout << "momentum error [MeV/c] = "; for (auto i : rsmom) cout << i << ' '; cout << endl;
  cout << "c2 function error = "; for (auto i : rsigma) cout << i << ' '; cout << endl;
  
  //define TVectorD for momentum, c2 function, error on momentum, error on c2 function 
  TVectorD cmom(rmom.size());
  for (size_t i = 0; i < rmom.size(); ++i) cmom[i] = rmom[i];
  TVectorD cc2(rc2.size());
  for (size_t i = 0; i < rc2.size(); ++i) cc2[i] = rc2[i];
  TVectorD csmom(rsmom.size());
  for (size_t i = 0; i < rsmom.size(); ++i) csmom[i] = rsmom[i];
  TVectorD csigma(rsigma.size());
  for (size_t i = 0; i < rsigma.size(); ++i) csigma[i] = rsigma[i];

  //root class for graph with x = cmom, y = cc2, dx = csmom, dy = csigma
  TGraphErrors *gr3 = new TGraphErrors(cmom, cc2, csmom, csigma);
  //define minimum momentum, maximum momentum, number of parameters of fit function
  double p_min = wmom[firstValid]; double p_max = pMax_ / 1000.; int npar = 2;
  cout << "p_min [GeV/c] = " << p_min << " p_max [GeV/c] = " << p_max << endl;
  //root class for fit to function 1 / (par[0] + par[1] / (x[0] * x[0]))
  TF1* fitfunc = new TF1("fitfunc", funzio, p_min, p_max, npar);
  //set limits to parameters alpha (0) in -0.5,1 and beta (1) in 0,5
  fitfunc->SetParLimits(0, -0.5, 1.);
  fitfunc->SetParLimits(1, -0., 5.);
  //execute fit using root class defined above
  gr3->Fit("fitfunc"); 

  //find optimal values of function parameters alpha, beta with associated errors
  double alpha, dalpha, beta, dbeta;
  alpha = fitfunc->GetParameter(0); dalpha = fitfunc->GetParError(0);
  beta = fitfunc->GetParameter(1); dbeta = fitfunc->GetParError(1); 

  //remove root objects used above for saving memory
  gr3->Delete(); fitfunc->Delete();
  //compute best p and error p 
  double p_best = sqrt(beta / (1 - alpha));
  double p_error = sqrt(pow(dbeta, 2.) / (4 * beta * (1 - alpha)) + beta * pow(dalpha, 2.) / (4 * pow((1 - alpha), 3.)));

  //return best p, error p 
  cout << "alpha = " << alpha << " and error on alpha = " << dalpha << endl;
  cout << "beta = " << beta << " and error on beta = " << dbeta << endl;
  cout << "best momentum = " << p_best << " and error on momentum = " << p_error << endl;
  if (p_best < p_min) cout << "best momentum " << p_best << " less than p_min " << p_min << endl;
  else if (p_best > p_max) cout << "best momentum " << p_best << " greater than p_max " << p_max << endl;
  return ScanResult(
    p_best, p_error, p_min, p_max,
    alpha, dalpha, beta, dbeta,
    rmom, rc2); }

//define drift origin coordinate [cm]
void TrajectoryMCSFitterICARUS::AnodeDistance(int cryo, int tpc, double x0) const { 
  //cryostat EAST (0), TPC EES (0) and EEN (1)
  if (cryo == 0 && tpc <= 1) x0 = -359.6;
  //cryostat EAST (0), TPC EWS (2) and EWN (3)
  if (cryo == 0 && tpc > 1) x0 = -60.8;
  //cryostat WEST (1), TPC WES (0) and WEN (1)
  if (cryo == 1 && tpc <= 1) x0 = 60.8;
  //cryostat WEST (1), TPC WWS (2) and WWN (3)
  if (cryo == 1 && tpc > 1) x0 = 359.6; 
}

//define rotation matrix depending on plane, tpc and cryostat (not used)
TMatrixD TrajectoryMCSFitterICARUS::ReferenceFrame(int plane) const { 
  //define rotation matrix 3x3 
  TMatrixD mat(3, 3);
  //populate rotation matrix for view 0 (induction 1)
  if (plane == 0) {
    mat(0, 0) = 1; mat(0, 1) = 0; mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = 1; mat(1, 2) = 0;
    mat(2, 0) = 0; mat(2, 1) = 0; mat(2, 2) = 1; }
  //populate rotation matrix for view 1 (induction 2)
  if (plane == 1) {
    mat(0, 0) = 1; mat(0, 1) = 0;              mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = cos(M_PI / 3);  mat(1, 2) = sin(M_PI / 3);
    mat(2, 0) = 0; mat(2, 1) = -sin(M_PI / 3); mat(2, 2) = cos(M_PI / 3); }
  //populate rotation matrix for view 2 (collection)
  if (plane == 2) {
    mat(0, 0) = 1; mat(0, 1) = 0;             mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = cos(M_PI / 3); mat(1, 2) = -sin(M_PI / 3);
    mat(2, 0) = 0; mat(2, 1) = sin(M_PI / 3); mat(2, 2) = cos(M_PI / 3); }
  return mat; 
}

//clean covariance matrix of jtail row and column
TMatrix TrajectoryMCSFitterICARUS::CleanCovariance(TMatrix cov, int jtail) const {
  //define size of input matrix cov
  int ncov = cov.GetNrows();
  //define new matrix covmod with reduced size
  TMatrix covmod(ncov - 1, ncov - 1);

  //fill matrix covmod with elements of matrix cov apart from row and column jtail
  for (int jrow = 0; jrow < ncov; jrow++) {
    for (int jcol = 0; jcol < ncov; jcol++) {
      if (jrow < jtail && jcol < jtail) covmod(jrow, jcol) = cov(jrow, jcol);
      if (jrow > jtail && jcol > jtail) covmod(jrow - 1, jcol - 1) = cov(jrow, jcol); } }

  //return modified matrix covmod without original row and column jtail
  return covmod; 
}

//check if some point is in plane = p
bool TrajectoryMCSFitterICARUS::isinplane(size_t index, unsigned int p) const {
  proxy::TrackPointData pd = pdata[index];
  art::Ptr<recob::Hit> hit = get<1>(pd);
  unsigned int plane = hit->WireID().Plane;
  if (plane == p) return true;
  else return false; 
}

//check if some point is in tpc = t
bool TrajectoryMCSFitterICARUS::isintpc(size_t index, unsigned int t) const {
  proxy::TrackPointData pd = pdata[index];
  art::Ptr<recob::Hit> hit = get<1>(pd);
  unsigned int tpc = hit->WireID().TPC;
  if (tpc == t) return true;
  else return false; 
}

//check tpc of last valid point
unsigned int TrajectoryMCSFitterICARUS::lasttpc(const recob::TrackTrajectory& traj) const {
  size_t index = traj.LastValidPoint();
  proxy::TrackPointData pd = pdata[index];
  art::Ptr<recob::Hit> hit = get<1>(pd);
  unsigned int tpc = hit->WireID().TPC;
  return tpc; 
}

//return 2d hit coordinates in plane = p and tpc = t
Vector_t TrajectoryMCSFitterICARUS::hit2d(const recob::TrackTrajectory& traj, size_t index, unsigned int p, unsigned int t) const {
  proxy::TrackPointData pd = pdata[index];
  art::Ptr<recob::Hit> hit = get<1>(pd);
  unsigned int plane = hit->WireID().Plane;
  unsigned int tpc = hit->WireID().TPC;
  if (plane == p && tpc == t) {

    /* Petrillo functions
    auto const& geom = *(lar::providerFrom<geo::Geometry>());
    geo::PlaneGeo const& geoplane = geom.Plane(hit->WireID());
    geo::Point_t myPoint = traj.LocationAtPoint(index);
    auto [wireCoord, driftDistance] = geoplane.DecomposePoint(myPoint);
    double x = driftDistance; 
    double y = wireCoord.Y(); 
    double z = wireCoord.X(); 
    */

    double x = hit->PeakTime() * 0.0622; 
    double y = hit->WireID().Wire * 0.3; 
    return Vector_t(x, y, 0); }
  else return Vector_t(); 
}

//return integral propto dedx
double TrajectoryMCSFitterICARUS::Integral(const recob::TrackTrajectory& traj, size_t index) const {
  proxy::TrackPointData pd = pdata[index];
  art::Ptr<recob::Hit> hit = get<1>(pd);
  return hit->Integral();
}

//find drift length of track [cm]
double TrajectoryMCSFitterICARUS::length1D(const recob::TrackTrajectory traj, unsigned int plane) const {
  //define index of first and last valid point
  size_t index = traj.FirstValidPoint();
  auto frstIndex = -1; 
  auto lastIndex = -1; 

  //in case of hits 2d, find drift length of track in chosen view and last tpc
  if (dimMode_ == 2) {
    unsigned int tpc = lasttpc(traj);
    //find index of first and last valid point
    while (index < traj.LastValidPoint()) {
      if (isinplane(index, plane) && isintpc(index, tpc)) {
        if (frstIndex == -1) frstIndex = index;
        lastIndex = index; }
      index = traj.NextValidPoint(index + 1); }
    //define first and last valid point
    Vector_t frstPoint = hit2d(traj, frstIndex, plane, tpc);
    Vector_t lastPoint = hit2d(traj, lastIndex, plane, tpc);
    //compute drift coordinate difference
    double dx = abs(frstPoint.X() - lastPoint.X());
    //return drift length 
    return dx; }
  
  //if not the case of hits 2d, find length of track along x direction
  else {
    //find index of first and last valid point
    frstIndex = traj.FirstValidPoint();
    lastIndex = traj.LastValidPoint();
    //define first and last valid point
    geo::Point_t frstPoint = traj.LocationAtPoint(frstIndex);
    geo::Point_t lastPoint = traj.LocationAtPoint(lastIndex);
    //compute x coordinate difference
    double dx = abs(frstPoint.X() - lastPoint.X());
    //return x length
    return dx; } 
}
  
//find drift-wire length of track [cm]
double TrajectoryMCSFitterICARUS::length2D(const recob::TrackTrajectory traj, unsigned int plane) const {
  //define index of first and last valid point
  size_t index = traj.FirstValidPoint();
  auto frstIndex = -1; 
  auto lastIndex = -1; 

  //in case of hits 2d, find drift-wire length of track in chosen view and last tpc
  if (dimMode_ == 2) {
    unsigned int tpc = lasttpc(traj);
    //find index of first and last valid point 
    while (index < traj.LastValidPoint()) {
      if (isinplane(index, plane) && isintpc(index, tpc)) {
        if (frstIndex == -1) frstIndex = index;
        lastIndex = index; }
      index = traj.NextValidPoint(index + 1); }
    //define first and last valid point
    Vector_t frstPoint = hit2d(traj, frstIndex, plane, tpc);
    Vector_t lastPoint = hit2d(traj, lastIndex, plane, tpc);
    //compute drift and wire coordinate differences
    double dx = abs(frstPoint.X() - lastPoint.X());
    double dy = abs(frstPoint.Y() - lastPoint.Y());
    //return drift-wire length 
    return sqrt(pow(dx, 2) + pow(dy, 2)); }
  
  //if not the case of hits 2d, find length of track in plane xy
  else {
    //find index of first and last valid point
    frstIndex = traj.FirstValidPoint();
    lastIndex = traj.LastValidPoint();
    //define first and last valid point
    geo::Point_t frstPoint = traj.LocationAtPoint(frstIndex);
    geo::Point_t lastPoint = traj.LocationAtPoint(lastIndex);
    //compute x and y coordinate differences
    double dx = abs(frstPoint.X() - lastPoint.X());
    double dy = abs(frstPoint.Y() - lastPoint.Y());
    //return xy length 
    return sqrt(pow(dx, 2) + pow(dy, 2)); } 
}
  
//find 3D length of track [cm]
double TrajectoryMCSFitterICARUS::length3D(const recob::TrackTrajectory traj, unsigned int plane) const {
  //define index of first and last valid point
  size_t index = traj.FirstValidPoint();
  auto frstIndex = -1; 
  auto lastIndex = -1; 

  //in case of hits 2d, find 3D length of track in chosen view and last tpc
  if (dimMode_ == 2) {
    unsigned int tpc = lasttpc(traj);
    //find index of first and last valid point
    while (index < traj.LastValidPoint()) {
      if (isinplane(index, plane) && isintpc(index, tpc)) {
        if (frstIndex == -1) frstIndex = index;
        lastIndex = index; }
      index = traj.NextValidPoint(index + 1); }
    //define first and last valid point 
    geo::Point_t frstPoint = traj.LocationAtPoint(frstIndex);
    geo::Point_t lastPoint = traj.LocationAtPoint(lastIndex);
    //compute x, y, z coordinate differences
    double dx = abs(frstPoint.X() - lastPoint.X());
    double dy = abs(frstPoint.Y() - lastPoint.Y());
    double dz = abs(frstPoint.Z() - lastPoint.Z());
    //return 3D length 
    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)); }
  
  //if not the case of hits 2d, find 3D length of track
  else {
    //find index of first and last valid point
    frstIndex = traj.FirstValidPoint();
    lastIndex = traj.LastValidPoint();
    //define first and last valid point 3d
    geo::Point_t frstPoint = traj.LocationAtPoint(frstIndex);
    geo::Point_t lastPoint = traj.LocationAtPoint(lastIndex);
    //compute x and y coordinate differences
    double dx = abs(frstPoint.X() - lastPoint.X());
    double dy = abs(frstPoint.Y() - lastPoint.Y());
    double dz = abs(frstPoint.Z() - lastPoint.Z());
    //return 3D length
    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)); } 
}

//compute directional factor for thetamcs in c2 function computation
double TrajectoryMCSFitterICARUS::ThetaMCSFactor(const recob::TrackTrajectory& traj, unsigned int plane) const {
  //recall length2D as length of track on chosen plane (I1, I2, C)
  double L2D = length2D(traj, plane);
  //recall length3D as length of track in three dimensions
  double L3D = length3D(traj, plane);

  //in case of hits 2d return factor 2d for thetamcs
  if (dimMode_ == 2) return L3D / L2D;
  //if not the case of 2d hits return factor 3d for thetamcs
  else return 1;
}

//compute expected angle (thetamcs) for a certain momentum
double TrajectoryMCSFitterICARUS::ThetaExpected(const recob::TrackTrajectory& traj, double p, double beta, double L) const {
  //define Highland parameters S2 [MeV], epsilon, X0 [cm]
  double S2 = 13.6; double eps = 0.038; double X0 = 14.;

  //compute theta MCS [mrad] just as claimed by Highland formula
  double thetamcs = S2 / (p * beta) * sqrt(L/X0) * (1 + eps * log(L/X0)) * ThetaMCSFactor(traj, planeMode_);
  return thetamcs;
}

//compute directional factor for thetaerr in c2 function computation
double TrajectoryMCSFitterICARUS::ThetaErrFactor(const recob::TrackTrajectory& traj, unsigned int plane) const {
  //recall length1D as length of track along drift direction
  double L1D = length1D(traj, plane);
  //recall length2D as length of track on chosen plane (I1, I2, C)
  double L2D = length2D(traj, plane);
  //recall length3D as length of track in three dimensions
  double L3D = length3D(traj, plane);

  //in case of hits 2d return factor 2d for thetaerr
  if (dimMode_ == 2) return sqrt(1 - pow(L1D, 2)/pow(L2D, 2));
  //if not the case of 2d hits return factor 3d for thetaerr
  else return sqrt(1 - pow(L1D, 2)/pow(L3D, 2));
}

//print delta3p ready for ntuples
double TrajectoryMCSFitterICARUS::PrintD3P() const {
  //use number of hits in different views as weight 
  vector<float> weight;
  weight.push_back(float(hits2dI1.size())); 
  weight.push_back(float(hits2dI2.size())); 
  weight.push_back(float(hits2dC.size()));
  
  //recall delta3p as computed in different views
  vector<float> delta3ps; 
  delta3ps.push_back(0.1 * d3pI1); 
  delta3ps.push_back(0.1 * d3pI2); 
  delta3ps.push_back(0.1 * d3pC);
  
  //compute weight average of delta3p
  float delta3pw = inner_product(delta3ps.begin(), delta3ps.end(), weight.begin(), 0.0f);
  delta3pw = delta3pw / accumulate(weight.begin(), weight.end(), 0.0f);
  
  //return and print delta3p for correct input parameters
  float delta3p;
  if (dimMode_ == 2) delta3p = delta3ps[planeMode_];
  else delta3p = delta3pw;
  cout << "delta3p [cm] = " << delta3p << endl; cout << " " << endl;
  return delta3p; 
}

//check if angles vector contains a minimum number of angles > 0 to avoid crash
void TrajectoryMCSFitterICARUS::ThetaCheck(vector<float> dthetaLin, vector<float> dthetaPoly, bool& checkLin, bool& checkPoly, unsigned int& firstseg, unsigned int& lastseg) const {
  bool first = false;
  auto countLin = 0;
  for (unsigned int j = 0; j < dthetaLin.size(); ++j) {
    if (dthetaLin[j] > 0.) {
      if (!first) {
        firstseg = j; first = true; }
      lastseg = j + 1;
      countLin++;
      if (countLin >= minNAngs_) checkLin = true; } 
    else {
      first = false;
      countLin = 0; } }

  auto countPoly = 0;
  for (unsigned int j = 0; j < dthetaPoly.size(); ++j) {
    if (dthetaPoly[j] > 0.) {
      countPoly++;
      if (countPoly >= minNAngs_ - 1) checkPoly = true; } 
    else countPoly = 0; } 
}

//geometrical check if track stops inside detector (used but info not added to ntuples)
bool TrajectoryMCSFitterICARUS::GeoStopCheck(const recob::TrackTrajectory& traj) const {
  size_t lastIndex = traj.LastValidPoint();
  geo::Point_t lastPoint = traj.LocationAtPoint(lastIndex); 
  double step = 20;

  double low_x; double high_x; 
  double last_x = lastPoint.X();
  if (last_x > 0) { low_x = 61.7; high_x = 358.73; }
  if (last_x < 0) { low_x = -358.73; high_x = -61.7; }
  cout << "last point x = " << last_x << " low x = " << low_x << " high x = " << high_x << endl;
  bool check_x = true; 
  if ((abs(last_x - low_x) < step) || (abs(last_x - high_x) < step)) {
    check_x = false;
    cout << "too near to borders in x direction!" << endl; }

  double low_y = -181.86; double high_y = 134.36; 
  double last_y = lastPoint.Y();
  cout << "last point y = " << last_y << " low y = " << low_y << " high y = " << high_y << endl;
  bool check_y = true; 
  if ((abs(last_y - low_y) < step) || (abs(last_y - high_y) < step)) {
    check_y = false;
    cout << "too near to borders in y direction!" << endl; }
  
  double low_z = -894.951; double high_z = 894.951; 
  double last_z = lastPoint.Z();
  cout << "last point z = " << last_z << " low z = " << low_z << " high z = " << high_z << endl;
  bool check_z = true; 
  if ((abs(last_z - low_z) < step) || (abs(last_z - high_z) < step)) {
    check_z = false;
    cout << "too near to borders in z direction!" << endl; }

  if (check_x && check_y && check_z) {
    cout << "stopping track found!" << endl;
    return true; }
  else {
    cout << "crossing track found!" << endl; 
    return false; }
}

//manual check to exclude hit belonging to specific indeces to avoid delta rays
bool TrajectoryMCSFitterICARUS::DeltaCheck(size_t index) const {
  unordered_set<size_t> excluded_indices = {0};
  if (excluded_indices.find(index) == excluded_indices.end()) return true;
  else return true;
}

//define cathode drift coordinates [cm]
double TrajectoryMCSFitterICARUS::CathodeDistance(unsigned int c) const { 
  //cryostat EAST (0)
  if (c == 0) return -210.2;
  //cryostat WEST (1)
  if (c == 1) return 210.2; 
  else return 0; 
}

//check to exclude hit near cathode
bool TrajectoryMCSFitterICARUS::CathodeCheck(const recob::TrackTrajectory& traj, size_t index) const {
  proxy::TrackPointData pd = pdata[index];
  art::Ptr<recob::Hit> hit = get<1>(pd);
  unsigned int c = hit->WireID().Cryostat;
  double cathode = CathodeDistance(c);
  double point = traj.LocationAtPoint(index).X();
  double tolerance = 15;
  if (point < cathode - tolerance || point > cathode + tolerance) return true;
  else return false;
}

vector<float> TrajectoryMCSFitterICARUS::CleanThetaLin(vector<float> dthetaLin) const {
  float threshold=0.1; //rad
 vector<float> dthetaLinClean;
  for (unsigned int j = 0; j < dthetaLin.size(); ++j) {
    if (dthetaLin[j] > threshold) {
      std::cout << "cleaning dthetaLin[j] = " << dthetaLin[j] << std::endl;
      dthetaLinClean.push_back(dthetaLin[j]);
}
  }
  return dthetaLinClean;
}
vector<float> TrajectoryMCSFitterICARUS::CleanThetaPoly(vector<float> dthetaPoly) const {
  float threshold=0.1;
 vector<float> dthetaPolyClean;
  for (unsigned int j = 0; j < dthetaPoly.size(); ++j) {
    std::cout << "cleaning dthetaPoly[j] = " << dthetaPoly[j] << std::endl;

    if (dthetaPoly[j] > threshold) {
      dthetaPolyClean.push_back(dthetaPoly[j]);
}
  }
return dthetaPolyClean;
}

/*--------------------------------------------------------------------------------------------*/
void TrajectoryMCSFitterICARUS::ProcessDeltaRays(const recob::TrackTrajectory& traj, int viewType, std::vector<int>& isDelta) const
{
    /// This is the algorithm to "tag" hits as delta ray candidates.
  /// The idea is to loop over all hits of the selected track
  /// and fill the hit field  "hit->isdRay" (BasicHits.h) with a value
  /// = 0   if the hit is NOT a delta ray candidate
  /// > 0   if the hit is     a delta ray candidate
  /*view_type = viewType; 
  hits = (&track);
  std::sort(hits->begin(),hits->end(),bWireLess());
  nHits = hits->size()*/

  //Create_ZY_arrays();
  //Fill_ZY_arrays();
  std::cout << "before tagging delta rays " << std::endl;
  TagDeltaRays(traj,viewType, isDelta);
  //TagCloseToDeltaRays();
  //TagOverlappingDeltaRays();

//  Delete_ZY_arrays();

  return ;
}
/*--------------------------------------------------------------------------------------------*/
void TrajectoryMCSFitterICARUS::Create_ZY_arrays()
{
}
/*--------------------------------------------------------------------------------------------*/
void TrajectoryMCSFitterICARUS::Fill_ZY_arrays(std::vector<double>& Z, std::vector<double>& Y)
{
  /*


  const float driftvelocity = 1.57;
  const float PitchCol=3.;     // Collection pitch fixed at 3 mm (irrelevant)

  AF::Hit2D* hit; 
  int NHits = nHits;
  for (int i=0;i<NHits;i++)
    {
      hit = (AF::Hit2D*) hits->at(i);
      
      Z[i] = hit->iWire*PitchCol;        // Wire coordinate [mm]
      Y[i] = hit->iDrift*driftvelocity;  // Drift time sample coor. [mm]
    }
  return;
  */
}
/*--------------------------------------------------------------------------------------------*/
void TrajectoryMCSFitterICARUS::TagDeltaRaysLocal()
{
  /// This routine tags hits which are delta-ray candidates.
  /// We work with hits ORDERED by increasing wire number.
  /// The tag is "hit->isdRay" and can have several values:
  ///  = 0  the hit is NO delta ray candidate
  ///  > 0  the hit is    delta ray candidate:
  ///  = 1  the hit is in a wire with more than one hit
  ///       and it is NOT the closest one to the track direction
  ///  = 2  the hit is far from next hit (good to tag short drays)
  ///  = 3  the hit is in the vecinity of a delta ray

  ///  At the end of this routine, DrDist if filled with:
  ///  a) distance from hit to track direction, for dray candidates
  ///  b) distance from hit to next hit, for non dray candidates
/*
  AF::Hit2D* hit;
  AF::Hit2D* hitf; 
  AF::Hit2D* hit2;
  int nHitsBack=0, nMaxBack=70;
  double x[nMaxBack],y[nMaxBack];
  int ifail, iWireprev=0;
  float xdr, ydr, distance;
  double sd, a0, a1;
  int imayor, nmd=20, ini, ifi;

  int NHits = nHits;  // number of hits linked to the track

  for (int i=0;i<NHits;i++)
    {
      hit = (AF::Hit2D*) hits->at(i);;  // hits ORDERED !! (by wire number)
      hit->isdRay=0;      // is not a delta-ray hit (default)
      //next line not only check if hit is empty, but also prevent from finding delta rays
      // in case the hit fitting has not been performed, which is not what I whant.
      // let's comment out the line for the moment
      //if(hit->Adctot <= 1.e-10) continue; // skip "empty" hits
      if (hit->iWire == iWireprev)  // wire with more than 1 hit -> dRay candidate
	{
	  hit->isdRay=1;  // d-ray candidate (flag also previous hit)
	  hitf = (AF::Hit2D*) hits->at(i-1);
	  hitf->isdRay = 1;
	}
      iWireprev = hit->iWire;
    }
  for (int i=0;i<NHits;i++)
    {
      hit = (AF::Hit2D*) hits->at(i);
      Dr[i] = hit->isdRay;
      DrDist[i]=0.;
    }

  for (int i=0;i<nMaxBack;i++){
    x[i]=0.;
    y[i]=0.;}

  for (int i=0;i<NHits;i++)
    {
      if(Dr[i])   // delta ray candidate found ...
	{
	  xdr=Z[i];
	  ydr=Y[i];
	  nHitsBack=0;
	  for (int j=(i-5);j>0;j--) // loop back ...
	    {
	      // take "no-delta ray" candidates only ...
	      if(Dr[j] == 0 && 
		 Y[j] != -999999.0 && Z[j] != -999999.0)
		{
		  x[nHitsBack]=Z[j];
		  y[nHitsBack]=Y[j];
		  nHitsBack++;
		}
	      if (nHitsBack==nMaxBack)
		break;
	    }
	  // Make a linear fit of the points
	  if (nHitsBack>3)
	    {
	      dlsqp1_(&nHitsBack, x, y, &a0, &a1, &sd, &ifail);
	      if(ifail==0){
		distance = distance_point_line(a0,a1,xdr,ydr);
		DrDist[i] = distance;
	      }
	    }
	}
    }
  for (int i=0;i<NHits;i++)
    {
      hit = (AF::Hit2D*) hits->at(i);
      if (Dr[i])
	{
	  imayor=0;
	  ini=(i+1)-nmd;
	  ifi=(i+1)+nmd;
	  if(ini<=0)       ini=1;
	  if(ifi>NHits)    ifi=NHits;
	  for(int j=ini;j<ifi;j++)
	    {
	      hit2 = (AF::Hit2D*) hits->at(j);
	      if(hit->iWire == hit2->iWire && i!=j)
		if(DrDist[i] > DrDist[j])
		  imayor = 1;
	    }
	  if(!imayor)
	    Dr[i] = 0;
	}
    }

  
  for (int i=0;i<(NHits-1);i++)         // loop forward
    if (!Dr[i] && Y[i]>0.)              // only for non d-ray candidates
      for(int k=(i+1);k<(NHits-1);k++)  // go forward searching for a non d-ray hit
	if(!Dr[k] && Y[k]>0.) {
	    DrDist[i]=sqrt((Z[k]-Z[i])*(Z[k]-Z[i])+(Y[k]-Y[i])*(Y[k]-Y[i]));
	    break; }
*/
  return;
}
/*--------------------------------------------------------------------------------------------*/
void TrajectoryMCSFitterICARUS::TagDeltaRays(const recob::TrackTrajectory& traj, int viewType, std::vector<int>& isDelta) const
{
  /// This routine tags hits which are delta-ray candidates.
  /// We work with hits ORDERED by increasing wire number.
  /// The tag is "hit->isdRay" and can have several values:
  ///  = 0  the hit is NO delta ray candidate
  ///  > 0  the hit is    delta ray candidate:
  ///  = 1  the hit is in a wire with more than one hit
  ///       and it is NOT the closest one to the track direction
  ///  = 2  the hit is far from next hit (good to tag short drays)
  ///  = 3  the hit is in the vecinity of a delta ray

  ///  At the end of this routine, DrDist if filled with:
  ///  a) distance from hit to track direction, for dray candidates
  ///  b) distance from hit to next hit, for non dray candidates

  std::vector<double> XP,YP; 
  std::vector<double> xp,yp; 
  TGraph* track;
  int ifail;
  float xdr, ydr, distance, mean, RMS;
  double a0, a1;
  int imenor, nmd=20, ini, ifi;
 std::vector<double> E,Dr,DrDist;
std::vector<recob::Hit> hits;
recob::Hit hit,hitf;
std::vector<int> isdRay;

  const float driftvelocity = 1.57;
  const float PitchCol=3.;     // Collection pitch fixed at 3 mm (irrelevant)
std::cout << " before filling vectors " << viewType << std::endl;


size_t index = traj.FirstValidPoint();
 while (index < traj.LastValidPoint()) {
    proxy::TrackPointData pd = pdata[index];
  art::Ptr<recob::Hit> hit = get<1>(pd);
    if (isinplane(index, viewType) && isintpc(index, lasttpc(traj))) {
      hits.push_back(*hit);
      YP.push_back(hit->WireID().Wire*PitchCol);        // Wire coordinate [mm]
      XP.push_back(hit->PeakTime()*driftvelocity);  // Drift time sample coor. [mm]
      }
    index = traj.NextValidPoint(index + 1); }
 
    
    int NHits=hits.size();
      for (int i=0;i<NHits;i++) 
      isdRay.push_back(0);          // it is not a delta-ray hit (default)

std::cout << " after filling vectors " << std::endl;
  
  for (int i=0;i<NHits;i++) {
      auto hit=hits.at(i);
      //next line not only check if hit is empty, but also prevent from finding delta rays
      // in case the hit fitting has not been performed, which is not what I whant.
      // let's comment out the line for the moment
      //if(hit->Adctot <= 1.e-10) continue; // skip "empty" hits
      std::cout << " hit " << i << " wire " << hit.WireID().Wire << std::endl;
      std::cout << " nhitsonwire " << HitsOnWire(hits,hit.WireID().Wire).size() << std::endl;
      std::vector<int> how=HitsOnWire(hits,hit.WireID().Wire);
      if (how.size()>1)  // wire with more than 1 hit -> dRay candidate
	{
    for(unsigned int jhow=0;jhow<how.size();jhow++) {
	  isdRay[how[jhow]]=1;  // d-ray candidate (flag also previous hit)
    if(viewType==2)  
    std::cout << " flagging delta i " << how[jhow] << " wire " << YP[how[jhow]]/PitchCol << " tick " << XP[how[jhow]]/driftvelocity << std::endl;
    }
	}
  
    }
    std::cout << " after filling how " << NHits << std::endl;
 
  for (int i=0;i<NHits;i++)
    {
      hit = hits.at(i);
      Dr.push_back(isdRay[i]);
      if(isdRay[i]) std::cout << " candidate dr " << i << std::endl;
      DrDist.push_back(0.);
      if(!Dr[i]){
	xp.push_back(XP[i]);
	yp.push_back(YP[i]);
      }
    }

if(yp.size()<2) return;

  track = new TGraph(yp.size(), &(yp[0]), &(xp[0]));
  ifail = track->Fit("pol1","0");
    std::cout << " after filling xp " << ifail <<  std::endl;

  if(ifail==0){
    TF1* pol1=(TF1*) track->GetFunction("pol1");
    double par[2];
    pol1->GetParameters(par);
    a0 = par[0];
    a1 = par[1];
    mean = 0.;
        std::cout << " before distance " << std::endl;

    for (int i=0;i<NHits;i++)
      {
	xdr=XP[i];
	ydr=YP[i];
      if(isdRay[i]) std::cout << " distance candidate dr " << i << std::endl;

	distance = distance_point_line(a0,a1,xdr,ydr);
	mean+=distance;
	DrDist[i] = distance;
      }
    mean /= (float)NHits;

    RMS = 0.;
    for (int i=0;i<NHits;i++)
    {
	RMS += (DrDist[i]-mean)*(DrDist[i]-mean);
    }
    RMS /= (float)NHits;
    RMS = sqrt(RMS);
  }
  else {
    std::cout<<"AF::FindDeltaRays::TagDeltaRays Error: cluster cannot be fitted to a straigh line. Delta Rays cannot be found"<<std::endl;
    return;
  }
  
  delete track;
  

  //temporarely tag with "2" single hits far from track in order to be distiguished from 
  // hits tagged with "1"

/*
  if(noisy[view_type-1]==1){
    plus=noisecut[view_type-1];
    plus *= RMS;
    //This step can be done only when the track (typically in Ind2) is immersed in a noisy region, unless you missidentify dRays..
    for (int i=0;i<NHits;i++){
       hit = hits.at(i);
//       if(view_type==3 && hit->iWire==2405)
// 	cout<<"0) "<<i<<" "<<Dr[i]<<" "<<DrDist[i]<<endl;

      if(!Dr[i] && DrDist[i] > mean + plus)
	Dr[i] = 2;

//       if(view_type==3 && hit->iWire==2405)
// 	    cout<<"1) "<<i<<" "<<Dr[i]<<" "<<DrDist[i]<<endl;
    }
  }
*/
    std::cout << " after deleting track " << std::endl;

    for (int i=0;i<NHits;i++)
    {
      hit = hits.at(i);
      if (Dr[i]==1)
	{
	  imenor=i;
	  ini=(i+1)-nmd;
	  ifi=(i+1)+nmd;
	  if(ini<=0)       ini=1;
	  if(ifi>NHits)    ifi=NHits;
	  for(int j=ini;j<ifi;j++)
	    {
	    auto  hit2 = hits.at(j);
	      if(hit.WireID().Wire == hit2.WireID().Wire && i!=j){
		if(DrDist[j] < DrDist[imenor])
		  imenor = j;
	      }
	    }
/*
	  if(noisy[view_type-1]==1){
	    if(DrDist[imenor] <= mean + plus) 
	      Dr[imenor] = 0;
	  }
	  else Dr[imenor] = 0;
    */
	}
    }
    std::cout << " after imenor " << std::endl;

  // re-tag single hits 2-->1
  for (int i=0;i<NHits;i++)
    if(Dr[i]==2)
      Dr[i] = 1;
  
  for (int i=0;i<(NHits-1);i++)         // loop forward
    if (!Dr[i] && XP[i]>0.)              // only for non d-ray candidates
      for(int k=(i+1);k<(NHits-1);k++)  // go forward searching for a non d-ray hit
	if(!Dr[k] && XP[k]>0.) {
	    DrDist[i]=sqrt((XP[k]-XP[i])*(XP[k]-XP[i])+(YP[k]-YP[i])*(YP[k]-YP[i]));
	    break; }

std::cout << " after drdist " << std::endl;

  for (int i=0;i<NHits;i++)
    {
      hit = hits.at(i);
      isdRay[i] = Dr[i];
      if(isdRay[i]) std::cout << " candidate delta ray " << i << std::endl;
      if(viewType==2&&isdRay[i] )  
      std::cout << " flagging final delta i " << i << " wire " << YP[i]/PitchCol << " tick " << XP[i]/driftvelocity << std::endl;
    }
    std::cout << " firstvalid " << traj.FirstValidPoint() << " lastvalid " << traj.LastValidPoint() << std::endl;
 index = traj.FirstValidPoint();
int count=0;
 while (index < traj.LastValidPoint()) {
    if (isinplane(index, viewType) && isintpc(index, lasttpc(traj))) {
      if(isdRay[count]) {isDelta[index]=1;
        std::cout << " flagging index " << index <<  std::endl;

      }
      count++;
      }
    index = traj.NextValidPoint(index + 1);
     }
 
    
  return;
  
}
/*--------------------------------------------------------------------------------------------*/
void TrajectoryMCSFitterICARUS::TagCloseToDeltaRays()
{
  /// This routine flags delta rays with tags nr.2 and 3

  /// It builds the distributions of the distances from each hit to the next hit
  /// "DrDist" for non already tagged d-ray hits
  /// (for Dr[i]=1 this is the distance to the straight line, but it does not matter)
  /// Compute the "mean" and the "RMS" and flag (Dr[i]=2) these hits that are far from
  /// the mean (one sigma) ... these hits are also d-ray candidates (very short) since
  /// they are far from the next hit.

  /// Finally, it flags (Dr[i]=3) also hits on the vecinity of the dray.
  /// All these hits are rejected when making segments for Mult.Scatt. analysis.

/*
  AF::Hit2D* hit;
  int auxi, auxj;
  int Nkeep, nflagback;
  float mean, RMS;

  int NHits = nHits;
  int *hIndex=NULL;

  hIndex = new int[NHits];

  for (int i=0;i<NHits;i++)
    hIndex[i] = i;

  // Order DrDist array in increasing order (by swaping)
  for (int i=0;i<(NHits-1);i++)
    for (int j=i+1;j<NHits;j++)
      if (DrDist[hIndex[i]] > DrDist[hIndex[j]])
	{
	  auxi = hIndex[i];
	  auxj = hIndex[j];
	  hIndex[i] = auxj;
	  hIndex[j] = auxi;
	}
  // Compute mean and RMS of the distributions of the distances
  Nkeep = 0;
  mean=0.;
  for (int i=0;i<NHits;i++)
    if (!Dr[hIndex[i]]){
      mean += DrDist[hIndex[i]];
      Nkeep++; }
  mean /= (float)Nkeep;
  RMS=0.;
  for (int i=0;i<NHits;i++)
    if (!Dr[hIndex[i]])
      RMS += (DrDist[hIndex[i]]-mean)*(DrDist[hIndex[i]]-mean);
  RMS /= (float)Nkeep;
  RMS = sqrt(RMS);

  // Flag this hit (Dr=2) if its distance to the next hit is higher than the mean
  // distance (taking all hits) + one sigma
  for (int i=0;i<NHits;i++)
    if(DrDist[hIndex[i]] > (mean+RMS))
      if (!Dr[hIndex[i]])               // if it is d-ray, it is already flaged
	Dr[hIndex[i]] = 2;
  
  // Flag also all the vecinity of the dray (the "nflagback" previous hits)
  // We reject all these hits when making segments for Mult.Scatt. analysis
  nflagback=3;
  for (int i=0;i<NHits;i++)
    if(Dr[i]==2)                          // hit close to d-ray
      for (int j=(i-1);j>=(i-nflagback);j--)
	if(!Dr[j])
	  Dr[j] = 3;


  for (int i=0;i<NHits;i++)
    {
      hit = (AF::Hit2D*) hits->at(i);
      hit->isdRay = Dr[i];
    }
  
  delete [] hIndex;
  return;
  */
}
/*--------------------------------------------------------------------------------------------*/
void TrajectoryMCSFitterICARUS::TagOverlappingDeltaRays()
{
  /*
  AF::Hit2D* hit;
  

  int NHits = nHits;  // number of hits linked to the track
  float mean, RMS, plus;
  float dE;
  int notDelta;
  std::vector<int> isdRayByCharge;
  std::vector<int> idxs;

  int wire,wireprev,wireprevprev,wirenext;

  mean = 0.;
  notDelta = 0;
  for (int j=0;j<NHits;j++)
    { 
      hit = (AF::Hit2D*) hits->at(j);  // hits ORDERED !! (by wire number)
      dE = view_type == T600::kViewColl? hit->GetEnergyDeposit() : hit->GetAmplitude();
      if(hit->isdRay==0) {
	isdRayByCharge.push_back(0);
	idxs.push_back(j);
	notDelta++;
	mean += dE; 
      }
    }

  //cout<<"idxs.size() "<<idxs.size()<<endl;
  if(idxs.size()>1){

    if(notDelta!=0){
      mean /= (float)notDelta;
    }

    RMS = 0.;
    for (int j=0;j<NHits;j++)
      {
	hit = (AF::Hit2D*) hits->at(j);  // hits ORDERED !! (by wire number)
	dE = view_type == T600::kViewColl ? hit->GetEnergyDeposit() : hit->GetAmplitude();
	if(hit->isdRay==0){
	  RMS += (dE-mean)*(dE-mean);
	}
      }
    if(notDelta!=0){
      RMS /= (float)notDelta;
    }
    RMS = sqrt(RMS);

    plus = threshold[view_type-1];
    plus *=RMS;


    for (int i=0;i<isdRayByCharge.size();i++)
      {
	hit = (AF::Hit2D*) hits->at(idxs[i]);  // hits ORDERED !! (by wire number)
	dE = view_type == T600::kViewColl? hit->GetEnergyDeposit() : hit->GetAmplitude();
      
	if(dE>(mean+plus)) isdRayByCharge[i]=5;
      }
  

    wire = ((AF::Hit2D*) hits->at(idxs[1]))->iWire;
    wireprev = ((AF::Hit2D*) hits->at(idxs[0]))->iWire;
    if(isdRayByCharge[1]==0 && isdRayByCharge[0]==5) {
      if(wire == wireprev+1)
	isdRayByCharge[0]=0;
    }

  
    for (int i=2;i<isdRayByCharge.size();i++)
      {
	wire = ((AF::Hit2D*) hits->at(idxs[i]))->iWire;
	wireprev = ((AF::Hit2D*) hits->at(idxs[i-1]))->iWire;
	wireprevprev = ((AF::Hit2D*) hits->at(idxs[i-2]))->iWire;
      
	if(isdRayByCharge[i]==0 && isdRayByCharge[i-1]==5 && isdRayByCharge[i-2]==0){
	  if((wire == wireprev+1) && (wireprev == wireprevprev+1)){
	    isdRayByCharge[i-1]=0;
	  }
	}
      }


    for (int i=1;i<isdRayByCharge.size()-1;i++)
      {
	wire = ((AF::Hit2D*) hits->at(idxs[i]))->iWire;
	wireprev = ((AF::Hit2D*) hits->at(idxs[i-1]))->iWire;
	wirenext = ((AF::Hit2D*) hits->at(idxs[i+1]))->iWire;
      
	if((isdRayByCharge[i]==5) && ((wire!=wireprev+1) && (wirenext!=wire+1)))
	  isdRayByCharge[i]=0;
	if((wire!=wireprev+1) && (wirenext==wire+1) && (isdRayByCharge[i]==5) && (isdRayByCharge[i+1]==0))
	  isdRayByCharge[i]=0;
	if((wire==wireprev+1) && (wirenext!=wire+1) && (isdRayByCharge[i]==5) && (isdRayByCharge[i-1]==0))
	  isdRayByCharge[i]=0;
      }


    wire = ((AF::Hit2D*) hits->at(idxs[isdRayByCharge.size()-1]))->iWire;
    wireprev = ((AF::Hit2D*) hits->at(idxs[isdRayByCharge.size()-2]))->iWire;
    if(isdRayByCharge[isdRayByCharge.size()-2]==0){
      if(wire == wireprev+1)
	isdRayByCharge[isdRayByCharge.size()-1]=0;
    }


    for (int i=0;i<isdRayByCharge.size();i++){
      hit = (AF::Hit2D*) hits->at(idxs[i]);
      if(hit->isdRay==0)hit->isdRay=isdRayByCharge[i];
    }

  }
  else{
    cout<<"AF::FindDeltaRays::TagOverlappingDeltaRays ---> No need to search for delta rays by Charge in this track: "<<endl;
    cout<<"all hits are tagged as outliers or track has only one hit"<<endl;
  }
  */
}
/*--------------------------------------------------------------------------------------------*/
std::vector<int> TrajectoryMCSFitterICARUS::HitsOnWire(std::vector<recob::Hit> hits,unsigned int iWire) const
// return the number of hits on a given wire
{
std::vector<int> how;

  for(unsigned int i=0;i<hits.size();i++)
    if(hits.at(i).WireID().Wire==iWire)
      how.push_back(i);
  return how;
}
/*--------------------------------------------------------------------------------------------*/
float TrajectoryMCSFitterICARUS::distance_point_line(const double a0, const double a1, const double xdr, const double ydr) const
{
/* Compute distance between a straight line "a1*x + a0" and a point (xdr,ydr)

   Inputs:
           a1 : slope                 of straight line.
           a0 : Y-axis crossing point of straight line.
           xdr,ydr : point coordinates.                                            
   Output:
           distance_point_line : distance point - straight line         */

  float distance;
  double vdir[2], vper[2], vmod, a1per, a0per, xcross, ycross;

  vdir[0]=1.;
  vdir[1]=atan(a1);
  vmod=sqrt(vdir[0]*vdir[0]+vdir[1]*vdir[1]);
  vdir[0] /= vmod;
  vdir[1] /= vmod;

  vper[0] = -vdir[1];
  vper[1] =  vdir[0];
  a1per = vper[1]/vper[0];

  a0per = ydr - a1per*xdr;
  xcross = (a0per-a0)/(a1-a1per);
  ycross = a1*xcross+a0;
  distance=sqrt((ydr-ycross)*(ydr-ycross)+(xdr-xcross)*(xdr-xcross));

  return distance;
}
