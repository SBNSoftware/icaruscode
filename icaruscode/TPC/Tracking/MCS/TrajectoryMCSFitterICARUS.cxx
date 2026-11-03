#include "TrajectoryMCSFitterICARUS.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardata/RecoBaseProxy/Track.h"

#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"

#include <numeric>

//toyMS stuff
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

using namespace std;
using namespace trkf;
using namespace recob::tracking;

//main function, return MCS momentum
recob::MCSFitResultGS TrajectoryMCSFitterICARUS::fitMcs(
  const recob::TrackTrajectory& traj, 
  int pid) const {
  //validity check
  bool check = false;

  //computation of sigma3p
  int nhits; double sigma3p; vector<double> delta3p;
  ComputeSigma3P(nhits, sigma3p, delta3p, wireSkip_);

  //define vector breakpoints with segmentation index
  vector<size_t> breakpoints; 
  //define vector seglens and cumseglens with segment lenghts [cm] and cumulative segment lengths [cm]
  vector<float> seglens; vector<float> cumseglens; 
  //define vector seghits and cumseghits with segment number of hits and cumulative segment number of hits
  vector<int> seghits; vector<int> cumseghits;
  vector<int> segsides;

  //break track into segments and populate vectors defined above
  breakTrajInSegments(traj, breakpoints, seglens, cumseglens, seghits, cumseghits, segsides);

  //check if number of segments is less than a certain value
  if (seglens.size() < minNSegs_) {
    if (verbose_ > 0) cout << "number of segments less than " << minNSegs_ << ", stopping fit" << endl;
    return recob::MCSFitResultGS(check); }
  //check if total number of hits is less than a certain value
  if (cumseghits.back() < minHits_) {
    if (verbose_ > 0) cout << "number of hits less than " << minHits_ << ", stopping fit" << endl;
    return recob::MCSFitResultGS(check); }

  //define trivial vector for check on directions
  auto trivial = Vector_t(0, 0, 0);

  //define vector dthetaLin with measured angles for linear fit
  vector<float> dthetaLin; dthetaLin.clear();
  //define vector psi0, psi1 with directions of adjacent segments
  Vector_t psi0; Vector_t psi1; 

  for (unsigned int p = 0; p < seglens.size(); p++) {
    if (verbose_ > 0) cout << "linear fit: iteration number " << p+1 << endl;
    //perform linear regression of current segment and memorize its direction in psi1
    linearRegression(traj, breakpoints[p], breakpoints[p+1], psi1);

    if (p == 0) {
      if (verbose_ > 0) cout << "first iteration, scattering angle impossible to compute" << endl; }
    if (p > 0) {
      //recall length [cm] of p-1, p segment
      double Ls0 = seglens[p-1]; 
      double Ls1 = seglens[p]; 
      if (verbose_ > 0) {
        cout << p << "-th segment length [cm] = " << Ls0 << endl;
        cout << p+1 << "-th segment length [cm] = " << Ls1 << endl; }
      
      //check if direction is backwards with respect to previous direction
      if (psi0.Dot(psi1) < 0.) psi1 *= -1.;
      if (verbose_ > 0) {
        cout << p << "-th measured direction = " << psi0 << endl;
        cout << p+1 << "-th measured direction = " << psi1 << endl; }

      //check if segment length is valid, return null if not valid
      if (Ls0 <= 0. || Ls1 <= 0.) {
        if (verbose_ > 0) {
          cout << "WARNING! invalid segment length found!" << endl; 
          cout << "measured angle [rad] between " << p << "-th and " << p+1 << "-th segment set as 0 " << endl; }
        dthetaLin.push_back(0.); }

      //check if segment direction is valid, return null if not valid
      else if (psi0 == trivial || psi1 == trivial) {
        if (verbose_ > 0) {
          cout << "WARNING! invalid segment direction found!" << endl;
          cout << "measured angle [rad] between " << p << "-th and " << p+1 << "-th segment set as 0 " << endl; }
        dthetaLin.push_back(0.); }

      else {
        //compute scalar product of adjacent segment directions, i.e. cosine of scattering angle
        double cosval = psi0.Dot(psi1);
        cosval = std::clamp(cosval, -1.0, 1.0);

        //compute scattering angle [rad] from scalar product
        double dt = acos(cosval);
        if (cosval == -1.0 || cosval == 1.0) dt = 0.001;

        //compute signed scattering angle
        //auto cross = psi0.Cross(psi1);
        //double sign = (cross.Dot(Vector_t(0, 0, 1)) >= 0) ? 1.0 : -1.0;
        //dt = dt * sign;

        //print scattering angle [rad] and add it to vector dthetaLin
        if (verbose_ > 0) cout << "measured angle [rad] between " << p << "-th and " << p+1 << "-th segment = " << dt << endl;
        if (cutMode_ == 0) dthetaLin.push_back(dt);
        else if (cutMode_ == 1) {
          if (p >= seglens.size() - cutAngles_) dthetaLin.push_back(0);
          else dthetaLin.push_back(dt); }
        else if (cutMode_ == 2) {
          if (p > cutAngles_) dthetaLin.push_back(0);
          else dthetaLin.push_back(dt); } } }
    psi0 = psi1; }

  //define vector dthetaPoly with measured angles for poligonal fit
  vector<float> dthetaPoly; dthetaPoly.clear(); 
  //define vector barycenters, bary with barycenter positions
  vector<Vector_t> barycenters; Vector_t bary; 

  for (unsigned int p = 0; p < seglens.size(); p++) {
    if (verbose_ > 0) cout << "poligonal fit: iteration number " << p+1 << endl;
    //find barycenter position of current segment and memorize it in bary
    findSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary);
    barycenters.push_back(bary);

    if (p == 0) {
      if (verbose_ > 0) cout << "first iteration, scattering angle impossible to compute" << endl; }
    if (p == 1) {
      if (verbose_ > 0) cout << "second iteration, scattering angle impossible to compute" << endl; }
    if (p > 1) {
      //recall length [cm] of p-2, p-1, p segment
      double Ls0 = seglens[p-2]; 
      double Ls1 = seglens[p-1]; 
      double Ls2 = seglens[p]; 
      if (verbose_ > 0) {
        cout << p-1 << "-th segment length [cm] = " << Ls0 << endl;
        cout << p << "-th segment length [cm] = " << Ls1 << endl;
        cout << p+1 << "-th segment length [cm] = " << Ls2 << endl; }

      //recall barycenter of p-2, p-1, p segment
      if (verbose_ > 0) {
        cout << p-1 << "-th barycenter = " << barycenters[p-2] << endl;
        cout << p << "-th barycenter = " << barycenters[p-1] << endl;
        cout << p+1 << "-th barycenter = " << barycenters[p] << endl; }

      //check if segment length is valid, return null if not valid
      if (Ls0 <= 0. || Ls1 <= 0. || Ls2 <= 0.) {
        if (verbose_ > 0) cout << "WARNING! invalid segment length found!" << endl; 
        dthetaPoly.push_back(0.); }

      //check if segment barycenter is valid, return null if not valid
      else if (barycenters[p-2] == trivial || barycenters[p-1] == trivial || barycenters[p] == trivial) {
        if (verbose_ > 0) {
          cout << "WARNING! invalid segment barycenter found!" << endl;
          cout << "measured angle [rad] between " << p-1 << "-th and " << p+1 << "-th segment set as 0 " << endl; }
        dthetaPoly.push_back(0.); }

      else {
        //compute direction between first (p-2) and second (p-1) barycenter, then normalize it
        Vector_t bcdir0 = (barycenters[p-1] - barycenters[p-2]).Unit();
        if (verbose_ > 0) cout << "direction between " << p-1 << "-th and " << p << "-th barycenter = " << bcdir0 << endl;

        //compute direction between second (p-1) and third (p) barycenter, then normalize it
        Vector_t bcdir1 = (barycenters[p] - barycenters[p-1]).Unit();
        if (verbose_ > 0) cout << "direction between " << p << "-th and " << p+1 << "-th barycenter = " << bcdir1 << endl;

        //compute scalar product of two adjacent directions, i.e. cosine of scattering angle
        double cosval = bcdir0.Dot(bcdir1);
        cosval = std::clamp(cosval, -1.0, 1.0);

        //compute scattering angle [rad] from scalar product
        double dt = acos(cosval);
        if (cosval == -1.0 || cosval == 1.0) dt = 0.001;

        //in case of 2D angle return signed scattering angle
        //auto cross = bcdir0.Cross(bcdir1);
        //double sign = (cross.Dot(ref) >= 0) ? 1.0 : -1.0;
        //dt = dt * sign;

        //print scattering angle [rad] and add it to vector dthetaPoly
        if (verbose_ > 0) cout << "measured angle [rad] between " << p-1 << "-th and " << p+1 << "-th segment = " << dt << endl;
        if (cutMode_ == 0) dthetaPoly.push_back(dt);
        else if (cutMode_ == 1) {
          if (p >= seglens.size() - cutAngles_) dthetaPoly.push_back(0);
          else dthetaPoly.push_back(dt); }
        else if (cutMode_ == 2) {
          if (p > cutAngles_) dthetaPoly.push_back(0);
          else dthetaPoly.push_back(dt); } } } }

  //print all vectors defined and populated at this point
  if (verbose_ > 0) {
    cout << "segment lengths [cm] = "; for (auto i : seglens) cout << i << ' '; cout << endl;
    cout << "cumulative segment lengths [cm] = "; for (auto i : cumseglens) cout << i << ' '; cout << endl;
    cout << "segment number of hits = "; for (auto i : seghits) cout << i << ' '; cout << endl;
    cout << "cumulative segment number of hits = "; for (auto i : cumseghits) cout << i << ' '; cout << endl;
    cout << "segment sides = "; for (auto i : segsides) cout << i << ' '; cout << endl;
    cout << "measured scattering angles, linear fit [rad] = "; for (auto i : dthetaLin) cout << i << ' '; cout << endl;
    cout << "measured scattering angles, poligonal fit [rad] = "; for (auto i : dthetaPoly) cout << i << ' '; cout << endl; }
  
  //check if there is minimum number of angles > 0
  unsigned int firstseg = 0; unsigned int lastseg = 0; 
  bool checkLin = false; bool checkPoly = false;
  ThetaCheck(dthetaLin, dthetaPoly, checkLin, checkPoly, firstseg, lastseg);
  if (!checkLin) {
    if (verbose_ > 0) {
      cout << "WARNING! number of scattering angles from linear fit is not greater or equal " << minNAngs_ << endl;
      cout << "fit stops here, return null" << endl; }
    return recob::MCSFitResultGS(check); }
  if (!checkPoly) {
    if (verbose_ > 0) {
      cout << "WARNING! number of scattering angles from poligonal fit is not greater or equal " << minNAngs_-1 << endl;
      cout << "fit stops here, return null" << endl; }
    return recob::MCSFitResultGS(check); }

  int nseg_used = lastseg - firstseg + 1;
  int nseg_tot = seglens.size();
  if (verbose_ > 0) cout << "nseg used = " << nseg_used << " nseg tot = " << nseg_tot << endl;

  //perform a c2 fit scan over all the computed angles
  float c2sigma = 0.05; 
  const ScanResult fitResult = C2Fit(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, c2sigma);

  //perform a likelihood fit scan over all the computed angles
  double best_p = 0.0; double error_p = 0.0;
  LikelihoodFit(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, best_p, error_p);

  //recall L1D, L2D, L3D functions
  double L1D = length1D(traj); 
  double L2D = length2D(traj);
  double L3D = length3D(traj); 

  //return pid, bestp and errp, sigma3p, length 1D 2D 3D, segment lens and number of hits, scattering linear-fit and poligonal angles
  check = true;
  return recob::MCSFitResultGS(
    pid, 
    fitResult.bestp, fitResult.errp, fitResult.minp, fitResult.maxp,
    fitResult.alpha, fitResult.dalpha, fitResult.beta, fitResult.dbeta,
    best_p, error_p,
    nhits, sigma3p, delta3p,
    L1D, L2D, L3D, 
    seglens, cumseglens,
    seghits, cumseghits, 
    segsides, 
    nseg_used, nseg_tot,
    dthetaLin, dthetaPoly,
    check); 
}

//break input trajectory into smaller pieces called segments 
void TrajectoryMCSFitterICARUS::breakTrajInSegments(
  const recob::TrackTrajectory& traj, 
  vector<size_t>& breakpoints, 
  vector<float>& seglens, 
  vector<float>& cumseglens, 
  vector<int>& seghits, 
  vector<int>& cumseghits, 
  vector<int>& segsides) const {
  //initialize vectors cumseglens and cumseghits (to use .back() later)
  cumseglens.push_back(0.); cumseghits.push_back(0);
  double seglen = 0.; int seghit = 0;

  //auto countsHit = [&](size_t i) { return IsInPlane(i) && IsInTPC(i, LastTPC(traj)) && CathodeCheck(traj, i); };
  auto countsHit = [&](size_t i) { return IsInPlane(i) && CathodeCheck(traj, i); };

  //determine index of first valid point of track, add it to vector breakpoints
  size_t index = traj.FirstValidPoint(); 
  breakpoints.push_back(index);

  //define position and side of first valid point
  auto pos0 = traj.LocationAtPoint(index);
  auto side0 = TPC(index);
  index = traj.NextValidPoint(index + 1);

  while (index != recob::TrackTrajectory::InvalidIndex) {
    //define position and side of current valid point, then update segment length and position of next valid point
    auto pos1 = traj.LocationAtPoint(index);
    auto side1 = TPC(index);
    seglen += (pos1 - pos0).R();
    pos0 = pos1;

    //break into a new segment if current segment length is track is changing side (cathode crossing)
    if (side1 != side0) {
      //add current valid point index to vector breakpoints
      breakpoints.push_back(index);

      //add current segment length [cm] and cumulative segment length [cm] to vector seglens, cumseglens
      seglens.push_back(0);
      cumseglens.push_back(cumseglens.back() + seglen);
      seglen = 0.;

      //add current segment number of hits and cumulative number of hits to vector seghits, cumseghits
      seghits.push_back(seghit);
      cumseghits.push_back(cumseghits.back() + seghit);
      seghit = countsHit(index) ? 1 : 0;

      //add current segment side to vector segsides
      segsides.push_back(side0);
      side0 = side1;
      index = traj.NextValidPoint(index + 1);
      continue;
    }

    //update number of hits
    if (countsHit(index)) seghit++;
    
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

      //add current segment side to vector segsides
      segsides.push_back(side0);
    }
    index = traj.NextValidPoint(index + 1);
  }

  //then add last segment
  if (seglen > 0.) {
    //add last valid point index to vector breakpoints
    breakpoints.push_back(traj.LastValidPoint() + 1);

    //add last segment length [cm] and total segment length [cm] to vector seglens, cumseglens
    if (seghit >= minHitsPerSegment_) seglens.push_back(seglen);
    else seglens.push_back(0.);
    cumseglens.push_back(cumseglens.back() + seglen);

    //add last segment number of hits and total number of hits to vector seghits, cumseghits
    seghits.push_back(seghit);
    cumseghits.push_back(cumseghits.back() + seghit);

    //add last segment side to vector segsides
    segsides.push_back(side0);
  }
  return;
}

//find barycenter position of input segment (useful for poligonal fit)
void TrajectoryMCSFitterICARUS::findSegmentBarycenter(
  const recob::TrackTrajectory& traj, 
  const size_t firstPoint, 
  const size_t lastPoint, 
  Vector_t& bary) const {
  //initalize number of points in segment
  int npoints = 0;

  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;
  
  size_t index = firstPoint;
  while (index < lastPoint) {
    if (IsInPlane(index) && CathodeCheck(traj, index)) {
      if ((fullTrackMode_ == 0) && IsInTPC(index, FirstTPC(traj))) {
        middlePointCalc.add(hit2D(index));
        npoints++; }
      else if ((fullTrackMode_ == 1) && IsInTPC(index, LastTPC(traj))) {
        middlePointCalc.add(hit2D(index));
        npoints++; }
      else if (fullTrackMode_ == 2) {
        middlePointCalc.add(hit2D(index));
        npoints++; } }
    index = traj.NextValidPoint(index + 1); }

  //check if number of points is greater than zero, otherwise return trivial vector
  if (npoints == 0) bary = Vector_t();
  //determine position of segment barycenter from position of valid points in the segment
  else bary = middlePointCalc.middlePoint();
}

//find average direction of trajectory between firstPoint and lastPoint
void TrajectoryMCSFitterICARUS::linearRegression(
  const recob::TrackTrajectory& traj, 
  const size_t firstPoint, 
  const size_t lastPoint, 
  Vector_t& pcdir) const {
  //initalize number of points in segment
  int npoints = 0;

  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;

  size_t index = firstPoint;
  while (index < lastPoint) {
    if (IsInPlane(index) && CathodeCheck(traj, index)) {
      if ((fullTrackMode_ == 0) && IsInTPC(index, FirstTPC(traj))) {
        middlePointCalc.add(hit2D(index));
        npoints++; }
      else if ((fullTrackMode_ == 1) && IsInTPC(index, LastTPC(traj))) {
        middlePointCalc.add(hit2D(index));
        npoints++; }
      else if (fullTrackMode_ == 2) {
        middlePointCalc.add(hit2D(index));
        npoints++; } } 
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
      if (IsInPlane(index) && CathodeCheck(traj, index)) {
        if ((fullTrackMode_ == 0) && IsInTPC(index, FirstTPC(traj))) {
          auto p = hit2D(index);
          const double xxw0 = p.X() - avgpos.X();
          const double yyw0 = p.Y() - avgpos.Y();
          const double zzw0 = p.Z() - avgpos.Z();
          m(0, 0) += xxw0 * xxw0 * norm; m(0, 1) += xxw0 * yyw0 * norm; m(0, 2) += xxw0 * zzw0 * norm;
          m(1, 0) += yyw0 * xxw0 * norm; m(1, 1) += yyw0 * yyw0 * norm; m(1, 2) += yyw0 * zzw0 * norm;
          m(2, 0) += zzw0 * xxw0 * norm; m(2, 1) += zzw0 * yyw0 * norm; m(2, 2) += zzw0 * zzw0 * norm; } 
        else if ((fullTrackMode_ == 1) && IsInTPC(index, LastTPC(traj))) {
          auto p = hit2D(index);
          const double xxw0 = p.X() - avgpos.X();
          const double yyw0 = p.Y() - avgpos.Y();
          const double zzw0 = p.Z() - avgpos.Z();
          m(0, 0) += xxw0 * xxw0 * norm; m(0, 1) += xxw0 * yyw0 * norm; m(0, 2) += xxw0 * zzw0 * norm;
          m(1, 0) += yyw0 * xxw0 * norm; m(1, 1) += yyw0 * yyw0 * norm; m(1, 2) += yyw0 * zzw0 * norm;
          m(2, 0) += zzw0 * xxw0 * norm; m(2, 1) += zzw0 * yyw0 * norm; m(2, 2) += zzw0 * zzw0 * norm; } 
        else if (fullTrackMode_ == 2) {
          auto p = hit2D(index);
          const double xxw0 = p.X() - avgpos.X();
          const double yyw0 = p.Y() - avgpos.Y();
          const double zzw0 = p.Z() - avgpos.Z();
          m(0, 0) += xxw0 * xxw0 * norm; m(0, 1) += xxw0 * yyw0 * norm; m(0, 2) += xxw0 * zzw0 * norm;
          m(1, 0) += yyw0 * xxw0 * norm; m(1, 1) += yyw0 * yyw0 * norm; m(1, 2) += yyw0 * zzw0 * norm;
          m(2, 0) += zzw0 * xxw0 * norm; m(2, 1) += zzw0 * yyw0 * norm; m(2, 2) += zzw0 * zzw0 * norm; } }

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
  }
}

//find most probable value of energy loss according to Landau distribution
double TrajectoryMCSFitterICARUS::energyLossLandau(
  const double mass, 
  const double energy) const {
  //define muon mass [GeV c-2], energy [GeV], beta, gamma
  double mass2 = mass * mass; 
  double energy2 = energy * energy;
  double beta2 = (energy2 - mass2) / energy2;  
  double gamma2 = 1 / (1. - beta2); 

  //define ionization potential I [MeV]
  double I = 188.E-6; 
  double I2 = I * I;

  //define density of liquid Argon rho, atomic number Z, mass number A, K [MeV cm-1]
  double rho = 1.396; double Z = 18; double A = 39.95; double K = 0.307075;
  double Const = K * rho * Z / A;

  //define electron mass [MeV c-2] 
  double me = 0.511; 

  //define Landau parameter j
  double j = 0.200;

  //return most probable value of energy loss [GeV] according to Landau distribution
  double xi = 0.5 * Const / beta2;
  double argument = 2 * me * beta2 * gamma2 * xi / I2;
  double MPV = xi * (log(argument) + j - beta2);
  return 0.001 * MPV;
}

//find mean energy loss according to Bethe-Bloch distribution
double TrajectoryMCSFitterICARUS::energyLossBetheBloch(
  const double mass, 
  const double energy) const {
  //define muon mass [GeV c-2], energy [GeV], beta, gamma
  double mass2 = mass * mass; 
  double energy2 = energy * energy;
  double beta2 = (energy2 - mass2) / energy2;  
  double gamma2 = 1 / (1. - beta2); 
  double gamma = sqrt(gamma2);

  //define ionization potential I [MeV]
  double I = 188.E-6; 
  double I2 = I * I;

  //define density of liquid Argon rho, atomic number Z, mass number A, K [MeV cm-1]
  double rho = 1.396; double Z = 18; double A = 39.95; double K = 0.307075;
  double Const = K * rho * Z / A;

  //define electron mass [MeV c-2] and mass ratio electron/muon
  double me = 0.511; 
  double mass_ratio = me / (1000 * mass);

  //define energy transfer Wmax [MeV]
  double Wmax = 2 * me * beta2 * gamma2;
  Wmax = Wmax / (1 + 2 * gamma * mass_ratio + mass_ratio * mass_ratio);

  //return mean value of energy loss [GeV] according to Bethe-Bloch distribution
  double argument = 2 * me * beta2 * gamma2 * Wmax / I2;
  double mean = (Const / beta2) * (0.5 * log(argument) - beta2);
  if (mean < 0.) return 0.;
  else return 0.001 * mean;
}

//find energy of track at a certain distance travelled
double TrajectoryMCSFitterICARUS::GetE(
  const double initial_E, 
  const double length_travelled, 
  const double mass) const {
  //divide length travelled in N pieces of 1 cm length
  unsigned int N = static_cast<int>(length_travelled);
  //define current energy [GeV]
  double current_E = initial_E;

  //compute energy according to Landau distribution (energy loss as MPV)
  if (eLossMode_ == 0) {
    for (unsigned int i = 0; i < N; ++i) {
      current_E -= energyLossLandau(mass, current_E);
      if (current_E <= mass) return 0.; }
    return current_E; }

  //compute energy according to MIP hypothesis (constant energy loss)
  else if (eLossMode_ == 1) {
    //define energy loss [GeV cm-1] for MIP muon in liquid Argon
    double kcal = 0.002105;
    current_E = initial_E - kcal * length_travelled;
    return current_E; }
  
  //compute energy according to Bethe-Bloch distribution (energy loss as mean)
  else if (eLossMode_ == 2) {
    for (unsigned int i = 0; i < N; ++i) {
      current_E -= energyLossBetheBloch(mass, current_E);
      if (current_E <= mass) return 0.; }
    return current_E; }

  else return 0.;
}

//compute sigma3p [mm] used for measurement error along drift direction
void TrajectoryMCSFitterICARUS::ComputeSigma3P(
  int& nhits, 
  double& sigma3p, 
  vector<double>& delta3p, 
  int wire_skip) const {
  double eps = 1e-6;

  vector<recob::Hit> hits;
  if (planeMode_ == 0) hits = hits2dI1;
  else if (planeMode_ == 1) hits = hits2dI2;
  else if (planeMode_ == 2) hits = hits2dC;

  nhits = hits.size(); sigma3p = 0.; delta3p.clear();
  if (nhits < 30) return; 

  double square_mean = 0.; double mean_square = 0.; int n = 0;
  double range = 0.3;

  wire_skip = wire_skip + 1;
  for (unsigned int i = 0; i < hits.size() - 2*wire_skip; i++) {
    recob::Hit h0 = hits.at(i);
    recob::Hit h1 = hits.at(i + wire_skip);
    recob::Hit h2 = hits.at(i + 2*wire_skip);

    double x0 = h0.PeakTime() * 0.0622; double y0 = h0.WireID().Wire * 0.3;
    double x1 = h1.PeakTime() * 0.0622; double y1 = h1.WireID().Wire * 0.3;
    double x2 = h2.PeakTime() * 0.0622; double y2 = h2.WireID().Wire * 0.3;

    if (abs(x2-x0) < eps || abs(x2-x1) < eps || abs(x1-x0) < eps || abs(y2-y0) < eps) continue;

    double xm = x0 + (x2-x0) * (y1-y0) / (y2-y0); 
    double d3p = x1 - xm;

    if (abs(d3p) < range) { 
      mean_square += d3p*d3p; square_mean += d3p; ++n; 
      delta3p.push_back(d3p); } }

  mean_square /= n; square_mean /= n; square_mean *= square_mean;
  sigma3p = sqrt(max(0.0, mean_square - square_mean));
  return; }

//build c2 function used in MCS fit
const double TrajectoryMCSFitterICARUS::C2Function(
  const recob::TrackTrajectory& traj, 
  std::vector<float>& seglens, 
  std::vector<float>& cumseglens, 
  std::vector<int>& seghits, 
  std::vector<float> dthetaLin, 
  std::vector<float> dthetaPoly, 
  double sigma3p, 
  int pid, 
  double input_p) const {
  //define theta_measured, theta_expected, theta_err that represent expected, measured and error angle
  double theta_measured, theta_expected, theta_err;

  //define ttallLin, ttallPoly 
  vector<float> ttall; vector<float> ttallLin; vector<float> ttallPoly;
  ttall.clear(); ttallLin.clear(); ttallPoly.clear();

  //check if there is minimum number of angles > 0
  unsigned int firstseg = 0; unsigned int lastseg = 0; 
  bool checkLin = false; bool checkPoly = false;
  ThetaCheck(dthetaLin, dthetaPoly, checkLin, checkPoly, firstseg, lastseg);
  bool check = checkLin && checkPoly;
  if (verbose_ > 1) cout << "first = " << firstseg << " last = " << lastseg << endl;
  if (!check) return 0;

  //compute number of valid segments 
  unsigned int nseg = lastseg - firstseg + 1;
  if (fullTrackMode_ == 0) {
    if (verbose_ > 1) cout << "nseg tot = " << seglens.size() << " nseg used = " << nseg << " minNSegs_ = " << minNSegs_ << endl;
    if (nseg < minNSegs_) return 0; }
  if (fullTrackMode_ == 1) {
    unsigned int nsegValid = 0;
    for (unsigned int j = firstseg; j <= lastseg && j < seglens.size(); ++j) if (seglens[j] > 0.) nsegValid++;
    if (verbose_ > 1) cout << "nseg tot = " << seglens.size() << " nseg = " << nseg << endl;
    if (verbose_ > 1) cout << "nseg tot = " << seglens.size() << " nseg used = " << nsegValid << " minNSegs_ = " << minNSegs_ << endl;
    if (nsegValid < minNSegs_) return 0; }

  vector<bool> keepLin(nseg-1, false), keepPoly(nseg-2, false);

  //define covariance matrices for both linear fit and poligonal
  TMatrixDSym mat(2*nseg-3);
  TMatrixDSym matFull(2*nseg-3);

  //define covariance matrices for linear fit
  TMatrixDSym matLin(nseg-1);
  TMatrixDSym matMCSLin(nseg-1);
  TMatrixDSym matErrLin(nseg-1);
  TMatrixDSym matOffDiagLin(nseg-1);

  //define covariance matrices for poligonal fit
  TMatrixDSym matPoly(nseg-2);
  TMatrixDSym matMCSPoly(nseg-2);
  TMatrixDSym matErrPoly(nseg-2);
  TMatrixDSym matOffDiagPoly(nseg-2);

  //define covariance matrix with mix terms
  TMatrixD matMCSMixed(nseg-1, nseg-2);

  //define propagation constant for measurement error associated to delta3p
  double kLin = sqrt(24.); double kPoly = sqrt(6.);

  //convert input momentum from [MeV/c] to [GeV/c]
  const double P = 0.001 * input_p;
  //compute total energy [GeV] from input momentum [GeV/c]
  const double PS = pow(P, 2);
  const double MS = pow(mass(pid), 2);
  const double E = sqrt(PS + MS);

  //notify computation of linear-fit angles for C2 function has started
  if (verbose_ > 1) cout << "starting computation linear-fit angles for C2 function ..." << endl;

  for (unsigned int jp = firstseg; jp < lastseg; jp++) {
    if (fullTrackMode_ == 2) { 
      if (dthetaLin[jp] <= 0.) continue; }
    //recall length [cm] of jp, jp+1 segment, return null if less or equal than zero
    double Ls0 = seglens[jp]; double Ls1 = seglens[jp+1];
    if (Ls0 <= 0. || Ls1 <= 0.) {
      if (verbose_ > 1) {
        cout << jp << "-th segment length [cm] = " << Ls0 << endl;
        cout << jp+1 << "-th segment length [cm] = " << Ls1 << endl;
        cout << "WARNING! invalid segment length found!" << endl;
        cout << "C2 function set as 0" << endl; }
      return 0; }

    //recall number of hits of jp, jp+1 segment, return null if less or equal than zero
    int nhit0 = seghits[jp]; int nhit1 = seghits[jp+1];
    if (nhit0 <= 0 || nhit1 <= 0) {
      if (verbose_ > 1) {
        cout << jp << "-th segment number of hit = " << nhit0 << endl;
        cout << jp+1 << "-th segment number of hit = " << nhit1 << endl;
        cout << "WARNING! invalid number of hits found!" << endl;
        cout << "C2 function set as 0" << endl; }
      return 0; }

    //compute energy [GeV] at jp, jp+1 segment, return null if less or equal than energy at rest
    const double E0 = GetE(E, cumseglens[jp], mass(pid));
    const double E1 = GetE(E, cumseglens[jp+1], mass(pid));
    const double E0S = E0 * E0; const double E1S = E1 * E1;
    if (E0S <= MS || E1S <= MS) {
      if (verbose_ > 1) {
        cout << jp << "-th segment energy [GeV] = " << E0 << endl;
        cout << jp+1 << "-th segment energy [GeV] = " << E1 << endl;
        cout << "WARNING! energy less than energy at rest!" << endl;
        cout << "C2 function set as 0" << endl; }
      return 0; }
    
    //compute momentum [GeV/c] and velocity of particle at jp, jp+1 segment
    const double P0 = sqrt(E0S - MS); 
    const double P1 = sqrt(E1S - MS); 
    const double beta0 = P0 / E0; 
    const double beta1 = P1 / E1; 

    //compute Highland term [mrad] just as claimed by Highland formula
    double var0 = ThetaExpected(traj, P0, beta0, Ls0);
    double var1 = ThetaExpected(traj, P1, beta1, Ls1);

    //print a few things about jp iteration
    if (verbose_ > 1) {
      cout << jp << "-th segment length [cm] = " << Ls0 << endl;
      cout << jp << "-th segment number of hit = " << nhit0 << endl;
      cout << jp << "-th segment energy [GeV] = " << E0 << endl;
      cout << jp << "-th segment momentum [GeV/c] = " << P0 << endl;
      cout << jp << "-th segment beta = v/c = " << beta0 << endl;
      cout << jp << "-th Highland term [mrad] = " << var0 << endl; }

    //print a few things about jp+1 iteration
    if (verbose_ > 1) {
      cout << jp+1 << "-th segment length [cm] = " << Ls1 << endl;
      cout << jp+1 << "-th segment number of hit = " << nhit1 << endl;
      cout << jp+1 << "-th segment energy [GeV] = " << E1 << endl;
      cout << jp+1 << "-th segment momentum [GeV/c] = " << P1 << endl;
      cout << jp+1 << "-th segment beta = v/c = " << beta1 << endl;
      cout << jp+1 << "-th Highland term [mrad] = " << var1 << endl; }
    
    //compute expected theta [mrad] from Highland terms
    var0 *= var0; var0 /= 3;
    var1 *= var1; var1 /= 3;
    theta_expected = sqrt(var0 + var1); 
    if (verbose_ > 1) cout << "expected angle [mrad] between " << jp << "-th and " << jp+1 << "-th segment = " << theta_expected << endl;

    //update ttall vector with measured scattering angle [mrad]
    theta_measured = 1000. * dthetaLin[jp];
    ttallLin.push_back(theta_measured);
    if (fullTrackMode_ == 2) { 
      keepLin[jp-firstseg] = true; }
    if (verbose_ > 1) cout << "measured angle [mrad] between " << jp << "-th and " << jp+1 << "-th segment = " << theta_measured << endl;

    //compute theta error [mrad] associated to delta3p
    theta_err = 1000 * (kLin * sigma3p)/(Ls0 * sqrt(float(nhit0))) * ThetaErrFactor(traj);
    if (verbose_ > 1) cout << "   error angle [mrad] between " << jp << "-th and " << jp+1 << "-th segment = " << theta_err << endl;

    //update matrix matMCSLin, matErrLin with on-diag MCS and error terms
    matMCSLin(jp-firstseg, jp-firstseg) = theta_expected * theta_expected;
    matErrLin(jp-firstseg, jp-firstseg) = theta_err * theta_err;

    //update matrix matOffDiagLin with off-diag MCS and error terms
    FillOffDiagMCSLin(matOffDiagLin, matMCSLin, jp-firstseg);
    FillOffDiagErrLin(matOffDiagLin, matErrLin, jp-firstseg); }

  //notify computation of polygonal angles for C2 function has started
  if (verbose_ > 1) cout << "starting computation polygonal angles for C2 function ..." << endl;

  for (unsigned int jp = firstseg; jp < lastseg - 1; jp++) {
    if (fullTrackMode_ == 2) { 
      if (dthetaPoly[jp] <= 0.) continue; }
    //recall length [cm] of jp, jp+1, jp+2 segment, return null if less or equal than zero
    double Ls0 = seglens[jp]; double Ls1 = seglens[jp+1]; double Ls2 = seglens[jp+2];
    if (Ls0 <= 0. || Ls1 <= 0. || Ls2 <= 0.) {
      if (verbose_ > 1) {
        cout << jp << "-th segment length [cm] = " << Ls0 << endl;
        cout << jp+1 << "-th segment length [cm] = " << Ls1 << endl;
        cout << jp+2 << "-th segment length [cm] = " << Ls2 << endl;
        cout << "WARNING! invalid segment length found!" << endl;
        cout << "C2 function set as 0" << endl; }
      return 0; }

    //recall number of hits of jp, jp+1, jp+2 segment, return null if less or equal than zero
    int nhit0 = seghits[jp]; int nhit1 = seghits[jp+1]; int nhit2 = seghits[jp+2];
    if (nhit0 <= 0 || nhit1 <= 0 || nhit2 <= 0) {
      if (verbose_ > 1) {
        cout << jp << "-th segment number of hit = " << nhit0 << endl;
        cout << jp+1 << "-th segment number of hit = " << nhit1 << endl;
        cout << jp+2 << "-th segment number of hit = " << nhit2 << endl;
        cout << "WARNING! invalid number of hits found!" << endl;
        cout << "C2 function set as 0" << endl; }
      return 0; }

    //compute energy [GeV] at jp, jp+1, jp+2 segment, return null if less or equal than energy at rest
    const double E0 = GetE(E, cumseglens[jp], mass(pid));
    const double E1 = GetE(E, cumseglens[jp+1], mass(pid));
    const double E2 = GetE(E, cumseglens[jp+2], mass(pid));
    const double E0S = E0 * E0; const double E1S = E1 * E1; const double E2S = E2 * E2;
    if (E0S <= MS || E1S <= MS || E2S <= MS) {
      if (verbose_ > 1) {
        cout << jp << "-th segment energy [GeV] = " << E0 << endl;
        cout << jp+1 << "-th segment energy [GeV] = " << E1 << endl;
        cout << jp+2 << "-th segment energy [GeV] = " << E2 << endl;
        cout << "WARNING! energy less than energy at rest!" << endl;
        cout << "C2 function set as 0" << endl; }
      return 0; }
    
    //compute momentum [GeV/c] and velocity of particle at jp, jp+1, jp+2 segment
    const double P0 = sqrt(E0S - MS); 
    const double P1 = sqrt(E1S - MS); 
    const double P2 = sqrt(E2S - MS); 
    const double beta0 = P0 / E0; 
    const double beta1 = P1 / E1; 
    const double beta2 = P2 / E2; 

    //compute Highland term [mrad] just as claimed by Highland formula
    double var0 = ThetaExpected(traj, P0, beta0, Ls0);
    double var1 = ThetaExpected(traj, P1, beta1, Ls1);
    double var2 = ThetaExpected(traj, P2, beta2, Ls2);

    //print a few things about jp iteration
    if (verbose_ > 1) {
      cout << jp << "-th segment length [cm] = " << Ls0 << endl;
      cout << jp << "-th segment number of hit = " << nhit0 << endl;
      cout << jp << "-th segment energy [GeV] = " << E0 << endl;
      cout << jp << "-th segment momentum [GeV/c] = " << P0 << endl;
      cout << jp << "-th segment beta = v/c = " << beta0 << endl;
      cout << jp << "-th Highland term [mrad] = " << var0 << endl; }

    //print a few things about jp+1 iteration
    if (verbose_ > 1) {
      cout << jp+1 << "-th segment length [cm] = " << Ls1 << endl;
      cout << jp+1 << "-th segment number of hit = " << nhit1 << endl;
      cout << jp+1 << "-th segment energy [GeV] = " << E1 << endl;
      cout << jp+1 << "-th segment momentum [GeV/c] = " << P1 << endl;
      cout << jp+1 << "-th segment beta = v/c = " << beta1 << endl;
      cout << jp+1 << "-th Highland term [mrad] = " << var1 << endl; }

    //print a few things about jp+2 iteration
    if (verbose_ > 1) {
      cout << jp+2 << "-th segment length [cm] = " << Ls2 << endl;
      cout << jp+2 << "-th segment number of hit = " << nhit2 << endl;
      cout << jp+2 << "-th segment energy [GeV] = " << E2 << endl;
      cout << jp+2 << "-th segment momentum [GeV/c] = " << P2 << endl;
      cout << jp+2 << "-th segment beta = v/c = " << beta2 << endl;
      cout << jp+2 << "-th Highland term [mrad] = " << var2 << endl; }
    
    //compute expected theta [mrad] from Highland terms
    var0 *= var0; var0 /= 12;
    var1 *= var1; var1 /= 3;
    var2 *= var2; var2 /= 12;
    theta_expected = sqrt(var0 + var1 + var2); 
    if (verbose_ > 1) cout << "expected angle [mrad] between " << jp << "-th and " << jp+2 << "-th segment = " << theta_expected << endl;

    //update ttall vector with measured scattering angle [mrad]
    theta_measured = 1000. * dthetaPoly[jp];
    ttallPoly.push_back(theta_measured);
    if (fullTrackMode_ == 2) { 
      keepPoly[jp-firstseg] = true; }
    if (verbose_ > 1) cout << "measured angle [mrad] between " << jp << "-th and " << jp+2 << "-th segment = " << theta_measured << endl;

    //compute theta error [mrad] associated to delta3p
    theta_err = 1000 * (kPoly * sigma3p)/(Ls0 * sqrt(float(nhit0))) * ThetaErrFactor(traj);
    if (verbose_ > 1) cout << "   error angle [mrad] between " << jp << "-th and " << jp+1 << "-th segment = " << theta_err << endl;

    //update matrix matMCSPoly, matErrPoly with on-diag MCS and error terms
    matMCSPoly(jp-firstseg, jp-firstseg) = theta_expected * theta_expected;
    matErrPoly(jp-firstseg, jp-firstseg) = theta_err * theta_err;

    //update matrix matOffDiagPoly with off-diag MCS and error terms
    FillOffDiagMCSPoly(matOffDiagPoly, matMCSPoly, jp-firstseg); 
    FillOffDiagErrPoly(matOffDiagPoly, matErrPoly, jp-firstseg); }

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

  auto compress = [](const TMatrixDSym& m, const vector<bool>& keep) {
    vector<int> idx;
    for (int i = 0; i < (int)keep.size(); ++i) 
      if (keep[i]) idx.push_back(i);
    TMatrixDSym out(idx.size());
    for (unsigned a = 0; a < idx.size(); ++a)
      for (unsigned b = 0; b < idx.size(); ++b) 
        out(a, b) = m(idx[a], idx[b]);
    return out; };
  vector<bool> keepFull = keepLin;
  keepFull.insert(keepFull.end(), keepPoly.begin(), keepPoly.end());

  //in case of both linear and poligonal fit, take matFull
  if (fitMode_ == 0) { 
    if (fullTrackMode_ == 0 || fullTrackMode_ == 1) {
      mat = matFull; }
    else if (fullTrackMode_ == 2) {
      TMatrixDSym matComp = compress(matFull, keepFull); 
      mat.ResizeTo(matComp.GetNrows(), matComp.GetNrows()); 
      mat = matComp; }
    ttall.insert(ttall.end(), ttallLin.begin(), ttallLin.end());
    ttall.insert(ttall.end(), ttallPoly.begin(), ttallPoly.end()); }
  //in case of only linear fit, take matLin
  if (fitMode_ == 1) { 
    if (fullTrackMode_ == 0 || fullTrackMode_ == 1) {
      mat.ResizeTo(matMCSLin.GetNrows(), matMCSLin.GetNrows());
      mat = matLin; }
    else if (fullTrackMode_ == 2) {
      TMatrixDSym matComp = compress(matLin, keepLin);
      mat.ResizeTo(matComp.GetNrows(), matComp.GetNrows()); 
      mat = matComp; }
    ttall = ttallLin; }
  //in case of only poligonal fit, take matPoly
  if (fitMode_ == 2) { 
    if (fullTrackMode_ == 0 || fullTrackMode_ == 1) {
      mat.ResizeTo(matMCSPoly.GetNrows(), matMCSPoly.GetNrows());
      mat = matPoly; }
    else if (fullTrackMode_ == 2) {
      TMatrixDSym matComp = compress(matPoly, keepPoly); 
      mat.ResizeTo(matComp.GetNrows(), matComp.GetNrows()); 
      mat = matComp; }
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
  //compute c2 function without trunc as average of mterms=
  if (ttall.empty()) return 0;
  vtcovmed /= (ttall.size());

  //define vector ttrunc that will be filled with accepted terms after cut
  vector<double> ttrunc;
  //define vector tails that will be filles with indexes of rejected terms
  vector<int> tails;
  //define ttrunctot as sum of accepted terms
  double ttrunctot = 0;
  //define vector atrunc that will be filled with ttall values of accepted terms
  vector<double> atrunc;

  //iterate over elements of vector terms
  double gaus = 3.;
  for (unsigned int jt = 0; jt < ttall.size(); jt++) {
    //check if element jt-th of vector terms is under threshold and observed angle is positive
    if (terms[jt] < vtcovmed * pow(gaus, 2) && ttall[jt] > 0.) {
      ttrunctot += terms[jt];
      atrunc.push_back(ttall[jt]);
      ttrunc.push_back(terms[jt]); } 
    else { 
      tails.push_back(jt); } }
  
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
  //compute c2 function with trunc as average of mtterms
  vmediomod /= (ttrunc.size());
  return vmediomod; 
}

//build likelihood function used in MCS fit
const double TrajectoryMCSFitterICARUS::LikelihoodFunction(
  const recob::TrackTrajectory& traj, 
  std::vector<float>& seglens, 
  std::vector<float>& cumseglens, 
  std::vector<int>& seghits, 
  std::vector<float> dthetaLin, 
  std::vector<float> dthetaPoly, 
  double sigma3p, 
  int pid, 
  double input_p, 
  bool verbose) const {
  //initalize log likelihood value
  double loglikelihood = 0.0;

  //check if there is minimum number of angles > 0
  unsigned int firstseg = 0; unsigned int lastseg = 0; 
  bool checkLin = false; bool checkPoly = false;
  ThetaCheck(dthetaLin, dthetaPoly, checkLin, checkPoly, firstseg, lastseg);
  bool check = checkLin && checkPoly;
  if (verbose_ > 1) cout << "first = " << firstseg << " last = " << lastseg << endl;
  if (!check) return numeric_limits<double>::max();

  //compute number of valid segments 
  unsigned int nseg = lastseg - firstseg + 1;
  if (fullTrackMode_ == 0) {
    if (nseg < minNSegs_) return numeric_limits<double>::max(); }
  if (fullTrackMode_ == 1) {
    unsigned int nsegValid = 0;
    for (unsigned int j = firstseg; j <= lastseg && j < seglens.size(); ++j) if (seglens[j] > 0.) nsegValid++;
    if (verbose_ > 1) cout << "nseg = " << nseg << " nsegValid = " << nsegValid << " minNSegs_ = " << minNSegs_ << endl;
    if (nsegValid < minNSegs_) return numeric_limits<double>::max(); }

  //define propagation constant for measurement error associated to delta3p
  double kLin = sqrt(24.); 

  //convert input momentum from [MeV/c] to [GeV/c]
  const double P = 0.001 * input_p;
  //compute total energy [GeV] from input momentum [GeV/c]
  const double PS = pow(P, 2);
  const double MS = pow(mass(pid), 2);
  const double E = sqrt(PS + MS);

  //notify computation of linear-fit angles for likelihood function has started
  if (verbose_ > 1) cout << "starting computation linear-fit angles for likelihood function ..." << endl; 

  for (unsigned int jp = firstseg; jp < lastseg; jp++) {
    if (fullTrackMode_ == 2) { 
      if (dthetaLin[jp] <= 0.) continue; }
    //recall length [cm] of jp, jp+1 segment, return null if less or equal than zero
    double Ls0 = seglens[jp]; double Ls1 = seglens[jp+1];
    if (Ls0 <= 0. || Ls1 <= 0.) {
      if (verbose_ > 1) {
        cout << jp << "-th segment length [cm] = " << Ls0 << endl;
        cout << jp+1 << "-th segment length [cm] = " << Ls1 << endl;
        cout << "WARNING! invalid segment length found!" << endl;
        cout << "C2 function set as 0" << endl; }
      return numeric_limits<double>::max(); }

    //recall number of hits of jp, jp+1 segment, return null if less or equal than zero
    int nhit0 = seghits[jp]; int nhit1 = seghits[jp+1];
    if (nhit0 <= 0 || nhit1 <= 0) {
      if (verbose_ > 1) {
        cout << jp << "-th segment number of hit = " << nhit0 << endl;
        cout << jp+1 << "-th segment number of hit = " << nhit1 << endl;
        cout << "WARNING! invalid number of hits found!" << endl;
        cout << "C2 function set as 0" << endl; }
      return numeric_limits<double>::max(); }

    //compute energy [GeV] at jp, jp+1 segment, return null if less or equal than energy at rest
    const double E0 = GetE(E, cumseglens[jp], mass(pid));
    const double E1 = GetE(E, cumseglens[jp+1], mass(pid));
    const double E0S = E0 * E0; const double E1S = E1 * E1;
    if (E0S <= MS || E1S <= MS) {
      if (verbose_ > 1) {
        cout << jp << "-th segment energy [GeV] = " << E0 << endl;
        cout << jp+1 << "-th segment energy [GeV] = " << E1 << endl;
        cout << "WARNING! energy less than energy at rest!" << endl;
        cout << "C2 function set as 0" << endl; }
      return numeric_limits<double>::max(); }
    
    //compute momentum [GeV/c] and velocity of particle at jp, jp+1 segment
    const double P0 = sqrt(E0S - MS); 
    const double P1 = sqrt(E1S - MS); 
    const double beta0 = P0 / E0; 
    const double beta1 = P1 / E1; 

    //compute Highland term [mrad] just as claimed by Highland formula
    double var0 = ThetaExpected(traj, P0, beta0, Ls0);
    double var1 = ThetaExpected(traj, P1, beta1, Ls1);

    //print a few things about jp iteration
    if (verbose_ > 1) {
      cout << jp << "-th segment length [cm] = " << Ls0 << endl;
      cout << jp << "-th segment number of hit = " << nhit0 << endl;
      cout << jp << "-th segment energy [GeV] = " << E0 << endl;
      cout << jp << "-th segment momentum [GeV/c] = " << P0 << endl;
      cout << jp << "-th segment beta = v/c = " << beta0 << endl;
      cout << jp << "-th Highland term [mrad] = " << var0 << endl; }

    //print a few things about jp+1 iteration
    if (verbose_ > 1) {
      cout << jp+1 << "-th segment length [cm] = " << Ls1 << endl;
      cout << jp+1 << "-th segment number of hit = " << nhit1 << endl;
      cout << jp+1 << "-th segment energy [GeV] = " << E1 << endl;
      cout << jp+1 << "-th segment momentum [GeV/c] = " << P1 << endl;
      cout << jp+1 << "-th segment beta = v/c = " << beta1 << endl;
      cout << jp+1 << "-th Highland term [mrad] = " << var1 << endl; }
    
    //compute expected theta [mrad] from Highland terms
    var0 *= var0; var0 /= 3;
    var1 *= var1; var1 /= 3;
    double theta_MCS = sqrt(var0 + var1); 
    if (verbose_ > 1) cout << "expected angle [mrad] between " << jp << "-th and " << jp+1 << "-th segment = " << theta_MCS << endl;

    //update ttall vector with measured scattering angle [mrad]
    double theta_measured = 1000. * dthetaLin[jp];
    if (verbose_ > 1) cout << "measured angle [mrad] between " << jp << "-th and " << jp+1 << "-th segment = " << theta_measured << endl;

    //compute theta error [mrad] associated to delta3p
    double theta_error = 1000 * (kLin * sigma3p)/(Ls0 * sqrt(float(nhit0))) * ThetaErrFactor(traj);
    if (verbose_ > 1) cout << "   error angle [mrad] between " << jp << "-th and " << jp+1 << "-th segment = " << theta_error << endl;

    //compute negative log likelihood
    double theta_expected = sqrt(theta_MCS * theta_MCS + theta_error * theta_error);
    //return infinity log likelihood if total MCS angle is null
    if (theta_expected == 0.0) {
      if (verbose_ > 1) cout << "WARNING! expected angle cannot be zero" << endl;
      return numeric_limits<double>::max(); } 

    //compute ratio between measured and expected scattering angles
    double theta_ratio = theta_measured / theta_expected;
    //compute log likelihood as sum of log, quadratic, fixed terms for every iteration
    loglikelihood += std::log(theta_expected) + 0.5 * theta_ratio * theta_ratio; }

  //return final value of log likelihood
  return loglikelihood; }

//perform fit of c2 function to determine MCS momentum
const TrajectoryMCSFitterICARUS::ScanResult TrajectoryMCSFitterICARUS::C2Fit(
  const recob::TrackTrajectory& traj, 
  std::vector<float>& seglens, 
  std::vector<float>& cumseglens, 
  std::vector<int>& seghits, 
  std::vector<float> dthetaLin, 
  std::vector<float> dthetaPoly, 
  double sigma3p, 
  int pid, 
  double sigma) const {
  if (verbose_ > 1) cout << "starting c2 function fit to find MCS momentum..." << endl; 
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
    if (verbose_ > 1) cout << "test momentum [MeV/c] = " << p_test << endl; 
    double c2 = C2Function(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, p_test);
    if (verbose_ > 1) cout << "c2 function = " << c2 << endl; 
    wmom[jMom] = p_test / 1000.; 
    wc2[jMom] = c2; 
    wsmom[jMom] = 1;
    wsigma[jMom] = sigma * c2; 
    if (firstValid < 0.01 && c2 > 0.01) firstValid = jMom;
    if (c2 > 0.) allZero = false;
    jMom++; }

  //return null is c2 function is null for all values of p test
  if (allZero) {
    if (verbose_ > 1) cout << "c2 function is zero everywhere! end fit" << endl;
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
    rsigma[jp - firstValid] = wsigma[jp];  }
  if (verbose_ > 1) {
    cout << "momentum [MeV/c] = "; for (auto i : rmom) cout << i << ' '; cout << endl;
    cout << "c2 function = "; for (auto i : rc2) cout << i << ' '; cout << endl;
    cout << "momentum error [MeV/c] = "; for (auto i : rsmom) cout << i << ' '; cout << endl;
    cout << "c2 function error = "; for (auto i : rsigma) cout << i << ' '; cout << endl; }
  
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
  if (verbose_ > 1) cout << "p_min [GeV/c] = " << p_min << " p_max [GeV/c] = " << p_max << endl;
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
  if (verbose_ > 1) {
    cout << "alpha = " << alpha << " and error on alpha = " << dalpha << endl;
    cout << "beta = " << beta << " and error on beta = " << dbeta << endl;
    cout << "best momentum [GeV/c] = " << p_best << " and error on momentum = " << p_error << endl; }
  auto dalpham = fabs(alpha + 0.5);
  auto dalphap = fabs(alpha - 1.);
  auto dbetam = fabs(beta - 0.);
  auto dbetap = fabs(beta - 5.);
  float lim = 0.001;

  if (dalpham < lim || dalphap < lim || dbetam < lim || dbetap < lim) {
    if (verbose_ > 1) {
      cout << "best momentum " << p_best << " but fit parameters at limit" << endl;
      cout << "computing C2NaiveFit" << endl; }
    return C2NaiveFit(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, sigma); }
  else if (p_best < p_min) {
    if (verbose_ > 1) {
      cout << "best momentum " << p_best << " less than p_min " << p_min << endl;
      cout << "computing C2NaiveFit" << endl; }
    return C2NaiveFit(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, sigma); }
  else if (p_best > p_max) {
    if (verbose_ > 1) {
      cout << "best momentum " << p_best << " greater than p_max " << p_max << endl;
      cout << "computing C2NaiveFit" << endl; }
    return C2NaiveFit(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, sigma); }

  return ScanResult(
    p_best, p_error, p_min, p_max,
    alpha, dalpha, beta, dbeta,
    rmom, rc2); }

//perform fit of likelihood function to determine MCS momentum
void TrajectoryMCSFitterICARUS::LikelihoodFit(
  const recob::TrackTrajectory& traj, 
  std::vector<float>& seglens, 
  std::vector<float>& cumseglens, 
  std::vector<int>& seghits, 
  std::vector<float> dthetaLin, 
  std::vector<float> dthetaPoly, 
  double sigma3p, 
  int pid, 
  double& best_p, 
  double& error_p) const {
  if (verbose_ > 1) cout << "starting likelihood function fit to find MCS momentum..." << endl;
  //initalize best value index, best log likelihood
  int best_idx = -1;
  double best_logL = std::numeric_limits<double>::max();
  double step = 10;
  //define vector vlogL, that will be filled with log likelihood values
  std::vector<float> vlogL;

  //iterate over p_test, compute log likelihood values
  for (double p_test = pMin_; p_test <= pMax_; p_test += step) {
    double logL = LikelihoodFunction(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, p_test, false);
    if (verbose_ > 1) cout << "test momentum for log likelihood [MeV/c] = " << p_test << " log likelihood = " << logL << endl; 
    if (logL < best_logL) {
      best_p = p_test;
      best_logL = logL;
      best_idx = vlogL.size(); }
    //add current loglikelihood to vector vlogL
    vlogL.push_back(logL); }

  //compute uncertainty from left side scan
  double lunc = -1.0;
  //check if best value index is not the first
  if (best_idx > 0) {
    //iterate for every index on the left of best value index
    for (int j = best_idx - 1; j >= 0; j--) {
      //define difference between best value loglikelood and some loglikelihood on the left
      double dLL = vlogL[j] - vlogL[best_idx];
      //define left uncertainty if this difference is less than 0.5
      if (dLL < 0.5) lunc = (best_idx - j) * step; 
      else break; } }

  //compute uncertainty from right side scan
  double runc = -1.0;
  //check if best value index is not the last
  if (best_idx < int(vlogL.size() - 1)) {  
    //iterate for every index on the right of best value index
    for (unsigned int j = best_idx + 1; j < vlogL.size(); j++) {
      //define difference between best value loglikelood and some loglikelihood on the right
      double dLL = vlogL[j] - vlogL[best_idx];
      //define right uncertainty if this difference is less than 0.5
      if (dLL < 0.5) runc = (j - best_idx) * step;
      else break; } }

  //print best value of momentum and uncertainty
  best_p = best_p/1000.; error_p = max(lunc,runc)/1000.;
  if (verbose_ > 1) cout << "best_p from likelihood [GeV/c] = " << best_p << "; uncertainty [GeV/c] = " << error_p << endl; 

  return; }

const TrajectoryMCSFitterICARUS::ScanResult TrajectoryMCSFitterICARUS::C2NaiveFit(
  const recob::TrackTrajectory& traj, 
  std::vector<float>& seglens, 
  std::vector<float>& cumseglens, 
  std::vector<int>& seghits, 
  std::vector<float> dthetaLin, 
  std::vector<float> dthetaPoly, 
  double sigma3p, 
  int pid, 
  double sigma) const {
  if (verbose_ > 1) cout << "starting c2 function fit to find MCS momentum..." << endl; 
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
    if (verbose_ > 1) cout << "test momentum [MeV/c] = " << p_test << endl; 
    double c2 = C2Function(traj, seglens, cumseglens, seghits, dthetaLin, dthetaPoly, sigma3p, pid, p_test);
    if (verbose_ > 1) cout << "c2 function = " << c2 << endl;
    wmom[jMom] = p_test / 1000.; 
    wc2[jMom] = c2; 
    if (wc2[jMom] > 0.) {
      if (verbose_ > 1) cout << "index = " << jMom << " p_test [GeV/c] = " << wmom[jMom] << " C2(p_test) = " << wc2[jMom] << endl; }
    wsmom[jMom] = 1;
    wsigma[jMom] = sigma * c2; 
    if (firstValid < 0.01 && c2 > 0.01) firstValid = jMom;
    if (c2 > 0.) allZero = false;
    jMom++; }

  //define minimum momentum and maximum momentum of fit 
  double p_min = wmom[firstValid]; double p_max = pMax_ / 1000.; 

  //return null is c2 function is null for all values of p test
  if (allZero) {
    if (verbose_ > 1) cout << "c2 function is zero everywhere! end fit" << endl;
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

  //search for first momentum whose c2 is greater than 1
  int over = -1;
  for (unsigned int jp = 0; jp < wmom.size(); jp++) {
    if (wc2[jp] > 1) {over = jp; break; } }
  auto p1 = wmom[over]; auto c2p1 = wc2[over];
  if (verbose_ > 1) cout << "over = " << over << " p1 [GeV/c] = " << p1 << " C2(p1) = " << c2p1 << endl;

  //search for first momentum whose c2 is less than 1
  int under = -1;
  for (unsigned int jp = wmom.size() - 1; jp >= 0; jp--) {
    if(wc2[jp] < 1) {under = jp; break; } }
  auto p2 = wmom[under]; auto c2p2 = wc2[under];
  if (verbose_ > 1) cout << "under = " << under << " p2 [GeV/c] = " << p2 << " C2(p2) = " << c2p2 << endl;
  
  //compute best p and error p with a naive fit
  double p_best; double p_error;
  if (under == -1) { 
    p_best = wmom[wmom.size()-1]; p_error = 0; }
  if (over == -1) { 
    p_best = wmom[0]; p_error = 0; }
  if (under >= 0 && over >= 0) {
    p_best = p2 + (p1-p2) * (1-c2p2) / (c2p1-c2p2);
    p_error = 0.5 * (p1-p2); }

  //return parameters alpha, beta with associated errors as null
  double alpha = 0.; double beta = 0.; 
  double dalpha = 0.; double dbeta = 0.; 

  //return best p, error p 
  if (verbose_ > 1) cout << "naive best momentum [GeV/c] = " << p_best << " and error on momentum = " << p_error << endl;
  if (p_best < p_min) {
    if (verbose_ > 1) cout << "naive best momentum " << p_best << " less than p_min " << p_min << endl; }
  else if (p_best > p_max) {
    if (verbose_ > 1) cout << "naive best momentum " << p_best << " greater than p_max " << p_max << endl; }
  return ScanResult(
    p_best, p_error, p_min, p_max,
    alpha, dalpha, beta, dbeta,
    rmom, rc2); }

//fill linear matrix with off-diag MCS terms
const void TrajectoryMCSFitterICARUS::FillOffDiagMCSLin(
  TMatrixDSym& matOD, 
  TMatrixD matdiag, 
  int jp) const {
  double whp = 0.174;
  if (jp >= 1) {
    matOD(jp, jp-1) += sqrt(matdiag(jp, jp) * matdiag(jp-1, jp-1)) * whp; 
    matOD(jp-1, jp) += sqrt(matdiag(jp, jp) * matdiag(jp-1, jp-1)) * whp; } 
}

//fill poligonal matrix with off-diag MCS terms
const void TrajectoryMCSFitterICARUS::FillOffDiagMCSPoly(
  TMatrixDSym& matOD, 
  TMatrixD matdiag, 
  int jp) const {
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
const void TrajectoryMCSFitterICARUS::FillOffDiagErrLin(
  TMatrixDSym& matOD, 
  TMatrixD matdiag, 
  int jp) const {
  double whp = -1./2.;
  if (jp >= 1) {
    matOD(jp, jp-1) += matdiag(jp, jp) * whp;
    matOD(jp-1, jp) += matdiag(jp, jp) * whp; } 
}

//fill poligonal matrix with off-diag error terms
const void TrajectoryMCSFitterICARUS::FillOffDiagErrPoly(
  TMatrixDSym& matOD, 
  TMatrixD matdiag, 
  int jp) const {
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
const void TrajectoryMCSFitterICARUS::FillMCSMixedTerms(
  TMatrixD& matMix, 
  TMatrixD matLin, 
  TMatrixD matPoly) const {
  double whp = 0.174;
  double whpp = 0.067;
  for (int jp = 0; jp < matPoly.GetNrows(); jp++) {
    matMix(jp, jp) += sqrt(matLin(jp, jp) * matPoly(jp, jp)) * whp;
    if (jp >= 1) matMix(jp-1, jp) += sqrt(matLin(jp, jp) * matPoly(jp-1, jp-1)) * whpp;
    if (jp < matPoly.GetNrows() - 1) matMix(jp+1, jp) += sqrt(matLin(jp, jp) * matPoly(jp+1, jp+1)) * whp; 
    if (jp < matPoly.GetNrows() - 2) matMix(jp, jp+2) += sqrt(matLin(jp, jp) * matPoly(jp+2, jp+2)) * whpp; } 
}

//clean covariance matrix of jtail row and column
TMatrix TrajectoryMCSFitterICARUS::CleanCovariance(
  TMatrix cov, 
  int jtail) const {
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

//gives the plane of index point: 0 for Induction-1, 1 for Induction-2, 2 for Collection
unsigned int TrajectoryMCSFitterICARUS::Plane(
  size_t index) const {
  return get<1>(pdata[index])->WireID().Plane;
}


//check if index point is in plane = p
bool TrajectoryMCSFitterICARUS::IsInPlane(
  size_t index) const {
  if (Plane(index) == planeMode_) return true;
  else return false; 
}

//gives the physical tpc of index point: 0 for TPC 0-1 (x < cathode), 1 for TPC 2-3 (x > cathode)
unsigned int TrajectoryMCSFitterICARUS::TPC(
  size_t index) const {
  return get<1>(pdata[index])->WireID().TPC / 2;
}

//check tpc of first valid point
unsigned int TrajectoryMCSFitterICARUS::FirstTPC(
  const recob::TrackTrajectory& traj) const {
  return TPC(traj.FirstValidPoint());
}

//check tpc of last valid point
unsigned int TrajectoryMCSFitterICARUS::LastTPC(
  const recob::TrackTrajectory& traj) const {
  return TPC(traj.LastValidPoint());
}

//check if index point is in tpc = t
bool TrajectoryMCSFitterICARUS::IsInTPC(
  size_t index, 
  unsigned int t) const {
  if (TPC(index) == t) return true;
  else return false; 
}

//gives 2D coordinate of index point: x = drift, y = wire
Vector_t TrajectoryMCSFitterICARUS::hit2D(
  size_t index) const {
  const recob::Hit& h = *get<1>(pdata[index]);
  const geo::WireID wid = h.WireID();
  if (wid.Plane != planeMode_) return Vector_t();

  //drift coordinate, corrected with t0
  const double shift = geom_->TPC(wid).DriftDir().X() * detProp_->DriftVelocity() * T0_;
  const double x = detProp_->ConvertTicksToX(h.PeakTime(), wid.Plane, wid.TPC, wid.Cryostat) + shift;

  //wire coordinate, with wire = n·c, n normal to wire, c wire center (absorbs offset 2536)
  const auto& plane = wr_->Plane(wid);
  const auto n = plane.GetIncreasingWireDirection();
  const auto c = wr_->Wire(wid).GetCenter();
  const double y = n.Y()*c.Y() + n.Z()*c.Z();
  return Vector_t(x, y, 0);
}

//find 1D (drift) length of track [cm]
double TrajectoryMCSFitterICARUS::length1D(
  const recob::TrackTrajectory traj) const {
  //define index of first and last valid point
  size_t index = traj.FirstValidPoint();
  auto frstIndex = -1; 
  auto lastIndex = -1; 

  //find drift length of track in chosen view
  while (index < traj.LastValidPoint()) {
    if (IsInPlane(index) && IsInTPC(index, LastTPC(traj))) {
      if (frstIndex == -1) frstIndex = index;
      lastIndex = index; }
    index = traj.NextValidPoint(index + 1); }
  //define first and last valid point
  Vector_t frstPoint = hit2D(frstIndex);
  Vector_t lastPoint = hit2D(lastIndex);
  //compute drift coordinate difference
  double dx = abs(frstPoint.X() - lastPoint.X());
  //return drift length 
  return dx;
}
  
//find 2D (drift-wire) length of track [cm]
double TrajectoryMCSFitterICARUS::length2D(
  const recob::TrackTrajectory traj) const {
  //define index of first and last valid point
  size_t index = traj.FirstValidPoint();
  auto frstIndex = -1; 
  auto lastIndex = -1; 

  //compute drift-wire length of track in chosen view 
  while (index < traj.LastValidPoint()) {
    if (IsInPlane(index) && IsInTPC(index, LastTPC(traj))) {
      if (frstIndex == -1) frstIndex = index;
      lastIndex = index; }
    index = traj.NextValidPoint(index + 1); }
  //define first and last valid point
  Vector_t frstPoint = hit2D(frstIndex);
  Vector_t lastPoint = hit2D(lastIndex);
  //compute drift and wire coordinate differences
  double dx = abs(frstPoint.X() - lastPoint.X());
  double dy = abs(frstPoint.Y() - lastPoint.Y());
  //return drift-wire length 
  return sqrt(pow(dx, 2) + pow(dy, 2)); 
}
  
//find 3D length of track [cm]
double TrajectoryMCSFitterICARUS::length3D(
  const recob::TrackTrajectory traj) const {
  //define index of first and last valid point
  size_t index = traj.FirstValidPoint();
  auto frstIndex = -1; 
  auto lastIndex = -1; 

  //find 3D length of track in chosen view 
  while (index < traj.LastValidPoint()) {
    if (IsInPlane(index) && IsInTPC(index, LastTPC(traj))) {
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
  return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2)); 
}

//compute directional factor for thetamcs in c2 function computation
double TrajectoryMCSFitterICARUS::ThetaMCSFactor(
  const recob::TrackTrajectory& traj) const {
  //recall length2D as length of track on chosen view (Ind1, Ind2, Coll)
  double L2D = length2D(traj);
  //recall length3D as length of track in three dimensions
  double L3D = length3D(traj);
  //return directional factor 1/cos(alpha) as ratio between L3D and L2D
  return L3D / L2D;
}

//compute expected angle (thetamcs) for a certain momentum
double TrajectoryMCSFitterICARUS::ThetaExpected(
  const recob::TrackTrajectory& traj, 
  double p, 
  double beta, 
  double L) const {
  //define Highland parameters S2 [MeV], epsilon, X0 [cm]
  double S2 = 13.6; double eps = 0.038; double X0 = 14.;

  //compute theta MCS [mrad] just as claimed by Highland formula
  double thetamcs = S2 / (p * beta) * sqrt(L/X0) * (1 + eps * log(L/X0)) * ThetaMCSFactor(traj);
  return thetamcs;
}

//compute directional factor for thetaerr in c2 function computation
double TrajectoryMCSFitterICARUS::ThetaErrFactor(
  const recob::TrackTrajectory& traj) const {
  //recall L1D as length of track along drift direction
  double L1D = length1D(traj);
  //recall L2D as length of track on chosen view (I1, I2, C)
  double L2D = length2D(traj);

  //return directional factor sin(beta) = sqrt(1 - cos(gamma)) with cos(gamma) as ratio between L1D and L2D
  return sqrt(1 - pow(L1D, 2)/pow(L2D, 2));
}

//check if angles vector contains a minimum number of angles > 0 to avoid crash
void TrajectoryMCSFitterICARUS::ThetaCheck(
  vector<float> dthetaLin, 
  vector<float> dthetaPoly, 
  bool& checkLin, 
  bool& checkPoly, 
  unsigned int& firstseg, 
  unsigned int& lastseg) const {
  bool first = false;
  auto countLin = 0;
  for (unsigned int j = 0; j < dthetaLin.size(); ++j) {
    if (dthetaLin[j] > 0.) {
      if (!first) { firstseg = j; first = true; }
      lastseg = j + 1;
      countLin++;
      if (countLin >= minNAngs_) checkLin = true; }
    else { 
      if (fullTrackMode_ == 0 || fullTrackMode_ == 1) { first = false; countLin = 0; } } }

  auto countPoly = 0;
  for (unsigned int j = 0; j < dthetaPoly.size(); ++j) {
    if (dthetaPoly[j] > 0.) {
      countPoly++;
      if (countPoly >= minNAngs_ - 1) checkPoly = true; }
    else { 
      if (fullTrackMode_ == 0 || fullTrackMode_ == 1) { countPoly = 0; } } }
}

//check to exclude hit near cathode
bool TrajectoryMCSFitterICARUS::CathodeCheck(
  const recob::TrackTrajectory& traj, 
  size_t index) const {
  art::Ptr<recob::Hit> hit = get<1>(pdata[index]);
  unsigned int c = hit->WireID().Cryostat;
  double cathode = geom_->TPC(geo::TPCID{c, 0}).GetCathodeCenter().X();
  double point = traj.LocationAtPoint(index).X();
  double tolerance = 15;
  if (point < cathode - tolerance || point > cathode + tolerance) return true;
  else return false;
}