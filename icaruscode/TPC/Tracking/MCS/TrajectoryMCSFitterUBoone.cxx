#include "TrajectoryMCSFitterUBoone.h"

#include "lardata/RecoObjects/TrackState.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

using namespace std;
using namespace trkf;
using namespace recob::tracking;
using PointFlags_t=recob::TrajectoryPointFlags;
using Flags_t=std::vector<PointFlags_t>;

recob::MCSFitResult TrajectoryMCSFitterUBoone::fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst) const {
  //notify algorithm has started
  cout << "starting uboone fitting... " << endl;

  //define vector breakpoints, that will be filled with integer indexes of trajectory segmentation
  vector<size_t> breakpoints;
  //define vector breakpointsgood, that will be filled with boolean indexes of trajectory segmentation
  vector<bool> breakpointsgood;
  //define vector segradlengths, that will be filled with segment lengths in units of radiation length
  vector<float> segradlengths;
  //define vector cumseglens, that will be filled with cumulative segment lengths
  vector<float> cumseglens;
  //divide trajectory in segments and fill vectors defined above
  breakTrajInSegments(traj, breakpoints, segradlengths, cumseglens, breakpointsgood, cutMode_, cutLength_);

  //return empty recob::MCSFitResult if there are less than 2 segments
  if (segradlengths.size() < 2) return recob::MCSFitResult();
  //define vector dtheta, that will be filled with scattering angles between adjacent segments
  vector<float> dtheta;
  //define vectors pcdir, that will be filled with directions of two adjacent segments
  Vector_t pcdir0;
  Vector_t pcdir1;

  //print cut mode
  if (cutMode_ == 0) {cout << "type of cut = no cut" << endl;}
  if (cutMode_ == 1) {cout << "type of cut = final cut" << endl;}
  if (cutMode_ == 2) {cout << "type of cut = initial cut" << endl;}
  cout << " " << endl;
  
  //iterate for every pair of adjacent segments
  for (unsigned int p = 0; p < segradlengths.size(); p++) {
    //print iteration number
    cout << "iteration number " << p << endl;
    
    //perform linear fit of segment and memorize its direction in pcdir1
    linearRegression(traj, breakpoints[p], breakpoints[p + 1], pcdir1);
    //ignore first iteration because it is used to compute direction of first segment
    if (p == 0) {cout << "first iteration, scattering angle impossible to compute" << endl;}

    //start to compute scattering angles from second iteration
    if (p > 0) {
      //ignore current iteration if segment lengths have anomalous values
      if (segradlengths[p] < -100. || segradlengths[p - 1] < -100.) {
        cout << "WARNING! invalid segment length found!" << endl;
        dtheta.push_back(-9.);
      } 
      //ignore current iteration if segmentation index is wrong
      if (!breakpointsgood[p]) {
        cout << "WARNING! breakpoint less than 25 cm from cathode!" << endl;
        dtheta.push_back(-3.);
      }
      else {
        //print positions of first and last point of segment
        auto startpos = traj.LocationAtPoint(breakpoints[p-1]);
        auto endpos = traj.LocationAtPoint(breakpoints[p]);
        cout << "segment starts at " << startpos << endl;
        cout << "segment ends at " << endpos << endl;
        //print current segment directions and lengths, and cumulative length
        cout << "initial direction, x = " << pcdir0.X() << " y = " << pcdir0.Y() << " z = " << pcdir0.Z() << endl;
        cout << "final direction, x = " << pcdir1.X() << " y = " << pcdir1.Y() << " z = " << pcdir1.Z() << endl;
        cout << p << "-th segment length [cm] = " << 14. * segradlengths[p - 1] << endl;
        cout << p + 1 << "-th segment length [cm] = " << 14. * segradlengths[p] << endl;
        cout << "cumulative length [cm] = " << cumseglens[p] << endl;

        //compute scalar product of two directions, i.e. cosine of scattering angle
        const double cosval = pcdir0.X() * pcdir1.X() + pcdir0.Y() * pcdir1.Y() + pcdir0.Z() * pcdir1.Z();
        //compute scattering angle from its cosine, in units of mrad
        double dt = 1000. * acos(cosval);
        //print current scattering angle between the two segments
        cout << p << "-th scattering angle [mrad] = " << dt << endl;

        //if no cut is performed
        if (cutMode_ == 0) {
          //fill vector dtheta with current scattering angle
          dtheta.push_back(dt);
        }

        //if cut at the end of the track is performed
        else if (cutMode_ == 1) {
          //ignore last 3 scattering angles
          if (p >= segradlengths.size() - 3) {
            cout << "scattering angle not considered for the fit" << endl;
            dtheta.push_back(-1.);
          }
          else {
            //fill vector dtheta with current scattering angle
            dtheta.push_back(dt);
          }
        }
        
        //if cut on the first part of the track is performed
        else if (cutMode_ == 2) {
          //ignores from 8th scattering angle
          //if (p > 7) {
          //  cout << "scattering angle not considered for the fit" << endl;
          //  dtheta.push_back(-2.);
          //} else {
          //fill vector dtheta with current scattering angle
          //  dtheta.push_back(dt);
          //}
          
          //ignore last 3 scattering angles
          if (p >= segradlengths.size() - 3) {
            cout << "scattering angle not considered for the fit" << endl;
            dtheta.push_back(-1.);
          }
          else {
            //fill vector dtheta with current scattering angle
            dtheta.push_back(dt);
          }
        }
      }
    }
    //update first segment direction for next iteration
    pcdir0 = pcdir1;
    cout << " " << endl;
  }

  //print vector of segment lengths
  cout << "segment lengths [cm] = ";
  for (auto i : segradlengths) {
    cout << 14. * i << ' ';
  } 
  cout << endl;

  //print vector of cumulative segment lengths
  cout << "cumulative segment length [cm] = ";
  for (auto i : cumseglens) {
    cout << i << ' ';
  }
  cout << endl;

  //remove hits from 1 of the 2 tpcs
  if (cutMode_ == 1) {
    //find last angle not computed due to proximity to cathode
    auto lastAngleNearToCat = std::find(dtheta.rbegin(), dtheta.rend(), -3.);
    
    //if no angles are not computed due to proximity to cathode, break; else, continue
    if (lastAngleNearToCat == dtheta.rend()) {
      dtheta=dtheta;
    } else {
      //get index of last angle not computed due to proximity to cathode
      int lastAngleNearToCatIndex = dtheta.size() - 1 - std::distance(dtheta.rbegin(), lastAngleNearToCat);

      //count how many angles there are after last angle not computed
      int positiveCount = 0;
      for (size_t i = lastAngleNearToCatIndex + 1; i < dtheta.size(); ++i) {
        if (dtheta[i] > 0.) {
          positiveCount++;
        } else {
          break;
        }
      }

      //verify that there are at least 7 angles after proximity to cathode
      if (positiveCount >= 7) {
        //remove from dtheta angles belonging to the other TPC (before proximity to cathode)
        for (int i = lastAngleNearToCatIndex - 1; i >= 0; --i) {
          if (dtheta[i] > 0.) {
            dtheta[i] = -4.;
          } else {
            continue; 
          }
        }
      } else {
        //remove from dtheta all angles
        for (float &val : dtheta) {
          if (val > 0.) {
            val = -4.;
          }
        }
      }
    }
  } else if (cutMode_ == 2) {
    //find first 7 consecutive angles correctly computed
    size_t startIndex = dtheta.size();
    int consecutiveCount = 0;
    for (size_t i = 0; i < dtheta.size(); ++i) {
      if (dtheta[i] > 0.) {
        consecutiveCount++;
        if (consecutiveCount == 7) {
          //get index of first of 7 consecutive angles correctly computed
          startIndex = i - 6;
          break;
        }
      } else {
        consecutiveCount = 0;
      }
    }
    
    //verify that there are 7 consecutive angles correctly computed
    if (consecutiveCount == 7) {
      for (size_t i = 0; i < dtheta.size(); ++i) {
        if (dtheta[i] > 0. && (i < startIndex || i >= startIndex + 7)) {
          dtheta[i] = -5.;
        }
      }
    } else {
      //remove from dtheta all angles
      for (float &val : dtheta) {
        if (val > 0.) {
          val = -5.;
        }
      }
    }
  }

  //print vector of scattering angles
  cout << "scattering angles [mrad] = ";
  for (auto i : dtheta) {
    cout << i << ' ';
  }
  cout << endl;

  //define vectors cumLen to memorize cumulative segment lengths, forward and backward
  vector<float> cumLenFwd;
  vector<float> cumLenBwd;
  //iterate for every segment to fill vectors cumLen
  for (unsigned int i = 0; i < cumseglens.size() - 2; i++) {
    //fill vector cumLenFwd with cumulative lengths from first to secondlast segment
    cumLenFwd.push_back(cumseglens[i]);
    //fill vector cumLenFwd with cumulative lengths from secondlast to first segment
    cumLenBwd.push_back(cumseglens.back() - cumseglens[i + 2]);
  }

  //perform a likelihood scan (forward; backward is ignored)
  const ScanResult fwdResult = doLikelihoodScan(dtheta, segradlengths, cumLenFwd, true, momDepConst, pid);

  //return result of MCS fit with interesting parameters
  return recob::MCSFitResult(pid, fwdResult.p, fwdResult.pUnc, fwdResult.logL, 0., 0., 0., segradlengths, dtheta);
}

void TrajectoryMCSFitterUBoone::breakTrajInSegments(const recob::TrackTrajectory& traj, vector<size_t>& breakpoints, vector<float>& segradlengths, vector<float>& cumseglens, vector<bool>& breakpointsgood, int cutMode, float cutLength) const {
  //set fiducial volume
  std::vector<geo::BoxBoundedGeo> fiducialVolumes;
  fiducialVolumes = setFiducialVolumes();
  //set excluded volumes
  std::vector<geo::BoxBoundedGeo> excludeVolumes;
  excludeVolumes = setExcludeVolumes();

  //define trajectory length
  double trajlen = traj.Length();
  //print trajectory length
  cout << "track length [cm] = " << trajlen << endl;

  //define segment length in function of trajectory length being less or greater than minimum length
  const double thisSegLen = (trajlen > (segLen_ * minNSegs_) ? segLen_ : trajlen / double(minNSegs_));
  //print required segment length, considered segment length, number of segments
  cout << "required segment length [cm] = " << segLen_ << endl;
  cout << "considered segment length [cm] = " << thisSegLen << endl;
  cout << "number of segments = " << max(minNSegs_, int(trajlen / segLen_)) << endl;

  //define inverse of Argon radiation length X0 = 14 cm
  constexpr double lar_radl_inv = 1. / 14.0;

  //first segment has zero cumulative length from previous segments
  cumseglens.push_back(0.);
  //initialize current segment length
  double thislen = 0.;

  //find first valid point index and add it to vector breakpoints
  auto nextValid = traj.FirstValidPoint();
  breakpoints.push_back(nextValid);
  //compute position of first valid point
  auto pos0 = traj.LocationAtPoint(nextValid);
  //determine if first valid point is inside fiducial volume and add it to vector breakpointsgood
  bool pos0good = isInVolume(fiducialVolumes, pos0) && !isInVolume(excludeVolumes, pos0);
  breakpointsgood.push_back(pos0good);
  
  //find next valid point index
  nextValid = traj.NextValidPoint(nextValid + 1);
  //initalize number of valid points in a segment
  int npoints = 0;
  
  //iterate over all valid points of trajectory
  while (nextValid != recob::TrackTrajectory::InvalidIndex) {
    //compute position of current valid point, then update segment length and position of next valid point
    auto pos1 = traj.LocationAtPoint(nextValid);
    thislen += ((pos1-pos0).R());
    pos0 = pos1;
    //increase number of valid points in a segment
    npoints++;

    //break into a new segment if current segment length is greater than required segment length
    if (thislen >= thisSegLen) {
      //determine if current valid point is inside fiducial volume and add it to vector breakpointsgood
      bool pos1good = isInVolume(fiducialVolumes, pos1) && !isInVolume(excludeVolumes, pos1);
      breakpointsgood.push_back(pos1good);
      //add current valid point index to vector breakpoints
      breakpoints.push_back(nextValid);

      //add segment length, in units of radiation length, to vector segradlengths
      if (npoints >= minHitsPerSegment_) {segradlengths.push_back(thislen * lar_radl_inv);}
      else {
        cout << "WARNING! invalid number of hits per segment!" << endl;
        continue;
      }

      //add cumulative segment length to vector cumseglens
      cumseglens.push_back(cumseglens.back() + thislen);

      //initialize segment length and number of valid points for next iteration
      thislen = 0.;
      npoints = 0;
    }

    //find next valid point index for next iteration
    nextValid = traj.NextValidPoint(nextValid + 1);
  }

  //then add last segment
  if (thislen > 0.) {
    //compute position of last valid point
    auto endpointpos = traj.LocationAtPoint(nextValid);
    //determine if last valid point is inside fiducial volume and add it to vector breakpointsgood
    bool endpointposgood = isInVolume(fiducialVolumes, endpointpos) && !isInVolume(excludeVolumes, endpointpos);
    breakpointsgood.push_back(endpointposgood);

    //add last valid point index to vector breakpoints
    breakpoints.push_back(traj.LastValidPoint() + 1);

    //add segment length, in units of radiation length, to vector segradlengths
    segradlengths.push_back(thislen * lar_radl_inv);

    //add cumulative segment length to vector cumseglens
    cumseglens.push_back(cumseglens.back() + thislen);
  }
  return;
}

const TrajectoryMCSFitterUBoone::ScanResult TrajectoryMCSFitterUBoone::doLikelihoodScan(std::vector<float>& dtheta, std::vector<float>& seg_nradlengths, std::vector<float>& cumLen, bool fwdFit, bool momDepConst, int pid) const {
  //initalize best value index as last, loglikelihood as infinity, momentum as -1
  int best_idx = -1;
  double best_logL = std::numeric_limits<double>::max();
  double best_p = -1.0;
  //define vector vlogL, that will be filled with loglikelihood values
  std::vector<float> vlogL;

  //iterate for different values of momentum to perform a likelihood scan
  for (double p_test = pMin_; p_test <= pMax_; p_test += pStep_) {
    //compute loglikelihood in function of momentum
    double logL = mcsLikelihood(p_test, angResol_, dtheta, seg_nradlengths, cumLen, fwdFit, momDepConst, pid);
    //substitute best value if current loglikelihood is minimum
    if (logL < best_logL) {
      best_p = p_test;
      best_logL = logL;
      best_idx = vlogL.size();
    }
    //add current loglikelihood to vector vlogL
    vlogL.push_back(logL);
  }
  
  //compute uncertainty from left side scan
  double lunc = -1.0;
  //check if best value index is not the first
  if (best_idx > 0) {
    //iterate for every index on the left of best value index
    for (int j = best_idx - 1; j >= 0; j--) {
      //define difference between best value loglikelood and some loglikelihood on the left
      double dLL = vlogL[j] - vlogL[best_idx];
      //define left uncertainty if this difference is less than 0.5
      if (dLL < 0.5) {lunc = (best_idx - j) * pStep_;} 
      else break;
    }
  }

  //compute uncertainty from right side scan
  double runc = -1.0;
  //check if best value index is not the last
  if (best_idx < int(vlogL.size() - 1)) {  
    //iterate for every index on the right of best value index
    for (unsigned int j = best_idx + 1; j < vlogL.size(); j++) {
      //define difference between best value loglikelood and some loglikelihood on the right
      double dLL = vlogL[j] - vlogL[best_idx];
      //define right uncertainty if this difference is less than 0.5
      if (dLL < 0.5) {runc = (j - best_idx) * pStep_;} 
      else break;
    }
  }

  //print best value of momentum and uncertainty
  cout << "best_p [GeV/c] = " << best_p << "; uncertainty [GeV/c] = " << max(lunc,runc) << endl;
  cout << " " << endl;

  //return best values of momentum, uncertainty and loglikelihood
  return ScanResult(best_p, std::max(lunc,runc), best_logL);
}

void TrajectoryMCSFitterUBoone::linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir) const {
  //initialize number of valid points in the segment
  int npoints = 0;
  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;

  //determine index of first valid point
  size_t nextValid = firstPoint;
  //iterate until last valid point of the segment
  while (nextValid < lastPoint) {
    //add position of current valid point to vector middlePointCalc
    middlePointCalc.add(traj.LocationAtPoint(nextValid));
    //determine index of next valid point
    nextValid = traj.NextValidPoint(nextValid + 1);
    //increase number of valid points in the segment
    npoints++;
  }

  //determine position of average point of segment from position of valid points in the segment
  const auto avgpos = middlePointCalc.middlePoint();
  //define normalization factor as inverse of number of valid points
  const double norm = 1. / double(npoints);
  //define covariance matrix as symmetric matrix 3x3
  TMatrixDSym m(3);

  //determine again index of first valid point
  nextValid = firstPoint;
  //iterate until last valid point of the segment
  while (nextValid < lastPoint) {
    //determine position of current valid point
    auto p = traj.LocationAtPoint(nextValid);
    //compute coordinate differences between current valid point and average point
    const double xxw0 = p.X() - avgpos.X();
    const double yyw0 = p.Y() - avgpos.Y();
    const double zzw0 = p.Z() - avgpos.Z();
    //update covariance matrix values with normalized values of coordinate differences above
    m(0, 0) += xxw0 * xxw0 * norm;
    m(0, 1) += xxw0 * yyw0 * norm;
    m(0, 2) += xxw0 * zzw0 * norm;
    m(1, 0) += yyw0 * xxw0 * norm;
    m(1, 1) += yyw0 * yyw0 * norm;
    m(1, 2) += yyw0 * zzw0 * norm;
    m(2, 0) += zzw0 * xxw0 * norm;
    m(2, 1) += zzw0 * yyw0 * norm;
    m(2, 2) += zzw0 * zzw0 * norm;
    //determine index of next valid point
    nextValid = traj.NextValidPoint(nextValid + 1);
  }

  //compute eigenvalues and eigenvectors of covariance matrix
  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //initialize maximum eigenvalue index and value
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  //iterate over eigenvalues
  for (int i = 1; i < 3; ++i) {
    //update maximum eigenvalue if current is greater
    if (eigenval(i) > maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  }

  //update segment direction with eigenvector associated to maximum eigenvalue
  pcdir = Vector_t(eigenvec(0, maxevalidx), eigenvec(1, maxevalidx), eigenvec(2, maxevalidx));

  //reverse segment direction if current direction is opposite to trajectory direction in the first valid point
  if (traj.DirectionAtPoint(firstPoint).Dot(pcdir) < 0.) pcdir *= -1.;
}

double TrajectoryMCSFitterUBoone::mcsLikelihood(double p, double theta0x, std::vector<float>& dthetaij, std::vector<float>& seg_nradl, std::vector<float>& cumLen, bool fwd, bool momDepConst, int pid) const {
  //determine beginning, end, incremental index for iteration if forward or backward 
  const int beg = (fwd ? 0 : (dthetaij.size() - 1));
  const int end = (fwd ? dthetaij.size() : -1);
  const int incr = (fwd ? +1 : -1);

  //determine particle mass from particle id, and compute squared mass
  const double m = mass(pid);
  const double m2 = m * m;
  //compute total energy of particle from E2=p2+m2
  const double Etot = sqrt(p * p + m2);

  //compute fixed term of loglikelihood
  double const fixedterm = 0.5 * log(2.0 * M_PI);

  //initalize loglikelihood value
  double result = 0;
  //iterate over trajectory segments
  for (int i = beg; i != end; i += incr) {
    //ignore current iteration if scattering angle has negative value
    if (dthetaij[i] < 0) {continue;}

    //determine energy in function of energy loss model, and compute squared energy
    const double Eij = GetE(Etot, cumLen[i], m);
    const double Eij2 = Eij * Eij;
    //return infinity loglikelihood if energy at current segment is less than energy at rest
    if (Eij2 <= m2) {
      result = numeric_limits<double>::max();
      break;
    }

    //compute momentum and velocity of particle at current segment
    const double pij = sqrt(Eij2 - m2);
    const double beta = sqrt(1. - (m2 / Eij2));

    //define Highland term (first one is tuned as explained in https://arxiv.org/abs/1703.06187)
    constexpr double tuned_HL_term1 = 11.0038;
    constexpr double HL_term2 = 0.038;
    //compute MCS angle from Highland formula, first Highland term might depend on momentum
    const double tH0 = ( (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) / (pij * beta) ) * (1.0 + HL_term2 * log(seg_nradl[i])) * sqrt(seg_nradl[i]);
    //compute total MCS angle as quadrature sum of Highland and angular resolution
    const double rms = sqrt(2.0 * (tH0 * tH0 + theta0x * theta0x));
    //return infinity loglikelihood if total MCS angle is null
    if (rms == 0.0) {
      cout << "WARNING! RMS cannot be zero" << endl;
      return numeric_limits<double>::max();
    } 

    //compute ratio between scattering and MCS angles
    const double arg = dthetaij[i]/rms;
    //compute loglikelihood as sum of log, quadratic, fixed terms for every iteration
    result += ( std::log( rms ) + 0.5 * arg * arg + fixedterm);
  }
  //return final value of loglikelihood
  return result;
}

double TrajectoryMCSFitterUBoone::energyLossLandau(const double mass2, const double e2, const double x) const {
  //return null energy if travelled distance is negative or null
  if (x <= 0.) return 0.;

  //define inverse of squared ionization potential 1/I2, in units of MeV-2 
  constexpr double Iinv2 = 1. / (188.E-6 * 188.E-6);
  //define product of density rho, atomic number Z, inverse of mass number 1/A, in units of mol cm-3
  constexpr double matConst = 1.4 * 18. / 40.;
  //define electron mass, in units of MeV/c2
  constexpr double me = 0.511;
  //define K constant, in units of MeV mol−1 cm2
  constexpr double kappa = 0.307075;
  //define Landau parameter j
  constexpr double j = 0.200;

  //compute squared velocity beta2 and squared Lorentz factor gamma2
  const double beta2 = (e2 - mass2) / e2;
  const double gamma2 = 1. / (1.0 - beta2);

  //return most probable value of energy, according to Landau distribution, in units of GeV
  const double epsilon = 0.5 * kappa * x * matConst / beta2;
  return 0.001 * epsilon * (log(2. * me * beta2 * gamma2 * epsilon * Iinv2) + j - beta2);
}

double TrajectoryMCSFitterUBoone::energyLossBetheBloch(const double mass, const double e2) const {
  //define inverse of ionization potential 1/I, in units of MeV-1
  constexpr double Iinv = 1. / 188.E-6;
  //define product of density rho, atomic number Z, inverse of mass number 1/A, in units of mol cm-3
  constexpr double matConst = 1.4 * 18. / 40.;
  //define electron mass, in units of MeV/c2
  constexpr double me = 0.511;
  //define K constant, in units of MeV mol−1 cm2
  constexpr double kappa = 0.307075;

  //compute squared velocity beta2 and squared Lorentz factor gamma2
  const double beta2 = (e2 - mass * mass) / e2;
  const double gamma2 = 1. / (1.0 - beta2);

  //compute ratio between electron mass and particle mass
  const double massRatio = me / mass;
  //compute argument of logaritmic term 
  const double argument = (2. * me * gamma2 * beta2 * Iinv) / sqrt(1 + 2 * sqrt(gamma2) * massRatio + massRatio * massRatio);
  
  //initialize energy loss, in units of MeV cm-1
  double dedx = kappa * matConst / beta2;
  //return null energy if particle mass is null or if argument of logaritmic term is not greater than exp(beta2)
  if (mass == 0.0) return (0.0);
  if (argument <= exp(beta2)) return (0.0);
  else {
    dedx *= (log(argument) - beta2);
    //return null energy if this energy loss is negative
    if (dedx < 0.) return (0.0);
  }

  //return energy loss according to Bethe Bloch, in units of GeV cm-1
  return 0.001 * dedx;
}

double TrajectoryMCSFitterUBoone::GetE(const double initial_E, const double length_travelled, const double m) const {
  //compute energy if particle is considered MIP with constant energy loss
  if (eLossMode_ == 1) {
    //define energy loss for MIP muon in liquid Argon, in units of GeV cm-1
    constexpr double kcal = 0.002105;
    //return energy as difference between initial energy and energy loss times travelled distance, in units of GeV
    return (initial_E - kcal * length_travelled);
  }

  //define step size as ratio between travelled distance and an arbitraty number (ex 10)
  const double step_size = length_travelled / nElossSteps_;
  //define initial energy, in units of GeV and squared mass, in units of (GeV/c2)2
  double current_E = initial_E;
  const double m2 = m * m;

  //compute energy if particle is considered with variable energy loss
  for (auto i = 0; i < nElossSteps_; ++i) {
    //compute energy according to Bethe Bloch, negleting density correction
    if (eLossMode_ == 2) {
      double dedx = energyLossBetheBloch(m, current_E * current_E);
      current_E -= (dedx * step_size);
    }

    //compute energy according to Landau distribution
    else {
      current_E -= energyLossLandau(m2, current_E * current_E, step_size);
    }

    //return null energy if current energy is less or equal than energy at rest
    if (current_E <= m) {return 0.;}
  }

  //return current energy of particle in function of travelled distance, in units of GeV
  return current_E;
}

std::vector<geo::BoxBoundedGeo> TrajectoryMCSFitterUBoone::setFiducialVolumes() const {
  std::vector<geo::BoxBoundedGeo> fiducialVolumes;
  std::vector<std::vector<geo::BoxBoundedGeo>> TPCVolumes;
  std::vector<geo::BoxBoundedGeo> ActiveVolumes;
  
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  for (auto const &cryoID: geometry->Iterate<geo::CryostatID>()) {
    std::vector<geo::BoxBoundedGeo> thisTPCVolumes;
    for (auto const& tpc: geometry->Iterate<geo::TPCGeo>(cryoID)) {
      thisTPCVolumes.push_back(tpc.ActiveBoundingBox());
    }
    TPCVolumes.push_back(std::move(thisTPCVolumes));
  }
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: TPCVolumes) {
    double xMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double xMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double yMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double yMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double zMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();
    double zMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();
    ActiveVolumes.emplace_back(xMin, xMax, yMin, yMax, zMin, zMax);
  }

  double fidInsetMinX = fiducialVolumeInsets_[0];
  double fidInsetMaxX = fiducialVolumeInsets_[1];
  double fidInsetMinY = fiducialVolumeInsets_[2];
  double fidInsetMaxY = fiducialVolumeInsets_[3];
  double fidInsetMinZ = fiducialVolumeInsets_[4];
  double fidInsetMaxZ = fiducialVolumeInsets_[5];

  for (const geo::BoxBoundedGeo &AV: ActiveVolumes) {
    fiducialVolumes.emplace_back(AV.MinX() + fidInsetMinX, AV.MaxX() - fidInsetMaxX,
                                 AV.MinY() + fidInsetMinY, AV.MaxY() - fidInsetMaxY,
                                 AV.MinZ() + fidInsetMinZ, AV.MaxZ() - fidInsetMaxZ);
  }
  return fiducialVolumes;
}

std::vector<geo::BoxBoundedGeo> TrajectoryMCSFitterUBoone::setExcludeVolumes() const {
  std::vector<geo::BoxBoundedGeo> excludeVolumes;
  if (excludeVolumes_.size()%6 != 0) {
    cout << "Error: excludeVolumes vector must have length multiple of 6, not excluding any regions" << endl;
    excludeVolumes.emplace_back(-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0);
    return excludeVolumes;
  }
  for (unsigned int i=0; i<excludeVolumes_.size()/6; i++) {
    excludeVolumes.emplace_back(excludeVolumes_[6*i+0], excludeVolumes_[6*i+1], excludeVolumes_[6*i+2], excludeVolumes_[6*i+3], excludeVolumes_[6*i+4], excludeVolumes_[6*i+5]);
  }
  return excludeVolumes;
}

bool TrajectoryMCSFitterUBoone::isInVolume(const std::vector<geo::BoxBoundedGeo> &volumes, const geo::Point_t &point) const {
  for (const geo::BoxBoundedGeo &volume: volumes) {
    if (point.X()>=volume.MinX() && point.X()<=volume.MaxX() &&
        point.Y()>=volume.MinY() && point.Y()<=volume.MaxY() &&
        point.Z()>=volume.MinZ() && point.Z()<=volume.MaxZ()) {
      return true;
    }
  }
  return false;
}
