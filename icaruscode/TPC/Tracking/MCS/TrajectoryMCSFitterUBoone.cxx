#include "TrajectoryMCSFitterUBoone.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

using namespace std;
using namespace trkf;
using namespace recob::tracking;
using PointFlags_t=recob::TrajectoryPointFlags;
using Flags_t=std::vector<PointFlags_t>;

recob::MCSFitResult TrajectoryMCSFitterUBoone::fitMcs( const recob::TrackTrajectory& traj, int pid) 
{
  std::cout << " fitting uboone " << std::endl;
  //cut the final portion of the track
  //float cutFinLen=40;
  for (size_t jp=0;jp<traj.NPoints();jp++) {
     // std::cout << " residual length " << traj.Length(jp) << std::endl;
 //     if(traj.Length(jp)<cutFinLen) traj.TrajectoryPoint(jp).SetIsValid(false);
}

  //
  // Break the trajectory in segments of length approximately equal to segLen_
  //
  vector<size_t> breakpoints;
  vector<float> segradlengths;
  vector<float> cumseglens;
breakpoints.clear();
segradlengths.clear();
cumseglens.clear();
//int cutMode=2;
breakpoints.clear();
segradlengths.clear();
cumseglens.clear();
 breakTrajInSegments(traj, breakpoints, segradlengths, cumseglens, cutMode_, cutLength_);
  //
  // Fit segment directions, and get 3D angles between them
  //
  if (segradlengths.size() < 2) return recob::MCSFitResult();
  vector<float> dtheta;
  Vector_t pcdir0;
  Vector_t pcdir1;
  for (unsigned int p = 0; p < segradlengths.size(); p++) {
    linearRegression(traj, breakpoints[p], breakpoints[p + 1], pcdir1);
    if (p > 0) {
      if (segradlengths[p] < -100. || segradlengths[p - 1] < -100.) { dtheta.push_back(-999.); }
      else {
        const double cosval =
          pcdir0.X() * pcdir1.X() + pcdir0.Y() * pcdir1.Y() + pcdir0.Z() * pcdir1.Z();
        //assert(std::abs(cosval)<=1);
        //units are mrad
        double dt = 1000. * acos(cosval); //should we try to use expansion for small angles?
        dtheta.push_back(dt);
      }
    }
    pcdir0 = pcdir1;
  }
  //
  // Perform likelihood scan in forward and backward directions
  //
  vector<float> cumLenFwd;
  vector<float> cumLenBwd;
  for (unsigned int i = 0; i < cumseglens.size() - 2; i++) {
    cumLenFwd.push_back(cumseglens[i]);
    cumLenBwd.push_back(cumseglens.back() - cumseglens[i + 2]);
  }
  double detAngResol = DetectorAngularResolution(std::abs(traj.StartDirection().Z()));
  std::cout << " before likelihood scan " << std::endl;
  const ScanResult fwdResult =
    doLikelihoodScan(dtheta, segradlengths, cumLenFwd, true, pid, detAngResol);
  //const ScanResult bwdResult =
  //i  doLikelihoodScan(dtheta, segradlengths, cumLenBwd, false, pid, detAngResol);
  //
  return recob::MCSFitResult(pid,
                             fwdResult.p,
                             fwdResult.pUnc,
                             fwdResult.logL,
                             0.,
                             0.,
                             0.,
                             segradlengths,
                             dtheta);
}

void TrajectoryMCSFitterUBoone::breakTrajInSegments(const recob::TrackTrajectory& traj,
                                              vector<size_t>& breakpoints,
                                              vector<float>& segradlengths,
						    vector<float>& cumseglens, int cutMode, float cutLength) const
{
  //
  // float finCutLen=50.;
  //  float iniUsedLen=100.;
   double trajlen=traj.Length();
    if(cutMode==1) trajlen = traj.Length()-cutLength;
    if(cutMode==2) trajlen = cutLength;
 
  art::ServiceHandle<geo::Geometry const> geom;
  auto const* _SCE = (applySCEcorr_ ? lar::providerFrom<spacecharge::SpaceChargeService>() : NULL);
  //

  const double thisSegLen =
    (trajlen > (segLen_ * minNSegs_) ? segLen_ : trajlen / double(minNSegs_));
   std::cout << "track with length=" << trajlen << " broken in nseg=" << std::max(minNSegs_,int(trajlen/segLen_)) << " of length=" << thisSegLen << " where segLen_=" << segLen_ << std::endl;
  //
  constexpr double lar_radl_inv = 1. / 14.0;
  cumseglens.push_back(0.); //first segment has zero cumulative length from previous segments
  double thislen = 0.;
  double totlen = 0.;
  auto nextValid = traj.FirstValidPoint();
  breakpoints.push_back(nextValid);

  auto pos0 = traj.LocationAtPoint(nextValid);
  if (applySCEcorr_) {
    geo::TPCID tpcid = geom->FindTPCAtPosition(pos0);
    geo::Vector_t pos0_offset(0., 0., 0.);
    if (tpcid.isValid) { pos0_offset = _SCE->GetCalPosOffsets(pos0, tpcid.TPC); }
    pos0.SetX(pos0.X() - pos0_offset.X());
    pos0.SetY(pos0.Y() + pos0_offset.Y());
    pos0.SetZ(pos0.Z() + pos0_offset.Z());
  }
  auto dir0 = traj.DirectionAtPoint(nextValid);
  nextValid = traj.NextValidPoint(nextValid + 1);
  int npoints = 0;
  while (nextValid != recob::TrackTrajectory::InvalidIndex) {
    // std::cout << " in while cycle totlen " << totlen <<  std::endl;
    auto pos1 = traj.LocationAtPoint(nextValid);
    // thislen += ( (pos1-pos0).R() );
    //  pos0=pos1;
    npoints++;
 bool condition=true;
  auto rrange = traj.Length(nextValid);
  float lengthFromStart=sqrt((traj.LocationAtPoint(nextValid)-traj.LocationAtPoint(traj.FirstValidPoint())).Mag2());
  if(cutMode==2) std::cout << " lengthFromStart "<< lengthFromStart << " cutLength " << cutLength <<  std::endl;
  if(cutMode==1) condition=( rrange>cutLength);
  if(cutMode==2) condition=(lengthFromStart<cutLength);
  //    std::cout << " nextValid "<< nextValid << " thislen " << thislen << "thisSegLen " << thisSegLen << " condition " << condition << std::endl;
    if (npoints == 0) dir0 = traj.DirectionAtPoint(nextValid);

    if (applySCEcorr_) {
      geo::TPCID tpcid = geom->FindTPCAtPosition(pos1);
      geo::Vector_t pos1_offset(0., 0., 0.);
      if (tpcid.isValid) { pos1_offset = _SCE->GetCalPosOffsets(pos1, tpcid.TPC); }
      pos1.SetX(pos1.X() - pos1_offset.X());
      pos1.SetY(pos1.Y() + pos1_offset.Y());
      pos1.SetZ(pos1.Z() + pos1_offset.Z());
    }
    //increments along the initial direction of the segment
    auto step = (pos1 - pos0).R();
    thislen += dir0.Dot(pos1 - pos0);
    totlen += step;
    pos0 = pos1;
    // //fixme: testing alternative approaches here
    // //test1: increments following scatters
    // auto step = (pos1-pos0).R();
    // thislen += step;
    // totlen += step;
    // pos0=pos1;
    // //test2: end-start distance along the initial direction of the segment
    // thislen = dir0.Dot(pos1-pos0);
    // totlen = (pos1-pos0).R();
    //
    npoints++;
    if (thislen >= (thisSegLen - segLenTolerance_)) {
     if(condition) {
      breakpoints.push_back(nextValid);
      if (npoints >= minHitsPerSegment_)
        segradlengths.push_back(thislen * lar_radl_inv);
      else
        segradlengths.push_back(-999.);
   

      cumseglens.push_back(totlen);
      std::cout << " adding segment " <<  cumseglens[cumseglens.size()-1] << " thislen " << thislen  << " thisSegLen " << thisSegLen << std::endl;
      thislen = 0.;
      npoints = 0;
    }
else break;
}
    nextValid = traj.NextValidPoint(nextValid + 1);
  }
  //then add last segment
  if (thislen > 0.) {
    std::cout << " last larrad1inv " << lar_radl_inv << std::endl;
    breakpoints.push_back(traj.LastValidPoint() + 1);
    segradlengths.push_back(thislen * lar_radl_inv);
    cumseglens.push_back(cumseglens.back() + thislen);
  }
  

  return;
}

const TrajectoryMCSFitterUBoone::ScanResult TrajectoryMCSFitterUBoone::doLikelihoodScan(
  std::vector<float>& dtheta,
  std::vector<float>& seg_nradlengths,
  std::vector<float>& cumLen,
  bool fwdFit,
  int pid,
  float pmin,
  float pmax,
  float pstep,
  float detAngResol) const
{
  int best_idx = -1;
  float best_logL = std::numeric_limits<float>::max();
  float best_p = -1.0;
  std::vector<float> vlogL;
  for (float p_test = pmin; p_test <= pmax; p_test += pstep) {
    float logL = mcsLikelihood(p_test, detAngResol, dtheta, seg_nradlengths, cumLen, fwdFit, pid);
    std::cout << " ptest " << p_test << " logL " << logL << " best_p " << best_p<< std::endl;
    if (logL < best_logL) {
      best_p = p_test;
      best_logL = logL;
      best_idx = vlogL.size();
    }
    vlogL.push_back(logL);
  }
   //
  //uncertainty from left side scan
  float lunc = -1.;
  if (best_idx > 0) {
    for (int j = best_idx - 1; j >= 0; j--) {
      float dLL = vlogL[j] - vlogL[best_idx];
      if (dLL >= 0.5) {
        lunc = (best_idx - j) * pstep;
        break;
      }
    }
  }
  //uncertainty from right side scan
  float runc = -1.;
  if (best_idx < int(vlogL.size() - 1)) {
    for (unsigned int j = best_idx + 1; j < vlogL.size(); j++) {
      float dLL = vlogL[j] - vlogL[best_idx];
      if (dLL >= 0.5) {
        runc = (j - best_idx) * pstep;
        break;
      }
    }
  }
 std::cout << " best_p " << best_p << std::endl;
  return ScanResult(best_p, std::max(lunc, runc), best_logL);
}

const TrajectoryMCSFitterUBoone::ScanResult TrajectoryMCSFitterUBoone::doLikelihoodScan(
  std::vector<float>& dtheta,
  std::vector<float>& seg_nradlengths,
  std::vector<float>& cumLen,
  bool fwdFit,
  int pid,
  float detAngResol) const
{

  //do a first, coarse scan
  const ScanResult& coarseRes = doLikelihoodScan(
    dtheta, seg_nradlengths, cumLen, fwdFit, pid, pMin_, pMax_, pStepCoarse_, detAngResol);

  float pmax = std::min(coarseRes.p + fineScanWindow_, pMax_);
  float pmin = std::max(coarseRes.p - fineScanWindow_, pMin_);
  if (coarseRes.pUnc < (std::numeric_limits<float>::max() - 1.)) {
    pmax = std::min(coarseRes.p + 2 * coarseRes.pUnc, pMax_);
    pmin = std::max(coarseRes.p - 2 * coarseRes.pUnc, pMin_);
  }

  //do the fine grained scan in a smaller region
  const ScanResult& refineRes =
    doLikelihoodScan(dtheta, seg_nradlengths, cumLen, fwdFit, pid, pmin, pmax, pStep_, detAngResol);

  return refineRes;
}

void TrajectoryMCSFitterUBoone::linearRegression(const recob::TrackTrajectory& traj,
                                           const size_t firstPoint,
                                           const size_t lastPoint,
                                           Vector_t& pcdir) const
{
  //
  art::ServiceHandle<geo::Geometry const> geom;
  auto const* _SCE = (applySCEcorr_ ? lar::providerFrom<spacecharge::SpaceChargeService>() : NULL);
  //
  int npoints = 0;
  geo::vect::MiddlePointAccumulator middlePointCalc;
  size_t nextValid = firstPoint;
  //fixme explore a max number of points to use for linear regression
  //while (nextValid<std::min(firstPoint+10,lastPoint)) {
  while (nextValid < lastPoint) {
    auto tempP = traj.LocationAtPoint(nextValid);
    if (applySCEcorr_) {
      geo::TPCID tpcid = geom->FindTPCAtPosition(tempP);
      geo::Vector_t tempP_offset(0., 0., 0.);
      if (tpcid.isValid) { tempP_offset = _SCE->GetCalPosOffsets(tempP, tpcid.TPC); }
      tempP.SetX(tempP.X() - tempP_offset.X());
      tempP.SetY(tempP.Y() + tempP_offset.Y());
      tempP.SetZ(tempP.Z() + tempP_offset.Z());
    }
    middlePointCalc.add(tempP);
    //middlePointCalc.add(traj.LocationAtPoint(nextValid));
    nextValid = traj.NextValidPoint(nextValid + 1);
    npoints++;
  }
  const auto avgpos = middlePointCalc.middlePoint();
  const double norm = 1. / double(npoints);
  //
  //assert(npoints>0);
  //
  TMatrixDSym m(3);
  nextValid = firstPoint;
  while (nextValid < lastPoint) {
    auto p = traj.LocationAtPoint(nextValid);
    if (applySCEcorr_) {
      geo::TPCID tpcid = geom->FindTPCAtPosition(p);
      geo::Vector_t p_offset(0., 0., 0.);
      if (tpcid.isValid) { p_offset = _SCE->GetCalPosOffsets(p, tpcid.TPC); }
      p.SetX(p.X() - p_offset.X());
      p.SetY(p.Y() + p_offset.Y());
      p.SetZ(p.Z() + p_offset.Z());
    }
    const double xxw0 = p.X() - avgpos.X();
    const double yyw0 = p.Y() - avgpos.Y();
    const double zzw0 = p.Z() - avgpos.Z();
    m(0, 0) += xxw0 * xxw0 * norm;
    m(0, 1) += xxw0 * yyw0 * norm;
    m(0, 2) += xxw0 * zzw0 * norm;
    m(1, 0) += yyw0 * xxw0 * norm;
    m(1, 1) += yyw0 * yyw0 * norm;
    m(1, 2) += yyw0 * zzw0 * norm;
    m(2, 0) += zzw0 * xxw0 * norm;
    m(2, 1) += zzw0 * yyw0 * norm;
    m(2, 2) += zzw0 * zzw0 * norm;
    nextValid = traj.NextValidPoint(nextValid + 1);
  }
  //
  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  for (int i = 1; i < 3; ++i) {
    if (eigenval(i) > maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  }
  //
  pcdir = Vector_t(eigenvec(0, maxevalidx), eigenvec(1, maxevalidx), eigenvec(2, maxevalidx));
  if (traj.DirectionAtPoint(firstPoint).Dot(pcdir) < 0.) pcdir *= -1.;
  //
}

double TrajectoryMCSFitterUBoone::mcsLikelihood(double p,
                                          double theta0x,
                                          std::vector<float>& dthetaij,
                                          std::vector<float>& seg_nradl,
                                          std::vector<float>& cumLen,
                                          bool fwd,
                                          int pid) const
{
  //
  const int beg = (fwd ? 0 : (dthetaij.size() - 1));
  const int end = (fwd ? dthetaij.size() : -1);
  const int incr = (fwd ? +1 : -1);
  //
  // bool print = false;
  //
  const double m = mass(pid);
  const double m2 = m * m;
  const double Etot = sqrt(p * p + m2); //Initial energy
  double Eij2 = 0.;
  //
  double const fixedterm = 0.5 * std::log(2.0 * M_PI);
  double result = 0;
  for (int i = beg; i != end; i += incr) {
    if (dthetaij[i] < 0) {
      //cout << "skip segment with too few points" << endl;
      continue;
    }
    //
    const double Eij = GetE(Etot, cumLen[i], m);
    Eij2 = Eij * Eij;
    //
    if (Eij2 <= m2) {
      result = std::numeric_limits<double>::max();
      break;
    }
    const double pij = sqrt(Eij2 - m2); //momentum at this segment
    const double beta = sqrt(1. - ((m2) / (pij * pij + m2)));
    constexpr double HighlandSecondTerm = 0.038;
    const double tH0 = (HighlandFirstTerm(pij) / (pij * beta)) *
                       (1.0 + HighlandSecondTerm * std::log(seg_nradl[i])) * sqrt(seg_nradl[i]);
    const double rms = sqrt(2.0 * (tH0 * tH0 + theta0x * theta0x));
    if (rms == 0.0) {
      std::cout << " Error : RMS cannot be zero ! " << std::endl;
      return std::numeric_limits<double>::max();
    }
    const double arg = dthetaij[i] / rms;
    result += (std::log(rms) + 0.5 * arg * arg + fixedterm);
  }
  return result;
}

double TrajectoryMCSFitterUBoone::energyLossLandau(const double mass2,
                                             const double e2,
                                             const double x) const
{
  //
  // eq. (33.11) in http://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (except density correction is ignored)
  //
  if (x <= 0.) return 0.;
  constexpr double Iinv2 = 1. / (188.E-6 * 188.E-6);
  constexpr double matConst = 1.4 * 18. / 40.; //density*Z/A
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  constexpr double j = 0.200;
  //
  const double beta2 = (e2 - mass2) / e2;
  const double gamma2 = 1. / (1.0 - beta2);
  const double epsilon = 0.5 * kappa * x * matConst / beta2;
  //
  return 0.001 * epsilon * (log(2. * me * beta2 * gamma2 * epsilon * Iinv2) + j - beta2);
}
//
double TrajectoryMCSFitterUBoone::energyLossBetheBloch(const double mass, const double e2) const
{
  // stolen, mostly, from GFMaterialEffects.
  constexpr double Iinv = 1. / 188.E-6;
  constexpr double matConst = 1.4 * 18. / 40.; //density*Z/A
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  //
  const double beta2 = (e2 - mass * mass) / e2;
  const double gamma2 = 1. / (1.0 - beta2);
  const double massRatio = me / mass;
  const double argument = (2. * me * gamma2 * beta2 * Iinv) *
                          std::sqrt(1 + 2 * std::sqrt(gamma2) * massRatio + massRatio * massRatio);
  //
  double dedx = kappa * matConst / beta2;
  //
  if (mass == 0.0) return (0.0);
  if (argument <= exp(beta2)) { dedx = 0.; }
  else {
    dedx *= (log(argument) - beta2) * 1.E-3; // Bethe-Bloch, converted to GeV/cm
    if (dedx < 0.) dedx = 0.;
  }
  return dedx;
}
//
double TrajectoryMCSFitterUBoone::GetE(const double initial_E,
                                 const double length_travelled,
                                 const double m) const
{
  //
  if (eLossMode_ == 1) {
    // ELoss mode: MIP (constant)
    constexpr double kcal = 0.002105;
    return (initial_E - kcal * length_travelled); //energy at this segment
  }
  //
  // Non constant energy loss distribution
  const double step_size = length_travelled / nElossSteps_;
  //
  double current_E = initial_E;
  const double m2 = m * m;
  //
  for (auto i = 0; i < nElossSteps_; ++i) {
    if (eLossMode_ == 2) {
      double dedx = energyLossBetheBloch(m, current_E);
      current_E -= (dedx * step_size);
    }
    else {
      // MPV of Landau energy loss distribution
      current_E -= energyLossLandau(m2, current_E * current_E, step_size);
    }
    if (current_E <= m) {
      // std::cout<<"WARNING: current_E less than mu mass. it is "<<current_E<<std::endl;
      return 0.;
    }
  }
  return current_E;
}
