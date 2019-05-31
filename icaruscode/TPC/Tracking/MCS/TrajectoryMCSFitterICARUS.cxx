#include "TrajectoryMCSFitterICARUS.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TH1.h"
#include "TFile.h"
#include "lardata/RecoBaseProxy/Track.h" //needed only if you do use the proxies


using namespace std;
using namespace trkf;
using namespace recob::tracking;

recob::MCSFitResult TrajectoryMCSFitterICARUS::fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst) const {

   //std::cout << " traj nhits " << traj.NHits() << std::endl;
   std::cout << " traj lenght " << traj.Length() << std::endl;
   GetOptimalSegLen(1000,traj.NPoints(),2,traj.Length());

std::cout << " D3p " << d3p << std::endl;
  //
  // Break the trajectory in segments of length approximately equal to segLen_
  //
  vector<size_t> breakpoints;
  vector<float> segradlengths;
  vector<float> cumseglens;
  breakTrajInSegments(traj, breakpoints, segradlengths, cumseglens);

std::cout << " n segments " << segradlengths.size() << std::endl;
  //
  // Fit segment directions, and get 3D angles between them
  //resi
  if (segradlengths.size()<2) return recob::MCSFitResult();
  vector<float> dtheta;
  Vector_t pcdir0;
  Vector_t pcdir1;
  for (unsigned int p = 0; p<segradlengths.size(); p++) {
    linearRegression(traj, breakpoints[p], breakpoints[p+1], pcdir1);
    if (p>0) {
      if (segradlengths[p]<-100. || segradlengths[p-1]<-100.) {
	dtheta.push_back(-999.);
      } else { 
	const double cosval = pcdir0.X()*pcdir1.X()+pcdir0.Y()*pcdir1.Y()+pcdir0.Z()*pcdir1.Z();
	//assert(std::abs(cosval)<=1);
	//units are mrad
	double dt = 1000.*acos(cosval);//should we try to use expansion for small angles?
	dtheta.push_back(dt);
std::cout << " linearfit angle " << dt << std::endl;    
  }
    }
    pcdir0 = pcdir1;
  }
 
 vector<Vector_t> barycenters;
 Vector_t bary;
vector<float> dthetaPoly;
for (unsigned int p = 0; p<segradlengths.size(); p++) {
    findSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary);
    barycenters.push_back(bary);
}
for (unsigned int p = 2; p<segradlengths.size(); p++) {
      if (segradlengths[p]<-100. || segradlengths[p-1]<-100. || segradlengths[p-2]<-100.) {
	dtheta.push_back(-999.);
      } else {
        Vector_t dbcp=barycenters[p]-barycenters[p-1];
        float norm=sqrt(dbcp.X()*dbcp.X()+dbcp.Y()*dbcp.Y()+dbcp.Z()*dbcp.Z());
        dbcp/=norm;
        std::cout << " dbcp " << dbcp << std::endl;
        Vector_t dbcm=barycenters[p-1]-barycenters[p-2];
         norm=sqrt(dbcm.X()*dbcm.X()+dbcm.Y()*dbcm.Y()+dbcm.Z()*dbcm.Z());
dbcm/=norm;
        std::cout << " dbcm " << dbcm << std::endl;
	const double cosval = dbcp.X()*dbcm.X()+dbcp.Y()*dbcm.Y()+dbcp.Z()*dbcm.Z();
	//assert(std::abs(cosval)<=1);
	//units are mrad
       std::cout << " cosval " << cosval << std::endl;	
       double dt = 1000.*acos(cosval);//should we try to use expansion for small angles?
	dthetaPoly.push_back(dt);
std::cout << " polygonal angle " << dt << std::endl;    

      }
  }
  //
  // Perform likelihood scan in forward and backward directions
  //

 // std::ofstream outDtheta("dtheta2gev.out",ios::app);
  //outDtheta << dtheta[0] << std::endl;
  vector<float> cumLenFwd;
  vector<float> cumLenBwd;
  for (unsigned int i = 0; i<cumseglens.size()-2; i++) {
    std::cout << " seglen " << cumseglens[i] << std::endl;
    cumLenFwd.push_back(cumseglens[i]);
    cumLenBwd.push_back(cumseglens.back()-cumseglens[i+2]);
  }
  const ScanResult fwdResult = doLikelihoodScan(dtheta, segradlengths, cumLenFwd, true,  momDepConst, pid);
  const ScanResult bwdResult = doLikelihoodScan(dtheta, segradlengths, cumLenBwd, false, momDepConst, pid);
  std::cout << " fwdResult " << fwdResult.p << " bwdResult " << bwdResult.p << std::endl;  
  const ScanResult fwdResultPoly = doLikelihoodScan(dthetaPoly, segradlengths, cumLenFwd, true,  momDepConst, pid);
  const ScanResult bwdResultPoly = doLikelihoodScan(dthetaPoly, segradlengths, cumLenBwd, false, momDepConst, pid);
  std::cout << " fwdResultPoly " << fwdResultPoly.p << " bwdResultPoly " << bwdResultPoly.p << std::endl;  
//
for(unsigned int j=0;j<segradlengths.size();j++)
 std::cout << " dtheta " << dtheta.at(j) << std::endl;
  return recob::MCSFitResult(pid,
			     fwdResult.p,fwdResult.pUnc,fwdResult.logL,
			     bwdResult.p,bwdResult.pUnc,bwdResult.logL,
			     segradlengths,dtheta);
}

void TrajectoryMCSFitterICARUS::breakTrajInSegments(const recob::TrackTrajectory& traj, vector<size_t>& breakpoints, vector<float>& segradlengths, vector<float>& cumseglens) const {
  //
  const double trajlen = traj.Length();
  const int nseg = std::max(minNSegs_,int(trajlen/segLen_));
  const double thisSegLen = trajlen/double(nseg);
  // std::cout << "track with length=" << trajlen << " broken in nseg=" << nseg << " of length=" << thisSegLen << " where segLen_=" << segLen_ << std::endl;
  //
  constexpr double lar_radl_inv = 1./14.0;
  cumseglens.push_back(0.);//first segment has zero cumulative length from previous segments
  double thislen = 0.;
  auto nextValid=traj.FirstValidPoint();
  breakpoints.push_back(nextValid);
  auto pos0 = traj.LocationAtPoint(nextValid);
  nextValid = traj.NextValidPoint(nextValid+1);
  int npoints = 0;
  while (nextValid!=recob::TrackTrajectory::InvalidIndex) {
    auto pos1 = traj.LocationAtPoint(nextValid);
    thislen += ( (pos1-pos0).R() );
    pos0=pos1;
    npoints++;
    if (thislen>=thisSegLen) {
      breakpoints.push_back(nextValid);
      if (npoints>=minHitsPerSegment_) segradlengths.push_back(thislen*lar_radl_inv);
      else segradlengths.push_back(-999.);
      cumseglens.push_back(cumseglens.back()+thislen);
      thislen = 0.;
      npoints = 0;
    }
    nextValid = traj.NextValidPoint(nextValid+1);
  }
  //then add last segment
  if (thislen>0.) {
    breakpoints.push_back(traj.LastValidPoint()+1);
    segradlengths.push_back(thislen*lar_radl_inv);
    cumseglens.push_back(cumseglens.back()+thislen);
  }
  return;
}

const TrajectoryMCSFitterICARUS::ScanResult TrajectoryMCSFitterICARUS::doLikelihoodScan(std::vector<float>& dtheta, std::vector<float>& seg_nradlengths, std::vector<float>& cumLen, bool fwdFit, bool momDepConst, int pid) const {
  int    best_idx  = -1;
  double best_logL = std::numeric_limits<double>::max();
  double best_p    = -1.0;
  std::vector<float> vlogL;
 for (double p_test = pMin_; p_test <= pMax_; p_test+=pStep_) {
    double logL = mcsLikelihood(p_test, angResol_, dtheta, seg_nradlengths, cumLen, fwdFit, momDepConst, pid);
    if (logL < best_logL) {
      best_p    = p_test;
      best_logL = logL;
      best_idx  = vlogL.size();
    }
// std::cout << " ptest " << p_test << " likeli " << logL << " bestp "<< best_p << std::endl;
//compute likelihood for MC momentum
  /*  double logL = mcsLikelihood(2., angResol_, dtheta, seg_nradlengths, cumLen, fwdFit, momDepConst, pid);
    if (logL < best_logL) {
      best_p    = 2.;
      best_logL = logL;
      best_idx  = vlogL.size();
    }*/
    vlogL.push_back(logL);
  }
  //
  //uncertainty from left side scan
  double lunc = -1.0;
  if (best_idx>0) {
    for (int j=best_idx-1;j>=0;j--) {
      double dLL = vlogL[j]-vlogL[best_idx];
      if ( dLL<0.5 ) {
	lunc = (best_idx-j)*pStep_;
      } else break;
    }
  }
  //uncertainty from right side scan
  double runc = -1.0;
  if (best_idx<int(vlogL.size()-1)) {  
    for (unsigned int j=best_idx+1;j<vlogL.size();j++) {
      double dLL = vlogL[j]-vlogL[best_idx];
      if ( dLL<0.5 ) {
	runc = (j-best_idx)*pStep_;
      } else break;
    }
  }
  return ScanResult(best_p, std::max(lunc,runc), best_logL);
}
void TrajectoryMCSFitterICARUS::findSegmentBarycenter(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& bary) const {
  int npoints = 0;
  geo::vect::MiddlePointAccumulator middlePointCalc;
  size_t nextValid = firstPoint;
  while (nextValid<lastPoint) {
    middlePointCalc.add(traj.LocationAtPoint(nextValid));
    nextValid = traj.NextValidPoint(nextValid+1);
    npoints++;
  }
  const auto avgpos = middlePointCalc.middlePoint();
  bary=avgpos;
  std::cout << " avgpos " << avgpos << std::endl;
}

void TrajectoryMCSFitterICARUS::linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir) const {
  //
  int npoints = 0;
  geo::vect::MiddlePointAccumulator middlePointCalc;
  size_t nextValid = firstPoint;
  while (nextValid<lastPoint) {
    middlePointCalc.add(traj.LocationAtPoint(nextValid));
    nextValid = traj.NextValidPoint(nextValid+1);
    npoints++;
  }
  const auto avgpos = middlePointCalc.middlePoint();
  const double norm = 1./double(npoints);
  //
  //assert(npoints>0);
  //
  TMatrixDSym m(3);
  nextValid = firstPoint;
  while (nextValid<lastPoint) {
    const auto p = traj.LocationAtPoint(nextValid);
    const double xxw0 = p.X()-avgpos.X();
    const double yyw0 = p.Y()-avgpos.Y();
    const double zzw0 = p.Z()-avgpos.Z();
    m(0, 0) += xxw0*xxw0*norm;
    m(0, 1) += xxw0*yyw0*norm;
    m(0, 2) += xxw0*zzw0*norm;
    m(1, 0) += yyw0*xxw0*norm;
    m(1, 1) += yyw0*yyw0*norm;
    m(1, 2) += yyw0*zzw0*norm;
    m(2, 0) += zzw0*xxw0*norm;
    m(2, 1) += zzw0*yyw0*norm;
    m(2, 2) += zzw0*zzw0*norm;
    nextValid = traj.NextValidPoint(nextValid+1);
  }
  //
  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  for (int i=1; i<3; ++i) {
    if (eigenval(i)>maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  } 
  //
  pcdir = Vector_t(eigenvec(0, maxevalidx),eigenvec(1, maxevalidx),eigenvec(2, maxevalidx));
  if (traj.DirectionAtPoint(firstPoint).Dot(pcdir)<0.) pcdir*=-1.;
  //
}

double TrajectoryMCSFitterICARUS::mcsLikelihood(double p, double theta0x, std::vector<float>& dthetaij, std::vector<float>& seg_nradl, std::vector<float>& cumLen, bool fwd, bool momDepConst, int pid) const {
  //
  const int beg  = (fwd ? 0 : (dthetaij.size()-1));
  const int end  = (fwd ? dthetaij.size() : -1);
  const int incr = (fwd ? +1 : -1);
  //
  bool print;
  if(p>1.29&&p<1.32)  print = true;//(p>1.999 && p<2.001);
  else print=false;
  //
  const double m = mass(pid);
  const double m2 = m*m;
  const double Etot = sqrt(p*p + m2);//Initial energy
  double Eij2 = 0.;
  //
  double const fixedterm = 0.5 * std::log( 2.0 * M_PI );
  double result = 0;
  for (int i = beg; i != end; i+=incr ) {
    if (dthetaij[i]<0) {
      //cout << "skip segment with too few points" << endl;
      continue;
    }
    //
    if (eLossMode_==1) {
      // ELoss mode: MIP (constant)
      constexpr double kcal = 0.002105;
      const double Eij = Etot - kcal*cumLen[i];//energy at this segment
      Eij2 = Eij*Eij;
    } else {
      // Non constant energy loss distribution
      const double Eij = GetE(Etot,cumLen[i],m);
      Eij2 = Eij*Eij;
    }
    //
    if ( Eij2 <= m2 ) {
      result = std::numeric_limits<double>::max();
      
      break;
    }
    const double pij = sqrt(Eij2 - m2);//momentum at this segment
    const double beta = sqrt( 1. - ((m2)/(pij*pij + m2)) );
    constexpr double tuned_HL_term1 = 11.0038; // https://arxiv.org/abs/1703.06187
    constexpr double HL_term2 = 0.038;
    const double tH0 = ( (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) / (pij*beta) ) * ( 1.0 + HL_term2 * std::log( seg_nradl[i] ) ) * sqrt( seg_nradl[i] );
    const double rms = sqrt( 2.0*( tH0 * tH0 + theta0x * theta0x ) );
    if (rms==0.0) {
      std::cout << " Error : RMS cannot be zero ! " << std::endl;
      return std::numeric_limits<double>::max();
    } 
    const double arg = dthetaij[i]/rms;
    result += ( std::log( rms ) + 0.5 * arg * arg + fixedterm);
    if (print && fwd==true) cout << "TrajectoryMCSFitterICARUS pij=" << pij << " dthetaij[i]=" << dthetaij[i] << " tH0=" << tH0 << " rms=" << rms << " prob=" << ( std::log( rms ) + 0.5 * arg * arg + fixedterm) << " const=" << (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) << " beta=" << beta << " red_length=" << seg_nradl[i] <<  " result " << result << endl;
  }
  //std::cout << " momentum " << p <<" likelihood " << result << std::endl; 
  return result;
}

double TrajectoryMCSFitterICARUS::energyLossLandau(const double mass2,const double e2, const double x) const {
  //
  // eq. (33.11) in http://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (except density correction is ignored)
  //
  if (x<=0.) return 0.;
  constexpr double Iinv2 = 1./(188.E-6*188.E-6);
  constexpr double matConst = 1.4*18./40.;//density*Z/A
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  constexpr double j = 0.200;
  //
  const double beta2 = (e2-mass2)/e2;
  const double gamma2 = 1./(1.0 - beta2);
  const double epsilon = 0.5*kappa*x*matConst/beta2;
  //
  return 0.001*epsilon*( log(2.*me*beta2*gamma2*epsilon*Iinv2) + j - beta2 );
}
//
double TrajectoryMCSFitterICARUS::energyLossBetheBloch(const double mass,const double e2) const {
  // stolen, mostly, from GFMaterialEffects.
  constexpr double Iinv = 1./188.E-6;
  constexpr double matConst = 1.4*18./40.;//density*Z/A
  constexpr double me = 0.511;
  constexpr double kappa = 0.307075;
  //
  const double beta2 = (e2-mass*mass)/e2;
  const double gamma2 = 1./(1.0 - beta2);
  const double massRatio = me/mass;
  const double argument = (2.*me*gamma2*beta2*Iinv) * std::sqrt(1+2*std::sqrt(gamma2)*massRatio + massRatio*massRatio);
  //
  double dedx = kappa*matConst/beta2;
  //
  if (mass==0.0) return(0.0);
  if (argument <= exp(beta2)) {
      dedx = 0.;
  } else{
    dedx *= (log(argument)-beta2)*1.E-3; // Bethe-Bloch, converted to GeV/cm
    if (dedx<0.) dedx = 0.;
  }
  return dedx;
}
//
double TrajectoryMCSFitterICARUS::GetE(const double initial_E, const double length_travelled, const double m) const {
  //
  const double step_size = length_travelled / nElossSteps_;
  //
  double current_E = initial_E;
  const double m2 = m*m;
  //
  for (auto i = 0; i < nElossSteps_; ++i) {
    if (eLossMode_==2) {
      double dedx = energyLossBetheBloch(m,current_E);
      current_E -= (dedx * step_size);
    } else {
      // MPV of Landau energy loss distribution
      current_E -= energyLossLandau(m2,current_E*current_E,step_size);
    }
    if ( current_E <= m ) {
      // std::cout<<"WARNING: current_E less than mu mass. it is "<<current_E<<std::endl;
      return 0.;
    }
  }
  return current_E;
}
double TrajectoryMCSFitterICARUS::GetOptimalSegLen(const double guess_p, const int n_points, const int plane, const double length_travelled) const {
  //
// check units of measurment! (energy, length...)
double initial_p=guess_p*-0.211*length_travelled;
double MIN_P=1000;
if(initial_p<MIN_P) initial_p=MIN_P;

double sigma=d3p; //to be replaced by computed D3P when ready!

//double cosa=1;
double sinb=1; //both to be replaced by angles when ready!
double DsTot=length_travelled; //to be checked!
double p0=13.6; //MeV! constant in MCS formula
double LogTerm=0.038; //log term in MCS formula
double X0=140.; //interaction length in mm

double xi=(sigma/DsTot)*(initial_p/2.)/p0*sqrt(X0/length_travelled)/sqrt(double(n_points))*sinb;
//double K=6*xi*xi;
double m=1/sqrt(6*xi);

//double sigmaMS=p0/(initial_p/2.)*sqrt(length_travelled/X0/m/cosa)/sqrt(2.)/cosa;
//double sigmaMeas=pow(m,1.5)*sqrt(double(6.))*sigma*sinb/cosa/DsTot/sqrt(double(n_points));

double alfa=1+LogTerm*log(DsTot/sinb/double(m)/X0);

double m_corr=(m==1)?m:sqrt(m/(m-1));

double xiCorr=sigma/DsTot*((initial_p/2.)/p0)*sqrt(X0/length_travelled)*sqrt((double)n_points)*sinb/alfa;

double m2=m_corr/sqrt(6.*xiCorr);
if(m==1) m2=10;
m2+=0.5;

return m2;

}
void TrajectoryMCSFitterICARUS::ComputeD3P()  {
{   
    double res;
    vector<double> h0;
    vector<double> alf;
    
    TH1D* hd3pv=new TH1D("hd3pv","hd3pv",100,-5.,5.);
    
    for (unsigned int j=0;j<hits2d.size()-2;j+=3) {
            
            double a;
            res=computeResidual(j,a);
            cout << " point "  << j << " residual " <<res << endl;
            //! for each triplet of consecutive hits, save absolute value of residuale
            if(abs(res)<5) { //mm
                h0.push_back(res);
                alf.push_back(a);
                hd3pv->Fill(res);
               
            }}
    
        bool writeHisto=true;
        if(writeHisto) {
       TFile *f = new TFile("d3phisto.root","UPDATE");
    hd3pv->Write();
        f->Close();
        f->Delete();
        }

    if(!h0.size())
       d3p=0.4;
    else
       d3p=hd3pv->GetRMS();
    
    //d3pvector=h0;
    //return d3p;
}
}
double TrajectoryMCSFitterICARUS::computeResidual(int i, double& alfa) const
{
    /*const auto p0 = traj.LocationAtPoint(i);
    const auto p1 = traj.LocationAtPoint(i+1);
    const auto p2 = traj.LocationAtPoint(i+2);

    auto x0=p0.X();  auto x1=p1.X();  auto x2=p2.X();
    auto y0=p0.Y();  auto y1=p1.Y();  auto y2=p2.Y();
   // auto z0=p0.Z();  auto z1=p1.Z();  auto z2=p2.Z();*/

recob::Hit h0=hits2d.at(i);
recob::Hit h1=hits2d.at(i+1);
recob::Hit h2=hits2d.at(i+2);

std::cout << " PeakTime " << h0.PeakTime() << std::endl;

float x0=h0.WireID().Wire*3; auto y0=h0.PeakTime()*0.622;
float x1=h1.WireID().Wire*3; auto y1=h1.PeakTime()*0.622;
float x2=h2.WireID().Wire*3; auto y2=h2.PeakTime()*0.622;

//std::cout << " x0 " << x0 << " x1 " << x1 << " x2 " << x2 << std::endl;
//std::cout << " y0 " << y0 << " y1 " << y1 << " y2 " << y2 << std::endl;
    double ym=y0+(x1-x0)/(x2-x0)*(y2-y0);
    if(abs(x2-x0)<0.001)
        return -999;
    if(abs(x2-x1)<0.001)
        return -999;
    if(abs(x1-x0)<0.001)
        return -999;
    double K=(x1-x0)/(x2-x0);
    
    alfa=1/sqrt(1+K*K+(1-K)*(1-K));
    double dy=y1-ym;
    double res=dy*(alfa);
   
    //cout << " y1 " << y1 << " ym " << ym << " alfa " << alfa << endl;
   // cout << " deltafit residual " << res << endl;
	   return res;
}
