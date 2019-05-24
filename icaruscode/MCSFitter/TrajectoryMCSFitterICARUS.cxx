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
 std::cout << " ptest " << p_test << " likeli " << logL << " bestp "<< best_p << std::endl;
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
      //result = std::numeric_limits<double>::max();
      
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
/*
double TrajectoryMCSFitterICARUS::RMSAngleMeasTrunc(recob::TrackTrajectory& tr, double nsigma, double& gau, int& ncut, double& ttave,bool write)
  //compute RMS of R distribution (see MCS ICARUS paper) i.e. ratio between measured 2D scattering angle and expcted theta_MS
{

 stringstream ssr;//create a stringstream
   ssr << run;//add number to the stream
   std::string sr=ssr.str();//return a string with the contents of the stream
 stringstream sse;//create a stringstream
   sse << evt;//add number to the stream
   std::string se=sse.str();//return a string with the contents of the stream
std::string name="MSdebug_"+sr+"_"+se+".out";
const char *cname = name.c_str();
ofstream debu(cname,ios::app);
cout << " DEBUG FILENAME MANAGER " << run << " " << evt << " write " << write << endl;
//nsigma=1000;
  int jerr=-1;
  wtcov.clear();
  wtnew.clear();
   chi2norm.clear();
verb=1; 

double thetams,thetaerr;
if(verb)
cout << " moduleangle d3pcu " << tr->D3PCU() << endl;
  vector<double> tt,ttall,tt0;
 double dstot,cos,beta,alfa;
 int n;

 if(tr->nPoints()<=2) {
 ncut=-1;
   return -1;
}
 int np=tr->nPointsUsed();
// np=tr->nPoints();
if(verb) {
cout << " tr npointsused " << np << endl;
 cout << " tr npointstot " << tr->nPoints() << endl;
}
EMatrix mat(2*np-3,2*np-3,0);
//EMatrix mat(np-2,np-2,0);
//EMatrix mat(np-1,np-1,0);
EMatrix materr(np-1,np-1,0);
EMatrix matpolyerr(np-2,np-2,0);
EMatrix matpolyms(np-2,np-2,0);
EMatrix matpoly(np-2,np-2,0);
EMatrix matcard(np-2,np-2,0);
EMatrix matmix(2*np-3,2*np-3,0);
 double ap0=0;
double apmedio=0;
double tpmedio=0;
 double cc=0;
int nbstop=0;


for(int jp=0;jp<tr->nPoints()-1;jp++) {
  if(tr->point(jp).isUsed()&&(!tr->point(jp).wIsOut())) {
    IcarusPoint& p=tr->point(jp);
      if((p.CP()>15)&&(p.DS()>0)) 
    
        nbstop++; 
}
}

//EMatrix mat(np-2,np-2,0);

for(int jp=1;jp<np;jp++) {
if(verb){
cout << " before isused " << jp << " nbstop " << nbstop << endl; 
cout << " point " << jp << " used " << tr->point(jp).isUsed() << endl;
}

cout << " point " << jp << " isused " << tr->point(jp).isUsed() << " isout " << tr->point(jp).isUsed() << endl;

  if(tr->point(jp).isUsed()&&(!tr->point(jp).wIsOut())) {
     IcarusPoint& p=tr->point(jp);
cout << " point " << jp << " cp " << p.CP() << " ds " << p.DS() << endl;
         if((p.CP()>0)&&(p.DS()>0)) {
             double c=p.CP();
   double d=p.DS();


  if(p.measurement().view()==2)
    cos=cosac;
  else
    cos=cosai;

double sigma0,sinb;
   int m;
IcarusPoint& pprev=tr->point(jp-1);
   double a=p.measVector()[2]-pprev.measVector()[2];
   double eseg=p.ESeg();


    beta=c/sqrt(c*c+MMU*MMU);

  //alfa=1;
  
   if(tr->D3PCU()>0) {
     sigma0=tr->D3PCU();
     sinb=tr->SinBC();
     dstot=tr->DsTotC();
dstot=4000*tr->CosAC();
     m=tr->nPoints();
     n=tr->NGrpC();
   }
   else {
      sigma0=tr->D3PIU()/sqrt(double(tr->NGrpI()));
      sinb=tr->SinBI();
     dstot=tr->DsTotI();
  m=tr->nPoints();
  
   n=tr->NGrpI();
   }
   if(!n)
     alfa=1;
   else
    alfa=1+0.038*log(dstot/m/140/cos);
    
double alfa1=1+0.038*log(dstot/m/140/cos);
double alfa2=1+0.038*log(dstot/n/140/sinb);
cout << " alfa1 " << alfa1 << " alfa2 " << alfa2 << endl;
   //sigma0=0.715;
//double washout=0.8585;
   double dxmedio=dstot/double(m);
    thetams=13.6/c/beta*sqrt(1./140./cos)*alfa/cos;

if(verb) {
cout << " eseg " << eseg << endl;
cout << " sigma0n " << sigma0/sqrt(n) << " n " << n << endl;
}
    thetaerr=sigma0*sqrt(12.);
  // thetaerr=0;
     double thetacorr=sqrt(thetams*thetams+thetaerr*thetaerr);
   //  ttall.push_back(acorr/thetams);
   // ttall.push_back(a/thetacorr);
   ttall.push_back(a);
   cout << jp << " Filling ttall " << ttall.size () << endl;
apmedio+=abs(a*a);
 ap0+=a;
tpmedio+=(thetams*thetams*d+thetaerr*thetaerr*6/d/d);
 cc+=(a*a)/(thetams*thetams*d+6*thetaerr*thetaerr/d/d);

	 }
 
 

	 tr->FillCovMatrixSegOnly(mat,jp,thetams*thetams,thetaerr*thetaerr,materr);
//	 cout << "addsegcov " << mat << " " << jp << endl;
	//tr->AddSegmentCovariance(mat,jp);

}
 }

cout << " end fit part ttall size " << ttall.size() << " np " << np << endl;

for(int jp=1;jp<np-1;jp++) {
if(verb){
cout << " before isused " << jp << " nbstop " << nbstop << endl; 
cout << " point " << jp << " used " << tr->point(jp).isUsed() << endl;
}

cout << " point " << jp << " isused " << tr->point(jp).isUsed() << " isout " << tr->point(jp).isUsed() << endl;

  if(tr->point(jp).isUsed()&&(!tr->point(jp).wIsOut())) {
     IcarusPoint& p=tr->point(jp);
cout << " point " << jp << " cp " << p.CP() << " ds " << p.DS() << endl;
         if((p.CP()>0)&&(p.DS()>0)) {
             double c=p.CP();
   double d=p.DS();


  if(p.measurement().view()==2)
    cos=cosac;
  else
    cos=cosai;

double sigma0,sinb;
   int m;
IcarusPoint& pprev=tr->point(jp-1);
   double a=p.AMeas();
   double eseg=p.ESeg();


    beta=c/sqrt(c*c+MMU*MMU);

  //alfa=1;
  
   if(tr->D3PCU()>0) {
     sigma0=tr->D3PCU();
     sinb=tr->SinBC();
     dstot=tr->DsTotC();
dstot=4000*tr->CosAC();
     m=tr->nPoints();
     n=tr->NGrpC();
   }
   else {
      sigma0=tr->D3PIU()/sqrt(double(tr->NGrpI()));
      sinb=tr->SinBI();
     dstot=tr->DsTotI();
  m=tr->nPoints();
  
   n=tr->NGrpI();
   }
   if(!n)
     alfa=1;
   else
    alfa=1+0.038*log(dstot/m/140/cos);
    

   //sigma0=0.715;
//double washout=0.8585;
   double dxmedio=dstot/double(m);
    //thetams=13.6/c/beta*sqrt(1./140./cos)*alfa/cos;
    thetams=13.6*sqrt(1./140./cos)*alfa/cos*0.7388;
if(verb) {
cout << " eseg " << eseg << endl;
cout << " sigma0n " << sigma0/sqrt(n) << " n " << n << endl;
}
    thetaerr=sigma0;
  // thetaerr=0;
     double thetacorr=sqrt(thetams*thetams+thetaerr*thetaerr);
   //  ttall.push_back(acorr/thetams);
   // ttall.push_back(a/thetacorr);
   ttall.push_back(a);
   cout << jp << " Filling ttall " << ttall.size () << endl;
apmedio+=abs(a*a);
 ap0+=a;
tpmedio+=(thetams*thetams*d+thetaerr*thetaerr*6/d/d);
 cc+=(a*a)/(thetams*thetams*d+6*thetaerr*thetaerr/d/d);

	 }
 
 
cout << " poly thetams " << thetams << " thetaerr " << thetaerr << endl;
	 tr->FillCovMatrix(matpoly,jp,thetams*thetams,thetaerr*thetaerr,matpolyms,matpolyerr,write);
tr->FillCovMixTerms(matmix,jp,np-2,thetams*thetams,thetaerr*thetaerr);
//	 cout << "addsegcov " << mat << " " << jp << endl;
	tr->AddSegmentCovariance(matpoly,jp,matcard);
//  cout << " after fillmatrix1  " << jp << " MAT " << mat << endl; 
for(int jm=0;jm<matpoly.num_row();jm++)
for(int jmm=0;jmm<matpoly.num_row();jmm++)
 mat[jm+np-1][jmm+np-1]=matpoly[jm][jmm];
}
 }



//mat+=matmix;

if(write)
debu << " direction matrix " << endl;
debu << materr << endl;
//cout << " MATERR " << materr << endl;
 EMatrix cov=mat;
 if(verb) cout << " cov " << cov << endl;
EVector vtall(ttall.size(),0);
for(int jv=0;jv<ttall.size();jv++)
  vtall[jv]=ttall[jv]; 
  cout << " both vtall " << vtall << endl;
 EMatrix invcov=cov.inverse(jerr);
 cout << " invcov " << cov.num_row() << endl;

 EVector vtcov=transpose(vtall)*cov.inverse(jerr)*vtall;
 if(verb) cout << " rotated vtall " << cov.inverse(jerr)*vtall << endl;
  EMatrix test1=(cov.inverse(jerr)*vtall);
for(int jv=0;jv<ttall.size();jv++)
  if(verb)cout << " rmsanglemeastrunc vtcov element "<< jv << " : " << transpose(vtall)[0][jv]*test1[jv][0]   << endl;

 vector<double> ttcov;
 for(int jv=0;jv<ttall.size();jv++)
   ttcov.push_back(transpose(vtall)[0][jv]*test1[jv][0]);

 double ttcovmed=0;
 for(int jv=0;jv<ttall.size();jv++)
   ttcovmed+=ttcov[jv];

 ttcovtot=ttcovmed;
 ttcovmed/=(ttall.size());

 
 
 // compute chi2 before truncation


int minhits=10;
//for(int kk=1;kk<nPoints();kk++) {
  vector<double> ttrunc;
vector<int> tails;
double ttrunctot=0;
vector<double> atrunc;
//nsigma=999;
for (int jt=0;jt<ttall.size();jt++) {
    if(abs(ttcov[jt])<ttcovmed*nsigma*nsigma) {
//   {
       ttrunctot+=ttcov[jt];
       atrunc.push_back(ttall[jt]);
       ttrunc.push_back(ttcov[jt]);
   if(verb)    std::cout << " keeping value " << ttcov[jt] << std::endl;
   }
    else {
if(verb) {cout << " ttcovmed " << ttcovmed << " nsigma " << nsigma << std::endl;
 std::cout << " discarding value " << ttcov[jt] << std::endl;}
      tails.push_back(jt);
    }
  }







 for (int jt=0;jt<np-2;jt++) {
int ist=0;
 for(int jta=0;jta<tails.size();jta++) {
  //if(tails[jta]+(np-1)==jt)
  // ist=1;
 // if(tails[jta]==-1&&jt==0)
 //  ist=1;
if(tails[jta]==jt-1)
   ist=1;
//if(tails[jta]==np-2&&jt==np-3)
//ist=1;
  }



cout << " poly jt " << jt << " ist " << ist << endl;
 }
 
  EVector vtrunc(atrunc.size(),0);
 for(int jv=0;jv<atrunc.size();jv++)
   vtrunc[jv]=atrunc[jv];

  EMatrix cov0=cov;
if(verb) {  cout << " Tails size " << tails.size() << " cov size " << cov.num_row() << endl;
  cout << " dirty covarience " << cov << endl;    
}

 if(tails.size())
    for(int jta=tails.size()-1;jta>=0;jta--) {
      cout << " jta " << jta << " tails " << tails[jta] << endl;
   
      int ncov=cov.num_row();
     
      EMatrix covmod=CleanCovariance(cov,tails[jta]); 
cout << " after cleaning cov3 " << cov.num_row() << " " << covmod.num_row() << endl;

      cov=covmod;     
   

}



 if(verb)cout << " Tails size " << tails.size() << " cleancov size " << cov.num_row() << endl;
 if(!cov.num_row())
   return -999;
 if(verb) cout << " cleaned covarience " << cov << endl;    

 invcov=cov.inverse(jerr);
 if(verb) { cout << " cleaned inverse covarience " << invcov.num_row() << endl;    
 cout << " vtrunc size " << vtrunc.num_row() << endl;
 cout << " vtrunc check " << vtrunc << endl;
}
 //EVector vtnew=transpose(vtrunc)*invcov*vtrunc;
 EMatrix test2=(cov.inverse(jerr)*vtrunc);
 vector<double> ttnew;
 for(int jv=0;jv<atrunc.size();jv++)
   ttnew.push_back(transpose(vtrunc)[0][jv]*test2[jv][0]);
for(int jv=0;jv<atrunc.size();jv++)
 if(verb) cout << " test2 element "<< jv << " : " << test2[jv][0]   << endl;

for(int jv=0;jv<atrunc.size();jv++)
 if(verb) cout << " vtrunc element "<< jv << " : " << transpose(vtrunc)[0][jv]*test2[jv][0]   << endl;

 double ttnewtot=0;
double ttnew1=0;
double ttnew2=0;
  for (int jt=0;jt<ttnew.size();jt++) 
         ttnewtot+=ttnew[jt];
 for (int jt=0;jt<ttnew.size()/2;jt++) 
         ttnew1+=ttnew[jt];
 for (int jt=ttnew.size()/2;jt<ttnew.size();jt++) 
         ttnew2+=ttnew[jt];
cout << " ttnewtot " << ttnewtot << endl;
cout << " ttnew1 " << ttnew1 << endl;
cout << " ttnew2 " << ttnew2 << endl;

int ndf=ttall.size()-tails.size();

double ttnewrms=ttnewtot/(double)ndf;
 
 ncut=ttall.size()-ttnew.size();
 if(verb) cout << " ncut fraction " << (double)ncut/ttall.size() << endl;

 double um=(double)(nsigma*nsigma/2.);

 double km=0.5;

 double corr=TMath::Gamma(km,um)-1./km/tgamma(km)*pow(um,km)*exp(-um);

 // ttnewrms/=corr;

 // for(int jp=0;jp<tr->nPoints()-2;jp++)
 //   chi2norm.push_back(Chi2Normal(jp,cov0,ttall));  

if(verb) cout << " ttcovmed " << ttcovmed << endl;
 for(int jp=0;jp<tr->nPoints()-2;jp++) {
   //  cout << " point " << jp << " chi2norm " << chi2norm[jp] << endl;
   // cout << " point " << jp << " ttcov " << ttcov[jp] << endl;
}

     
GaussCorrection* gauss=new GaussCorrection();
double gcorr=gauss->GaussFactor(ttcov.size(),1./sqrt(3.));
if(verb) {cout << " gauss correction " << gcorr << endl;
cout << " chi2 correction " << corr << endl;
}
//tnewrms/=gcorr;



 gau=gcorr;

 wtcov=ttcov;

 ecorr=corr;
if(verb)cout << " rms after gausscorr " << ttnewrms << endl;
//exit(44); 
     return ttnewrms;
 
 
   
}

double TrajectoryMCSFitterICARUS::mcsCorrelatedLikelihood(double p, double theta0x, std::vector<float>& dthetaij, std::vector<float>& seg_nradl, std::vector<float>& cumLen, bool fwd, bool momDepConst, int pid) const {
  //
  const int beg  = (fwd ? 0 : (dthetaij.size()-1));
  const int end  = (fwd ? dthetaij.size() : -1);
  const int incr = (fwd ? +1 : -1);
  //
  // bool print = false;//(p>1.999 && p<2.001);
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
    // if (print && fwd==true) cout << "TrajectoryMCSFitterICARUS pij=" << pij << " dthetaij[i]=" << dthetaij[i] << " tH0=" << tH0 << " rms=" << rms << " prob=" << ( std::log( rms ) + 0.5 * arg * arg + fixedterm) << " const=" << (momDepConst ? MomentumDependentConstant(pij) : tuned_HL_term1) << " beta=" << beta << " red_length=" << seg_nradl[i] << endl;
  }
  return result;
}
*/
