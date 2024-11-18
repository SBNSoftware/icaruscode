
#include "TrajectoryMCSFitterICARUS.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "lardata/RecoBaseProxy/Track.h" //needed only if you do use the proxies


using namespace std;
using namespace trkf;
using namespace recob::tracking;

recob::MCSFitResult TrajectoryMCSFitterICARUS::fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst) const {

  // std::cout << " traj nhits " << traj.NPoints() << std::endl;
   //std::cout << " traj lenght " << traj.Length() << std::endl;
float minlen=40;
if(traj.Length()<minlen) {
  vector<float> dum;
  return recob::MCSFitResult(pid,
			    0,0,0,
0,0,0,
			     dum,dum);
}


 GetOptimalSegLen(traj,1000,traj.NPoints(),2,traj.Length());

 // std::cout << " D3p " << d3p << std::endl;
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

  std::cout << " n segments " << segradlengths.size() << std::endl;
  //
  // Fit segment directions, and get 3D angles between them
  //resi
  if (segradlengths.size()<2) return recob::MCSFitResult();
  vector<float> dtheta;
  Vector_t pcdir0;
  Vector_t pcdir1;
    std::cout << " segradl size " << segradlengths.size() << " cumseg size " << cumseglens.size() <<  std::endl;

  for (unsigned int p = 0; p<segradlengths.size(); p++) {
    //linearRegression2D(traj, breakpoints[p], breakpoints[p+1], pcdir1); //2d
    linearRegression(traj, breakpoints[p], breakpoints[p+1], pcdir1); //3d
    if (p>0) {
      cout << " p " << p << " segradlength " << segradlengths[p] << endl;
      if (segradlengths[p]<-100. || segradlengths[p-1]<-100.
     || (cutMode_==2 && (p == 1 ||  p > 7))
      ) 
      {
         dtheta.push_back(-999.); 
         std::cout << " excluding dtheta " << std::endl;
 
      } 
      else { 
                 std::cout << " including dtheta " << std::endl;

	const double cosval = pcdir0.X()*pcdir1.X()+pcdir0.Y()*pcdir1.Y()+pcdir0.Z()*pcdir1.Z();
	//assert(std::abs(cosval)<=1);
	//units are mrad
	double dt = acos(cosval);//should we try to use expansion for small angles?
	dtheta.push_back(dt);
  //if(dt>3.) {
    cout << " pcdir0 "<< pcdir0 << endl;
    cout << " pcdir1 "<< pcdir1 << endl;

  //}
  //if(dt>3.) 
  std::cout << " cosval " << cosval << " linearfit angle " << dt <<  " dtheta size " << dtheta.size() << std::endl;
  }
    }
    pcdir0 = pcdir1;
  }
 
 vector<Vector_t> barycenters;
 Vector_t bary;
vector<float> dthetaPoly;
for (unsigned int p = 0; p<segradlengths.size(); p++) {
    cout << " before finding 2d barycenters for poly " << endl;
  //  find2DSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary); //2d
   findSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary); //3d
    barycenters.push_back(bary);
}
for (unsigned int p = 2; p<segradlengths.size(); p++) {
      if (segradlengths[p]<-100. || segradlengths[p-1]<-100. || segradlengths[p-2]<-100.
//      ||(cutMode_==2 && p > 7)
      ) {
	dthetaPoly.push_back(-999.);
      } else {
        Vector_t dbcp=barycenters[p]-barycenters[p-1];
        float norm=sqrt(dbcp.X()*dbcp.X()+dbcp.Y()*dbcp.Y()+dbcp.Z()*dbcp.Z());
        dbcp/=norm;
        cout << " barycenters " << barycenters[p-1] << " " << barycenters[p-2] << endl;
       std::cout << " dbcp " << dbcp << std::endl;
        Vector_t dbcm=barycenters[p-1]-barycenters[p-2];
         norm=sqrt(dbcm.X()*dbcm.X()+dbcm.Y()*dbcm.Y()+dbcm.Z()*dbcm.Z());
dbcm/=norm;
     std::cout << " dbcm " << dbcm << std::endl;
	const double cosval = dbcp.X()*dbcm.X()+dbcp.Y()*dbcm.Y()+dbcp.Z()*dbcm.Z();
	//assert(std::abs(cosval)<=1);
	//units are mrad
      std::cout << " cosval " << cosval << std::endl;
       double dt = acos(cosval);//should we try to use expansion for small angles?
	dthetaPoly.push_back(dt);
    std::cout << " polygonal angle " << dt*1000. << std::endl;

      }
  }
  //
  //
 std::cout << " before computing c2 " << std::endl;
 std::vector<float> ttall600; ttall600.clear();

 double c2= C2Function(traj,cumseglens,breakpoints,dtheta,dthetaPoly,ttall600,600.);
 std::cout << " c2 " << c2 << std::endl;
 //exit(11);
//return recob::MCSFitResult(pid,0.,0.,0.,0.,0.,0.,segradlengths,dtheta);
const ScanResult fwdResult = C2Fit(dtheta,dthetaPoly, segradlengths, cumseglens, breakpoints,true,  1., pid,0.05,traj);
dtheta.insert( dtheta.end(), dthetaPoly.begin(), dthetaPoly.end() );
  return recob::MCSFitResult(pid, fwdResult.p,fwdResult.pUnc,fwdResult.logL,fwdResult.p,fwdResult.pUnc,fwdResult.logL,ttall600,dtheta);
}
void TrajectoryMCSFitterICARUS::breakTrajInSegments(const recob::TrackTrajectory& traj, vector<size_t>& breakpoints, vector<float>& segradlengths, vector<float>& cumseglens, int cutMode, float cutLength) const {
  //
//cout << " breakpoints size " <<  breakpoints.size() << std::endl;
  // float finCutLen=50.;
  // float iniUsedLen=100.;
   double trajlen=traj.Length();
    if(cutMode==1) trajlen = traj.Length()-cutLength;
    if(cutMode==2) trajlen = cutLength;
  const int nseg = std::max(minNSegs_,int(trajlen/segLen_));
  const double thisSegLen = trajlen/double(nseg);
   //std::cout << "track with length=" << trajlen << " broken in nseg=" << nseg << " of length=" << thisSegLen << " where segLen_=" << segLen_ << std::endl;
  //
  constexpr double lar_radl_inv = 1./14.0;
  cumseglens.push_back(0.);//first segment has zero cumulative length from previous segments
  double thislen = 0.;
  auto nextValid=traj.FirstValidPoint();
  breakpoints.push_back(nextValid);
  auto pos0 = traj.LocationAtPoint(nextValid);



  nextValid = traj.NextValidPoint(nextValid+1);
  int npoints = 0;
 

  while (nextValid != recob::TrackTrajectory::InvalidIndex) {
std::cout << " in while cycle " << nextValid << std::endl;
    auto pos1 = traj.LocationAtPoint(nextValid);
    thislen += ( (pos1-pos0).R() );
    pos0=pos1;
    npoints++;
 bool condition=true;
  auto rrange = traj.Length(nextValid);
  float lengthFromStart=sqrt((traj.LocationAtPoint(nextValid)-traj.LocationAtPoint(traj.FirstValidPoint())).Mag2());
  // if(cutMode==2) std::cout << " lengthFromStart "<< lengthFromStart << " iniUsedLen " << cutLength <<  std::endl;
  if(cutMode==1) condition=( rrange>cutLength);
  if(cutMode==2) condition=(lengthFromStart<cutLength);
      std::cout << " nextValid "<< nextValid << " thislen " << thislen << "thisSegLen " << thisSegLen << " condition " << condition << std::endl;
    if (thislen>=thisSegLen) {
     if(condition) {
     std::cout << " adding breakpoint " <<  nextValid << std::endl;

      breakpoints.push_back(nextValid);
          std::cout << " break size " << breakpoints.size() << " last break " << breakpoints[breakpoints.size()-1] << std::endl;

      if (npoints>=minHitsPerSegment_) segradlengths.push_back(thislen*lar_radl_inv);
      else segradlengths.push_back(-999.);
      if(segradlengths[segradlengths.size()-1]==-999.) 
       cout << " adding weird segradlength npoints " << npoints << endl;
      cumseglens.push_back(cumseglens.back()+thislen);
     std::cout << " break size " << breakpoints.size() << " last cumseglen " << cumseglens[cumseglens.size()-1] << " last segradlen " << segradlengths[segradlengths.size()-1] <<std::endl;

      thislen = 0.;
      npoints = 0;
}
else {
  if (thislen>0.) {
    //then add last segment
    breakpoints.push_back(nextValid);
      std::cout << " break size " << breakpoints.size() << " very last break " << breakpoints[breakpoints.size()-1] << std::endl;

    segradlengths.push_back(thislen*lar_radl_inv);
    cumseglens.push_back(cumseglens.back()+thislen);
    std::cout << " break size " << breakpoints.size() << " very last cumseglen " << cumseglens[cumseglens.size()-1] << " last segradlen " << segradlengths[segradlengths.size()-1] <<std::endl;

  }
break;
}
    }
    nextValid = traj.NextValidPoint(nextValid+1);
  

  }
//exit(33);
  return;
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
  std::cout << " 3dbary avgpos " << avgpos << std::endl;
}

void TrajectoryMCSFitterICARUS::find2DSegmentBarycenter(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& bary) const {
  int npoints = 0;
  float wsum=0; float ssum=0;
  size_t nextValid = firstPoint;
  std::vector<recob::Hit> v;
  //unsigned int iniTPC=0;
  // Get track collection proxy and parallel mcs fit data (associated hits loaded by default)
  // Note: if tracks were produced from a TrackTrajectory collection you could access the original trajectories adding ',proxy::withOriginalTrajectory()' to the list of arguments
cout << " pdata size " << pdata.size() << endl;
  while (nextValid<lastPoint) {
   proxy::TrackPointData pd=pdata[nextValid];
   auto hit=std::get<1>(pd);
  if(hit->WireID().Plane==2) {
 unsigned int tpc=hit->WireID().TPC;
 unsigned int cryo=hit->WireID().Cryostat;
  unsigned int plane=hit->WireID().Plane;
TVectorD proj(3);
proj(0)=hit->WireID().Wire*3; proj(1)=(hit->PeakTime()*0.622); proj[2]=0;
TVectorD xyz=ReferenceFrame(plane,tpc,cryo)*proj;
double ori=DriftOrigin(plane,tpc,cryo);
cout << " reference frame xyz " << xyz(0) << " " << xyz(1) << " " << xyz(2) << endl;
wsum+=xyz(0);
double absx=xyz(1)+ori;
ssum+=absx;
   cout << "2dbary wireid " << hit->WireID() << " peaktime " << hit->PeakTime() << " xyz1 " << xyz[1] <<" absx " << absx << endl;
    npoints++;
  }
      nextValid = traj.NextValidPoint(nextValid+1);

  }
  const auto wmed = float(wsum/npoints);
  const auto smed = float(ssum/npoints);
 
 //Vector_t bary;
 bary.SetXYZ(wmed,smed,0);
std::cout << " wmed " << wmed << " smed " << smed << std::endl;
}

void TrajectoryMCSFitterICARUS::linearRegression2D(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir) const {
  //
  int npoints = 0;
  geo::vect::MiddlePointAccumulator middlePointCalc;
  size_t nextValid = firstPoint;
  while (nextValid<lastPoint) {
    middlePointCalc.add(traj.LocationAtPoint(nextValid));
    nextValid = traj.NextValidPoint(nextValid+1);
    npoints++;
  }

  Vector_t avgpos;
   cout << " before finding 2d barycenters for fit " << endl;
  //find2DSegmentBarycenter(traj,firstPoint,lastPoint,avgpos); //2d
  findSegmentBarycenter(traj,firstPoint,lastPoint,avgpos); //3d
  const double norm = 1./double(npoints);
  //
  //assert(npoints>0);
  //
  TMatrixDSym m(2);
  nextValid = firstPoint;
  while (nextValid<lastPoint) {
    const auto p = traj.LocationAtPoint(nextValid);
    const double xxw0 = p.X()-avgpos.X();
    const double yyw0 = p.Y()-avgpos.Y();
    m(0, 0) += xxw0*xxw0*norm;
    m(0, 1) += xxw0*yyw0*norm;
    m(1, 0) += yyw0*xxw0*norm;
    m(1, 1) += yyw0*yyw0*norm;
  
    nextValid = traj.NextValidPoint(nextValid+1);
  }
  //
  const TMatrixDSymEigen me(m);
  const auto& eigenval = me.GetEigenValues();
  const auto& eigenvec = me.GetEigenVectors();
  //
  int maxevalidx = 0;
  double maxeval = eigenval(0);
  for (int i=1; i<2; ++i) {
    if (eigenval(i)>maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  } 
  //
  pcdir = Vector_t(eigenvec(0, maxevalidx),eigenvec(1, maxevalidx),0.);
  if (traj.DirectionAtPoint(firstPoint).Dot(pcdir)<0.) pcdir*=-1.;
  //
  cout << " pcdir " << pcdir.X() << " " << pcdir.Y() << endl;
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
      double dedx = 1000.*energyLossBetheBloch(m,current_E);
      current_E -= (dedx * step_size);
    } else {
      // MPV of Landau energy loss distribution
      current_E -= 1000.*energyLossLandau(m2,current_E*current_E,step_size);
    }
    if ( current_E <= m ) {
      // std::cout<<"WARNING: current_E less than mu mass. it is "<<current_E<<std::endl;
      return 0.;
    }
  }
  return current_E;
}
double TrajectoryMCSFitterICARUS::GetOptimalSegLen(const recob::TrackTrajectory& tr,const double guess_p, const int n_points, const int plane, const double length_travelled) const {
  //
// check units of measurment! (energy, length...)
double initial_p=guess_p*-0.211*length_travelled;
double MIN_P=1000;
if(initial_p<MIN_P) initial_p=MIN_P;

double sigma=d3p; //to be replaced by computed D3P when ready!

//double cosa=1;
double sinb=1; //both to be replaced by angles when ready!
 double DsTot;
 if(cutMode()==0) DsTot=sqrt((tr.LocationAtPoint(tr.LastValidPoint())-tr.LocationAtPoint(0)).Mag2());
if(cutMode()==2) DsTot=cutLength();
 if(cutMode()==1) DsTot=sqrt((tr.LocationAtPoint(tr.LastValidPoint())-tr.LocationAtPoint(0)).Mag2())-cutLength();
 TVector3 start=tr.LocationAtPoint<TVector3>(0);
  TVector3 end=tr.LocationAtPoint<TVector3>(tr.LastValidPoint());
  TVector3 avdir=end-start;
  TVector3 projColl(0,0,0);
  projColl[0]=1; projColl[1]=sqrt(3.)/2.; projColl[2]=0.5;
  //  double cos=avdir*projColl;
  //double DsTot=length_travelled; //to be checked!
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
void TrajectoryMCSFitterICARUS::ComputeD3P(int plane)  {
{   
    double res;
    vector<double> h0;
    vector<double> alf;
    
    std::vector<recob::Hit> hits;
    if(plane==0) hits=hits2dI1;
    if(plane==1) hits=hits2dI2;
    if(plane==2) hits=hits2dC;

    TH1D* hd3pv=new TH1D("hd3pv","hd3pv",100,-5.,5.);
if(!hits.size()) return;    

    for (unsigned int j=0;j<hits.size()-2;j+=3) {
            
            double a;
            res=computeResidual(j,a,hits);
            //cout << " point "  << j << " residual " <<res << endl;
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
double TrajectoryMCSFitterICARUS::computeResidual(int i, double& alfa, std::vector<recob::Hit> hits2d) const
{


recob::Hit h0=hits2d.at(i);
recob::Hit h1=hits2d.at(i+1);
recob::Hit h2=hits2d.at(i+2);

//std::cout << " PeakTime " << h0.PeakTime() << std::endl;

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

const double TrajectoryMCSFitterICARUS::C2Function(const recob::TrackTrajectory& tr, std::vector<float> cumseglens, std::vector<long unsigned int> breakpoints, std::vector<float> dtheta,std::vector<float> dthetaPoly,std::vector<float>& ttall, double p0) const
  //compute RMS of R distribution (see MCS ICARUS paper) i.e. ratio between measured 2D scattering angle and expcted theta_MS
{
  std::vector<double> wtcov,wtnew,chi2norm;
//nsigma=1000;
  wtcov.clear();
  wtnew.clear();
   chi2norm.clear();
   //int verb=1; 
   std::cout << " after clears " << std::endl;
double thetams,thetaerr;
  vector<double> tt,tt0;
 double dstot,beta,alfa;
 int n;
 int firstseg=1;
 unsigned int lastseg=cumseglens.size()-1;
 for(unsigned int js=0;js<cumseglens.size()-1;js++) {
  if(dtheta[js]!=-999&&js)
   {
    firstseg=js;
    break;
   }
 }
  for(unsigned int js=cumseglens.size()-2;js>=0;js--) {
  if(dtheta[js]!=-999)
   {
    lastseg=js;
    break;
   }
 }
  
 unsigned int nseg=lastseg-firstseg+1;
 // int ncut;
 if(nseg<3) {
   // ncut=-1;
   return 0;
}
 std::cout << " before matrices " << nseg << std::endl;

TMatrixDSym mat(2*nseg-3);

TMatrixDSym materr(nseg-1);
TMatrixDSym matpolyerr(nseg-2);
TMatrixDSym matpolyms(nseg-2);
TMatrixDSym matpoly(nseg-2);
TMatrixDSym matcard(nseg-2);
TMatrixDSym matmix(2*nseg-3);
 double ap0=0;
double apmedio=0;
double tpmedio=0;
 double cc=0;
 int nbstop=0;
 double MMU=105.6;
 
 int np=tr.NPoints();
  TVector3 start=tr.LocationAtPoint<TVector3>(0);
  TVector3 end=tr.LocationAtPoint<TVector3>(tr.LastValidPoint());
  TVector3 avdir=end-start;
  TVector3 projColl(0,0,0);
  projColl[0]=1; projColl[1]=sqrt(3.)/2.; projColl[2]=0.5;
 //double  cos=avdir*projColl;
 double sinb=avdir*projColl;
 double cos=(collLength()/10./(tr.Length())); //2d
  cos=1; //3d
 sinb=1;
 cout << " sinb " << sinb << endl;
 //  double dstot;

  
 cout << " nseg " << nseg <<  " dtheta size " << dtheta.size() << " firstseg " << firstseg << endl;
for(unsigned int jp=firstseg;jp<lastseg;jp++) {
  // if(tr.point(jp).isUsed()&&(!tr.point(jp).wIsOut()))
  //{
  float cp=0;
  float ds=0;

      ds=sqrt((tr.LocationAtPoint(breakpoints[jp])-tr.LocationAtPoint(breakpoints[jp-1])).Mag2());
 
  const double Etot = sqrt(p0*p0 + MMU*MMU);//Initial energy
  double Eij2 = 0.;
  //double result = 0;
    if (eLossMode_==1) {
      // ELoss mode: MIP (constant)
      constexpr double kcal = 2.105;
      const double Eij = Etot - kcal*cumseglens[jp];//energy at this segment
      Eij2 = Eij*Eij;
    } else {
      // Non constant energy loss distribution
      const double Eij = GetE(Etot,cumseglens[jp],MMU);
      Eij2 = Eij*Eij;
    }
    //
    if ( Eij2 <= MMU*MMU ) {
      std::cout << " ptest " << p0 << " hitting Bragg peak, returning 0 " << std::endl;
      return 0;
    }
    const double pij = sqrt(Eij2 - MMU*MMU);//momentum at this segment
      cp=pij;
      if(cp>15.&&ds>0)    
        nbstop++; 


  double a=dtheta[jp];
 
  cout<< " filling a: jp " << jp << " jp-1 " << jp-1 << " a " << a << endl;
  double eseg=0;

  double MMU=105.6;
   beta=cp/sqrt(cp*cp+MMU*MMU);

  alfa=1;
  double sigma0;
  
   //  if(tr->D3PCU()>0) {
  //ComputeD3P();
  sigma0=d3p;
  cout << " eseg " << eseg << " sigma0 " <<sigma0 << endl;

 int np=tr.NPoints();
 if(cutMode()==0) dstot=sqrt((tr.LocationAtPoint(tr.LastValidPoint())-tr.LocationAtPoint(0)).Mag2());
if(cutMode()==2) dstot=cutLength();
 if(cutMode()==1) dstot=sqrt((tr.LocationAtPoint(tr.LastValidPoint())-tr.LocationAtPoint(0)).Mag2())-cutLength();
 
 int nsegtot=cumseglens.size()-1;
  n=(int)(np/nsegtot);


   if(!n)
     alfa=1;
   else
    alfa=1+0.038*log(dstot/nsegtot/140./cos);
    
alfa=1.;

   //sigma0=0.715;
double washout=0.8585;
   double dxmedio=dstot/double(nsegtot)*10.; //mm
   cout << " weird dxmedio " << dxmedio << " dstot " << dstot << " nsegtot " << nsegtot << " cutmode " << cutMode() << endl;
   dxmedio=140.;
   //washout=1;
thetams=13.6/cp/beta*sqrt(1./140./cos)*alfa/cos*washout*sqrt(dxmedio);

double thetams0=13.6/cp;

    std::cout << " poly beta " << beta << " alfa " << alfa << " cp " << cp << " dxmedio " << dxmedio << std::endl;
    std::cout << " complex thetams  " << thetams << " basic thetams " << thetams0 << std::endl;

   float collPointsRatio=float(hits2dC.size())/float(tr.NPoints());
   float avPointsSeg=float(breakpoints[breakpoints.size()-1] )/breakpoints.size();
float avCollPointsSeg=collPointsRatio*float(breakpoints[breakpoints.size()-1] )/breakpoints.size();

cout << " hits2d " << hits2dC.size() << " n points " << tr.NPoints() << endl;
cout << " last breakpoint " << breakpoints[breakpoints.size()-1] << " n seg " << cumseglens.size() << endl;
cout << " collpointsratio " << collPointsRatio << endl;
cout << " avcollpointsseg " << avCollPointsSeg << endl;

sinb=cosTrackDrift(tr); //3d
//sinb=1; //2d
    thetaerr=sigma0*sqrt(24.)/(dxmedio)/sqrt(float(avCollPointsSeg))/sinb; //2d
    thetaerr=sigma0*sqrt(24.)/(dxmedio)/sqrt(float(avPointsSeg))/sinb; //3d
  // thetaerr=0;
    double thetacorr=sqrt(thetams*thetams+thetaerr*thetaerr);
    // ttall.push_back(acorr/thetams);
     std::cout << " filling ttall fit " << a*1000. << " " << thetams*1000. << " " << thetaerr*1000.  << std::endl;
    ttall.push_back(a/thetacorr);
  
   std::cout << " all that stuff " << thetams << thetaerr << dstot << cos << beta << alfa << cp << ds<< n << ap0 << apmedio << tpmedio << cc<< washout << dxmedio << endl;
 
   FillCovMatrixSegOnly(tr,mat,jp,thetams*thetams,thetaerr*thetaerr,materr,breakpoints);
   std::cout << " before addsegcov nseg " << nseg << std::endl;

  // AddSegmentCovariance(tr,mat,jp);
   std::cout << " after addsegcov " << jp << std::endl;
 //}
 }

cout << " firstseg " << firstseg << " lastseg " << lastseg << endl;

for(unsigned int jp=firstseg;jp<lastseg-1;jp++) {
  std::cout << jp << " breakpoint " << breakpoints[jp] << " momentum " << tr.MomentumAtPoint(breakpoints[jp]) << std::endl;
  std::cout << "elossmode " << eLossMode_ << std::endl;
  float cp=0;  float ds=0;

      ds=sqrt((tr.LocationAtPoint(breakpoints[jp])-tr.LocationAtPoint(breakpoints[jp-1])).Mag2());
 
  const double Etot = sqrt(p0*p0 + MMU*MMU);//Initial energy
  double Eij2 = 0.;
  //double result = 0;
    if (eLossMode_==1) {
      // ELoss mode: MIP (constant)
      constexpr double kcal = 0.002105;
      const double Eij = Etot - kcal*cumseglens[jp];//energy at this segment
      Eij2 = Eij*Eij;
    } else {
      // Non constant energy loss distribution
      const double Eij = GetE(Etot,cumseglens[jp],MMU);
      Eij2 = Eij*Eij;
    }
    //
    if ( Eij2 <= MMU*MMU ) {
      return 0;
    }
    const double pij = sqrt(Eij2 - MMU*MMU);//momentum at this segment
      cp=pij;
      if(cp>15&&ds>0)    
        nbstop++; 


  double a=dthetaPoly[jp-1];

  cout<< a << endl;
  double eseg=0;

  double MMU=105.6;
   beta=cp/sqrt(cp*cp+MMU*MMU);

  alfa=1;
  double sigma0;

  //ComputeD3P();
  sigma0=d3p;
  cout << " eseg " << eseg << " sigma0 " << sigma0 << endl;
  sinb=1;
  
  cos=(collLength()/10./(tr.Length())); //2d
  cos=1; //3d
  //std::cout << " colllength " << collLength() << " trlength " << tr.Length() << " cosine " << cos << std::endl;
 

 if(cutMode()==0) dstot=sqrt((tr.LocationAtPoint(tr.LastValidPoint())-tr.LocationAtPoint(0)).Mag2());
if(cutMode()==2) dstot=cutLength();
 if(cutMode()==1) dstot=sqrt((tr.LocationAtPoint(tr.LastValidPoint())-tr.LocationAtPoint(0)).Mag2())-cutLength();
 
  n=(int)(np/cumseglens.size()-1);
int nsegtot=cumseglens.size()-1;

   if(!n)
     alfa=1;
   else
    alfa=1+0.038*log(dstot/nsegtot/140/cos);

  alfa=1.;
std::cout << " checking n " << n << " alfa " << alfa << " dstot " << dstot << " cos " << cos <<std::endl;
   //sigma0=0.715;
double washout=0.7388;
   double dxmedio=dstot/double(nsegtot)*10.;//mm
   cout << " weird dxmedio " << dxmedio << endl;
   dxmedio=140.;
  // thetaerr=0;
  //  double thetacorr=sqrt(thetams*thetams+thetaerr*thetaerr);
    // ttall.push_back(acorr/thetams);
   //  ttall.push_back(a/thetacorr);
  
   std::cout << " all that stuff again " << thetams << thetaerr << dstot << cos << beta << alfa << cp << ds<< n << ap0 << apmedio << tpmedio << cc<< washout << dxmedio << endl;

   // sigma0=0.715;
    //washout=1;
   thetams=13.6/cp/beta*sqrt(1./140./cos)*alfa/cos*washout*sqrt(dxmedio);
   cout << " thetams poly cp " << cp << " beta " << beta << " alfa " << alfa << " dxmedio " << dxmedio << endl;
   
   float collPointsRatio=float(hits2dC.size())/float(tr.NPoints());
   float avCollPointsSeg=collPointsRatio*float(breakpoints[breakpoints.size()-1] )/breakpoints.size();
   float avPointsSeg=float(breakpoints[breakpoints.size()-1] )/breakpoints.size();


cout << " hits2d " << hits2dC.size() << " n points " << tr.NPoints() << endl;
cout << " last breakpoint " << breakpoints[breakpoints.size()-1] << " n seg " << cumseglens.size() << endl;
cout << " collpointsratio " << collPointsRatio << endl;
cout << " avcollpointsseg " << avCollPointsSeg << endl;
sinb=cosTrackDrift(tr);  //3d
//sinb=1; //2d
    thetaerr=sigma0*sqrt(6.)/(dxmedio)/sqrt(float(avCollPointsSeg))/sinb; //2d
    thetaerr=sigma0*sqrt(6.)/(dxmedio)/sqrt(float(avPointsSeg))/sinb; //3d
cout << " thetaerr " << thetaerr << endl;
cout << " sigma0" << sigma0 << " dxmedio " << dxmedio << " sqrt " << sqrt(float(avCollPointsSeg)) << " sinb " << sinb << endl;
    cout << " npnsegtot ratio " << np/nsegtot <<" np " << np << " nsegtot " << nsegtot << endl;
  // thetaerr=0;
   double thetacorr=sqrt(thetams*thetams+thetaerr*thetaerr);
   std::cout << " beta " << beta << " alfa " << alfa << " cp " << cp << " dstot " << dstot << std::endl;

     std::cout << " filling ttall poly " << a*1000. << " " << thetams*1000. << " " << thetaerr*1000.  << std::endl;
   ttall.push_back(a/thetacorr);
    
     //apmedio+=abs(a*a);
     // ap0+=a;
     //tpmedio+=(thetams*thetams*d+thetaerr*thetaerr*6/d/d);
     // cc+=(a*a)/(thetams*thetams*d+6*thetaerr*thetaerr/d/d);

     	 
 
 
  //cout << " poly thetams " << thetams << " thetaerr " << thetaerr << endl;
	 cout << " before covma "  << jp << endl;
   //FillCovMatrix(tr,matpoly,jp,thetams*thetams,thetaerr*thetaerr,matpolyms,matpolyerr,breakpoints,firstseg);
	 cout << "before mixx "  << jp << endl;
   //FillCovMixTerms(tr,matmix,jp,nseg-2,thetams*thetams,thetaerr*thetaerr);
	 cout << "after mixx "  << jp << endl;
	 //AddSegmentCovariance(tr,matpoly,jp);
	 std::cout << " after fillmatrix1  " << jp << std::endl; 
for(int jm=0;jm<matpoly.GetNrows();jm++)
for(int jmm=0;jmm<matpoly.GetNrows();jmm++)
  mat(jm+nseg-1,jmm+nseg-1)=matpoly(jm,jmm);
//if(jp==np-2)
 std::cout << " after matpoly  " << jp << " np " << nseg << std::endl; 
}

//mat+=matmix;

TMatrixD cov(TMatrixD::kUnit,mat);

TMatrixD vtall(ttall.size(),1);
for(unsigned int jv=0;jv<ttall.size();jv++)
  vtall(jv,0)=ttall[jv]; 

TMatrixD invcov=cov.Invert();
TMatrixD vtcov=invcov*vtall;
TMatrixD tvtall(TMatrixD::kTransposed,vtall); 

 double vtcovmed=0;
 std::vector<double> terms;
 for(unsigned int jv=0;jv<ttall.size();jv++) {
 double term=tvtall(0,jv)*vtcov(jv,0);
 terms.push_back(term);
 vtcovmed+=term;
 }
 vtcovmed/=(ttall.size());
 
 
 
 // compute chi2 before truncation

vector<double> ttrunc;
vector<int> tails;
double ttrunctot=0;
vector<double> atrunc;
//nsigma=999;
 int nsigma=3;
for (unsigned int jt=0;jt<ttall.size();jt++) {
  cout << " threshold " << vtcovmed*nsigma*nsigma << " term " << terms[jt] << endl;
  if(terms[jt]<vtcovmed*nsigma*nsigma) {
//   {
    ttrunctot+=terms[jt];
    atrunc.push_back(ttall[jt]);
    ttrunc.push_back(terms[jt]);
    std::cout << " keeping value " << jt << std::endl;
   }
    else {
      tails.push_back(jt);
    }
  }
   
  TMatrixD vtrunc(atrunc.size(),1);
 for(unsigned int jv=0;jv<atrunc.size();jv++)
   vtrunc(jv,0)=atrunc[jv];

std::cout << " after vtrunc tails size " << tails.size() << std::endl;

   std::vector<TMatrix> covs;
   covs.push_back(cov);
 if(tails.size()) {
    for(int jta=tails.size()-1;jta>=0;jta--) {   
    cout << " cleaning covariance jta " << jta << endl; 
     TMatrix covtemp=CleanCovariance(covs[covs.size()-1],tails[jta]); 
     cout << " after cleaning jta " << jta << endl;
covs.push_back(covtemp);

}
}
TMatrix covcut=covs[covs.size()-1];
std::cout << " after covcut " << std::endl;
 if(!covcut.GetNrows())
   return -999;

std::cout << " before invcovmod " << std::endl;
 TMatrixD invcovmod=covcut.Invert();
std::cout << " vtrunc " << std::endl;
vtrunc.Print();
 TMatrixD vtcovmod=invcovmod*vtrunc;
std::cout << "  vtcovmod " << std::endl;
vtcovmod.Print();
TMatrixD tvtrunc(TMatrixD::kTransposed,vtrunc); 
std::cout << " after tvtrunc " << std::endl;

 double vmediomod=0;
 std::vector<double> tterms;
 for(unsigned int jv=0;jv<ttrunc.size();jv++) {
  cout << " jv " << jv << endl;
 double tterm=tvtrunc(0,jv)*vtcovmod(jv,0);
 cout << " adding term " << tterm << endl;
 tterms.push_back(tterm);
 vmediomod+=tterm;
 }
 std::cout << " after vmediomod " << std::endl;

 vmediomod/=(ttrunc.size());

 double vtrunctot=0;

  for ( int jt=0;jt<vtrunc.GetNrows();jt++) 
    vtrunctot+=vtrunc(jt,0);
 std::cout << " after vtrunctot " << std::endl;

//int ndf=ttall.size()-tails.size();

double ttnewrms=vmediomod;
 std::cout << " ttnewrms " << ttnewrms << std::endl; 
 //exit(11);
 // ncut=ttall.size()-vtnew.GetNrows();

      
//GaussCorrection* gaus=new GaussCorrection();
//double gcorr=gaus->GaussFactor(ttcov.size(),1./sqrt(3.));
//gau=gcorr;

 //gau=1.;
 
//    return ttnewrms;
    

  return ttnewrms;
}
/////////////////////////////////

const void TrajectoryMCSFitterICARUS::FillCovMatrixSegOnly(recob::TrackTrajectory tr, TMatrixDSym mat,unsigned int jp,double sms,double serr,TMatrixDSym materr,std::vector<long unsigned int> breaks) const
{
   /*
  std::cout << " begin segonly " << std::endl;
 unsigned int ip=jp-1;
auto pos0 = tr.LocationAtPoint(ip);
auto pos1 = tr.LocationAtPoint(jp);
double dx= ( (pos1-pos0).R() );
auto pos2 = tr.LocationAtPoint(jp+1);
double dxp= ( (pos2-pos1).R() );
 double dxmed=(dxp+dx)/2;
if(dx<0.00001)
  dxmed=dxp;
  c0=sms*dxmed+serr;


 double s0=1/double(breaks[jp]-breaks[ip]);
double sp=1/double(breaks[jp+1]-breaks[jp]);
std::cout << " segonly s0 " << s0 << " sp " << sp << std::endl;


// double c0=1;
 //double cp;//,cm,cpp,cmm;


 double wash1=0.174;
 double wash0=0.859;

// double serrmod=serr*sp;

 //c0=sms*dx*wash0+serr/dx/dx*s0+serr*sp/dx/dx;
 //c0=sms*dxmed;
 
  std::cout << " before filling mat "  << std::endl;
  
  mat(jp,jp)=c0;
  

   std::cout << " after filling mat "  << std::endl;
   
  //materr[ip][ip]=serr/dx/dx*s0+serr*sp/dx/dx;
//cout << " diagonal error " << serr/dx/dx*s0+serr*sp/dx/dx+tcard << " diagonal tcard " << tcard << endl;
 std::cout << " filling covsegonly MS " << jp << "," << jp << " : " << sms*dx*wash0 << std::endl; 
// cout << " filling covsegonly MS " << jp << "," << jp << " : " << sms*dx*wash0 << endl; 
//cout << " filling covsegonly err " << jp << "," << jp << " : " << serr/dx/dx*s0+serr*sp/dx/dx<< endl; 

//cout << " filling covsegonly tcard " << jp << "," << jp << " : " << sqrt(tcard) << " dx " << dx << endl; 

//cout << " np " << np << endl;
 if(jp<tr.NPoints()-1) {


double cp=-serr*sp/dx/dxp+sms*wash1*dx;

 std::cout << " offdiagonal error " << serr*sp/dx/dxp << std::endl;
 mat(ip,ip+1)=cp;
if(ip==0) {
  //cout << "segonly serrsp " << serr*sp <<  " sp " << sp << " dx " << dx << " dxp " << dxp <<endl;
  //cout << "segonly product "<< -serr*sp/dx/dxp << endl;
}
  }


//SYMMETRIZE
  if(jp<tr.NPoints()) {
    //mat[ip+1][ip]=mat[ip][ip+1];
}
//debu << mat << endl;
  std::cout << " end segonly " << std::endl;


 double dxm,dxmm,dxp,dxpp;

  int ip=jp-1;
 dxm=point(jp-1).DS();
 dxp=point(jp).DS();
 // dxm=dxmedio;
 //dxp=dxmedio;
 cout <<" point " << jp << " dxm " << dxm << " dxp " << dxp << endl;


 double c0,cp,cm,cpp,cmm;
 double dxmed=(dxp+dxm)/2;
 if(dxm<0.00001)
  dxmed=dxp;
  c0=sms*dxmed+serr;
 //c0=sms*dxmed;
 cout << " serr " << serr << " dxp " << dxp << " dxm " << dxm << " dxmed " << dxmed << endl;
 cout << " MS term " << sms*dxmed << " measurement term " << 2*serr << endl;
 cout << "measMS ratio" << 2*serr/(sms*dxmed) << endl;
  mat[ip][ip]=c0;

  cout << " filling covseg matrix " << jp << "," << jp << " : " << c0 << endl; 
cout << " np " << np << endl;
  if(jp<np+1) {
    //dxp=point(ip+1).DS();
dxpp=point(jp+1).DS();
 //dxpp=dxmedio;
double errp=point(jp+1).SMSeg(); 
double errm=point(jp).SMSeg(); 
 cp=-abs(errm*errp);
  cout << " serr " << serr << " sqrt " << sqrt(serr) << " errp " << errp << endl;
 cout << " cp " << cp << endl;
 mat[ip][ip+1]=cp;
  cout << " filling covseg matrix " << jp << "," << jp+1 << " : " << cp << endl; 
  }

 if(jp>1) {
 dxmm=point(jp-1).DS();
 //dxmm=dxmedio;
 // if(jp==21||jp==22) {
 //cout << " check- element " << jp << " point " << ip << " dxp " << dxp <<" dxm " << dxm << " dxmm " << dxmm << endl;
 //cout << " term1 " << 2/dxm/dxm << " term2 " << 1/dxm/dxp << " term 3 " << 1/dxm/dxmm << endl;
 // }
double errm=point(jp-1).SMSeg(); 
double errp=point(jp).SMSeg(); 
 cm=-abs(errp*errm);
   cout << " filling covseg matrix " << jp << "," << jp-1 << " : " << cm << endl; 
   mat[ip][ip-1]=cm;
 }
 */
}
 
///////////////////////////////////////////////////////////////////
const void TrajectoryMCSFitterICARUS::AddSegmentCovariance(recob::TrackTrajectory tr,TMatrixDSym mat, int jp) const
{
int ip=jp-1;
int np=mat.GetNrows();
TMatrixDSym matbd(np);
 std::cout << " after creating matbd " << std::endl;

  TVector3 start=tr.LocationAtPoint<TVector3>(0);
  TVector3 end=tr.LocationAtPoint<TVector3>(tr.LastValidPoint());
  TVector3 avdir=end-start;
  TVector3 projColl(0,0,0);
  projColl[0]=1; projColl[1]=sqrt(3.)/2.; projColl[2]=0.5;

double t0p=0;
double t0m=0;

double tppp=0;

//cout << " jm " << jm << " nboarDX " << point(jm).NBoarDX() << endl;
//cout << " jm+1 " << jm+1 << " nboarDX " << point(jm+1).NBoards() << endl;
//for(int jb1=0;jb1<point(jm).NBoards();jb1++)
// cout << " jm nseg " <<  point(jm).NBSeg(jb1) << endl;
//for(int jb2=0;jb2<point(jm+1).NBoards();jb2++)
// cout << " jm+1 nseg " <<  point(jm+1).NBSeg(jb2) << endl;

 double dxpp=-1;
auto pos0 = tr.LocationAtPoint(jp-1);
auto pos1 = tr.LocationAtPoint(jp);
auto pos2 = tr.LocationAtPoint(jp+1);
double dxm= ( (pos1-pos0).R() );
double dxp= ( (pos2-pos1).R() );
//auto pos3=pos2;
 if(jp<mat.GetNrows()) {
auto pos3 = tr.LocationAtPoint(jp+2);
 dxpp= ( (pos3-pos2).R() );
 }
double sm,sp,s0;
//sm=point(jm-1).Sin2B();
//s0=point(jm).Sin2B();
//sp=point(jm+1).Sin2B();
 sm=1; s0=1; sp=1;

 //double t00=point(jm).ESeg();
 //double tpp=point(jm+1).ESeg();
 //double tmm=point(jm-1).ESeg();
 double t00=0; double tpp=0; double tmm=0;

  
//cout << "jm-1: " << jm-1 << " t00 " << t00 << " tpp " << tpp << " tmm " << tmm << endl;
//cout << "diag jm-1: " << jm-1 << " e00 " << point(jm).ESeg()<< " epp " << point(jm+1).ESeg() << " tmm " << point(jm-1).ESeg() << endl;
double t1,t2,t3;
t1=t00*t00*s0*s0*(1/dxp+1/dxm)*(1/dxp+1/dxm);
t2=tpp*tpp*sp*sp/dxp/dxp;
t3=tmm*tmm*sm*sm/dxm/dxm;

//if(jm==3) {
cout << " t1 " << t1 << " t2 " << t2 << " t3 " << t3 << " ip " << ip << endl;
//cout << " DIAG0 " << t00*t00	 << " DIAG1 " << tpp*tpp << " DIAG2 " << tmm*tmm << endl;
//cout << " ANG0 " << (1/dxp+1/dxm)*s0 << " ANG1 " << 1/dxp/dxp*s0*s0 << " ANG2 " << 1/dxm/dxm*s0*s0 << endl;
//cout << " s0 " << point(jm).Sin2B() << " dxm " << dxm << " dxp " << dxp << endl;
//}


matbd(ip,ip)=(t1+t2+t3);
matbd(ip,ip)+=(-2*t0p*t0p*s0*sp/dxp*(1/dxp+1/dxm)-2*t0m*t0m*s0*sm/dxm*(1/dxp+1/dxm));
 cout << " t4 " << (-2*t0p*t0p*s0*sp/dxp*(1/dxp+1/dxm)-2*t0m*t0m*s0*sm/dxm*(1/dxp+1/dxm)) << endl;

 if(jp<mat.GetNrows()){
   //double spp=point(jm+2).Sin2B();
   /// double spp=1;
   cout << " ip " << ip << " ip+1 " << ip+1 <<" matbd size " << matbd.GetNrows() << " " << matbd.GetNcols() << endl;
  matbd(ip,ip+1)=-t00*t00*s0*s0/dxp*(1/dxm+1/dxp)-tpp*tpp*sp*sp/dxp*(1/dxpp+1/dxp);
cout << " board off " << ip << " " << ip+1 << endl;
matbd(ip,ip+1)+=(tppp*tppp*sp*s0/dxpp/dxp+t0m*t0m*s0*sp/dxp/dxm+t0p*t0p*s0*sp*(1/dxp/dxp+(1/dxp+1/dxm)*(1/dxp+1/dxpp)));
//cout << " board off term 1 " << tppp*tppp*sp*s0/dxpp/dxp << " term2 " << t0m*t0m*s0*sp/dxp/dxm << " term3 " << t0p*t0p*s0*sp*(1/dxp/dxp+(1/dxp+1/dxm)*(1/dxpp+1/dxp)) << endl;
}
 if(jp<mat.GetNrows()-1) {
   cout << " plus2 ip " << ip << " jp " << jp <<" matbd size " << matbd.GetNrows() << endl;
  matbd(ip,ip+2)=tpp*tpp*sp*sp/dxp/dxpp;
  //double spp=point(jp+2).Sin2B();
  double spp=1;
auto pos3 = tr.LocationAtPoint(jp+2);
auto pos4 = tr.LocationAtPoint(jp+3);
double dxppp= ( (pos4-pos3).R() );
 matbd(ip,ip+2)+=-t0p*t0p*s0*spp/dxpp*(1/dxp+1/dxm)-tppp*tppp*s0*spp/dxp*(1/dxppp+1/dxpp);
 //cout << " board offoff " << ip << " " << ip+2 <<  " term1 " << t0p*t0p*s0*sp/dxp*(1/dxp+1/dxm) << " t0p2 " << t0p*t0p << " angle " << s0*sp/dxp*(1/dxp+1/dxm) << endl;
}
 if(jp<mat.GetNrows()-2) {
   //double spp=point(jm+2).Sin2B();
   double spp=1;
 matbd(ip,ip+2)=tppp*tppp*sp*spp/dxp/dxpp;

}
cout << " before symmetrization " << endl;
 if(jp<mat.GetNrows())
matbd(ip+1,ip)=matbd(ip,ip+1);
 if(jp<mat.GetNrows()-1) 
matbd(ip+2,ip)=matbd(ip,ip+2);
 if(jp<mat.GetNrows()-2) 
matbd(ip+3,ip)=matbd(ip,ip+3);

mat+=matbd;

}
const void TrajectoryMCSFitterICARUS::FillCovMatrix(recob::TrackTrajectory tr,TMatrixDSym mat,int jp,double sms,double serr, TMatrixDSym matms, TMatrixDSym materr,std::vector<long unsigned int> breakpoints, int firstseg) const
{

  int matsize=mat.GetNrows();
  //cout << " before filling MATMS " << matms << endl;
 //cout << " sms " << sms << " serr " << serr <<  " np " << np << endl;

   //cout << " point " << jp << " DS " << point(jp).DS() << endl;

 //double dxmedio=76.7179;

  double MMU=105.6;
auto pos0 = tr.LocationAtPoint(jp-1);
auto pos1 = tr.LocationAtPoint(jp);
auto pos2 = tr.LocationAtPoint(jp+1);
double dxm= ( (pos1-pos0).R() );
double dxp= ( (pos2-pos1).R() );
  int ip=jp-firstseg;

double pp,ppp,p0;
p0=tr.MomentumAtPoint(jp);
//pm=tr.MomentumAtPoint(jp-1);
pp=tr.MomentumAtPoint(jp+1);

 double bp,b0;
 //bm=pm/sqrt(pm*pm+MMU*MMU);
 b0=p0/sqrt(p0*p0+MMU*MMU);
bp=pp/sqrt(pp*pp+MMU*MMU);

 double em,ep,e0;
 em=breakpoints[jp-1]-breakpoints[jp-2];
 e0=breakpoints[jp]-breakpoints[jp-1];
 ep=breakpoints[jp+1]-breakpoints[jp];

 double sm,sp,s0;
 sm=1; s0=1; sp=1;
 //sm=point(jp-1).Sin2B();
 //s0=point(jp).Sin2B();
 //sp=point(jp+1).Sin2B();

 
 // dxm=dxmedio;
 //dxp=dxmedio;
 std::cout << " element " << jp << " point " << ip << " dxm " << dxm << " dxp " << dxp << std::endl;

 double c0,cp,cpp;
// c0=sms/p0/p0/b0/b0*(dxp+dxm)/2+2*serr/e0/e0*s0*s0*(1/dxp/dxp+1/dxm/dxm+1/dxp/dxm);
c0=sms/p0/p0/b0/b0*(dxp+dxm)/2+serr*(s0*s0/e0/e0*(1/dxp+1/dxm)*(1/dxp+1/dxm)+sp*sp/ep/ep/dxp/dxp+sm*sm/em/em/dxm/dxm);
 std::cout << " c0 " << c0 << " ip " << ip << std::endl;
 mat(0,0)=c0;
 std::cout << " c0 " << c0 << " ip " << ip << std::endl;
 matms(ip,ip)=sms/p0/p0/b0/b0*(dxp+dxm)/2;
 std::cout << " c0 " << c0 << " ip " << ip << std::endl;
 materr(ip,ip)=serr*(s0*s0/e0/e0*(1/dxp+1/dxm)*(1/dxp+1/dxm)+sp*sp/ep/ep/dxp/dxp+sm*sm/em/em/dxm/dxm);
//point(jp).setMMR(serr*(s0*s0/e0/e0*(1/dxp+1/dxm)*(1/dxp+1/dxm)+sp*sp/ep/ep/dxp/dxp+sm*sm/em/em/dxm/dxm)/(sms*(dxp+dxm)/2));
 std::cout << " before 1 " << std::endl;
  if(jp<matsize-1) {
    //dxp=point(ip+1).DS();
auto pos3 = tr.LocationAtPoint(jp+2);
double dxpp= ( (pos3-pos2).R() );
 //dxpp=dxmedio;

//double p0p=sqrt(p0*pp)
 double smsp=sms/p0/pp/b0/bp;
double errp=serr*(s0*s0/e0/e0/dxp*(1/dxm+1/dxp)+sp*sp/ep/ep/dxp*(1/dxpp+1/dxp));

  cp=-errp+smsp*dxp*0.393;
  mat(jp,jp+1)=cp;

matms(jp,jp+1)=smsp*(dxp)*0.393;
materr(jp,jp+1)=-errp;

  }
  std::cout << " before -1 matsize " << matsize << " jp " << jp << " ip " << ip << std::endl;
  if(jp<matsize-2) {
auto pos3 = tr.LocationAtPoint(jp+2);
double dxpp= ( (pos3-pos2).R() );
//ppp=point(jp+2).CP();
ppp=tr.MomentumAtPoint(jp+2);
double bpp=ppp/sqrt(ppp*ppp+MMU*MMU);
//double epp=breakpoints[jp+2]-breakpoints[jp+1];
 //spp=point(jp+2).Sin2B();
//double spp=1;
double smspp=sms/p0/ppp/bpp/b0;
double errpp=serr*sp*sp/ep/ep;
 // dxpp=dxmedio;
  cpp=errpp/dxp/dxpp;
//double smspp=sms*p0/sqrt(ppp*p0);
 mat(jp,jp+2)=cpp+smspp*(dxp+dxpp)/2*0.0148;
matms(jp,jp+2)=smspp*(dxp+dxpp)/2*0.0148;
materr(jp,jp+2)=cpp;
 
  }
  std::cout << " before symmetrize " << std::endl;

//SYMMETRIZE
 if(jp<matsize-1) {
mat(jp+1,jp)=mat(jp,jp+1);
matms(jp+1,jp)=matms(jp,jp+1);
materr(jp+1,jp)=materr(jp,jp+1);
}
if(jp<matsize-2) {
mat(jp+2,jp)=mat(jp,jp+2);
matms(jp+2,jp)=matms(jp,jp+2);
materr(jp+2,jp)=materr(jp,jp+2);
}

}
 ///////////////////////////////////////////////////////////
const void TrajectoryMCSFitterICARUS::FillCovMixTerms(recob::TrackTrajectory tr, TMatrixDSym mat,int jp,int ns,double sms,double serr) const
{
//int np=mat.GetNrows();
 cout << " sms " << sms << " serr " << serr << endl;

 //double dxmedio=76.7179;
double p0=tr.MomentumAtPoint(jp);


 double dxm,dxp;

  int ip=jp-1;
 

auto pos0 = tr.LocationAtPoint(jp-1);
auto pos1 = tr.LocationAtPoint(jp);
auto pos2 = tr.LocationAtPoint(jp+1);
 dxm= ( (pos1-pos0).R() );
 dxp= ( (pos2-pos1).R() );



 // dxm=dxmedio;
 //dxp=dxmedio;
 cout << " element " << jp << " point " << ip << " dxm " << dxm << " dxp " << dxp << endl;

double wash0=0.735;
double wash1=0.055;
double washfit=0.859;
 double c0;
 if(ip<ns-1) {

double  c0=sms*(dxp+dxm)/2/p0/p0*washfit*wash0;
 std::cout << " mix term a " << ip  << " " << ip+ns-2 << " value " << c0 << std::endl; 
  std::cout << " jp " << jp  << " ip " << ip << " ns " << ns << std::endl; 

 mat(ip,ip+ns-2)=c0;
 //mat[ip+ns-2][ip]=c0;
 //cout << " mix term " << ip  << " " << ip+ns-2 << " value " << c0 << endl; 
 //cout << " mix term " <<  ip  << " " << ip+ns-2 << " sms " << sms << " dx " << (dxp+dxm)/2 << endl;
 }
  if(ip<ns-2) {
    std::cout << " filling covmix matrix b" << jp << "," << jp << " : " << c0 << std::endl; 
     cout << " mix term b " << ip  << " " << ip+ns-1 << " value " << c0 << endl; 
      std::cout << " jp " << jp  << " ip " << ip << " ns " << ns << std::endl; 

 c0=sms*(dxp+dxm)/2/p0/p0*washfit*wash0;

  mat(ip,ip+ns-2)=c0;
 
}
if(ip<ns-3) {
   cout << " filling covmix matrix c " << jp << "," << jp << " : " << c0 << endl; 
   c0=sms*(dxp+dxm)/2/p0/p0*washfit*wash1;

   mat(ip,ip+ns-3)=c0;
   //mat[ip+ns][ip]=c0;
   cout << " mix term c " << ip  << " " << ip+ns << " value " << c0 << endl; 
   // cout << " mix term " << ip+ns  << " " << ip << " value " << c0 << endl; 
 }
 if(ip<ns) {
   cout << " filling covmix matrix d " << jp << "," << jp << " : " << c0 << endl; 
   c0=sms*(dxp+dxm)/2/p0/p0*washfit*wash1;

   mat(ip,ip+ns)=c0;
   //mat[ip+ns-3][ip]=c0;
   cout << " mix term d " << ip  << " " << ip+ns-3 << " value " << c0 << endl; 
 }


}
/**********************************************************/
const TrajectoryMCSFitterICARUS::ScanResult TrajectoryMCSFitterICARUS::C2Fit(std::vector<float>& dtheta,std::vector<float>& dthetaPoly, std::vector<float>& seg_nradlengths, std::vector<float>& cumLen,std::vector<size_t>& breaks, bool fwdFit, bool momDepConst, int pid, float sigma, const recob::TrackTrajectory& traj) const
  //fit C2(p) function with expected dependency (see ICARUS MCS paper)
{
  
//if(!fwdFit) exit(22);
int nMom=ceil((pMax_-pMin_)/pStep_)+1;

 TVectorD wmom(nMom);
 TVectorD wsmom(nMom);
 TVectorD wsigma(nMom);
 TVectorD wc2(nMom);
std::cout << " begin c2fit " << nMom << std::endl;

int jMom=0;
std::cout << " before loop " << std::endl;
double firstValid=0;
for (double p_test = pMin_; p_test <= pMax_; p_test+=pStep_) {
//for (double p_test = 500; p_test <= 500; p_test+=pStep_) {
std::cout << " before c2function " << std::endl;
std::vector<float> ttall; ttall.clear();
   // double c2 = C2Function(p_test, angResol_, dtheta, seg_nradlengths, cumLen, breaks,fwdFit, momDepConst, pid, 0,traj);
   double c2= C2Function(traj,cumLen,breaks,dtheta,dthetaPoly,ttall,p_test); 
std::cout << " p_test " << p_test << " c2 " << c2 << std::endl;
std::cout << " jmom " << jMom << " nMom " << nMom << std::endl;
    wmom[jMom]=p_test/1000.;
 wc2[jMom]=c2; wsigma[jMom]=sigma*c2; wsmom[jMom]=1;
   if(firstValid<0.01&&c2>0.01) firstValid=jMom;
    jMom++;
}
std::cout << " after momentum loop " << jMom << std::endl;
TVectorD rmom(nMom-firstValid);
 TVectorD rsmom(nMom-firstValid);
 TVectorD rsigma(nMom-firstValid);
 TVectorD rc2(nMom-firstValid);
for(int jp=firstValid;jp<jMom;jp++) {
rmom[jp-firstValid]=wmom[jp];
rsmom[jp-firstValid]=wsmom[jp];
rsigma[jp-firstValid]=wsigma[jp];
rc2[jp-firstValid]=wc2[jp];
}
 
const TVectorD cmom=rmom;
const TVectorD cc2=rc2;
const TVectorD csmom=rsmom;
const TVectorD csigma=rsigma;

cout << " firstValid " << firstValid << " wmom " << wmom[firstValid] << endl;
cmom.Print();
cc2.Print();


 TGraphErrors *gr3 = new TGraphErrors(cmom,cc2,csmom,csigma);
 TF1* fitfunc= new TF1("fitfunc",funzio,wmom[firstValid],pMax_,2);
fitfunc->SetParLimits(0,-0.5,1.);
fitfunc->SetParLimits(1,-0.,5.);

      gr3->Fit("fitfunc"); 

    double fixedTerm,momTerm;
    double eFixedTerm,eMomTerm;
//fit function parameters and errors (see paper)
            fixedTerm=fitfunc->GetParameter(0);
            momTerm=fitfunc->GetParameter(1); 
	    eFixedTerm=fitfunc->GetParError(0);
            eMomTerm=fitfunc->GetParError(1); 

 

  gr3->Delete();
  fitfunc->Delete();
  double best_p=sqrt(momTerm/(1-fixedTerm));
  double error_p=sqrt(pow(eMomTerm,2.)/(4*momTerm*(1-fixedTerm))+momTerm*pow(eFixedTerm,2.)/(4*pow((1-fixedTerm),3.)));
  std::cout << " end c2fit best_p" << best_p << std::endl;
  if(best_p<.001) { //naive interpolation
  int over=-1;
  if(rc2[0]>1) over=0;
 for(int jp=1;jp<nMom-firstValid;jp++) {
if(rc2[jp]>1.&&rc2[jp-1]<1.) over=jp;
 }
 if(over>0&&over<nMom-firstValid) best_p=rmom[over-1]+(rmom[over]-rmom[over-1])/(rc2[over]-rc2[over-1]);
 if(over==0) best_p=rmom[0]; //underflow
if(over==-1) best_p=rmom[nMom-firstValid-1]; //overflow
  std::cout << " end c2fit naive best_p" << best_p << std::endl;

  }

  
 // exit(11);
  
  //float best_p=0; float error_p=0;
  return ScanResult(best_p, error_p, 0.);
}

 
double TrajectoryMCSFitterICARUS::collLength() const
{
recob::Hit h0=hits2dC[0];
int jf=hits2dC.size()-1;
int NS0=h0.WireID().TPC;
//float x0=h0.WireID().Wire*3; auto y0=h0.PeakTime()*0.622;
int jb=jf;
for(unsigned int jh=0;jh<hits2dC.size();jh++) {
  int NS=hits2dC[jh].WireID().TPC;
  if(NS!=NS0)
   jb=jh-1;
}
recob::Hit hf=hits2dC[jf];
recob::Hit hb=hits2dC[jb];
recob::Hit hbp=hits2dC[jb+1];

int w0=h0.WireID().Wire;
int wf=hf.WireID().Wire;
int wb=hb.WireID().Wire;
float dxf=(wb-w0)*3.; 
float dyf=(hb.PeakTime()-h0.PeakTime())*0.622; 
std::cout << " dxf " << dxf << " dyf " << dyf << std::endl;
std::cout << " jf " << jf << " jb " << jb << std::endl;
std::cout << " wf " <<wf << " wb " << wb<< " w0 " << w0<<std::endl;
if(jb!=jf) {
float dxfp=(hf.WireID().Wire-hbp.WireID().Wire)*3.; 
float dyfp=(hf.PeakTime()-hbp.PeakTime())*0.622;  
dxf+=dxfp;
dyf+=dyfp;
std::cout << " dxfp " << dxfp << " dyfp " << dyfp << std::endl;

}
std::cout << " corrected dxf " << dxf << " dyf " << dyf << std::endl;

return sqrt(dxf*dxf+dyf*dyf);
}
double TrajectoryMCSFitterICARUS::collWireLength()
{
recob::Hit h0=hits2dC[0];
recob::Hit hf=hits2dC[hits2dC.size()-1];
float x0=h0.WireID().Wire*3; 
float xf=hf.WireID().Wire*3; 
return xf-x0;
}
 double TrajectoryMCSFitterICARUS::cosTrackDrift( recob::TrackTrajectory tr) const
{
auto start=tr.Start();  
auto end=tr.End();
auto l3d=sqrt((end-start).Mag2()); 
auto l1d=(end-start).X();
return l1d/l3d;

}
/******************************************************************/
TMatrixD TrajectoryMCSFitterICARUS::ReferenceFrame(int plane, int tpc, int cryo) const  { 
  //! Initialize reference frames relative to wire planes (views)
cout << " begin reference frame " << plane << " " << tpc << endl;
 TVector vec1(3);TVector vec2(3);TVector vec3(3);

//const double PITCH=2.99;
const double PIG=3.1416;

if(plane==0) {
  //view 0 (induction 1)
  vec1(0)=1;vec1(1)=vec1[2]=0;
  vec2(1)=1;vec2(0)=vec2[2]=0;
  vec3[2]=1;vec3(0)=vec3(1)=0;
 
}

if(plane==1) {
   //view 1 (induction 2)
  if(tpc>1) { 

  vec1(0)=-cos(PIG/3);vec1(1)=0;vec1[2]=sin(PIG/3);
  vec2(1)=-1;vec2(0)=vec2[2]=0;
  vec3(0)=sin(PIG/3);vec3(1)=0;vec3[2]=cos(PIG/3);
  }
  else {

  vec1(0)=cos(PIG/3);vec1(1)=0;vec1[2]=sin(PIG/3);
  vec2(1)=1;vec2(0)=vec2[2]=0;
  vec3(0)=sin(PIG/3);vec3(1)=0;vec3[2]=-cos(PIG/3);
  }
}
  if(plane==2) {
  //view 2 (collection)
    if(tpc>1) { 
 
 vec1(0)=cos(PIG/3); vec1(1)=0;vec1[2]=sin(PIG/3);
  vec2(1)=1;vec2(0)=vec2[2]=0;
  vec3(0)=sin(PIG/3);vec3(1)=0;vec3[2]=-cos(PIG/3);
 
  
  }
  else {

  vec1(0)=-cos(PIG/3);vec1(1)=0;vec1[2]=sin(PIG/3);
  vec2(1)=1;vec2(0)=vec2[2]=0;
  vec3(0)=sin(PIG/3);vec3(1)=0;vec3[2]=cos(PIG/3);
 
  }
  } 
/*SISTEMARE DIVERSE TPC LOGICHE IN CASO DI INDUZIONE 1!


  //view 3 (half of induction 1)
 
  vec1(0)=-1;vec1(1)=vec1[2]=0;
  vec2(1)=1;vec2(0)=vec2[2]=0;
  vec3[2]=-1;vec3(0)=vec3(1)=0;

  ori(0)=2112*PITCH;  

ReferenceFrame rf3(ori,vec1,vec2,vec3);

  rf3.setOrigin(ori); 
  rf3.setRefMatrix();

  addRefFrame(rf3);
*/
TMatrixD mat(3,3);
//cout << " reference frame mat before fill " << endl;
for(int j=0;j<3;j++) {
mat(0,j)=vec1(j);
mat(1,j)=vec2(j);
mat(2,j)=vec3(j);
}
//cout << " reference frame mat after fill " << endl;
//mat.Print();
return mat;
  }
  /******************************************************************/
double TrajectoryMCSFitterICARUS::DriftOrigin(int plane, int tpc, int cryo) const  { 
double x0;
if(cryo==1&&tpc>1) x0=3596.3;
if(cryo==0&&tpc<=1) x0=-3596.3;
if(cryo==1&&tpc<=1) x0=608.0;
if(cryo==0&&tpc>1) x0=-608.0;
  return x0;
}
TMatrix TrajectoryMCSFitterICARUS::CleanCovariance(TMatrix  cov, int jtail) const
{
  int ncov=cov.GetNrows();
  cout << " clean ncov " << ncov << endl;
  cout << " clean jtail " << jtail << endl;

  TMatrix covmod(ncov-1,ncov-1);

  for(int jrow=0;jrow<ncov;jrow++) {
    for(int jcol=0;jcol<ncov;jcol++) {
      if(jrow<jtail&&jcol<jtail) {
	cout << " writingcov " << jrow << " " << jcol << endl;
	covmod(jrow,jcol)=cov(jrow,jcol);
	}
          if(jrow>jtail&&jcol>jtail){
	    cout << " movingcov " << jrow << " " << jcol << endl;
      	covmod(jrow-1,jcol-1)=cov(jrow,jcol);
	}
    }}

cout << " cleaned matrix " << endl;
covmod.Print();
return covmod;
 
}