#include "TrajectoryMCSFitterICARUS.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardata/RecoBaseProxy/Track.h"

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
recob::MCSFitResult TrajectoryMCSFitterICARUS::fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst) const {
  //notify algorithm has started
  cout << "starting icarus fit..." << endl;
  cout << " " << endl;

  //print cut mode and specify in what it consist
  cout << "cut mode = " << cutMode_ << endl;
  if (cutMode_ == 0) cout << "this means we are not cutting track" << endl;
  else if (cutMode_ == 1) cout << "this means we are cutting last 3 angles of track" << endl;
  else if (cutMode_ == 2) cout << "this means we are cutting nothing but first 7 angles of track" << endl;
  cout << " " << endl;

  //print dim mode and specify in what it consist
  cout << "dim mode = " << dimMode_ << endl;
  if (dimMode_ == 2) cout << "this means we are fitting 2D hits" << endl;
  else if (dimMode_ == 3) cout << "this means we are fitting 3D hits" << endl;
  cout << " " << endl;

  //print plane mode and specify in what it consist
  cout << "plane mode = " << planeMode_ << endl;
  if (planeMode_ == 0) cout << "this means 2D hits are in Induction-1 view" << endl;
  else if (planeMode_ == 1) cout << "this means 2D hits are in Induction-2 view" << endl;
  else if (planeMode_ == 2) cout << "this means 2D hits are in Collection view" << endl;
  cout << " " << endl;

  //return empty recob::MCSFitResult if track length is less than a certain value (ex 150 cm)
  float minlen = 150;
  if (traj.Length() < minlen) {
    cout << "nominal track length less than " << minlen << " cm, stopping fit" << endl;
    return recob::MCSFitResult();
  }
  float norma = sqrt((traj.LocationAtPoint(traj.LastValidPoint()) - traj.LocationAtPoint(traj.FirstValidPoint())).Mag2());
  if (norma < minlen) {
    cout << "effective track length less than " << minlen << " cm, stopping fit" << endl;
    return recob::MCSFitResult();
  };

  //print 3d hits
  cout << "printing 3D hits ..." << endl;
  auto hit3D = traj.FirstValidPoint();
  while (hit3D != recob::TrackTrajectory::InvalidIndex) {
    proxy::TrackPointData pd = pdata[hit3D];
    auto hit = get<1>(pd);
    auto tpc = hit->WireID().TPC; auto cryo = hit->WireID().Cryostat; auto plane = hit->WireID().Plane;
    cout << "TPC = " << tpc << " cryostat = " << cryo << " plane = " << plane << " 3D hit position = " << traj.LocationAtPoint(hit3D) << endl;
    hit3D = traj.NextValidPoint(hit3D + 1);
  }
  cout << " " << endl;

  //print info of first hit (TPC, cryostat, plane)
  auto hit0 = traj.FirstValidPoint();
  proxy::TrackPointData pd0 = pdata[hit0];
  auto info0 = get<1>(pd0);
  auto tpc0 = info0->WireID().TPC; auto cryo0 = info0->WireID().Cryostat; auto plane0 = info0->WireID().Plane;
  cout << "first hit TPC = " << tpc0 << " cryostat = " << cryo0 << " and plane = " << plane0 << endl;

  //print info of last hit (TPC, cryostat, plane)
  auto hit1 = traj.LastValidPoint();
  proxy::TrackPointData pd1 = pdata[hit1];
  auto info1 = get<1>(pd1);
  auto tpc1 = info1->WireID().TPC; auto cryo1 = info1->WireID().Cryostat; auto plane1 = info1->WireID().Plane;
  cout << "last hit TPC = " << tpc1 << " cryostat = " << cryo1 << " and plane = " << plane1 << endl;
  cout << " " << endl;

  //check if track is crossing z = 0
  auto coord0 = traj.LocationAtPoint(hit0); auto coord1 = traj.LocationAtPoint(hit1);
  double z0 = coord0.Z(); double z1 = coord1.Z();
  if (z0 * z1 < 0.) {
    cout << "track is crossing z = 0, stopping fit" << endl;
    return recob::MCSFitResult();
  }

  //print number of hits in the same TPC of first and last hit, for each plane
  vector<int> n_tots; vector<int> n_firsts; vector<int> n_lasts; 
  for (unsigned int p = 0; p <= 2; ++p) {
    int n_tot = 0; int n_first = 0; int n_last = 0;
    cout << "as regards hits belonging to plane " << p << endl;
    auto hit = traj.FirstValidPoint();
    while (hit != recob::TrackTrajectory::InvalidIndex) {
      proxy::TrackPointData pd = pdata[hit];
      auto info = get<1>(pd);
      auto tpc = info->WireID().TPC; auto plane = info->WireID().Plane;
      if (plane == p) {
        ++n_tot;
        if (tpc == tpc0) ++n_first;
        else if (tpc == tpc1) ++n_last;
      }
      hit = traj.NextValidPoint(hit + 1);
    }
    n_tots.push_back(n_tot); n_firsts.push_back(n_first); n_lasts.push_back(n_last);
    cout << "hits in plane " << p << " in same TPC of first hit = " << n_first << " and ratio = " << float(n_first)/float(n_tot) << endl;
    cout << "hits in plane " << p << " in same TPC of last hit = " << n_last << " and ratio = " << float(n_last)/float(n_tot) << endl;
    cout << "total hits in plane " << p << " = " << n_tot << " and total ratio = " << float(n_first+n_last)/float(n_tot) << endl;
  }
  cout << " " << endl;

  //return empty recob::MCSFitResult if there are less than a certain number of hits 2D (ex 30)
  int n_min2D = 30;
  if ((dimMode_ == 2) && (n_tots[planeMode_] < n_min2D)) {
    cout << "number of hits in 2D, plane " << planeMode_ << " less than " << n_min2D << ", stopping fit" << endl;
    return recob::MCSFitResult();
  } 

  //return empty recob::MCSFitResult if there are less than a certain number of hits 3D (ex 50)
  int n_min3D = 50; int n_sum = n_tots[0] + n_tots[1] + n_tots[2];
  if ((dimMode_ == 3) && (n_sum < n_min3D)) {
    cout << "number of hits in 3D less than " << n_min3D << ", stopping fit" << endl;
    return recob::MCSFitResult();
  }

  //define vectors breakpoints with indexes of segmentation, segradlengths with segment lengths in units of X0, cumseglens with cumulative segment lengths
  vector<size_t> breakpoints; vector<float> segradlengths; vector<float> cumseglens;
  //divide trajectory in segments and fill vectors defined above
  breakTrajInSegments(traj, breakpoints, segradlengths, cumseglens, cutMode_, cutLength_);
  //return empty recob::MCSFitResult if there are less than a certain number of segments
  if (segradlengths.size() < minNSegs_) {
    cout << "number of segments less than " << minNSegs_ << ", stopping fit" << endl;
    return recob::MCSFitResult();
  }

  //notify linear fit has started
  cout << "starting computing scattering angles with linear fit..." << endl;
  cout << " " << endl;

  //define vectors dtheta with scattering angles, pcdir with directions of adjacent segments
  vector<float> dtheta_3D; Vector_t pcdir0_3D; Vector_t pcdir1_3D; 
  vector<float> dtheta_2D; Vector_t pcdir0_2D; Vector_t pcdir1_2D;
  dtheta_3D.clear(); dtheta_2D.clear();

  //iterate for every pair of adjacent segments
  for (unsigned int p = 0; p < segradlengths.size(); p++) {
    //print iteration number
    cout << "linear fit: iteration number " << p << endl;

    //perform linear fit of segment and memorize its direction in pcdir1
    linearRegression2D(traj, breakpoints[p], breakpoints[p+1], pcdir1_2D, tpc0, planeMode_);
    linearRegression(traj, breakpoints[p], breakpoints[p + 1], pcdir1_3D);

    //ignore first iteration because it is used to compute direction of first segment
    if (p == 0) {cout << "first iteration, scattering angle impossible to compute" << endl;}

    //start to compute scattering angles from second iteration
    if (p > 0) {
      //ignore current iteration if segment lengths have anomalous values
      if (segradlengths[p] < -100. || segradlengths[p - 1] < -100.) {
        cout << "WARNING! invalid segment length found!" << endl;
        dtheta_3D.push_back(-9.); dtheta_2D.push_back(-9.);
      }
      
      else {
        //print current segment directions
        cout << p << "-th 3D direction, x = " << pcdir0_3D.X() << " y = " << pcdir0_3D.Y() << " z = " << pcdir0_3D.Z() << endl;
        cout << p + 1 << "-th 3D direction, x = " << pcdir1_3D.X() << " y = " << pcdir1_3D.Y() << " z = " << pcdir1_3D.Z() << endl;

        //compute scalar product of two directions, i.e. cosine of scattering angle
	      double cosval = pcdir0_3D.X() * pcdir1_3D.X() + pcdir0_3D.Y() * pcdir1_3D.Y() + pcdir0_3D.Z() * pcdir1_3D.Z();
        if (cosval < -1.0) cosval = -1.0;
        if (cosval > 1.0) cosval = 1.0;
        //compute scattering angle from its cosine, in units of rad
        double dt_3D = acos(cosval);
        //print current scattering angle between the two segments
        cout << p << "-th 3D scattering angle [rad] = " << dt_3D << endl;
  
        //if no cut is performed
        if (cutMode_ == 0) {
          //fill vector dtheta with current scattering angle
          cout << "3D scattering angle considered for the fit" << endl;
          dtheta_3D.push_back(dt_3D);
        }

        //if cut at the end of the track is performed
        else if (cutMode_ == 1) {
          //ignore last 3 scattering angles
          if (p >= segradlengths.size() - 3) {
            cout << "3D scattering angle not considered for the fit" << endl;
            dtheta_3D.push_back(-1.);
          } else {
            //fill vector dtheta with current scattering angle
            cout << "3D scattering angle considered for the fit" << endl;
            dtheta_3D.push_back(dt_3D);
          }
        }
        
        //if cut on the first part of the track is performed
        else if (cutMode_ == 2) {
          //ignores other than first 7 scattering angles
          if (p > 7) {
            cout << "3D scattering angle not considered for the fit" << endl;
            dtheta_3D.push_back(-2.);
          } else {
            //fill vector dtheta with current scattering angle
            cout << "3D scattering angle considered for the fit" << endl;
            dtheta_3D.push_back(dt_3D);
          }
        }

        //print current segment directions
        cout << p << "-th 2D direction, x = " << pcdir0_2D.X() << " y = " << pcdir0_2D.Y() << endl;
        cout << p + 1 << "-th 2D direction, x = " << pcdir1_2D.X() << " y = " << pcdir1_2D.Y() << endl;

        auto pcdir_zero = Vector_t(1, 0, 0);
        if (pcdir0_2D == pcdir_zero || pcdir1_2D == pcdir_zero) {
          cout << "hits belong to wrong TPC! 2D barycenter impossible to compute" << endl;
          dtheta_2D.push_back(-8.);
        } 
        else {
          //compute scalar product of two directions, i.e. cosine of scattering angle
          double cosval_2D = pcdir0_2D.X() * pcdir1_2D.X() + pcdir0_2D.Y() * pcdir1_2D.Y();
          if (cosval_2D < -1.0) cosval_2D = -1.0;
          if (cosval_2D > 1.0) cosval_2D = 1.0;
          cout << "cosval_2D = " << cosval_2D << endl;
          //compute scattering angle from its cosine, in units of rad
          double dt_2D = acos(cosval_2D);
          //print current scattering angle between the two segments
          cout << p << "-th 2D scattering angle [rad] = " << dt_2D << endl;

          //if no cut is performed
          if (cutMode_ == 0) {
            //fill vector dtheta with current scattering angle
            cout << "2D scattering angle considered for the fit" << endl;
            dtheta_2D.push_back(dt_2D);
          }

          //if cut at the end of the track is performed
          else if (cutMode_ == 1) {
            //ignore last 3 scattering angles
            if (p >= segradlengths.size() - 3) {
              cout << "2D scattering angle not considered for the fit" << endl;
              dtheta_2D.push_back(-1.);
            } else {
              //fill vector dtheta with current scattering angle
              cout << "2D scattering angle considered for the fit" << endl;
              dtheta_2D.push_back(dt_2D);
            }
          }
          
          //if cut on the first part of the track is performed
          else if (cutMode_ == 2) {
            //ignores other than first 6 scattering angles
            if (p > 7) {
              cout << "2D scattering angle not considered for the fit" << endl;
              dtheta_2D.push_back(-2.);
            } else {
              //fill vector dtheta with current scattering angle
              cout << "2D scattering angle considered for the fit" << endl;
              dtheta_2D.push_back(dt_2D);
            }
          }
        }
      }
    }
    //update first segment direction for next iteration
    pcdir0_3D = pcdir1_3D; pcdir0_2D = pcdir1_2D;
    cout << " " << endl;
  }

  //define vectors barycenters with barycenter positions, dthetaPoly with scattering angles
  vector<Vector_t> barycenters_3D; Vector_t bary_3D; vector<float> dthetaPoly_3D;
  vector<Vector_t> barycenters_2D; Vector_t bary_2D; vector<float> dthetaPoly_2D;
  dthetaPoly_3D.clear(); dthetaPoly_2D.clear();

  //iterate for every triplet of adjacent segments
  for (unsigned int p = 0; p < segradlengths.size(); p++) {
    //print iteration number
    cout << "poligonal fit: iteration number " << p << endl;

    //find barycenter of segment and memorize its position in bary
    find2DSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary_2D, tpc0, planeMode_);
    barycenters_2D.push_back(bary_2D);
    findSegmentBarycenter(traj, breakpoints[p], breakpoints[p+1], bary_3D);
    barycenters_3D.push_back(bary_3D);

    //ignore first iteration because it is used to compute barycenter of first segment
    if (p == 0) {cout << "first iteration, scattering angle impossible to compute" << endl;}
    //ignore second iteration because it is used to compute barycenter of second segment
    if (p == 1) {cout << "second iteration, scattering angle impossible to compute" << endl;}

    //start to compute scattering angles from third iteration
    if (p > 1) {
      //ignore current iteration if segment lengths have anomalous values
      if (segradlengths[p] < -100. || segradlengths[p - 1] < -100. || segradlengths[p - 2] < -100.) {
        cout << "WARNING! invalid segment length found!" << endl;
        dthetaPoly_3D.push_back(-9.);
        dthetaPoly_2D.push_back(-9.);
      } 

      else {
        //print positions of first and last point of segment
        cout << p - 1 << "-th 3D barycenter = " << barycenters_3D[p - 2] << endl;
        cout << p << "-th 3D barycenter = " << barycenters_3D[p - 1] << endl;
        cout << p + 1 << "-th 3D barycenter = " << barycenters_3D[p] << endl;

        //define delta barycenter plus, direction between first and second barycenter, normalized
        Vector_t dbcp_3D = barycenters_3D[p] - barycenters_3D[p - 1];
        float normp_3D = sqrt(dbcp_3D.X() * dbcp_3D.X() + dbcp_3D.Y() * dbcp_3D.Y() + dbcp_3D.Z() * dbcp_3D.Z());
        dbcp_3D /= normp_3D;
        cout << "direction between " << p - 1 << "-th and " << p << "-th 3D barycenter = " << dbcp_3D << endl;

        //define delta barycenter minus, direction between second and third barycenter, normalized
        Vector_t dbcm_3D = barycenters_3D[p - 1] - barycenters_3D[p - 2];
        float normm_3D = sqrt(dbcm_3D.X() * dbcm_3D.X() + dbcm_3D.Y() * dbcm_3D.Y() + dbcm_3D.Z() * dbcm_3D.Z());
        dbcm_3D /= normm_3D;
        cout << "direction between " << p << "-th and " << p + 1 << "-th 3D barycenter = " << dbcm_3D << endl;

        //compute scalar product of two directions, i.e. cosine of scattering angle
        double cosval = dbcp_3D.X() * dbcm_3D.X() + dbcp_3D.Y() * dbcm_3D.Y() + dbcp_3D.Z() * dbcm_3D.Z();
        if (cosval < -1.0) cosval = -1.0;
        if (cosval > 1.0) cosval = 1.0;
        //compute scattering angle from its cosine, in units of rad
        double dt_3D = acos(cosval);
        //print current scattering angle between the two segments
        cout << p - 1 << "-th 3D scattering angle [rad] = " << dt_3D << endl;

        //if no cut is performed
        if (cutMode_ == 0) {
          //fill vector dtheta with current scattering angle
          cout << "3D scattering angle considered for the fit" << endl;
          dthetaPoly_3D.push_back(dt_3D);
        }

        //if cut at the end of the track is performed
        else if (cutMode_ == 1) {
          //ignore last 3 scattering angles
          if (p >= segradlengths.size() - 3) {
            cout << "3D scattering angle not considered for the fit" << endl;
            dthetaPoly_3D.push_back(-11.);
          } else {
            //fill vector dtheta with current scattering angle
            cout << "3D scattering angle considered for the fit" << endl;
            dthetaPoly_3D.push_back(dt_3D);
          }
        }
        
        //if cut on the first part of the track is performed
        else if (cutMode_ == 2) {
          //ignores other than first 6 scattering angles
          if (p > 7) {
            cout << "3D scattering angle not considered for the fit" << endl;
            dthetaPoly_3D.push_back(-22.);
          } else {
            //fill vector dtheta with current scattering angle
            cout << "3D scattering angle considered for the fit" << endl;
            dthetaPoly_3D.push_back(dt_3D);
          }
        }

        //print positions of first and last point of segment
        cout << p - 1 << "-th 2D barycenter = " << barycenters_2D[p - 2] << endl;
        cout << p << "-th 2D barycenter = " << barycenters_2D[p - 1] << endl;
        cout << p + 1 << "-th 2D barycenter = " << barycenters_2D[p] << endl;

        //check on 2D barycenter
        auto bary_zero = Vector_t(0, 0, 0);
        if (barycenters_2D[p - 2] == bary_zero || barycenters_2D[p - 1] == bary_zero || barycenters_2D[p] == bary_zero) {
          cout << "some hits belong to wrong TPC! 2D barycenter impossible to compute" << endl;
          dthetaPoly_2D.push_back(-88.);
        } 
        else {
          //define delta barycenter plus, direction between first and second barycenter, normalized
          Vector_t dbcp_2D = barycenters_2D[p] - barycenters_2D[p - 1];
          float normp_2D = sqrt(dbcp_2D.X() * dbcp_2D.X() + dbcp_2D.Y() * dbcp_2D.Y());
          dbcp_2D /= normp_2D;
          cout << "direction between " << p - 1 << "-th and " << p << "-th 2D barycenter = " << dbcp_2D << endl;

          //define delta barycenter minus, direction between second and third barycenter, normalized
          Vector_t dbcm_2D = barycenters_2D[p - 1] - barycenters_2D[p - 2];
          float normm_2D = sqrt(dbcm_2D.X() * dbcm_2D.X() + dbcm_2D.Y() * dbcm_2D.Y());
          dbcm_2D /= normm_2D;
          cout << "direction between " << p << "-th and " << p + 1 << "-th 2D barycenter = " << dbcm_2D << endl;

          //compute scalar product of two directions, i.e. cosine of scattering angle
          double cosval_2D = dbcp_2D.X() * dbcm_2D.X() + dbcp_2D.Y() * dbcm_2D.Y();
          if (cosval_2D < -1.0) cosval_2D = -1.0;
          if (cosval_2D > 1.0) cosval_2D = 1.0;
          cout << "cosval_2D = " << cosval_2D << endl;
          //compute scattering angle from its cosine, in units of rad
          double dt_2D = acos(cosval_2D);
          //print current scattering angle between the two segments
          cout << p - 1 << "-th 2D scattering angle [rad] = " << dt_2D << endl;

          //if no cut is performed
          if (cutMode_ == 0) {
            //fill vector dtheta with current scattering angle
            cout << "2D scattering angle considered for the fit" << endl;
            dthetaPoly_2D.push_back(dt_2D);
          }

          //if cut at the end of the track is performed
          else if (cutMode_ == 1) {
            //ignore last 3 scattering angles
            if (p >= segradlengths.size() - 3) {
              cout << "2D scattering angle not considered for the fit" << endl;
              dthetaPoly_2D.push_back(-11.);
            } else {
              //fill vector dtheta with current scattering angle
              cout << "2D scattering angle considered for the fit" << endl;
              dthetaPoly_2D.push_back(dt_2D);
            }
          }
          
          //if cut on the first part of the track is performed
          else if (cutMode_ == 2) {
            //ignores other than first 6 scattering angles
            if (p > 7) {
              cout << "2D scattering angle not considered for the fit" << endl;
              dthetaPoly_2D.push_back(-22.);
            } else {
              //fill vector dtheta with current scattering angle
              cout << "2D scattering angle considered for the fit" << endl;
              dthetaPoly_2D.push_back(dt_2D);
            }
          }
        }
      }
    }
    cout << " " << endl;
  }

  //print vector segradlengths with segment lengths, in units of cm
  cout << "segment lengths [cm] = ";
  for (auto i : segradlengths) {cout << i << ' ';}
  cout << endl;

  //print vector cumseglens with cumulative segment lengths, in units of cm
  cout << "cumulative segment length [cm] = ";
  for (auto i : cumseglens) {cout << i << ' ';}
  cout << endl;

  //print vector dtheta with scattering angles for linear fit, in units of rad
  cout << "3D scattering angles, linear fit [rad] = ";
  for (auto i : dtheta_3D) {cout << i << ' ';}
  cout << endl;

  //print vector dtheta_2D with 2D scattering angles for linear fit, in units of rad
  cout << "2D scattering angles, linear fit [rad] = ";
  for (auto i : dtheta_2D) {cout << i << ' ';}
  cout << endl;

  //print vector dtheta with scattering angles for poligonal fit, in units of rad
  cout << "3D scattering angles, poligonal fit [rad] = ";
  for (auto i : dthetaPoly_3D) {cout << i << ' ';}
  cout << endl;

  //print vector dtheta with scattering angles for poligonal fit, in units of rad
  cout << "2D scattering angles, poligonal fit [rad] = ";
  for (auto i : dthetaPoly_2D) {cout << i << ' ';}
  cout << endl;
  cout << " " << endl;

  //notify performing c2 fit has started
  cout << "starting performing c2 fit..." << endl;
  cout << " " << endl;

  //compute number of hits in 2D, for each plane
  vector<float> weight;
  weight.push_back(float(hits2dI1.size())); 
  weight.push_back(float(hits2dI2.size())); 
  weight.push_back(float(hits2dC.size()));
  //compute delta3p for each plane
  vector<float> delta3ps; 
  delta3ps.push_back(0.1 * d3pI1); 
  delta3ps.push_back(0.1 * d3pI2); 
  delta3ps.push_back(0.1 * d3pC);
  //compute weight average of delta3p using number of hits in each plane as weight
  float delta3pw = inner_product(delta3ps.begin(), delta3ps.end(), weight.begin(), 0.0f);
  delta3pw = delta3pw / accumulate(weight.begin(), weight.end(), 0.0f);
  //return delta3p in function of fhicl parameters, in units of cm
  float delta3p;
  if (dimMode_ == 2) delta3p = delta3ps[planeMode_];
  else if (dimMode_ == 3) delta3p = delta3pw;

  //define vectors dthetaLin, dthetaPoly with scattering angles in linear and poly case
  vector<float> dthetaLin; 
  vector<float> dthetaPoly;
  if (dimMode_ == 2) {dthetaLin = dtheta_2D; dthetaPoly = dthetaPoly_2D;}
  else if (dimMode_ == 3) {dthetaLin = dtheta_3D; dthetaPoly = dthetaPoly_3D;}
  else return recob::MCSFitResult();
  //define vector dtheta and populate it with both linear and polygonal scattering angles
  vector<float> dtheta; 
  dtheta.insert(dtheta.end(), dthetaLin.begin(), dthetaLin.end());
  dtheta.insert(dtheta.end(), dthetaPoly.begin(), dthetaPoly.end());

  //perform a c2fit scan over all the computed angles
  float c2sigma = 0.05; 
  const ScanResult fwdResult = C2Fit(dthetaLin, dthetaPoly, segradlengths, cumseglens, breakpoints, pid, c2sigma, traj);
  return recob::MCSFitResult(pid, fwdResult.p, fwdResult.pUnc, delta3p, 0, 0, 0, segradlengths, dtheta);
}

//break input trajectory into smaller pieces called segments 
void TrajectoryMCSFitterICARUS::breakTrajInSegments(const recob::TrackTrajectory& traj, vector<size_t>& breakpoints, vector<float>& segradlengths, vector<float>& cumseglens, int cutMode, float cutLength) const {
  //define trajectory length
  double trajlen = traj.Length();
  //print trajectory length
  cout << "track length [cm] = " << trajlen << endl;

  //define segment length in function of trajectory length being less or greater than minimum length
  const double thisSegLen = (trajlen > (segLen_ * minNSegs_) ? segLen_ : trajlen / double(minNSegs_));
  //print required segment length, considered segment length, number of segments
  cout << "required segment length [cm] = " << segLen_ << endl;
  cout << "considered segment length [cm] = " << thisSegLen << endl;
  unsigned int nss = trajlen / segLen_;
  cout << "number of segments = " << max(minNSegs_, nss) << endl;

  //define inverse of Argon radiation length X0 = 14 cm
  constexpr double lar_radl_inv = 1./14.0;

  //first segment has zero cumulative length from previous segments
  cumseglens.push_back(0.);
  //initialize current segment length
  double thislen = 0.;

  //find first valid point index and add it to vector breakpoints
  auto nextValid = traj.FirstValidPoint();
  breakpoints.push_back(nextValid);
  //compute position of first valid point
  auto pos0 = traj.LocationAtPoint(nextValid);

  //find next valid point index
  nextValid = traj.NextValidPoint(nextValid + 1);
  //initalize number of valid points in a segment
  int npoints = 0;

  //iterate over all valid points of trajectory
  while (nextValid != recob::TrackTrajectory::InvalidIndex) {
    //compute position of current valid point, then update segment length and position of next valid point
    auto pos1 = traj.LocationAtPoint(nextValid);
    thislen += (pos1 - pos0).R();
    pos0 = pos1;
    //increase number of valid points in a segment
    npoints++;
    
    //break into a new segment if current segment length is greater than required segment length
    if (thislen >= thisSegLen) {
      //add current valid point index to vector breakpoints
      breakpoints.push_back(nextValid);

      //add segment length, in units of radiation length, to vector segradlengths
      if (npoints >= minHitsPerSegment_) segradlengths.push_back(thislen * lar_radl_inv);
      else {
        cout << "WARNING! invalid number of hits per segment!" << endl;
        segradlengths.push_back(-999.);
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
    //add last valid point index to vector breakpoints
    breakpoints.push_back(traj.LastValidPoint() + 1);

    //add segment length, in units of radiation length, to vector segradlengths
    segradlengths.push_back(thislen * lar_radl_inv);

    //add cumulative segment length to vector cumseglens
    cumseglens.push_back(cumseglens.back() + thislen);
  }
  return;
}

//find barycenter of input segment (useful for poligonal fit)
void TrajectoryMCSFitterICARUS::findSegmentBarycenter(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& bary) const {
  //initalize number of points in segment
  int npoints = 0;
  //define an accumulator with position of valid points in the segment
  geo::vect::MiddlePointAccumulator middlePointCalc;
  
  //determine index of first valid point of the segment
  size_t nextValid = firstPoint;
  //iterate until last valid point of the segment
  while (nextValid < lastPoint) {
    //add position of current valid point to vector middlePointCalc
    middlePointCalc.add(traj.LocationAtPoint(nextValid));
    //find next valid point index for next iteration
    nextValid = traj.NextValidPoint(nextValid + 1);
    //increase number of valid points in a segment
    npoints++;
  }

  //determine position of segment barycenter from position of valid points in the segment
  const auto avgpos = middlePointCalc.middlePoint();
  bary = avgpos;
}

//find average direction of trajectory between firstPoint and lastPoint
void TrajectoryMCSFitterICARUS::linearRegression(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir) const {
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
    //find next valid point index for next iteration
    nextValid = traj.NextValidPoint(nextValid + 1);
    //increase number of valid points in the segment
    npoints++;
  }

  //determine position of segment barycenter from position of valid points in the segment
  const auto avgpos = middlePointCalc.middlePoint();
  //define normalization factor as inverse of number of valid points
  const double norm = 1./double(npoints);
  //define covariance matrix m as symmetric matrix 3x3
  TMatrixDSym m(3);

  //determine again index of first valid point
  nextValid = firstPoint;
  //iterate again until last valid point of the segment
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

//find barycenter of input segment (useful for poligonal fit) on some plane
void TrajectoryMCSFitterICARUS::find2DSegmentBarycenter(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& bary2D, double t, unsigned int p) const {
  //initalize number of points in segment
  int npoints = 0;
  //define cumulative sum wsum, ssum over two coordinates w, s
  float sum_x = 0; float sum_y = 0; 
  float kx = 0.622 * 0.1; float ky = 3 * 0.1; 

  //determine index of first valid point of the segment
  size_t nextValid = firstPoint;
  //iterate until last valid point of the segment
  while (nextValid < lastPoint) {
    //get data information of current hit
    proxy::TrackPointData pd = pdata[nextValid];
    auto hit = get<1>(pd);
    unsigned int tpc = hit->WireID().TPC;
    unsigned int plane = hit->WireID().Plane;
    //select only those hits belonging to input plane and tpc of last hit of track
    if (plane == p && tpc == t) {
      //define 2D coordinates new_x (along drift direction) and new_y (rotation of vertical direction)
      auto new_x = hit->PeakTime(); auto new_y = hit->WireID().Wire;
      new_x = kx * new_x; new_y = ky * new_y;
      //update cumulative sums over new_x and new_y coordinates
      sum_x += new_x; sum_y += new_y;
      //increase number of valid points in a segment
      npoints++;
    }
    //find next valid point index for next iteration
    nextValid = traj.NextValidPoint(nextValid + 1);
  }

  //return barycenter coordinates if there are hits, otherwise return null
  if (npoints > 0) {
    //compute barycenter coordinates as average position over xprime and yprime coordinates
    const auto avg_x = float(sum_x / npoints);
    const auto avg_y = float(sum_y / npoints);
    bary2D.SetXYZ(avg_x, avg_y, 0);
  } else {
    bary2D.SetXYZ(0, 0, 0);
  }
}

//find average direction of trajectory between firstPoint and lastPoint on some plane
void TrajectoryMCSFitterICARUS::linearRegression2D(const recob::TrackTrajectory& traj, const size_t firstPoint, const size_t lastPoint, Vector_t& pcdir2D, double t, unsigned int p) const {
  //initalize number of points in segment
  int npoints = 0;
  //define cumulative sum wsum, ssum over two coordinates w, s
  float sum_x = 0; float sum_y = 0; 
  float kx = 0.622 * 0.1; float ky = 3 * 0.1; 

  //determine index of first valid point of the segment
  size_t nextValid = firstPoint;
  //iterate until last valid point of the segment
  while (nextValid < lastPoint) {
    //get data information of current hit
    proxy::TrackPointData pd = pdata[nextValid];
    auto hit = get<1>(pd);
    unsigned int tpc = hit->WireID().TPC;
    unsigned int plane = hit->WireID().Plane;
    //select only those hits belonging to input plane and tpc of last hit of track
    if (plane == p && tpc == t) {
      //define 2D coordinates new_x (along drift direction) and new_y (rotation of vertical direction)
      auto new_x = hit->PeakTime(); auto new_y = hit->WireID().Wire;
      new_x = kx * new_x; new_y = ky * new_y;
      //update cumulative sums over new_x and new_y coordinates
      sum_x += new_x; sum_y += new_y;
      //increase number of valid points in a segment
      npoints++;
    }
    //find next valid point index for next iteration
    nextValid = traj.NextValidPoint(nextValid + 1);
  }

  //compute average position over xprime and yprime coordinates
  const auto avg_x = float(sum_x / npoints);
  const auto avg_y = float(sum_y / npoints);
  //define normalization factor as inverse of number of valid points
  const double norm = 1./double(npoints);
  //define covariance matrix m as symmetric matrix 2x2
  TMatrixDSym m(2);

  //determine again index of first valid point
  nextValid = firstPoint;
  //iterate again until last valid point of the segment
  while (nextValid < lastPoint) {
    //get data information of current hit
    proxy::TrackPointData pd = pdata[nextValid];
    auto hit = get<1>(pd);
    unsigned int tpc = hit->WireID().TPC;
    unsigned int plane = hit->WireID().Plane;

    //select only those hits belonging to input plane and tpc of last hit of track
    if (plane == p && tpc == t) {
      //define 2D coordinates new_x (along drift direction) and new_y (rotation of vertical direction)
      auto new_x = kx * hit->PeakTime();
      auto new_y = ky * hit->WireID().Wire;
      //compute coordinate differences between current valid point and average point
      const double xxw0 = new_x - avg_x;
      const double yyw0 = new_y - avg_y;
      //update covariance matrix values with normalized values of coordinate differences above
      m(0, 0) += xxw0 * xxw0 * norm;
      m(0, 1) += xxw0 * yyw0 * norm;
      m(1, 0) += yyw0 * xxw0 * norm;
      m(1, 1) += yyw0 * yyw0 * norm;
    }
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
  for (int i = 1; i < 2; ++i) {
    //update maximum eigenvalue if current is greater
    if (eigenval(i) > maxeval) {
      maxevalidx = i;
      maxeval = eigenval(i);
    }
  } 
  
  //update segment direction with eigenvector associated to maximum eigenvalue
  pcdir2D = Vector_t(eigenvec(0, maxevalidx), eigenvec(1, maxevalidx), 0);
}

//find most probable value of energy loss according to Landau distribution
double TrajectoryMCSFitterICARUS::energyLossLandau(const double mass2, const double e2, const double x) const {
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

//find mean energy loss according to BetheBloch distribution
double TrajectoryMCSFitterICARUS::energyLossBetheBloch(const double mass,const double e2) const {
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

//find energy of track at a certain distance travelled
double TrajectoryMCSFitterICARUS::GetE(const double initial_E, const double length_travelled, const double m) const {
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
    //compute energy according to Bethe Bloch, negleting density correction, in units of GeV
    if (eLossMode_ == 2) {
      double dedx = energyLossBetheBloch(m, current_E * current_E);
      current_E -= (dedx * step_size);
    }

    //compute energy according to Landau distribution, in units of GeV
    else {
      current_E -= energyLossLandau(m2, current_E * current_E, step_size);
    }

    //return null energy if current energy is less or equal than energy at rest
    if (current_E <= m) {return 0.;}
  }

  //return current energy of particle in function of travelled distance, in units of GeV
  return current_E;
}

//function not used yet
double TrajectoryMCSFitterICARUS::GetOptimalSegLen(const recob::TrackTrajectory& tr, const double guess_p, const int n_points, const int plane, const double length_travelled) const {
  //define energy loss for MIP muon in liquid Argon, in units of MeV mm-1
  constexpr double kcal = 0.2105;
  //define momentum as difference between guess and energy loss times travelled distance, in units of MeV
  double initial_p = guess_p - kcal * length_travelled;
  double MIN_P = 1000;
  if(initial_p < MIN_P) initial_p = MIN_P;

  //define measurement error as delta3p
  //double sigma = d3pC;

  //define two angles a by cos(a) and b by sin(b)
  //double cosa = 1;
  //double sinb = 1; 

  //define total distance travelled (check if it is equal to length_travelled)
  //double DsTot;
  //if (cutMode() == 0) DsTot = sqrt((tr.LocationAtPoint(tr.LastValidPoint()) - tr.LocationAtPoint(0)).Mag2());
  //if (cutMode() == 1) DsTot = sqrt((tr.LocationAtPoint(tr.LastValidPoint()) - tr.LocationAtPoint(0)).Mag2()) - cutLength();
  //if (cutMode() == 2) DsTot = cutLength();

  //define start point, end point, average direction and projection to collection plane
  TVector3 start = tr.LocationAtPoint<TVector3>(0);
  TVector3 end = tr.LocationAtPoint<TVector3>(tr.LastValidPoint());
  TVector3 avdir = end - start;
  TVector3 projColl(0, 0, 0);
  projColl[0] = 1; projColl[1] = sqrt(3.) / 2.; projColl[2] = 0.5;
  //double cos = avdir * projColl;

  //define Highland parameters S2 in units of MeV, epsilon, X0 in units of mm
  //double S2 = 13.6; double epsilon = 0.038; double X0 = 140.;

  //lots of calculation which cannot be understood (Filippo aiutami a capirle)
  //double xi = (sigma / DsTot) * (initial_p / 2.) / S2 * sqrt(X0 / length_travelled) / sqrt(double(n_points)) * sinb;
  //double K = 6 * xi * xi;
  //double m = 1 / sqrt(6 * xi);
  //double sigmaMS = S2 / (initial_p / 2.) * sqrt(length_travelled / X0 / m / cosa) / sqrt(2.) / cosa;
  //double sigmaMeas = pow(m, 1.5) * sqrt(double(6.)) * sigma * sinb / cosa / DsTot / sqrt(double(n_points));
  //double alfa = 1 + epsilon * log(DsTot / sinb / double(m) / X0);
  //double m_corr = (m == 1)? m : sqrt(m / (m - 1));
  //double xiCorr = sigma / DsTot * ((initial_p / 2.) / S2) * sqrt(X0 / length_travelled) * sqrt((double)n_points) * sinb/ alfa;
  //double m2 = m_corr / sqrt(6. * xiCorr);
  //if (m == 1) m2 = 10;
  //m2 += 0.5;
  return 14.;
}

//compute delta3p used for measurement error along drift direction
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

  //return d3p = 0.4 as default value if there are no valid residuals
  if (!h0.size()) d3p = 0.4;
  //otherwise return d3p as rms of residuals saved in histogram 
  else d3p = hd3pv->GetRMS();
  
  //fill d3p with right quantity depending on the plane
  if (plane == 0) d3pI1 = d3p;
  if (plane == 1) d3pI2 = d3p;
  if (plane == 2) d3pC = d3p;
}

//compute delta3p in 3D used for measurement error along drift direction
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

  //return d3p = 0.4 as default value if there are no valid residuals
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

//compute 3D residual used in delta3p computation
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
const double TrajectoryMCSFitterICARUS::C2Function(const recob::TrackTrajectory& tr, vector<float> cumseglens, vector<long unsigned int> breakpoints, vector<float> dthetaLin, vector<float> dthetaPoly, double p0, const double m) const {
  //define thetaexp, thetams, thetaerr that represent expected, measured and error angle
  double thetaexp, thetams, thetaerr;

  //define ttallLin, ttallPoly 
  vector<float> ttall; vector<float> ttallLin; vector<float> ttallPoly;
  ttall.clear(); ttallLin.clear(); ttallPoly.clear();

  //define first segment index and check if it corresponds to positive angle
  unsigned int firstseg = 0; bool check = false;
  for (unsigned int js = 0; js < dthetaLin.size(); js++) {
    if (dthetaLin[js] > 0. && dthetaPoly[js] > 0.) {
      cout << "first linear angle = " << dthetaLin[js] << " first polygonal angle = " << dthetaPoly[js] << " first index = " << js << endl;
      firstseg = js; check = true;
      break;
    }
  }
  if (!check) return 0;

  //define last segment index and check if it corresponds to positive angle
  unsigned int lastseg = cumseglens.size() - 2;
  for (unsigned int js = cumseglens.size() - 3; js >= 0; js--) {
    if (dthetaLin[js] > 0.) {
      cout << "last linear angle = " << dthetaLin[js] << " last polygonal angle = " << dthetaPoly[js-1] << " last index = " << js << endl;
      lastseg = js+1;
      break;
    }
  }

  //compute number of segments as difference (+1) between first and last segment indexes
  unsigned int nseg = lastseg - firstseg + 1;
  unsigned int nsegtot = cumseglens.size() - 1;
  //check if number of segments is greater than a certain value (ex 3) if not return null
  if (nseg < 7) return 0;

  //compute segment length for current segment, in units of cm
  double Ls = tr.Length() / double(nsegtot); Ls = 14.;

  //define Highland parameters S2 in units of MeV, epsilon, X0 in units of cm
  double S2 = 13.6; double epsilon = 0.038; double X0 = 14.;
  cout << "Highland parameter S2 [MeV] = " << S2 << endl;
  cout << "Highland parameter epsilon = " << epsilon << endl;
  cout << "Argon radiation length X0 [cm] = " << X0 << endl;
  cout << "used segment length [cm] = " << Ls << endl;
  
  //define alfa that represent 1 + epsilon * log (X / X0) term in Highland formula
  double alfa = 1 + epsilon * log(Ls / X0);
  cout << "bracket term in Highland formula 1 + epsilon * log (X / X0) = " << alfa << endl;
  cout << " " << endl;

  //define covariance matrix mat, sum of MCS and measurement
  TMatrixDSym mat(2 * nseg - 3);
  TMatrixDSym matFull(2 * nseg - 3);

  //define covariance matrix matLin (linear part), sum of MCS and measurement
  TMatrixDSym matLin(nseg - 1);
  //define covariance matrix matMCSLin, with MCS of linear part
  TMatrixDSym matMCSLin(nseg - 1);
  //define covariance matrix matErrLin, with measurement of linear part
  TMatrixDSym matErrLin(nseg - 1);

  //define covariance matrix matPoly (polygonal part), sum of MCS and measurement
  TMatrixDSym matPoly(nseg - 2);
  //define covariance matrix matMCSPoly, with MCS of polygonal part
  TMatrixDSym matMCSPoly(nseg - 2);
  //define covariance matrix matErrPoly, with measurement of polygonal part
  TMatrixDSym matErrPoly(nseg - 2);

  //define matrix matmix with mix terms
  TMatrixDSym matmix(2 * nseg - 3);
 
  //define number of points in the track
  int np = tr.NPoints();

  //compute ratio between 2D points and 3D points of a track, for each plane
  vector<float> ratios2D3D; 
  ratios2D3D.push_back(float(hits2dI1.size()) / float(np));
  ratios2D3D.push_back(float(hits2dI2.size()) / float(np));
  ratios2D3D.push_back(float(hits2dC.size()) / float(np));

  //define average number of 3D points in a segment
  float avg3Dperseg = float(np)/float(nsegtot);
  cout << "average number of 3D points in a segment = " << avg3Dperseg << endl;

  //define average number of 2D points in a segment
  float avg2Dperseg = ratios2D3D[planeMode_] * avg3Dperseg;
  cout << "average number of 2D points in a segment = " << avg2Dperseg << endl;
  cout << " " << endl;

  //compute measurement error delta3p, in units of mm
  vector<float> d3ps; double sigma0;
  d3ps.push_back(d3pI1); d3ps.push_back(d3pI2); d3ps.push_back(d3pC); 
  cout << "delta3p in induction1 view [mm] = " << d3ps[0] << " weight = " << ratios2D3D[0] << endl;
  cout << "delta3p in induction2 view [mm] = " << d3ps[1] << " weight = " << ratios2D3D[1] << endl;
  cout << "delta3p in collection view [mm] = " << d3ps[2] << " weight = " << ratios2D3D[2] << endl;
  float d3p = inner_product(d3ps.begin(), d3ps.end(), ratios2D3D.begin(), 0.0f);
  d3p = d3p / accumulate(ratios2D3D.begin(), ratios2D3D.end(), 0.0f);
  cout << "delta3p weight average in all views [mm] = " << d3p << endl;
  cout << " " << endl;

  //define sigma0 from delta3p, in units of cm
  if (dimMode_ == 2) sigma0 = 0.1 * d3ps[planeMode_];
  else if (dimMode_ == 3) sigma0 = 0.1 * d3p;
  else sigma0 = 0.1;

  //notify computation of scattering angles for linear fit has started
  cout << "starting computation of scattering angles for linear fit ..." << endl;
  cout << " " << endl;

  //convert input p0 with units of MeV/c to p1 with units of GeV/c
  const double p1 = 0.001 * p0;
  cout << "momentum [GeV/c] = " << p1 << endl;
  cout << "mass [GeV/c2] = " << m << endl;
  //compute total energy from momentum, in units of GeV
  const double p2 = p1 * p1;
  const double m2 = m * m;
  const double Etot = sqrt(p2 + m2);
  cout << "energy [GeV] = " << Etot << endl;
  cout << " " << endl;

  for (unsigned int jp = firstseg; jp < lastseg; jp++) {
    //print iteration number
    cout << "linear fit: segment number " << jp << endl;
    
    //compute current energy at current length, in units of MeV
    const double Eij = GetE(Etot, cumseglens[jp], m);
    cout << "current energy [GeV] = " << Eij << endl;
    const double Eij2 = Eij * Eij;
    //check if energy at current segment is greater than energy at rest, if not return null
    if (Eij2 <= m2) {
      cout << "energy at current segment is less than energy at rest" << endl;
      return 0;
    }
    
    //compute momentum and velocity of particle at current segment
    const double pij = sqrt(Eij2 - m2);
    const double beta = pij / Eij;
    cout << "current momentum [GeV/c] = " << pij << " and current beta = " << beta << endl;

    //update ttall vector with measured scattering angle (linear fit)
    thetaexp = 1000. * dthetaLin[jp];
    ttallLin.push_back(thetaexp);
    cout << "current theta MCS measured [mrad] = " << thetaexp << endl;

    //compute theta MCS just as claimed by Highland formula
    if (dimMode_ == 2) thetams = (1 / sqrt(2)) * S2 / (pij * beta) * sqrt(Ls / X0) * alfa;
    if (dimMode_ == 3) thetams = S2 / (pij * beta) * sqrt(Ls / X0) * alfa;
    cout << "current theta MCS expected from Highland [mrad] = " << thetams << endl;

    //define measurement error associated to delta3p (2D and 3D)
    if (dimMode_ == 2) thetaerr = 1000. * sigma0 * sqrt(24.) / Ls / sqrt(float(avg2Dperseg));
    if (dimMode_ == 3) thetaerr = sqrt(2) * 1000. * sigma0 * sqrt(24.) / Ls / sqrt(float(avg3Dperseg));
    cout << "current theta MCS error from delta3p [mrad] = " << abs(thetaerr) << endl;

    //update matrix matMCSLin, matErrLin with MCS and measurement error squared angle
    matMCSLin(jp-firstseg, jp-firstseg) = thetams * thetams;
    matErrLin(jp-firstseg, jp-firstseg) = thetaerr * thetaerr;
  
    //FillCovMatrixSegOnly(tr, mat, jp, thetams * thetams, thetaerr * thetaerr, materr, breakpoints);
    //AddSegmentCovariance(tr, mat, jp);
    cout << " " << endl;
  }

  for (unsigned int jp = firstseg; jp < lastseg - 1; jp++) {
    //print iteration number
    cout << "poligonal fit: segment number " << jp << endl;
    
    //compute current energy at current length, in units of MeV
    const double Eij = GetE(Etot, cumseglens[jp], m);
    cout << "current energy [GeV] = " << Eij << endl;
    const double Eij2 = Eij * Eij;
    //check if energy at current segment is greater than energy at rest, if not return null
    if (Eij2 <= m2) {
      cout << "energy at current segment is less than energy at rest" << endl;
      return 0;
    }
    
    //compute momentum and velocity of particle at current segment
    const double pij = sqrt(Eij2 - m2);
    const double beta = pij / Eij;
    cout << "current momentum [GeV/c] = " << pij << " and current beta = " << beta << endl;

    //update ttall vector with measured scattering angle (linear fit)
    thetaexp = 1000. * dthetaPoly[jp];
    ttallPoly.push_back(thetaexp);
    cout << "current theta MCS measured [mrad] = " << thetaexp << endl;

    //compute theta MCS just as claimed by Highland formula
    if (dimMode_ == 2) thetams = (1 / sqrt(2)) * S2 / (pij * beta) * sqrt(Ls / X0) * alfa;
    if (dimMode_ == 3) thetams = S2 / (pij * beta) * sqrt(Ls / X0) * alfa;
    cout << "current theta MCS expected from Highland [mrad] = " << thetams << endl;

    //define measurement error associated to delta3p
    if (dimMode_ == 2) thetaerr = 1000. * sigma0 * sqrt(6.) / Ls / sqrt(float(avg2Dperseg));
    if (dimMode_ == 3) thetaerr = sqrt(2) * 1000. * sigma0 * sqrt(6.) / Ls / sqrt(float(avg3Dperseg));
    cout << "current theta MCS error from delta3p [mrad] = " << abs(thetaerr) << endl;

    //update matrix matMCSPoly, matErrPoly with MCS and measurement error squared angle
    matMCSPoly(jp-firstseg, jp-firstseg) = thetams * thetams;
    matErrPoly(jp-firstseg, jp-firstseg) = thetaerr * thetaerr;
  
    //FillCovMatrix(tr, matpoly, jp, thetams * thetams, thetaerr * thetaerr, matpolyms, matpolyerr, breakpoints, firstseg);
    //FillCovMixTerms(tr, matmix, jp, nseg - 2, thetams * thetams, thetaerr * thetaerr);
    //AddSegmentCovariance(tr, matpoly, jp);
    cout << " " << endl;
  }

  //sum MCS and error matrices, linear fit
  for (int jm = 0; jm < matMCSLin.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matMCSLin.GetNrows(); jmm++) {  
      matLin(jm, jmm) = matMCSLin(jm, jmm) + matErrLin(jm, jmm); } }
  cout << "matrix matLin, sum of MCS and error terms for linear fit = " << endl; matLin.Print();
  cout << "dimension of matrix matLin = " << matLin.GetNrows() << endl;

  //sum MCS and error matrices, polygonal fit
  for (int jm = 0; jm < matMCSPoly.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matMCSPoly.GetNrows(); jmm++) {
      matPoly(jm, jmm) = matMCSPoly(jm, jmm) + matErrPoly(jm, jmm); } }
  cout << "matrix matPoly, sum of MCS and error terms for polygonal fit = " << endl; matPoly.Print();
  cout << "dimension of matrix matPoly = " << matPoly.GetNrows() << endl;

  //merge linear and polygonal matrices matLin, matPoly into matFull
  for (int jm = 0; jm < matLin.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matLin.GetNrows(); jmm++) {
      matFull(jm, jmm) = matLin(jm, jmm); } }
  for (int jm = 0; jm < matPoly.GetNrows(); jm++) {
    for (int jmm = 0; jmm < matPoly.GetNrows(); jmm++) {
      matFull(jm + matLin.GetNrows(), jmm + matLin.GetNrows()) = matPoly(jm, jmm); } }
  cout << "matrix matFull, combination of matLin and matPoly = " << endl; matFull.Print();
  cout << "dimension of matrix matFull = " << matFull.GetNrows() << endl;

  //in case of both linear and poly matrices, take matFull
  if (fitMode_ == 0) {
    mat = matFull;
    ttall.insert(ttall.end(), ttallLin.begin(), ttallLin.end());
    ttall.insert(ttall.end(), ttallPoly.begin(), ttallPoly.end()); }
  //in case of only linear matrix, take matLin
  if (fitMode_ == 1) {
    mat.ResizeTo(matMCSLin.GetNrows(), matMCSLin.GetNrows());
    mat = matLin; 
    ttall = ttallLin; }
  //in case of only polygonal matrix, take matPoly
  if (fitMode_ == 2) {
    mat.ResizeTo(matMCSPoly.GetNrows(), matMCSPoly.GetNrows());
    mat = matPoly;
    ttall = ttallPoly; }

  //print matrix mat which depends on fitMode chosen 
  cout << "fit mode = " << fitMode_ << endl;
  cout << "matrix mat = " << endl; mat.Print();
  cout << "dimension of matrix mat = " << mat.GetNrows() << endl;
  
  //initialize matrix cov as a copy of matrix mat
  TMatrixD cov = mat;
  cout << "matrix cov, copy of matrix mat = " << endl;
  cov.Print();

  //define column matrix vtall with same size of vector ttall, and fill it with elements of vector ttall
  TMatrixD vtall(ttall.size(), 1);
  for (unsigned int jv = 0; jv < ttall.size(); jv++) vtall(jv, 0) = ttall[jv]; 
  cout << "column matrix vtall, with elements of vector ttall = " << endl;
  vtall.Print();

  //define row matrix tvtall as the transpose of column matrix vtall
  TMatrixD tvtall(TMatrixD::kTransposed, vtall); 
  cout << "row matrix tvtall, transpose of column matrix vtall = " << endl;
  tvtall.Print();

  //define matrix invcov as the inverse of matrix cov
  TMatrixD invcov = cov.Invert();
  cout << "matrix invcov, inverse of matrix cov = " << endl;
  invcov.Print();

  //define column matrix vtcov as matrix product of invcov and vtall
  TMatrixD vtcov = invcov * vtall;
  cout << "matrix vtcov, product of matrix invcov and column matrix vtall = " << endl;
  vtcov.Print();

  //initialize vtcovmed and vector terms
  double vtcovmed = 0;
  vector<double> terms;
  //iterate over elements of row matrix tvtall and column matrix vtcov
  for (unsigned int jv = 0; jv < ttall.size(); jv++) {
    double term = tvtall(0, jv) * vtcov(jv, 0);
    //add jv-th product between element of tvtall and element of vtcov to terms
    terms.push_back(term);
    //update cumulative sum of products between elements of tvtall and elements of vtcov
    vtcovmed += term;
  }
  TMatrixD mterms(terms.size(), 1);
  for (unsigned int jv = 0; jv < terms.size(); jv++) mterms(jv, 0) = terms[jv]; 
  cout << "column matrix mterms, product of matrix tvtall and matrix vtcov = " << endl;
  mterms.Print();
  //divide vtcovmed over number of ttall elements
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
  for (unsigned int jt = 0; jt < ttall.size(); jt++) {
    //check if element jt-th of vector terms is under threshold, if not add to vector tails
    if(terms[jt] < vtcovmed * 9.) {
      ttrunctot += terms[jt];
      atrunc.push_back(ttall[jt]);
      ttrunc.push_back(terms[jt]);
    } else {
      tails.push_back(jt);
    }
  }
  
  //define column matrix vtrunc and fill it with ttall values of accepted terms, then print it
  TMatrixD vtrunc(atrunc.size(), 1);
  for (unsigned int jv = 0; jv < atrunc.size(); jv++) vtrunc(jv, 0) = atrunc[jv];
  cout << "column matrix vtrunc, with accepted elements of vector ttall = " << endl;
  vtrunc.Print();

  //define row matrix tvtrunc as the transpose of column matrix vtrunc
  TMatrixD tvtrunc(TMatrixD::kTransposed, vtrunc); 
  cout << "row matrix tvtrunc, transpose of column matrix vtrunc = " << endl;
  tvtrunc.Print();

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
      covs.push_back(covtemp);
    }
  }

  //define matrix covcut which is final covariance matrix after cut
  TMatrix covcut = covs[covs.size() - 1];
  cout << "matrix covcut, final covariance matrix after cut = " << endl;
  covcut.Print();
  //check if this matrix is not empty and eventually return null
  if (!covcut.GetNrows()) return -999;

  //define matrix invcovmod as the copy of matrix covcut
  TMatrixD invcovmod = covcut;
  cout << "matrix invcovmod, copy of matrix covcut = " << endl;
  invcovmod.Print();

  //define matrix vtcovmod as product of matrix invcovmod and column matrix vtrunc, then print it
  TMatrixD vtcovmod = invcovmod * vtrunc;
  cout << "matrix vtcovmod, product of matrix invcovmod and column matrix vtrunc = " << endl;
  vtcovmod.Print();

  //initialize vmediomod and vector tterms
  double vmediomod = 0;
  vector<double> tterms;
  //iterate over elements of row matrix tvtrunc and column matrix vtcovmod
  for (unsigned int jv = 0; jv < ttrunc.size(); jv++) {
    double tterm = tvtrunc(0, jv) * vtcovmod(jv, 0);
    //add jv-th product between element of tvtrunc and element of vtcovmod to terms
    tterms.push_back(tterm);
    //update cumulative sum of products between elements of tvtrunc and elements of vtcovmod
    vmediomod += tterm;
  }
  TMatrixD mtterms(tterms.size(), 1);
  for (unsigned int jv = 0; jv < tterms.size(); jv++) mtterms(jv, 0) = tterms[jv]; 
  cout << "column matrix mtterms, product of matrix tvtrunc and matrix vtcovmod = " << endl;
  mtterms.Print();
  //divide vtmediomod over number of ttrunc elements, then return it
  vmediomod /= (ttrunc.size());
  cout << "correctly computed c2 function " << endl;
  return vmediomod;
}

//function not used yet
const void TrajectoryMCSFitterICARUS::FillCovMatrixSegOnly(recob::TrackTrajectory tr, TMatrixDSym mat, unsigned int jp, double sms, double serr, TMatrixDSym materr, vector<long unsigned int> breaks) const {
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

//add segment covariance to input matrix
const void TrajectoryMCSFitterICARUS::AddSegmentCovariance(recob::TrackTrajectory tr, TMatrixDSym mat, int jp) const {
  //define index previous to input index
  int ip = jp - 1;
  //define dimension of input matrix
  int np = mat.GetNrows();
  //define matrix matbd with same dimension of input matrix
  TMatrixDSym matbd(np);

  //define start point, end point, average direction and projection to collection plane
  TVector3 start = tr.LocationAtPoint<TVector3>(0);
  TVector3 end = tr.LocationAtPoint<TVector3>(tr.LastValidPoint());
  TVector3 avdir = end - start;
  TVector3 projColl(0, 0, 0);
  projColl[0] = 1; projColl[1] = sqrt(3.)/2.; projColl[2] = 0.5;

  //initialize some variables for propagation error
  double t0m = 0; double t0p = 0;
  double tmm = 0; double t00 = 0; double tpp = 0;
  //tmm = point(jm - 1).ESeg(); t00 = point(jm).ESeg(); tpp = point(jm + 1).ESeg();
  double sm = 1; double s0 = 1; double sp = 1;
  //sm = point(jm - 1).Sin2B(); s0 = point(jm).Sin2B(); sp = point(jm + 1).Sin2B();
  double tppp = 0;
  double spp = 1;

  //define segment before and after input point (segmentation point) 
  auto pos0 = tr.LocationAtPoint(jp - 1);
  auto pos1 = tr.LocationAtPoint(jp);
  auto pos2 = tr.LocationAtPoint(jp + 1);
  double dxm = (pos1 - pos0).R();
  double dxp = (pos2 - pos1).R();
  //define also segment after segment after input point 
  double dxpp = -1;
  if (jp < mat.GetNrows()) {
    auto pos3 = tr.LocationAtPoint(jp + 2);
    dxpp = (pos3 - pos2).R();
  }

  //compute contribution terms to covariance matrix
  double t1 = t00 * t00 * s0 * s0 * (1 / dxp + 1 / dxm) * (1 / dxp + 1 / dxm);
  double t2 = tpp * tpp * sp * sp / dxp / dxp;
  double t3 = tmm * tmm * sm * sm / dxm / dxm;
  matbd(ip, ip) = (t1 + t2 + t3);
  matbd(ip, ip) += (-2 * t0p * t0p * s0 * sp / dxp * (1 / dxp + 1 / dxm) - 2 * t0m * t0m * s0 * sm / dxm * (1 / dxp + 1 / dxm));

  //add contribution to non-diagonal terms of covariance matrix
  if (jp < np){
    //spp = point(jm + 2).Sin2B();
    matbd(ip, ip + 1) = - t00 * t00 * s0 * s0 / dxp * (1 / dxm + 1 / dxp) - tpp * tpp * sp * sp / dxp * (1 / dxpp + 1 / dxp);
    matbd(ip, ip + 1) += (tppp * tppp * sp * s0 / dxpp / dxp + t0m * t0m * s0 * sp / dxp / dxm + t0p * t0p * s0 * sp * (1 / dxp / dxp + (1 / dxp + 1 / dxm) * (1 / dxp + 1 / dxpp)));
  }
  if (jp < np - 1) {
    matbd(ip, ip + 2) = tpp * tpp * sp * sp / dxp / dxpp;
    //spp = point(jp + 2).Sin2B();
    auto pos3 = tr.LocationAtPoint(jp + 2);
    auto pos4 = tr.LocationAtPoint(jp + 3);
    double dxppp = (pos4 - pos3).R();
    matbd(ip, ip + 2) += - t0p * t0p * s0 * spp / dxpp * (1 / dxp + 1 / dxm) - tppp * tppp * s0 * spp / dxp * (1 / dxppp + 1 / dxpp);
  }
  if (jp < np - 2) {
    //spp = point(jm + 2).Sin2B();
    matbd(ip, ip + 2) = tppp * tppp * sp * spp / dxp / dxpp;
  }

  //guarantee matrix matbd stay symmetric
  if (jp < np) matbd(ip + 1, ip) = matbd(ip, ip + 1);
  if (jp < np - 1) matbd(ip + 2, ip) = matbd(ip, ip + 2);
  if (jp < np - 2) matbd(ip + 3, ip) = matbd(ip, ip + 3);

  //return sum of input matrix mat and new matrix matbd
  mat += matbd;
}

//fill covariance matrix
const void TrajectoryMCSFitterICARUS::FillCovMatrix(recob::TrackTrajectory tr, TMatrixDSym mat, int jp, double sms, double serr, TMatrixDSym matms, TMatrixDSym materr, vector<long unsigned int> breakpoints, int firstseg) const {
  //define index previous to firstseg
  int ip = jp - firstseg;
  //define dimension of input matrix
  int matsize = mat.GetNrows();
  //double dxmedio = 76.7179;

  //define muon mass in MeV/c2
  double MMU = 105.6;

  //define segment before and after input point (segmentation point) 
  auto pos0 = tr.LocationAtPoint(jp - 1);
  auto pos1 = tr.LocationAtPoint(jp);
  auto pos2 = tr.LocationAtPoint(jp + 1);
  double dxm = (pos1 - pos0).R();
  double dxp = (pos2 - pos1).R();

  //define momentum at input point, before input point, after input point
  //double pm = tr.MomentumAtPoint(jp - 1);
  double p0 = tr.MomentumAtPoint(jp);
  double pp = tr.MomentumAtPoint(jp + 1);

  //define beta=v/c at input point, before input point, after input point
  //double bm = pm / sqrt(pm * pm + MMU * MMU);
  double b0 = p0 / sqrt(p0 * p0 + MMU * MMU);
  double bp = pp / sqrt(pp * pp + MMU * MMU);

  //define index difference at input point, before input point, after input point
  double em, ep, e0;
  em = breakpoints[jp - 1] - breakpoints[jp - 2];
  e0 = breakpoints[jp] - breakpoints[jp - 1];
  ep = breakpoints[jp + 1] - breakpoints[jp];

  //define angular coefficient at input point, before input point, after input point
  double sm = 1; double s0 = 1; double sp = 1;
  //sm = point(jp - 1).Sin2B();
  //s0 = point(jp).Sin2B();
  //sp = point(jp + 1).Sin2B();

  //define covariance contribution due to MCS and error
  double c0;
  //c0 = sms / p0 / p0 / b0 / b0 * (dxp + dxm) / 2 + 2 * serr * s0 * s0 / e0 / e0 * (1 / dxp / dxp + 1 / dxm / dxm + 1 / dxp / dxm);
  c0 = sms / p0 / p0 / b0 / b0 * (dxp + dxm) / 2 + serr * (s0 * s0 / e0 / e0 * (1 / dxp + 1 / dxm) * (1 / dxp + 1 / dxm) + sp * sp / ep / ep / dxp / dxp + sm * sm / em / em / dxm / dxm);
  
  //update input matrix mat, matms, materr
  mat(0, 0) = c0;
  matms(ip, ip) = sms / p0 / p0 / b0 / b0 * (dxp + dxm) / 2;
  materr(ip, ip) = serr * (s0 * s0 / e0 / e0 * (1 / dxp + 1 / dxm) * (1 / dxp + 1 / dxm) + sp * sp / ep / ep / dxp / dxp + sm * sm / em / em / dxm / dxm);

  //define covariance contribution non-diagonal between point jp and jp+1
  if (jp < matsize - 1) {
    //dxp = point(ip + 1).DS();
    auto pos3 = tr.LocationAtPoint(jp + 2);
    double dxpp = (pos3 - pos2).R();
    //dxpp = dxmedio;

    //double p0p = sqrt(p0 * pp)
    double smsp = sms / p0 / pp / b0 / bp;
    double errp = serr * (s0 * s0 / e0 / e0 / dxp * (1 / dxm + 1 / dxp) + sp * sp / ep / ep / dxp * (1 / dxpp + 1 / dxp));

    mat(jp, jp + 1) = -errp + smsp * dxp * 0.393;
    matms(jp, jp + 1) = smsp * dxp * 0.393;
    materr(jp, jp + 1) = -errp;
  }

  //define covariance contribution non-diagonal between point jp and jp+2
  if (jp < matsize - 2) {
    auto pos3 = tr.LocationAtPoint(jp + 2);
    double dxpp = (pos3 - pos2).R();
    double ppp = tr.MomentumAtPoint(jp + 2);
    double bpp = ppp / sqrt(ppp * ppp + MMU * MMU);
    //double epp = breakpoints[jp + 2] - breakpoints[jp + 1];
    //spp = point(jp + 2).Sin2B();
    //double spp = 1;
    double smspp = sms / p0 / ppp / bpp / b0;
    double errpp = serr * sp * sp / ep / ep;
    //dxpp = dxmedio;
    double cpp = errpp / dxp / dxpp;
    //double smspp = sms * p0 / sqrt(ppp * p0);
    mat(jp, jp + 2) = cpp + smspp * (dxp + dxpp) / 2 * 0.0148;
    matms(jp, jp + 2) = smspp * (dxp + dxpp) / 2 * 0.0148;
    materr(jp, jp + 2) = cpp;
  }

  //guarantee input matrix mat, matms, materr are symmetric
  if (jp < matsize - 1) {
    mat(jp + 1, jp) = mat(jp, jp + 1);
    matms(jp + 1, jp) = matms(jp, jp + 1);
    materr(jp + 1, jp) = materr(jp, jp + 1);
  }
  if (jp < matsize - 2) {
    mat(jp + 2, jp) = mat(jp, jp + 2);
    matms(jp + 2, jp) = matms(jp, jp + 2);
    materr(jp + 2, jp) = materr(jp, jp + 2);
  }
}

//fill covariance matrix mix terms
const void TrajectoryMCSFitterICARUS::FillCovMixTerms(recob::TrackTrajectory tr, TMatrixDSym mat, int jp, int ns, double sms, double serr) const {
  //define index previous to input index
  int ip = jp - 1;
  //define dimension of input matrix
  //int np = mat.GetNrows();
  //double dxmedio = 76.7179;

  //define segment before and after input point (segmentation point) 
  auto pos0 = tr.LocationAtPoint(jp - 1);
  auto pos1 = tr.LocationAtPoint(jp);
  auto pos2 = tr.LocationAtPoint(jp + 1);
  double dxm = (pos1 - pos0).R();
  double dxp = (pos2 - pos1).R();
  
  //define momentum at input index
  double p0 = tr.MomentumAtPoint(jp);

  //define wash coefficient
  double wash0 = 0.735;
  double wash1 = 0.055;
  double washfit = 0.859;
  double c0;

  //update mix terms in input matrix mat
  if (ip < ns - 1) {
    c0 = sms * (dxp + dxm) / 2 / p0 / p0 * washfit * wash0;
    mat(ip, ip + ns - 2) = c0;
    //mat[ip + ns - 2][ip] = c0;
  }
  if (ip < ns - 2) {
    c0 = sms * (dxp + dxm) / 2 / p0 / p0 * washfit * wash0;
    mat(ip, ip + ns - 2) = c0;
  }
  if (ip < ns - 3) {
    c0 = sms * (dxp + dxm) / 2 / p0 / p0 * washfit * wash1;
    mat(ip, ip + ns - 3) = c0;
    //mat[ip + ns][ip] = c0;
  }
  if (ip < ns) {
    c0 = sms * (dxp + dxm) / 2 / p0 / p0 * washfit * wash1;
    mat(ip, ip + ns) = c0;
    //mat[ip + ns - 3][ip] = c0;
  }
}

//perform fit of c2 function to determine MCS momentum
const TrajectoryMCSFitterICARUS::ScanResult TrajectoryMCSFitterICARUS::C2Fit(vector<float>& dthetaLin, vector<float>& dthetaPoly, vector<float>& seg_nradlengths, vector<float>& cumLen, vector<size_t>& breaks, int pid, float sigma, const recob::TrackTrajectory& traj) const {
  //define number of steps
  int nMom = ceil((pMax_ - pMin_) / pStep_) + 1;
  //define vector momentum, c2 function, error on momentum, error on c2 function 
  TVectorD wmom(nMom);
  TVectorD wc2(nMom);
  TVectorD wsmom(nMom);
  TVectorD wsigma(nMom);

  //iterate over p_test, compute c2 function and populate vectors wmom, wc2, wsmom, wsigma
  int jMom = 0; int firstValid = 0; bool allZero = true;
  for (double p_test = pMin_; p_test <= pMax_; p_test += pStep_) {
    cout << "test momentum [MeV/c] = " << p_test << endl; cout << " " << endl;
    double c2 = C2Function(traj, cumLen, breaks, dthetaLin, dthetaPoly, p_test, mass(pid));
    cout << "c2 function = " << c2 << endl; cout << " " << endl;
    wmom[jMom] = p_test / 1000.; 
    wc2[jMom] = c2; 
    wsmom[jMom] = 1;
    wsigma[jMom] = sigma * c2; 
    if (firstValid < 0.01 && c2 > 0.01) firstValid = jMom;
    if (c2 > 0.) allZero = false;
    jMom++;
    cout << "firstValid = " << firstValid << " jMom = " << jMom << " nMom = " << nMom << endl;
  }

  //return null is c2 function is null for all values of p test
  if (allZero) {
    cout << "c2 function is zero everywhere! end fit" << endl;
    return ScanResult(0, 0, 0);
}

  //define new vectors for momentum, c2 function, error on momentum, error on c2 function only for c2!=0
  TVectorD rmom(nMom - firstValid);
  TVectorD rsmom(nMom - firstValid);
  TVectorD rsigma(nMom - firstValid);
  TVectorD rc2(nMom - firstValid);
  for (int jp = firstValid; jp < jMom; jp++) {
    rmom[jp - firstValid] = wmom[jp];
    rc2[jp - firstValid] = wc2[jp];
    rsmom[jp - firstValid] = wsmom[jp];
    rsigma[jp - firstValid] = wsigma[jp];
  }
  
  //define once again new vectors for momentum, c2 function, error on momentum, error on c2 function 
  const TVectorD cmom = rmom;
  const TVectorD cc2 = rc2;
  const TVectorD csmom = rsmom;
  const TVectorD csigma = rsigma;

  //root class for graph with x = cmom, y = cc2, dx = csmom, dy = csigma
  TGraphErrors *gr3 = new TGraphErrors(cmom, cc2, csmom, csigma);
  //define minimum momentum, maximum momentum, number of parameters of fit function
  double pmin = wmom[firstValid]; double pmax = pMax_; int npar = 2;
  //root class for fit to function 1 / (par[0] + par[1] / (x[0] * x[0]))
  TF1* fitfunc = new TF1("fitfunc", funzio, pmin, pmax, npar);
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
  double best_p = sqrt(beta / (1 - alpha));
  double error_p = sqrt(pow(dbeta, 2.) / (4 * beta * (1 - alpha)) + beta * pow(dalpha, 2.) / (4 * pow((1 - alpha), 3.)));

  //in case of best p being too low, search where c2 function change sign
  if (best_p < 0.001) {
    int over = -1;
    if (rc2[0] > 1) over = 0;
    for (int jp = 1; jp < nMom - firstValid; jp++) {
      if (rc2[jp] > 1. && rc2[jp - 1] < 1.) over = jp;
    }
    if (over > 0 && over < nMom - firstValid) best_p = rmom[over - 1] + (rmom[over] - rmom[over - 1]) / (rc2[over] - rc2[over - 1]);
    if (over == 0) best_p = rmom[0]; //underflow
    if (over == -1) best_p = rmom[nMom - firstValid - 1]; //overflow

  }

  //return best p, error p 
  cout << "best momentum = " << best_p << " and error momentum = " << error_p << endl;
  return ScanResult(best_p, error_p, 0.);
}

//find length on collection plane
double TrajectoryMCSFitterICARUS::collLength() const {
  //define first hit on collection plane
  recob::Hit h0 = hits2dC[0];
  //define index associated to last hit on collection plane
  int jf = hits2dC.size() - 1;
  //define number of TPC associated to first hit
  int NS0 = h0.WireID().TPC;

  //float x0 = h0.WireID().Wire * 3; auto y0 = h0.PeakTime() * 0.622;

  //find index of last hit in the same TPC of first hit
  int jb = jf;
  for (unsigned int jh = 0; jh < hits2dC.size(); jh++) {
    //define number of TPC associated to current hit
    int NS = hits2dC[jh].WireID().TPC;
    if (NS != NS0) jb = jh - 1;
  }

  //define last hit - last hit in the same TPC of first hit
  recob::Hit hb = hits2dC[jb];
  //define afterlast hit - hit after last hit in the same TPC, in case there is one
  //note there is supposed to be one in any case, since we are currently looking to cathode-crossing tracks
  recob::Hit hp = hits2dC[jb + 1];
  //define very last hit - last hit, regardless of TPC
  recob::Hit hf = hits2dC[jf];

  //define wire ID of first hit, last hit and very last hit
  int w0 = h0.WireID().Wire;
  int wb = hb.WireID().Wire;
  int wp = hp.WireID().Wire;
  int wf = hf.WireID().Wire;
  //compute distance on w and s coordinates between last and first hit, in the same TPC
  float dw = (wb - w0) * 3.;
  cout << "distance between first and last hit, on w coordinate = " << dw << endl;
  float ds = (hb.PeakTime() - h0.PeakTime()) * 0.622;
  cout << "distance between first and last hit, on s coordinate = " << ds << endl;
  
  //in case there is hit after last hit in the same TPC, compute correction to distance
  if (jb != jf) {
    float dw_corr = (wf - wp) * 3.; 
    float ds_corr = (hf.PeakTime() - hp.PeakTime()) * 0.622;  
    dw += dw_corr;
    ds += ds_corr;
  }

  //print and save corrected distance between first and very last hit
  cout << "corrected distance between first and last hit, on w coordinate = " << dw << endl;
  cout << "corrected distance between first and last hit, on s coordinate = " << ds << endl;
  cout << "corrected modulus distance between first and last hit = " << sqrt(dw * dw + ds * ds) << endl;
  return sqrt(dw * dw + ds * ds);
}

//find projection on wire coordinate of length on collection plane
double TrajectoryMCSFitterICARUS::collWireLength() {
  //define first and last hits on collection plane
  recob::Hit h0 = hits2dC[0];
  recob::Hit hf = hits2dC[hits2dC.size() - 1];
  //compute wire coordinate on collection plane of first and last hit
  float x0 = h0.WireID().Wire * 3; 
  float xf = hf.WireID().Wire * 3; 
  //return distance on collection plane between first and last hit
  return xf - x0;
}

//find cosine of angle between 1D drift length and 3D length track
double TrajectoryMCSFitterICARUS::cosTrackDrift(recob::TrackTrajectory tr) const {
  //define start point and end point of track
  auto start = tr.Start();
  auto end = tr.End();
  //compute 3D and 1D length of track (1D over drift direction)
  auto l3d = sqrt((end - start).Mag2()); 
  auto l1d = (end - start).X();
  //return cosine of angle between 1D drift length track and 3D length track
  return l1d / l3d;
}

//define rotation matrix depending on plane, tpc and cryostat
TMatrixD TrajectoryMCSFitterICARUS::ReferenceFrame(int plane) const { 
  //define rotation matrix 3x3 
  TMatrixD mat(3, 3);

  //populate rotation matrix for view 0 (induction 1)
  if (plane == 0) {
    mat(0, 0) = 1; mat(0, 1) = 0; mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = 1; mat(1, 2) = 0;
    mat(2, 0) = 0; mat(2, 1) = 0; mat(2, 2) = 1;
  }

  //populate rotation matrix for view 1 (induction 2)
  if (plane == 1) {
    mat(0, 0) = 1; mat(0, 1) = 0;              mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = cos(M_PI / 3);  mat(1, 2) = sin(M_PI / 3);
    mat(2, 0) = 0; mat(2, 1) = -sin(M_PI / 3); mat(2, 2) = cos(M_PI / 3);
  }

  //populate rotation matrix for view 2 (collection)
  if (plane == 2) {
    mat(0, 0) = 1; mat(0, 1) = 0;             mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = cos(M_PI / 3); mat(1, 2) = -sin(M_PI / 3);
    mat(2, 0) = 0; mat(2, 1) = sin(M_PI / 3); mat(2, 2) = cos(M_PI / 3);
  }

  return mat;
}

//define drift origin coordinate 
double TrajectoryMCSFitterICARUS::DriftOrigin(int plane, int tpc, int cryo) const { 
  //initialize drift coordinate of anode
  double x0;

  //cryostat EAST (0), tpc EES (0) and EEN (1)
  if (cryo == 0 && tpc <= 1) x0 = -3596.3;
  //cryostat EAST (0), tpc EWS (2) and EWN (3)
  if (cryo == 0 && tpc > 1) x0 = -608.0;
  //cryostat WEST (1), tpc WES (0) and WEN (1)
  if (cryo == 1 && tpc <= 1) x0 = 608.0;
  //cryostat WEST (1), tpc WWS (2) and WWN (3)
  if (cryo == 1 && tpc > 1) x0 = 3596.3;

  //return drift coordinate of anode
  return x0;
}

//define cathode drift coordinates 
void TrajectoryMCSFitterICARUS::CathodeDistance(int cryo, double x0, double x1) const { 
  //cryostat EAST (0)
  if (cryo == 0) { x0 = -235.; x1 = -185.; }
  //cryostat WEST (1)
  if (cryo == 1) { x0 = 185.; x1 = 235.; }
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
      if (jrow < jtail && jcol < jtail) {
	      covmod(jrow, jcol) = cov(jrow, jcol);
	    }
      if(jrow > jtail && jcol > jtail) {
      	covmod(jrow - 1, jcol - 1) = cov(jrow, jcol);
	    }
    }
  }

  //return modified matrix covmod without original row and column jtail
  return covmod;
}