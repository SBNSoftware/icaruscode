#ifndef TrajectoryMCSFitterICARUS_H
#define TrajectoryMCSFitterICARUS_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/TrackState.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "lardata/RecoBaseProxy/Track.h" //needed only if you do use the proxies

namespace trkf {
  /**
   * @file  larreco/RecoAlg/TrajectoryMCSFitterICARUS.h
   * @class trkf::TrajectoryMCSFitterICARUS
   *
   * @brief Class for Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory
   *
   * Class for Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory.
   *
   * Inputs are: a Track or Trajectory, and various fit parameters (pIdHypothesis, minNumSegments, segmentLength, pMin, pMax, pStep)
   *
   * Outputs are: a recob::MCSFitResult, containing:
   *   resulting momentum, momentum uncertainty, and best likelihood value (both for fwd and bwd fit);
   *   vector of segment (radiation) lengths, vector of scattering angles, and PID hypothesis used in the fit.
   *
   * For configuration options see TrajectoryMCSFitterICARUS#Configs
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */
  class TrajectoryMCSFitterICARUS {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> pIdHypothesis {
        Name("pIdHypothesis"),
	      Comment("Default particle ID hypothesis to be used in the fit when not specified."),
	      13
      };
      fhicl::Atom<int> minNumSegments {
        Name("minNumSegments"),
	      Comment("Minimum number of segments the track is split into."),
	      7
      };
      fhicl::Atom<int> minNumAngles {
        Name("minNumAngles"),
	      Comment("Minimum number of angles the track should have to perform fit correctly."),
	      3
      };
      fhicl::Atom<double> segmentLength {
        Name("segmentLength"),
	      Comment("Nominal length of track segments used in the fit."),
	      14.
      };
      fhicl::Atom<int> minHitsPerSegment {
        Name("minHitsPerSegment"),
	      Comment("Exclude segments with less hits than this value."),
	      3
      };
      fhicl::Atom<int> minHits {
        Name("minHits"),
	      Comment("Exclude tracks with less hits than this value."),
	      30
      };
      fhicl::Atom<int> nElossSteps {
        Name("nElossSteps"),
	      Comment("Number of steps for computing energy loss uptream to current segment."),
	      10
      };
      fhicl::Atom<int> eLossMode {
        Name("eLossMode"),
	      Comment("0 (default) for MPV Landau, 1 for MIP, 2 for Bethe-Bloch."),
	      0
      };
      fhicl::Atom<double> pMin {
        Name("pMin"),
	      Comment("Minimum momentum value in likelihood scan, in units of MeV/c."),
	      100
      };
      fhicl::Atom<double> pMax {
        Name("pMax"),
	      Comment("Maximum momentum value in likelihood scan, in units of MeV/c."),
	      1500
      };
      fhicl::Atom<double> pStep {
        Name("pStep"),
	      Comment("Step in momentum value in likelihood scan, in units of MeV/c."),
	      10
      };
      fhicl::Atom<unsigned int> cutMode{
        Name("cutMode"),
        Comment("0 for full track, 1 for cutting final X angles, 2 for keeping initial X angles. X is given by cutAngles."),
        0
      };
      fhicl::Atom<unsigned int> cutAngles{
        Name("cutAngles"),
        Comment("Number of angles to cut for MCS fitting."),
        0
      };
      fhicl::Atom<unsigned int> dimMode{
        Name("dimMode"),
        Comment("2 for 2D, 3 for 3D. Flag of dimension mode for ICARUS fit."),
        0
      };
      fhicl::Atom<unsigned int> planeMode{
        Name("planeMode"),
        Comment("0 for Induction-1, 1 for Induction-2, 2 for Collection view. Flag of plane mode for ICARUS 2D fit."),
        0
      };
      fhicl::Atom<unsigned int> fitMode{
        Name("fitMode"),
        Comment("0 for both linear and polygonal angles, 1 for only linear angles, 2 for only polygonal angles. Flag of fit mode for ICARUS fit."),
        0
      };
    };
    using Parameters = fhicl::Table<Config>;

    TrajectoryMCSFitterICARUS(
      int pIdHyp, 
      int minNSegs, 
      int minNAngs,
      double segLen, 
      int minHitsPerSegment, 
      int minHits, 
      int nElossSteps, 
      int eLossMode, 
      double pMin, 
      double pMax, 
      double pStep, 
      unsigned int cutMode, 
      unsigned int cutAngles, 
      unsigned int dimMode, 
      unsigned int planeMode, 
      unsigned int fitMode) {
        pIdHyp_ = pIdHyp;
        minNSegs_ = minNSegs;
        minNAngs_ = minNAngs;
        segLen_ = segLen;
        minHitsPerSegment_ = minHitsPerSegment;
        minHits_ = minHits;
        nElossSteps_ = nElossSteps;
        eLossMode_ = eLossMode;
        pMin_ = pMin;
        pMax_ = pMax;
        pStep_ = pStep;
        cutMode_ = cutMode;
        cutAngles_ = cutAngles;
        dimMode_ = dimMode;
        planeMode_ = planeMode;
        fitMode_ = fitMode;
    }
    explicit TrajectoryMCSFitterICARUS(const Parameters & p):TrajectoryMCSFitterICARUS(
      p().pIdHypothesis(),
      p().minNumSegments(),
      p().minNumAngles(),
      p().segmentLength(),
      p().minHitsPerSegment(),
      p().minHits(),
      p().nElossSteps(),
      p().eLossMode(),
      p().pMin(),
      p().pMax(),
      p().pStep(),
      p().cutMode(),
      p().cutAngles(),
      p().dimMode(),
      p().planeMode(),
      p().fitMode()){}
    
    recob::MCSFitResult fitMcs(
      const recob::TrackTrajectory& traj, 
      bool momDepConst = true) const { 
        return fitMcs(traj, pIdHyp_, momDepConst); }

    recob::MCSFitResult fitMcs(
      const recob::Track& track, 
      bool momDepConst = true) const { 
        return fitMcs(track, pIdHyp_, momDepConst); }

    recob::MCSFitResult fitMcs(
      const recob::Trajectory& traj, 
      bool momDepConst = true) const { 
        return fitMcs(traj, pIdHyp_, momDepConst); }

    recob::MCSFitResult fitMcs(
      const recob::TrackTrajectory& traj, 
      int pid, 
      bool momDepConst = true) const;

    recob::MCSFitResult fitMcs(
      const recob::Track& track, 
      int pid, 
      bool momDepConst = true) const { 
        return fitMcs(track.Trajectory(), pid, momDepConst); }

    recob::MCSFitResult fitMcs(
      const recob::Trajectory& traj, 
      int pid, 
      bool momDepConst = true) const {
        recob::TrackTrajectory::Flags_t flags(traj.NPoints());
        const recob::TrackTrajectory tt(traj, std::move(flags));
        return fitMcs(tt, pid, momDepConst);
    }

    void breakTrajInSegments(
      const recob::TrackTrajectory& traj, 
      std::vector<size_t>& breakpoints, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<int>& cumseghits) const;

    void findSegmentBarycenter(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& bary) const;

    void linearRegression(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& pcdir) const;

    void find2DSegmentBarycenter(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& bary2D) const;

    void linearRegression2D(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& pcdir2D) const;

    double GetOptimalSegLen(
      const recob::TrackTrajectory& traj, 
      const double guess_p, 
      const int n_points, 
      const int plane, 
      const double length_travelled) const;

    double computeResidual(
      int i, 
      double& alfa, 
      std::vector<recob::Hit> h) const;

    double computeResidual3D(
      recob::TrackTrajectory traj, 
      int i, 
      double& alfa) const;

    void ComputeD3P(
      int plane);

    float ComputeD3P3D(
      const recob::TrackTrajectory& tr) const;

    struct ScanResult {
      public:
        ScanResult(
          double ap, 
          double apUnc, 
          double alogL): p(ap), pUnc(apUnc), logL(alogL) {}
        double p, pUnc, logL; };

    void set2DHitsC(std::vector<recob::Hit> h) {hits2dC = h;}
    void set2DHitsI2(std::vector<recob::Hit> h) {hits2dI2 = h;}
    void set2DHitsI1(std::vector<recob::Hit> h) {hits2dI1 = h;}
    void setPointData(std::vector<proxy::TrackPointData> h) {pdata = h;}
    void setCRTShift(float s) {CRTshift = s;}

    double mass (int pid) const {
      if (abs(pid) == 13) return mumass;
      if (abs(pid) == 211) return pimass;
      if (abs(pid) == 321) return kmass;
      if (abs(pid) == 2212) return pmass;
      return util::kBogusD; }

    unsigned int cutMode() const { 
      return cutMode_; }
    unsigned int cutAngles() const { 
      return cutAngles_; }
    unsigned int dimMode() const { 
      return dimMode_; }
    unsigned int planeMode() const { 
      return planeMode_; }
    unsigned int fitMode() const { 
      return fitMode_; }

    double energyLossLandau(
      const double mass2, 
      const double e2, 
      const double x) const;

    double energyLossBetheBloch(
      const double mass, 
      const double e2) const;

    double GetE(
      const double initial_E, 
      const double length_travelled, 
      const double mass) const;

    const double C2Function(
      const recob::TrackTrajectory& traj, 
      std::vector<size_t>& breakpoints, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<int>& cumseghits, 
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      int pid, 
      double p0) const;
    
    const void FillOffDiagMCSLin(
      TMatrixDSym& mat, 
      TMatrixD matdiag, 
      int jp) const;

    const void FillOffDiagMCSPoly(
      TMatrixDSym& mat, 
      TMatrixD matdiag, 
      int jp) const;

    const void FillOffDiagErrLin(
      TMatrixDSym& mat, 
      TMatrixD matdiag, 
      int jp) const;

    const void FillOffDiagErrPoly(
      TMatrixDSym& mat, 
      TMatrixD matdiag, 
      int jp) const;  

    const void FillMCSMixedTerms(
      TMatrixD& matMix, 
      TMatrixD matLin, 
      TMatrixD matPoly) const; 
    
    static Double_t funzio(
      Double_t *x, 
      Double_t *par) { return 1. / (par[0] + par[1] / x[0] / x[0]); } 

    const ScanResult C2Fit(
      const recob::TrackTrajectory& traj, 
      std::vector<size_t>& breakpoints, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<int>& cumseghits, 
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      int pid, 
      double sigma) const;

    void AnodeDistance(
      int cryo, 
      int tpc, 
      double x0) const;

    void CathodeDistance(
      int cryo, 
      double x0) const;

    TMatrixD ReferenceFrame(
      int plane) const;

    TMatrix CleanCovariance(
      TMatrix cov, 
      int jtail) const;

    bool isinplane(
      size_t index, 
      unsigned int p) const;

    bool isintpc(
      size_t index, 
      unsigned int t) const;

    unsigned int lasttpc(
      const recob::TrackTrajectory& traj) const;

    Vector_t hit2d(
      const recob::TrackTrajectory& traj, 
      size_t index, 
      unsigned int p, 
      unsigned int t) const;

    double length1D(
      const recob::TrackTrajectory traj, 
      unsigned int plane) const;

    double length2D(
      const recob::TrackTrajectory traj, 
      unsigned int plane) const;

    double length3D(
      const recob::TrackTrajectory traj, 
      unsigned int plane) const;

    double ThetaMCSFactor(
      const recob::TrackTrajectory& traj, 
      unsigned int plane) const;

    double ThetaErrFactor(
      const recob::TrackTrajectory& traj, 
      unsigned int plane) const;

    double PrintD3P() const;

    void ThetaCheck(
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      bool& checkLin, 
      bool& checkPoly, 
      unsigned int& firstseg, 
      unsigned int& lastseg) const;

    void GeoStopCheck(
      const recob::TrackTrajectory& traj) const;

  private:
    int pIdHyp_;
    unsigned int minNSegs_;
    int minNAngs_;
    double segLen_;
    int minHitsPerSegment_;
    int minHits_;
    int nElossSteps_;
    int eLossMode_;
    double pMin_;
    double pMax_;
    double pStep_;
    unsigned int cutMode_;
    unsigned int cutAngles_;
    unsigned int dimMode_;
    unsigned int planeMode_;
    unsigned int fitMode_;
    std::vector<recob::Hit> hits2dC;
    std::vector<recob::Hit> hits2dI2;
    std::vector<recob::Hit> hits2dI1;
    std::vector<proxy::TrackPointData> pdata;
    float d3pC; float d3pI1; float d3pI2; 
    float CRTshift; }; }

#endif