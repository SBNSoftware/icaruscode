#ifndef TrajectoryMCSFitterICARUS_H
#define TrajectoryMCSFitterICARUS_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h"
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
   * @brief Class for C2-function fit of Multiple Coulomb Scattering angles between segments within a Track - or Trajectory.
   *
   * Input are a Track - or Trajectory - and various parameters: pIdHypothesis, minNumSegments, minNumAngles, segmentLength, minHitsPerSegment, minHits, nElossSteps, eLossMode, pMin, pMax, pStep, cutMode, cutAngles, dimMode, planeMode, fitMode.
   *
   * Output is a recob::MCSFitResultGS, containing:
   *   resulting momentum, momentum uncertainty, and best likelihood value (both for fwd and bwd fit);
   *   vector of segment (radiation) lengths, vector of scattering angles, and PID hypothesis used in the fit.
   *
   * Note that the comulative segment length is what is used to compute the energy loss, but the segment length is actually slightly different, so the output can be used to reproduce the original results but they will not be identical (but very close).
   *
   * @author  G. Chiello (University of Pisa)
   * @date    2025
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
	      13 };
      fhicl::Atom<int> minNumSegments {
        Name("minNumSegments"),
	      Comment("Minimum number of segments the track is split into."),
	      7 };
      fhicl::Atom<int> minNumAngles {
        Name("minNumAngles"),
	      Comment("Minimum number of angles the track should have to perform fit correctly."),
	      3 };
      fhicl::Atom<double> segmentLength {
        Name("segmentLength"),
	      Comment("Nominal length of track segments used in the fit."),
	      14. };
      fhicl::Atom<int> minHitsPerSegment {
        Name("minHitsPerSegment"),
	      Comment("Exclude segments with less hits than this value."),
	      3 };
      fhicl::Atom<int> minHits {
        Name("minHits"),
	      Comment("Exclude tracks with less hits than this value."),
	      30 };
      fhicl::Atom<int> nElossSteps {
        Name("nElossSteps"),
	      Comment("Number of steps for computing energy loss uptream to current segment."),
	      10 };
      fhicl::Atom<int> eLossMode {
        Name("eLossMode"),
	      Comment("0 (default) for MPV Landau, 1 for MIP, 2 for Bethe-Bloch."),
	      0 };
      fhicl::Atom<double> pMin {
        Name("pMin"),
	      Comment("Minimum momentum value [MeV/c] in C2-function fit."),
	      100 };
      fhicl::Atom<double> pMax {
        Name("pMax"),
	      Comment("Maximum momentum value [MeV/c] in C2-function fit."),
	      1500 };
      fhicl::Atom<double> pStep {
        Name("pStep"),
	      Comment("Step in momentum value [MeV/c] in C2-function fit."),
	      10 };
      fhicl::Atom<unsigned int> cutMode{
        Name("cutMode"),
        Comment("0 for full track, 1 for cutting final X angles, 2 for keeping initial X angles. X is given by cutAngles."),
        0 };
      fhicl::Atom<unsigned int> cutAngles{
        Name("cutAngles"),
        Comment("Number of angles to cut for MCS fit."),
        0 };
      fhicl::Atom<unsigned int> dimMode{
        Name("dimMode"),
        Comment("2 for 2D, 3 for 3D. Flag of dimension mode for ICARUS fit."),
        0 };
      fhicl::Atom<unsigned int> planeMode{
        Name("planeMode"),
        Comment("0 for Induction-1, 1 for Induction-2, 2 for Collection view. Flag of plane mode for ICARUS 2D fit."),
        0 };
      fhicl::Atom<unsigned int> fitMode{
        Name("fitMode"),
        Comment("0 for both linear and polygonal angles, 1 for only linear angles, 2 for only polygonal angles. Flag of fit mode for ICARUS fit."),
        0 };
      fhicl::Atom<bool> removeDeltas{
        Name("removeDeltas"),
        Comment("If true, remove identified deltas from MCS calculation."),
        0 }; };

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
      unsigned int fitMode,
      bool removeDeltas) {
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
        removeDeltas_ = removeDeltas;}

    explicit TrajectoryMCSFitterICARUS(
      const Parameters & p) : TrajectoryMCSFitterICARUS(
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
        p().fitMode(),
        p().removeDeltas()) {}
    
    recob::MCSFitResultGS fitMcs(
      const recob::TrackTrajectory& traj, 
      bool momDepConst = true) const { 
        return fitMcs(traj, pIdHyp_, momDepConst); }

    recob::MCSFitResultGS fitMcs(
      const recob::Track& track, 
      bool momDepConst = true) const { 
        return fitMcs(track, pIdHyp_, momDepConst); }

    recob::MCSFitResultGS fitMcs(
      const recob::Trajectory& traj, 
      bool momDepConst = true) const { 
        return fitMcs(traj, pIdHyp_, momDepConst); }

    recob::MCSFitResultGS fitMcs(
      const recob::TrackTrajectory& traj, 
      int pid, 
      bool momDepConst = true) const;

    recob::MCSFitResultGS fitMcs(
      const recob::Track& track, 
      int pid, 
      bool momDepConst = true) const { 
        return fitMcs(track.Trajectory(), pid, momDepConst); }

    recob::MCSFitResultGS fitMcs(
      const recob::Trajectory& traj, 
      int pid, 
      bool momDepConst = true) const {
        recob::TrackTrajectory::Flags_t flags(traj.NPoints());
        const recob::TrackTrajectory tt(traj, std::move(flags));
        return fitMcs(tt, pid, momDepConst); }

    void breakTrajInSegments(
      const recob::TrackTrajectory& traj, 
      std::vector<size_t>& breakpoints, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<int>& cumseghits,
      std::vector<int>& isdi) const;

    void findSegmentBarycenter(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& bary,
      std::vector<int>& isdi) const;

    void linearRegression(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& pcdir,
      std::vector<int>& isdi) const;

    void find2DSegmentBarycenter(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& bary2D,
      std::vector<int>& isdi) const;

    void linearRegression2D(
      const recob::TrackTrajectory& traj, 
      const size_t firstPoint, 
      const size_t lastPoint, 
      recob::tracking::Vector_t& pcdir2D,
      std::vector<int>& isdi) const;

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
        double bestp, errp, minp, maxp;
        double alpha, dalpha, beta, dbeta;
        std::vector<float> testp, c2function;
        std::vector<int> tailssize;
        ScanResult() = default;
        ScanResult(
          double best, double err, double min, double max,
          double a, double da, double b, double db,
          std::vector<float> test, std::vector<float> c2, std::vector<int> ts) :
            bestp(best), errp(err), minp(min), maxp(max),
            alpha(a), dalpha(da), beta(b), dbeta(db),
            testp(test), c2function(c2), tailssize(ts) {} };

    void set2DHitsC(
      std::vector<recob::Hit> h) { hits2dC = h; }

    void set2DHitsI2(
      std::vector<recob::Hit> h) { hits2dI2 = h; }

    void set2DHitsI1(
      std::vector<recob::Hit> h) { hits2dI1 = h; }

    void setPointData(
      std::vector<proxy::TrackPointData> h) { pdata = h; }

    void setCRTShift(
      float s) { CRTshift = s; }

    void setRangeP(
      float r) { rangeP = r; }

    double mass (
      int pid) const {
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

    bool removeDeltas() const { 
      return removeDeltas_; }

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
      double p0,
      int& tailssize) const;
    
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

    double Integral(
      const recob::TrackTrajectory& traj, 
      size_t index) const;

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

    double ThetaExpected(
      const recob::TrackTrajectory& traj, 
      double p,
      double beta,
      double L) const;

    double PrintD3P() const;

    void ThetaCheck(
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      bool& checkLin, 
      bool& checkPoly, 
      unsigned int& firstseg, 
      unsigned int& lastseg) const;

    bool GeoStopCheck(
      const recob::TrackTrajectory& traj) const;
    
    bool DeltaCheck(
      size_t index) const;

    double CathodeDistance(
      unsigned int c) const;
    
    bool CathodeCheck(
      const recob::TrackTrajectory& traj,
      size_t index) const;

    std::vector<float> CleanThetaLin(
      std::vector<float> dthetaLin) const;

    std::vector<float> CleanThetaPoly(
      std::vector<float> dthetaPoly) const;

    float C2PRange() { return c2prange; }

    int TailsPRange() { return tailsprange; }

    float distance_point_line(
      const double a0, 
      const double a1, 
      const double xdr, 
      const double ydr) const;
    
    void ProcessDeltaRays(
      const recob::TrackTrajectory& traj, 
      int viewType, 
      std::vector<int>& isd, 
      std::vector<int>& isdi) const;
    
    void Create_ZY_arrays();
    
    void Delete_ZY_arrays();
    
    void Fill_ZY_arrays(
      std::vector<double>& Z, 
      std::vector<double>& Y);
    
    void TagOverlappingDeltaRays(
      const recob::TrackTrajectory& traj,
      int view_type, 
      std::vector<int>& isDelta,
      std::vector<int>& isdi, 
      float threshold) const;
    
    void TagDeltaRaysLocal();
    
    void TagDeltaRays(
      const recob::TrackTrajectory& traj, 
      int viewType, 
      std::vector<int>& isDelta,
      std::vector<int>& isdi) const;
    
    void TagDeltaRaysCylinder(
      const recob::TrackTrajectory& traj, 
      int viewType, 
      std::vector<int>& isDelta,
      std::vector<int>& isdi) const;
    
    void TagCloseToDeltaRays();
    
    std::vector<int> HitsOnWire(
      std::vector<recob::Hit> hits,
      unsigned int iWire) const;

    void runToyMS() const;

    void linearRegressionToy(
      std::vector<double> xx, 
      std::vector<double> yy, 
      std::vector<double> zz, 
      Vector_t& pcdir) const ;

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
    bool removeDeltas_;
    std::vector<recob::Hit> hits2dC;
    std::vector<recob::Hit> hits2dI2;
    std::vector<recob::Hit> hits2dI1;
    std::vector<proxy::TrackPointData> pdata;
    float d3pC; float d3pI1; float d3pI2; 
    float CRTshift;
    float rangeP;
    float c2prange; 
    int tailsprange; }; }

#endif