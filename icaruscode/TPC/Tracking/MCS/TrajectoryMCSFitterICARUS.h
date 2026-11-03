#ifndef TrajectoryMCSFitterICARUS_H
#define TrajectoryMCSFitterICARUS_H

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/RecoObjects/TrackState.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "lardata/RecoBaseProxy/Track.h" //needed only if you do use the proxies

namespace trkf {
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
	      5 };
      fhicl::Atom<int> minNumAngles {
        Name("minNumAngles"),
	      Comment("Minimum number of angles the track must have to perform fit correctly."),
	      4 };
      fhicl::Atom<double> segmentLength {
        Name("segmentLength"),
	      Comment("Nominal length of track segments used in the fit."),
	      10. };
      fhicl::Atom<int> minHitsPerSegment {
        Name("minHitsPerSegment"),
	      Comment("Exclude segments with less hits than this value."),
	      3 };
      fhicl::Atom<int> minHits {
        Name("minHits"),
	      Comment("Exclude tracks with less hits than this value."),
	      30 };
      fhicl::Atom<double> minLength {
        Name("minLength"),
	      Comment("Exclude tracks with length less than this value."),
	      100 };
      fhicl::Atom<int> eLossMode {
        Name("eLossMode"),
	      Comment("0 for MPV Landau, 1 for MIP, 2 for Bethe-Bloch."),
	      2 };
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
      fhicl::Atom<unsigned int> planeMode{
        Name("planeMode"),
        Comment("0 for Induction-1, 1 for Induction-2, 2 for Collection view."),
        0 };
      fhicl::Atom<unsigned int> fitMode{
        Name("fitMode"),
        Comment("0 for both linear and polygonal angles, 1 for only linear angles, 2 for only polygonal angles. Only for C2 function computation."),
        0 };
      fhicl::Atom<unsigned int> fullTrackMode{
        Name("fullTrackMode"),
        Comment("0 for portion of track only in the first TPC, 1 for portion of track only in the last TPC, 2 for full track."),
        0 };
      fhicl::Atom<unsigned int> selectionMode {
        Name("selectionMode"),
        Comment("0 for cosmic stopping tracks. Mode for tracks not stopping and not cosmic is not available yet."),
        0 };
      fhicl::Atom<unsigned int> wireSkip{
        Name("wireSkip"),
        Comment("Number of wires to be skipped in sigma3p computation."),
        0 };
      fhicl::Atom<unsigned int> verbose{
        Name("verbose"),
        Comment("0 for not printing out anything, 1 for printing out something, 2 for printing out everything."),
        0 }; };

    using Parameters = fhicl::Table<Config>;

    TrajectoryMCSFitterICARUS(
      int pIdHyp, 
      int minNSegs, 
      int minNAngs,
      double segLen, 
      int minHitsPerSegment, 
      int minHits, 
      double minLength,
      int eLossMode, 
      double pMin, 
      double pMax, 
      double pStep, 
      unsigned int cutMode, 
      unsigned int cutAngles, 
      unsigned int planeMode, 
      unsigned int fitMode,
      unsigned int fullTrackMode,
      unsigned int selectionMode,
      unsigned int wireSkip,
      unsigned int verbose) {
        pIdHyp_ = pIdHyp;
        minNSegs_ = minNSegs;
        minNAngs_ = minNAngs;
        segLen_ = segLen;
        minHitsPerSegment_ = minHitsPerSegment;
        minHits_ = minHits;
        minLength_ = minLength;
        eLossMode_ = eLossMode;
        pMin_ = pMin;
        pMax_ = pMax;
        pStep_ = pStep;
        cutMode_ = cutMode;
        cutAngles_ = cutAngles;
        planeMode_ = planeMode;
        fitMode_ = fitMode; 
        fullTrackMode_ = fullTrackMode;
        selectionMode_ = selectionMode; 
        wireSkip_ = wireSkip; 
        verbose_ = verbose; }

    explicit TrajectoryMCSFitterICARUS(
      const Parameters & p) : TrajectoryMCSFitterICARUS(
        p().pIdHypothesis(),
        p().minNumSegments(),
        p().minNumAngles(),
        p().segmentLength(),
        p().minHitsPerSegment(),
        p().minHits(),
        p().minLength(),
        p().eLossMode(),
        p().pMin(),
        p().pMax(),
        p().pStep(),
        p().cutMode(),
        p().cutAngles(),
        p().planeMode(),
        p().fitMode(),
        p().fullTrackMode(),
        p().selectionMode(),
        p().wireSkip(),
        p().verbose()) {}
    
    recob::MCSFitResultGS fitMcs(
      const recob::TrackTrajectory& traj) const { 
        return fitMcs(traj, pIdHyp_); }

    recob::MCSFitResultGS fitMcs(
      const recob::Track& track) const { 
        return fitMcs(track, pIdHyp_); }

    recob::MCSFitResultGS fitMcs(
      const recob::Trajectory& traj) const { 
        return fitMcs(traj, pIdHyp_); }

    recob::MCSFitResultGS fitMcs(
      const recob::TrackTrajectory& traj, 
      int pid) const;

    recob::MCSFitResultGS fitMcs(
      const recob::Track& track, 
      int pid) const { 
        return fitMcs(track.Trajectory(), pid); }

    recob::MCSFitResultGS fitMcs(
      const recob::Trajectory& traj, 
      int pid) const {
        recob::TrackTrajectory::Flags_t flags(traj.NPoints());
        const recob::TrackTrajectory tt(traj, std::move(flags));
        return fitMcs(tt, pid); }

    struct ScanResult {
      public:
        double bestp, errp, minp, maxp;
        double alpha, dalpha, beta, dbeta;
        ScanResult() = default;
        ScanResult(
          double best, double err, double min, double max,
          double a, double da, double b, double db,
          std::vector<float> test, std::vector<float> c2) :
            bestp(best), errp(err), minp(min), maxp(max),
            alpha(a), dalpha(da), beta(b), dbeta(db) {} };

    // functions in MCS fitter

    void breakTrajInSegments(
      const recob::TrackTrajectory& traj, 
      std::vector<size_t>& breakpoints, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<int>& cumseghits,
      std::vector<int>& segsides) const;

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

    double energyLossLandau(
      const double mass, 
      const double energy) const;

    double energyLossBetheBloch(
      const double mass, 
      const double energy) const;

    double GetE(
      const double initial_E, 
      const double length_travelled, 
      const double mass) const;

    void ComputeSigma3P(
      int& nhits,
      double& sigma3p, 
      std::vector<double>& delta3p, 
      int wire_skip) const;
    
    static Double_t funzio(
      Double_t *x, 
      Double_t *par) { return 1. / (par[0] + par[1] / x[0] / x[0]); } 

    const double C2Function(
      const recob::TrackTrajectory& traj, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      double sigma3p, 
      int pid, 
      double p) const;

    const ScanResult C2Fit(
      const recob::TrackTrajectory& traj, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      double sigma3p, 
      int pid, 
      double sigma) const;

    const ScanResult C2NaiveFit(
      const recob::TrackTrajectory& traj, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      double sigma3p, 
      int pid, 
      double sigma) const;

    const double LikelihoodFunction(
      const recob::TrackTrajectory& traj, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      double sigma3p, 
      int pid, 
      double p0,
      bool verbose) const;

    void LikelihoodFit(
      const recob::TrackTrajectory& traj, 
      std::vector<float>& seglens, 
      std::vector<float>& cumseglens, 
      std::vector<int>& seghits, 
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      double sigma3p, 
      int pid, 
      double& best_p,
      double& error_p) const;
    
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

    TMatrix CleanCovariance(
      TMatrix cov, 
      int jtail) const;

    unsigned int Plane(
      size_t index) const;

    bool IsInPlane(
      size_t index) const;

    unsigned int TPC(
      size_t index) const;

    bool IsInTPC(
      size_t index, 
      unsigned int t) const;

    unsigned int FirstTPC(
      const recob::TrackTrajectory& traj) const;

    unsigned int LastTPC(
      const recob::TrackTrajectory& traj) const;

    Vector_t hit2D(
      size_t index) const;

    double length1D(
      const recob::TrackTrajectory traj) const;

    double length2D(
      const recob::TrackTrajectory traj) const;

    double length3D(
      const recob::TrackTrajectory traj) const;

    double ThetaMCSFactor(
      const recob::TrackTrajectory& traj) const;

    double ThetaErrFactor(
      const recob::TrackTrajectory& traj) const;

    double ThetaExpected(
      const recob::TrackTrajectory& traj, 
      double p,
      double beta,
      double L) const;

    void ThetaCheck(
      std::vector<float> dthetaLin, 
      std::vector<float> dthetaPoly, 
      bool& checkLin, 
      bool& checkPoly, 
      unsigned int& firstseg, 
      unsigned int& lastseg) const;
    
    bool CathodeCheck(
      const recob::TrackTrajectory& traj,
      size_t index) const;

    // functions in MCS producer

    void set2DHitsC(
      std::vector<recob::Hit> h) { hits2dC = h; }

    void set2DHitsI2(
      std::vector<recob::Hit> h) { hits2dI2 = h; }

    void set2DHitsI1(
      std::vector<recob::Hit> h) { hits2dI1 = h; }

    void setPointData(
      std::vector<proxy::TrackPointData> h) { pdata = h; }

    void setWR(
      geo::WireReadoutGeom const& wr) { wr_ = &wr; }
      
    void setDetProp(
      detinfo::DetectorPropertiesData const& detProp) { detProp_ = &detProp; }

    void setGeom(
      geo::GeometryCore const* geom) { geom_ = geom; }

    void setPFPT0(
      float t0) { PFPT0_ = t0; }
      
    void setCRTT0(
      float t0) { CRTT0_ = t0; }

    void setT0(
      float t0) { T0_ = t0; }

    void sethasPFPT0(
      bool has) { hasPFPT0_ = has; }
      
    void sethasCRTT0(
      bool has) { hasCRTT0_ = has; }
      
    void sethasT0(
      bool has) { hasT0_ = has; }

    void setRangeP(
      float r) { rangeP = r; }

    double mass (
      int pid) const {
        if (abs(pid) == 13) return mumass;
        if (abs(pid) == 211) return pimass;
        if (abs(pid) == 321) return kmass;
        if (abs(pid) == 2212) return pmass;
        return util::kBogusD; }

    double minLength() const {
      return minLength_; }

    unsigned int cutMode() const { 
      return cutMode_; }

    unsigned int cutAngles() const { 
      return cutAngles_; }

    unsigned int planeMode() const { 
      return planeMode_; }

    unsigned int fitMode() const { 
      return fitMode_; }

    unsigned int fullTrackMode() const {
      return fullTrackMode_; }

    unsigned int selectionMode() const { 
      return selectionMode_; }

    unsigned int wireSkip() const { 
      return wireSkip_; }

    unsigned int verbose() const { 
      return verbose_; }

  private:
    int pIdHyp_;
    unsigned int minNSegs_;
    int minNAngs_;
    double segLen_;
    int minHitsPerSegment_;
    int minHits_;
    double minLength_;
    int eLossMode_;
    double pMin_;
    double pMax_;
    double pStep_;
    unsigned int cutMode_;
    unsigned int cutAngles_;
    unsigned int planeMode_;
    unsigned int fitMode_;
    unsigned int fullTrackMode_;
    unsigned int selectionMode_;
    unsigned int wireSkip_;
    unsigned int verbose_;
    std::vector<recob::Hit> hits2dC;
    std::vector<recob::Hit> hits2dI2;
    std::vector<recob::Hit> hits2dI1;
    std::vector<proxy::TrackPointData> pdata;
    geo::WireReadoutGeom const* wr_;
    geo::GeometryCore const* geom_;
    detinfo::DetectorPropertiesData const* detProp_;
    float PFPT0_; float CRTT0_; float T0_;
    bool hasPFPT0_; bool hasCRTT0_; bool hasT0_;
    float rangeP; }; }

#endif