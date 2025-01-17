#ifndef TrajectoryMCSFitterUBoone_H
#define TrajectoryMCSFitterUBoone_H

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "canvas/Persistency/Common/Ptr.h"

#include "lardata/RecoObjects/TrackState.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "lardata/RecoBaseProxy/Track.h" //needed only if you do use the proxies

#include <array>
#include <utility>
#include <vector>

namespace trkf {
  /**
   * @file  larreco/RecoAlg/TrajectoryMCSFitterUBoone.h
   * @class trkf::TrajectoryMCSFitterUBoone
   *
   * @brief Class for Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory
   *
   * Class for Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory.
   *
   * Inputs are: a Track or Trajectory, and various fit parameters (pIdHypothesis, minNumSegments, segmentLength, pMin, pMax, pStep, angResol)
   *
   * Outputs are: a recob::MCSFitResult, containing:
   *   resulting momentum, momentum uncertainty, and best likelihood value (both for fwd and bwd fit);
   *   vector of comulative segment (radiation) lengths, vector of scattering angles, and PID hypothesis used in the fit.
   *   Note that the comulative segment length is what is used to compute the energy loss, but the segment length is actually slightly different,
   *   so the output can be used to reproduce the original results but they will not be identical (but very close).
   *
   * For configuration options see TrajectoryMCSFitterUBoone#Config
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */
  class TrajectoryMCSFitterUBoone {
    //
  public:
    //
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> pIdHypothesis{
        Name("pIdHypothesis"),
        Comment("Default particle Id Hypothesis to be used in the fit when not specified."),
        13};
      fhicl::Atom<unsigned int> minNumSegments{
        Name("minNumSegments"),
        Comment("Minimum number of segments the track is split into."),
        3};
      fhicl::Atom<double> segmentLength{
        Name("segmentLength"),
        Comment("Nominal length of track segments used in the fit."),
        14.};
      fhicl::Atom<int> minHitsPerSegment{
        Name("minHitsPerSegment"),
        Comment("Exclude segments with less hits than this value."),
        2};
      fhicl::Atom<int> nElossSteps{
        Name("nElossSteps"),
        Comment("Number of steps for computing energy loss uptream to current segment."),
        10};
      fhicl::Atom<int> eLossMode{
        Name("eLossMode"),
        Comment("Default is MPV Landau. Choose 1 for MIP (constant); 2 for Bethe-Bloch."),
        0};
      fhicl::Atom<double> pMin{Name("pMin"),
                               Comment("Minimum momentum value in likelihood scan."),
                               0.01};
      fhicl::Atom<double> pMax{Name("pMax"),
                               Comment("Maximum momentum value in likelihood scan."),
                               7.50};
      fhicl::Atom<double> pStep{Name("pStep"),
                                Comment("Step in momentum value in fine grained likelihood scan."),
                                0.01};
      fhicl::Atom<double> angResol {
        Name("angResol"),
        Comment("Angular resolution parameter used in modified Highland formula. Unit is mrad."),
        10.0};
      fhicl::Sequence<double> fiducialVolumeInsets {
        Name("fiducialVolumeInsets"),
        Comment("insets from the TPC boundaries to exclude from the fit")};
      fhicl::Sequence<double> excludeVolumes {
        Name("excludeVolumes"),
        Comment("regions to exclude")};
      fhicl::Atom<int> cutMode{
        Name("cutMode"),
        Comment("Flag of track cutting mode. 0=full track, 1=cut final XX cm, 2=keep initial XX cm. XX is given by cutLength variables"),
        false};
      fhicl::Atom<float> cutLength{
        Name("cutLength"),
        Comment("initial or final cutting length to use for MCS"),
        false};
    };
    using Parameters = fhicl::Table<Config>;
    //
    TrajectoryMCSFitterUBoone(int pIdHyp,
                              unsigned int minNSegs,
                              double segLen,
                              int minHitsPerSegment,
                              int nElossSteps,
                              int eLossMode,
                              double pMin,
                              double pMax,
                              double pStep,
                              double angResol, 
                              std::vector<double>fiducialVolumeInsets, 
                              std::vector<double> excludeVolumes,
                              int cutMode,
                              float cutLength) {
      pIdHyp_ = pIdHyp;
      minNSegs_ = minNSegs;
      segLen_ = segLen;
      minHitsPerSegment_ = minHitsPerSegment;
      nElossSteps_ = nElossSteps;
      eLossMode_ = eLossMode;
      pMin_ = pMin;
      pMax_ = pMax;
      pStep_ = pStep;
      angResol_ = angResol;
      fiducialVolumeInsets_ = fiducialVolumeInsets;
      if(fiducialVolumeInsets_.size() != 6) fiducialVolumeInsets_ = std::vector<double>(6, 0.);
      excludeVolumes_ = excludeVolumes;
      cutMode_ = cutMode;
      cutLength_ = cutLength;
    }
    explicit TrajectoryMCSFitterUBoone(const Parameters& p) 
    : TrajectoryMCSFitterUBoone(p().pIdHypothesis(),
                                p().minNumSegments(),
                                p().segmentLength(),
                                p().minHitsPerSegment(),
                                p().nElossSteps(),
                                p().eLossMode(),
                                p().pMin(),
                                p().pMax(),
                                p().pStep(),
                                p().angResol(),
                                p().fiducialVolumeInsets(),
                                p().excludeVolumes(),
                                p().cutMode(),
                                p().cutLength()) {}
    //
    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj, bool momDepConst = true) const { return fitMcs(traj,pIdHyp_,momDepConst); }
    recob::MCSFitResult fitMcs(const recob::Track& track,          bool momDepConst = true) const { return fitMcs(track,pIdHyp_,momDepConst); }
    recob::MCSFitResult fitMcs(const recob::Trajectory& traj,      bool momDepConst = true) const { return fitMcs(traj,pIdHyp_,momDepConst); }
    //
    recob::MCSFitResult fitMcs(const recob::TrackTrajectory& traj, int pid, bool momDepConst = true) const;
    recob::MCSFitResult fitMcs(const recob::Track& track,          int pid, bool momDepConst = true) const { return fitMcs(track.Trajectory(),pid,momDepConst); }
    recob::MCSFitResult fitMcs(const recob::Trajectory& traj,      int pid, bool momDepConst = true) const {
      recob::TrackTrajectory::Flags_t flags(traj.NPoints());
      const recob::TrackTrajectory tt(traj,std::move(flags));
      return fitMcs(tt,pid,momDepConst);
    }
    //
    void breakTrajInSegments(const recob::TrackTrajectory& traj,
                             std::vector<size_t>& breakpoints,
                             std::vector<float>& segradlengths,
                             std::vector<float>& cumseglens, 
                             std::vector<bool>& breakpointsgood,
                             int cutMode, 
                             float cutLength) const;
    void linearRegression(const recob::TrackTrajectory& traj,
                          const size_t firstPoint,
                          const size_t lastPoint,
                          recob::tracking::Vector_t& pcdir) const;
    double mcsLikelihood(double p,
                         double theta0x,
                         std::vector<float>& dthetaij,
                         std::vector<float>& seg_nradl,
                         std::vector<float>& cumLen,
                         bool fwd,
                         bool momDepConst,
                         int pid) const;
    //
    struct ScanResult {
      public:
        ScanResult(double ap, double apUnc, double alogL) : p(ap), pUnc(apUnc), logL(alogL) {}
        double p, pUnc, logL;
    };
    //
    void set2DHitsC(std::vector<recob::Hit> h) {hits2dC = h;}
    void set2DHitsI2(std::vector<recob::Hit> h) {hits2dI2 = h;}
    void set2DHitsI1(std::vector<recob::Hit> h) {hits2dI1 = h;}
    void setPointData(std::vector<proxy::TrackPointData> h) {pdata = h;}
    //
    const ScanResult doLikelihoodScan(std::vector<float>& dtheta,
                                      std::vector<float>& seg_nradlengths,
                                      std::vector<float>& cumLen,
                                      bool fwdFit,
                                      bool momDepConst,
                                      int pid) const;
    //
    inline double MomentumDependentConstant(const double p) const {
      //these are from https://arxiv.org/abs/1703.06187
      constexpr double a = 0.1049;
      constexpr double c = 11.0038;
      return (a / (p * p)) + c;
    }
    double mass(int pid) const
    {
      if (abs(pid) == 13) { return mumass; }
      if (abs(pid) == 211) { return pimass; }
      if (abs(pid) == 321) { return kmass; }
      if (abs(pid) == 2212) { return pmass; }
      return util::kBogusD;
    }
    double energyLossBetheBloch(const double mass, const double e2) const;
    double energyLossLandau(const double mass2, const double E2, const double x) const;
    //
    double GetE(const double initial_E, const double length_travelled, const double mass) const;
    //
    unsigned int minNSegs() const { return minNSegs_; }
    double segLen() const { return segLen_; }
    //
    std::vector<geo::BoxBoundedGeo> setFiducialVolumes() const;
    std::vector<geo::BoxBoundedGeo> setExcludeVolumes() const;
    bool isInVolume(const std::vector<geo::BoxBoundedGeo> &volumes, const geo::Point_t &point) const;

    bool isintpc(
      size_t index, 
      unsigned int t) const;

    unsigned int lasttpc(
      const recob::TrackTrajectory& traj) const;
      
    //
    int cutMode() const { return cutMode_; }
    float cutLength() const { return cutLength_; }
    //
  private:
    int pIdHyp_;
    unsigned int minNSegs_;
    double segLen_;
    int minHitsPerSegment_;
    int nElossSteps_;
    int eLossMode_;
    double pMin_;
    double pMax_;
    double pStep_;
    double angResol_;
    std::vector<double> fiducialVolumeInsets_;
    std::vector<double> excludeVolumes_;
    int cutMode_;
    float cutLength_;
    std::vector<recob::Hit> hits2dC;
    std::vector<recob::Hit> hits2dI2;
    std::vector<recob::Hit> hits2dI1;
    std::vector<proxy::TrackPointData> pdata;
  };
}

#endif