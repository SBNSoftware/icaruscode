/**
 * @file   icaruscode/PMT/Calibration/TrackTools/StoppingTrackSelector.h
 * @brief  Selects stopping-muon tracks and matches them in time to OpFlashes.
 * @author Matteo Vicenzi
 *
 * The selection is not reimplemented here: this class builds the `sbn::TrackInfo`
 * that `sbn::ITCSSelectionTool` expects and hands it to the very same art tool that
 * is used for the calibration ntuples (`sbncode/Calibration/TrackCaloSkimmerSelectStoppingTrack_tool.cc`).
 * This is configured with `@local::crthittagged_stopping_selection` from `icarus_trackcalo_skimmer.fcl`
 * so the cut values cannot drift from what is used in the calibration ntuples.
 *
 *  Only the `sbn::TrackInfo` members the stopping tool actually reads are filled:
 * `dir.y`, `end`, `hit_{min,max}_time_p2_tpc{E,W}`, `hits2[].{oncalo,rr,dqdx}` and
 * `whicht0`. Everything else keeps its default.
 */

#ifndef ICARUSCODE_PMT_CALIBRATION_TRACKTOOLS_STOPPINGTRACKSELECTOR_H
#define ICARUSCODE_PMT_CALIBRATION_TRACKTOOLS_STOPPINGTRACKSELECTOR_H

// framework
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

// SBN
#include "sbnobj/Common/Calibration/TrackCaloSkimmerObj.h"
#include "sbncode/Calibration/ITCSSelectionTool.h"

#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace icarus {

  // ---------------------------------------------------------------------------
  /**
   * @brief One selected stopping track, with its matched flash if there is one.
   * 
   * All times are in **microseconds**. `anab::T0::Time()` is in nanoseconds and is
   * converted once, inside `StoppingTrackSelector::select()`.
   */
  struct TrackFlashMatch {

    static constexpr double NoTime = -99999.0;
    static constexpr float  NoPos  = -999999.0f;

    // track
    int      trackID   = -1;      ///< `recob::Track::ID()`.
    unsigned cryostat  = 999;     ///< From the first track hit, not from the label.
    int      whichT0   = -1;      ///< TrackCaloSkimmer convention: 0 Pandora, 2 CRTHit.
    bool     selected  = false;   ///< Result of `sbn::ITCSSelectionTool::DoSelect()`.
    double   trackT0   = NoTime;  ///< Best track time [us].
    float    startX = NoPos, startY = NoPos, startZ = NoPos; /// Track start, drift-corrected [cm].
    float    endX = NoPos, endY = NoPos, endZ = NoPos; ///< Track end, drift-corrected [cm].
    float    dirY = NoPos;        ///< Start direction, y component.
    float    centroidY = NoPos, centroidZ = NoPos; /// Charge-weighted centroid of the track's plane-2 hits [cm].
    float    length = NoPos;      ///< Track length [cm].
    float    medianEnddQdx = NoPos; ///< Median dQ/dx of the last hits: what the cut uses.

    // matched flash
    int    flashID   = -1;      ///< Index in the OpFlash collection; -1 if unmatched.
    double flashTime = NoTime;  ///< `recob::OpFlash::Time()` [us].
    double deltaT    = NoTime;  ///< `flashTime - trackT0` [us].
    float  flashPE   = NoPos;
    float  flashY    = NoPos;
    float  flashZ    = NoPos; 
    float  radius    = NoPos;   ///< flashZY - centroidZY [cm]

    /// Whether a flash was matched to this track.
    bool hasFlash() const { return flashID >= 0; }

  }; // struct TrackFlashMatch


  // ---------------------------------------------------------------------------
  /**
   * @brief Runs the calibration stopping-track selection, then matches in time.
   *
   * One instance covers every cryostat; `select()` is called once per cryostat
   * index, which indexes into the configured label sequences.
   *
   */
  class StoppingTrackSelector {

  public:

    explicit StoppingTrackSelector(fhicl::ParameterSet const& pset);

    /// Number of configured cryostats.
    std::size_t nCryostats() const { return fTrackLabels.size(); }

    /// Selects stopping tracks in one cryostat and matches them to flashes.
    std::vector<TrackFlashMatch> select(art::Event const& e, std::size_t iCryo);

  private:

    /// What `buildTrackInfo()` returns alongside the `sbn::TrackInfo` itself:
    /// quantities the selection tool does not read but the match does.
    struct TrackExtras {
      double t0_ns     = std::numeric_limits<double>::signaling_NaN();
      float  centroidY = TrackFlashMatch::NoPos;
      float  centroidZ = TrackFlashMatch::NoPos;
    };

    /// Per-event associations, built once in `select()` and reused for every track.
    struct EventAssns;

    /// Fills the `sbn::TrackInfo` members the stopping tool reads.
    /// Returns false if the track cannot be evaluated at all.
    bool buildTrackInfo(art::Event const& e, EventAssns const& assns,
                        std::size_t iCryo,
                        art::Ptr<recob::PFParticle> const& pfp,
                        art::Ptr<recob::Track> const& track,
                        sbn::TrackInfo& info, TrackExtras& extra) const;

    /// Picks the best flash by `MatchMetric`, or leaves `m.flashID` at -1.
    void matchFlash(std::vector<recob::OpFlash> const& flashes,
                    TrackFlashMatch& m) const;

    /// Median dQ/dx over plane-2 hits with `oncalo && rr < MediandQdxRRMax`.
    /// Same quantity the stopping tool cuts on; diagnostic only.
    float medianEnddQdx(std::vector<sbn::TrackHitInfo> const& hits) const;

    /// Whether a point is inside any cryostat active volume.
    bool pointInActiveVolume(geo::Point_t const& p) const;

    // --- configuration -------------------------------------------------------
    std::vector<art::InputTag> fTrackLabels;   ///< `recob::Track`, one per cryostat.
    std::vector<art::InputTag> fPFPLabels;     ///< `recob::PFParticle`, one per cryostat.
    std::vector<art::InputTag> fPFPT0Labels;   ///< Pandora `anab::T0`, one per cryostat.
    std::vector<art::InputTag> fCaloLabels;    ///< `anab::Calorimetry`, one per cryostat.
    std::vector<art::InputTag> fFlashLabels;   ///< `recob::OpFlash`, one per cryostat.
    art::InputTag fCRTHitT0Label;              ///< CRT hit tagging, both cryostats.

    bool  fAllowShowerLikePFPs;   ///< Keep PFParticles with `PdgCode() == 11`.
    bool  fIncludeCRTHitTagging;  ///< Consider the CRT-hit T0 at all.
    bool  fIncludeTopCRT;         ///< Accept `CRTHitT0TaggingInfo::Sys == 0`.
    bool  fIncludeSideCRT;        ///< Accept `CRTHitT0TaggingInfo::Sys == 1`.

    /// CRT-hit distance vetoes, one value per cryostat. They are not the same
    /// in the two cryostats: `icarus_trackcalo_skimmer.fcl:120-134` 
    std::vector<double> fTopCRTDistanceCutStopping;
    std::vector<double> fTopCRTDistanceCutPassing;
    std::vector<double> fSideCRTDistanceCutStopping;
    std::vector<double> fSideCRTDistanceCutPassing;

    double fMediandQdxRRMax;      ///< [cm]; mirrors the tool, diagnostic only.

    bool   fUseTimeWindow;        /// // whether the window rejects flashes, off by default
    double fMatchWindowLow;       ///< [us] minimum `flashTime - trackT0`.
    double fMatchWindowHigh;      ///< [us] maximum `flashTime - trackT0`.
    double fMinFlashPE;
    double fMaxFlashTrackZ;       ///< [cm] `|flash.ZCenter() - centroidZ|`.
    std::string fMatchMetric;     ///< "MinRadius" or "MinDeltaT".
    bool   fRequireFlashMatch;    ///< Drop selected tracks with no flash.

    std::string fLogCategory;

    // --- state ---------------------------------------------------------------
    std::unique_ptr<sbn::ITCSSelectionTool> fSelTool; ///< Built once; holds state.
    std::vector<geo::BoxBoundedGeo> fActiveVolumes;   ///< One per cryostat.

  }; // class StoppingTrackSelector

} // namespace icarus

#endif // ICARUSCODE_PMT_CALIBRATION_TRACKTOOLS_STOPPINGTRACKSELECTOR_H
