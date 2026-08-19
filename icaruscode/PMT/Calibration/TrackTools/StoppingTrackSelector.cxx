/**
 * @file   icaruscode/PMT/Calibration/TrackTools/StoppingTrackSelector.cxx
 * @brief  Implementation of icarus::StoppingTrackSelector.
 * @author Matteo Vicenzi
 *
 * Everything that fills `sbn::TrackInfo` mirrors `sbncode/Calibration/TrackCaloSkimmer_module.cc` 
 * line by line, so that the tool sees exactly what it sees in the calibration ntuple job. 
 * 
 * Cross-references to that file are given as `TrackCaloSkimmer_module.cc:<line>`.
 */

#include "icaruscode/PMT/Calibration/TrackTools/StoppingTrackSelector.h"

// framework
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

// SBN
#include "sbnobj/Common/CRT/CRTHitT0TaggingInfo.hh"

#include <algorithm>
#include <cmath>

// -----------------------------------------------------------------------------
namespace {

  /// Minimum time of plane-2 hits in TPC E (`TPCE`) or W, ICARUS convention.
  /// Verbatim from `TrackCaloSkimmer_module.cc:488-508`: the -1 initial value is
  /// load bearing, the tool reads a negative value as "no hits here, skip the cut".
  float hitMinTime(std::vector<sbn::TrackHitInfo> const& hits, bool TPCE)
  {
    double min = -1;
    for (sbn::TrackHitInfo const& h : hits) {
      bool const hit_is_TPCE = h.h.tpc <= 1; // ICARUS: TPC 0,1 are E within the cryostat
      if (h.oncalo && hit_is_TPCE == TPCE) {
        if (min < 0. || h.h.time < min) min = h.h.time;
      }
    }
    return min;
  }

  /// As `hitMinTime()`, for the maximum. `TrackCaloSkimmer_module.cc:512-532`.
  float hitMaxTime(std::vector<sbn::TrackHitInfo> const& hits, bool TPCE)
  {
    double max = -1;
    for (sbn::TrackHitInfo const& h : hits) {
      bool const hit_is_TPCE = h.h.tpc <= 1;
      if (h.oncalo && hit_is_TPCE == TPCE) {
        if (max < 0. || h.h.time > max) max = h.h.time;
      }
    }
    return max;
  }

} // local namespace


// -----------------------------------------------------------------------------
icarus::StoppingTrackSelector::StoppingTrackSelector(fhicl::ParameterSet const& pset)
  : fTrackLabels     { pset.get<std::vector<art::InputTag>>("TrackLabels") }
  , fPFPLabels       { pset.get<std::vector<art::InputTag>>("PFPLabels") }
  , fPFPT0Labels     { pset.get<std::vector<art::InputTag>>("PFPT0Labels") }
  , fCaloLabels      { pset.get<std::vector<art::InputTag>>("CaloLabels") }
  , fFlashLabels     { pset.get<std::vector<art::InputTag>>("FlashLabels") }
  , fCRTHitT0Label   { pset.get<art::InputTag>("CRTHitT0Label", art::InputTag{}) }
  // every default below is the TrackCaloSkimmer default (TrackCaloSkimmer_module.cc:88-95),
  // deliberately: a default that differs from the skimmer's is a silent selection bias.
  , fAllowShowerLikePFPs { pset.get<bool>("AllowShowerLikePFPs", false) }
  , fIncludeCRTHitTagging{ pset.get<bool>("IncludeCRTHitTagging", false) }
  , fIncludeTopCRT       { pset.get<bool>("IncludeTopCRT", false) }
  , fIncludeSideCRT      { pset.get<bool>("IncludeSideCRT", false) }
  , fTopCRTDistanceCutStopping {
      pset.get<std::vector<double>>("TopCRTDistanceCut_stopping",
        std::vector<double>(fTrackLabels.size(), 100.)) }
  , fTopCRTDistanceCutPassing {
      pset.get<std::vector<double>>("TopCRTDistanceCut_throughgoing",
        std::vector<double>(fTrackLabels.size(), 100.)) }
  , fSideCRTDistanceCutStopping {
      pset.get<std::vector<double>>("SideCRTDistanceCut_stopping",
        std::vector<double>(fTrackLabels.size(), 100.)) }
  , fSideCRTDistanceCutPassing {
      pset.get<std::vector<double>>("SideCRTDistanceCut_throughgoing",
        std::vector<double>(fTrackLabels.size(), 100.)) }
  , fMediandQdxRRMax { pset.get<fhicl::ParameterSet>("SelectionTool")
                         .get<double>("MediandQdxRRMax", 5.) }
  , fUseTimeWindow   { pset.get<bool>("UseTimeWindow", false) }
  , fMatchWindowLow  { pset.get<double>("MatchWindowLow", -30.) }
  , fMatchWindowHigh { pset.get<double>("MatchWindowHigh", 30.) }
  , fMinFlashPE      { pset.get<double>("MinFlashPE", 0.) }
  , fMaxFlashTrackZ  { pset.get<double>("MaxFlashTrackZ", 1e9) }
  , fMatchMetric     { pset.get<std::string>("MatchMetric", "MinRadius") }
  , fRequireFlashMatch { pset.get<bool>("RequireFlashMatch", true) }
  , fLogCategory     { pset.get<std::string>("LogCategory", "StoppingTrackSelector") }
{

  // every per-cryostat sequence must have the same length
  std::size_t const n = fTrackLabels.size();
  if (fPFPLabels.size() != n || fPFPT0Labels.size() != n
      || fCaloLabels.size() != n || fFlashLabels.size() != n)
  {
    throw art::Exception(art::errors::Configuration)
      << "StoppingTrackSelector: TrackLabels, PFPLabels, PFPT0Labels, CaloLabels and"
         " FlashLabels must all have the same size (one entry per cryostat).\n";
  }

  // the CRT distance vetoes are per cryostat too: the two ICARUS cryostats are
  // configured with different values upstream, so a scalar cannot express them
  if (fTopCRTDistanceCutStopping.size()  != n || fTopCRTDistanceCutPassing.size()  != n
      || fSideCRTDistanceCutStopping.size() != n || fSideCRTDistanceCutPassing.size() != n)
  {
    throw art::Exception(art::errors::Configuration)
      << "StoppingTrackSelector: the four {Top,Side}CRTDistanceCut_* sequences must"
         " each have one entry per cryostat (" << n << " expected).\n";
  }

  if (fMatchMetric != "MinRadius" && fMatchMetric != "MinDeltaT") {
    throw art::Exception(art::errors::Configuration)
      << "StoppingTrackSelector: unknown MatchMetric '" << fMatchMetric
      << "'; expected 'MinRadius' or 'MinDeltaT'.\n";
  }

  // the selection tool itself: built once, it carries a mutable prescale counter
  fSelTool = art::make_tool<sbn::ITCSSelectionTool>(
    pset.get<fhicl::ParameterSet>("SelectionTool"));

  // active volumes, one per cryostat, built exactly as the stopping tool does
  // (TrackCaloSkimmerSelectStoppingTrack_tool.cc:79-107). Used for the CRT-hit
  // distance veto, which picks its cut depending on whether the track stops.
  geo::GeometryCore const* geometry = lar::providerFrom<geo::Geometry>();
  for (auto const& cryoid : geometry->Iterate<geo::CryostatID>()) {
    std::vector<geo::BoxBoundedGeo> tpcs;
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryoid)) {
      tpcs.push_back(TPC.ActiveBoundingBox());
    }
    if (tpcs.empty()) continue;

    auto const byMinX = [](auto const& l, auto const& r) { return l.MinX() < r.MinX(); };
    auto const byMinY = [](auto const& l, auto const& r) { return l.MinY() < r.MinY(); };
    auto const byMinZ = [](auto const& l, auto const& r) { return l.MinZ() < r.MinZ(); };
    auto const byMaxX = [](auto const& l, auto const& r) { return l.MaxX() < r.MaxX(); };
    auto const byMaxY = [](auto const& l, auto const& r) { return l.MaxY() < r.MaxY(); };
    auto const byMaxZ = [](auto const& l, auto const& r) { return l.MaxZ() < r.MaxZ(); };

    fActiveVolumes.emplace_back(
      std::min_element(tpcs.begin(), tpcs.end(), byMinX)->MinX(),
      std::max_element(tpcs.begin(), tpcs.end(), byMaxX)->MaxX(),
      std::min_element(tpcs.begin(), tpcs.end(), byMinY)->MinY(),
      std::max_element(tpcs.begin(), tpcs.end(), byMaxY)->MaxY(),
      std::min_element(tpcs.begin(), tpcs.end(), byMinZ)->MinZ(),
      std::max_element(tpcs.begin(), tpcs.end(), byMaxZ)->MaxZ());
  }

} // icarus::StoppingTrackSelector::StoppingTrackSelector()


// -----------------------------------------------------------------------------
bool icarus::StoppingTrackSelector::pointInActiveVolume(geo::Point_t const& p) const
{
  for (geo::BoxBoundedGeo const& av : fActiveVolumes) {
    if (av.ContainsPosition(p)) return true;
  }
  return false;
}


// -----------------------------------------------------------------------------
float icarus::StoppingTrackSelector::medianEnddQdx
  (std::vector<sbn::TrackHitInfo> const& hits) const
{
  // same collection the tool medians over (TrackCaloSkimmerSelectStoppingTrack_tool.cc:139-158)
  std::vector<float> dqdx;
  for (sbn::TrackHitInfo const& h : hits) {
    if (h.oncalo && h.rr < fMediandQdxRRMax) dqdx.push_back(h.dqdx);
  }
  if (dqdx.empty()) return TrackFlashMatch::NoPos;

  std::sort(dqdx.begin(), dqdx.end());
  std::size_t const n = dqdx.size();
  return (n % 2 == 1) ? dqdx[n / 2] : 0.5f * (dqdx[n / 2 - 1] + dqdx[n / 2]);
}


// -----------------------------------------------------------------------------
/// Associations are expensive to build, so they are made once per event and
/// per cryostat rather than once per track.
struct icarus::StoppingTrackSelector::EventAssns {

  EventAssns(art::Event const& e,
             art::Handle<std::vector<recob::PFParticle>> const& pfpHandle,
             art::Handle<std::vector<recob::Track>> const& trackHandle,
             art::InputTag const& trackTag, art::InputTag const& caloTag,
             art::InputTag const& pfpT0Tag, art::InputTag const& crtT0Tag)
    // same tags the skimmer uses (TrackCaloSkimmer_module.cc:245-247)
    : hits   { trackHandle, e, trackTag }
    , calo   { trackHandle, e, caloTag }
    , t0PFP  { pfpHandle,   e, pfpT0Tag }
    , t0CRT  { trackHandle, e, crtT0Tag }
    , crtTag { trackHandle, e, crtT0Tag }
    {}

  art::FindManyP<recob::Hit, recob::TrackHitMeta> hits;
  art::FindManyP<anab::Calorimetry>               calo;
  art::FindManyP<anab::T0>                        t0PFP;   ///< on the PFParticle
  art::FindManyP<anab::T0>                        t0CRT;   ///< on the track
  art::FindManyP<sbn::crt::CRTHitT0TaggingInfo>   crtTag;  ///< on the track
};


// -----------------------------------------------------------------------------
bool icarus::StoppingTrackSelector::buildTrackInfo(
  art::Event const& e, EventAssns const& assns, std::size_t iCryo,
  art::Ptr<recob::PFParticle> const& pfp,
  art::Ptr<recob::Track> const& track,
  sbn::TrackInfo& info, TrackExtras& extra) const
{

  if (!assns.hits.isValid()) return false;

  auto const& hits = assns.hits.at(track.key());
  auto const& thms = assns.hits.data(track.key());
  if (hits.empty()) return false;

  info.id       = track->ID();
  info.length   = track->Length();
  info.cryostat = hits.front()->WireID().Cryostat; // :1211-1213

  info.dir.x = track->StartDirection().X();
  info.dir.y = track->StartDirection().Y();
  info.dir.z = track->StartDirection().Z();

  info.start.y = track->Start().Y();
  info.start.z = track->Start().Z();
  info.end.y   = track->End().Y();
  info.end.z   = track->End().Z();
  // NB: end.x is set below, per T0 branch. With no T0 it stays a signaling NaN,
  // which is what makes the fiducial-X test fail rather than pass by accident.

  // --- T0, following TrackCaloSkimmer_module.cc:352-430 -----------------------
  double t0Pandora = std::numeric_limits<double>::signaling_NaN();
  bool hasT0Pandora = false;
  if (assns.t0PFP.isValid() && !assns.t0PFP.at(pfp.key()).empty()) {
    t0Pandora = assns.t0PFP.at(pfp.key()).at(0)->Time(); // [ns]
    hasT0Pandora = true;
  }

  double t0CRTHit = std::numeric_limits<double>::signaling_NaN();
  bool hasT0CRTHit = false;
  if (!hasT0Pandora && fIncludeCRTHitTagging && !fCRTHitT0Label.empty()) {

    if (assns.t0CRT.isValid() && assns.crtTag.isValid()
        && !assns.t0CRT.at(track.key()).empty() && !assns.crtTag.at(track.key()).empty())
    {
      sbn::crt::CRTHitT0TaggingInfo const& tag = *assns.crtTag.at(track.key()).at(0);
      double const time = assns.t0CRT.at(track.key()).at(0)->Time(); // [ns]

      // wall veto (:387)
      bool const sysRejected =
        (tag.Sys == 0 && !fIncludeTopCRT) || (tag.Sys == 1 && !fIncludeSideCRT);

      // distance veto: the cut depends on whether the track stops, and the test
      // point is built from Start() even though the skimmer calls it `end` (:390-397)
      auto const& dprop =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
      geo::GeometryCore const* geometry = lar::providerFrom<geo::Geometry>();

      geo::Point_t testPoint { track->Start().X(), track->Start().Y(), track->Start().Z() };
      int const driftDir = geometry->TPC(hits.at(0)->WireID()).DriftDir().X();
      double const driftv = dprop.DriftVelocity();
      testPoint.SetX(testPoint.X() + time * driftDir * driftv * 1e-3);

      bool const trackIsStopping = pointInActiveVolume(testPoint);
      double const distCut = (tag.Sys == 0)
        ? (trackIsStopping ? fTopCRTDistanceCutStopping[iCryo]
                           : fTopCRTDistanceCutPassing[iCryo])
        : (trackIsStopping ? fSideCRTDistanceCutStopping[iCryo]
                           : fSideCRTDistanceCutPassing[iCryo]);

      if (!sysRejected && !(tag.Distance > distCut)) {
        t0CRTHit = time;
        hasT0CRTHit = true;
      }
    }
  }

  if (hasT0Pandora) {
    info.whicht0 = 0;
    info.start.x = track->Start().X();  // no drift shift with a Pandora T0 (:1176-1180)
    info.end.x   = track->End().X();
    extra.t0_ns = t0Pandora;
  }
  else if (hasT0CRTHit) {
    info.whicht0 = 2;
    auto const& dprop =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
    geo::GeometryCore const* geometry = lar::providerFrom<geo::Geometry>();
    int const driftDir = geometry->TPC(hits.at(0)->WireID()).DriftDir().X();
    double const driftv = dprop.DriftVelocity(dprop.Efield(), dprop.Temperature());
    double const xshift = driftDir * driftv * t0CRTHit * 1e-3; // :1195-1200
    info.start.x   = track->Start().X() + xshift;
    info.end.x     = track->End().X() + xshift;
    info.xShiftCRT = xshift;
    extra.t0_ns = t0CRTHit;
  }
  else {
    info.whicht0 = -1;  // the tool rejects this under RequireT0: true
    extra.t0_ns = std::numeric_limits<double>::signaling_NaN();
  }

  // --- plane-2 hits, following MakeHit (:1380-1500) ---------------------------
  static std::vector<art::Ptr<anab::Calorimetry>> const emptyCalo;
  auto const& calos = assns.calo.isValid() ? assns.calo.at(track.key()) : emptyCalo;

  // charge-weighted centroid, accumulated over the same plane-2 hits. Weight is
  // the raw hit integral, as in TPCPMTBarycenterMatchProducer; only hits with a
  // valid trajectory point contribute, since the others have no 3D position.
  double sumW = 0., sumWY = 0., sumWZ = 0.;

  for (std::size_t i = 0; i < hits.size(); ++i) {

    if (hits[i]->WireID().Plane != 2) continue;

    sbn::TrackHitInfo hinfo;
    hinfo.h.time  = hits[i]->PeakTime();
    hinfo.h.tpc   = hits[i]->WireID().TPC;
    hinfo.h.plane = hits[i]->WireID().Plane;
    hinfo.h.wire  = hits[i]->WireID().Wire;
    hinfo.h.id    = static_cast<int>(hits[i].key());

    // a hit off the fitted trajectory never gets calorimetry attached (:1462-1468).
    // Skipping this test would mark hits `oncalo` that the ntuple marks false and
    // would shift the end-of-track median dQ/dx, hence the selection itself.
    bool const badhit =
      (thms[i]->Index() == std::numeric_limits<unsigned int>::max())
      || (!track->HasValidPoint(thms[i]->Index()));
    hinfo.ontraj = !badhit;

    if (!badhit) {

      geo::Point_t const loc = track->LocationAtPoint(thms[i]->Index());
      double const w = hits[i]->Integral();
      if (w > 0.) { sumW += w; sumWY += w * loc.Y(); sumWZ += w * loc.Z(); }

      for (art::Ptr<anab::Calorimetry> const& c : calos) {
        if (c->PlaneID().Plane != hinfo.h.plane) continue;
        // found the plane: now find this hit within it
        for (std::size_t iCalo = 0; iCalo < c->dQdx().size(); ++iCalo) {
          if (c->TpIndices()[iCalo] == hits[i].key()) {
            hinfo.oncalo = true;
            hinfo.pitch  = c->TrkPitchVec()[iCalo];
            hinfo.dqdx   = c->dQdx()[iCalo];
            hinfo.rr     = c->ResidualRange()[iCalo];
            break;
          }
        }
        break; // only the first calorimetry object on the plane, as in the skimmer
      }
    }

    info.hits2.push_back(hinfo);
  }

  // all four must be set: the tool skips a cut when the value is negative, so
  // leaving one at its default would apply the cut to a TPC with no hits
  if (sumW > 0.) {
    extra.centroidY = sumWY / sumW;
    extra.centroidZ = sumWZ / sumW;
  }

  info.hit_min_time_p2_tpcE = hitMinTime(info.hits2, true);
  info.hit_max_time_p2_tpcE = hitMaxTime(info.hits2, true);
  info.hit_min_time_p2_tpcW = hitMinTime(info.hits2, false);
  info.hit_max_time_p2_tpcW = hitMaxTime(info.hits2, false);

  return true;

} // icarus::StoppingTrackSelector::buildTrackInfo()


// -----------------------------------------------------------------------------
void icarus::StoppingTrackSelector::matchFlash
  (std::vector<recob::OpFlash> const& flashes, TrackFlashMatch& m) const
{

  // Default ranking is "MinRadius": the barycenter distance between the flash and
  // the track's charge centroid, the same quantity TPCPMTBarycenterMatchProducer
  // ranks on. Time is recorded but does not rank by default.
  bool const hasCentroid =
    (m.centroidY != TrackFlashMatch::NoPos) && (m.centroidZ != TrackFlashMatch::NoPos);
  bool const rankOnRadius = (fMatchMetric == "MinRadius") && hasCentroid;

  if (!hasCentroid) {
    mf::LogWarning(fLogCategory)
      << "track " << m.trackID << " has no charge centroid; matching on |dt| instead";
  }

  double best = std::numeric_limits<double>::max();

  for (std::size_t i = 0; i < flashes.size(); ++i) {

    recob::OpFlash const& flash = flashes[i];

    double const dt = flash.Time() - m.trackT0;    // both [us]
    if (fUseTimeWindow && (dt < fMatchWindowLow || dt > fMatchWindowHigh)) continue;
    if (flash.TotalPE() < fMinFlashPE) continue;
    if (hasCentroid
        && std::abs(flash.ZCenter() - m.centroidZ) > fMaxFlashTrackZ) continue;

    double const radius = hasCentroid
      ? std::hypot(flash.YCenter() - m.centroidY, flash.ZCenter() - m.centroidZ)
      : TrackFlashMatch::NoPos;

    double const metric = rankOnRadius ? radius : std::abs(dt);
    if (metric >= best) continue;

    best        = metric;
    m.flashID   = static_cast<int>(i);
    m.flashTime = flash.Time();
    m.deltaT    = dt;
    m.radius    = radius;
    m.flashPE   = flash.TotalPE();
    m.flashY    = flash.YCenter();
    m.flashZ    = flash.ZCenter();
  }

} // icarus::StoppingTrackSelector::matchFlash()


// -----------------------------------------------------------------------------
std::vector<icarus::TrackFlashMatch> icarus::StoppingTrackSelector::select
  (art::Event const& e, std::size_t iCryo)
{

  std::vector<TrackFlashMatch> matches;

  if (iCryo >= fTrackLabels.size()) {
    throw art::Exception(art::errors::LogicError)
      << "StoppingTrackSelector::select(): cryostat index " << iCryo
      << " out of range (" << fTrackLabels.size() << " configured).\n";
  }

  auto const pfpHandle = e.getHandle<std::vector<recob::PFParticle>>(fPFPLabels[iCryo]);
  if (!pfpHandle.isValid()) {
    mf::LogWarning(fLogCategory)
      << "No recob::PFParticle with label '" << fPFPLabels[iCryo].encode() << "'";
    return matches;
  }

  auto const trackHandle = e.getHandle<std::vector<recob::Track>>(fTrackLabels[iCryo]);
  if (!trackHandle.isValid()) {
    mf::LogWarning(fLogCategory)
      << "No recob::Track with label '" << fTrackLabels[iCryo].encode() << "'";
    return matches;
  }

  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  art::FindManyP<recob::Track> fmTracks(pfpHandle, e, fTrackLabels[iCryo]);
  if (!fmTracks.isValid()) {
    mf::LogWarning(fLogCategory)
      << "No PFParticle -> Track association with label '"
      << fTrackLabels[iCryo].encode() << "'";
    return matches;
  }

  // built once per event and cryostat, not once per track
  EventAssns const assns { e, pfpHandle, trackHandle, fTrackLabels[iCryo],
                           fCaloLabels[iCryo], fPFPT0Labels[iCryo], fCRTHitT0Label };

  // an absent flash collection is not fatal: every track simply stays unmatched
  static std::vector<recob::OpFlash> const noFlashes;
  auto const flashHandle = e.getHandle<std::vector<recob::OpFlash>>(fFlashLabels[iCryo]);
  std::vector<recob::OpFlash> const& flashes =
    flashHandle.isValid() ? *flashHandle : noFlashes;
  if (!flashHandle.isValid()) {
    mf::LogWarning(fLogCategory)
      << "No recob::OpFlash with label '" << fFlashLabels[iCryo].encode() << "'";
  }

  for (art::Ptr<recob::PFParticle> const& pfp : pfps) {

    // the skimmer's own pre-filters (TrackCaloSkimmer_module.cc:316-325)
    if (pfp->PdgCode() == 11 && !fAllowShowerLikePFPs) continue;

    auto const& tracks = fmTracks.at(pfp.key());
    if (tracks.size() != 1) continue;
    art::Ptr<recob::Track> const& track = tracks.at(0);

    sbn::TrackInfo info;
    TrackExtras extra;
    if (!buildTrackInfo(e, assns, iCryo, pfp, track, info, extra)) continue;

    // the actual selection same art tool the calibration ntuples use
    if (!fSelTool->DoSelect(info)) continue;

    TrackFlashMatch m;
    m.trackID       = info.id;
    m.cryostat      = static_cast<unsigned>(info.cryostat);
    m.whichT0       = info.whicht0;
    m.selected      = true;
    m.trackT0       = extra.t0_ns / 1000.;   // ns -> us, once and only here
    m.centroidY     = extra.centroidY;
    m.centroidZ     = extra.centroidZ;
    m.startX        = info.start.x;
    m.startY        = info.start.y;
    m.startZ        = info.start.z;
    m.endX          = info.end.x;
    m.endY          = info.end.y;
    m.endZ          = info.end.z;
    m.dirY          = info.dir.y;
    m.length        = info.length;
    m.medianEnddQdx = medianEnddQdx(info.hits2);

    matchFlash(flashes, m);
    if (fRequireFlashMatch && !m.hasFlash()) continue;

    matches.push_back(m);

  } // for PFParticles

  mf::LogTrace(fLogCategory)
    << e.id() << " cryostat " << iCryo << ": " << matches.size() << " selected tracks";

  return matches;

} // icarus::StoppingTrackSelector::select()
