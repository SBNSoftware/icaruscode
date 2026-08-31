// std includes
#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>

// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace wiremod {

// WireMod ThetaXW replicated from WireModUtility::ThetaXW so reco-track angles
// are computed with exactly the same convention as truth angles from WireModifierXXW.
inline double ThetaXW(double dxdr, double dydr, double dzdr, double planeAngle)
{
  const double s = std::sin(planeAngle - 0.5 * M_PI);
  const double c = std::cos(planeAngle - 0.5 * M_PI);
  const double cosG = std::abs(dydr * s + dzdr * c);
  return (cosG > 1e-9) ? std::abs(std::atan(dxdr / cosG)) : 0.5 * M_PI;
}

// Map type: channel → (peakTime, plane, reco_thxw_deg)
// Used to carry pre-computed reco angles keyed by the hit they belong to.
using RecoThxwMap = std::unordered_multimap<
                      raw::ChannelID_t,
                      std::tuple<float, int, float>>;

static RecoThxwMap buildRecoThxwMap(art::Event const& evt,
                                    std::vector<art::InputTag> const& tags,
                                    geo::WireReadoutGeom const* wr)
{
  RecoThxwMap m;
  for (auto const& tag : tags) {
    auto h = evt.getHandle<std::vector<recob::Track>>(tag);
    if (!h.isValid()) {
      mf::LogWarning("WireModHitDiagnostic")
        << "Track collection not found: " << tag; continue;
    }
    art::FindManyP<recob::Hit, recob::TrackHitMeta> hfm(h, evt, tag);
    if (!hfm.isValid()) continue;
    for (size_t iT = 0; iT < h->size(); ++iT) {
      auto const& trk   = (*h)[iT];
      auto const& hits  = hfm.at(iT);
      auto const& metas = hfm.data(iT);
      for (size_t iH = 0; iH < hits.size(); ++iH) {
        size_t pt = metas[iH]->Index();
        if (pt == std::numeric_limits<unsigned int>::max() || !trk.HasValidPoint(pt))
          continue;
        auto dir = trk.DirectionAtPoint(pt);
        if (std::abs(dir.X()) + std::abs(dir.Y()) + std::abs(dir.Z()) < 1e-6) continue;
        auto const& wid = hits[iH]->WireID();
        double thetaZ   = wr->Plane(wid.planeID()).ThetaZ();
        float  thxw_deg = static_cast<float>(
                            ThetaXW(dir.X(), dir.Y(), dir.Z(), thetaZ) * TMath::RadToDeg());
        m.emplace(hits[iH]->Channel(),
                  std::make_tuple(hits[iH]->PeakTime(), static_cast<int>(wid.Plane), thxw_deg));
      }
    }
  }
  return m;
}

// Returns -999 if no match within ±matchWindow ticks.
static float lookupRecoThxw(RecoThxwMap const& m,
                             raw::ChannelID_t ch, float peakTime, int plane,
                             float matchWindow)
{
  auto range = m.equal_range(ch);
  float best = -999.f, bestDt = matchWindow + 1.f;
  for (auto it = range.first; it != range.second; ++it) {
    float tp; int pl; float thxw;
    std::tie(tp, pl, thxw) = it->second;
    if (pl != plane) continue;
    float dt = std::abs(tp - peakTime);
    if (dt < bestDt) { bestDt = dt; best = thxw; }
  }
  return best;
}

class WireModHitDiagnostic : public art::EDAnalyzer
{
public:
  explicit WireModHitDiagnostic(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& evt) override;

private:
  std::vector<art::InputTag> fBeforeHitLabels;   // dummyGaushit2dTPC*
  std::vector<art::InputTag> fAfterHitLabels;    // gaushit2dTPC*
  std::vector<art::InputTag> fScaleInfoLabels;   // WireModTPC* (scale products)
  std::vector<art::InputTag> fTrackLabels;       // pandoraTrackGaus* on after-hits
  std::vector<art::InputTag> fDummyTrackLabels;  // pandoraTrackGaus* on before-hits (optional)
  float fMatchWindowTicks;
  bool  fOnlyMC;

  const geo::WireReadoutGeom* fWireReadout =
    &(art::ServiceHandle<geo::WireReadout const>()->Get());

  // ThetaXW 2D: truth vs reco (after-hit tracks) and truth vs dummy-reco (before-hit tracks)
  TH2F* fThxwTrueVsReco            = nullptr;
  TH2F* fThxwTrueVsReco_p[3]      = {nullptr, nullptr, nullptr};
  TH1F* fThxwResid                 = nullptr;
  TH2F* fThxwTrueVsDummyReco       = nullptr;
  TH2F* fThxwTrueVsDummyReco_p[3] = {nullptr, nullptr, nullptr};
  TH1F* fThxwDummyResid            = nullptr;

  // HitMatchTree: one entry per before-hit
  TTree*  fHitMatchTree;
  Int_t   tm_run, tm_subrun, tm_event;
  Int_t   tm_channel, tm_cryo, tm_tpc, tm_plane, tm_wire;
  Float_t tm_before_peak_time, tm_before_peak_amp, tm_before_integral, tm_before_rms;
  Int_t   tm_before_start_tick, tm_before_end_tick;
  Float_t tm_before_gof;
  Int_t   tm_before_dof, tm_before_n_multi;
  Int_t   tm_is_mc;
  Float_t tm_scale_q, tm_scale_sigma;
  Float_t tm_truth_x;
  Float_t tm_theta_xw;              // truth ThetaXW from WireMod [deg]
  Float_t tm_reco_theta_xw;         // reco ThetaXW from after-hit track [deg]
  Float_t tm_dummy_reco_theta_xw;   // reco ThetaXW from before-hit track [deg]
  Float_t tm_truth_dirx, tm_truth_diry, tm_truth_dirz;
  Float_t tm_pitch;
  Float_t tm_dqdx_before;
  Float_t tm_truth_e, tm_truth_nelec;
  Int_t   tm_is_matched, tm_n_after_matches, tm_n_before_matches;
  Float_t tm_after_peak_time, tm_after_peak_amp, tm_after_integral, tm_after_rms;
  Int_t   tm_after_start_tick, tm_after_end_tick;
  Float_t tm_after_gof;
  Int_t   tm_after_dof, tm_after_n_multi;
  Float_t tm_ratio_integral, tm_ratio_rms;

  // NewHitTree: one entry per after-hit without a before-hit match
  TTree*  fNewHitTree;
  Int_t   tn_run, tn_subrun, tn_event;
  Int_t   tn_channel, tn_cryo, tn_tpc, tn_plane, tn_wire;
  Float_t tn_peak_time, tn_peak_amp, tn_integral, tn_rms;
  Int_t   tn_start_tick, tn_end_tick;
  Float_t tn_gof;
  Int_t   tn_dof;
};

//----------------------------------------------------------------------------
WireModHitDiagnostic::WireModHitDiagnostic(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}
  , fMatchWindowTicks{pset.get<float>("MatchWindowTicks", 10.f)}
  , fOnlyMC{pset.get<bool>("OnlyMC", false)}
{
  for (auto const& s : pset.get<std::vector<std::string>>("BeforeHitLabels"))
    fBeforeHitLabels.emplace_back(s);
  for (auto const& s : pset.get<std::vector<std::string>>("AfterHitLabels"))
    fAfterHitLabels.emplace_back(s);
  for (auto const& s : pset.get<std::vector<std::string>>("ScaleInfoLabels",
                                                           std::vector<std::string>{}))
    fScaleInfoLabels.emplace_back(s);
  for (auto const& s : pset.get<std::vector<std::string>>("TrackLabels",
                                                           std::vector<std::string>{}))
    fTrackLabels.emplace_back(s);
  for (auto const& s : pset.get<std::vector<std::string>>("DummyTrackLabels",
                                                           std::vector<std::string>{}))
    fDummyTrackLabels.emplace_back(s);

  art::ServiceHandle<art::TFileService> tfs;

  const int nb = 90; const double lo = 0., hi = 90.;

  // POST-WireMod reco vs truth
  fThxwTrueVsReco = tfs->make<TH2F>("ThxwTrueVsReco",
      "Truth vs POST-WireMod reco #theta_{XW} (all planes);"
      "truth #theta_{XW} [deg];reco #theta_{XW} [deg]",
      nb, lo, hi, nb, lo, hi);
  for (int p = 0; p < 3; ++p)
    fThxwTrueVsReco_p[p] = tfs->make<TH2F>(Form("ThxwTrueVsReco_p%d", p),
        Form("Truth vs POST-WireMod reco #theta_{XW} (plane %d);"
             "truth #theta_{XW} [deg];reco #theta_{XW} [deg]", p),
        nb, lo, hi, nb, lo, hi);
  fThxwResid = tfs->make<TH1F>("ThxwResid",
      "POST-WireMod reco - truth #theta_{XW};"
      "#Delta#theta_{XW} [deg];hits", 200, -45., 45.);

  // PRE-WireMod reco vs truth (only filled when DummyTrackLabels are provided)
  fThxwTrueVsDummyReco = tfs->make<TH2F>("ThxwTrueVsDummyReco",
      "Truth vs PRE-WireMod reco #theta_{XW} (all planes);"
      "truth #theta_{XW} [deg];reco #theta_{XW} [deg]",
      nb, lo, hi, nb, lo, hi);
  for (int p = 0; p < 3; ++p)
    fThxwTrueVsDummyReco_p[p] = tfs->make<TH2F>(Form("ThxwTrueVsDummyReco_p%d", p),
        Form("Truth vs PRE-WireMod reco #theta_{XW} (plane %d);"
             "truth #theta_{XW} [deg];reco #theta_{XW} [deg]", p),
        nb, lo, hi, nb, lo, hi);
  fThxwDummyResid = tfs->make<TH1F>("ThxwDummyResid",
      "PRE-WireMod reco - truth #theta_{XW};"
      "#Delta#theta_{XW} [deg];hits", 200, -45., 45.);

  fHitMatchTree = tfs->make<TTree>("HitMatchTree",
                                   "WireMod before-after hit comparison");
  fHitMatchTree->Branch("run",                    &tm_run,                    "run/I");
  fHitMatchTree->Branch("subrun",                 &tm_subrun,                 "subrun/I");
  fHitMatchTree->Branch("event",                  &tm_event,                  "event/I");
  fHitMatchTree->Branch("channel",                &tm_channel,                "channel/I");
  fHitMatchTree->Branch("cryo",                   &tm_cryo,                   "cryo/I");
  fHitMatchTree->Branch("tpc",                    &tm_tpc,                    "tpc/I");
  fHitMatchTree->Branch("plane",                  &tm_plane,                  "plane/I");
  fHitMatchTree->Branch("wire",                   &tm_wire,                   "wire/I");
  fHitMatchTree->Branch("before_peak_time",       &tm_before_peak_time,       "before_peak_time/F");
  fHitMatchTree->Branch("before_peak_amp",        &tm_before_peak_amp,        "before_peak_amp/F");
  fHitMatchTree->Branch("before_integral",        &tm_before_integral,        "before_integral/F");
  fHitMatchTree->Branch("before_rms",             &tm_before_rms,             "before_rms/F");
  fHitMatchTree->Branch("before_start_tick",      &tm_before_start_tick,      "before_start_tick/I");
  fHitMatchTree->Branch("before_end_tick",        &tm_before_end_tick,        "before_end_tick/I");
  fHitMatchTree->Branch("before_gof",             &tm_before_gof,             "before_gof/F");
  fHitMatchTree->Branch("before_dof",             &tm_before_dof,             "before_dof/I");
  fHitMatchTree->Branch("before_n_multi",         &tm_before_n_multi,         "before_n_multi/I");
  fHitMatchTree->Branch("is_mc",                  &tm_is_mc,                  "is_mc/I");
  fHitMatchTree->Branch("scale_q",                &tm_scale_q,                "scale_q/F");
  fHitMatchTree->Branch("scale_sigma",            &tm_scale_sigma,            "scale_sigma/F");
  fHitMatchTree->Branch("truth_x",                &tm_truth_x,                "truth_x/F");
  // ThetaXW angles: all three perspectives
  fHitMatchTree->Branch("theta_xw",              &tm_theta_xw,              "theta_xw/F");
  fHitMatchTree->Branch("reco_theta_xw",         &tm_reco_theta_xw,         "reco_theta_xw/F");
  fHitMatchTree->Branch("dummy_reco_theta_xw",   &tm_dummy_reco_theta_xw,   "dummy_reco_theta_xw/F");
  fHitMatchTree->Branch("truth_dirx",             &tm_truth_dirx,             "truth_dirx/F");
  fHitMatchTree->Branch("truth_diry",             &tm_truth_diry,             "truth_diry/F");
  fHitMatchTree->Branch("truth_dirz",             &tm_truth_dirz,             "truth_dirz/F");
  fHitMatchTree->Branch("pitch",                  &tm_pitch,                  "pitch/F");
  fHitMatchTree->Branch("dqdx_before",            &tm_dqdx_before,            "dqdx_before/F");
  fHitMatchTree->Branch("truth_e",                &tm_truth_e,                "truth_e/F");
  fHitMatchTree->Branch("truth_nelec",            &tm_truth_nelec,            "truth_nelec/F");
  fHitMatchTree->Branch("is_matched",             &tm_is_matched,             "is_matched/I");
  fHitMatchTree->Branch("n_after_matches",        &tm_n_after_matches,        "n_after_matches/I");
  fHitMatchTree->Branch("n_before_matches",       &tm_n_before_matches,       "n_before_matches/I");
  fHitMatchTree->Branch("after_peak_time",        &tm_after_peak_time,        "after_peak_time/F");
  fHitMatchTree->Branch("after_peak_amp",         &tm_after_peak_amp,         "after_peak_amp/F");
  fHitMatchTree->Branch("after_integral",         &tm_after_integral,         "after_integral/F");
  fHitMatchTree->Branch("after_rms",              &tm_after_rms,              "after_rms/F");
  fHitMatchTree->Branch("after_start_tick",       &tm_after_start_tick,       "after_start_tick/I");
  fHitMatchTree->Branch("after_end_tick",         &tm_after_end_tick,         "after_end_tick/I");
  fHitMatchTree->Branch("after_gof",              &tm_after_gof,              "after_gof/F");
  fHitMatchTree->Branch("after_dof",              &tm_after_dof,              "after_dof/I");
  fHitMatchTree->Branch("after_n_multi",          &tm_after_n_multi,          "after_n_multi/I");
  fHitMatchTree->Branch("ratio_integral",         &tm_ratio_integral,         "ratio_integral/F");
  fHitMatchTree->Branch("ratio_rms",              &tm_ratio_rms,              "ratio_rms/F");

  fNewHitTree = tfs->make<TTree>("NewHitTree",
                                 "After-mod hits with no before-mod counterpart");
  fNewHitTree->Branch("run",        &tn_run,        "run/I");
  fNewHitTree->Branch("subrun",     &tn_subrun,     "subrun/I");
  fNewHitTree->Branch("event",      &tn_event,      "event/I");
  fNewHitTree->Branch("channel",    &tn_channel,    "channel/I");
  fNewHitTree->Branch("cryo",       &tn_cryo,       "cryo/I");
  fNewHitTree->Branch("tpc",        &tn_tpc,        "tpc/I");
  fNewHitTree->Branch("plane",      &tn_plane,      "plane/I");
  fNewHitTree->Branch("wire",       &tn_wire,       "wire/I");
  fNewHitTree->Branch("peak_time",  &tn_peak_time,  "peak_time/F");
  fNewHitTree->Branch("peak_amp",   &tn_peak_amp,   "peak_amp/F");
  fNewHitTree->Branch("integral",   &tn_integral,   "integral/F");
  fNewHitTree->Branch("rms",        &tn_rms,        "rms/F");
  fNewHitTree->Branch("start_tick", &tn_start_tick, "start_tick/I");
  fNewHitTree->Branch("end_tick",   &tn_end_tick,   "end_tick/I");
  fNewHitTree->Branch("gof",        &tn_gof,        "gof/F");
  fNewHitTree->Branch("dof",        &tn_dof,        "dof/I");
}

// global TPC index: 0=EE (cryo0,tpcset0) 1=EW (cryo0,tpcset1)
//                   2=WE (cryo1,tpcset0) 3=WW (cryo1,tpcset1)
int globalTPC(raw::ChannelID_t ch, geo::WireReadoutGeom const* wr)
{
  auto rop = wr->ChannelToROP(ch);
  return static_cast<int>(rop.Cryostat) * 2 + static_cast<int>(rop.TPCset);
}

//----------------------------------------------------------------------------
void WireModHitDiagnostic::analyze(art::Event const& evt)
{
  tm_run    = tn_run    = static_cast<Int_t>(evt.run());
  tm_subrun = tn_subrun = static_cast<Int_t>(evt.subRun());
  tm_event  = tn_event  = static_cast<Int_t>(evt.event());

  RecoThxwMap recoThxwByChannel      = buildRecoThxwMap(evt, fTrackLabels,      fWireReadout);
  RecoThxwMap dummyRecoThxwByChannel = buildRecoThxwMap(evt, fDummyTrackLabels, fWireReadout);

  std::vector<recob::Hit const*> afterHits;
  for (auto const& tag : fAfterHitLabels) {
    auto h = evt.getHandle<std::vector<recob::Hit>>(tag);
    if (!h.isValid()) {
      mf::LogWarning("WireModHitDiagnostic")
        << "After-hit collection not found: " << tag; continue;
    }
    for (auto const& hit : *h) afterHits.push_back(&hit);
  }

  std::unordered_multimap<raw::ChannelID_t, size_t> afterByChannel;
  afterByChannel.reserve(afterHits.size() * 2);
  for (size_t i = 0; i < afterHits.size(); ++i)
    afterByChannel.emplace(afterHits[i]->Channel(), i);

  std::vector<int> afterNBeforeMatches(afterHits.size(), 0);
  RecoThxwMap truthThxwByChannel;

  for (size_t j = 0; j < fBeforeHitLabels.size(); ++j)
  {
    auto bHandle = evt.getHandle<std::vector<recob::Hit>>(fBeforeHitLabels[j]);
    if (!bHandle.isValid()) {
      mf::LogWarning("WireModHitDiagnostic")
        << "Before-hit collection not found: " << fBeforeHitLabels[j]; continue;
    }
    auto const& localHits = *bHandle;

    // scale-info parallel vectors for this TPC (indexed 1:1 with localHits)
    std::vector<float> const* pScaleQ     = nullptr;
    std::vector<float> const* pScaleSigma = nullptr;
    std::vector<float> const* pTruthX     = nullptr;
    std::vector<float> const* pThetaXW    = nullptr;
    std::vector<int>   const* pIsMC       = nullptr;
    std::vector<float> const* pDirX       = nullptr;
    std::vector<float> const* pDirY       = nullptr;
    std::vector<float> const* pDirZ       = nullptr;
    std::vector<float> const* pPitch      = nullptr;
    std::vector<float> const* pDQdX       = nullptr;
    std::vector<float> const* pTruthE     = nullptr;
    std::vector<float> const* pTruthNelec = nullptr;

    if (j < fScaleInfoLabels.size()) {
      auto const& siLabel = fScaleInfoLabels[j];
      auto sqH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "scaleQ"));
      auto ssH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "scaleSigma"));
      auto txH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "truthX"));
      auto thH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "thetaXW"));
      auto imH  = evt.getHandle<std::vector<int>>  (art::InputTag(siLabel.label(), "isMC"));
      auto dxH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "dirX"));
      auto dyH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "dirY"));
      auto dzH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "dirZ"));
      auto ptH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "pitch"));
      auto dqH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "dQdX"));
      auto teH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "truthE"));
      auto tnH  = evt.getHandle<std::vector<float>>(art::InputTag(siLabel.label(), "truthNelec"));
      if (sqH.isValid())  pScaleQ     = &(*sqH);
      if (ssH.isValid())  pScaleSigma = &(*ssH);
      if (txH.isValid())  pTruthX     = &(*txH);
      if (thH.isValid())  pThetaXW    = &(*thH);
      if (imH.isValid())  pIsMC       = &(*imH);
      if (dxH.isValid())  pDirX       = &(*dxH);
      if (dyH.isValid())  pDirY       = &(*dyH);
      if (dzH.isValid())  pDirZ       = &(*dzH);
      if (ptH.isValid())  pPitch      = &(*ptH);
      if (dqH.isValid())  pDQdX       = &(*dqH);
      if (teH.isValid())  pTruthE     = &(*teH);
      if (tnH.isValid())  pTruthNelec = &(*tnH);
    }

    // pre-pass: count before-hits per after-hit so n_before_matches is ready for tree fill.
    // Window on PeakTime (not StartTick/EndTick) to avoid multi-Gaussian group mismatches.
    for (size_t i = 0; i < localHits.size(); ++i) {
      auto const& bh = localHits[i];
      if (fOnlyMC) {
        int isMC = (pIsMC && i < pIsMC->size()) ? (*pIsMC)[i] : -1;
        if (isMC != 1) continue;
      }
      float winLo = bh.PeakTime() - fMatchWindowTicks;
      float winHi = bh.PeakTime() + fMatchWindowTicks;
      auto range = afterByChannel.equal_range(bh.Channel());
      for (auto it = range.first; it != range.second; ++it) {
        if (afterHits[it->second]->PeakTime() >= winLo &&
            afterHits[it->second]->PeakTime() <= winHi)
          ++afterNBeforeMatches[it->second];
      }
    }

    for (size_t i = 0; i < localHits.size(); ++i)
    {
      auto const& bh  = localHits[i];
      auto const& wid = bh.WireID();

      tm_channel           = static_cast<Int_t>(bh.Channel());
      tm_cryo              = static_cast<Int_t>(wid.Cryostat);
      tm_tpc               = globalTPC(bh.Channel(), fWireReadout);
      tm_plane             = static_cast<Int_t>(wid.Plane);
      tm_wire              = static_cast<Int_t>(wid.Wire);
      tm_before_peak_time  = bh.PeakTime();
      tm_before_peak_amp   = bh.PeakAmplitude();
      tm_before_integral   = bh.Integral();
      tm_before_rms        = bh.RMS();
      tm_before_start_tick = bh.StartTick();
      tm_before_end_tick   = bh.EndTick();
      tm_before_gof        = bh.GoodnessOfFit();
      tm_before_dof        = bh.DegreesOfFreedom();
      tm_before_n_multi    = bh.Multiplicity();

      tm_is_mc             = (pIsMC       && i < pIsMC->size())       ? (*pIsMC)[i]       : -1;
      tm_scale_q           = (pScaleQ     && i < pScaleQ->size())     ? (*pScaleQ)[i]     : -999.f;
      tm_scale_sigma       = (pScaleSigma && i < pScaleSigma->size()) ? (*pScaleSigma)[i] : -999.f;
      tm_truth_x           = (pTruthX     && i < pTruthX->size())     ? (*pTruthX)[i]     : -999.f;
      tm_theta_xw          = (pThetaXW    && i < pThetaXW->size())    ? (*pThetaXW)[i]    : -999.f;
      tm_truth_dirx        = (pDirX       && i < pDirX->size())       ? (*pDirX)[i]       : -999.f;
      tm_truth_diry        = (pDirY       && i < pDirY->size())       ? (*pDirY)[i]       : -999.f;
      tm_truth_dirz        = (pDirZ       && i < pDirZ->size())       ? (*pDirZ)[i]       : -999.f;
      tm_pitch             = (pPitch      && i < pPitch->size())       ? (*pPitch)[i]      : -999.f;
      tm_dqdx_before       = (pDQdX       && i < pDQdX->size())       ? (*pDQdX)[i]       : -999.f;
      tm_truth_e           = (pTruthE     && i < pTruthE->size())     ? (*pTruthE)[i]     : -1.f;
      tm_truth_nelec       = (pTruthNelec && i < pTruthNelec->size()) ? (*pTruthNelec)[i] : -1.f;

      tm_dummy_reco_theta_xw = lookupRecoThxw(dummyRecoThxwByChannel,
                                               bh.Channel(), bh.PeakTime(), tm_plane,
                                               fMatchWindowTicks);

      if (tm_theta_xw > -900.f)
        truthThxwByChannel.emplace(bh.Channel(),
          std::make_tuple(bh.PeakTime(), tm_plane, tm_theta_xw));

      if (fOnlyMC && tm_is_mc != 1) continue;

      float winLo = bh.PeakTime() - fMatchWindowTicks;
      float winHi = bh.PeakTime() + fMatchWindowTicks;

      size_t bestIdx  = std::numeric_limits<size_t>::max();
      float  bestDist = std::numeric_limits<float>::max();
      int    nMatch   = 0;

      auto range = afterByChannel.equal_range(bh.Channel());
      for (auto it = range.first; it != range.second; ++it) {
        float aPeak = afterHits[it->second]->PeakTime();
        if (aPeak < winLo || aPeak > winHi) continue;
        ++nMatch;
        float dist = std::abs(aPeak - bh.PeakTime());
        if (dist < bestDist) { bestDist = dist; bestIdx = it->second; }
      }

      tm_n_after_matches  = nMatch;
      tm_is_matched       = (bestIdx != std::numeric_limits<size_t>::max()) ? 1 : 0;
      tm_n_before_matches = (bestIdx != std::numeric_limits<size_t>::max())
                              ? afterNBeforeMatches[bestIdx] : 0;

      if (tm_is_matched) {
        auto const* ah  = afterHits[bestIdx];
        tm_after_peak_time  = ah->PeakTime();
        tm_after_peak_amp   = ah->PeakAmplitude();
        tm_after_integral   = ah->Integral();
        tm_after_rms        = ah->RMS();
        tm_after_start_tick = ah->StartTick();
        tm_after_end_tick   = ah->EndTick();
        tm_after_gof        = ah->GoodnessOfFit();
        tm_after_dof        = ah->DegreesOfFreedom();
        tm_after_n_multi    = ah->Multiplicity();
        tm_ratio_integral   = (bh.Integral() > 0.f)
                                ? ah->Integral() / bh.Integral() : -999.f;
        tm_ratio_rms        = (bh.RMS() > 0.f)
                                ? ah->RMS() / bh.RMS() : -999.f;

        tm_reco_theta_xw = lookupRecoThxw(recoThxwByChannel,
                                           ah->Channel(), ah->PeakTime(), tm_plane,
                                           fMatchWindowTicks);
      } else {
        tm_after_peak_time  = -999.f; tm_after_peak_amp  = -999.f;
        tm_after_integral   = -999.f; tm_after_rms       = -999.f;
        tm_after_start_tick = -999;   tm_after_end_tick   = -999;
        tm_after_gof        = -999.f; tm_after_dof        = -999;
        tm_after_n_multi    = -999;
        tm_ratio_integral   = -999.f; tm_ratio_rms        = -999.f;
        tm_reco_theta_xw    = -999.f;
      }

      fHitMatchTree->Fill();
    }
  }

  for (size_t i = 0; i < afterHits.size(); ++i) {
    if (afterNBeforeMatches[i] > 0) continue;
    auto const* ah  = afterHits[i];
    auto const& wid = ah->WireID();
    tn_channel    = static_cast<Int_t>(ah->Channel());
    tn_cryo       = static_cast<Int_t>(wid.Cryostat);
    tn_tpc        = globalTPC(ah->Channel(), fWireReadout);
    tn_plane      = static_cast<Int_t>(wid.Plane);
    tn_wire       = static_cast<Int_t>(wid.Wire);
    tn_peak_time  = ah->PeakTime();
    tn_peak_amp   = ah->PeakAmplitude();
    tn_integral   = ah->Integral();
    tn_rms        = ah->RMS();
    tn_start_tick = ah->StartTick();
    tn_end_tick   = ah->EndTick();
    tn_gof        = ah->GoodnessOfFit();
    tn_dof        = ah->DegreesOfFreedom();
    fNewHitTree->Fill();
  }

  for (auto const& tag : fTrackLabels) {
    auto h = evt.getHandle<std::vector<recob::Track>>(tag);
    if (!h.isValid()) continue;
    art::FindManyP<recob::Hit, recob::TrackHitMeta> hfm(h, evt, tag);
    if (!hfm.isValid()) continue;
    for (size_t iT = 0; iT < h->size(); ++iT) {
      auto const& trk   = (*h)[iT];
      auto const& hits  = hfm.at(iT);
      auto const& metas = hfm.data(iT);
      for (size_t iH = 0; iH < hits.size(); ++iH) {
        size_t pt = metas[iH]->Index();
        if (pt == std::numeric_limits<unsigned int>::max() || !trk.HasValidPoint(pt)) continue;
        auto dir = trk.DirectionAtPoint(pt);
        if (std::abs(dir.X()) + std::abs(dir.Y()) + std::abs(dir.Z()) < 1e-6) continue;
        auto const& wid    = hits[iH]->WireID();
        int         plane  = static_cast<int>(wid.Plane);
        double thetaZ      = fWireReadout->Plane(wid.planeID()).ThetaZ();
        float recoThxw     = static_cast<float>(
                               ThetaXW(dir.X(), dir.Y(), dir.Z(), thetaZ) * TMath::RadToDeg());
        float trueThxw = lookupRecoThxw(truthThxwByChannel,
                                         hits[iH]->Channel(), hits[iH]->PeakTime(), plane,
                                         fMatchWindowTicks);
        if (trueThxw < 0.f) continue;
        fThxwTrueVsReco->Fill(trueThxw, recoThxw);
        if (plane >= 0 && plane < 3) fThxwTrueVsReco_p[plane]->Fill(trueThxw, recoThxw);
        fThxwResid->Fill(recoThxw - trueThxw);
      }
    }
  }

  for (auto const& tag : fDummyTrackLabels) {
    auto h = evt.getHandle<std::vector<recob::Track>>(tag);
    if (!h.isValid()) continue;
    art::FindManyP<recob::Hit, recob::TrackHitMeta> hfm(h, evt, tag);
    if (!hfm.isValid()) continue;
    for (size_t iT = 0; iT < h->size(); ++iT) {
      auto const& trk   = (*h)[iT];
      auto const& hits  = hfm.at(iT);
      auto const& metas = hfm.data(iT);
      for (size_t iH = 0; iH < hits.size(); ++iH) {
        size_t pt = metas[iH]->Index();
        if (pt == std::numeric_limits<unsigned int>::max() || !trk.HasValidPoint(pt)) continue;
        auto dir = trk.DirectionAtPoint(pt);
        if (std::abs(dir.X()) + std::abs(dir.Y()) + std::abs(dir.Z()) < 1e-6) continue;
        auto const& wid    = hits[iH]->WireID();
        int         plane  = static_cast<int>(wid.Plane);
        double thetaZ      = fWireReadout->Plane(wid.planeID()).ThetaZ();
        float dummyRecoThxw = static_cast<float>(
                                ThetaXW(dir.X(), dir.Y(), dir.Z(), thetaZ) * TMath::RadToDeg());
        float trueThxw = lookupRecoThxw(truthThxwByChannel,
                                         hits[iH]->Channel(), hits[iH]->PeakTime(), plane,
                                         fMatchWindowTicks);
        if (trueThxw < 0.f) continue;
        fThxwTrueVsDummyReco->Fill(trueThxw, dummyRecoThxw);
        if (plane >= 0 && plane < 3) fThxwTrueVsDummyReco_p[plane]->Fill(trueThxw, dummyRecoThxw);
        fThxwDummyResid->Fill(dummyRecoThxw - trueThxw);
      }
    }
  }
}

DEFINE_ART_MODULE(WireModHitDiagnostic)

} // namespace wiremod
