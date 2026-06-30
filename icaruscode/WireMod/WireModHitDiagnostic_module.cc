// std includes
#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>

// ROOT includes
#include "TTree.h"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace wiremod {

class WireModHitDiagnostic : public art::EDAnalyzer
{
public:
  explicit WireModHitDiagnostic(fhicl::ParameterSet const& pset);
  void analyze(art::Event const& evt) override;

private:
  std::vector<art::InputTag> fBeforeHitLabels;  // dummyGaushit2dTPC*
  std::vector<art::InputTag> fAfterHitLabels;   // gaushit2dTPC*
  std::vector<art::InputTag> fScaleInfoLabels;  // WireModTPC* (produce scale products)
  float fMatchWindowTicks;
  bool  fOnlyMC; // if true, skip overlay hits (is_mc==0) when filling HitMatchTree

  const geo::WireReadoutGeom* fWireReadout =
    &(art::ServiceHandle<geo::WireReadout const>()->Get());

  // HitMatchTree: one entry per before-hit
  TTree*  fHitMatchTree;
  Int_t   tm_run, tm_subrun, tm_event;
  Int_t   tm_channel, tm_cryo, tm_tpc, tm_plane, tm_wire;
  Float_t tm_before_peak_time, tm_before_peak_amp, tm_before_integral, tm_before_rms;
  Int_t   tm_before_start_tick, tm_before_end_tick;
  Float_t tm_before_gof;
  Int_t   tm_before_dof, tm_before_n_multi;
  // MC / overlay flag from WireModifierXXW BackTracker
  Int_t   tm_is_mc;      // 1 = MC hit, 0 = overlay data hit
  // scale info from WireModifierXXW (-999 if overlay / no EDep match)
  Float_t tm_scale_q;
  Float_t tm_scale_sigma;
  Float_t tm_truth_x;     // truth X coordinate used for graph lookup (cm)
  Float_t tm_theta_xw;    // ThetaXW angle used for graph lookup (degrees)
  Float_t tm_truth_dirx;  // truth dxdr
  Float_t tm_truth_diry;  // truth dydr
  Float_t tm_truth_dirz;  // truth dzdr
  Float_t tm_pitch;       // wirePitch / cosG [cm]
  Float_t tm_dqdx_before; // integral_before / pitch [ADC/cm], no lifetime correction
  // match result
  Int_t   tm_is_matched;
  Int_t   tm_n_after_matches;
  Float_t tm_after_peak_time, tm_after_peak_amp, tm_after_integral, tm_after_rms;
  Int_t   tm_after_start_tick, tm_after_end_tick;
  Float_t tm_after_gof;
  Int_t   tm_after_dof;
  Float_t tm_ratio_integral;
  Float_t tm_ratio_rms;

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

  art::ServiceHandle<art::TFileService> tfs;

  fHitMatchTree = tfs->make<TTree>("HitMatchTree",
                                   "WireMod before-after hit comparison");
  fHitMatchTree->Branch("run",               &tm_run,               "run/I");
  fHitMatchTree->Branch("subrun",            &tm_subrun,            "subrun/I");
  fHitMatchTree->Branch("event",             &tm_event,             "event/I");
  fHitMatchTree->Branch("channel",           &tm_channel,           "channel/I");
  fHitMatchTree->Branch("cryo",              &tm_cryo,              "cryo/I");
  fHitMatchTree->Branch("tpc",               &tm_tpc,               "tpc/I");
  fHitMatchTree->Branch("plane",             &tm_plane,             "plane/I");
  fHitMatchTree->Branch("wire",              &tm_wire,              "wire/I");
  fHitMatchTree->Branch("before_peak_time",  &tm_before_peak_time,  "before_peak_time/F");
  fHitMatchTree->Branch("before_peak_amp",   &tm_before_peak_amp,   "before_peak_amp/F");
  fHitMatchTree->Branch("before_integral",   &tm_before_integral,   "before_integral/F");
  fHitMatchTree->Branch("before_rms",        &tm_before_rms,        "before_rms/F");
  fHitMatchTree->Branch("before_start_tick", &tm_before_start_tick, "before_start_tick/I");
  fHitMatchTree->Branch("before_end_tick",   &tm_before_end_tick,   "before_end_tick/I");
  fHitMatchTree->Branch("before_gof",        &tm_before_gof,        "before_gof/F");
  fHitMatchTree->Branch("before_dof",        &tm_before_dof,        "before_dof/I");
  fHitMatchTree->Branch("before_n_multi",    &tm_before_n_multi,    "before_n_multi/I");
  fHitMatchTree->Branch("is_mc",             &tm_is_mc,             "is_mc/I");
  fHitMatchTree->Branch("scale_q",           &tm_scale_q,           "scale_q/F");
  fHitMatchTree->Branch("scale_sigma",       &tm_scale_sigma,       "scale_sigma/F");
  fHitMatchTree->Branch("truth_x",           &tm_truth_x,           "truth_x/F");
  fHitMatchTree->Branch("theta_xw",          &tm_theta_xw,          "theta_xw/F");
  fHitMatchTree->Branch("truth_dirx",        &tm_truth_dirx,        "truth_dirx/F");
  fHitMatchTree->Branch("truth_diry",        &tm_truth_diry,        "truth_diry/F");
  fHitMatchTree->Branch("truth_dirz",        &tm_truth_dirz,        "truth_dirz/F");
  fHitMatchTree->Branch("pitch",             &tm_pitch,             "pitch/F");
  fHitMatchTree->Branch("dqdx_before",       &tm_dqdx_before,       "dqdx_before/F");
  fHitMatchTree->Branch("is_matched",        &tm_is_matched,        "is_matched/I");
  fHitMatchTree->Branch("n_after_matches",   &tm_n_after_matches,   "n_after_matches/I");
  fHitMatchTree->Branch("after_peak_time",   &tm_after_peak_time,   "after_peak_time/F");
  fHitMatchTree->Branch("after_peak_amp",    &tm_after_peak_amp,    "after_peak_amp/F");
  fHitMatchTree->Branch("after_integral",    &tm_after_integral,    "after_integral/F");
  fHitMatchTree->Branch("after_rms",         &tm_after_rms,         "after_rms/F");
  fHitMatchTree->Branch("after_start_tick",  &tm_after_start_tick,  "after_start_tick/I");
  fHitMatchTree->Branch("after_end_tick",    &tm_after_end_tick,    "after_end_tick/I");
  fHitMatchTree->Branch("after_gof",         &tm_after_gof,         "after_gof/F");
  fHitMatchTree->Branch("after_dof",         &tm_after_dof,         "after_dof/I");
  fHitMatchTree->Branch("ratio_integral",    &tm_ratio_integral,    "ratio_integral/F");
  fHitMatchTree->Branch("ratio_rms",         &tm_ratio_rms,         "ratio_rms/F");

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

  // Gather all after-hits once (searched by channel for every before-hit)
  std::vector<recob::Hit const*> afterHits;
  for (auto const& tag : fAfterHitLabels)
  {
    auto h = evt.getHandle<std::vector<recob::Hit>>(tag);
    if (!h.isValid())
    {
      mf::LogWarning("WireModHitDiagnostic")
        << "After-hit collection not found: " << tag;
      continue;
    }
    for (auto const& hit : *h)
      afterHits.push_back(&hit);
  }

  // Build channel -> [after-hit indices] index
  std::unordered_multimap<raw::ChannelID_t, size_t> afterByChannel;
  afterByChannel.reserve(afterHits.size() * 2);
  for (size_t i = 0; i < afterHits.size(); ++i)
    afterByChannel.emplace(afterHits[i]->Channel(), i);

  std::vector<bool> afterMatched(afterHits.size(), false);

  // Process before-hit collections one at a time to keep the hit index
  // aligned with the scale-info parallel vectors from WireModifierXXW.
  for (size_t j = 0; j < fBeforeHitLabels.size(); ++j)
  {
    auto bHandle = evt.getHandle<std::vector<recob::Hit>>(fBeforeHitLabels[j]);
    if (!bHandle.isValid())
    {
      mf::LogWarning("WireModHitDiagnostic")
        << "Before-hit collection not found: " << fBeforeHitLabels[j];
      continue;
    }
    auto const& localHits = *bHandle;

    // Try to load scale-info and isMC parallel vectors for this TPC
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

    if (j < fScaleInfoLabels.size())
    {
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

      // MC flag and scale info (parallel vectors indexed by hit position in this collection)
      tm_is_mc        = (pIsMC       && i < pIsMC->size())       ? (*pIsMC)[i]       : -1;
      tm_scale_q      = (pScaleQ     && i < pScaleQ->size())     ? (*pScaleQ)[i]     : -999.f;
      tm_scale_sigma  = (pScaleSigma && i < pScaleSigma->size()) ? (*pScaleSigma)[i] : -999.f;
      tm_truth_x      = (pTruthX     && i < pTruthX->size())     ? (*pTruthX)[i]     : -999.f;
      tm_theta_xw     = (pThetaXW    && i < pThetaXW->size())    ? (*pThetaXW)[i]    : -999.f;
      tm_truth_dirx   = (pDirX       && i < pDirX->size())       ? (*pDirX)[i]       : -999.f;
      tm_truth_diry   = (pDirY       && i < pDirY->size())       ? (*pDirY)[i]       : -999.f;
      tm_truth_dirz   = (pDirZ       && i < pDirZ->size())       ? (*pDirZ)[i]       : -999.f;
      tm_pitch        = (pPitch      && i < pPitch->size())      ? (*pPitch)[i]      : -999.f;
      tm_dqdx_before  = (pDQdX       && i < pDQdX->size())       ? (*pDQdX)[i]       : -999.f;

      // if OnlyMC is set, skip overlay hits
      if (fOnlyMC && tm_is_mc != 1) continue;

      // Search for matching after-hit on the same channel
      float winLo  = static_cast<float>(bh.StartTick()) - fMatchWindowTicks;
      float winHi  = static_cast<float>(bh.EndTick())   + fMatchWindowTicks;

      size_t bestIdx  = std::numeric_limits<size_t>::max();
      float  bestDist = std::numeric_limits<float>::max();
      int    nMatch   = 0;

      auto range = afterByChannel.equal_range(bh.Channel());
      for (auto it = range.first; it != range.second; ++it)
      {
        float aPeak = afterHits[it->second]->PeakTime();
        if (aPeak < winLo || aPeak > winHi) continue;
        ++nMatch;
        float dist = std::abs(aPeak - bh.PeakTime());
        if (dist < bestDist) { bestDist = dist; bestIdx = it->second; }
      }

      tm_n_after_matches = nMatch;
      tm_is_matched      = (bestIdx != std::numeric_limits<size_t>::max()) ? 1 : 0;

      if (tm_is_matched)
      {
        for (auto it = range.first; it != range.second; ++it)
        {
          float aPeak = afterHits[it->second]->PeakTime();
          if (aPeak >= winLo && aPeak <= winHi)
            afterMatched[it->second] = true;
        }
        auto const* ah   = afterHits[bestIdx];
        tm_after_peak_time  = ah->PeakTime();
        tm_after_peak_amp   = ah->PeakAmplitude();
        tm_after_integral   = ah->Integral();
        tm_after_rms        = ah->RMS();
        tm_after_start_tick = ah->StartTick();
        tm_after_end_tick   = ah->EndTick();
        tm_after_gof        = ah->GoodnessOfFit();
        tm_after_dof        = ah->DegreesOfFreedom();
        tm_ratio_integral   = (bh.Integral() > 0.f)
                                ? ah->Integral() / bh.Integral() : -999.f;
        tm_ratio_rms        = (bh.RMS() > 0.f)
                                ? ah->RMS() / bh.RMS() : -999.f;
      }
      else
      {
        tm_after_peak_time  = -999.f;
        tm_after_peak_amp   = -999.f;
        tm_after_integral   = -999.f;
        tm_after_rms        = -999.f;
        tm_after_start_tick = -999;
        tm_after_end_tick   = -999;
        tm_after_gof        = -999.f;
        tm_after_dof        = -999;
        tm_ratio_integral   = -999.f;
        tm_ratio_rms        = -999.f;
      }

      fHitMatchTree->Fill();
    }
  }

  // NewHitTree: after-hits not claimed by any before-hit
  for (size_t i = 0; i < afterHits.size(); ++i)
  {
    if (afterMatched[i]) continue;

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
}

DEFINE_ART_MODULE(WireModHitDiagnostic)

} // namespace wiremod
