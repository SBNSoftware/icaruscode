////////////////////////////////////////////////////////////////////////
// Class:       VerticalTrackFlashWaveformAna
// Plugin Type: analyzer
// File:        VerticalTrackFlashWaveformAna_module.cc
//
// Select nearly-vertical cosmic muons inside a Z-slice, match them to a
// recob::OpFlash by YZ-barycenter proximity, and dump for every PMT in
// the flash the OpHit info (integral/amplitude/time/pe) together with
// the single raw::OpDetWaveform covering the flash time on that channel.
//
// Output tree feeds the PMT gain calibration on vertical cosmic muons.
//
// mailto:mvicenzi@bnl.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "icaruscode/Timing/IPMTTimingCorrectionService.h"
#include "icaruscode/Timing/PMTTimingCorrections.h"

#include "TTree.h"

#include <cmath>
#include <limits>
#include <map>
#include <vector>

namespace pmtcalib {
  class VerticalTrackFlashWaveformAna;
}

class pmtcalib::VerticalTrackFlashWaveformAna : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Sequence<art::InputTag> TrackLabels{
        Name("TrackLabels"),
        Comment("recob::Track labels, one per cryostat (E, W)")};

    fhicl::Sequence<art::InputTag> OpFlashLabels{
        Name("OpFlashLabels"),
        Comment("recob::OpFlash labels, one per cryostat (E, W) in the same order as TrackLabels")};

    fhicl::Sequence<art::InputTag> PandoraLabels{
        Name("PandoraLabels"),
        Comment("Pandora producer labels (one per cryostat, E/W), used to fetch the"
                " hit->SpacePoint association for the track charge barycenter")};

    fhicl::Atom<art::InputTag> OpDetWaveformLabel{
        Name("OpDetWaveformLabel"),
        Comment("Single raw::OpDetWaveform collection covering both cryostats")};

    fhicl::Atom<art::InputTag> SimPhotonsLabel{
        Name("SimPhotonsLabel"),
        Comment("Single sim::SimPhoton collection covering both cryostats"),
        ""};

    fhicl::Atom<double> MinAbsDirY{
        Name("MinAbsDirY"),
        Comment("Minimum |dir.Y()| for the track to be considered vertical"),
        0.90};

    fhicl::Atom<double> ZMin{
        Name("ZMin"),
        Comment("Minimum Z (cm): both track endpoints must have Z >= ZMin"),
        -500.0};

    fhicl::Atom<double> ZMax{
        Name("ZMax"),
        Comment("Maximum Z (cm): both track endpoints must have Z <= ZMax"),
        500.0};

    fhicl::Atom<double> MinTrackLength{
        Name("MinTrackLength"),
        Comment("Minimum track length in cm"),
        100.0};

    fhicl::Atom<bool> RequireDownwards{
        Name("RequireDownwards"),
        Comment("If true, only keep tracks going downward (dir.Y() < 0 after optional flip)"),
        true};

    fhicl::Atom<double> MaxDeltaX{
        Name("MaxDeltaX"),
        Comment("Max |X_end - X_start| (cm) between track endpoints"),
        20.0};

    fhicl::Atom<double> MaxDeltaZ{
        Name("MaxDeltaZ"),
        Comment("Max |Z_end - Z_start| (cm) between track endpoints"),
        20.0};

    fhicl::Atom<double> BarycenterMaxDist{
        Name("BarycenterMaxDist"),
        Comment("Max YZ distance (cm) between track midpoint and flash barycenter"),
        30.0};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit VerticalTrackFlashWaveformAna(Parameters const& config);

  VerticalTrackFlashWaveformAna(VerticalTrackFlashWaveformAna const&) = delete;
  VerticalTrackFlashWaveformAna(VerticalTrackFlashWaveformAna&&) = delete;
  VerticalTrackFlashWaveformAna& operator=(VerticalTrackFlashWaveformAna const&) = delete;
  VerticalTrackFlashWaveformAna& operator=(VerticalTrackFlashWaveformAna&&) = delete;

  void beginJob() override;
  void analyze(art::Event const& e) override;

private:
  // --- configuration ---
  std::vector<art::InputTag> fTrackLabels;
  std::vector<art::InputTag> fOpFlashLabels;
  std::vector<art::InputTag> fPandoraLabels;
  art::InputTag fOpDetWaveformLabel;
  art::InputTag fSimPhotonsLabel;
  double fMinAbsDirY;
  double fZMin;
  double fZMax;
  double fMinTrackLength;
  bool fRequireDownwards;
  double fMaxDeltaX;
  double fMaxDeltaZ;
  double fBarycenterMaxDist;

  // --- detector timing (optical tick in us) ---
  double fOpticalTickUs = 0.002; // ICARUS PMT digitizer: 500 MHz -> 2 ns

  // --- PMT timing corrections service (DB-provided) ---
  // OpHit/OpFlash times already include these corrections, but raw
  // OpDetWaveform::TimeStamp() does not. We add the same correction to the
  // waveform timestamp when checking coverage of a flash time.
  icarusDB::PMTTimingCorrections const* fTimingCorrections = nullptr;

  // --- output tree ---
  TTree* fTree = nullptr;
  int fNSelTracks = 0;

  // event
  unsigned int fRun = 0;
  unsigned int fEvent = 0;
  int fCryo = -1;

  // track
  double fTrkStartX = 0., fTrkStartY = 0., fTrkStartZ = 0.;
  double fTrkEndX = 0., fTrkEndY = 0., fTrkEndZ = 0.;
  double fTrkDirX = 0., fTrkDirY = 0., fTrkDirZ = 0.;
  double fTrkLength = 0.;
  double fTrkBaryX = 0., fTrkBaryY = 0., fTrkBaryZ = 0.;

  // flash
  double fFlashTime = 0.;
  double fFlashTotalPE = 0.;
  double fFlashY = 0., fFlashZ = 0.;
  double fFlashWidthY = 0., fFlashWidthZ = 0.;
  int fFlashNPMTs = 0;
  double fFlashBarycenterDist = 0.;

  // per-PMT scalars (one tree row per PMT/waveform)
  int fPMTChannel = -1;
  double fPMTOpHitIntegral = 0.;
  double fPMTOpHitAmplitude = 0.;
  double fPMTOpHitPeakTime = 0.;
  double fPMTOpHitStartTime = 0.;
  double fPMTOpHitPE = 0.;
  std::vector<short> fPMTWaveform;
  double fPMTWaveformTimestamp = 0.;
  double fPMTPosX = 0., fPMTPosY = 0., fPMTPosZ = 0.;
  double fPMTBarycenterDist = 0.;

  // simphotons
  int fNSimPhotons = 0;
  std::vector<double> fSimPhotonTimes;
  std::vector<double> fSimPhotonInitX;
  std::vector<double> fSimPhotonInitY;
  std::vector<double> fSimPhotonInitZ;

  // --- helpers (stubs; to be filled in next iteration) ---

  /// Initialise branches on the output tree.
  void setupTree();

  /// Clear all per-entry buffers.
  void resetEntry();

  /// Select a nearly-vertical track inside the Z-slice.
  /// Returns true and fills start/end/dir/length into members if selected.
  bool selectVerticalTrack(recob::Track const& track);

  /// Find the OpFlash in `flashes` whose YZ-barycenter is closest to
  /// (trkY, trkZ); return its index (or -1 if none within BarycenterMaxDist).
  /// On success sets fFlashBarycenterDist.
  int matchFlashByBarycenter(std::vector<art::Ptr<recob::OpFlash>> const& flashes,
                             double trkY, double trkZ);

  /// Total timing correction (us) to apply to a raw waveform TimeStamp so
  /// that it lives in the same reference frame as OpHit/OpFlash times.
  /// Corrections are additive; the sum matches what OpHitTimingCorrection
  /// applies to OpHits.
  double getTimingCorrection(raw::Channel_t channel) const;

  /// Scan the full waveform collection and return the raw::OpDetWaveform on
  /// `channel` whose time window [TimeStamp+corr, TimeStamp+corr+size*tick)
  /// contains `flashTime`. `correction` is the per-channel correction already
  /// extracted once by the caller. Returns nullptr if none found.
  raw::OpDetWaveform const* findCoveringWaveform(
      std::vector<raw::OpDetWaveform> const& wfs,
      raw::Channel_t channel,
      double flashTime,
      double correction) const;

  /// Fill per-PMT vectors for the matched flash by looping the waveform
  /// collection once per PMT channel in the flash.
  void fillFlashPMTs(
      recob::OpFlash const& flash,
      std::vector<art::Ptr<recob::OpHit>> const& ophits,
      std::vector<raw::OpDetWaveform> const& wfs,
      std::vector<sim::SimPhotons> const& simphotonsCollection,
      detinfo::DetectorTimings const& timings);
};

// =====================================================================

pmtcalib::VerticalTrackFlashWaveformAna::VerticalTrackFlashWaveformAna(Parameters const& config)
  : art::EDAnalyzer{config}
  , fTrackLabels{config().TrackLabels()}
  , fOpFlashLabels{config().OpFlashLabels()}
  , fPandoraLabels{config().PandoraLabels()}
  , fOpDetWaveformLabel{config().OpDetWaveformLabel()}
  , fSimPhotonsLabel{config().SimPhotonsLabel()}
  , fMinAbsDirY{config().MinAbsDirY()}
  , fZMin{config().ZMin()}
  , fZMax{config().ZMax()}
  , fMinTrackLength{config().MinTrackLength()}
  , fRequireDownwards{config().RequireDownwards()}
  , fMaxDeltaX{config().MaxDeltaX()}
  , fMaxDeltaZ{config().MaxDeltaZ()}
  , fBarycenterMaxDist{config().BarycenterMaxDist()}
  , fTimingCorrections{lar::providerFrom<icarusDB::IPMTTimingCorrectionService const>()}
{
  if (fTrackLabels.size() != fOpFlashLabels.size()
      || fTrackLabels.size() != fPandoraLabels.size()) {
    throw art::Exception(art::errors::Configuration)
        << "TrackLabels, OpFlashLabels and PandoraLabels must all have the same"
           " size (one per cryostat).";
  }
}

void pmtcalib::VerticalTrackFlashWaveformAna::beginJob()
{
  setupTree();
  fNSelTracks = 0;
}

void pmtcalib::VerticalTrackFlashWaveformAna::analyze(art::Event const& e)
{
  fRun = e.run();
  fEvent = e.event();

  // Refresh the optical tick period from DetectorClocksService each event
  // (cheap and keeps data/MC config consistent).
  auto const clocks =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  fOpticalTickUs = clocks.OpticalClock().TickPeriod();

  // Single collections spanning both cryostats.
  auto const& wfs =
      *e.getValidHandle<std::vector<raw::OpDetWaveform>>(fOpDetWaveformLabel);

  auto const timings = detinfo::makeDetectorTimings(clocks);

  static const std::vector<sim::SimPhotons> emptySimPhotons;
  auto const & simphotonsCollection = (!fSimPhotonsLabel.empty())
                                      ? e.getProduct<std::vector<sim::SimPhotons>>(fSimPhotonsLabel)
                                      : emptySimPhotons;

  for (std::size_t icryo = 0; icryo < fTrackLabels.size(); ++icryo) {
    auto const& trackHandle =
        e.getValidHandle<std::vector<recob::Track>>(fTrackLabels[icryo]);
    auto const& tracks = *trackHandle;
    auto const& flashHandle =
        e.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashLabels[icryo]);
    auto const& flashes = *flashHandle;

    if (flashes.empty()) continue;

    // Promote OpFlash vector to Ptrs (for the matcher signature) and fetch
    // the flash->OpHit associations once per cryostat.
    std::vector<art::Ptr<recob::OpFlash>> flashPtrs;
    flashPtrs.reserve(flashes.size());
    for (std::size_t i = 0; i < flashes.size(); ++i)
      flashPtrs.emplace_back(flashHandle, i);

    art::FindManyP<recob::OpHit> ophitsPtr(flashHandle, e, fOpFlashLabels[icryo]);

    // Track -> Hit association, from the track producer (pattern borrowed
    // from sbncode TrackCaloSkimmer_module).
    art::FindManyP<recob::Hit> fmtrkHits(trackHandle, e, fTrackLabels[icryo]);

    std::vector<art::Ptr<recob::Hit>> const emptyHitVector;

    for (std::size_t itrk = 0; itrk < tracks.size(); ++itrk) {
      resetEntry();
      if (!selectVerticalTrack(tracks[itrk])) continue;

      art::Ptr<recob::Track> const trkPtr(trackHandle, itrk);
      std::vector<art::Ptr<recob::Hit>> const& trkHits =
          fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;

      // Charge-weighted track barycenter using hit Integral as weight and
      // the associated SpacePoint(s) for 3D position — same chain as
      // sbncode TrackCaloSkimmer: Hit -> SpacePoint assns live under the
      // pandora (PFP) producer label.
      art::FindManyP<recob::SpacePoint> fmtrkHitSPs(trkHits, e, fPandoraLabels[icryo]);

      double sumCharge = 0.;
      double sumX = 0., sumY = 0., sumZ = 0.;
      for (std::size_t k = 0; k < trkHits.size(); ++k) {
        auto const& sps = fmtrkHitSPs.at(k);
        if (sps.empty()) continue;
        double const q = trkHits[k]->Integral();
        auto const* xyz = sps.front()->XYZ();
        sumCharge += q;
        sumX += xyz[0] * q;
        sumY += xyz[1] * q;
        sumZ += xyz[2] * q;
      }
      if (sumCharge <= 0.) continue;
      fTrkBaryX = sumX / sumCharge;
      fTrkBaryY = sumY / sumCharge;
      fTrkBaryZ = sumZ / sumCharge;
  
      int const iflash = matchFlashByBarycenter(flashPtrs, fTrkBaryY, fTrkBaryZ);
      if (iflash < 0) continue;

      auto const& flash = flashes[iflash];
      fFlashTime = flash.Time();
      fFlashTotalPE = flash.TotalPE();
      fFlashY = flash.YCenter();
      fFlashZ = flash.ZCenter();
      fFlashWidthY = flash.YWidth();
      fFlashWidthZ = flash.ZWidth();
  
      mf::LogInfo("VerticalTrackFlashWaveformAna") << "Vertical track matched to flash: " 
        << "\n TrackYDir: " << fTrkDirY
        << "\n TrackX: " << fTrkStartX << " " << fTrkBaryX << " " << fTrkEndX
        << "\n TrackY: " << fTrkStartY << " " << fTrkBaryY << " " << fTrkEndY
        << "\n TrackZ: " << fTrkStartZ << " " << fTrkBaryZ << " " << fTrkEndZ
        << "\n TrkLength: " << fTrkLength
        << "\n FlashTime: " << fFlashTime
        << "\n FlashY: " << fFlashY
        << "\n FlashZ: " << fFlashZ 
        << "\n FlashPE: " << fFlashTotalPE;

      fNSelTracks++;
      fCryo = static_cast<int>(icryo);

      auto const& ophits = ophitsPtr.at(iflash);
      fillFlashPMTs(flash, ophits, wfs, simphotonsCollection, timings);

    }
  }
 
  mf::LogInfo("VerticalTrackFlashWaveformAna") << "Total selected tracks: " << fNSelTracks; 
}

// --- helper stubs (bodies to be added in the next step) -----------------

void pmtcalib::VerticalTrackFlashWaveformAna::setupTree()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("vertical", "Vertical muon / flash / PMT waveform tree");

  fTree->Branch("run", &fRun);
  fTree->Branch("event", &fEvent);
  fTree->Branch("cryo", &fCryo);

  fTree->Branch("trk_start_x", &fTrkStartX);
  fTree->Branch("trk_start_y", &fTrkStartY);
  fTree->Branch("trk_start_z", &fTrkStartZ);
  fTree->Branch("trk_end_x", &fTrkEndX);
  fTree->Branch("trk_end_y", &fTrkEndY);
  fTree->Branch("trk_end_z", &fTrkEndZ);
  fTree->Branch("trk_dir_x", &fTrkDirX);
  fTree->Branch("trk_dir_y", &fTrkDirY);
  fTree->Branch("trk_dir_z", &fTrkDirZ);
  fTree->Branch("trk_length", &fTrkLength);
  fTree->Branch("trk_barycenter_x", &fTrkBaryX);
  fTree->Branch("trk_barycenter_y", &fTrkBaryY);
  fTree->Branch("trk_barycenter_z", &fTrkBaryZ);

  fTree->Branch("flash_time", &fFlashTime);
  fTree->Branch("flash_total_pe", &fFlashTotalPE);
  fTree->Branch("flash_y", &fFlashY);
  fTree->Branch("flash_z", &fFlashZ);
  fTree->Branch("flash_width_y", &fFlashWidthY);
  fTree->Branch("flash_width_z", &fFlashWidthZ);
  fTree->Branch("flash_n_pmts", &fFlashNPMTs);
  fTree->Branch("flash_barycenter_dist", &fFlashBarycenterDist);

  fTree->Branch("pmt_channel", &fPMTChannel, "pmt_channel/I");
  fTree->Branch("pmt_ophit_integral", &fPMTOpHitIntegral, "pmt_ophit_integral/D");
  fTree->Branch("pmt_ophit_amplitude", &fPMTOpHitAmplitude, "pmt_ophit_amplitude/D");
  fTree->Branch("pmt_ophit_peak_time", &fPMTOpHitPeakTime, "pmt_ophit_peak_time/D");
  fTree->Branch("pmt_ophit_start_time", &fPMTOpHitStartTime, "pmt_ophit_start_time/D");
  fTree->Branch("pmt_ophit_pe", &fPMTOpHitPE, "pmt_ophit_pe/D");
  fTree->Branch("pmt_waveform", &fPMTWaveform);
  fTree->Branch("pmt_waveform_timestamp", &fPMTWaveformTimestamp, "pmt_waveform_timestamp/D");
  fTree->Branch("pmt_pos_x", &fPMTPosX, "pmt_pos_x/D");
  fTree->Branch("pmt_pos_y", &fPMTPosY, "pmt_pos_y/D");
  fTree->Branch("pmt_pos_z", &fPMTPosZ, "pmt_pos_z/D");
  fTree->Branch("pmt_barycenter_dist", &fPMTBarycenterDist, "pmt_barycenter_dist/D");

  if(!fSimPhotonsLabel.empty())
  {
    fTree->Branch("n_sim_photons", &fNSimPhotons, "n_sim_photons/I");
    fTree->Branch("sim_photon_times", &fSimPhotonTimes);
    fTree->Branch("sim_photon_init_x", &fSimPhotonInitX);
    fTree->Branch("sim_photon_init_y", &fSimPhotonInitY);
    fTree->Branch("sim_photon_init_z", &fSimPhotonInitZ);
  }

}

void pmtcalib::VerticalTrackFlashWaveformAna::resetEntry()
{
  fCryo = -1;
  fTrkStartX = fTrkStartY = fTrkStartZ = 0.;
  fTrkEndX = fTrkEndY = fTrkEndZ = 0.;
  fTrkDirX = fTrkDirY = fTrkDirZ = 0.;
  fTrkLength = 0.;
  fTrkBaryX = fTrkBaryY = fTrkBaryZ = 0.;

  fFlashTime = fFlashTotalPE = 0.;
  fFlashY = fFlashZ = fFlashWidthY = fFlashWidthZ = 0.;
  fFlashNPMTs = 0;
  fFlashBarycenterDist = 0.;

  fPMTChannel = -1;
  fPMTOpHitIntegral = 0.;
  fPMTOpHitAmplitude = 0.;
  fPMTOpHitPeakTime = 0.;
  fPMTOpHitStartTime = 0.;
  fPMTOpHitPE = 0.;
  fPMTWaveform.clear();
  fPMTWaveformTimestamp = 0.;
  fPMTPosX = fPMTPosY = fPMTPosZ = 0.;

  fNSimPhotons = 0;
  fSimPhotonTimes.clear();
  fSimPhotonInitX.clear();
  fSimPhotonInitY.clear();
  fSimPhotonInitZ.clear();
}

bool pmtcalib::VerticalTrackFlashWaveformAna::selectVerticalTrack(recob::Track const& track)
{
  
  if (track.Length() < fMinTrackLength) return false;

  auto const start = track.Start();
  auto const end = track.End();

  if (start.Z() < fZMin || start.Z() > fZMax) return false;
  if (end.Z() < fZMin || end.Z() > fZMax) return false;
  if (std::abs(end.X() - start.X()) > fMaxDeltaX) return false;
  if (std::abs(end.Z() - start.Z()) > fMaxDeltaZ) return false;

  auto dir = track.StartDirection();
  // Flip to downgoing if requested (same idiom as TimeTrackTreeStorageCRT).
  if (fRequireDownwards && dir.Y() > 0.0) dir = -track.EndDirection();
  if (fRequireDownwards && dir.Y() >= 0.0) return false;
  if (std::abs(dir.Y()) < fMinAbsDirY) return false;

  fTrkStartX = start.X(); fTrkStartY = start.Y(); fTrkStartZ = start.Z();
  fTrkEndX   = end.X();   fTrkEndY   = end.Y();   fTrkEndZ   = end.Z();
  fTrkDirX   = dir.X();   fTrkDirY   = dir.Y();   fTrkDirZ   = dir.Z();
  fTrkLength = track.Length();

  return true;
}

int pmtcalib::VerticalTrackFlashWaveformAna::matchFlashByBarycenter(
    std::vector<art::Ptr<recob::OpFlash>> const& flashes,
    double trkY, double trkZ)
{
  int best = -1;
  double bestDist = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < flashes.size(); ++i) {
    double const dy = flashes[i]->YCenter() - trkY;
    double const dz = flashes[i]->ZCenter() - trkZ;
    double const d = std::sqrt(dy * dy + dz * dz);
    if (d < bestDist) { bestDist = d; best = static_cast<int>(i); }
  }
  if (best < 0 || bestDist > fBarycenterMaxDist) return -1;
  fFlashBarycenterDist = bestDist;
  return best;
}

double pmtcalib::VerticalTrackFlashWaveformAna::getTimingCorrection(
    raw::Channel_t channel) const
{
  // Same additive combination applied to OpHits in OpHitTimingCorrection.
  return fTimingCorrections->getLaserCorrections(channel)
       + fTimingCorrections->getCosmicsCorrections(channel);
}

raw::OpDetWaveform const*
pmtcalib::VerticalTrackFlashWaveformAna::findCoveringWaveform(
    std::vector<raw::OpDetWaveform> const& wfs,
    raw::Channel_t channel,
    double flashTime,
    double correction) const
{
  // Linear scan (performance is not a concern). Both OpDetWaveform::TimeStamp()
  // and OpFlash::Time() are in microseconds. The raw timestamp is un-corrected,
  // while flashTime is already corrected, so shift the waveform window by the
  // per-channel correction before the containment check.
  for (auto const& wf : wfs) {
    if (wf.ChannelNumber() != channel) continue;
    double const t0 = wf.TimeStamp() + correction;
    double const t1 = t0 + wf.size() * fOpticalTickUs;
    if (flashTime >= t0 && flashTime < t1) return &wf;
  }
  return nullptr;
}

void pmtcalib::VerticalTrackFlashWaveformAna::fillFlashPMTs(
    recob::OpFlash const& flash,
    std::vector<art::Ptr<recob::OpHit>> const& ophits,
    std::vector<raw::OpDetWaveform> const& wfs,
    std::vector<sim::SimPhotons> const& simphotonsCollection,
    detinfo::DetectorTimings const& timings)
{
  auto const& channelMap = art::ServiceHandle<geo::WireReadout const>()->Get();

  // A flash can have more than one OpHit on the same channel. Following
  // ICARUSFlashAssAna, keep only the earliest-in-time ophit per channel and
  // record its values for that single hit.
  std::map<raw::Channel_t, recob::OpHit const*> firstHit;
  for (auto const& hit : ophits) {
    auto [it, inserted] = firstHit.try_emplace(hit->OpChannel(), hit.get());
    if (!inserted && hit->StartTime() < it->second->StartTime())
      it->second = hit.get();
  }

  fFlashNPMTs = static_cast<int>(firstHit.size());

  for (auto const& [ch, hitPtr] : firstHit) {
    recob::OpHit const& hit = *hitPtr;

    // Extract the per-channel correction once, then use it to shift the
    // raw waveform timestamp into the same frame as the flash time.
    double const correction = getTimingCorrection(ch);

    raw::OpDetWaveform const* wf = findCoveringWaveform(wfs, ch, flash.AbsTime(), correction);
    if (!wf) {
      mf::LogWarning("VerticalTrackFlashWaveformAna")
          << "No raw::OpDetWaveform covering flash AbsTime=" << flash.AbsTime()
          << " us on channel " << ch;
      continue;
    }

    auto const pmtPos = channelMap.OpDetGeoFromOpChannel(ch).GetCenter();

    fPMTChannel = static_cast<int>(ch);
    fPMTOpHitIntegral = hit.Area();
    fPMTOpHitAmplitude = hit.Amplitude();
    fPMTOpHitPeakTime = hit.PeakTime();
    fPMTOpHitStartTime = hit.StartTime();
    fPMTOpHitPE = hit.PE();
    fPMTWaveform = wf->Waveform();
    fPMTWaveformTimestamp = wf->TimeStamp() + getTimingCorrection(fPMTChannel);
    fPMTPosX = pmtPos.X();
    fPMTPosY = pmtPos.Y();
    fPMTPosZ = pmtPos.Z();
    fPMTBarycenterDist = std::hypot( fPMTPosX - fTrkBaryX, fPMTPosY - fTrkBaryY, fPMTPosZ - fTrkBaryZ );

    // SimPhotons belonging to this channel and falling inside the current
    // waveform window. Collection is empty when no SimPhotons product is
    // available (data, or MC without photon propagation).
    fNSimPhotons = 0;
    fSimPhotonTimes.clear();
    fSimPhotonInitX.clear();
    fSimPhotonInitY.clear();
    fSimPhotonInitZ.clear();

    double const wfstart = fPMTWaveformTimestamp;
    double const wfend   = wfstart + wf->size() * fOpticalTickUs;

    for (auto const& simphotons : simphotonsCollection) {

      if (simphotons.OpChannel() != fPMTChannel) continue;

      for (auto const& ph : simphotons) {
        detinfo::timescales::simulation_time const photonTime{ ph.Time };
        double const t = timings.toElectronicsTime(photonTime).value();

        if (t < wfstart || t > wfend) continue;

        fSimPhotonTimes.push_back(t);
        fSimPhotonInitX.push_back(ph.InitialPosition.X());
        fSimPhotonInitY.push_back(ph.InitialPosition.Y());
        fSimPhotonInitZ.push_back(ph.InitialPosition.Z());
      }
    }
    fNSimPhotons = static_cast<int>(fSimPhotonTimes.size());

    fTree->Fill();
  }
}

DEFINE_ART_MODULE(pmtcalib::VerticalTrackFlashWaveformAna)
