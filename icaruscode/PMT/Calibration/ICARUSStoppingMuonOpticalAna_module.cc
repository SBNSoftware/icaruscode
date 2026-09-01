////////////////////////////////////////////////////////////////////////
// Class:       ICARUSStoppingMuonOpticalAna
// Plugin Type: analyzer
// File:        ICARUSStoppingMuonOpticalAna_module.cc
//
// Optical dump written only for events containing a stopping muon matched to a
// flash. Flashes and OpHits are then kept for the whole readout; only the
// waveform samples are restricted to the fragments covering a matched flash.
//
// The stopping-track selection is not implemented here: it is delegated to
// icarus::StoppingTrackSelector, which runs the same art tool that is used
// in the calibration ntuples.
//
// Branch names and types are deliberately identical to ICARUSOpticalDebug
// wherever the two write the same quantity, so that a notebook or a comparison
// script works unchanged against either output.
//
// mailto: mvicenzi@bnl.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/DelegatedParameter.h"

#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larana/OpticalDetector/IPedAlgoMakerTool.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "icaruscode/Timing/PMTTimingCorrections.h"
#include "icaruscode/Timing/IPMTTimingCorrectionService.h"
#include "icaruscode/PMT/Calibration/TrackTools/StoppingTrackSelector.h"

#include "TTree.h"

#include <algorithm>
#include <array>
#include <cmath>         // std::abs(double)
#include <cstddef>
#include <cstdlib>       // std::abs(int)
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>       // std::accumulate
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace icarus { class ICARUSStoppingMuonOpticalAna; }

class icarus::ICARUSStoppingMuonOpticalAna : public art::EDAnalyzer {

public:

  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    // --- optical inputs: same keys as ICARUSOpticalDebug ---------------------

    fhicl::Sequence<art::InputTag> OpHitLabels{
        Name("OpHitLabels"),
        Comment("Tags for the recob::OpHit data products")};

    fhicl::Sequence<art::InputTag> FlashLabels{
        Name("FlashLabels"),
        Comment("Tags for the recob::OpFlash data products, one per cryostat.")};

    fhicl::Atom<bool> SaveWaveforms{
        Name("SaveWaveforms"),
        Comment("Set to save the raw::OpDetWaveforms")};

    fhicl::Sequence<art::InputTag> OpDetWaveformLabels{
        Name("OpDetWaveformLabels"),
        Comment("Tags for the raw::OpDetWaveform data products")};

    fhicl::Atom<float> OpHitThresholdADC{
        Name("OpHitThresholdADC"),
        Comment("Threshold in ADC for an OpHit to be considered")};

    fhicl::Atom<art::InputTag> MCParticleTruthLabel{
        Name("MCParticleTruthLabel"),
        Comment("Tag for simb::MCParticle truth info; empty on data"),
        art::InputTag{}};

    fhicl::Atom<art::InputTag> TriggerLabel{
        Name("TriggerLabel"),
        Comment("Unshifted trigger information.")};

    fhicl::Sequence<art::InputTag> SimPhotonsLabels{
        Name("SimPhotonsLabels"),
        Comment("sim::SimPhotons collections to explore")};

    fhicl::DelegatedParameter PedAlgoPset{
        Name("PedAlgoRollingMeanMaker"),
        Comment("parameters of the pedestal extraction algorithm.")};

    // --- selection --------------------------------------------

    fhicl::DelegatedParameter Selector{
        Name("Selector"),
        Comment("configuration of icarus::StoppingTrackSelector")};

    fhicl::Atom<bool> SkipUnmatchedEvents{
        Name("SkipUnmatchedEvents"),
        Comment("Write no optical tree for an event in which no stopping track was matched to a flash."),
        true};

  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit ICARUSStoppingMuonOpticalAna(Parameters const& config);

  ICARUSStoppingMuonOpticalAna(ICARUSStoppingMuonOpticalAna const&) = delete;
  ICARUSStoppingMuonOpticalAna(ICARUSStoppingMuonOpticalAna&&) = delete;
  ICARUSStoppingMuonOpticalAna& operator=(ICARUSStoppingMuonOpticalAna const&) = delete;
  ICARUSStoppingMuonOpticalAna& operator=(ICARUSStoppingMuonOpticalAna&&) = delete;

  void beginJob() override;
  void analyze(art::Event const& e) override;

private:

  /// Matched flash times [us, absolute], indexed by cryostat (East=0, West=1),
  /// which is also the index into FlashLabels. One entry per matched flash.
  using FlashTimes = std::vector<std::vector<double>>;

  // --- helpers, thin wrappers over the geometry services --------------------
  geo::CryostatID::CryostatID_t getCryostatByChannel(int channel) const;
  int getSideByChannel(int channel) const;
  std::array<double, 3> getChannelXYZ(int channel) const;
  double getTimingCorrection(int channel) const;

  // --- per-event steps ------------------------------------------------------
  void fillTrackMatchTree(std::vector<std::vector<TrackFlashMatch>> const& matches);
  /// Writes every flash and returns the times of the matched ones.
  FlashTimes fillFlashes(art::Event const& e,
                         std::vector<std::vector<TrackFlashMatch>> const& matches);
  /// Writes one row per fragment; samples only for those covering a matched flash.
  void fillWaveforms(art::Event const& e, FlashTimes const& flashTimes);
  /// Writes every out-of-flash OpHit: they are the cheap record of all the light.
  void fillUnmatchedOpHits(art::Event const& e);
  void fillMCTruth(art::Event const& e);

  /// Fills the in-flash OpHit vectors.
  void fillFlashOpHits(std::vector<art::Ptr<recob::OpHit>> const& ophits);

  // --- configuration --------------------------------------------------------
  std::vector<art::InputTag> fOpHitLabels;
  std::vector<art::InputTag> fFlashLabels;
  std::vector<art::InputTag> fOpDetWaveformLabels;
  art::InputTag fMCParticleLabel;
  art::InputTag fTriggerLabel;
  std::vector<art::InputTag> fSimPhotonsLabels;
  bool const fSaveWaveforms;
  float const fOpHitThresholdADC;

  bool const fSkipUnmatchedEvents;
  double const fOpticalTickPeriod;

  // --- services and state ---------------------------------------------------
  geo::GeometryCore const* fGeom;
  geo::WireReadoutGeom const* fChannelMapAlg;
  icarusDB::PMTTimingCorrections const& fPMTTimingCorrectionsService;
  std::unique_ptr<pmtana::PMTPedestalBase> fPedAlgo;
  StoppingTrackSelector fSelector;

  // --- trees ----------------------------------------------------------------
  TTree* fOpFlashTree = nullptr;
  std::vector<TTree*> fOpHitTrees;
  std::vector<TTree*> fOpDetWaveformTrees;
  TTree* fMCParticleTree = nullptr;
  TTree* fTrackMatchTree = nullptr;

  // --- branch buffers -------------------------------------------------------
  // common
  int m_run;
  int m_subrun;
  int m_event;
  int m_timestamp;

  // flash tree (names and types as in ICARUSOpticalDebug)
  int m_flash_id;
  int m_flash_cryo;
  int m_multiplicity;
  int m_multiplicity_left;
  int m_multiplicity_right;
  float m_sum_pe;
  float m_flash_time;
  float m_flash_y;
  float m_flash_z;
  float m_flash_width_y;
  float m_flash_width_z;
  int m_flash_nhits;
  bool m_track_matched;
  std::vector<int> m_channel_id;
  std::vector<float> m_hit_x, m_hit_y, m_hit_z;
  std::vector<float>  m_hit_start_time, m_hit_rise_time, m_hit_peak_time, m_hit_width;
  std::vector<float> m_hit_timing_corr, m_hit_area, m_hit_pe, m_hit_amplitude;

  // out-of-flash OpHit tree
  int m_channel;
  float m_x, m_y, m_z;
  float m_integral;
  float m_amplitude;
  float m_start_time;
  float m_peak_time;
  float m_rise_time;
  float m_timing_corr;
  float m_width;
  float m_abs_start_time;
  float m_pe;
  float m_fast_to_total;

  // waveform tree
  int m_wfchannel;
  float m_wfstart;           ///< raw, uncorrected TimeStamp() [us]
  float m_wf_timing_corr;    ///< per-channel PMT timing correction [us]
  bool  m_wf_saved;          ///< whether `wf`/`bs` hold the samples for this row
  unsigned int m_nticks;
  std::vector<short> m_wf;
  std::vector<float> m_bs;

  // truth tree
  int track_gen_pdg;
  float track_gen_time, track_gen_E;
  /// simb::MCParticle::EndT() of the primary muon [ns]. Beware: this is the
  /// stop for mu- (it matches the atomic-cascade daughters) but the decay
  /// for mu+ (it matches the decay daughters).
  float track_end_time;
  float time_shift; //< Per-event offset applied to the waveforms and the SimPhotons [ns].
  float track_gen_x, track_gen_y, track_gen_z;
  int michel_gen_pdg;
  int michel_process;   ///< 0 free decay, 1 atomic cascade, 2 bound decay, -1 none
  float michel_gen_time, michel_gen_E;
  float michel_gen_x, michel_gen_y, michel_gen_z;
  bool michel_in_av;
  bool track_end_in_av;

  // track-match tree
  int t_track_id;
  int t_cryo;
  int t_whicht0;
  int t_flash_id;
  float t_t0;
  float t_start_x, t_start_y, t_start_z;
  float t_end_x, t_end_y, t_end_z;
  float t_dir_y;
  float t_length;
  float t_median_end_dqdx;
  float t_centroid_y, t_centroid_z;
  float t_flash_time, t_flash_pe, t_flash_y, t_flash_z;
  float t_dt;
  float t_radius;
  float t_trange_low, t_trange_high;   ///< drift-allowed interval [us], no margin

}; // class icarus::ICARUSStoppingMuonOpticalAna


// -----------------------------------------------------------------------------
icarus::ICARUSStoppingMuonOpticalAna::ICARUSStoppingMuonOpticalAna
  (Parameters const& config)
  : EDAnalyzer(config)
  , fOpHitLabels(config().OpHitLabels())
  , fFlashLabels(config().FlashLabels())
  , fOpDetWaveformLabels(config().OpDetWaveformLabels())
  , fMCParticleLabel(config().MCParticleTruthLabel())
  , fTriggerLabel(config().TriggerLabel())
  , fSimPhotonsLabels(config().SimPhotonsLabels())
  , fSaveWaveforms(config().SaveWaveforms())
  , fOpHitThresholdADC(config().OpHitThresholdADC())
  , fSkipUnmatchedEvents(config().SkipUnmatchedEvents())
  , fOpticalTickPeriod(
      art::ServiceHandle<detinfo::DetectorClocksService const>()
        ->DataForJob().OpticalClock().TickPeriod())
  , fGeom(lar::providerFrom<geo::Geometry>())
  , fChannelMapAlg(&art::ServiceHandle<geo::WireReadout const>()->Get())
  , fPMTTimingCorrectionsService(
      *(lar::providerFrom<icarusDB::IPMTTimingCorrectionService const>()))
  , fPedAlgo(art::make_tool<opdet::IPedAlgoMakerTool>(
      config().PedAlgoPset.get<fhicl::ParameterSet>())->makeAlgo())
  , fSelector(config().Selector.get<fhicl::ParameterSet>())
{

  // The selector indexes TrackFlashMatch from the flash product it was
  // given, and this module indexes flash_id into the one it was given. 
  // Those have to be the same products:
  //
  //  flash_labels: [ "opflashCryoE", "opflashCryoW" ]
  //   ...
  //   FlashLabels:          @local::flash_labels
  //   Selector.FlashLabels: @local::flash_labels
  //
  // This check makes a divergence a job-start error
  auto const selectorFlashLabels =
    config().Selector.get<fhicl::ParameterSet>()
      .get<std::vector<art::InputTag>>("FlashLabels", {});

  if (selectorFlashLabels != fFlashLabels) {
    throw art::Exception(art::errors::Configuration)
      << "ICARUSStoppingMuonOpticalAna: FlashLabels and Selector.FlashLabels must"
         " be the same list; point both at one @local:: definition.\n";
  }

} // constructor


// -----------------------------------------------------------------------------
// geometry helpers: thin wrappers over the services, as in ICARUSOpticalDebug
// -----------------------------------------------------------------------------
geo::CryostatID::CryostatID_t
icarus::ICARUSStoppingMuonOpticalAna::getCryostatByChannel(int channel) const
{
  return fChannelMapAlg->OpDetGeoFromOpChannel(channel).ID().Cryostat;
}


int icarus::ICARUSStoppingMuonOpticalAna::getSideByChannel(int channel) const
{
  // channels run east to west, north to south: [0:89] and [180:269] are the east
  // wall of their cryostat (0), [90:179] and [270:359] the west wall (1)
  return (channel / 90) % 2;
}


std::array<double, 3>
icarus::ICARUSStoppingMuonOpticalAna::getChannelXYZ(int channel) const
{
  auto const PMTxyz = fChannelMapAlg->OpDetGeoFromOpChannel(channel).GetCenter();
  return std::array<double, 3>{ PMTxyz.X(), PMTxyz.Y(), PMTxyz.Z() };
}


double icarus::ICARUSStoppingMuonOpticalAna::getTimingCorrection(int channel) const
{
  return fPMTTimingCorrectionsService.getLaserCorrections(channel)
       + fPMTTimingCorrectionsService.getCosmicsCorrections(channel);
}


// -----------------------------------------------------------------------------
void icarus::ICARUSStoppingMuonOpticalAna::beginJob()
{
  art::ServiceHandle<art::TFileService const> tfs;

  // --- out-of-flash OpHits, one tree per label ------------------------------
  for (art::InputTag const& label : fOpHitLabels) {
    std::string const name = label.label() + "_ttree";
    std::string const info = "Out-of-flash recob::OpHit with label " + label.label();
    TTree* ttree = tfs->make<TTree>(name.c_str(), info.c_str());
    ttree->Branch("run", &m_run, "run/I");
    ttree->Branch("subrun", &m_subrun, "subrun/I");
    ttree->Branch("event", &m_event, "event/I");
    ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
    ttree->Branch("channel_id", &m_channel, "channel_id/I");
    ttree->Branch("integral", &m_integral, "integral/F");
    ttree->Branch("amplitude", &m_amplitude, "amplitude/F");
    ttree->Branch("start_time", &m_start_time, "start_time/F");
    ttree->Branch("peak_time", &m_peak_time, "peak_time/F");
    ttree->Branch("rise_time", &m_rise_time, "rise_time/F");
    ttree->Branch("abs_start_time", &m_abs_start_time, "abs_start_time/F");
    ttree->Branch("timing_corr", &m_timing_corr, "timing_corr/F");
    ttree->Branch("pe", &m_pe, "pe/F");
    ttree->Branch("width", &m_width, "width/F");
    ttree->Branch("x", &m_x, "x/F");
    ttree->Branch("y", &m_y, "y/F");
    ttree->Branch("z", &m_z, "z/F");
    ttree->Branch("fast_to_total", &m_fast_to_total, "fast_to_total/F");
    fOpHitTrees.push_back(ttree);
  }

  // --- flashes and their OpHits ---------------------------------------------
  // One tree, not one per cryostat: `flash_cryo` already distinguishes them, and
  // splitting only forces the analysis to read and concatenate two trees.
  {
    fOpFlashTree = tfs->make<TTree>("flash_tree", "recob::OpFlash and their OpHits");
    TTree* ttree = fOpFlashTree;
    ttree->Branch("run", &m_run, "run/I");
    ttree->Branch("subrun", &m_subrun, "subrun/I");
    ttree->Branch("event", &m_event, "event/I");
    ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
    ttree->Branch("flash_id", &m_flash_id, "flash_id/I");
    ttree->Branch("flash_cryo", &m_flash_cryo, "flash_cryo/I");
    ttree->Branch("multiplicity", &m_multiplicity, "multiplicity/I");
    ttree->Branch("multiplicity_right", &m_multiplicity_right, "multiplicity_right/I");
    ttree->Branch("multiplicity_left", &m_multiplicity_left, "multiplicity_left/I");
    ttree->Branch("sum_pe", &m_sum_pe, "sum_pe/F");
    ttree->Branch("flash_time", &m_flash_time, "flash_time/F");
    ttree->Branch("flash_y", &m_flash_y, "flash_y/F");
    ttree->Branch("flash_z", &m_flash_z, "flash_z/F");
    ttree->Branch("flash_width_y", &m_flash_width_y, "flash_width_y/F");
    ttree->Branch("flash_width_z", &m_flash_width_z, "flash_width_z/F");
    ttree->Branch("flash_nhits", &m_flash_nhits, "flash_nhits/I");
    ttree->Branch("track_matched", &m_track_matched, "track_matched/O");
    ttree->Branch("channel_id", &m_channel_id);
    ttree->Branch("hit_x", &m_hit_x);
    ttree->Branch("hit_y", &m_hit_y);
    ttree->Branch("hit_z", &m_hit_z);
    ttree->Branch("hit_start_time", &m_hit_start_time);
    ttree->Branch("hit_rise_time", &m_hit_rise_time);
    ttree->Branch("hit_peak_time", &m_hit_peak_time);
    ttree->Branch("hit_width", &m_hit_width);
    ttree->Branch("hit_timing_corr", &m_hit_timing_corr);
    ttree->Branch("hit_area", &m_hit_area);
    ttree->Branch("hit_pe", &m_hit_pe);
    ttree->Branch("hit_amplitude", &m_hit_amplitude);
  }

  // --- waveforms, one tree per label ----------------------------------------
  if (fSaveWaveforms) {
    for (art::InputTag const& label : fOpDetWaveformLabels) {
      std::string const name = label.label() + label.instance() + "_wfttree";
      std::string const info = "raw::OpDetWaveform with label " + label.label();
      TTree* ttree = tfs->make<TTree>(name.c_str(), info.c_str());
      ttree->Branch("run", &m_run, "run/I");
      ttree->Branch("subrun", &m_subrun, "subrun/I");
      ttree->Branch("event", &m_event, "event/I");
      ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
      ttree->Branch("channel_id", &m_wfchannel, "channel_id/I");
      ttree->Branch("nticks", &m_nticks, "nticks/i");
      ttree->Branch("wf_start", &m_wfstart, "wf_start/F");
      ttree->Branch("timing_corr", &m_wf_timing_corr, "timing_corr/F");
      ttree->Branch("saved", &m_wf_saved, "saved/O"); // waveform saved?
      ttree->Branch("bs", &m_bs);
      ttree->Branch("wf", &m_wf);
      fOpDetWaveformTrees.push_back(ttree);
    }
  }

  // --- one row per selected track -------------------------------------------
  fTrackMatchTree = tfs->make<TTree>("trackmatch_tree",
    "stopping tracks flash-matched");
  fTrackMatchTree->Branch("run", &m_run, "run/I");
  fTrackMatchTree->Branch("subrun", &m_subrun, "subrun/I");
  fTrackMatchTree->Branch("event", &m_event, "event/I");
  fTrackMatchTree->Branch("timestamp", &m_timestamp, "timestamp/I");
  fTrackMatchTree->Branch("track_id", &t_track_id, "track_id/I");
  fTrackMatchTree->Branch("cryo", &t_cryo, "cryo/I");
  fTrackMatchTree->Branch("whicht0", &t_whicht0, "whicht0/I");
  fTrackMatchTree->Branch("flash_id", &t_flash_id, "flash_id/I");
  fTrackMatchTree->Branch("t0", &t_t0, "t0/F");
  fTrackMatchTree->Branch("start_x", &t_start_x, "start_x/F");
  fTrackMatchTree->Branch("start_y", &t_start_y, "start_y/F");
  fTrackMatchTree->Branch("start_z", &t_start_z, "start_z/F");
  fTrackMatchTree->Branch("end_x", &t_end_x, "end_x/F");
  fTrackMatchTree->Branch("end_y", &t_end_y, "end_y/F");
  fTrackMatchTree->Branch("end_z", &t_end_z, "end_z/F");
  fTrackMatchTree->Branch("dir_y", &t_dir_y, "dir_y/F");
  fTrackMatchTree->Branch("length", &t_length, "length/F");
  fTrackMatchTree->Branch("median_end_dqdx", &t_median_end_dqdx, "median_end_dqdx/F");
  fTrackMatchTree->Branch("centroid_y", &t_centroid_y, "centroid_y/F");
  fTrackMatchTree->Branch("centroid_z", &t_centroid_z, "centroid_z/F");
  fTrackMatchTree->Branch("flash_time", &t_flash_time, "flash_time/F");
  fTrackMatchTree->Branch("flash_pe", &t_flash_pe, "flash_pe/F");
  fTrackMatchTree->Branch("flash_y", &t_flash_y, "flash_y/F");
  fTrackMatchTree->Branch("flash_z", &t_flash_z, "flash_z/F");
  fTrackMatchTree->Branch("dt", &t_dt, "dt/F");
  fTrackMatchTree->Branch("radius", &t_radius, "radius/F");
  fTrackMatchTree->Branch("trange_low", &t_trange_low, "trange_low/F");
  fTrackMatchTree->Branch("trange_high", &t_trange_high, "trange_high/F");

  // --- truth ----------------------------------------------------------------
  if (!fMCParticleLabel.empty()) {
    std::string const name =
      fMCParticleLabel.label() + fMCParticleLabel.instance() + "_MCParticle";
    fMCParticleTree = tfs->make<TTree>(name.c_str(),
      ("simb::MCParticle truth: " + fMCParticleLabel.label()).c_str());
    fMCParticleTree->Branch("run", &m_run, "run/I");
    fMCParticleTree->Branch("subrun", &m_subrun, "subrun/I");
    fMCParticleTree->Branch("event", &m_event, "event/I");
    fMCParticleTree->Branch("track_gen_pdg", &track_gen_pdg, "track_gen_pdg/I");
    fMCParticleTree->Branch("track_gen_time", &track_gen_time, "track_gen_time/F");
    fMCParticleTree->Branch("track_gen_E", &track_gen_E, "track_gen_E/F");
    fMCParticleTree->Branch("track_gen_x", &track_gen_x, "track_gen_x/F");
    fMCParticleTree->Branch("track_gen_y", &track_gen_y, "track_gen_y/F");
    fMCParticleTree->Branch("track_gen_z", &track_gen_z, "track_gen_z/F");
    fMCParticleTree->Branch("track_end_time", &track_end_time, "track_end_time/F");
    fMCParticleTree->Branch("time_shift", &time_shift, "time_shift/F");
    fMCParticleTree->Branch("track_end_in_av", &track_end_in_av, "track_end_in_av/O");
    fMCParticleTree->Branch("michel_gen_pdg", &michel_gen_pdg, "michel_gen_pdg/I");
    fMCParticleTree->Branch("michel_process", &michel_process, "michel_process/I");
    fMCParticleTree->Branch("michel_gen_time", &michel_gen_time, "michel_gen_time/F");
    fMCParticleTree->Branch("michel_gen_E", &michel_gen_E, "michel_gen_E/F");
    fMCParticleTree->Branch("michel_gen_x", &michel_gen_x, "michel_gen_x/F");
    fMCParticleTree->Branch("michel_gen_y", &michel_gen_y, "michel_gen_y/F");
    fMCParticleTree->Branch("michel_gen_z", &michel_gen_z, "michel_gen_z/F");
    fMCParticleTree->Branch("michel_in_av", &michel_in_av, "michel_in_av/O");
  }

} // beginJob()


// -----------------------------------------------------------------------------
void icarus::ICARUSStoppingMuonOpticalAna::fillTrackMatchTree
  (std::vector<std::vector<TrackFlashMatch>> const& matches)
{
  for (std::size_t iCryo = 0; iCryo < matches.size(); ++iCryo) {
    for (TrackFlashMatch const& m : matches[iCryo]) {
      t_track_id        = m.trackID;
      t_cryo            = static_cast<int>(m.cryostat);
      t_whicht0         = m.whichT0;
      t_flash_id        = m.flashID;
      t_t0              = m.trackT0;
      t_start_x         = m.startX;
      t_start_y         = m.startY;
      t_start_z         = m.startZ;
      t_end_x           = m.endX;
      t_end_y           = m.endY;
      t_end_z           = m.endZ;
      t_dir_y           = m.dirY;
      t_length          = m.length;
      t_median_end_dqdx = m.medianEnddQdx;
      t_centroid_y      = m.centroidY;
      t_centroid_z      = m.centroidZ;
      t_flash_time      = m.flashTime;
      t_flash_pe        = m.flashPE;
      t_flash_y         = m.flashY;
      t_flash_z         = m.flashZ;
      t_dt              = m.deltaT;
      t_radius          = m.radius;
      t_trange_low      = m.trangeLow;
      t_trange_high     = m.trangeHigh;
      fTrackMatchTree->Fill();
    }
  }
}


// -----------------------------------------------------------------------------
void icarus::ICARUSStoppingMuonOpticalAna::fillFlashOpHits
  (std::vector<art::Ptr<recob::OpHit>> const& ophits)
{
  m_channel_id.clear();
  m_hit_x.clear(); m_hit_y.clear(); m_hit_z.clear();
  m_hit_start_time.clear(); m_hit_rise_time.clear(); m_hit_peak_time.clear();
  m_hit_width.clear(); m_hit_timing_corr.clear();
  m_hit_area.clear(); m_hit_pe.clear(); m_hit_amplitude.clear();

  std::unordered_map<int, float> sumpe_map;

  for (art::Ptr<recob::OpHit> const& ophit : ophits) {

    if (ophit->Amplitude() < fOpHitThresholdADC) continue;

    int const ch = ophit->OpChannel();
    m_channel_id.push_back(ch);

    auto const pos = getChannelXYZ(ch);
    m_hit_x.push_back(static_cast<float>(pos[0]));
    m_hit_y.push_back(static_cast<float>(pos[1]));
    m_hit_z.push_back(static_cast<float>(pos[2]));

    // RiseTime() is relative to StartTime(): make it absolute, as in the reference
    m_hit_start_time.push_back(ophit->StartTime());
    m_hit_peak_time.push_back(ophit->PeakTime());
    m_hit_rise_time.push_back(ophit->StartTime() + ophit->RiseTime());
    m_hit_timing_corr.push_back(static_cast<float>(getTimingCorrection(ch)));
    m_hit_width.push_back(ophit->Width());

    m_hit_area.push_back(static_cast<float>(ophit->Area()));           // ADC x tick
    m_hit_amplitude.push_back(static_cast<float>(ophit->Amplitude())); // ADC
    m_hit_pe.push_back(static_cast<float>(ophit->PE()));

    sumpe_map[ch] += ophit->PE();
  }

  m_multiplicity_left = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0,
    [this](int v, std::pair<int const, float> const& p)
    { return getSideByChannel(p.first) == 0 ? v + 1 : v; });

  m_multiplicity_right = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0,
    [this](int v, std::pair<int const, float> const& p)
    { return getSideByChannel(p.first) == 1 ? v + 1 : v; });

  m_multiplicity = m_multiplicity_left + m_multiplicity_right;
  m_flash_nhits  = static_cast<int>(m_channel_id.size());
}


// -----------------------------------------------------------------------------
/// Writes every flash, matched or not, and returns the absolute times of the
/// matched ones -- those select which waveform fragments are kept.
///
/// All flashes are written on purpose. The Michel decays 0-8 us after the muon
/// and often reconstructs as its own OpFlash, which might not be matched.
icarus::ICARUSStoppingMuonOpticalAna::FlashTimes
icarus::ICARUSStoppingMuonOpticalAna::fillFlashes
  (art::Event const& e, std::vector<std::vector<TrackFlashMatch>> const& matches)
{
  FlashTimes flashTimes(fFlashLabels.size());

  for (std::size_t iCryo = 0; iCryo < fFlashLabels.size(); ++iCryo) {

    // TrackFlashMatch::flashID indexes the flash product of this same cryostat
    std::set<int> matchedHere;
    for (TrackFlashMatch const& m : matches[iCryo])
      if (m.hasFlash()) matchedHere.insert(m.flashID);

    art::InputTag const& label = fFlashLabels[iCryo];
    auto const flashHandle = e.getHandle<std::vector<recob::OpFlash>>(label);
    if (!flashHandle.isValid()) {
      mf::LogWarning("ICARUSStoppingMuonOpticalAna")
        << "No recob::OpFlash with label '" << label.encode() << "'";
      continue;
    }

    // a selected track pointing at an index this product does not have means the
    // selector and this module are not reading the same collection
    for (int const idx : matchedHere) {
      if (idx < 0 || static_cast<std::size_t>(idx) >= flashHandle->size()) {
        throw art::Exception(art::errors::LogicError)
          << "Matched flash index " << idx << " is out of range for '"
          << label.encode() << "' (" << flashHandle->size() << " flashes). The"
             " selector and this module are reading different products.\n";
      }
    }

    art::FindManyP<recob::OpHit> ophitsPtr(flashHandle, e, label);

    for (std::size_t idx = 0; idx < flashHandle->size(); ++idx) {

      recob::OpFlash const& flash = (*flashHandle)[idx];

      m_flash_id       = static_cast<int>(idx);
      m_flash_cryo     = static_cast<int>(iCryo);
      m_flash_time     = flash.Time();
      m_sum_pe         = flash.TotalPE();
      m_flash_y        = flash.YCenter();
      m_flash_z        = flash.ZCenter();
      m_flash_width_y  = flash.YWidth();
      m_flash_width_z  = flash.ZWidth();
      m_track_matched  = (matchedHere.count(static_cast<int>(idx)) > 0);

      // absolute time: waveform timestamps live on that scale
      if (m_track_matched) flashTimes[iCryo].push_back(flash.AbsTime());

      fillFlashOpHits(ophitsPtr.at(idx));

      fOpFlashTree->Fill();
    }
  }

  return flashTimes;
}


// -----------------------------------------------------------------------------
/// Writes every OpHit that is not in a flash, over the whole readout and both
/// cryostats -- no time restriction.
///
/// OpHits are the cheap summary of all the light: they are reconstructed from
/// every fragment the readout produced, so keeping them all gives the light
/// record across the full 3 ms for ~1.2 MB/event. That is what covers the ~3% of
/// Michels decaying past the covering fragment's ~7.7 us: they have no waveform
/// samples, but they are still here as hits.
void icarus::ICARUSStoppingMuonOpticalAna::fillUnmatchedOpHits
  (art::Event const& e)
{
  for (std::size_t iLabel = 0; iLabel < fOpHitLabels.size(); ++iLabel) {

    art::InputTag const& label = fOpHitLabels[iLabel];
    auto const ophitHandle = e.getHandle<std::vector<recob::OpHit>>(label);
    if (!ophitHandle.isValid() || ophitHandle->empty()) {
      mf::LogWarning("ICARUSStoppingMuonOpticalAna")
        << "No recob::OpHit in collection with label '" << label.encode() << "'";
      continue;
    }

    unsigned kept = 0;

    // one association per cryostat: a hit counts as flash-matched only against
    // the flash collection of its own cryostat
    for (std::size_t iCryo = 0; iCryo < fFlashLabels.size(); ++iCryo) {

      art::FindOneP<recob::OpFlash> flashPtr(ophitHandle, e, fFlashLabels[iCryo]);
      if (!flashPtr.isValid()) continue;   // no flash product for this cryostat

      for (std::size_t idx = 0; idx < ophitHandle->size(); ++idx) {

        if (!flashPtr.at(idx).isNull()) continue;   // in a flash: not our business

        recob::OpHit const& ophit = (*ophitHandle)[idx];

        if (ophit.Amplitude() < fOpHitThresholdADC) continue;

        // the OpHit product spans both cryostats, so the channel still has to be
        // resolved: otherwise a hit of the other cryostat looks unmatched purely
        // because we asked the wrong association
        if (getCryostatByChannel(ophit.OpChannel()) != iCryo) continue;

        m_channel        = ophit.OpChannel();
        auto const pos   = getChannelXYZ(ophit.OpChannel());
        m_x              = static_cast<float>(pos[0]);
        m_y              = static_cast<float>(pos[1]);
        m_z              = static_cast<float>(pos[2]);
        m_integral       = static_cast<float>(ophit.Area());
        m_amplitude      = static_cast<float>(ophit.Amplitude());
        m_start_time     = ophit.StartTime();
        m_peak_time      = ophit.PeakTime();
        m_rise_time      = ophit.StartTime() + ophit.RiseTime();
        m_width          = ophit.Width();
        m_abs_start_time = ophit.PeakTimeAbs() + (m_start_time - m_peak_time);
        m_timing_corr    = static_cast<float>(getTimingCorrection(ophit.OpChannel()));
        m_pe             = static_cast<float>(ophit.PE());
        m_fast_to_total  = ophit.FastToTotal();

        fOpHitTrees[iLabel]->Fill();
        ++kept;
      }
    }

    mf::LogInfo("ICARUSStoppingMuonOpticalAna")
      << "'" << label.encode() << "': " << kept
      << " out-of-flash OpHits written";
  }
}


// -----------------------------------------------------------------------------
/// One row per fragment. `wf`/`bs` are filled only for the fragment covering a
/// matched flash; `saved` says which, and `nticks` is the true length either way.
void icarus::ICARUSStoppingMuonOpticalAna::fillWaveforms
  (art::Event const& e, FlashTimes const& flashTimes)
{
  if (!fSaveWaveforms) return;

  for (std::size_t iLabel = 0; iLabel < fOpDetWaveformLabels.size(); ++iLabel) {

    art::InputTag const& label = fOpDetWaveformLabels[iLabel];
    auto const waveHandle = e.getHandle<std::vector<raw::OpDetWaveform>>(label);
    if (!waveHandle.isValid() || waveHandle->empty()) {
      mf::LogWarning("ICARUSStoppingMuonOpticalAna")
        << "No raw::OpDetWaveform in collection with label '" << label.encode() << "'";
      continue;
    }

    std::size_t kept = 0;

    for (raw::OpDetWaveform const& wave : *waveHandle) {

      int const channel = wave.ChannelNumber();
      std::size_t const iCryo = getCryostatByChannel(channel);

      // Not a time-scale conversion: TimeStamp() is already absolute. 
      // This is the per-channel DB calibration, which OpHits and OpFlashes carry and a raw
      // waveform does not.
      double const correction = getTimingCorrection(channel);
      double const wStart = wave.TimeStamp() + correction;
      double const wEnd   = wStart + wave.Waveform().size() * fOpticalTickPeriod;

      bool covers = false;
      if (iCryo < flashTimes.size()) {
        for (double const t : flashTimes[iCryo]) {
          if (t >= wStart && t < wEnd) { covers = true; break; }
        }
      }

      m_wfchannel      = channel;
      m_nticks         = wave.Waveform().size();
      m_wfstart        = wave.TimeStamp();   // raw, as in ICARUSOpticalDebug
      m_wf_timing_corr = static_cast<float>(correction);
      m_wf_saved       = covers;

      if (covers) {
        m_wf = wave.Waveform();

        fPedAlgo->Evaluate(wave.Waveform());
        auto const& mean = fPedAlgo->Mean();
        m_bs.assign(mean.begin(), mean.end());   // double -> float

        ++kept;
      }
      else {
        // metadata-only row: nticks still reports the true fragment length, so
        // `saved` (or wf.empty()) is what distinguishes the two
        m_wf.clear();
        m_bs.clear();
      }

      fOpDetWaveformTrees[iLabel]->Fill();
    }

    mf::LogInfo("ICARUSStoppingMuonOpticalAna")
      << "'" << label.encode() << "': samples kept for " << kept << " of "
      << waveHandle->size() << " waveform fragments (metadata for all)";
  }
}


// -----------------------------------------------------------------------------
void icarus::ICARUSStoppingMuonOpticalAna::fillMCTruth(art::Event const& e)
{
  if (fMCParticleLabel.empty() || fMCParticleTree == nullptr) return;

  auto const particleHandle = e.getHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);
  if (!particleHandle.isValid() || particleHandle->empty()) {
    mf::LogWarning("ICARUSStoppingMuonOpticalAna")
      << "No simb::MCParticle with label '" << fMCParticleLabel.encode() << "'";
    return;
  }

  // The 'shifted' producer (sbncode AdjustSimForTrigger) adds a per-event offset
  // to the OpDetWaveforms and the SimPhotons, but NOT to these MCParticles, so the
  // times printed below are not on the waveform clock. The offset is
  //   shift = clockData.TriggerTime() - emulatedTrigger.TriggerTime()
  // and it changes event to event, so it cannot be folded in as a constant. 
  // To put a truth time on the waveform clock:
  //   t_rel = G4ToElecTime(T_mcparticle + shift_ns) - TriggerTime()
  double timeShift_ns = 0.;    // G4 times -> the optical clock; 0 if unavailable

  if (!fTriggerLabel.empty()) {

    auto const triggerHandle = e.getHandle<std::vector<raw::Trigger>>(fTriggerLabel);

    if (!triggerHandle) {
      std::cout << "TimeShift: no raw::Trigger with label '"
                << fTriggerLabel.encode() << "'" << std::endl;
    }
    else if (triggerHandle->size() != 1) {
      std::cout << "TimeShift: " << triggerHandle->size()
                << " raw::Trigger in this event, expected exactly 1" << std::endl;
    }
    else {
      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
      double const emulated = triggerHandle->front().TriggerTime();   // us
      double const nominal  = clockData.TriggerTime();                // us
      double const shift    = nominal - emulated;                     // us
      timeShift_ns = shift * 1000.;
      time_shift   = static_cast<float>(timeShift_ns);
      std::cout << "TimeShift: nominal=" << nominal
                << " us  emulated=" << emulated
                << " us  shift=" << shift << " us (" << shift * 1000. << " ns)"
                << "  G4ToElecTime(0)=" << clockData.G4ToElecTime(0.) << " us"
                << std::endl;
    }
  }

  // reset both blocks every event: otherwise an event with no primary muon
  // reports the previous event's muon
  track_gen_pdg = -1;
  track_gen_time = -1.f; track_gen_E = -1.f;
  track_gen_x = -1.f; track_gen_y = -1.f; track_gen_z = -1.f;
  track_end_time = -1.f;
  track_end_in_av = false;
  michel_gen_pdg = -1;
  michel_process = -1;
  michel_gen_time = -1.f; michel_gen_E = -1.f;
  michel_gen_x = -1.f; michel_gen_y = -1.f; michel_gen_z = -1.f;
  michel_in_av = false;

  auto const inActiveVolume = [this](double x, double y, double z) {
    geo::Point_t const p { x, y, z };
    for (auto const& cryo : fGeom->Iterate<geo::CryostatGeo>()) {
      for (auto const& tpc : fGeom->Iterate<geo::TPCGeo>(cryo.ID())) {
        if (tpc.ActiveBoundingBox().ContainsPosition(p)) return true;
      }
    }
    return false;
  };

  // first pass: find the primary muon
  int primaryMuonID = -1;
  simb::MCParticle const* primaryMuon = nullptr;
  unsigned nPrimaryMuons = 0;

  for (simb::MCParticle const& par : *particleHandle) {
    if (std::abs(par.PdgCode()) != 13 || par.Mother() != 0) continue;
    ++nPrimaryMuons;
    if (primaryMuon) continue;   // keep the first, do not let the last overwrite
    primaryMuonID = par.TrackId();
    primaryMuon   = &par;
  }

  if (nPrimaryMuons > 1) {
    mf::LogWarning("ICARUSStoppingMuonOpticalAna")
      << nPrimaryMuons << " primary muons in this event; the truth tree describes"
         " the first one only";
  }

  if (primaryMuon) {
    track_gen_pdg  = primaryMuon->PdgCode();
    track_gen_time = primaryMuon->T();     // ns
    track_gen_E    = primaryMuon->E();     // GeV
    track_gen_x    = primaryMuon->Vx();
    track_gen_y    = primaryMuon->Vy();
    track_gen_z    = primaryMuon->Vz();
    track_end_time = primaryMuon->EndT();   // ns
    track_end_in_av =
      inActiveVolume(primaryMuon->EndX(), primaryMuon->EndY(), primaryMuon->EndZ());
  }

  // second pass: which electron daughter, if any, is the Michel.
  //
  // Geant4 reports mu- bound decay and mu- nuclear capture under one process name,
  // "muMinusCaptureAtRest", so Process() cannot tell them apart, while mu+ 
  // univocally decays as "Decay".

  // For mu-: only the decay emits an anti-nu_e (capture is mu- + p -> n + nu_mu). 
  // The atomic cascade runs in both branches, so a bound-decay mu- has 
  // several e- daughters; the Michel is the most energetic one, 
  // tens of MeV against sub-MeV Auger electrons.
  if (primaryMuonID >= 0) {

    constexpr int    kPdgAntiNuE      = -12;
    constexpr double kElectronMassGeV = 0.000510998946;
    constexpr double kMichelMinKEGeV  = 0.002;   // 2 MeV: far above the Auger cascade,
                                                 // far below the 30-53 MeV Michel

    // the anti-nu_e must be a daughter of this muon
    bool neutrinosStored = false;
    bool hasAntiNuE      = false;
    for (simb::MCParticle const& par : *particleHandle) {
      int const absPdg = std::abs(par.PdgCode());
      if (absPdg == 12 || absPdg == 14 || absPdg == 16) neutrinosStored = true;
      if (par.PdgCode() == kPdgAntiNuE && par.Mother() == primaryMuonID) hasAntiNuE = true;
    }

    // the Michel candidate: the most energetic electron daughter of the muon
    // coming from either decay or capture.
    simb::MCParticle const* hardestElectron = nullptr;
    for (simb::MCParticle const& par : *particleHandle) {

      std::cout << "Particle: TrackID=" << par.TrackId() << " PDG=" << par.PdgCode() << " Mother=" << par.Mother()
                << " Process=" << par.Process() << " EndProcess=" << par.EndProcess()
                << " E=" << par.E()
                << " T=" << par.T() << " EndT=" << par.EndT()
                << " nTraj=" << par.NumberTrajectoryPoints() << std::endl;

      if (par.Mother() != primaryMuonID) continue;
      if (std::abs(par.PdgCode()) != 11) continue;
      if (par.Process() != "Decay" && par.Process() != "muMinusCaptureAtRest") continue;
      if (!hardestElectron || par.E() > hardestElectron->E()) hardestElectron = &par;
    }

    if (hardestElectron) {

      int process = -1;
      if (hardestElectron->Process() == "Decay") {
        process = 0;                                  // mu+ (or a mu- decaying in flight)
      }
      else if (neutrinosStored) {                     // if neutrinos are not stored, use them!
        process = hasAntiNuE ? 2 : 1;                 // bound decay vs nuclear capture
      }
      else {
        double const ke = hardestElectron->E() - kElectronMassGeV;
        process = (ke >= kMichelMinKEGeV) ? 2 : 1;    // fallback: energy on its own
      }

      // on the capture branch this is the hardest Auger electron, not a Michel:
      // michel_process == 1 flags that
      michel_process  = process;
      michel_gen_pdg  = hardestElectron->PdgCode();
      michel_gen_time = hardestElectron->T();   // ns
      michel_gen_E    = hardestElectron->E();   // GeV
      michel_gen_x    = hardestElectron->Vx();
      michel_gen_y    = hardestElectron->Vy();
      michel_gen_z    = hardestElectron->Vz();
      michel_in_av    = inActiveVolume(hardestElectron->Vx(), hardestElectron->Vy(),
                                       hardestElectron->Vz());
    }
  }

  // --- SimPhotons: per-track photon map, binned in time ------------------------
  //
  // sim::OnePhoton::MotherTrackID stamps every detected photon with the G4 track
  // that made it, which is what SimPhotonsLite drops and what makes a light-based
  // truth label possible. Attribution rather than a time cut is required: the
  // muon's own slow component is still arriving underneath any delayed peak.
  //
  // Bins are 100 ns, indexed relative to the muon disappearance on the *photon*
  // clock (the photons carry the AdjustSimForTrigger shift, the MCParticles do
  // not). Bin 0 therefore holds the delayed peak, and is used as the window for
  //   S = non-muon photons in bin 0   (at the disappearance, non-muon == delayed)
  //   B = muon photons in bin 0       (the slow tail the peak has to stand on)
  if (!fSimPhotonsLabels.empty() && primaryMuon) {

    std::unordered_map<int, simb::MCParticle const*> byID;
    std::unordered_map<int, std::vector<int>> children;
    for (simb::MCParticle const& par : *particleHandle) {
      byID[par.TrackId()] = &par;
      children[par.Mother()].push_back(par.TrackId());
    }

    // disappearance = latest direct daughter. Correct for both charges: for mu-
    // it takes the capture/decay products over the atomic cascade, and for mu+
    // every daughter sits at the decay anyway.
    double tDisappear = std::numeric_limits<double>::lowest();
    for (int const id : children[primaryMuonID]) {
      tDisappear = std::max(tDisappear, byID[id]->T());
    }
    double const tDisShifted = tDisappear + timeShift_ns;

    constexpr double kBinNs    = 100.;
    constexpr int    kBinFirst = -10;
    constexpr int    kBinLast  =  30;

    for (art::InputTag const& tag : fSimPhotonsLabels) {

      auto const photonHandle = e.getHandle<std::vector<sim::SimPhotons>>(tag);
      if (!photonHandle) {
        std::cout << "SimPhotons '" << tag.encode() << "': NOT PRESENT" << std::endl;
        continue;
      }

      std::map<int, std::size_t> total;                  // trackId -> photons
      std::map<int, std::map<int, std::size_t>> binned;   // trackId -> bin -> photons
      std::size_t nPhotons = 0, nZero = 0, nUnresolved = 0, nOutside = 0;

      for (sim::SimPhotons const& channel : *photonHandle) {
        for (sim::OnePhoton const& photon : channel) {

          ++nPhotons;
          if (photon.MotherTrackID == 0) ++nZero;
          int const mother = std::abs(photon.MotherTrackID);  // EM daughters are signed
          if (byID.count(mother) == 0) ++nUnresolved;
          ++total[mother];

          int const bin = static_cast<int>(
            std::floor((photon.Time - tDisShifted) / kBinNs));
          if (bin >= kBinFirst && bin <= kBinLast) ++binned[mother][bin];
          else ++nOutside;
        }
      }

      std::cout << "SimPhotons '" << tag.encode() << "': " << photonHandle->size()
                << " channels, " << nPhotons << " photons"
                << "  (MotherTrackID zero=" << nZero
                << " unresolved=" << nUnresolved
                << " distinct=" << total.size() << ")\n"
                << "  disappearance " << tDisappear << " ns (G4) -> "
                << tDisShifted << " ns (photon clock); bins are 100 ns from there,"
                   " " << nOutside << " photons outside ["
                << kBinFirst << "," << kBinLast << "]\n"
                << "  trackId      pdg        T[ns]    total   bins (index:photons)"
                << std::endl;

      std::vector<std::pair<std::size_t, int>> order;
      order.reserve(total.size());
      for (auto const& entry : total) order.emplace_back(entry.second, entry.first);
      std::sort(order.rbegin(), order.rend());

      for (auto const& [count, trackId] : order) {
        auto const it = byID.find(trackId);
        std::cout << std::setw(9) << trackId
                  << std::setw(9) << (it == byID.end() ? 0 : it->second->PdgCode())
                  << std::setw(13) << (it == byID.end() ? -1. : it->second->T())
                  << std::setw(9) << count << "   ";
        auto const bit = binned.find(trackId);
        if (bit != binned.end()) {
          for (auto const& [bin, n] : bit->second) std::cout << bin << ":" << n << " ";
        }
        std::cout << std::endl;
      }

      // S/B in bin 0, the 100 ns window opening at the disappearance
      std::size_t S = 0, B = 0;
      for (auto const& [trackId, bins] : binned) {
        auto const it = bins.find(0);
        if (it == bins.end()) continue;
        if (trackId == primaryMuonID) B += it->second;
        else                          S += it->second;
      }
      std::cout << "  bin 0 = [" << tDisShifted << ", " << tDisShifted + kBinNs
                << ") ns:  S(non-muon)=" << S << "  B(muon)=" << B << "  S/B=";
      if (B > 0) std::cout << (double(S) / B);
      else       std::cout << "n/a (no muon light in the window)";
      std::cout << std::endl;
    }
  }

  fMCParticleTree->Fill();
}


// -----------------------------------------------------------------------------
void icarus::ICARUSStoppingMuonOpticalAna::analyze(art::Event const& e)
{

  m_run       = e.id().run();
  m_subrun    = e.id().subRun();
  m_event     = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second

  // --- selection, before anything is written --------------------------------
  std::vector<std::vector<TrackFlashMatch>> matches(fSelector.nCryostats());
  bool anyMatched = false;

  for (std::size_t iCryo = 0; iCryo < fSelector.nCryostats(); ++iCryo) {
    matches[iCryo] = fSelector.select(e, iCryo);
    for (TrackFlashMatch const& m : matches[iCryo])
      if (m.hasFlash()) anyMatched = true;
  }

  // if no tracks were matched at all, and you require matches
  // do not save anything from this event
  if (fSkipUnmatchedEvents && !anyMatched){
    mf::LogWarning("ICARUSStoppingMuonOpticalAna") 
      << "Zero track matches in run " << m_run << " subrun " << m_subrun << " event " << m_event << "!";
    return;
  }

  // save the matches in the file
  fillTrackMatchTree(matches);

  // --- optical dump ---------------------------------------------------------
  // Flashes and OpHits are written in full: they are the cheap record of all the
  // light in the event. Only the waveform samples are restricted, to the fragments 
  // covering a flash matched with a stopping muon.
  FlashTimes const flashTimes = fillFlashes(e, matches);
  fillWaveforms(e, flashTimes);
  fillUnmatchedOpHits(e);
  fillMCTruth(e);

} // analyze()


DEFINE_ART_MODULE(icarus::ICARUSStoppingMuonOpticalAna)
