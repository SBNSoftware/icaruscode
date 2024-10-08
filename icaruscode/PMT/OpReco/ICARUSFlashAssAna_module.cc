////////////////////////////////////////////////////////////////////////
// Class:       ICARUSFlashAssAna
// Plugin Type: analyzer (art v3_06_03)
// File:        ICARUSFlashAssAna_module.cc
//
// Generated at Tue Jun 29 13:43:54 2021 by Andrea Scarpelli using cetskelgen
// from cetlib version v3_11_01.
//
// Module that dumps the association between Flashes and OpHit.
// These trees make up the optical information in the calibration ntuples.
//
// mailto:ascarpel@bnl.gov, mvicenzi@bnl.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "art_root_io/TFileService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/IcarusObj/PMTBeamSignal.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

#include "TTree.h"

#include <vector>
#include <map>
#include <numeric> // std::accumulate
#include <limits>
#include <cstddef>

namespace opana
{
  class ICARUSFlashAssAna;
}

class opana::ICARUSFlashAssAna : public art::EDAnalyzer
{

public:
  struct Config
  {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> TriggerLabel{
        Name("TriggerLabel"),
        Comment("Label for the Trigger fragment label")};

    fhicl::Atom<bool> DumpWaveformsInfo{
        Name("DumpWaveformsInfo"),
        Comment("Set the option to save some aggregated waveform information")};

    fhicl::Atom<bool> SaveRawWaveforms{
        Name("SaveRawWaveforms"),
        Comment("Set to save the full raw::OpDetWaveforms")};

    fhicl::Atom<bool> UseSharedBaseline{
        Name("UseSharedBaseline"),
        Comment("Set the option to use icarus::WaveformBaseline")};

    fhicl::Sequence<art::InputTag> OpDetWaveformLabels{
        Name("OpDetWaveformLabels"),
        Comment("Tags for the raw::OpDetWaveform data products")};

    fhicl::Sequence<art::InputTag> BaselineLabels{
        Name("BaselineLabels"),
        Comment("Tags for the icarus::WaveformBaseline data products")};

    fhicl::Sequence<art::InputTag> OpHitLabels{
        Name("OpHitLabels"),
        Comment("Tags for the recob::OpHit data products")};

    fhicl::Sequence<art::InputTag> FlashLabels{
        Name("FlashLabels"),
        Comment("Tags for the recob::Flashe data products")};

    fhicl::Atom<art::InputTag> RWMLabel{
        Name("RWMLabel"),
        Comment("Tag for the RWM std::vector<icarus::timing::PMTBeamSignal> data product")};

    fhicl::Atom<float> PEOpHitThreshold{
        Name("PEOpHitThreshold"),
        Comment("Threshold in PE for an OpHit to be considered in the information calculated for a flash")};

  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit ICARUSFlashAssAna(Parameters const &config);

  ICARUSFlashAssAna(ICARUSFlashAssAna const &) = delete;
  ICARUSFlashAssAna(ICARUSFlashAssAna &&) = delete;
  ICARUSFlashAssAna &operator=(ICARUSFlashAssAna const &) = delete;
  ICARUSFlashAssAna &operator=(ICARUSFlashAssAna &&) = delete;

  void analyze(art::Event const &e) override;
  void beginJob() override;

  /// Compute median of waveform for baseline (fallback option)
  template <typename T>
  T Median(std::vector<T> data) const;

  /// Return cryostat from PMT channel_id
  geo::CryostatID::CryostatID_t getCryostatByChannel(int channel);

  /// Return wall from PMT channel_id
  int getSideByChannel(const int channel);

  /// Process OpHits in the absence of flashes
  void processOpHits(art::Event const &e, unsigned int cryo);

  /// Process Ophits in the presence of flashes
  void processOpHitsFlash(std::vector<art::Ptr<recob::OpHit>> const &ophits,
                          int &multiplicity_left, int &multiplicity_right,
                          float &sum_pe_left, float &sum_pe_right,
                          std::vector<double> &pmt_start_time,
                          std::vector<double> &pmt_rise_time,
                          std::vector<double> &pmt_start_rwm_time,
                          std::vector<double> &pmt_pe,
                          std::vector<double> &pmt_amplitude,
                          TTree *ophittree);

  /// Return RWM-relative time from a trigger-relative time
  float getRWMRelativeTime(int channel, float t);

  /// Return the RWM-relative flash interaction time
  float getFlashBunchTime(std::vector<double> pmt_start_time_rwm,
                          std::vector<double> pmt_rise_time);

private:
  //----------
  // Input parameters

  art::InputTag fTriggerLabel;
  bool fSaveWaveformInfo;
  bool fSaveRawWaveforms;
  bool fUseSharedBaseline;
  std::vector<art::InputTag> fOpDetWaveformLabels;
  std::vector<art::InputTag> fBaselineLabels;
  std::vector<art::InputTag> fOpHitLabels;
  std::vector<art::InputTag> fFlashLabels;
  art::InputTag fRWMLabel;
  float fPEOpHitThreshold;
  bool fDebug;

  //----------
  // Output trees

  TTree *fEventTree;
  std::vector<TTree *> fOpDetWaveformTrees;
  std::vector<TTree *> fOpFlashTrees;
  std::vector<TTree *> fOpHitTrees;
  std::vector<TTree *> fOpHitFlashTrees;

  //----------------
  // Output variables

  // Common
  int m_run;
  int m_event;
  int m_timestamp;

  // Event/trigger tree
  float m_beam_gate_start = -99999;
  float m_beam_gate_width = -99999;
  int m_beam_type = -1;
  int m_trigger_type = -1;
  unsigned int m_gate_type;
  std::string m_gate_name;
  uint64_t m_trigger_timestamp;
  uint64_t m_gate_start_timestamp;
  uint64_t m_trigger_gate_diff;
  uint64_t lvdsCryoE[2];
  uint64_t lvdsCryoW[2];
  uint16_t addersCryoE[2];
  uint16_t addersCryoW[2];

  // Flash trees
  int m_flash_id;
  int m_multiplicity;
  int m_multiplicity_left;
  int m_multiplicity_right;
  float m_sum_pe;
  float m_sum_pe_left;
  float m_sum_pe_right;
  float m_flash_time;
  float m_flash_time_rwm;
  float m_flash_y;
  float m_flash_width_y;
  float m_flash_z;
  float m_flash_width_z;
  std::vector<double> m_pmt_start_time;
  std::vector<double> m_pmt_rise_time;
  std::vector<double> m_pmt_start_time_rwm;
  std::vector<double> m_pmt_pe;
  std::vector<double> m_pmt_amplitude;

  // Ophit trees
  int m_channel_id;
  float m_integral;  // in ADC x tick
  float m_amplitude; // in ADC
  float m_start_time;
  float m_peak_time;
  float m_rise_time;
  float m_width;
  float m_abs_start_time;
  float m_start_time_rwm;
  float m_peak_time_rwm;
  float m_pe;
  float m_fast_to_total;

  // Waveform trees
  float m_wf_start;
  short m_baseline;
  short m_chargesum;
  int m_nticks;
  std::vector<short> m_wf;

  // Geometry tree
  std::vector<float> m_pmt_x;
  std::vector<float> m_pmt_y;
  std::vector<float> m_pmt_z;

  //----------
  // Support variables/products

  geo::GeometryCore const *fGeom;
  std::vector<icarus::timing::PMTBeamSignal> fRWMTimes;
};

// ----------------------------------------------------------------------------

opana::ICARUSFlashAssAna::ICARUSFlashAssAna(Parameters const &config)
    : EDAnalyzer(config),
      fTriggerLabel(config().TriggerLabel()),
      fSaveWaveformInfo(config().DumpWaveformsInfo()),
      fSaveRawWaveforms(config().SaveRawWaveforms()),
      fUseSharedBaseline(config().UseSharedBaseline()),
      fOpDetWaveformLabels(config().OpDetWaveformLabels()),
      fBaselineLabels(config().BaselineLabels()),
      fOpHitLabels(config().OpHitLabels()),
      fFlashLabels(config().FlashLabels()),
      fRWMLabel(config().RWMLabel()),
      fPEOpHitThreshold(config().PEOpHitThreshold()),
      fGeom(lar::providerFrom<geo::Geometry>())
{
}

// ----------------------------------------------------------------------------

void opana::ICARUSFlashAssAna::beginJob()
{

  art::ServiceHandle<art::TFileService const> tfs;

  // Setting up the GEOMETRY tree
  // Channel id corresponds to vector index
  TTree *fGeoTree = tfs->make<TTree>("geotree", "geometry information");
  fGeoTree->Branch("pmt_x", &m_pmt_x);
  fGeoTree->Branch("pmt_y", &m_pmt_y);
  fGeoTree->Branch("pmt_z", &m_pmt_z);

  for (std::size_t opch = 0; opch < fGeom->NOpChannels(); ++opch)
  {

    auto const PMTxyz = fGeom->OpDetGeoFromOpChannel(opch).GetCenter();
    m_pmt_x.push_back(PMTxyz.X());
    m_pmt_y.push_back(PMTxyz.Y());
    m_pmt_z.push_back(PMTxyz.Z());
  }

  fGeoTree->Fill();

  // Settinp up the EVENT tree
  // Trigger information and LVDS status
  fEventTree = tfs->make<TTree>("eventstree", "higher level information on the event");
  fEventTree->Branch("run", &m_run, "run/I");
  fEventTree->Branch("event", &m_event, "event/I");
  fEventTree->Branch("timestamp", &m_timestamp, "timestamp/I");
  fEventTree->Branch("beam_gate_start", &m_beam_gate_start, "beam_gate_start/F");
  fEventTree->Branch("beam_gate_width", &m_beam_gate_width, "beam_gate_width/F");
  fEventTree->Branch("beam_type", &m_beam_type, "beam_type/I");
  fEventTree->Branch("gate_type", &m_gate_type, "gate_type/b");
  fEventTree->Branch("gate_name", &m_gate_name);
  fEventTree->Branch("trigger_type", &m_trigger_type, "trigger_type/I");
  fEventTree->Branch("trigger_timestamp", &m_trigger_timestamp, "trigger_timestamp/l");
  fEventTree->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");
  fEventTree->Branch("trigger_gate_diff", &m_trigger_gate_diff, "trigger_gate_diff/l");
  fEventTree->Branch("lvdsCryoE", &lvdsCryoE, "lvdsCryoE[2]/l");
  fEventTree->Branch("lvdsCryoW", &lvdsCryoW, "lvdsCryoW[2]/l");
  fEventTree->Branch("addersCryoE", &addersCryoE, "addersCryoE[2]/s");
  fEventTree->Branch("addersCryoW", &addersCryoW, "addersCryoW[2]/s");

  // Setting up the WAVEFORM trees (one per product label)
  // This tree will hold some aggregated optical waveform information
  // The flag must be enabled to have the information saved
  if (!fOpDetWaveformLabels.empty() && fSaveWaveformInfo)
  {

    for (auto const &label : fOpDetWaveformLabels)
    {

      std::string name = label.label() + label.instance() + "wfttree";
      std::string info = "TTree with aggregated optical waveform information with label: " + label.label();

      TTree *ttree = tfs->make<TTree>(name.c_str(), info.c_str());
      ttree->Branch("run", &m_run, "run/I");
      ttree->Branch("event", &m_event, "event/I");
      ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
      ttree->Branch("channel_id", &m_channel_id, "channel_id/I");
      ttree->Branch("wf_start", &m_wf_start, "wf_start/F");
      ttree->Branch("baseline", &m_baseline, "baseline/s");
      ttree->Branch("chargesum", &m_chargesum, "chargesum/s");
      ttree->Branch("nticks", &m_nticks, "nticks/I");
      if (fSaveRawWaveforms)
        ttree->Branch("wf", &m_wf);

      fOpDetWaveformTrees.push_back(ttree);
    }
  }

  // Setting up the OPHIT trees (one per cryostat)
  // This ttree will hold the ophit information when a flash is not found in the event
  // NB: information of the optical hits in events where flashes are present are lost
  for (auto const &label : fOpHitLabels)
  {

    std::string name = label.label() + "_ttree";
    std::string info = "TTree for the recob::OpHit objects with label " + label.label() + " in events without flashes.";

    TTree *ttree = tfs->make<TTree>(name.c_str(), info.c_str());
    ttree->Branch("run", &m_run, "run/I");
    ttree->Branch("event", &m_event, "event/I");
    ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
    ttree->Branch("channel_id", &m_channel_id, "channel_id/I");
    ttree->Branch("integral", &m_integral, "integral/F");
    ttree->Branch("amplitude", &m_amplitude, "amplitude/F");
    ttree->Branch("start_time", &m_start_time, "start_time/F");
    ttree->Branch("peak_time", &m_peak_time, "peak_time/F");
    ttree->Branch("rise_time", &m_rise_time, "rise_time/F");
    ttree->Branch("abs_start_time", &m_abs_start_time, "abs_start_time/F");
    ttree->Branch("start_time_rwm", &m_start_time_rwm, "start_time_rwm/F");
    ttree->Branch("peak_time_rwm", &m_peak_time_rwm, "peak_time_rwm/F");
    ttree->Branch("pe", &m_pe, "pe/F");
    ttree->Branch("width", &m_width, "width/F");
    ttree->Branch("fast_to_total", &m_fast_to_total, "fast_to_total/F");

    fOpHitTrees.push_back(ttree);
  }

  // Setting up the OPFLASH/OPHITS trees (one per cryostat)
  // These ttrees hold the information for the ophits and the flashes
  // NB: information of the optical hits is stored differently when flashes are found
  if (!fFlashLabels.empty())
  {

    for (auto const &label : fFlashLabels)
    {

      // TTree for the flash in a given cryostat
      std::string name = label.label() + "_flashtree";
      std::string info = "TTree for the recob::Flashes with label " + label.label();

      TTree *ttree = tfs->make<TTree>(name.c_str(), info.c_str());
      ttree->Branch("run", &m_run, "run/I");
      ttree->Branch("event", &m_event, "event/I");
      ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
      ttree->Branch("flash_id", &m_flash_id, "flash_id/I");
      ttree->Branch("multiplicity", &m_multiplicity, "multiplicity/I");
      ttree->Branch("multiplicity_right", &m_multiplicity_right, "multiplicity_right/I");
      ttree->Branch("multiplicity_left", &m_multiplicity_left, "multiplicity_left/I");
      ttree->Branch("sum_pe", &m_sum_pe, "sum_pe/F");
      ttree->Branch("sum_pe_right", &m_sum_pe_right, "sum_pe_right/F");
      ttree->Branch("sum_pe_left", &m_sum_pe_left, "sum_pe_left/F");
      ttree->Branch("flash_time", &m_flash_time, "flash_time/F");
      ttree->Branch("flash_time_rwm", &m_flash_time_rwm, "flash_time_rwm/F");
      ttree->Branch("flash_y", &m_flash_y, "flash_y/F");
      ttree->Branch("flash_width_y", &m_flash_width_y, "flash_width_y/F");
      ttree->Branch("flash_z", &m_flash_z, "flash_z/F");
      ttree->Branch("flash_width_z", &m_flash_width_z, "flash_width_z/F");
      ttree->Branch("pmt_x", &m_pmt_x);
      ttree->Branch("pmt_y", &m_pmt_y);
      ttree->Branch("pmt_z", &m_pmt_z);
      ttree->Branch("time_pmt", &m_pmt_start_time);
      ttree->Branch("time_pmt_rwm", &m_pmt_start_time_rwm);
      ttree->Branch("pe_pmt", &m_pmt_pe);
      ttree->Branch("amplitude_pmt", &m_pmt_amplitude);

      fOpFlashTrees.push_back(ttree);

      // Now the ttree for the OpHit associated in the flash
      name = label.label() + "_ophittree";
      info = "Three for the recob::OpHit associated with an OpHitFlash" + label.label();

      TTree *ophittree = tfs->make<TTree>(name.c_str(), info.c_str());
      ophittree->Branch("run", &m_run, "run/I");
      ophittree->Branch("event", &m_event, "event/I");
      ophittree->Branch("timestamp", &m_timestamp, "timestamp/I");
      ophittree->Branch("flash_id", &m_flash_id, "flash_id/I");
      ophittree->Branch("channel_id", &m_channel_id, "channel_id/I");
      ophittree->Branch("integral", &m_integral, "integral/F");
      ophittree->Branch("amplitude", &m_amplitude, "amplitude/F");
      ophittree->Branch("start_time", &m_start_time, "start_time/F");
      ophittree->Branch("peak_time", &m_peak_time, "peak_time/F");
      ophittree->Branch("rise_time", &m_rise_time, "rise_time/F");
      ophittree->Branch("abs_start_time", &m_abs_start_time, "abs_start_time/F");
      ophittree->Branch("start_time_rwm", &m_start_time_rwm, "start_time_rwm/F");
      ophittree->Branch("peak_time_rwm", &m_peak_time_rwm, "peak_time_rwm/F");
      ophittree->Branch("pe", &m_pe, "pe/F");
      ophittree->Branch("width", &m_width, "width/F");
      ophittree->Branch("fast_to_total", &m_fast_to_total, "fast_to_total/F");

      fOpHitFlashTrees.push_back(ophittree);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename T>
T opana::ICARUSFlashAssAna::Median(std::vector<T> data) const
{

  std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
  return data[data.size() / 2];
}

// ----------------------------------------------------------------------------

geo::CryostatID::CryostatID_t opana::ICARUSFlashAssAna::getCryostatByChannel(int channel)
{

  const geo::OpDetGeo &opdetgeo = fGeom->OpDetGeoFromOpChannel(channel);
  geo::CryostatID::CryostatID_t cid = opdetgeo.ID().Cryostat;
  return cid;
}

// ----------------------------------------------------------------------------

int opana::ICARUSFlashAssAna::getSideByChannel(const int channel)
{

  /*
  Channels are numbered from east to west, from North (cryo side) to South (beam side)
  We look in the opposide direction wrt to the beam direction South->North:

  - Left is the east wall of each cryostat;
  - Right is the west side of each cryostat;
  - [ 0:89 ] and [180:269] are on the left,
    the return value of the function is 0;
  - [ 90-179 ] and [ 270:359 ] are on the right,
    the return value of the function is 1;
  */

  int side = channel / 90; // always round down
  return side % 2;
}

// ----------------------------------------------------------------------------

float opana::ICARUSFlashAssAna::getRWMRelativeTime(int channel, float t)
{

  if (fRWMTimes.empty())
    return icarus::timing::NoTime;

  auto rwm = fRWMTimes.at(channel);
  if (!rwm.isValid())
  {
    mf::LogTrace("ICARUSFlashAssAna") << "No RWM signal for channel " << channel << " "
                                      << "(Crate " << rwm.crate << ", Board " << rwm.digitizerLabel
                                      << ", SpecialChannel " << rwm.specialChannel << ")"
                                      << " in event " << m_event << " gate " << m_gate_name;
    return icarus::timing::NoTime;
  }

  float rwm_trigger = rwm.startTime; // rwm time w.r.t. trigger time [us]
  return (t - rwm_trigger);
}

// ----------------------------------------------------------------------------

float opana::ICARUSFlashAssAna::getFlashBunchTime(std::vector<double> pmt_start_time_rwm,
                                                  std::vector<double> pmt_rise_time)
{

  float tfirst_left = std::numeric_limits<float>::max();
  float tfirst_right = std::numeric_limits<float>::max();

  // if no RWM info available, all pmt_start_time_rwm are invalid
  // return icarus::timing::NoTime as well for the flash
  if (fRWMTimes.empty())
    return icarus::timing::NoTime;

  int nleft = 0;
  int nright = 0;
  for (std::size_t i = 0; i < pmt_start_time_rwm.size(); i++)
  {

    int side = getSideByChannel(i);
    float t = pmt_start_time_rwm[i] + pmt_rise_time[i]; // rise time w.r.t. rwm

    // exclude channels that have no signals
    // these have time == 0
    if (t > -1e-10 && t < 1e-10)
      continue;

    // if any RWM copy is missing (therefore missing for an entire PMT crate),
    // it might not be possible to use the first hits (they might not have a RMW time)
    // so return icarus::timing::NoTime as in other bad cases
    if (!fRWMTimes[i].isValid())
      return icarus::timing::NoTime;

    // count hits separetely on the two walls
    if (side == 0)
    {
      nleft++;
      if (t < tfirst_left)
        tfirst_left = t;
    }
    else if (side == 1)
    {
      nright++;
      if (t < tfirst_right)
        tfirst_right = t;
    }
  }

  // if there are no hits in one of the walls...
  if (nleft < 1 || nright < 1)
  {
    mf::LogWarning("ICARUSFlashAssAna") << "Flash " << m_flash_id << " doesn't have hits on both walls!"
                                        << "Left: " << nleft << " t " << tfirst_left << " "
                                        << "Right: " << nright << " t " << tfirst_right;
    // return what we have...
    return (tfirst_left < tfirst_right) ? tfirst_left : tfirst_right;
  }

  return (tfirst_left + tfirst_right) / 2.;
}

// ----------------------------------------------------------------------------

void opana::ICARUSFlashAssAna::processOpHits(art::Event const &e, unsigned int cryo)
{

  // if no OpHits have been selected at all (and no flashes as well!)
  if (fOpHitLabels.empty())
  {
    mf::LogError("ICARUSFlashAssAna") << "No recob::OpHit labels selected.";
    return;
  }

  for (std::size_t iOpHitLabel = 0; iOpHitLabel < fOpHitLabels.size(); iOpHitLabel++)
  {

    auto const label = fOpHitLabels[iOpHitLabel];
    auto const &ophits = e.getProduct<std::vector<recob::OpHit>>(label);

    // we want our ophits to be valid and not empty
    if (ophits.empty())
    {
      mf::LogError("ICARUSFlashAssAna") << "Invalid recob::OpHit with label '" << label.encode() << "'";
      continue;
    }

    for (auto const &ophit : ophits)
    {

      const int channel_id = ophit.OpChannel();
      if (getCryostatByChannel(channel_id) != cryo)
      {
        continue;
      }

      m_channel_id = channel_id;
      m_integral = ophit.Area();       // in ADC x tick
      m_amplitude = ophit.Amplitude(); // in ADC
      m_width = ophit.Width();
      m_pe = ophit.PE();
      m_fast_to_total = ophit.FastToTotal();

      // save times: start, peak, rise
      // rise is relative to start
      m_start_time = ophit.StartTime();
      m_peak_time = ophit.PeakTime();
      m_rise_time = ophit.RiseTime();
      m_start_time_rwm = getRWMRelativeTime(channel_id, m_start_time);
      m_peak_time_rwm = getRWMRelativeTime(channel_id, m_peak_time);
      m_abs_start_time = ophit.PeakTimeAbs() + (m_start_time - m_peak_time);

      fOpHitTrees[iOpHitLabel]->Fill();
    }
  }
}

// ----------------------------------------------------------------------------

void opana::ICARUSFlashAssAna::processOpHitsFlash(std::vector<art::Ptr<recob::OpHit>> const &ophits,
                                                  int &multiplicity_left, int &multiplicity_right,
                                                  float &sum_pe_left, float &sum_pe_right,
                                                  std::vector<double> &pmt_start_time,
                                                  std::vector<double> &pmt_rise_time,
                                                  std::vector<double> &pmt_start_time_rwm,
                                                  std::vector<double> &pmt_pe,
                                                  std::vector<double> &pmt_amplitude,
                                                  TTree *ophittree)
{

  // we calculate the total charge clustered in the flash per channel taking part to the flash
  // at the same time we store the times of the first ophit in each channel and its amplitude
  // same loop is used to saved OPHITS info

  std::unordered_map<int, float> sumpe_map;
  for (auto const ophit : ophits)
  {

    if (ophit->PE() < fPEOpHitThreshold)
    {
      continue;
    }

    const int channel_id = ophit->OpChannel();
    sumpe_map[channel_id] += ophit->PE();

    m_channel_id = channel_id;
    m_integral = ophit->Area();       // in ADC x tick
    m_amplitude = ophit->Amplitude(); // in ADC
    m_width = ophit->Width();
    m_pe = ophit->PE();
    m_fast_to_total = ophit->FastToTotal();

    // save times: start, peak, rise
    // rise is relative to start
    m_start_time = ophit->StartTime();
    m_peak_time = ophit->PeakTime();
    m_rise_time = ophit->RiseTime();
    m_start_time_rwm = getRWMRelativeTime(channel_id, m_start_time);
    m_peak_time_rwm = getRWMRelativeTime(channel_id, m_peak_time);
    m_abs_start_time = ophit->PeakTimeAbs() + (m_start_time - m_peak_time);

    pmt_pe[channel_id] += ophit->PE();

    // select the first ophit (by time) in each channel
    if ((pmt_start_time[channel_id] == 0) || (pmt_start_time[channel_id] > m_start_time))
    {
      pmt_start_time[channel_id] = m_start_time;
      pmt_rise_time[channel_id] = m_rise_time;
      pmt_start_time_rwm[channel_id] = m_start_time_rwm;
      pmt_amplitude[channel_id] = m_amplitude;
    }

    ophittree->Fill();
  }

  m_multiplicity_left = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0,
                                        [&](int value, const std::map<int, float>::value_type &p)
                                        { return getSideByChannel(p.first) == 0 ? ++value : value; });

  m_multiplicity_right = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0,
                                         [&](int value, const std::map<int, float>::value_type &p)
                                         { return getSideByChannel(p.first) == 1 ? ++value : value; });

  m_sum_pe_left = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0.0,
                                  [&](float value, const std::map<int, float>::value_type &p)
                                  { return getSideByChannel(p.first) == 0 ? value + p.second : value; });

  m_sum_pe_right = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0.0,
                                   [&](float value, const std::map<int, float>::value_type &p)
                                   { return getSideByChannel(p.first) == 1 ? value + p.second : value; });
}

// ----------------------------------------------------------------------------

void opana::ICARUSFlashAssAna::analyze(art::Event const &e)
{

  // Collect global event metadata
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second

  // -----
  // TRIGGER INFO
  // We work out the trigger information here

  if (!fTriggerLabel.empty())
  {

    // Beam gate information
    auto const &beamgateInfo = e.getProduct<std::vector<sim::BeamGateInfo>>(fTriggerLabel);

    if (!beamgateInfo.empty())
    {
      for (auto const &beamgate : beamgateInfo)
      {
        m_beam_gate_start = beamgate.Start();
        m_beam_gate_width = beamgate.Width();
        m_beam_type = beamgate.BeamType();
      }
    }
    else
    {
      mf::LogError("ICARUSFlashAssAna") << "No sim::BeamGateInfo associated to label: " << fTriggerLabel.label() << "\n";
    }

    // Now trigger information
    auto const &extraInfo = e.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel);

    if (extraInfo.isValid())
    {

      sbn::triggerSource bit = extraInfo.sourceType;
      m_gate_type = (unsigned int)bit;
      m_gate_name = bitName(bit);

      m_trigger_type = value(extraInfo.triggerType);
      m_trigger_timestamp = extraInfo.triggerTimestamp;
      m_gate_start_timestamp = extraInfo.beamGateTimestamp;
      m_trigger_gate_diff = m_trigger_timestamp - m_gate_start_timestamp;

      // majority lvds info
      lvdsCryoE[0] = extraInfo.cryostats[0].LVDSstatus[0];
      lvdsCryoE[1] = extraInfo.cryostats[0].LVDSstatus[1];
      lvdsCryoW[0] = extraInfo.cryostats[1].LVDSstatus[0];
      lvdsCryoW[1] = extraInfo.cryostats[1].LVDSstatus[1];

      // adders lvds info
      addersCryoE[0] = extraInfo.cryostats[0].sectorStatus[0];
      addersCryoE[1] = extraInfo.cryostats[0].sectorStatus[1];
      addersCryoW[0] = extraInfo.cryostats[1].sectorStatus[0];
      addersCryoW[1] = extraInfo.cryostats[1].sectorStatus[1];
    }
    else
    {
      mf::LogError("ICARUSFlashAssAna") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "!";
    }
  }

  // -----
  // RWM INFO
  // We work out the RWM information here
  // it might be empty if offbeam or missing, bu that's okay!

  if (!fRWMLabel.empty())
  {

    fRWMTimes = e.getProduct<std::vector<icarus::timing::PMTBeamSignal>>(fRWMLabel);
    if (fRWMTimes.empty())
      mf::LogTrace("ICARUSFlashAssAna") << "Data product std::vector<icarus::timing::PMTBeamSignal> for '"
                                        << fRWMLabel.label() << "' is empty in " << m_gate_name << " event!";
  }

  // -----
  // WAVEFORM INFO
  // Now we work on the waveforms if we are allowed to
  // Full raw waveforms are dumped only if option is set explicitly!
  // Waveforms are a complex business: if running at stage0, they're all available
  // if running at stage1, only on-beam waveforms are available but it's a different product!

  if (!fOpDetWaveformLabels.empty() && fSaveWaveformInfo)
  {

    for (std::size_t i = 0; i < fOpDetWaveformLabels.size(); i++)
    {

      auto const wflabel = fOpDetWaveformLabels[i];
      auto const bslabel = fBaselineLabels[i];
      auto const &waveforms = e.getProduct<std::vector<raw::OpDetWaveform>>(wflabel);
      auto const &baselines = e.getProduct<std::vector<icarus::WaveformBaseline>>(bslabel);

      // waveforms data product ("daqPMT") is dropped after stage0, and only "daqPMTonbeam" is kept
      // check if collection is present and print warning otherwise (requires change in fhicl config)
      if (!waveforms.empty())
      {

        std::size_t idx = 0;
        for (auto const &wave : waveforms)
        {

          m_channel_id = wave.ChannelNumber();
          m_nticks = wave.Waveform().size();
          m_wf_start = wave.TimeStamp();

          // try using icarus::WaveformBaseline, default to median if not
          // doesn't work with daqPMTonbeam as order gets messed up, so check if size is the same
          // FIXME: move to use art:Assn which can work both for daqPMT and daqPMTonbeam?
          m_baseline = (fUseSharedBaseline && (baselines.size() == waveforms.size())) ? baselines.at(idx).baseline() : Median(wave.Waveform());

          m_chargesum = std::accumulate(wave.Waveform().begin(), wave.Waveform().end(), 0,
                                        [&](short x, short y)
                                        { return ((m_baseline - x) + (m_baseline - y)); });

          // if required, save the full waveforms as well
          if (fSaveRawWaveforms)
            m_wf = wave.Waveform();

          fOpDetWaveformTrees[i]->Fill();
          idx++;
        }
      }
      else
      {
        mf::LogWarning("ICARUSFlashAssAna") << "Data product std::vector<raw::OpDetWaveform> for " << wflabel.label()
                                            << " is missing! Was it dropped before this stage?";
      }
    }
  }

  // -----
  // FLASHES/OPHITS INFO
  // Now we take care of the flashes:
  // we separate the case where we have a flash and the case where we don't have a flash
  // the difference is in how ophits are stored...

  if (!fFlashLabels.empty())
  {

    // hold the cryostat information
    std::vector<unsigned int> cids;

    for (std::size_t iFlashLabel = 0; iFlashLabel < fFlashLabels.size(); iFlashLabel++)
    {

      auto const label = fFlashLabels[iFlashLabel];
      auto const &flash_handle = e.getValidHandle<std::vector<recob::OpFlash>>(label);
      auto const &flashes = *flash_handle;

      // want our flashes to be valid and not empty
      if (flashes.empty())
      {
        mf::LogWarning("ICARUSFlashAssAna")
            << "No recob::OpFlash in collection with label '" << label.encode() << "'";
      }
      else
      {

        art::FindManyP<recob::OpHit> ophitsPtr(flash_handle, e, label);

        std::size_t idx = 0;
        for (auto const &flash : flashes)
        {

          m_pmt_pe.resize(360);
          m_pmt_start_time.resize(360);
          m_pmt_rise_time.resize(360);
          m_pmt_start_time_rwm.resize(360);
          m_pmt_amplitude.resize(360);

          m_flash_id = idx;
          m_flash_time = flash.Time();
          m_sum_pe = flash.TotalPE();
          m_flash_y = flash.YCenter();
          m_flash_width_y = flash.YWidth();
          m_flash_z = flash.ZCenter();
          m_flash_width_z = flash.ZWidth();

          auto const &ophits = ophitsPtr.at(idx);

          // we keep track of the cryistats where the flashes are found;
          geo::CryostatID::CryostatID_t cid = getCryostatByChannel(ophits.front()->OpChannel());
          auto const found = std::find(cids.begin(), cids.end(), cid);
          if (found != cids.end())
          {
            cids.push_back(cid);
          }

          // get the multiplicity, the number of PE per side, and the RWM-relative flash time
          // also store the first ophits for every channel in the flash
          processOpHitsFlash(ophits,
                             m_multiplicity_left, m_multiplicity_right,
                             m_sum_pe_left, m_sum_pe_right, m_pmt_start_time,
                             m_pmt_rise_time, m_pmt_start_time_rwm, m_pmt_pe, m_pmt_amplitude,
                             fOpHitFlashTrees[iFlashLabel]);

          m_multiplicity = m_multiplicity_left + m_multiplicity_right;

          // get the flash interaction time w.r.t. RWM
          // this is currently the mean between the first ophits on opposite walls
          m_flash_time_rwm = getFlashBunchTime(m_pmt_start_time_rwm, m_pmt_rise_time);

          fOpFlashTrees[iFlashLabel]->Fill();

          m_pmt_pe.clear();
          m_pmt_start_time.clear();
          m_pmt_rise_time.clear();
          m_pmt_start_time_rwm.clear();

          idx++;
        }
      }
    }

    // If the flashes did not cover all three cryostats..
    // ..well, we save the ophits on what is missing
    for (unsigned int cid = 0; cid < fGeom->Ncryostats(); cid++)
    {

      auto const found = std::find(cids.begin(), cids.end(), cid);
      if (found == cids.end())
      {
        processOpHits(e, cid);
      }
    }
  }

  else
  {
    mf::LogError("ICARUSFlashAssAna") << "No recob::OpFlash labels selected";

    // we save the ophits anyways even in absence of flashes
    for (unsigned int cid = 0; cid < fGeom->Ncryostats(); cid++)
    {
      processOpHits(e, cid);
    }
  }

  fEventTree->Fill();
}

DEFINE_ART_MODULE(opana::ICARUSFlashAssAna)
