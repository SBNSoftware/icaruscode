////////////////////////////////////////////////////////////////////////
// Class:       ICARUSOpticalDebug
// Plugin Type: analyzer (art v3_06_03)
// File:        ICARUSOpticalDebug_module.cc
//
// Module that dumps OpFlashes and their OpHit content.
//
// mailto:mvicenzi@bnl.gov
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

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "icaruscode/Timing/PMTTimingCorrections.h"
#include "icaruscode/Timing/IPMTTimingCorrectionService.h"

#include "TTree.h"

#include <array>
#include <vector>
#include <map>
#include <numeric> // std::accumulate
#include <limits>
#include <cstddef>

namespace opana
{
  class ICARUSOpticalDebug;
}

class opana::ICARUSOpticalDebug : public art::EDAnalyzer
{

public:
  struct Config
  {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Sequence<art::InputTag> OpHitLabels{
        Name("OpHitLabels"),
        Comment("Tags for the recob::OpHit data products")};

    fhicl::Sequence<art::InputTag> FlashLabels{
        Name("FlashLabels"),
        Comment("Tags for the recob::Flash data products")};

    fhicl::Atom<float> OpHitThresholdADC{
        Name("OpHitThresholdADC"),
        Comment("Threshold in ADC for an OpHit to be considered")};

    fhicl::Atom<std::string> LogCategory{
        fhicl::Name("LogCategory"),
        fhicl::Comment("category tag used for messages to message facility"),
        "OpHitTimingCorrection" // default
    };

  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit ICARUSOpticalDebug(Parameters const &config);

  ICARUSOpticalDebug(ICARUSOpticalDebug const &) = delete;
  ICARUSOpticalDebug(ICARUSOpticalDebug &&) = delete;
  ICARUSOpticalDebug &operator=(ICARUSOpticalDebug const &) = delete;
  ICARUSOpticalDebug &operator=(ICARUSOpticalDebug &&) = delete;

  void analyze(art::Event const &e) override;
  void beginJob() override;

  /// Return cryostat from PMT channel_id
  geo::CryostatID::CryostatID_t getCryostatByChannel(int channel);

  /// Return wall from PMT channel_id
  int getSideByChannel(const int channel);

  /// Return PMT position from channel_id
  std::array<double, 3> getChannelXYZ(const int channel);

  /// Return timing correction(s) per channel id
  double getTimingCorrection(const int channel);

  /// Process OpHits in the absence of flashes
  void processOpHits(art::Event const &e, unsigned int cryo);

  /// Process Ophits in the presence of flashes
  void processOpHitsFlash(std::vector<art::Ptr<recob::OpHit>> const &ophits,
                          int &multiplicity_left, int &multiplicity_right,
                          std::vector<int> &channel_id,
                          std::vector<double> &hit_x,
                          std::vector<double> &hit_y,
                          std::vector<double> &hit_z,
                          std::vector<float> &hit_start_time,
                          std::vector<float> &hit_rise_time,
                          std::vector<float> &hit_peak_time,
                          std::vector<double> &hit_timing_corr,
                          std::vector<double> &hit_area,
                          std::vector<double> &hit_pe,
                          std::vector<double> &hit_amplitude);

private:
  //----------
  // Input parameters
  //
  std::vector<art::InputTag> fOpHitLabels;
  std::vector<art::InputTag> fFlashLabels;

  //----------
  // Output trees

  std::vector<TTree *> fOpFlashTrees; // flash + matched hits
  std::vector<TTree *> fOpHitTrees;   // unmatched hits

  //----------------
  // Output variables

  // Common
  int m_run;
  int m_event;
  int m_timestamp;

  // Flash + hits tree
  int m_flash_id;
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
  std::vector<int> m_channel_id;
  std::vector<double> m_hit_x;
  std::vector<double> m_hit_y;
  std::vector<double> m_hit_z;
  std::vector<float> m_hit_start_time;
  std::vector<float> m_hit_rise_time;
  std::vector<float> m_hit_peak_time;
  std::vector<double> m_hit_timing_corr;
  std::vector<double> m_hit_area;
  std::vector<double> m_hit_pe;
  std::vector<double> m_hit_amplitude;

  // Ophit trees
  int m_channel;
  double m_x;
  double m_y;
  double m_z;
  double m_integral;  // in ADC x tick
  double m_amplitude; // in ADC
  float m_start_time;
  float m_peak_time;
  float m_rise_time;
  double m_timing_corr;
  float m_width;
  float m_abs_start_time;
  double m_pe;
  float m_fast_to_total;

  //----------
  // Support variables/products

  geo::GeometryCore const *fGeom;
  geo::WireReadoutGeom const *fChannelMapAlg;
  float const fOpHitThresholdADC;
  icarusDB::PMTTimingCorrections const &fPMTTimingCorrectionsService;
};

// ----------------------------------------------------------------------------

opana::ICARUSOpticalDebug::ICARUSOpticalDebug(Parameters const &config)
    : EDAnalyzer(config),
      fOpHitLabels(config().OpHitLabels()),
      fFlashLabels(config().FlashLabels()),
      fGeom(lar::providerFrom<geo::Geometry>()),
      fChannelMapAlg(&art::ServiceHandle<geo::WireReadout const>()->Get()),
      fOpHitThresholdADC(config().OpHitThresholdADC()),
      fPMTTimingCorrectionsService(*(lar::providerFrom<icarusDB::IPMTTimingCorrectionService const>()))
{
}

// ----------------------------------------------------------------------------

std::array<double, 3> opana::ICARUSOpticalDebug::getChannelXYZ(const int channel)
{
  auto const PMTxyz = fChannelMapAlg->OpDetGeoFromOpChannel(channel).GetCenter();
  return std::array<double, 3>{PMTxyz.X(), PMTxyz.Y(), PMTxyz.Z()};
}

// ----------------------------------------------------------------------------

double opana::ICARUSOpticalDebug::getTimingCorrection(const int channel)
{
  double laserTimeCorrection = fPMTTimingCorrectionsService.getLaserCorrections(channel);
  double cosmicsCorrection = fPMTTimingCorrectionsService.getCosmicsCorrections(channel);
  return laserTimeCorrection + cosmicsCorrection;
}

// ----------------------------------------------------------------------------


void opana::ICARUSOpticalDebug::beginJob()
{
  art::ServiceHandle<art::TFileService const> tfs;

  // Setting up the OpHit trees (one per cryostat)
  // This ttree will hold the unmatched ophit information
  for (auto const &label : fOpHitLabels)
  {
    std::string name = label.label() + "_ttree";
    std::string info = "TTree for unmatched recob::OpHit with label " + label.label();
    TTree *ttree = tfs->make<TTree>(name.c_str(), info.c_str());
    ttree->Branch("run", &m_run, "run/I");
    ttree->Branch("event", &m_event, "event/I");
    ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
    ttree->Branch("channel_id", &m_channel, "channel_id/I");
    ttree->Branch("integral", &m_integral, "integral/D");
    ttree->Branch("amplitude", &m_amplitude, "amplitude/D");
    ttree->Branch("start_time", &m_start_time, "start_time/F");
    ttree->Branch("peak_time", &m_peak_time, "peak_time/F");
    ttree->Branch("rise_time", &m_rise_time, "rise_time/F");
    ttree->Branch("abs_start_time", &m_abs_start_time, "abs_start_time/F");
    ttree->Branch("timing_corr", &m_timing_corr, "timing_corr/D");
    ttree->Branch("pe", &m_pe, "pe/D");
    ttree->Branch("width", &m_width, "width/F");
    ttree->Branch("x", &m_x, "x/D");
    ttree->Branch("y", &m_y, "y/D");
    ttree->Branch("z", &m_z, "z/D");
    ttree->Branch("fast_to_total", &m_fast_to_total, "fast_to_total/F");

    fOpHitTrees.push_back(ttree);
  }

  // Setting up the OPFLASH/OPHITS trees (one per cryostat)
  // These ttrees hold the information for matched ophits and their flashes
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
      ttree->Branch("flash_time", &m_flash_time, "flash_time/F");
      ttree->Branch("flash_y", &m_flash_y, "flash_y/F");
      ttree->Branch("flash_z", &m_flash_z, "flash_z/F");
      ttree->Branch("flash_width_y", &m_flash_width_y, "flash_width_y/F");
      ttree->Branch("flash_width_z", &m_flash_width_z, "flash_width_z/F");
      ttree->Branch("flash_nhits", &m_flash_nhits, "flash_nhits/I");
      ttree->Branch("channel_id", &m_channel_id);
      ttree->Branch("hit_x", &m_hit_x);
      ttree->Branch("hit_y", &m_hit_y);
      ttree->Branch("hit_z", &m_hit_z);
      ttree->Branch("hit_start_time", &m_hit_start_time);
      ttree->Branch("hit_rise_time", &m_hit_rise_time);
      ttree->Branch("hit_peak_time", &m_hit_peak_time);
      ttree->Branch("hit_timing_corr", &m_hit_timing_corr);
      ttree->Branch("hit_area", &m_hit_area);
      ttree->Branch("hit_pe", &m_hit_pe);
      ttree->Branch("hit_amplitude", &m_hit_amplitude);

      fOpFlashTrees.push_back(ttree);
    }
  }
}

// ----------------------------------------------------------------------------

geo::CryostatID::CryostatID_t opana::ICARUSOpticalDebug::getCryostatByChannel(int channel)
{
  const geo::OpDetGeo &opdetgeo = fChannelMapAlg->OpDetGeoFromOpChannel(channel);
  geo::CryostatID::CryostatID_t cid = opdetgeo.ID().Cryostat;
  return cid;
}

// ----------------------------------------------------------------------------

int opana::ICARUSOpticalDebug::getSideByChannel(const int channel)
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

void opana::ICARUSOpticalDebug::processOpHitsFlash(std::vector<art::Ptr<recob::OpHit>> const &ophits,
                                                 int &multiplicity_left, int &multiplicity_right,
                                                 std::vector<int> &channel_id,
                                                 std::vector<double> &hit_x,
                                                 std::vector<double> &hit_y,
                                                 std::vector<double> &hit_z,
                                                 std::vector<float> &hit_start_time,
                                                 std::vector<float> &hit_rise_time,
                                                 std::vector<float> &hit_peak_time,
                                                 std::vector<double> &hit_timing_corr,
                                                 std::vector<double> &hit_area,
                                                 std::vector<double> &hit_pe,
                                                 std::vector<double> &hit_amplitude)
{

  std::unordered_map<int, float> sumpe_map;
  for (auto const ophit : ophits)
  {

    // check if above threshold to be printed
    if (ophit->Amplitude() < fOpHitThresholdADC)
      continue;

    const int ch = ophit->OpChannel();
    channel_id.push_back(ch);

    auto pos = getChannelXYZ(ch);
    hit_x.push_back(pos[0]);
    hit_y.push_back(pos[1]);
    hit_z.push_back(pos[2]);

    // save times: start, peak, rise
    // rise is relative to start, make it absolute
    hit_start_time.push_back(ophit->StartTime());
    hit_peak_time.push_back(ophit->PeakTime());
    hit_rise_time.push_back(ophit->StartTime() + ophit->RiseTime());
    hit_timing_corr.push_back( getTimingCorrection(ch) );

    hit_area.push_back(ophit->Area());           // in ADC x tick
    hit_amplitude.push_back(ophit->Amplitude()); // in ADC
    hit_pe.push_back(ophit->PE());

    sumpe_map[ch] += ophit->PE();
  }

  m_multiplicity_left = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0,
                                        [&](int value, const std::map<int, float>::value_type &p)
                                        { return getSideByChannel(p.first) == 0 ? ++value : value; });

  m_multiplicity_right = std::accumulate(sumpe_map.begin(), sumpe_map.end(), 0,
                                         [&](int value, const std::map<int, float>::value_type &p)
                                         { return getSideByChannel(p.first) == 1 ? ++value : value; });
}

// ----------------------------------------------------------------------------

void opana::ICARUSOpticalDebug::analyze(art::Event const &e)
{

  // Collect global event metadata
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second

  // hold the cryostat information
  std::vector<unsigned int> cids;

  // -----
  // FLASHES/OPHITS INFO
  // Now we take care of the flashes:

  if (!fFlashLabels.empty())
  {

    for (std::size_t iFlashLabel = 0; iFlashLabel < fFlashLabels.size(); iFlashLabel++)
    {
      auto const label = fFlashLabels[iFlashLabel];
      auto const &flash_handle = e.getValidHandle<std::vector<recob::OpFlash>>(label);
      auto const &flashes = *flash_handle;

      // want our flashes to be valid and not empty
      if (flashes.empty())
      {
        mf::LogWarning("ICARUSOpticalDebug")
            << "No recob::OpFlash in collection with label '" << label.encode() << "'";
      }
      else
      {
        art::FindManyP<recob::OpHit> ophitsPtr(flash_handle, e, label);
        std::size_t idx = 0;
        for (auto const &flash : flashes)
        {
          m_channel_id.clear();
          m_hit_x.clear();
          m_hit_y.clear();
          m_hit_z.clear();
          m_hit_start_time.clear();
          m_hit_rise_time.clear();
          m_hit_peak_time.clear();
          m_hit_timing_corr.clear();
          m_hit_area.clear();
          m_hit_pe.clear();
          m_hit_amplitude.clear();

          m_flash_id = idx;
          m_flash_time = flash.Time();
          m_sum_pe = flash.TotalPE();
          m_flash_z = flash.ZCenter();
          m_flash_y = flash.YCenter();
          m_flash_width_y = flash.YWidth();
          m_flash_width_z = flash.ZWidth();

          auto const &ophits = ophitsPtr.at(idx);

          // we keep track of the cryostats where the flashes are found;
          geo::CryostatID::CryostatID_t cid = getCryostatByChannel(ophits.front()->OpChannel());
          auto const found = std::find(cids.begin(), cids.end(), cid);
          if (found == cids.end()) // if seen for the first time
            cids.push_back(cid);   // save it!

          // also store the matched ophits in the flash
          processOpHitsFlash(ophits,
                             m_multiplicity_left, m_multiplicity_right,
                             m_channel_id, m_hit_x, m_hit_y, m_hit_z,
                             m_hit_start_time, m_hit_rise_time, m_hit_peak_time,
                             m_hit_timing_corr,
                             m_hit_area, m_hit_pe, m_hit_amplitude);

          m_multiplicity = m_multiplicity_left + m_multiplicity_right;
          m_flash_nhits = m_channel_id.size();

          mf::LogInfo("ICARUSOpticalDebug") << "Cryo " << cid << ", flashID " << m_flash_id
                                          << ", nophits " << m_flash_nhits;

          fOpFlashTrees[iFlashLabel]->Fill();
          idx++;
        } // flash loop
      }
    } // both flash products

    // -----
    // UNMATCHED OPHITS INFO
    // Now we take care of the hits that didn't match with any flash.
    // But only from cryostats that had a flash --> if no flashes at all, don't save anything

    if (fOpHitLabels.empty())
    {
      mf::LogError("ICARUSOpticalDebug") << "No recob::OpHit labels selected.";
      return;
    }

    for (std::size_t iOpHitLabel = 0; iOpHitLabel < fOpHitLabels.size(); iOpHitLabel++)
    {

      auto const label = fOpHitLabels[iOpHitLabel];
      auto const &ophits_handle = e.getValidHandle<std::vector<recob::OpHit>>(label);
      auto const &ophits = *ophits_handle;

      // we want our ophits to be valid and not empty
      if (ophits.empty())
      {
        mf::LogError("ICARUSOpticalDebug") << "Invalid recob::OpHit with label '" << label.encode() << "'";
        continue;
      }

      unsigned int matched = 0;
      unsigned int unmatched = 0;

      // for the cryos that saw flashes
      for (auto cid : cids)
      {
        auto const assn_label = fFlashLabels[cid];
        art::FindOneP<recob::OpFlash> flashPtr(ophits_handle, e, assn_label);

        // loop over this hits --> if the single hit has an associated flash, skip it
        for (std::size_t idx = 0; idx < ophits.size(); ++idx)
        {
          art::Ptr<recob::OpFlash> flash_ptr = flashPtr.at(idx);

          // Check if this OpHit has an associated OpFlash
          if (flash_ptr.isNull()) // if unmatched, store it
          {
            auto const &ophit = ophits[idx];

            // check if above threshold to be printed
            if (ophit.Amplitude() < fOpHitThresholdADC)
              continue;

            m_channel = ophit.OpChannel();

            auto pos = getChannelXYZ(ophit.OpChannel());
            m_x = pos[0];
            m_y = pos[1];
            m_z = pos[2];
            m_integral = ophit.Area();
            m_amplitude = ophit.Amplitude();
            m_start_time = ophit.StartTime();
            m_peak_time = ophit.PeakTime();
            m_rise_time = ophit.StartTime() + ophit.RiseTime();
            m_width = ophit.Width();
            m_abs_start_time = ophit.PeakTimeAbs() + (m_start_time - m_peak_time);
            m_timing_corr = getTimingCorrection(ophit.OpChannel());
            m_pe = ophit.PE();
            m_fast_to_total = ophit.FastToTotal();

            fOpHitTrees[iOpHitLabel]->Fill();
            unmatched++;
          }
          else
          { // if matched, already stored with its flash
            matched++;
          }
        } // for each ophit
      } // for each cryo/flashlabel

      mf::LogInfo("ICARUSOpticalDebug")
          << "Total flash-matched ophits: " << matched << "\n"
          << "Total unmatched ophits: " << unmatched;

    } // for ophit product
  }
  else
  {
    mf::LogError("ICARUSOpticalDebug") << "No recob::OpFlash labels selected";
  }
}

DEFINE_ART_MODULE(opana::ICARUSOpticalDebug)
