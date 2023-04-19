////////////////////////////////////////////////////////////////////////
// Class:       FilterCRTPMTMatching
// Plugin Type: analyzer (Unknown Unknown)
// File:        FilterCRTPMTMatching_module.cc
//
// Generated at Thu Feb 17 13:31:00 2022 by Biswaranjan Behera using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Data product includes

#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// C++ includes
#include <map>
#include <numeric>
#include <vector>

// ROOT includes
#include "TTree.h"
#include "TVector3.h"

using std::map;
using std::vector;

namespace icarus::crt {
class FilterCRTPMTMatching;
struct CRTPMT {
  double tof;    ///< Time difference between CRT Hit and optical flash [ns]
  double distance;    ///< Distance between CRT Hit and optical flash centroid [cm]
  art::Ptr<sbn::crt::CRTHit> CRTHit;
};
enum hasCRTHit {
  /**
   * Type of the matching.
   *
   * Secret decoder ring:
   *  * `0`: No CRT match
   *  * `1`: matched with Top CRT hit before optical flash
   *  * `2`: matched with Side CRT hit before optical flash
   *  * `3`: entering from Top and exiting from side
   *  * `4`: matched with a Top CRT after the optical flash
   *  * `5`: matched with a Side CRT after the optical flash
   *  * `6`: matched with multiple Top CRT hits before the optical flash
   *  * `7`: matched with multiple Top CRT hits before the optical flash and
   * more then 1 side CRT hits after the optical flash
   *  * `8`: all the other cases
   *
   */
  noMatch = 0,
  enTop = 1,
  enSide = 2,
  enTop_exSide = 3,
  exTop = 4,
  exSide = 5,
  enTop_mult = 6,
  enTop_exSide_mult = 7,
  others = 9
};
struct CRTMatches {
  // vector of pairs where first is the matched Time of Flight and second is the
  // matched CRTHit
  std::vector<CRTPMT> entering;
  std::vector<CRTPMT> exiting;
  hasCRTHit matchType;
};
struct MatchedCRT {
  geo::Point_t CRTHitPos;
  double CRTPMTTimeDiff_ns;
  double CRTTime_us;
  int CRTSys;
  int CRTRegion;
};
struct FlashType {
  geo::Point_t FlashPos;
  double FlashTime_us;
  double FlashGateTime_ns;
  bool inBeam;
  bool inGate;
  hasCRTHit Classification;
  std::vector<MatchedCRT> CRTmatches;
};
struct EventCRTPMT {
  bool CosmicLike;
  std::vector<FlashType> inGateFlashes;
};
}  // namespace icarus::crt

using namespace icarus::crt;

bool flashInTime(double flashTime, int gateType, double gateDiff,
                 double gateWidth) {
  //   {1, 2.2},  // BNB
  //   {2, 10.1}, // NuMI
  //   {3, 2.2},  // BNB Offbeam
  //   {4, 10.1}, // NuMI Offbeam

  // As a reminder, I will leave here the commented part of the vetoOffset, in
  // case something changes in the future
  /*double vetoOffset = 3.5;*/
  double activeGate = gateWidth /*- vetoOffset*/;

  double relFlashTime = flashTime + gateDiff / 1000. /*- vetoOffset*/;
  mf::LogDebug("FilterCRTPMTMatching FlashInTime")
      << "Gate Diff " << gateDiff / 1000 << " Ftime+gateDiff "
      << flashTime + gateDiff / 1000. << " " << activeGate;

  return ((relFlashTime > 0) && (relFlashTime < activeGate));
}

icarus::crt::CRTMatches CRTHitmatched(
    double flashTime, geo::Point_t const& flashpos,
    std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, double interval) {
  std::vector<icarus::crt::CRTPMT> enteringCRTHits;
  std::vector<icarus::crt::CRTPMT> exitingCRTHits;
  hasCRTHit MatchType;
  int topen = 0, topex = 0, sideen = 0, sideex = 0;
  for (auto const& crtHit : crtHits) {
    double tof = crtHit->ts1_ns - flashTime * 1e3;
    double distance =
        (flashpos - geo::Point_t{crtHit->x_pos, crtHit->y_pos, crtHit->z_pos})
            .R();
    if (abs(tof) >= interval) continue;
    if (tof < 0) {
      if (crtHit->plane > 36)
        sideen++;
      else
        topen++;
      CRTPMT m_match = {tof, distance, crtHit};
      enteringCRTHits.push_back(m_match);
    } else if (tof >= 0) {
      if (crtHit->plane > 36)
        sideex++;
      else
        topex++;
      CRTPMT m_match = {tof, distance, crtHit};
      exitingCRTHits.push_back(m_match);
    }
  }
  if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0)
    MatchType = noMatch;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0)
    MatchType = enTop;
  else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0)
    MatchType = enSide;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1)
    MatchType = enTop_exSide;
  else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0)
    MatchType = exTop;
  else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1)
    MatchType = exSide;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0)
    MatchType = enTop_mult;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1)
    MatchType = enTop_exSide_mult;
  else
    MatchType = others;

  return {std::move(enteringCRTHits), std::move(exitingCRTHits), MatchType};
}

class icarus::crt::FilterCRTPMTMatching : public art::EDFilter {

/**
 * @brief Rejects events with only incoming particles in time with the beam.
 *
 * Here is a long and detailed explanation on what are the choices.
 */

 public:
  using CRTHit = sbn::crt::CRTHit;

  explicit FilterCRTPMTMatching(fhicl::ParameterSet const& p);

  // Required functions.
  //void getTriggerConf(art::Run const&);
  void beginRun(art::Run const& run);
  bool filter(art::Event&) override;

 private:
  // Declare member data here.

  static bool HitCompare(const art::Ptr<CRTHit>& h1,
                         const art::Ptr<CRTHit>& h2);
  void ClearVecs();
  std::vector<art::InputTag> fFlashLabels;

  art::InputTag fOpFlashModuleLabel0;
  art::InputTag fOpFlashModuleLabel1;
  art::InputTag fOpFlashModuleLabel2;
  art::InputTag fOpFlashModuleLabel3;
  art::InputTag fCrtHitModuleLabel;
  art::InputTag fTriggerLabel;
  art::InputTag fTriggerConfigurationLabel;

  icarus::TriggerConfiguration fTriggerConfiguration;

  // double                     fFlashTimeCut;

  int fEvent;   ///< number of the event being processed.
  int fRun;     ///< number of the run being processed.
  int fSubRun;  ///< number of the sub-run being processed.

  // add trigger data product vars
  unsigned int m_gate_type;
  std::string m_gate_name;
  uint64_t m_trigger_timestamp;
  uint64_t m_gate_start_timestamp;
  uint64_t m_trigger_gate_diff;
  uint64_t m_gate_width;

  std::string fFilterLevel;  // Filter level, default is loose
  int fnOpHitToTrigger;  // Number of OpHit above threshold to mimic the trigger
  double fTimeOfFlightInterval;  // Time of Flight interval to find the match
  bool fOutputTree;              // Output tree or not
  int fPMTADCThresh;
  std::vector<int> fTopBefore;
  std::vector<int> fTopAfter;
  std::vector<int> fSideBefore;
  std::vector<int> fSideAfter;
  double fBNBBeamGateMin;
  double fBNBBeamGateMax;
  double fBNBinBeamMin;
  double fBNBinBeamMax;
  double fNuMIBeamGateMin;
  double fNuMIBeamGateMax;
  double fNuMIinBeamMin;
  double fNuMIinBeamMax;

  TTree* fMatchTree;

  int fFiltered;  ///< Was the event filtered out?

  std::vector<int> fHitsInTime;  ///< Classification of the Hits in Time

  // matchTree vars

  double fOpFlashTime_us;
  double fOpFlashPos_X;
  double fOpFlashPos_Y;
  double fOpFlashPos_Z;
  int fClassification;
  bool finGate;
  bool finBeam;
  double fFlashBeamTime_ns;
  vector<double> fCRTHitPos_X;
  vector<double> fCRTHitPos_Y;
  vector<double> fCRTHitPos_Z;
  vector<int> fCRTHitRegion;
  vector<int> fCRTHitSystem;
  vector<double> fCRTHitTime_us;
  vector<double> fCRTFlashTime_ns;

  geo::GeometryCore const* fGeometryService;  ///< pointer to Geometry provider.
};

icarus::crt::FilterCRTPMTMatching::FilterCRTPMTMatching(
    fhicl::ParameterSet const& p)
    : EDFilter{p},
      fOpFlashModuleLabel0(
          p.get<art::InputTag>("OpFlashModuleLabel0")),
      fOpFlashModuleLabel1(
          p.get<art::InputTag>("OpFlashModuleLabel1")),
      fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit")),
      fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")),
      fTriggerConfigurationLabel(
          p.get<art::InputTag>("TriggerConfiguration", "triggerconfig")),
      fFilterLevel(p.get<std::string>("FilterLevel", "loose")),
      fnOpHitToTrigger(p.get<int>("nOpHitToTrigger")),
      fTimeOfFlightInterval(p.get<double>("TimeOfFlightInterval")),
      fPMTADCThresh(p.get<int>("PMTADCThresh")),
      fBNBBeamGateMin(p.get<double>("BNBBeamGateMin")),
      fBNBBeamGateMax(p.get<double>("BNBBeamGateMax")),
      fBNBinBeamMin(p.get<double>("BNBinBeamMin")),
      fBNBinBeamMax(p.get<double>("BNBinBeamMax")),
      fNuMIBeamGateMin(p.get<double>("NuMIBeamGateMin")),
      fNuMIBeamGateMax(p.get<double>("NuMIBeamGateMax")),
      fNuMIinBeamMin(p.get<double>("NuMIinBeamMin")),
      fNuMIinBeamMax(p.get<double>("NuMIinBeamMax")),
      fGeometryService(lar::providerFrom<geo::Geometry>()) {
  fFlashLabels = { fOpFlashModuleLabel0, fOpFlashModuleLabel1 };
  art::ServiceHandle<art::TFileService> tfs;
  if (fOutputTree == true) {
    fMatchTree =
        tfs->make<TTree>("matchTree", "CRTHit - OpHit/Flash matching analysis");
    fMatchTree->Branch("event", &fEvent, "event/I");
    fMatchTree->Branch("run", &fRun, "run/I");
    fMatchTree->Branch("subrun", &fSubRun, "subrun/I");
    fMatchTree->Branch("Filter", &fFiltered);
    fMatchTree->Branch("FlashClassification", &fClassification);
    fMatchTree->Branch("FlashBarycenter_X", &fOpFlashPos_X);
    fMatchTree->Branch("FlashBarycenter_Y", &fOpFlashPos_Y);
    fMatchTree->Branch("FlashBarycenter_Z", &fOpFlashPos_Z);
    fMatchTree->Branch("FlashTime_us", &fOpFlashTime_us);
    fMatchTree->Branch("inGate", &finGate);
    fMatchTree->Branch("inBeam", &finBeam);
    fMatchTree->Branch("FlashGateTime", &fFlashBeamTime_ns);
    fMatchTree->Branch("CRTHit_X", &fCRTHitPos_X);
    fMatchTree->Branch("CRTHit_Y", &fCRTHitPos_Y);
    fMatchTree->Branch("CRTHit_Z", &fCRTHitPos_Z);
    fMatchTree->Branch("CRTHitTime_us", &fCRTHitTime_us);
    fMatchTree->Branch("CRTFlashTimeDiff_ns", &fCRTFlashTime_ns);
    fMatchTree->Branch("CRTHitRegion", &fCRTHitRegion);
    fMatchTree->Branch("CRTHitDetector", &fCRTHitSystem);

    fMatchTree->Branch("gate_type", &m_gate_type, "gate_type/b");
    fMatchTree->Branch("gate_name", &m_gate_name);
    fMatchTree->Branch("trigger_timestamp", &m_trigger_timestamp,
                       "trigger_timestamp/l");
    fMatchTree->Branch("gate_start_timestamp", &m_gate_start_timestamp,
                       "gate_start_timestamp/l");
    fMatchTree->Branch("trigger_gate_diff", &m_trigger_gate_diff,
                       "trigger_gate_diff/l");
  }
}

//void icarus::crt::FilterCRTPMTMatching::getTriggerConf(art::Run const& r) {
//  fTriggerConfiguration =
//      r.getProduct<icarus::TriggerConfiguration>(fTriggerConfigurationLabel);
//}

void icarus::crt::FilterCRTPMTMatching::beginRun(art::Run const& r) {
  fTriggerConfiguration =
      r.getProduct<icarus::TriggerConfiguration>(fTriggerConfigurationLabel);
}

bool icarus::crt::FilterCRTPMTMatching::filter(art::Event& e) {
  mf::LogDebug("FilterCRTPMTMatching: ") << "beginning analyis";

  // Start by fetching some basic event information for our n-tuple.
  fEvent = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();

  ClearVecs();
  // add trigger info
  auto const& triggerInfo = e.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel);
  sbn::triggerSource bit = triggerInfo.sourceType;
  m_gate_type = (unsigned int)bit;
  m_gate_name = bitName(bit);
  m_trigger_timestamp = triggerInfo.triggerTimestamp;
  m_gate_start_timestamp = triggerInfo.beamGateTimestamp;
  m_trigger_gate_diff = triggerInfo.triggerTimestamp - triggerInfo.beamGateTimestamp;
  m_gate_width = fTriggerConfiguration.getGateWidth(m_gate_type);

  // CRTHits
  art::Handle<std::vector<CRTHit>> crtHitListHandle;
  std::vector<art::Ptr<CRTHit>> crtHitList;
  if (e.getByLabel(fCrtHitModuleLabel, crtHitListHandle))
    art::fill_ptr_vector(crtHitList, crtHitListHandle);
  if ((fFilterLevel != "loose") && (fFilterLevel != "medium") && (fFilterLevel != "tight"))
    throw art::Exception{ art::errors::Configuration } << "Invalid CRT/PMT filter level: '" << fFilterLevel << "'\n";

  mf::LogInfo("FilterCRTPMTMatching::FilteringLevel ") << fFilterLevel;
  std::vector<FlashType> thisEventFlashes;
  for (art::InputTag const& flashLabel : fFlashLabels) {
    auto const flashHandle =
        e.getHandle<std::vector<recob::OpFlash>>(flashLabel);
    art::FindMany<recob::OpHit> findManyHits(flashHandle, e, flashLabel);

    for (auto const& [iflash, flash] : util::enumerate(*flashHandle)) {
      hasCRTHit eventType = others;
      double tflash = flash.Time();
      vector<recob::OpHit const*> const& hits = findManyHits.at(iflash);
      int nPMTsTriggering = 0;
      double firstTime = 999999;
      geo::vect::MiddlePointAccumulator flashCentroid;
      double ampsum = 0, t_m = 0;
      std::vector<double> fHitX, fHitY, fHitZ, fHitT, fHitA;
      for (auto const& hit : hits) {
        if (hit->Amplitude() > fPMTADCThresh) nPMTsTriggering++;
        if (firstTime > hit->PeakTime()) firstTime = hit->PeakTime();
        geo::Point_t const pos =
            fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel())
                .GetCenter();
        double amp = hit->Amplitude();
        ampsum += amp;
        fHitX.push_back(pos.X());
        fHitY.push_back(pos.Y());
        fHitZ.push_back(pos.Z());
        fHitT.push_back(hit->PeakTime());
        fHitA.push_back(amp);
        flashCentroid.add(pos, amp);
        t_m = t_m + hit->PeakTime();
      }
      geo::Point_t flash_pos = flashCentroid.middlePoint();
      t_m = t_m / nPMTsTriggering;
      if (nPMTsTriggering < fnOpHitToTrigger) {
        continue;
      }
      bool inTime = flashInTime(firstTime, m_gate_type, m_trigger_gate_diff,
                                m_gate_width);
      double fThisRelGateTime = m_trigger_gate_diff + tflash * 1e3;
      bool fThisInTime_gate = false;
      bool fThisInTime_beam = false;
      if (m_gate_type == 1 || m_gate_type == 3) {  // BNB OffBeamBNB
        if (fThisRelGateTime > fBNBBeamGateMin &&
            fThisRelGateTime < fBNBBeamGateMax)
          fThisInTime_gate = true;
        if (fThisRelGateTime > fBNBinBeamMin &&
            fThisRelGateTime < fBNBinBeamMax)
          fThisInTime_beam = true;
      }
      if (m_gate_type == 2 || m_gate_type == 4) {  // NuMI OffBeamNuMI
        if (fThisRelGateTime > fNuMIBeamGateMin &&
            fThisRelGateTime < fNuMIBeamGateMax)
          fThisInTime_gate = true;
        if (fThisRelGateTime > fNuMIinBeamMin &&
            fThisRelGateTime < fNuMIinBeamMax)
          fThisInTime_beam = true;
      }
      inTime = fThisInTime_gate;
      icarus::crt::CRTMatches CRTmatches = CRTHitmatched(
          firstTime, flash_pos, crtHitList, fTimeOfFlightInterval);
      int TopEn = 0, TopEx = 0, SideEn = 0, SideEx = 0;
      auto nCRTHits = CRTmatches.entering.size() + CRTmatches.exiting.size();
      std::vector<MatchedCRT> thisFlashCRTmatches;
      eventType = CRTmatches.matchType;
      if (nCRTHits > 0) {
        for (auto const& entering : CRTmatches.entering) {
          vector<double> CRTpos {entering.CRTHit->x_pos,
                                 entering.CRTHit->y_pos,
                                 entering.CRTHit->z_pos};
          geo::Point_t thisCRTpos {entering.CRTHit->x_pos,
                                   entering.CRTHit->y_pos,
                                   entering.CRTHit->z_pos};
          double CRTtime = entering.CRTHit->ts1_ns / 1e3;
          int CRTRegion = entering.CRTHit->plane;
          int CRTSys = 0;
          if (CRTRegion >= 36)
            CRTSys =
                1;  // Very lazy way to determine if the Hit is a Top or a Side.
                    // Will update it when bottom CRT will be availble.
          double CRTTof_opflash = CRTtime - tflash;
          std::vector<int> HitFebs;
          for (auto crts : entering.CRTHit->feb_id) {
            HitFebs.emplace_back((int)crts);
          }
          if (CRTSys == 0) TopEn++;
          if (CRTSys == 1) SideEn++;
	  MatchedCRT thisCRTMatch = { /* .CRTHitPos = */ thisCRTpos, // C++20: restore initializers
                                     /* .CRTPMTTimeDiff_ns = */ CRTTof_opflash,
                                     /* .CRTTime_us = */ CRTtime,
                                     /* .CRTSys = */ CRTSys,
                                     /* .CRTRegion = */ CRTRegion};
          thisFlashCRTmatches.push_back(thisCRTMatch);
        }
        for (auto const& exiting : CRTmatches.exiting) {
          vector<double> CRTpos {exiting.CRTHit->x_pos,
				 exiting.CRTHit->y_pos,
                                 exiting.CRTHit->z_pos};
          geo::Point_t thisCRTpos {exiting.CRTHit->x_pos,
                                   exiting.CRTHit->y_pos,
                                   exiting.CRTHit->z_pos};
          double CRTtime = exiting.CRTHit->ts1_ns / 1e3;
          int CRTRegion = exiting.CRTHit->plane;
          int CRTSys = 0;
          if (CRTRegion >= 36) CRTSys = 1;
          double CRTTof_opflash = CRTtime - tflash;
          std::vector<int> HitsFebs;
          for (auto crts : exiting.CRTHit->feb_id) {
            HitsFebs.emplace_back((int)crts);
          }
          if (CRTSys == 0) TopEx++;
          if (CRTSys == 1) SideEx++;
          MatchedCRT thisCRTMatch = { /* .CRTHitPos = */ thisCRTpos, // C++20: restore initializers
                                     /* .CRTPMTTimeDiff_ns = */ CRTTof_opflash,
                                     /* .CRTTime_us = */ CRTtime,
                                     /* .CRTSys = */ CRTSys,
                                     /* .CRTRegion = */ CRTRegion};
          thisFlashCRTmatches.push_back(thisCRTMatch);
        }
      }
      FlashType thisFlashType = { /* .FlashPos = */ flash_pos, // C++20: restore initializers
                                 /* .FlashTime_us = */ tflash,
                                 /* .FlashGateTime_ns = */ fThisRelGateTime,
                                 /* .inBeam = */ fThisInTime_beam,
                                 /* .inGate = */ inTime,
                                 /* .Classification = */ eventType,
                                 /* .CRTmatches = */ thisFlashCRTmatches};
      if (inTime == true) thisEventFlashes.push_back(thisFlashType);
    }
  }
  bool hasOnlyCosmics = false;
  if (fFilterLevel == "loose")
    hasOnlyCosmics = false;  // By default, with loose Filtering, everything is
                         // "intersting", nothing tagged as clear cosmic
  else if (fFilterLevel == "medium") {
    hasOnlyCosmics = true;
    for (const auto& h : thisEventFlashes) {
      bool isCosmic = false;
      if (h.Classification == enTop || h.Classification == enTop_exSide ||
          h.Classification == enTop_mult ||
          h.Classification == enTop_exSide_mult)
        isCosmic = true;
      // With Medium filter, everything (inTime) which is associated with Top
      // CRT Hit before the Flash is filtered as clear cosmic
      else
        isCosmic = false;
      hasOnlyCosmics = hasOnlyCosmics && isCosmic;
    }
  } else if (fFilterLevel == "tight") {
    // With Tight filter, everything (inTime) associated with a CRT Hit before
    // the Flash is filtered as clear cosmic
    hasOnlyCosmics = true;
    for (const auto& h : thisEventFlashes) {
      bool isCosmic = false;
      if (h.Classification != noMatch || h.Classification != exTop ||
          h.Classification != exSide)
        isCosmic = true;
      else
        isCosmic = false;
      hasOnlyCosmics = hasOnlyCosmics && isCosmic;
    }
  }
  fFiltered = hasOnlyCosmics;
  EventCRTPMT thisEvent = {/* .Filter = */ hasOnlyCosmics, // C++20: restore initializers
                           /* .inGateFlashes = */ thisEventFlashes};
  for (const auto& f : thisEventFlashes) {
    fClassification = f.Classification;
    fOpFlashPos_X = f.FlashPos.X();
    fOpFlashPos_Y = f.FlashPos.Y();
    fOpFlashPos_Z = f.FlashPos.Z();
    fOpFlashTime_us = f.FlashTime_us;
    finGate = f.inGate;
    finBeam = f.inBeam;
    fFlashBeamTime_ns = f.FlashGateTime_ns;
    for (const auto& crt : f.CRTmatches) {
      fCRTHitPos_X.push_back(crt.CRTHitPos.X());
      fCRTHitPos_Y.push_back(crt.CRTHitPos.Y());
      fCRTHitPos_Z.push_back(crt.CRTHitPos.Z());
      fCRTHitTime_us.push_back(crt.CRTTime_us);
      fCRTFlashTime_ns.push_back(crt.CRTPMTTimeDiff_ns);
      fCRTHitRegion.push_back(crt.CRTRegion);
      fCRTHitSystem.push_back(crt.CRTSys);
    }
    fMatchTree->Fill();
    ClearVecs();
  }
  if (thisEventFlashes.empty()) {
    fClassification = 9;
    fOpFlashPos_X = 0;
    fOpFlashPos_Y = 0;
    fOpFlashPos_Z = 0;
    fOpFlashTime_us = 0;
    finGate = false;
    finBeam = false;
    fFlashBeamTime_ns = 0;
    fMatchTree->Fill();
  }
  return !hasOnlyCosmics;
}

void icarus::crt::FilterCRTPMTMatching::ClearVecs() {
  // matchTree
  fCRTHitPos_X.clear();
  fCRTHitPos_Y.clear();
  fCRTHitPos_Z.clear();
  fCRTHitTime_us.clear();
  fCRTFlashTime_ns.clear();
  fCRTHitRegion.clear();
  fCRTHitSystem.clear();
}

DEFINE_ART_MODULE(icarus::crt::FilterCRTPMTMatching)
