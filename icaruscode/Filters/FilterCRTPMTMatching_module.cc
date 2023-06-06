/**
  * @file icaruscode/Filters/FilterCRTPMTMatching_module.cc
  * @author Francesco Poppi, mail: poppi@bo.infn.it
*/


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
//#include "icaruscode/IcarusObj/CRTPMTMatching.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"

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
}

class icarus::crt::FilterCRTPMTMatching : public art::EDFilter {

/**
 * @brief Rejects events with only incoming particles in time with the beam.
 *
 * This filtering module is based on CRT the CRT PMT matching.
 * It starts by considering only the flashes which are within the beam gate
 * and it tries to match them with one or more CRT Hits using a configurable
 * time interval, which, at the time (19/04/2023) is chosen as +/-100 ns.
 * If the flash is matched, the relative time will determine if the mu candidate
 * is entering or exiting the TPC. Depending on the amound and relative time of 
 * the match, a classification is provided to the flash.
 * The filtering module has 3 options: loose (everything goes through), medium
 * (removes events where the flash is matched with one Top CRT hit before the
 * flash) and tight (removes all the events where the flashes are matched with
 * CRT hits before the Flash).
 * If multiple flashes are in the beam gate, the filtering logic applies to the 
 * AND of the flashes classification.
 * For questions and maintenance ask Francesco Poppi.
 */

 public:
  using CRTHit = sbn::crt::CRTHit;
  using CRTPMTMatching = sbn::crt::CRTPMTMatching;
  explicit FilterCRTPMTMatching(fhicl::ParameterSet const& p);
  // Required functions.
  bool beginRun(art::Run& run) override;
  bool filter(art::Event&) override;

 private:
  // Declare member data here.

  void ClearVecs();

  art::InputTag fCrtPmtModuleLabel;
  art::InputTag fTriggerLabel;
  art::InputTag fTriggerConfigurationLabel;

  icarus::TriggerConfiguration fTriggerConfiguration;

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
  bool fOutputTree;              // Output tree or not
  bool fSpillOnly;	///< Consider only flashes in Spill

  TTree* fMatchTree = nullptr;

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

};

icarus::crt::FilterCRTPMTMatching::FilterCRTPMTMatching(
    fhicl::ParameterSet const& p)
    : EDFilter{p},
      fCrtPmtModuleLabel(p.get<art::InputTag>("CrtPmtModuleLabel")),
      fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")),
      fTriggerConfigurationLabel(
          p.get<art::InputTag>("TriggerConfiguration", "triggerconfig")),
      fFilterLevel(p.get<std::string>("FilterLevel", "loose")),
      fOutputTree(p.get<bool>("MakeMatchTree", true)),
      fSpillOnly(p.get<bool>("SpillOnly", true)){
  if (fOutputTree) {
    art::ServiceHandle<art::TFileService> tfs;
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

bool icarus::crt::FilterCRTPMTMatching::beginRun(art::Run& r) {
  fTriggerConfiguration =
      r.getProduct<icarus::TriggerConfiguration>(fTriggerConfigurationLabel);
  return true;
}

bool icarus::crt::FilterCRTPMTMatching::filter(art::Event& e) {
  mf::LogDebug("FilterCRTPMTMatching: ") << "beginning analyis";

  fEvent = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();
  ClearVecs();

  auto const& triggerInfo = e.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel);
  sbn::triggerSource bit = triggerInfo.sourceType;
  m_gate_type = (unsigned int)bit;
  m_gate_name = bitName(bit);
  m_trigger_timestamp = triggerInfo.triggerTimestamp;
  m_gate_start_timestamp = triggerInfo.beamGateTimestamp;
  m_trigger_gate_diff = triggerInfo.triggerTimestamp - triggerInfo.beamGateTimestamp;
  m_gate_width = fTriggerConfiguration.getGateWidth(m_gate_type);

  auto const& crtpmtMatches = e.getProduct<std::vector<CRTPMTMatching>>(fCrtPmtModuleLabel);

  if ((fFilterLevel != "loose") && (fFilterLevel != "medium") && (fFilterLevel != "tight"))
    throw art::Exception{ art::errors::Configuration } << "Invalid CRT/PMT filter level: '" << fFilterLevel << "'\n";

  mf::LogInfo("FilterCRTPMTMatching::FilteringLevel ") << fFilterLevel;
  std::vector<CRTPMTMatching> EventFlashes;
  for(auto const & crtpmt : crtpmtMatches) {
    if(fSpillOnly==true){
	    if(crtpmt.flashInBeam==1) EventFlashes.push_back(crtpmt);
    }
    else {
	    if(crtpmt.flashInGate==1) EventFlashes.push_back(crtpmt);
    }
  } 

  bool hasOnlyCosmics = false;
  if (fFilterLevel == "loose") {
    hasOnlyCosmics = false;  // By default, with loose Filtering, everything is
                             // "intersting", nothing tagged as clear cosmic
    if(fSpillOnly==true && EventFlashes.empty()) hasOnlyCosmics = true;
  }
  else if (fFilterLevel == "medium") {
    hasOnlyCosmics = true;
    for (const auto& f : EventFlashes) {
      bool isCosmic = false;
      if (f.flashClassification == sbn::crt::MatchType::enTop || f.flashClassification == sbn::crt::MatchType::enTop_exSide ||
          f.flashClassification == sbn::crt::MatchType::enTop_mult ||
          f.flashClassification == sbn::crt::MatchType::enTop_exSide_mult) {
        isCosmic = true;
	}
      	// With Medium filter, everything (inTime) which is associated with Top
        // CRT Hit before the Flash is filtered as clear cosmic
      else
     	isCosmic = false;
      hasOnlyCosmics = hasOnlyCosmics && isCosmic;
    }
  }
  else if (fFilterLevel == "tight") {
    hasOnlyCosmics = true;
    for (const auto& f : EventFlashes) {
      bool isCosmic = false;
      if ((f.flashClassification != sbn::crt::MatchType::noMatch) && (f.flashClassification != sbn::crt::MatchType::exTop) &&
          (f.flashClassification != sbn::crt::MatchType::exSide))
        isCosmic = true;
        // With Medium filter, everything (inTime) which is associated with Top
        // CRT Hit before the Flash is filtered as clear cosmic
      else
         isCosmic = false;
      hasOnlyCosmics = hasOnlyCosmics && isCosmic;
    }
  }
  fFiltered = hasOnlyCosmics;
  if (fMatchTree) {
    for(const auto& f : EventFlashes){
      fClassification = (int) f.flashClassification;
      fOpFlashPos_X = f.flashPosition.X();
      fOpFlashPos_Y = f.flashPosition.Y();
      fOpFlashPos_Z = f.flashPosition.Z();
      fOpFlashTime_us = f.flashTime;
      finGate = f.flashInGate;
      finBeam = f.flashInBeam;
      fFlashBeamTime_ns = f.flashGateTime*1e3;
      for (const auto& crt : f.matchedCRTHits){
	fCRTHitPos_X.push_back(crt.position.X());
	fCRTHitPos_Y.push_back(crt.position.Y());
	fCRTHitPos_Z.push_back(crt.position.Z());
        fCRTHitTime_us.push_back(crt.time);
        fCRTFlashTime_ns.push_back(crt.PMTTimeDiff*1e3);
        fCRTHitRegion.push_back(crt.region);
        fCRTHitSystem.push_back(crt.sys);
      }
      fMatchTree->Fill();
      ClearVecs();
    }
    if (EventFlashes.empty()) {
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
