////////////////////////////////////////////////////////////////////////
// Class:       CRTPMTMatchingAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTPMTMatchingAna_module.cc
//
// Generated at Thu Feb 17 13:31:00 2022 by Biswaranjan Behera using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
//#include "larsim/MCCheater/PhotonBackTrackerService.h"

// Data product includes
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
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

namespace icarus {
namespace crt {
class CRTPMTMatchingAna;
struct CRTPMT {
  double tof;
  double distance;
  art::Ptr<sbn::crt::CRTHit> CRTHit;
};
enum hasCRTHit {
  noMatch = 0,       // No CRT match
  enTop = 1,         // matched with Top CRT hit before optical flash
  enSide = 2,        // matched with Side CRT hit before optical flash
  enTop_exSide = 3,  // entering from Top and exiting from side
  exTop = 4,         // matched with a Top CRT after the optical flash
  exSide = 5,        // matched with a Side CRT after the optical flash
  enTop_mult =
      6,  // matched with multiple Top CRT hits before the optical flash
  enTop_exSide_mult =
      7,  //  matched with multiple Top CRT hits before the optical flash and
          //  more then 1 side CRT hits after the optical flash
  others = 9  // all the other cases
};
struct MatchedCRT {
  // vector of pairs where first is the matched Time of Flight and second is the
  // matched CRTHit
  std::vector<CRTPMT> entering;
  std::vector<CRTPMT> exiting;
  hasCRTHit matchType;
};
struct TrajPoint {
  double X;
  double Y;
  double Z;
  double T;
};
}  // namespace crt
}  // namespace icarus

using namespace icarus::crt;

bool flashInTime(double const& flashTime, int gateType, double gateDiff,
                 double gateWidth) {
  // std::map<unsigned int, double> gateLength{
  //  {1, 2.7}, //BNB
  //  {2, 10.6}, //NuMI
  //  {3, 2.7}, //BNB Offbeam
  //  {4, 10.6}, //NuMI Offbeam
  //};

  // double vetoOffset = 3.5;
  double activeGate = gateWidth;

  double relFlashTime = flashTime + gateDiff / 1000.;
  mf::LogInfo("CRTPMTMatching::flashInTime")
      << "Gate Diff " << gateDiff / 1000 << " Ftime+gateDiff "
      << flashTime + gateDiff / 1000. << " " << activeGate;
  // std::cout<<"Gate Diff "<<gateDiff/1000<<"  "<<flashTime<<" Ftime+gateDiff
  // "<<flashTime + gateDiff/1000.<<" "<<activeGate<<std::endl;
  return ((relFlashTime > 0) && (relFlashTime < activeGate));
}

icarus::crt::MatchedCRT CRTHitmatched(
    const double& flashTime, geo::Point_t const& flashpos,
    std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, const double& interval) {
  std::vector<icarus::crt::CRTPMT> enteringCRTHits;
  std::vector<icarus::crt::CRTPMT> exitingCRTHits;

  hasCRTHit MatchType;
  int topen = 0, topex = 0, sideen = 0, sideex = 0;
  for (auto const crtHit : crtHits) {
    double tof = crtHit->ts1_ns / 1e3 - flashTime;
    // double distance = std::hypot(flashpos[0]-crtHit->x_pos,
    // flashpos[1]-crtHit->y_pos, flashpos[2]-crtHit->z_pos);
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

class icarus::crt::CRTPMTMatchingAna : public art::EDAnalyzer {
 public:
  using CRTHit = sbn::crt::CRTHit;

  explicit CRTPMTMatchingAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTPMTMatchingAna(CRTPMTMatchingAna const&) = delete;
  CRTPMTMatchingAna(CRTPMTMatchingAna&&) = delete;
  CRTPMTMatchingAna& operator=(CRTPMTMatchingAna const&) = delete;
  CRTPMTMatchingAna& operator=(CRTPMTMatchingAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

 private:
  // Declare member data here.

  bool HitCompare(const art::Ptr<CRTHit>& h1, const art::Ptr<CRTHit>& h2);
  void ClearVecs();

  // art::InputTag fOpHitModuleLabel;
  art::InputTag fOpFlashModuleLabel0;
  art::InputTag fOpFlashModuleLabel1;
  art::InputTag fOpFlashModuleLabel2;
  art::InputTag fOpFlashModuleLabel3;
  art::InputTag fCrtHitModuleLabel;
  art::InputTag fTriggerLabel;
  art::InputTag fTriggerConfiguration;
  // tart::InputTag fCrtTrackModuleLabel;

  int fEvent;   ///< number of the event being processed
  int fRun;     ///< number of the run being processed
  int fSubRun;  ///< number of the sub-run being processed

  // add trigger data product vars
  unsigned int m_gate_type;
  std::string m_gate_name;
  uint64_t m_trigger_timestamp;
  uint64_t m_gate_start_timestamp;
  uint64_t m_trigger_gate_diff;
  // uint64_t m_gate_crt_diff;
  uint64_t m_gate_width;

  double fCoinWindow;
  double fOpDelay;
  double fCrtDelay;
  int fPMTADCThresh;
  int fFlashPeThresh;
  int fHitPeThresh;
  double fFlashVelocity;
  double fFlashZOffset;
  double fHitVelocityMax;
  double fHitVelocityMin;
  int fTopBefore;
  int fTopAfter;
  int fSideBefore;
  int fSideAfter;
  double fBNBBeamGateMin;
  double fBNBBeamGateMax;
  double fBNBinBeamMin;
  double fBNBinBeamMax;

  double fNuMIBeamGateMin;
  double fNuMIBeamGateMax;
  double fNuMIinBeamMin;
  double fNuMIinBeamMax;
  //  CRTBackTracker* bt;
  CRTCommonUtils* crtutil;

  std::vector<art::InputTag> fFlashLabels;

  TTree* fMatchTree;

  // matchTree vars

  vector<double> fOpFlashPos;  // Position of the optical Flash Barycenter
  double fOpFlashTime_us;  // Time of the optical Flash w.r.t. Global Trigger
  double fOpFlashTimehit_us;
  double fOpFlashTimeAbs;
  bool fInTime;  // Was the OpFlash in beam spill time?
  bool fInTime_gate;
  bool fInTime_beam;
  double fOpFlashPE;       // Total PEs of the optical Flash
  vector<double> fOpHitX;  // X position of the OP hit
  vector<double> fOpHitY;  // Y position of the OP hit
  vector<double> fOpHitZ;  // Z position of the OP hit
  vector<double> fOpHitT;  // T of the OP hit
  vector<double> fOpHitA;  // Amplitude (ADC) of the OP hit
  int fNCRTmatch;
  vector<vector<double>> fMatchedCRTpos;  // Position of the matched CRTs
  // vector<double[3]>       fMatchedCRTpos_err; // Position error of the
  // matched CRTs
  vector<double>
      fMatchedCRTtime_us;  // Time of the matched CRT Hits w.r.t. Global Trigger
  vector<ULong64_t>
      fMatchedCRTtime_abs;        // Time of the matched CRT Hits in Unix time
  vector<int> fMatchedCRTregion;  // Region of the matched CRT Hits
  vector<vector<int>> fMatchedCRTmodID;
  vector<int> fMatchedCRTsys;           // Subsystem of the matched CRT Hit
  vector<double> fMatchedCRTamplitude;  // Amplitude of the matched CRT Hit
  vector<int> fDirection;  // Was the matched CRT before or after the flash?
                           // entering/exiting
  vector<double>
      fDistance;  // Distance between matched CRT and light barycenter
  vector<double>
      fTofOpHit;  // Time of Flight between matched CRT and first Optical Hit
  vector<double>
      fTofOpFlash;  // Time of Flight between matched CRT and Optical Flash
  vector<double> fVelocity;  // Assuming the correct match, evaluate the speed
                             // of the particle
  vector<double>
      fCRTGateDiff;  // Difference between CRTHit and BeamGate opening
  int fEventType;    // Was triggered the event?
  double fRelGateTime;
  /*int                    fNCrt;         //number of CRT hits per event
  vector<vector<double>> fCrtXYZT;      //Matched CRT hit x,y,z,t [cm/ns]
  vector<vector<double>> fCrtXYZErr;    //Matched CRT hit x,y,z,t error [cm/ns]
  vector<int>            fCrtRegion;    //Matched CRT hit region code
  vector<double>         fCrtPE;        //total PE's for Matched CRT hit
  vector<double>         fTofHit;       //CRT - PMT TOF using OpHit
  double         fTofFlash;     //CRT - PMT TOF using OpFlash
  double         fTofFlashHit;
  vector<double>         fDistHit;      //CRT - PMT distance [cm]
  double         fDistFlash;    //CRT - flash barycenter distance [cm]
  double         fDistFlashHit;
  vector<double>         fTofPeHit;     //total PE for matched OpHit
  double         fTofPeFlash;   //total PE for matched OpFlash
  vector<double>         fTofPeFlashHit;

  vector<vector<double>> fTofXYZTHit;   //position/time [cm/ns] for matched
  OpHit vector<double> fTofXYZTFlash; //position/time [cm/ns] for matched
  OpFlash vector<vector<double>> fTofXYZTFlashHit; vector<int> fTofTpcHit; //TPC
  for matched OpHit int            fTofTpcFlash;  //TPC for matched OpFlash
  vector<bool>           fMatchHit;     //was CRT hit matched to OpHit?
  bool           fMatchFlash;   //was CRT hit matched to OpFlash?
  bool           fInTime;   //was the OpFlash in beam spill time?
  vector<bool>           fTrackFilt;    //excluded by track filter */

  geo::GeometryCore const* fGeometryService;  ///< pointer to Geometry provider

  TTree* mSelectionTree;

  int mEvent;         ///< number of the event being processed
  int mRun;           ///< number of the run being processed
  int mSubRun;        ///< number of the sub-run being processed
  bool mEventFilter;  // Event should be filtered?
  int mEventType;     // EventType
  bool mInBeam;
};

icarus::crt::CRTPMTMatchingAna::CRTPMTMatchingAna(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // ,
      //,fOpHitModuleLabel(p.get<art::InputTag>("OpHitModuleLabel","ophit"))
      ,
      fOpFlashModuleLabel0(
          p.get<art::InputTag>("OpFlashModuleLabel0", "opflashTPC0")),
      fOpFlashModuleLabel1(
          p.get<art::InputTag>("OpFlashModuleLabel1", "opflashTPC1"))
      //  ,fOpFlashModuleLabel2(p.get<art::InputTag>("OpFlashModuleLabel2","opflashTPC2"))
      //,fOpFlashModuleLabel3(p.get<art::InputTag>("OpFlashModuleLabel3","opflashTPC3"))
      ,
      fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit")),
      fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")),
      fTriggerConfiguration(
          p.get<art::InputTag>("TriggerConfiguration", "triggerconfig")),
      fCoinWindow(p.get<double>("CoincidenceWindow", 60.0)),
      fOpDelay(p.get<double>("OpDelay", 55.1)),
      fCrtDelay(p.get<double>("CrtDelay", 1.6e6)),
      fPMTADCThresh(p.get<int>("PMTADCThresh", 400)),
      fFlashPeThresh(p.get<int>("FlashPeThresh", 9000)),
      fHitPeThresh(p.get<int>("HitPeThresh", 700)),
      fFlashVelocity(p.get<double>("FlashVelocityThresh", -40.)),
      fFlashZOffset(p.get<double>("FlashZOffset", 0.)),
      fHitVelocityMax(p.get<double>("HitVelocityMax", 20.)),
      fHitVelocityMin(p.get<double>("HitVelocityMin", 1.)),
      fBNBBeamGateMin(p.get<double>("BNBBeamGateMin", -550)),
      fBNBBeamGateMax(p.get<double>("BNBBeamGateMax", 2300)),
      fBNBinBeamMin(p.get<double>("BNBinBeamMin", -300)),
      fBNBinBeamMax(p.get<double>("BNBinBeamMax", 1300)),
      fNuMIBeamGateMin(p.get<double>("NuMIBeamGateMin", -550)),
      fNuMIBeamGateMax(p.get<double>("NuMIBeamGateMax", 10000)),
      fNuMIinBeamMin(p.get<double>("NuMIinBeamMin", -300)),
      fNuMIinBeamMax(p.get<double>("NuMIinBeamMax", 9100))
      // bt(new CRTBackTracker(p.get<fhicl::ParameterSet>("CRTBackTrack"))),
      ,
      crtutil(new CRTCommonUtils())
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  fFlashLabels.push_back(fOpFlashModuleLabel0);
  fFlashLabels.push_back(fOpFlashModuleLabel1);
  // fFlashLabels[2] = fOpFlashModuleLabel2;
  // fFlashLabels[3] = fOpFlashModuleLabel3;

  // Get a pointer to the geometry service provider.
  fGeometryService = lar::providerFrom<geo::Geometry>();

  art::ServiceHandle<art::TFileService> tfs;

  fMatchTree =
      tfs->make<TTree>("matchTree", "CRTHit - OpHit/Flash matching analysis");

  fMatchTree->Branch("event", &fEvent, "event/I");
  fMatchTree->Branch("run", &fRun, "run/I");
  fMatchTree->Branch("subrun", &fSubRun, "subrun/I");
  fMatchTree->Branch("opFlash_pos", &fOpFlashPos);
  fMatchTree->Branch("opFlash_time_us", &fOpFlashTime_us);
  fMatchTree->Branch("opFlashAbsTime", &fOpFlashTimeAbs);
  fMatchTree->Branch("opFlash_time_firsthit_us", &fOpFlashTimehit_us);
  fMatchTree->Branch("opHit_x", &fOpHitX);
  fMatchTree->Branch("opHit_y", &fOpHitY);
  fMatchTree->Branch("opHit_z", &fOpHitZ);
  fMatchTree->Branch("opHit_t", &fOpHitT);
  fMatchTree->Branch("opHit_amplitude", &fOpHitA);
  fMatchTree->Branch("inTime", &fInTime);
  fMatchTree->Branch("inTime_beam", &fInTime_beam);
  fMatchTree->Branch("inTime_gate", &fInTime_gate);
  fMatchTree->Branch("CRT_matches", &fNCRTmatch);
  fMatchTree->Branch("CRT_pos", &fMatchedCRTpos);
  fMatchTree->Branch("CRT_Time_us", &fMatchedCRTtime_us);
  fMatchTree->Branch("CRT_absTime", &fMatchedCRTtime_abs);
  fMatchTree->Branch("CRT_region", &fMatchedCRTregion);
  fMatchTree->Branch("CRT_FEB", &fMatchedCRTmodID);
  fMatchTree->Branch("CRT_subsys", &fMatchedCRTsys);
  fMatchTree->Branch("CRT_pes", &fMatchedCRTamplitude);
  fMatchTree->Branch("direction", &fDirection);
  fMatchTree->Branch("distance", &fDistance);
  fMatchTree->Branch("TOF_opflash", &fTofOpFlash);
  fMatchTree->Branch("TOF_ophit", &fTofOpHit);
  fMatchTree->Branch("TopHitsBefore", &fTopBefore);
  fMatchTree->Branch("TopHitsAfter", &fTopAfter);
  fMatchTree->Branch("SideHitsBefore", &fSideBefore);
  fMatchTree->Branch("SideHitsAfter", &fSideAfter);
  fMatchTree->Branch("CRT_gate_diff", &fCRTGateDiff);
  fMatchTree->Branch("eventType", &fEventType);
  fMatchTree->Branch("RelativeGateTime", &fRelGateTime);
  fMatchTree->Branch("gate_type", &m_gate_type, "gate_type/b");
  fMatchTree->Branch("gate_name", &m_gate_name);
  fMatchTree->Branch("trigger_timestamp", &m_trigger_timestamp,
                     "trigger_timestamp/l");
  fMatchTree->Branch("gate_start_timestamp", &m_gate_start_timestamp,
                     "gate_start_timestamp/l");
  fMatchTree->Branch("trigger_gate_diff", &m_trigger_gate_diff,
                     "trigger_gate_diff/l");
  // fMatchTree->Branch("gate_crt_diff",&m_gate_crt_diff, "gate_crt_diff/l");

  mSelectionTree = tfs->make<TTree>(
      "FilterTree", "Event Filter based on CRT Hit - PMT matches");
  mSelectionTree->Branch("event", &fEvent, "event/I");
  mSelectionTree->Branch("run", &fRun, "run/I");
  mSelectionTree->Branch("subrun", &fSubRun, "subrun/I");
  mSelectionTree->Branch("EventFilter", &mEventFilter);
  mSelectionTree->Branch("EventType", &mEventType);
  mSelectionTree->Branch("InBeam", &mInBeam);
}

void icarus::crt::CRTPMTMatchingAna::analyze(art::Event const& e) {
  // Implementation of required member function here.
  /*
  geo::CryostatGeo const& cryo0 = fGeometryService->Cryostat(0);
  geo::CryostatGeo const& cryo1 = fGeometryService->Cryostat(1);
  geo::TPCGeo const& tpc00 = cryo0.TPC(0);
  geo::TPCGeo const& tpc01 = cryo0.TPC(1);
  geo::TPCGeo const& tpc10 = cryo1.TPC(0);
  geo::TPCGeo const& tpc11 = cryo1.TPC(1);
  */
  MF_LOG_DEBUG("CRTPMTMatching: ") << "beginning analyis" << '\n';
  std::cout << "START" << std::endl;
  // Start by fetching some basic event information for our n-tuple.
  fEvent = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();

  mEvent = e.id().event();
  mRun = e.run();
  mSubRun = e.subRun();

  ClearVecs();

  // add trigger info
  if (!fTriggerLabel.empty()) {
    art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
    e.getByLabel(fTriggerLabel, trigger_handle);
    if (trigger_handle.isValid()) {
      sbn::triggerSource bit = trigger_handle->sourceType;
      m_gate_type = (unsigned int)bit;
      m_gate_name = bitName(bit);
      m_trigger_timestamp = trigger_handle->triggerTimestamp;
      m_gate_start_timestamp = trigger_handle->beamGateTimestamp;
      m_trigger_gate_diff =
          trigger_handle->triggerTimestamp - trigger_handle->beamGateTimestamp;
      // Read Beam Gate Size
      art::Run const& run = e.getRun();
      auto const& trigger_configuration =
          run.getProduct<icarus::TriggerConfiguration>(fTriggerConfiguration);
      m_gate_width = trigger_configuration.getGateWidth(m_gate_type);
    } else {
      mf::LogError("CRTPMTMatching:")
          << "No sbn::ExtraTriggerInfo associated to label: "
          << fTriggerLabel.encode() << "\n";
    }
  } else {
    // std::cout  << "Trigger Data product " << fTriggerLabel.label() << " not
    // found!\n" ;
  }

  // OpHits
  // art::Handle< std::vector<recob::OpHit> > opHitListHandle;
  // std::vector< art::Ptr<recob::OpHit> >    opHitList;
  // if( e.getByLabel(fOpHitModuleLabel,opHitListHandle) )
  //  art::fill_ptr_vector(opHitList, opHitListHandle);

  // OpFlash
  std::array<art::Handle<std::vector<recob::OpFlash>>, 2U> flashHandles;
  for (int i = 0; i < 2; i++) {
    flashHandles[i] = e.getHandle<std::vector<recob::OpFlash>>(fFlashLabels[i]);
  }
  /*
  map<int, art::Handle< std::vector<recob::OpFlash> > > flashHandles;
  std::map<int,std::vector< art::Ptr<recob::OpFlash> >> opFlashLists;
  //  fNFlash = 0;
  for(int i=0; i<2; i++) {
    if( e.getByLabel(fFlashLabels[i],flashHandles[i]) )
      art::fill_ptr_vector(opFlashLists[i], flashHandles[i]);
  }*/

  // CRTHits
  art::Handle<std::vector<CRTHit>> crtHitListHandle;
  std::vector<art::Ptr<CRTHit>> crtHitList;
  if (e.getByLabel(fCrtHitModuleLabel, crtHitListHandle))
    art::fill_ptr_vector(crtHitList, crtHitListHandle);

  // fNCrt = crtHitList.size();

  bool Filter = false;
  hasCRTHit Type = others;
  // for(auto const& flashList : opFlashLists){
  for (art::InputTag const& flashLabel : fFlashLabels) {
    auto const flashHandle =
        e.getHandle<std::vector<recob::OpFlash>>(flashLabel);
    art::FindMany<recob::OpHit> findManyHits(flashHandle, e, flashLabel);

    for (auto const& [iflash, flash] : util::enumerate(*flashHandle)) {
      // auto const& flash = flashList.second[iflash];
      hasCRTHit eventType = others;
      double tflash = flash.Time();
      double tAbsflash = flash.AbsTime();
      vector<recob::OpHit const*> const& hits = findManyHits.at(iflash);
      int nPMTsTriggering = 0;
      double firstTime = 999999;
      geo::vect::MiddlePointAccumulator flashCentroid;
      // double flash_pos[3]={0,0,0};
      double ampsum = 0, t_m = 0;
      for (auto const& hit : hits) {
        if (hit->Amplitude() > fPMTADCThresh) nPMTsTriggering++;
        if (firstTime > hit->PeakTime()) firstTime = hit->PeakTime();
        geo::Point_t const pos =
            fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel())
                .GetCenter();
        double amp = hit->Amplitude();
        ampsum += amp;
        fOpHitX.push_back(pos.X());
        fOpHitY.push_back(pos.Y());
        fOpHitZ.push_back(pos.Z());
        fOpHitT.push_back(hit->PeakTime());
        fOpHitA.push_back(amp);
        /*flash_pos[0]=flash_pos[0]+pos.X()*amp;
        flash_pos[1]=flash_pos[1]+pos.Y()*amp;
        flash_pos[2]=flash_pos[2]+pos.Z()*amp;*/
        flashCentroid.add(pos, amp);
        t_m = t_m + hit->PeakTime();
      }
      /*flash_pos[0]=flash_pos[0]/ampsum;
      flash_pos[1]=flash_pos[1]/ampsum;
      flash_pos[2]=flash_pos[2]/ampsum;*/
      // auto flash_pos = flashCentroid.middlePointAs<double[3]>();
      geo::Point_t flash_pos = flashCentroid.middlePoint();
      t_m = t_m / nPMTsTriggering;
      if (nPMTsTriggering < 5) {
        continue;
      }
      double gateDiff =
          (Long64_t)m_trigger_timestamp - (Long64_t)m_gate_start_timestamp;
      bool inTime = flashInTime(firstTime, m_gate_type, gateDiff, m_gate_width);
      fRelGateTime = gateDiff + (tAbsflash - 1500) * 1e3;
      fInTime_gate = false;
      fInTime_beam = false;
      if (m_gate_type == 1 || m_gate_type == 3) {  // BNB OffBeamBNB
        if (fRelGateTime > fBNBBeamGateMin && fRelGateTime < fBNBBeamGateMax)
          fInTime_gate = true;
        if (fRelGateTime > fBNBinBeamMin && fRelGateTime < fBNBinBeamMax)
          fInTime_beam = true;
      }
      if (m_gate_type == 2 || m_gate_type == 4) {  // NuMI OffBeamNuMI
        if (fRelGateTime > fNuMIBeamGateMin && fRelGateTime < fNuMIBeamGateMax)
          fInTime_gate = true;
        if (fRelGateTime > fNuMIinBeamMin && fRelGateTime < fNuMIinBeamMax)
          fInTime_beam = true;
      }
      inTime = fInTime_gate;
      // Now get the CRT. Search the CRT Hits within -100 from the flash time
      //  for the future a differentiation between Top and Side and some
      //  considerations based on the proximity of the light are necessary
      icarus::crt::MatchedCRT CRTmatches =
          CRTHitmatched(firstTime, flash_pos, crtHitList, 0.1);
      int TopEn = 0, TopEx = 0, SideEn = 0, SideEx = 0;

      auto nCRTHits = CRTmatches.entering.size() + CRTmatches.exiting.size();
      double peflash = flash.TotalPE();
      if (nCRTHits >= 0) {
        eventType = CRTmatches.matchType;
        for (auto& entering : CRTmatches.entering) {
          vector<double> CRTpos = {entering.CRTHit->x_pos,
                                   entering.CRTHit->y_pos,
                                   entering.CRTHit->z_pos};
          double CRTtime = entering.CRTHit->ts1_ns / 1e3;
          ULong64_t CRTAbsTime = entering.CRTHit->ts0_s;
          int CRTRegion = entering.CRTHit->plane;
          int CRTSys = 0;
          if (CRTRegion >= 36)
            CRTSys =
                1;  // Very lazy way to determine if the Hit is a Top or a Side.
                    // Will update it when bottom CRT will be availble
          double CRTPe = entering.CRTHit->peshit;
          int CRTDirection = 0;  // 0=entering, 1=exiting
          double CRTDistance = entering.distance;
          double CRTTof_ophit = entering.tof;
          double CRTTof_opflash = CRTtime - tflash;
          double CRTGateDiff = CRTAbsTime - m_gate_start_timestamp;
          std::vector<int> HitFebs;
          for (auto& crts : entering.CRTHit->feb_id) {
            HitFebs.emplace_back((int)crts);
          }
          if (CRTSys == 0) TopEn++;
          if (CRTSys == 1) SideEn++;
          fMatchedCRTmodID.emplace_back(HitFebs);
          fMatchedCRTpos.emplace_back(CRTpos);
          fMatchedCRTtime_us.emplace_back(CRTtime);
          fTofOpFlash.emplace_back(CRTTof_opflash);
          fMatchedCRTtime_abs.emplace_back(CRTAbsTime);
          fMatchedCRTregion.emplace_back(CRTRegion);
          fMatchedCRTsys.emplace_back(CRTSys);
          fMatchedCRTamplitude.emplace_back(CRTPe);
          fDirection.emplace_back(CRTDirection);
          fDistance.emplace_back(CRTDistance);
          fTofOpHit.emplace_back(CRTTof_ophit);
          fCRTGateDiff.emplace_back(CRTGateDiff);
          HitFebs.clear();
          // Fill
        }
        for (auto& exiting : CRTmatches.exiting) {
          vector<double> CRTpos = {exiting.CRTHit->x_pos, exiting.CRTHit->y_pos,
                                   exiting.CRTHit->z_pos};
          double CRTtime = exiting.CRTHit->ts1_ns / 1e3;
          ULong64_t CRTAbsTime = exiting.CRTHit->ts0_s;
          int CRTRegion = exiting.CRTHit->plane;
          int CRTSys = 0;
          if (CRTRegion >= 36) CRTSys = 1;
          double CRTPe = exiting.CRTHit->peshit;
          int CRTDirection = 1;  // 0=entering, 1=exiting
          double CRTDistance = exiting.distance;
          double CRTTof_ophit = exiting.tof;
          double CRTTof_opflash = CRTtime - tflash;
          double CRTGateDiff = CRTAbsTime - m_gate_start_timestamp;
          std::vector<int> HitFebs;
          for (auto& crts : exiting.CRTHit->feb_id) {
            HitFebs.emplace_back((int)crts);
          }
          if (CRTSys == 0) TopEx++;
          if (CRTSys == 1) SideEx++;
          fMatchedCRTmodID.emplace_back(HitFebs);
          fMatchedCRTpos.push_back(std::move(CRTpos));
          fMatchedCRTtime_us.emplace_back(CRTtime);
          fMatchedCRTtime_abs.emplace_back(CRTAbsTime);
          fMatchedCRTregion.emplace_back(CRTRegion);
          fMatchedCRTsys.emplace_back(CRTSys);
          fMatchedCRTamplitude.emplace_back(CRTPe);
          fDirection.emplace_back(CRTDirection);
          fDistance.emplace_back(CRTDistance);
          fTofOpFlash.emplace_back(CRTTof_opflash);
          fTofOpHit.emplace_back(CRTTof_ophit);
          fCRTGateDiff.emplace_back(CRTGateDiff);
          // Fill
        }
        if ((eventType == enTop || eventType == enTop_exSide ||
             eventType == enTop_mult || eventType == enTop_exSide_mult) &&
            inTime == true)
          Filter = true;
        if (inTime == true) Type = eventType;
      }  // if matched
      fTopBefore = TopEn;
      fTopAfter = TopEx;
      fSideBefore = SideEn;
      fSideAfter = SideEx;
      fOpFlashPE = peflash;
      fOpFlashPos = {(double)flash_pos.X(), (double)flash_pos.Y(),
                     (double)flash_pos.Z()};
      fOpFlashTime_us = tflash;
      fOpFlashTimeAbs = tAbsflash;
      // fRelGateTime=gateDiff+(tAbsflash-1500)*1e3;
      //  std::cout<<"Gate diff "<<gateDiff<<" relTime
      //  "<<fRelGateTime<<std::endl;
      fOpFlashTimehit_us = firstTime;
      fInTime = inTime;
      fEventType = eventType;
      fNCRTmatch = nCRTHits;
      fMatchTree->Fill();
      fOpHitX.clear();
      fOpHitY.clear();
      fOpHitZ.clear();
      fOpHitT.clear();
      fOpHitA.clear();
      fOpFlashPos.clear();
      fMatchedCRTmodID.clear();
      fMatchedCRTpos.clear();
      fMatchedCRTtime_us.clear();
      fMatchedCRTtime_abs.clear();
      fMatchedCRTregion.clear();
      fMatchedCRTsys.clear();
      fMatchedCRTamplitude.clear();
      fDirection.clear();
      fDistance.clear();
      fTofOpHit.clear();
      fTofOpFlash.clear();
      fVelocity.clear();
      fCRTGateDiff.clear();
    }  // for Flash
  }
  mEventType = Type;
  mEventFilter = Filter;
  mInBeam = fInTime_beam;
  mSelectionTree->Fill();
}

void icarus::crt::CRTPMTMatchingAna::ClearVecs() {
  // matchTree
  fOpHitX.clear();
  fOpHitY.clear();
  fOpHitZ.clear();
  fOpHitT.clear();
  fOpHitA.clear();
  fOpFlashPos.clear();
  fMatchedCRTmodID.clear();
  fMatchedCRTpos.clear();
  fMatchedCRTtime_us.clear();
  fMatchedCRTtime_abs.clear();
  fMatchedCRTregion.clear();
  fMatchedCRTsys.clear();
  fMatchedCRTamplitude.clear();
  fDirection.clear();
  fDistance.clear();
  fDistance.clear();
  fTofOpHit.clear();
  fVelocity.clear();
  fCRTGateDiff.clear();
  /*fCrtXYZT.clear();
  fCrtXYZErr.clear();
  fCrtRegion.clear();
  fCrtPE.clear();
  fTofHit.clear();
  fTofFlash.clear();
  fTofFlashHit.clear();
  fDistHit.clear();
  fDistFlash.clear();
  fDistFlashHit.clear();
  fTofPeHit.clear();
  fTofPeFlash.clear();
  fTofPeFlashHit.clear();
  fTofXYZTHit.clear();
  fTofXYZTFlash.clear();
  fTofXYZTFlashHit.clear();
  fTofTpcHit.clear();
  fTofTpcFlash.clear();
  fMatchHit.clear();
  fMatchFlash.clear();*/
}

DEFINE_ART_MODULE(icarus::crt::CRTPMTMatchingAna)
