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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
//#include "larsim/MCCheater/PhotonBackTrackerService.h"

//Data product includes
#include "nusimdata/SimulationBase/MCGeneratorInfo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

//C++ includes
#include <vector>
#include <map>

//ROOT includes
#include "TTree.h"
#include "TVector3.h"

using std::vector;
using std::map;


namespace icarus {
  namespace crt {
    class CRTPMTMatchingAna;
  }
}

using namespace icarus::crt;

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

  art::InputTag fOpHitModuleLabel;
  art::InputTag fOpFlashModuleLabel0;
  art::InputTag fOpFlashModuleLabel1;
  art::InputTag fOpFlashModuleLabel2;
  art::InputTag fOpFlashModuleLabel3;
  art::InputTag fCrtHitModuleLabel;
  art::InputTag fTriggerLabel;
  //tart::InputTag fCrtTrackModuleLabel;

  int fEvent;        ///< number of the event being processed
  int fRun;          ///< number of the run being processed
  int fSubRun;       ///< number of the sub-run being processed

  //add trigger data product vars
  unsigned int m_gate_type;
  std::string m_gate_name;
  uint64_t m_trigger_timestamp;
  uint64_t m_gate_start_timestamp;
  uint64_t m_trigger_gate_diff;
  uint64_t m_gate_crt_diff;

  double        fCoinWindow;
  double        fOpDelay;
  double        fCrtDelay;
  int           fFlashPeThresh;
  int           fHitPeThresh;
  double        fFlashVelocity;
  double        fFlashZOffset;
  double        fHitVelocityMax;
  double        fHitVelocityMin;

  //  CRTBackTracker* bt;
  CRTCommonUtils *crtutil;

  map<int,art::InputTag> fFlashLabels;

  TTree* fMatchTree;

  //matchTree vars
  int                    fNCrt;         //number of CRT hits per event
  vector<vector<double>> fCrtXYZT;      //CRT hit x,y,z,t [cm/ns]
  vector<vector<double>> fCrtXYZErr;    //CRT hit x,y,z,t error [cm/ns]
  vector<int>            fCrtRegion;    //CRT hit region code
  vector<double>         fCrtPE;        //total PE's for CRT hit
  vector<double>         fTofHit;       //CRT - PMT TOF using OpHit
  vector<double>         fTofFlash;     //CRT - PMT TOF using OpFlash
  vector<double>         fTofFlashHit;
  vector<double>         fDistHit;      //CRT - PMT distance [cm]
  vector<double>         fDistFlash;    //CRT - flash barycenter distance [cm]
  vector<double>         fDistFlashHit;
  vector<double>         fTofPeHit;     //total PE for matched OpHit
  vector<double>         fTofPeFlash;   //total PE for matched OpFlash
  vector<double>         fTofPeFlashHit;

  vector<vector<double>> fTofXYZTHit;   //position/time [cm/ns] for matched OpHit
  vector<vector<double>> fTofXYZTFlash; //position/time [cm/ns] for matched OpFlash
  vector<vector<double>> fTofXYZTFlashHit;
  vector<int>            fTofTpcHit;    //TPC for matched OpHit
  vector<int>            fTofTpcFlash;  //TPC for matched OpFlash
  vector<bool>           fMatchHit;     //was CRT hit matched to OpHit?
  vector<bool>           fMatchFlash;   //was CRT hit matched to OpFlash?
  vector<bool>           fTrackFilt;    //excluded by track filter

  geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
};


icarus::crt::CRTPMTMatchingAna::CRTPMTMatchingAna(fhicl::ParameterSet const& p)
  :EDAnalyzer{p}  // ,
  ,fOpHitModuleLabel(p.get<art::InputTag>("OpHitModuleLabel","ophit"))
  ,fOpFlashModuleLabel0(p.get<art::InputTag>("OpFlashModuleLabel0","opflashTPC0"))
  ,fOpFlashModuleLabel1(p.get<art::InputTag>("OpFlashModuleLabel1","opflashTPC1"))
  //  ,fOpFlashModuleLabel2(p.get<art::InputTag>("OpFlashModuleLabel2","opflashTPC2"))
  //,fOpFlashModuleLabel3(p.get<art::InputTag>("OpFlashModuleLabel3","opflashTPC3"))
  ,fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel","crthit"))
  ,fTriggerLabel(p.get<art::InputTag>("TriggerLabel","daqTrigger"))
  //  ,fCrtTrackModuleLabel(p.get<art::InputTag>("CrtTrackModuleLabel","crttrack"))
  ,fCoinWindow(p.get<double>("CoincidenceWindow",60.0))
  ,fOpDelay(p.get<double>("OpDelay",55.1))
  ,fCrtDelay(p.get<double>("CrtDelay",1.6e6))
  ,fFlashPeThresh(p.get<int>("FlashPeThresh",9000))
  ,fHitPeThresh(p.get<int>("HitPeThresh",700))
  ,fFlashVelocity(p.get<double>("FlashVelocityThresh",-40.))
  ,fFlashZOffset(p.get<double>("FlashZOffset",0.))
  ,fHitVelocityMax(p.get<double>("HitVelocityMax",20.))
  ,fHitVelocityMin(p.get<double>("HitVelocityMin",1.))
   // bt(new CRTBackTracker(p.get<fhicl::ParameterSet>("CRTBackTrack"))),
  ,crtutil(new CRTCommonUtils())
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fFlashLabels[0] = fOpFlashModuleLabel0;
  fFlashLabels[1] = fOpFlashModuleLabel1;
  fFlashLabels[2] = fOpFlashModuleLabel2;
  fFlashLabels[3] = fOpFlashModuleLabel3;

  // Get a pointer to the geometry service provider.
  fGeometryService = lar::providerFrom<geo::Geometry>();

  art::ServiceHandle<art::TFileService> tfs;

  fMatchTree = tfs->make<TTree>("matchTree","CRTHit - OpHit/Flash matching analysis");

  fMatchTree->Branch("event",        &fEvent,          "event/I");
  fMatchTree->Branch("run",          &fRun,            "run/I");
  fMatchTree->Branch("subrun",       &fSubRun,         "subrun/I");
  fMatchTree->Branch("ncrt",         &fNCrt,      "nCrtHit/I");
  fMatchTree->Branch("crtXYZT",      &fCrtXYZT);
  fMatchTree->Branch("crtXYZErr",    &fCrtXYZErr);
  fMatchTree->Branch("crtPE",        &fCrtPE);
  fMatchTree->Branch("crtRegion",    &fCrtRegion);
  fMatchTree->Branch("tofHit",       &fTofHit);
  fMatchTree->Branch("tofFlash",     &fTofFlash);
  fMatchTree->Branch("tofFlashHit",  &fTofFlashHit);
  fMatchTree->Branch("peHit",        &fTofPeHit);
  fMatchTree->Branch("peFlash",      &fTofPeFlash);
  fMatchTree->Branch("peFlashHit",   &fTofPeFlashHit);
  fMatchTree->Branch("xyztHit",      &fTofXYZTHit);
  fMatchTree->Branch("xyztFlash",    &fTofXYZTFlash);
  fMatchTree->Branch("xyztFlashHit", &fTofXYZTFlashHit);
  fMatchTree->Branch("distHit",      &fDistHit);
  fMatchTree->Branch("distFlash",    &fDistFlash);
  fMatchTree->Branch("distFlashHit", &fDistFlashHit);
  fMatchTree->Branch("tpcHit",       &fTofTpcHit);
  fMatchTree->Branch("tpcFlash",     &fTofTpcFlash);
  fMatchTree->Branch("gate_type", &m_gate_type, "gate_type/b");
  fMatchTree->Branch("gate_name", &m_gate_name);
  fMatchTree->Branch("trigger_timestamp", &m_trigger_timestamp, "trigger_timestamp/l");
  fMatchTree->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");
  fMatchTree->Branch("trigger_gate_diff", &m_trigger_gate_diff, "trigger_gate_diff/l");
  fMatchTree->Branch("gate_crt_diff",&m_gate_crt_diff, "gate_crt_diff/l");
}

void icarus::crt::CRTPMTMatchingAna::analyze(art::Event const& e)
{
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

  // Start by fetching some basic event information for our n-tuple.
  fEvent  = e.id().event();
  fRun    = e.run();
  fSubRun = e.subRun();

  ClearVecs();

  //add trigger info
  if( !fTriggerLabel.empty() ) {

    art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
    e.getByLabel( fTriggerLabel, trigger_handle );
    if( trigger_handle.isValid() ) {
      sbn::triggerSource bit = trigger_handle->sourceType;
      m_gate_type = (unsigned int)bit;
      m_gate_name = bitName(bit);
      m_trigger_timestamp = trigger_handle->triggerTimestamp;
      m_gate_start_timestamp =  trigger_handle->beamGateTimestamp;
      m_trigger_gate_diff = trigger_handle->triggerTimestamp - trigger_handle->beamGateTimestamp;

    }
    else{
      mf::LogError("CRTPMTMatching:") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
    }
  }
  else {
    std::cout  << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
  }


  //OpHits
  art::Handle< std::vector<recob::OpHit> > opHitListHandle;
  std::vector< art::Ptr<recob::OpHit> >    opHitList;
  if( e.getByLabel(fOpHitModuleLabel,opHitListHandle) )
    art::fill_ptr_vector(opHitList, opHitListHandle);


  //OpFlash
  map<int, art::Handle< std::vector<recob::OpFlash> > > flashHandles;
  std::map<int,std::vector< art::Ptr<recob::OpFlash> >> opFlashLists;
  //  fNFlash = 0;
  for(int i=0; i<2; i++) {
    if( e.getByLabel(fFlashLabels[i],flashHandles[i]) )
      art::fill_ptr_vector(opFlashLists[i], flashHandles[i]);
  }


  //CRTHits
  art::Handle< std::vector<CRTHit> > crtHitListHandle;
  std::vector< art::Ptr<CRTHit> >    crtHitList;
  if( e.getByLabel(fCrtHitModuleLabel,crtHitListHandle))
    art::fill_ptr_vector(crtHitList, crtHitListHandle);

  fNCrt = crtHitList.size();

  for(auto const& crt : crtHitList){
    vector<double> xyzt, xyzerr;
    TVector3 rcrt(crt->x_pos,crt->y_pos,crt->z_pos);
    m_gate_crt_diff = m_gate_start_timestamp - crt->ts0_ns;
    xyzt.push_back(rcrt.X());
    xyzt.push_back(rcrt.Y());
    xyzt.push_back(rcrt.Z());

    double tcrt = double(m_gate_start_timestamp - crt->ts0_ns)/1e3;
    tcrt = -tcrt+1e6;

    //    double tcrt = (int32_t)crt->ts0_ns;// - fCrtDelay;
    //uint64_t tcrt = crt->ts0_ns;

    //    std::cout <<"T0: " << crt->ts0_ns << " , tcrt : " << tcrt << std::endl;
    xyzt.push_back(tcrt);
    fCrtXYZT.push_back(xyzt);

    xyzerr.push_back(crt->x_err);
    xyzerr.push_back(crt->y_err);
    xyzerr.push_back(crt->z_err);
    fCrtXYZErr.push_back(xyzerr);

    fCrtPE.push_back(crt->peshit);
    fCrtRegion.push_back(crtutil->AuxDetRegionNameToNum(crt->tagger));

    // -- flash match --
    int matchtpc = -1;
    double tdiff = DBL_MAX, rdiff=DBL_MAX, peflash=DBL_MAX;
    bool matched=false;
    xyzt.clear();
    double flashHitT = DBL_MAX, flashHitPE=DBL_MAX, flashHitDiff=DBL_MAX;
    vector<double> flashHitxyzt;

    for(auto const& flashList : opFlashLists) {

      art::FindManyP<recob::OpHit> findManyHits(flashHandles[flashList.first], e, fFlashLabels[flashList.first]);

      for(size_t iflash=0; iflash<flashList.second.size(); iflash++) {

	auto const& flash = flashList.second[iflash];
	if(flash->TotalPE()<fFlashPeThresh){
	  continue;
	}

	double tflash = flash->Time();//*1e3;//-fOpDelay;

	TVector3 rflash(0,flash->YCenter(),flash->ZCenter());
	TVector3 vdiff = rcrt-rflash;
	//std::cout << "flash time: "<<flash->Time() << " absdiff : "<< abs(tcrt-tflash)<< std::endl;
	if(abs(tcrt-tflash)<abs(tdiff)) {
	  peflash = flash->TotalPE();
	  tdiff = tcrt-tflash;
	  rdiff = vdiff.Mag();
	  xyzt.clear();
	  xyzt.push_back(rflash.X());
	  xyzt.push_back(rflash.Y());
	  xyzt.push_back(rflash.Z());
	  xyzt.push_back(tflash);
	  matched = true;
	  matchtpc = flashList.first;

	  vector<art::Ptr<recob::OpHit>> hits = findManyHits.at(iflash);
	  for(auto const& hit : hits) {
	    double tPmt = hit->PeakTime();//*1.e3;//-fOpDelay;
	    if( tPmt < flashHitT) {
	      flashHitT = tPmt;
	      flashHitPE = hit->PE();
	      //FlashHit position/time
	      double pos[3];
	      fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter(pos);
	      flashHitxyzt.clear();
	      for(int i=0; i<3; i++) flashHitxyzt.push_back(pos[i]);
	      flashHitxyzt.push_back(flashHitT);

	      //FlashHit distance
	      TVector3 rflashHit(pos[0],pos[1],pos[2]);
	      TVector3 vdiffHit = rcrt-rflashHit;
	      flashHitDiff = vdiffHit.Mag();
	    }
	  }//loop over flash hits
	}//if minimum tdiff
      }//for OpFlash in this flash list
    }//for flash lists
    if(!matched) {
      peflash = DBL_MAX;
      for(int i=0; i<4; i++) xyzt.push_back(DBL_MAX);
    }

    fMatchFlash.push_back(matched);
    fTofFlash.push_back(tdiff);
    fTofPeFlash.push_back(peflash);
    fTofXYZTFlash.push_back(xyzt);
    fDistFlash.push_back(rdiff);
    fTofTpcFlash.push_back(matchtpc);
    fTofFlashHit.push_back(tcrt-flashHitT);
    fTofPeFlashHit.push_back(flashHitPE);
    fTofXYZTFlashHit.push_back(flashHitxyzt);
    fDistFlashHit.push_back(flashHitDiff);

    // -- match OpHits to CRTHits --
    tdiff = DBL_MAX;
    matched = false;
    peflash = DBL_MAX;
    rdiff = DBL_MAX;
    xyzt.clear();

    double pemax = 0.;
    for(auto const& hit : opHitList) {
      double thit = hit->PeakTime();//*1e3-fOpDelay;
      if(hit->PE()<fHitPeThresh){
	continue;
      }


      //if(abs(tcrt-thit)<abs(tdiff)) {
      if(abs(tcrt-thit)<fCoinWindow && hit->PE()>pemax) {
	pemax = hit->PE();

	//hitXYZT
	double pos[3];
	fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter(pos);

	//distHit
	TVector3 rhit (pos[0],pos[1],pos[2]);
	TVector3 vdiff = rcrt-rhit;
	rdiff = vdiff.Mag();
	peflash = hit->PE();
	tdiff = tcrt-thit;
	xyzt.clear();
	for(int i=0; i<3; i++) xyzt.push_back(pos[i]);
	xyzt.push_back(thit);

	matched = true;
	std::cout << "thit: "<<thit << " , tcrt: " << tcrt << " , tdiff " << tdiff  << std::endl;
      }//if min tof
    }//for OpHits
    if(!matched) {
      tdiff = DBL_MAX;
      peflash = DBL_MAX;
      rdiff = DBL_MAX;
      xyzt.clear();
      for(int i=0; i<4; i++) xyzt.push_back(DBL_MAX);
    }
    fMatchHit.push_back(matched);
    fTofHit.push_back(tdiff);
    fTofPeHit.push_back(peflash);
    fTofXYZTHit.push_back(xyzt);
    fDistHit.push_back(rdiff);
    fTofTpcHit.push_back(matchtpc);
  }//for CRTHits

  fMatchTree->Fill();
}

void icarus::crt::CRTPMTMatchingAna::ClearVecs()
{
  //matchTree
  fCrtXYZT.clear();
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
  fMatchFlash.clear();
}

bool icarus::crt::CRTPMTMatchingAna::HitCompare(const art::Ptr<CRTHit>& hit1, const art::Ptr<CRTHit>& hit2) {

  if(hit1->ts1_ns != hit2->ts1_ns) return false;
  if(hit1->plane  != hit2->plane) return false;
  if(hit1->x_pos  != hit2->x_pos) return false;
  if(hit1->y_pos  != hit2->y_pos) return false;
  if(hit1->z_pos  != hit2->z_pos) return false;
  if(hit1->x_err  != hit2->x_err) return false;
  if(hit1->y_err  != hit2->y_err) return false;
  if(hit1->z_err  != hit2->z_err) return false;
  if(hit1->tagger != hit2->tagger) return false;

  return true;
}

DEFINE_ART_MODULE(icarus::crt::CRTPMTMatchingAna)
