////////////////////////////////////////////////////////////////////////
// Class:       CrtOpHitMatchAnalysis
// Plugin Type: analyzer (art v3_05_00)
// File:        CrtOpHitMatchAnalysis_module.cc
//
// Generated at Tue Apr 14 19:49:18 2020 by Christopher Hilgenberg using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

//Framework includes
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

//LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


//Data product includes
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"

//C++ includes
#include <vector>

//ROOT includes
#include "TTree.h"

using std::vector;

class CrtOpHitMatchAnalysis;

class CrtOpHitMatchAnalysis : public art::EDAnalyzer {
public:
  explicit CrtOpHitMatchAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CrtOpHitMatchAnalysis(CrtOpHitMatchAnalysis const&) = delete;
  CrtOpHitMatchAnalysis(CrtOpHitMatchAnalysis&&) = delete;
  CrtOpHitMatchAnalysis& operator=(CrtOpHitMatchAnalysis const&) = delete;
  CrtOpHitMatchAnalysis& operator=(CrtOpHitMatchAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  const detinfo::DetectorClocks* fClock;

  art::InputTag fOpHitModuleLabel;
  art::InputTag fOpFlashModuleLabel0;
  art::InputTag fOpFlashModuleLabel1;
  art::InputTag fCrtHitModuleLabel;
  double        fCoinWindow;
  double        fOpDelay;

  TTree* fTree;

  int            fNCrtHit;
  int            fNOpHit;
  int            fNOpFlash;
  vector<double> fTCrt;
  vector<double> fTOp;
  vector<double> fTFlash;
  vector<double> fPeFlash;
  vector<double> fTofHit;
  vector<double> fTofFlash;
  vector<double> fTofPe;
  double         fTTrig;

  
};


CrtOpHitMatchAnalysis::CrtOpHitMatchAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
    fOpHitModuleLabel(p.get<art::InputTag>("OpHitModuleLabel","ophit")),
    fOpFlashModuleLabel0(p.get<art::InputTag>("OpFlashModuleLabel0","opflashTPC0")),
    fOpFlashModuleLabel1(p.get<art::InputTag>("OpFlashModuleLabel1","opflashTPC1")),
    fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel","crtsimhit")),
    fCoinWindow(p.get<double>("CoincidenceWindow",60.0)),
    fOpDelay(p.get<double>("OpDelay",55.1)),
    fNCrtHit(0),
    fNOpHit(0),
    fNOpFlash(0),
    fTCrt({}),
    fTOp({}),
    fTFlash({}),
    fPeFlash({}),
    fTofHit({}),
    fTofFlash({}),
    fTofPe({}),
    fTTrig(DBL_MAX) {
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("CrtOpHitTree","CRTHit - OpHit matching analysis");
  fTree->Branch("nOpHit",   &fNOpHit,   "nOpHit/I");
  fTree->Branch("nOpFlash", &fNOpFlash, "nOpFlash/I");
  fTree->Branch("nCrtHit",  &fNCrtHit,  "nCrtHit/I");
  fTree->Branch("tCrt",     &fTCrt);
  fTree->Branch("tOp",      &fTOp);
  fTree->Branch("tFlash",   &fTFlash);
  fTree->Branch("peFlash",  &fPeFlash);
  fTree->Branch("tofHit",   &fTofHit);
  fTree->Branch("tofFlash", &fTofFlash);
  fTree->Branch("tofPe",    &fTofPe);
  fTree->Branch("tTTrig",   &fTTrig,    "tTrig/D");

  fClock = lar::providerFrom<detinfo::DetectorClocksService>();
}

void CrtOpHitMatchAnalysis::analyze(art::Event const& e)
{
  fTTrig = fClock->TriggerTime();

  fTCrt.clear();
  fTOp.clear();
  fTFlash.clear();
  fPeFlash.clear();
  fTofHit.clear();
  fTofFlash.clear();
  fTofPe.clear();

  //OpHits
  art::Handle< std::vector<recob::OpHit> > opHitListHandle;
  std::vector< art::Ptr<recob::OpHit> >    opHitList;
  if( e.getByLabel(fOpHitModuleLabel,opHitListHandle) )
      art::fill_ptr_vector(opHitList, opHitListHandle);

  fNOpHit = opHitList.size();

  for(auto const& ophit : opHitList){
        fTOp.push_back(ophit->PeakTime()*1e3-fOpDelay);
  }

  //OpFlash
  art::Handle< std::vector<recob::OpFlash> > opFlashListHandle0;
  art::Handle< std::vector<recob::OpFlash> > opFlashListHandle1;
  std::vector< art::Ptr<recob::OpFlash> >    opFlashList;
  if( e.getByLabel(fOpFlashModuleLabel0,opFlashListHandle0) )
      art::fill_ptr_vector(opFlashList, opFlashListHandle0);
  if( e.getByLabel(fOpFlashModuleLabel1,opFlashListHandle1) )
      art::fill_ptr_vector(opFlashList, opFlashListHandle1);

  fNOpFlash = opFlashList.size();

  for(auto const& flash : opFlashList){
        fTFlash.push_back(flash->Time()*1e3-fOpDelay);
        fPeFlash.push_back(flash->TotalPE());
  }

  //CRTHits
  art::Handle< std::vector<icarus::crt::CRTHit> > crtHitListHandle;
  std::vector< art::Ptr<icarus::crt::CRTHit> >    crtHitList;
  if( e.getByLabel(fCrtHitModuleLabel,crtHitListHandle))
      art::fill_ptr_vector(crtHitList, crtHitListHandle);

  fNCrtHit = crtHitList.size();

  for(int icrt=0; icrt<fNCrtHit; icrt++){

      auto const& crthit = crtHitList[icrt];
      fTCrt.push_back((int32_t)crthit->ts0_ns + 1.1e6);

      double diff = DBL_MAX;
      int imatch=-1;
      for(int iflash=0; iflash<fNOpFlash; iflash++){
          if(abs(fTCrt.back()-fTFlash[iflash])<abs(diff)) {
		diff = fTCrt.back()-fTFlash[iflash];
                imatch = iflash;
          }
      }
      if(imatch>-1){
        fTofFlash.push_back(diff);
        fTofPe.push_back(fPeFlash[imatch]);
      }

      diff = DBL_MAX;
      imatch=-1;
      for(int iophit=0; iophit<fNOpHit; iophit++) {
          if(abs(fTCrt.back()-fTOp[iophit])<abs(diff)) {
                diff = fTCrt.back()-fTOp[iophit];
                imatch = iophit;
          }
      }
      if(imatch<0) continue;
      fTofHit.push_back(diff);
      //fTofPe.push_back(fPeFlash[imatch])'
  }

  fTree->Fill();
}

void CrtOpHitMatchAnalysis::beginJob()
{
}

void CrtOpHitMatchAnalysis::endJob()
{
}

DEFINE_ART_MODULE(CrtOpHitMatchAnalysis)
