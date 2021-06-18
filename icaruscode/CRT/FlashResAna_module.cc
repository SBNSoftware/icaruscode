////////////////////////////////////////////////////////////////////////
// Class:       FlashResAna
// Plugin Type: analyzer (art v3_05_01)
// File:        FlashResAna_module.cc
//
// Generated at Sun Jun 28 11:06:51 2020 by Christopher Hilgenberg using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

//framework includes
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

//larsoft includes
//#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimPhotons.h"

//C++ includes
#include <vector>
#include <map>

//ROOT includes
#include <TTree.h>
#include <TH1F.h>

namespace icarus {
    class FlashResAna;
}

using std::vector;
using std::map;

class icarus::FlashResAna : public art::EDAnalyzer {
public:
  explicit FlashResAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlashResAna(FlashResAna const&) = delete;
  FlashResAna(FlashResAna&&) = delete;
  FlashResAna& operator=(FlashResAna const&) = delete;
  FlashResAna& operator=(FlashResAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  art::InputTag fGenLabel;
  art::InputTag fPhotLabel;
  art::InputTag fFlashLabel0;
  art::InputTag fFlashLabel1;
  art::InputTag fFlashLabel2;
  art::InputTag fFlashLabel3;
  art::InputTag fHitLabel;

  map<int,art::InputTag> fFlashLabels;

  TTree* fTree;
  int    fNFlash;
  bool   fCC;
  double fFlashPE;
  int    fFlashTPC;
  double fFlashWidth[4];
  int    fMaxChan;
  double fMaxPE;
  double fNuE;
  double fNuXYZT[4];
  double fFlashXYZT[4];
  double fFlashXYZTMin[4];
  double fFlashXYZTMax[4];
  double fDelta[4];
  double fDeltaPmt;
  bool   fInFrame;
  int    fNPmtFlash;
  double fHitXYZT[4];
  double fHitDist;
  double fHitTof;
  int    fNHit;
  double fHitPE;
  double fTPhot;

  TH1F* fDeltaTFlash;

  unsigned int GetMaxChan(vector<double> const& pes);

  geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider

};


icarus::FlashResAna::FlashResAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fGenLabel(p.get<art::InputTag>("GenLabel","generator")),
    fPhotLabel(p.get<art::InputTag>("PhotLabel","largeant")),
    fFlashLabel0(p.get<art::InputTag>("FlashLabel0","opflashTPC0")),
    fFlashLabel1(p.get<art::InputTag>("FlashLabel1","opflashTPC1")),
    fFlashLabel2(p.get<art::InputTag>("FlashLabel2","opflashTPC2")),
    fFlashLabel3(p.get<art::InputTag>("FlashLabel3","opflashTPC3")),
    fHitLabel(p.get<art::InputTag>("HitLabel","ophit"))
{
    fFlashLabels[0] = fFlashLabel0;
    fFlashLabels[1] = fFlashLabel1;
    fFlashLabels[2] = fFlashLabel2;
    fFlashLabels[3] = fFlashLabel3;

    art::ServiceHandle<art::TFileService> tfs;

    fTree = tfs->make<TTree>("anatree","OpFlash resolution");
    fTree->Branch("nFlash",    &fNFlash,   "nFlash/I");
    fTree->Branch("nuE",       &fNuE,      "nuE/D");
    fTree->Branch("cc",        &fCC,       "cc/O");
    fTree->Branch("nuXYZT",    fNuXYZT,    "nuXYZT[4]/D");
    fTree->Branch("flashXYZT", fFlashXYZT, "flashXYZT[4]/D");
    fTree->Branch("flashXYZTMin", fFlashXYZTMin, "flashXYZTMin[4]/D");
    fTree->Branch("flashXYZTMax", fFlashXYZTMax, "flashXYZTMax[4]/D");
    fTree->Branch("flashWidth",fFlashWidth,"flashWidth[4]/D");
    fTree->Branch("flashPE",   &fFlashPE,  "flashPE/D");
    fTree->Branch("flashTPC",  &fFlashTPC, "flashTPC/I");
    fTree->Branch("maxChan",   &fMaxChan,  "maxChan/i");
    fTree->Branch("maxPE",     &fMaxPE,    "maxPE/D");
    fTree->Branch("delta",     fDelta,     "delta[4]/D");
    fTree->Branch("inFrame",   &fInFrame,  "inFrame/O");
    fTree->Branch("deltaPmt",  &fDeltaPmt, "deltaPmt/D");
    fTree->Branch("nPmtFlash", &fNPmtFlash,"nPmtFlash/I");
    fTree->Branch("hitXYZT",   fHitXYZT,   "hittXYZT[4]/D");
    fTree->Branch("hitDist",   &fHitDist,  "hitDist/D");
    fTree->Branch("hitTof",    &fHitTof,   "hitTof/D");
    fTree->Branch("nHit",      &fNHit,     "nHit/I");
    fTree->Branch("hitPE",     &fHitPE,    "hitPE/D"); 
    fTree->Branch("tPhot",     &fTPhot,    "tPhot/D");

    fDeltaTFlash = tfs->make<TH1F>("dtflash","delay between flashes;#Delta t [ns];",25000,0,500000);

    fGeometryService = lar::providerFrom<geo::Geometry>();

}

void icarus::FlashResAna::analyze(art::Event const& e)
{
    //initialize for each event
    fNFlash   = 0;
    fCC       = false;
    fFlashPE  = 0.;
    fFlashTPC = -1;
    fMaxChan  = -1;
    fMaxPE    = 0.;
    fNuE      = 0.;
    fDeltaPmt = 0.;
    fInFrame  = false;;
    fNPmtFlash= -1;
    double dt = DBL_MAX;
    geo::CryostatGeo const& cryo0 = fGeometryService->Cryostat(0);
 
    for(size_t i=0; i<4; i++) {
        fNuXYZT[i] = dt;
        fFlashXYZT[i] = dt;
        fFlashXYZTMin[i] = dt;
        fFlashXYZTMax[i] = dt;
        fFlashWidth[i] = dt;
        fDelta[i] = dt;
    }


    //neutrino truth
    auto const& mctruths = //vector of MCTruths from GENIE
      *e.getValidHandle<vector<simb::MCTruth>>(fGenLabel);

    auto const& nu = mctruths[0].GetNeutrino();
    fNuE = nu.Nu().E();
    if(nu.CCNC()==0)
        fCC = true;

    const TLorentzVector nuxyzt = nu.Nu().Position(0);
    fNuXYZT[0] = nuxyzt.X();
    fNuXYZT[1] = nuxyzt.Y();
    fNuXYZT[2] = nuxyzt.Z();
    fNuXYZT[3] = nuxyzt.T();

    //SimPhot
    art::Handle< std::vector<sim::SimPhotons>> photHandle;
    std::vector< art::Ptr<sim::SimPhotons> >   photList;
    if( e.getByLabel(fPhotLabel, photHandle) )
        art::fill_ptr_vector(photList, photHandle);

    double tphotmin = DBL_MAX;
    for(auto const& phot : photList){
        for(size_t i=0; i<phot->size(); i++) {
            if ((phot->at(i)).Time < tphotmin)
                tphotmin = (phot->at(i)).Time;
        }
    }
    fTPhot = tphotmin;

    //OpFlash
    map<int,art::Handle< std::vector<recob::OpFlash> > > flashHandles;
    map<int,vector< art::Ptr<recob::OpFlash> >> flashLists;
    if( e.getByLabel(fFlashLabel0,flashHandles[0]) )
        art::fill_ptr_vector(flashLists[0], flashHandles[0]);
    if( e.getByLabel(fFlashLabel1,flashHandles[1]) )
        art::fill_ptr_vector(flashLists[1], flashHandles[1]);
    if( e.getByLabel(fFlashLabel2,flashHandles[2]) )
        art::fill_ptr_vector(flashLists[2], flashHandles[2]);
    if( e.getByLabel(fFlashLabel3,flashHandles[3]) )
        art::fill_ptr_vector(flashLists[3], flashHandles[3]);

    for(auto flashList : flashLists) {

        int tpc = flashList.first;
        size_t nflashtpc = flashList.second.size();
        if(nflashtpc==0) continue;
        fNFlash += nflashtpc;

        std::sort(flashList.second.begin(),flashList.second.end(),
                  [](art::Ptr<recob::OpFlash> f1, art::Ptr<recob::OpFlash> f2)
                     {return f1->Time()<f2->Time();}
                 );

        art::FindManyP<recob::OpHit> findManyHits(
                          flashHandles[tpc], e, fFlashLabels[tpc]);

        //auto const& flash0 = flashList.second[0];
        for(size_t iflash=0; iflash < nflashtpc; iflash++) {
            auto const& flash0 = flashList.second[iflash]; //) {
            if(abs(fNuXYZT[3]-flash0->Time()*1.0e3)<abs(dt)) {
                dt             = flash0->Time()*1.0e3 - fNuXYZT[3];
                fFlashXYZT[0]  = 0.;
                fFlashXYZT[1]  = flash0->YCenter();
                fFlashXYZT[2]  = flash0->ZCenter();
                fFlashXYZT[3]  = flash0->Time()*1.0e3;
                fFlashWidth[0] = 0.;
                fFlashWidth[1] = flash0->YWidth();
                fFlashWidth[2] = flash0->ZWidth();
                fFlashWidth[3] = flash0->TimeWidth();
                fFlashPE       = flash0->TotalPE();
                fFlashTPC      = tpc;
                auto const& pes= flash0->PEs();
                fNPmtFlash     = pes.size();
                fMaxChan       = GetMaxChan(pes);
                fMaxPE         = flash0->PE(fMaxChan);
                fInFrame       = flash0->InBeamFrame();
                for(size_t i=0; i<4; i++)
                    fDelta[i] = fFlashXYZT[i] - fNuXYZT[i];

                double tPmtMax=-DBL_MAX, tPmtMin=DBL_MAX;
                vector<art::Ptr<recob::OpHit>> hits = findManyHits.at(iflash);
                for(auto const& hit : hits) {
                    if(hit->OpChannel() == fMaxChan){
                        fDeltaPmt = hit->PeakTime()*1.0e3 - fNuXYZT[3];
                    }
                    if(hit->PeakTime()<tPmtMin){
                        tPmtMin = hit->PeakTime();
                        geo::OpDetGeo const& opDet = cryo0.OpDet(hit->OpChannel());
                        double pos[3];
                        opDet.GetCenter(pos);
                        for(int i=0; i<3; i++) fFlashXYZTMin[i] = pos[i];
                        fFlashXYZTMin[3] = tPmtMin;
                    }
                    if(hit->PeakTime()>tPmtMax){
                        tPmtMax = hit->PeakTime();
                        geo::OpDetGeo const& opDet = cryo0.OpDet(hit->OpChannel());
                        double pos[3];
                        opDet.GetCenter(pos);
                        for(int i=0; i<3; i++) fFlashXYZTMax[i] = pos[i];
                        fFlashXYZTMax[3] = tPmtMax;
                    }
                }
            }
        }

        if(nflashtpc>1)
        for(size_t i=0; i<nflashtpc-1; i++) {
              auto const& flash = flashList.second.at(i);
              auto const& flashnext = flashList.second.at(i+1);
              fDeltaTFlash->Fill((flashnext->Time()-flash->Time())*1.0e3);
        }
    }//for flash lists

    //OpHits
    art::Handle< std::vector<recob::OpHit> > opHitHandle;
    std::vector< art::Ptr<recob::OpHit> >    opHitList;

    if( e.getByLabel(fHitLabel,opHitHandle) ) {

        art::fill_ptr_vector(opHitList, opHitHandle);
        fNHit = opHitList.size();
        if(fNHit==0){
            std::cout << "empty OpHit vector!" << std::endl;
            for(int i=0; i<4; i++) fHitXYZT[i] = DBL_MAX;
            fHitDist = DBL_MAX;
            fHitTof = DBL_MAX;
            fHitPE = DBL_MAX; 
        }

        else{ 
            std::sort(opHitList.begin(),opHitList.end(),
                    [](art::Ptr<recob::OpHit> h1, art::Ptr<recob::OpHit> h2)
                       {return h1->PeakTime()<h2->PeakTime();}
                   );
            if(opHitList.size()>1 && opHitList[0]->PeakTime()>opHitList.back()->PeakTime())
                std::cout << "OpHits time sort failed!" << std::endl;

            geo::OpDetGeo const& opDet = cryo0.OpDet(opHitList[0]->OpChannel());
            double pos[3];
            opDet.GetCenter(pos);
            for(int i=0; i<3; i++) fHitXYZT[i] = pos[i];
            fHitXYZT[3] = opHitList[0]->PeakTime()*1.0e3; //time in ns
            fHitDist = 0.;
            for(int i=0; i<3; i++) 
                fHitDist+=pow(fHitXYZT[i]-fNuXYZT[i],2);
            fHitDist = sqrt(fHitDist);
            fHitTof  = fHitXYZT[3] - fNuXYZT[3];
            fHitPE = opHitList[0]->PE();
        }
    }    
    else{
        std::cout << "no OpHit handle found!" << std::endl;
    }

    fTree->Fill();

}

unsigned int icarus::FlashResAna::GetMaxChan(vector<double> const& pes){

    size_t imax = 0;
    double max = 0.;
    for(size_t i=0; i<pes.size(); i++){
        if(pes[i]>max) {
            max = pes[i];
            imax = i;
        }
    }

    return imax;
}

DEFINE_ART_MODULE(icarus::FlashResAna)
