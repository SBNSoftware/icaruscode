////////////////////////////////////////////////////////////////////////
// Class:       CRTAutoVeto
// Plugin Type: analyzer (art v3_05_01)
// File:        CRTAutoVeto_module.cc
//
// Generated at Sat Jun 13 10:33:50 2020 by Christopher Hilgenberg using cetskelgen
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

//larsoft includes
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

//C++ includes
#include <vector>
#include <map>
#include <string>
#include <utility>

//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

//local includes
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"


namespace icarus {
 namespace crt {
    class CRTAutoVeto;
 }
}

using namespace icarus::crt;
using std::vector;
using std::map;

class icarus::crt::CRTAutoVeto : public art::EDAnalyzer {

 public:
    explicit CRTAutoVeto(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
  
    // Plugins should not be copied or assigned.
    CRTAutoVeto(CRTAutoVeto const&) = delete;
    CRTAutoVeto(CRTAutoVeto&&) = delete;
    CRTAutoVeto& operator=(CRTAutoVeto const&) = delete;
    CRTAutoVeto& operator=(CRTAutoVeto&&) = delete;
  
    // Required functions.
    void analyze(art::Event const& e) override;

 private:

    //config vars
    art::InputTag fGenLabel;
    art::InputTag fSimLabel;
    art::InputTag fCRTHitLabel;
    //define fiducial volume
    double fFidXOut;
    double fFidXIn;
    double fFidYTop;
    double fFidYBot;
    double fFidZUp;
    double fFidZDown;

    //histos
    TH2F* fNuVertexXZHist;
    TH2F* fNuVertexXYHist;
    TH2F* fNuVertexYZHist;
    TH2F* fNuVertexXZAVHist;
    TH2F* fNuVertexXYAVHist; 
    TH2F* fNuVertexYZAVHist; 
    TH2F* fNuVertexXZFVHist; 
    TH2F* fNuVertexXYFVHist; 
    TH2F* fNuVertexYZFVHist; 
                  
    TH1F* fNuE; 
    TH1F* fNuEAV;
    TH1F* fNuEAVNotFV;
    TH1F* fNuEFV;
    TH1F* fNuEVeto;
    TH1F* fNuEVetoAV;
    TH1F* fNuEVetoAVNotFV;
    TH1F* fNuEVetoFV;

    TH1F* fCRTReg;
    TH1F* fCRTRegAV;
    TH1F* fCRTRegAVNotFV;
    TH1F* fCRTRegFV;
    TH1F* fNCRTReg;
    TH1F* fNCRTRegAV;
    TH1F* fNCRTRegAVNotFV;
    TH1F* fNCRTRegFV;

    //tree
    TTree* fTree;
    int   fEvent;
    int   fNumNu;
    int   fNuPDG;
    float fNuE;
    bool  fNuCC;
    bool  fNuAV;
    bool  fNuFV;
    int                   fNumHit;  //number of CRTHits
    vector<int>           fHitRegs; //region code of each hit
    vector<vector<float>> fHitXYZT; //position/time of each hit vector<{x,y,z,t}>

   
};//end class 


//constructor
CRTAutoVeto::CRTAutoVeto(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
    fGenLabel(p.get<art::InputTag>("GenLabel","generator")),
    fSimLabel(p.get<art::InputTag>("SimLabel","largeant")),
    fCRTHitLabel(p.get<art::InputTag>("CRTHitLabel","crthit")),
    fFidXOut(p.get<double>("FiducialXOuter", 25.0)),
    fFidXIn(p.get<double>("FiducialXInner", 0.0)),
    fFidYTop(p.get<double>("FiducialYTop", 25.0)),
    fFidYBot(p.get<double>("FiducialYBottom",25.0)),
    fFidZUp(p.get<double>("FiducialZUpstream",30.0)),
    fFidZDown(p.get<double>("FiducialZDownstream",50.0))
{

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
  
    //2D nu vertex pos
    fNuVertexXZHist = tfs->make<TH2F>("nu_vtx_XZ", "#nu vertex: XZ; Z [cm]; X [cm]",
                                     100, -1000, 1000, 100, -450, 450);
    fNuVertexXYHist = tfs->make<TH2F>("nu_vtx_XY", "#nu vertex: XY; X [cm]; Y [cm]",
                                     100, -450, 450, 100, -250, 250);
    fNuVertexYZHist = tfs->make<TH2F>("nu_vtx_YZ", "#nu vertex: YZ; Z [cm]; Y[cm]",
                               100, -1000, 1000, 100, -250, 250);
    fNuVertexXZAVHist = tfs->make<TH2F>("nu_vtx_AV_XZ", "#nu vertex: XZ; Z [cm]; X [cm]",
                               100, -1000, 1000, 100, -450, 450);
    fNuVertexXYAVHist = tfs->make<TH2F>("nu_vtx_AV_XY", "#nu vertex: XY; X [cm]; Y [cm]",
                               100, -450, 450, 100, -250, 250);
    fNuVertexYZAVHist = tfs->make<TH2F>("nu_vtx_AV_YZ", "#nu vertex: YZ; Z [cm]; Y[cm]",
                               100, -1000, 1000, 100, -250, 250);
    fNuVertexXZFVHist = tfs->make<TH2F>("nu_vtx_FV_XZ", "#nu vertex: XZ; Z [cm]; X [cm]",
                               100, -1000, 1000, 100, -450, 450);
    fNuVertexXYFVHist = tfs->make<TH2F>("nu_vtx_FV_XY", "#nu vertex: XY; X [cm]; Y [cm]",
                               100, -450, 450, 100, -250, 250);
    fNuVertexYZFVHist = tfs->make<TH2F>("nu_vtx_FV_YZ", "#nu vertex: YZ; Z [cm]; Y[cm]",
                               100, -1000, 1000, 100, -250, 250);
  
    //nu info
    fNuE            = tfs->make<TH1F>("nuE",           "#nu energy; E_{#nu} [GeV]",100,0,10);
    fNuEAV          = tfs->make<TH1F>("nuEAV",         "#nu energy; E_{#nu} [GeV]",100,0,10);
    fNuEAVNotFV     = tfs->make<TH1F>("nuEAVNotFV",    "#nu energy; E_{#nu} [GeV]",100,0,10);
    fNuEFV          = tfs->make<TH1F>("nuEFV",         "#nu energy; E_{#nu} [GeV]",100,0,10);
    fNuEVeto        = tfs->make<TH1F>("nuEVeto",       "#nu energy; E_{#nu} [GeV]",100,0,10);
    fNuEVetoAV      = tfs->make<TH1F>("nuEVetoAV",     "#nu energy; E_{#nu} [GeV]",100,0,10);
    fNuEVetoAVNotFV = tfs->make<TH1F>("nuEVetoAVNotFV","#nu energy; E_{#nu} [GeV]",100,0,10);
    fNuEVetoFV      = tfs->make<TH1F>("nuEVetoFV",     "#nu energy; E_{#nu} [GeV]",100,0,10);
  
    //CRT info
    fCRTReg         = tfs->make<TH1F>("crtReg",        "CRT Region; region ID",23,29,52);
    fCRTRegAV       = tfs->make<TH1F>("crtRegAV",      "CRT Region; region ID",23,29,52);
    fCRTRegAVNotFV  = tfs->make<TH1F>("crtRegAVNotFV" ,"CRT Region; region ID",23,29,52);
    fCRTRegFV       = tfs->make<TH1F>("crtRegFV",      "CRT Region; region ID",23,29,52);
    fNCRTReg        = tfs->make<TH1F>("nCRTReg",       "No. CRT Regions Hit; no. regions",4,0,4);
    fNCRTRegAV      = tfs->make<TH1F>("nCRTRegAV",     "No. CRT Regions Hit; no. regions",4,0,4);
    fNCRTRegAVNotFV = tfs->make<TH1F>("nCRTRegAVNotFV","No. CRT Regions Hit; no. regions",4,0,4);
    fNCRTRegFV      = tfs->make<TH1F>("nCRTRegFV",     "No. CRT Regions Hit; no. regions",4,0,4);

    //Tree
    fTree  = tfs->make<TTree>("VetoTree", "auto-veto info");
    fTree->Branch("event",   &fEvent,  "event/I");
    fTree->Branch("numNu",   &fNumNu,  "numNu/I");
    fTree->Branch("nuPDG",   &fNuPDG,  "nuPDG/I");
    fTree->Branch("nuE",     &fNuE,    "nuE/F");
    fTree->Branch("cc",      &fNuCC,   "cc/O");
    fTree->Branch("av",      &fNuAV,   "av/O");
    fTree->Branch("fv",      &fNuFV,   "fv/O");
    fTree->Branch("numHit",  &fNumHit, "numHit/I");
    fTree->Branch("hitReg",  &fHitRegs);
    fTree->Branch("hitXYZT", &fHitXYZT);

}

void CRTAutoVeto::analyze(art::Event const& ev)
{
    fEvent = ev.id().event();
  
    auto const& mctruths = //vector of MCTruths
      *ev.getValidHandle<vector<simb::MCTruth>>(fGenLabel);
  
    auto const& crthits = //vector of CRTHits
      *ev.getValidHandle<vector<CRTHit>>(fCRTHitLabel);
  
    fNumNu = mctruths.size();
    fNumHit = crthits.size();
  
    //loop over neutrinos
    for(auto const& mctruth : mctruths) {
    
        auto const& nu = mctruth.GetNeutrino();
  
        fNuPDG = nu.Nu().PDGCode();
        fNuE   = nu.Nu().E();
  
        fNuCC = false;
        if(nu.CCNC()==0)
            fNuCC = true;
  
  
  
    }//for nus

    CRTCommonUtils util(); 
    fHitRegs.clear();
    fHitXYZT.clear();
  
    //loop over CRTHits
    for(auto const& hit : crthits) {
  
       fHitRegs.push_back(util.AuxDetRegionNameToNum(hit.tagger)); 
       fHitXYZT.push_back({hit.x_pos,hit.y_pos,hit.z_pos,hit.ts0_ns});
  
    }

    fTree->Fill();

}//analyze

DEFINE_ART_MODULE(CRTAutoVeto)
