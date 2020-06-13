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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

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

//namespace {

//}//local namespace

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

    //utilities
    CRTBackTracker* bt;
    CRTCommonUtils* crtutil;
    geo::GeometryCore const* fGeoService;

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
    fFidZDown(p.get<double>("FiducialZDownstream",50.0)),
    bt(new CRTBackTracker(p.get<fhicl::ParameterSet>("CRTBackTrack"))),
    crtutil(new CRTCommonUtils())
{
    fGeoService = lar::providerFrom<geo::Geometry>();

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

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
  
        fNuPDG = nu.Nu().PdgCode();
        fNuE   = nu.Nu().E();
  
        fNuCC = false;
        if(nu.CCNC()==0)
            fNuCC = true;
  
  
  
    }//for nus

    fHitRegs.clear();
    fHitXYZT.clear();
  
    //loop over CRTHits
    for(auto const& hit : crthits) {
  
       fHitRegs.push_back(crtutil->AuxDetRegionNameToNum(hit.tagger)); 
       fHitXYZT.push_back({hit.x_pos,hit.y_pos,hit.z_pos,(float)hit.ts0_ns});
  
    }

    fTree->Fill();

}//analyze

DEFINE_ART_MODULE(CRTAutoVeto)
