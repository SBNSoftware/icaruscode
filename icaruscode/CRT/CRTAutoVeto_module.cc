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
#include "art/Framework/Services/Registry/ServiceHandle.h"

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
#include <TVector3.h>

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

    bool IsAV(TVector3 const& point);
    bool IsFV(TVector3 const& point);

    //config vars
    art::InputTag fGenLabel;
    art::InputTag fSimLabel;
    art::InputTag fCRTHitLabel;

    //define fiducial volume
    double fFidXOut; //distance from outer drift LAr face, outer TPCs [cm]
    double fFidXIn; //distance from outer drift LAr face, inner TPCs [cm]
    double fFidYTop; //distance from top active LAr face    [cm]
    double fFidYBot; //distance from bottom active LAr face [cm]
    double fFidZUp;  //distance from upstream active LAr face [cm]
    double fFidZDown; //distance from downstream active LAr face [cm]

    //utilities
    CRTBackTracker* bt;          //CRTBackTracker extracts truth info from CRT hits
    CRTCommonUtils* crtutil;     //CRT related utilities 
    geo::GeometryCore const* fGeoService; //access geometry

    //tree
    TTree* fTree; //pointer to analysis tree
    int   fEvent; //art event number
    int   fNumNu; //number of neutrinos in this event
    int   fNuPDG; //neutrino PDG code
    float fNuE;  //nu energy [GeV]
    bool  fNuCC; //is the nu interaction CC?
    bool  fNuAV; //was nu vertex in active volume?
    bool  fNuFV; //was nu vertex in fiducial volume?
    int                   fNumHit;  //number of CRTHits
    vector<int>           fHitRegs; //region code of each hit
    vector<vector<float>> fHitXYZT; //position/time of each hit vector<{x,y,z,t}>
    vector<vector<int>>   fHitPDGs; //PDG code of all particles generating a CRTHit
    vector<int>           fNHitPDGs; //number of particles contributing to CRT hit
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
    fTree->Branch("hitPDG",  &fHitPDGs);
    fTree->Branch("nHitPDG", &fNHitPDGs);
}

void CRTAutoVeto::analyze(art::Event const& ev)
{
  
    auto const& mctruths = //vector of MCTruths from GENIE
      *ev.getValidHandle<vector<simb::MCTruth>>(fGenLabel);

    auto const& simparticles = //vector of MCParticles from G4
      *ev.getValidHandle<vector<simb::MCParticle>>(fSimLabel);
  
    auto const& crthits = //vector of CRTHits
      *ev.getValidHandle<vector<CRTHit>>(fCRTHitLabel);

    map< int, const simb::MCParticle*> particleMap;

    fEvent = ev.id().event();  
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
  
        const TVector3 point = mctruth.GetNeutrino().Nu().Position(0).Vect();
        fNuAV = IsAV(point);
        fNuFV = false;
        if(fNuAV)
            fNuFV = IsFV(point);  
  
    }//for nus

    //loop over G4 tracks
    for(auto const& particle : simparticles){

        particleMap[particle.TrackId()] = &particle;

    }//G4 tracks

    fHitRegs.clear();
    fHitXYZT.clear();
    fHitPDGs.clear();
    fNHitPDGs.clear();

    //loop over CRTHits
    for(auto const& hit : crthits) {
  
       fHitRegs.push_back(crtutil->AuxDetRegionNameToNum(hit.tagger)); 
       fHitXYZT.push_back({hit.x_pos,hit.y_pos,hit.z_pos,(float)hit.ts0_ns});
       fHitPDGs.push_back({});
 
       int npdg=0;
       for(const int id :  bt->AllTrueIds(ev,hit)){
           if(id<0) continue;
           if(particleMap.find(id)==particleMap.end()) continue;
           fHitPDGs.back().push_back(particleMap[id]->PdgCode());
           npdg++;
       }
       fNHitPDGs.push_back(npdg);

    }

    fTree->Fill();

}//analyze

bool CRTAutoVeto::IsAV(TVector3 const& point){

    geo::CryostatGeo const& cryo0 = fGeoService->Cryostat(0);
    geo::CryostatGeo const& cryo1 = fGeoService->Cryostat(1);
    geo::TPCGeo const& tpc00 = cryo0.TPC(0);
    geo::TPCGeo const& tpc01 = cryo0.TPC(1);
    geo::TPCGeo const& tpc10 = cryo1.TPC(0);
    geo::TPCGeo const& tpc11 = cryo1.TPC(1);

    if(tpc00.ContainsPosition(point) ||
       tpc01.ContainsPosition(point) ||
       tpc10.ContainsPosition(point) ||
       tpc11.ContainsPosition(point) )
        return true; 

    return false;
}

bool CRTAutoVeto::IsFV(TVector3 const& point) {

    geo::CryostatGeo const& cryo0 = fGeoService->Cryostat(0);
    geo::CryostatGeo const& cryo1 = fGeoService->Cryostat(1);
    geo::TPCGeo const& tpc00 = cryo0.TPC(0);
    geo::TPCGeo const& tpc01 = cryo0.TPC(1);
    geo::TPCGeo const& tpc10 = cryo1.TPC(0);
    geo::TPCGeo const& tpc11 = cryo1.TPC(1);

    if(tpc00.ContainsPosition(point)
        && tpc00.InFiducialX(point.X(),fFidXOut,fFidXIn)
        && tpc00.InFiducialY(point.Y(),fFidYBot,fFidYTop)
        && tpc00.InFiducialZ(point.Z(),fFidZUp,fFidZDown) )
        return true;
    if(tpc01.ContainsPosition(point)
        && tpc01.InFiducialX(point.X(),fFidXOut,fFidXIn)
        && tpc01.InFiducialY(point.Y(),fFidYBot,fFidYTop)
        && tpc01.InFiducialZ(point.Z(),fFidZUp,fFidZDown) )
        return true;
    if(tpc10.ContainsPosition(point)
        && tpc10.InFiducialX(point.X(),fFidXOut,fFidXIn)
        && tpc10.InFiducialY(point.Y(),fFidYBot,fFidYTop)
        && tpc10.InFiducialZ(point.Z(),fFidZUp,fFidZDown) )
        return true;
    if(tpc11.ContainsPosition(point)
        && tpc11.InFiducialX(point.X(),fFidXIn,fFidXOut)
        && tpc11.InFiducialY(point.Y(),fFidYBot,fFidYTop)
        && tpc11.InFiducialZ(point.Z(),fFidZUp,fFidZDown) )
        return true;

    return false;
}

DEFINE_ART_MODULE(CRTAutoVeto)
