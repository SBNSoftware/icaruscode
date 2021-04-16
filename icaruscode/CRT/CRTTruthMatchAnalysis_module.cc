////////////////////////////////////////////////////////////////////////
// Class:       CRTTruthMatchAnalysis
// Plugin Type: analyzer (art v3_05_00)
// File:        CRTTruthMatchAnalysis_module.cc
//
// Generated at Fri Apr 24 14:45:31 2020 by Christopher Hilgenberg using cetskelgen
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

//LarSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

//icaruscode includes
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

//C++ includes
#include <vector>
#include <map>
#include <string>
#include <utility>

//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TVector3.h>

namespace icarus{ 
 namespace crt {
    class CRTTruthMatchAnalysis;
 }
}

using namespace icarus::crt;
using std::vector;
using std::map;

class icarus::crt::CRTTruthMatchAnalysis : public art::EDAnalyzer {

 public:
    using CRTHit = sbn::crt::CRTHit;
   
    explicit CRTTruthMatchAnalysis(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTTruthMatchAnalysis(CRTTruthMatchAnalysis const&) = delete;
    CRTTruthMatchAnalysis(CRTTruthMatchAnalysis&&) = delete;
    CRTTruthMatchAnalysis& operator=(CRTTruthMatchAnalysis const&) = delete;
    CRTTruthMatchAnalysis& operator=(CRTTruthMatchAnalysis&&) = delete;

    void beginJob() override;
    void analyze(art::Event const& e) override;

 private:

    float GetDeltaR(art::Ptr<CRTHit> h1, art::Ptr<CRTHit> h2);
    int GetDeltaT(art::Ptr<CRTHit> h1, art::Ptr<CRTHit> h2);
    TVector3 HitToPosVec(art::Ptr<CRTHit> hit);

    art::InputTag fSimulationLabel;
    art::InputTag fCRTTrueHitLabel;
    art::InputTag fCRTDataLabel;
    art::InputTag fCRTSimHitLabel;

    CRTBackTracker bt;
    CRTCommonUtils* fCrtutils;

    const float fPosmax;
    const float fPosbinning;
    const float fTimemax;
    const float fTbinning;

    const bool fVerbose;

    TH1F* fNTrueHitMu; //number of CRTTrueHits per muon track
    TH1F* fNDataMu;    //number of CRTData per muon track
    TH1F* fNSimHitMu;   //number of CRTSimHits per muon track
    TH1F* fNTrueSimDiffMu; //difference in number of true - sim hits

    TH1F* fNTrueHitMu_top;
    TH1F* fNTrueHitMu_side;
    TH1F* fNTrueHitMu_bot;
//     TH1F* fNDataMu_top;
//     TH1F* fNDataMu_side;
//     TH1F* fNDataMu_bot;
    TH1F* fNSimHitMu_top;
    TH1F* fNSimHitMu_side;
    TH1F* fNSimHitMu_bot;

    TH1F* fNRegTrue;
    TH1F* fNRegSim;

    TH1F* fFebsTrue;
    TH1F* fFebsSim;

    TH1F *fRres_top, *fXres_top, *fYres_top, *fZres_top, *fTres_top;
    TH1F *fRres_rim, *fXres_rim, *fYres_rim, *fZres_rim, *fTres_rim;
    TH1F *fRres_lat, *fXres_lat, *fYres_lat, *fZres_lat, *fTres_lat;
    TH1F *fRres_nor, *fXres_nor, *fYres_nor, *fZres_nor, *fTres_nor;
    TH1F *fRres_sth, *fXres_sth, *fYres_sth, *fZres_sth, *fTres_sth;
    TH1F *fRres_bot, *fXres_bot, *fYres_bot, *fZres_bot, *fTres_bot;

    //size_t        fNViews;
    //vector<TH2F*> fChargeViews;
    //vector<TGraph*> fTrueHitPoints;
    //vector<TGraph*> fSimHitPoints;

};//class definition

//constructor
CRTTruthMatchAnalysis::CRTTruthMatchAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSimulationLabel(p.get<art::InputTag>("SimulationLabel","largeant")),
    fCRTTrueHitLabel(p.get<art::InputTag>("CRTTrueHitLabel","crttruehit")),
    fCRTDataLabel(p.get<art::InputTag>("CRTDataLabel","crtdaq")),
    fCRTSimHitLabel(p.get<art::InputTag>("CRTSimHitLabel","crthit")),
    bt(p.get<fhicl::ParameterSet>("CRTBackTrack")),
    fCrtutils(new CRTCommonUtils()),
    fPosmax(p.get<float>("PosMax",1000.)),
    fPosbinning(p.get<float>("PosBinning",1.)),
    fTimemax(p.get<float>("TimeMax",30.)),
    fTbinning(p.get<float>("Tbinning",1.)),
    fVerbose(p.get<bool>("Verbose",false))
{}

void CRTTruthMatchAnalysis::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    fNTrueHitMu      = tfs->make<TH1F>("hNTrueHitMu",     "No. True Hits per #mu",      10,0,10);
    fNDataMu         = tfs->make<TH1F>("hNDataMu",        "No. FEB Triggers per #mu",   10,0,10);
    fNSimHitMu       = tfs->make<TH1F>("hNSimHitMu",      "No. Sim Hits per #mu",       10,0,10);
    fNTrueSimDiffMu  = tfs->make<TH1F>("hNTrueSimDiffMu", "No. True-Sim Hits per #mu",  14,-6,8); 

    fNTrueHitMu_top  = tfs->make<TH1F>("hNTrueHitMu_top", "No. True Hits per #mu",      10,0,10);
    fNTrueHitMu_side = tfs->make<TH1F>("hNTrueHitMu_side","No. True Hits per #mu",      10,0,10);
    fNTrueHitMu_bot  = tfs->make<TH1F>("hNTrueHitMu_bot", "No. True Hits per #mu",      10,0,10);
    fNSimHitMu_top   = tfs->make<TH1F>("hNSimHitMu_top",  "No. Sim Hits per #mu",       10,0,10);
    fNSimHitMu_side  = tfs->make<TH1F>("hNSimHitMu_side", "No. Sim Hits per #mu",       10,0,10);
    fNSimHitMu_bot   = tfs->make<TH1F>("hNSimHitMu_bot",  "No. Sim Hits per #mu",       10,0,10);

    fNRegTrue        = tfs->make<TH1F>("hNRegTrue",       "No. True Hits per Region",   10,0,10);
    fNRegSim         = tfs->make<TH1F>("hNRegSim",        "No. Sim Hits per Region",    10,0,10); 

    fFebsTrue        = tfs->make<TH1F>("hFebsTrue",       "mac5s in true hits",         300,0,300);
    fFebsSim         = tfs->make<TH1F>("hFebsSim",        "mac5s in sim hits",          300,0,300);

    const float absmax = fPosmax+0.5*fPosbinning;
    const float absmin = -1*absmax;
    const int nbinsabs = 2*absmax/fPosbinning;

    const float rmax = fPosmax;
    const int nbinsr = rmax/fPosbinning;

    const float tmax = fTimemax+0.5*fTbinning;
    const float tmin = -1*tmax;
    const int nbinst = 2*tmax/fTbinning;

    art::TFileDirectory resdir = tfs->mkdir("resolution");
    fRres_top        = resdir.make<TH1F>("hRres_top",       "CRTHit #sigma_{r}: Top",       nbinsr,0,rmax);
    fRres_rim        = resdir.make<TH1F>("hRres_rim",       "CRTHit #sigma_{r}: Rim",       nbinsr,0,rmax);
    fRres_bot        = resdir.make<TH1F>("hRres_bot",       "CRTHit #sigma_{r}: Bottom",    nbinsr,0,rmax);
    fRres_lat        = resdir.make<TH1F>("hRres_lat",       "CRTHit #sigma_{r}: East/West", nbinsr,0,rmax);
    fRres_nor        = resdir.make<TH1F>("hRres_nor",       "CRTHit #sigma_{r}: North",     nbinsr,0,rmax);
    fRres_sth        = resdir.make<TH1F>("hRres_sth",       "CRTHit #sigma_{r}: South",     nbinsr,0,rmax);

    fXres_top        = resdir.make<TH1F>("hXres_top",       "CRTHit #sigma_{x}: Top",       nbinsabs,absmin,absmax);
    fXres_rim        = resdir.make<TH1F>("hXres_rim",       "CRTHit #sigma_{x}: Rim",       nbinsabs,absmin,absmax);
    fXres_bot        = resdir.make<TH1F>("hXres_bot",       "CRTHit #sigma_{x}: Bottom",    nbinsabs,absmin,absmax);
    fXres_lat        = resdir.make<TH1F>("hXres_lat",       "CRTHit #sigma_{x}: East/West", nbinsabs,absmin,absmax);
    fXres_nor        = resdir.make<TH1F>("hXres_nor",       "CRTHit #sigma_{x}: North",     nbinsabs,absmin,absmax);
    fXres_sth        = resdir.make<TH1F>("hXres_sth",       "CRTHit #sigma_{x}: South",     nbinsabs,absmin,absmax);

    fYres_top        = resdir.make<TH1F>("hYres_top",       "CRTHit #sigma_{y}: Top",       nbinsabs,absmin,absmax);
    fYres_rim        = resdir.make<TH1F>("hYres_rim",       "CRTHit #sigma_{y}: Rim",       nbinsabs,absmin,absmax);
    fYres_bot        = resdir.make<TH1F>("hYres_bot",       "CRTHit #sigma_{y}: Bottom",    nbinsabs,absmin,absmax);
    fYres_lat        = resdir.make<TH1F>("hYres_lat",       "CRTHit #sigma_{y}: East/West", nbinsabs,absmin,absmax);
    fYres_nor        = resdir.make<TH1F>("hYres_nor",       "CRTHit #sigma_{y}: North",     nbinsabs,absmin,absmax);
    fYres_sth        = resdir.make<TH1F>("hYres_sth",       "CRTHit #sigma_{y}: South",     nbinsabs,absmin,absmax);

    fZres_top        = resdir.make<TH1F>("hZres_top",       "CRTHit #sigma_{z}: Top",       nbinsabs,absmin,absmax);
    fZres_rim        = resdir.make<TH1F>("hZres_rim",       "CRTHit #sigma_{z}: Rim",       nbinsabs,absmin,absmax);
    fZres_bot        = resdir.make<TH1F>("hZres_bot",       "CRTHit #sigma_{z}: Bottom",    nbinsabs,absmin,absmax);
    fZres_lat        = resdir.make<TH1F>("hZres_lat",       "CRTHit #sigma_{z}: East/West", nbinsabs,absmin,absmax);
    fZres_nor        = resdir.make<TH1F>("hZres_nor",       "CRTHit #sigma_{z}: North",     nbinsabs,absmin,absmax);
    fZres_sth        = resdir.make<TH1F>("hZres_sth",       "CRTHit #sigma_{z}: South",     nbinsabs,absmin,absmax);

    fTres_top        = resdir.make<TH1F>("hTres_top",       "CRTHit #sigma_{t}: Top",       nbinst,tmin,tmax);
    fTres_rim        = resdir.make<TH1F>("hTres_rim",       "CRTHit #sigma_{t}: Rim",       nbinst,tmin,tmax);
    fTres_bot        = resdir.make<TH1F>("hTres_bot",       "CRTHit #sigma_{t}: Bottom",    nbinst,tmin,tmax);
    fTres_lat        = resdir.make<TH1F>("hTres_lat",       "CRTHit #sigma_{t}: East/West", nbinst,tmin,tmax);
    fTres_nor        = resdir.make<TH1F>("hTres_nor",       "CRTHit #sigma_{t}: North",     nbinst,tmin,tmax);
    fTres_sth        = resdir.make<TH1F>("hTres_sth",       "CRTHit #sigma_{t}: South",     nbinst,tmin,tmax);

    /*art::TFileDirectory viewdir = tfs->mkdir("chargeview");

    fNViews = 0;
    for(int i=0; i<100; i++) {
        std::string gnametrue = "gview_true"+std::to_string(i);
        std::string gnamesim  = "gview_sim"+std::to_string(i);
        std::string hname = "hview"+std::to_string(i);
        std::string gtitletrue = "Charge View True Point "+std::to_string(i);
        std::string gtitlesim  = "Charge View Sim Point "+std::to_string(i);
        std::string htitle = "Charge View " + std::to_string(i);
        fChargeViews.push_back(viewdir.make<TH2F>(hname.c_str(), htitle.c_str(),16,0,16,16,0,16));
        fChargeViews.back()->GetXaxis()->SetTitle("X-layer channel");
        fChargeViews.back()->GetYaxis()->SetTitle("Y-layer channel");
        fTrueHitPoints.push_back(viewdir.makeAndRegister<TGraph>(gnametrue.c_str(),gtitletrue.c_str(),1));
        fSimHitPoints.push_back(viewdir.makeAndRegister<TGraph>(gnamesim.c_str(),gtitlesim.c_str(),1));
        fTrueHitPoints.back()->SetMarkerStyle(47);
        fTrueHitPoints.back()->SetMarkerSize(10);
        fTrueHitPoints.back()->SetMarkerColor(kRed);
        fSimHitPoints.back()->SetMarkerStyle(8);
        fSimHitPoints.back()->SetMarkerSize(10);
        fSimHitPoints.back()->SetMarkerColor(kMagenta);
    }*/

    //style options
    fNTrueHitMu->SetLineWidth(2);
    fNTrueHitMu_top->SetLineWidth(2);
    fNTrueHitMu_side->SetLineWidth(2);
    fNTrueHitMu_bot->SetLineWidth(2);
    fNSimHitMu->SetLineWidth(2);
    fNSimHitMu_top->SetLineWidth(2);
    fNSimHitMu_side->SetLineWidth(2);
    fNSimHitMu_bot->SetLineWidth(2);

    fNTrueHitMu->SetLineColor(kBlack);
    fNTrueHitMu_top->SetLineColor(kBlue);
    fNTrueHitMu_side->SetLineColor(kRed);
    fNTrueHitMu_bot->SetLineColor(kGreen);
    fNSimHitMu->SetLineColor(kBlack);
    fNSimHitMu_top->SetLineColor(kBlue);
    fNSimHitMu_side->SetLineColor(kRed);
    fNSimHitMu_bot->SetLineColor(kGreen);

    fRres_top->SetLineWidth(2);
    fRres_rim->SetLineWidth(2);
    fRres_bot->SetLineWidth(2);
    fRres_lat->SetLineWidth(2);
    fRres_nor->SetLineWidth(2);
    fRres_sth->SetLineWidth(2);

    fXres_top->SetLineWidth(2);
    fXres_rim->SetLineWidth(2);
    fXres_bot->SetLineWidth(2);
    fXres_lat->SetLineWidth(2);
    fXres_nor->SetLineWidth(2);
    fXres_sth->SetLineWidth(2);

    fYres_top->SetLineWidth(2);
    fYres_rim->SetLineWidth(2);
    fYres_bot->SetLineWidth(2);
    fYres_lat->SetLineWidth(2);
    fYres_nor->SetLineWidth(2);
    fYres_sth->SetLineWidth(2);

    fZres_top->SetLineWidth(2);
    fZres_rim->SetLineWidth(2);
    fZres_bot->SetLineWidth(2);
    fZres_lat->SetLineWidth(2);
    fZres_nor->SetLineWidth(2);
    fZres_sth->SetLineWidth(2);

    fTres_top->SetLineWidth(2);
    fTres_rim->SetLineWidth(2);
    fTres_bot->SetLineWidth(2);
    fTres_lat->SetLineWidth(2);
    fTres_nor->SetLineWidth(2);
    fTres_sth->SetLineWidth(2);

    fRres_top->GetXaxis()->SetTitle("#||{#vec{r_{true}} - #vec{r_{reco}}} [cm]");
    fRres_rim->GetXaxis()->SetTitle("#||{#vec{r_{true}} - #vec{r_{reco}}} [cm]");
    fRres_bot->GetXaxis()->SetTitle("#||{#vec{r_{true}} - #vec{r_{reco}}} [cm]");
    fRres_lat->GetXaxis()->SetTitle("#||{#vec{r_{true}} - #vec{r_{reco}}} [cm]");
    fRres_nor->GetXaxis()->SetTitle("#||{#vec{r_{true}} - #vec{r_{reco}}} [cm]");
    fRres_sth->GetXaxis()->SetTitle("#||{#vec{r_{true}} - #vec{r_{reco}}} [cm]");

    fXres_top->GetXaxis()->SetTitle("#||{x_{true} - x_{reco}} [cm]");
    fXres_rim->GetXaxis()->SetTitle("#||{x_{true} - x_{reco}} [cm]");
    fXres_bot->GetXaxis()->SetTitle("#||{x_{true} - x_{reco}} [cm]");
    fXres_lat->GetXaxis()->SetTitle("#||{x_{true} - x_{reco}} [cm]");
    fXres_nor->GetXaxis()->SetTitle("#||{x_{true} - x_{reco}} [cm]");
    fXres_sth->GetXaxis()->SetTitle("#||{x_{true} - x_{reco}} [cm]");

    fYres_top->GetXaxis()->SetTitle("#||{y_{true} - y_{reco}} [cm]");
    fYres_rim->GetXaxis()->SetTitle("#||{y_{true} - y_{reco}} [cm]");
    fYres_bot->GetXaxis()->SetTitle("#||{y_{true} - y_{reco}} [cm]");
    fYres_lat->GetXaxis()->SetTitle("#||{y_{true} - y_{reco}} [cm]");
    fYres_nor->GetXaxis()->SetTitle("#||{y_{true} - y_{reco}} [cm]");
    fYres_sth->GetXaxis()->SetTitle("#||{y_{true} - y_{reco}} [cm]");

    fZres_top->GetXaxis()->SetTitle("#||{z_{true} - z_{reco}} [cm]");
    fZres_rim->GetXaxis()->SetTitle("#||{z_{true} - z_{reco}} [cm]");
    fZres_bot->GetXaxis()->SetTitle("#||{z_{true} - z_{reco}} [cm]");
    fZres_lat->GetXaxis()->SetTitle("#||{z_{true} - z_{reco}} [cm]");
    fZres_nor->GetXaxis()->SetTitle("#||{z_{true} - z_{reco}} [cm]");
    fZres_sth->GetXaxis()->SetTitle("#||{z_{true} - z_{reco}} [cm]");

    fTres_top->GetXaxis()->SetTitle("#||{t_{true} - t_{reco}} [ns]");
    fTres_rim->GetXaxis()->SetTitle("#||{t_{true} - t_{reco}} [ns]");
    fTres_bot->GetXaxis()->SetTitle("#||{t_{true} - t_{reco}} [ns]");
    fTres_lat->GetXaxis()->SetTitle("#||{t_{true} - t_{reco}} [ns]");
    fTres_nor->GetXaxis()->SetTitle("#||{t_{true} - t_{reco}} [ns]");
    fTres_sth->GetXaxis()->SetTitle("#||{t_{true} - t_{reco}} [ns]");


}

//main body
void CRTTruthMatchAnalysis::analyze(art::Event const& e)
{

    size_t nmu=0, nmu_dep=0, ntrue=0, ndata=0, nsim=0;
    size_t nmiss_true = 0, nmiss_data = 0, nmiss_sim = 0;

    //Intitialize CRTBackTracker
    //bt.Initialize(e);

    //MCParticles
    art::Handle< vector<simb::MCParticle> > mcHandle;
    map< int, const simb::MCParticle*> particleMap;
    if (!e.getByLabel(fSimulationLabel, mcHandle))
    {
        // If we have no MCParticles at all in an event, then we're in
        // big trouble. Throw an exception.
        throw cet::exception("CRTTruthMatchAnalysis")
          << " No simb::MCParticle objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    //fill particleMap (trackID->MCParticle object)
    for(auto const& particle : (*mcHandle)) {
        particleMap[particle.TrackId()] = &particle;
    }
    //mf::LogPrint("CRTTruthMatchAnalysis")
    //  << "found " << particleMap.size() << " MCParticles";

    //AuxDetSimChannels
    art::Handle< vector<sim::AuxDetSimChannel> > adscHandle;
    vector< art::Ptr<sim::AuxDetSimChannel> > adscList;
    if( e.getByLabel(fSimulationLabel,adscHandle) )
        art::fill_ptr_vector(adscList,adscHandle);
    else
        throw cet::exception("CRTTruthMatchAnalysis")
          << " No sim::AuxDetSimChannel objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;

    vector<int> muIds; //list of muon trackIds passing through CRT
    for(auto const& adsc : adscList) {
        for(auto const& ide : adsc->AuxDetIDEs()) {
            if(particleMap.find(ide.trackID)==particleMap.end())
                continue;
            if(abs(particleMap[ide.trackID]->PdgCode())!=13) //muons only
                continue;
            muIds.push_back(ide.trackID);
            nmu_dep++;
        }
    }
    //remove repeat IDs
    std::sort(muIds.begin(),muIds.end());
    muIds.erase(std::unique(muIds.begin(),muIds.end()),muIds.end());
    nmu = muIds.size();
    //mf::LogInfo("CRTTruthMatchAnalysis")
    //   << "found " << nmu << " muons passing through CRT"
    //   << " with " << nmu_dep << " associated energy deposits ("
    //   << 1.0*nmu_dep/nmu << ")";

    //CRTTrueHits
    art::Handle< vector<CRTHit> > trueHitHandle;
    vector< art::Ptr<CRTHit> > trueHitList;
    map<int,vector<art::Ptr<CRTHit>>> muToTrueHits; 

    if( e.getByLabel(fCRTTrueHitLabel,trueHitHandle) ) {
        art::fill_ptr_vector(trueHitList,trueHitHandle);

        for(auto const& hit : trueHitList) {
            for(const int id: bt.AllTrueIds(e,*hit)) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToTrueHits[id].push_back(hit);
            }//for trackIDs
        }//for CRTTrueHits

        ntrue = muToTrueHits.size();
    }//if CRTTrueHits
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTTrueHits found";

    //CRTData
    art::Handle< vector<CRTData> > dataHandle;
    vector< art::Ptr<CRTData> > dataList;
    map<int,vector<art::Ptr<CRTData>>> muToData;
    if( e.getByLabel(fCRTDataLabel,dataHandle)) {
        art::fill_ptr_vector(dataList,dataHandle);
        //mf::LogPrint("CRTTruthMatchAnalysis")
        //  << "found " << dataList.size() << " CRTData";

        for(auto const& data : dataList) {
            for(const int id: bt.AllTrueIds(e,*data)) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToData[id].push_back(data);
                ndata++;
            }//for trackIDs
        }//for data

        for(const int trk : muIds) {
            if(muToData.find(trk)==muToData.end()){
                nmiss_data++;
                fNDataMu->Fill(0);
            }
            else
                fNDataMu->Fill(muToData[trk].size());
        }

    }//if data
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTData found";

    //CRTSimHits
    art::Handle< vector<CRTHit> > simHitHandle;
    vector< art::Ptr<CRTHit> > simHitList;
    map<int,vector<art::Ptr<CRTHit>>> muToSimHits;
    if( e.getByLabel(fCRTSimHitLabel,simHitHandle) ) {
        art::fill_ptr_vector(simHitList,simHitHandle);

        for(auto const& hit : simHitList) {
            for(const int id: bt.AllTrueIds(e,*hit)) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToSimHits[id].push_back(hit);
                nsim++;
            }//for trackIDs
        }//for CRTSimHits

        if(fVerbose)
            mf::LogPrint("CRTTruthMatchAnalysis")
              << "  found " << nsim << " muon-associated CRTSimHits ("
              << 100.0*nsim/simHitList.size() << " %)";
    }//if CRTSimHits
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTSimHits found";

    //====tallying matches==========
    // loop over muon IDs
    size_t nmatch=0, nmiss=0, nmiss_track=0, nmiss_feb=0;
    size_t nmatch_top=0;
    for(const int trk : muIds) {
        //counters
        int ntruehit=0, nsimhit=0;
        int ntrue_top=0, ntrue_side=0, ntrue_bot=0;
        int nsim_top=0, nsim_side=0, nsim_bot=0;
        map<string,size_t> nregtrue, nregsim;

        bool simhasmu = true;
        //True Hits
        if(muToTrueHits.find(trk)!=muToTrueHits.end()){
            ntruehit = muToTrueHits[trk].size();
            for(auto const& hit : muToTrueHits[trk]) {
                //mac addresses
                for(auto const& id : hit->feb_id)
                    fFebsTrue->Fill((int)id);
                //region
                string truereg = hit->tagger;
                if(nregtrue.find(truereg)==nregtrue.end())
                    nregtrue[truereg] = 1;
                else 
                    nregtrue[truereg]++;
                //type
                char truetype = fCrtutils->GetRegTypeFromRegName(truereg);
                if(truetype=='c') ntrue_top++;
                if(truetype=='m') ntrue_side++;
                if(truetype=='d') ntrue_bot++;

                //Sim Hits
                if(muToSimHits.find(trk)!=muToSimHits.end()) {
                    for(auto const& simhit : muToSimHits[trk]) {
                        std::string simreg = simhit->tagger; //region
                        if(simreg != truereg) continue; //make same region
                        if(truetype == 'c' || truetype=='d'){  
                            if(simhit->feb_id[0] == hit->feb_id[0]){ //same mac5

                                if(simreg=="Top") {
                                    fRres_top->Fill(GetDeltaR(hit,simhit));
                                    fXres_top->Fill(hit->x_pos - simhit->x_pos);
                                    fYres_top->Fill(hit->y_pos - simhit->y_pos);
                                    fZres_top->Fill(hit->z_pos - simhit->z_pos);
                                    fTres_top->Fill(GetDeltaT(hit,simhit));
                                    nmatch++;
                                    /*if(nmatch_top<10&&fNViews<100) {
                                        //need way to map global back to local coords

                                        TVector3 truepos = 
                                          fCrtutils->WorldToModuleCoords(HitToPosVec(hit),fCrtutils->MacToAuxDetID(hit->feb_id[0],0));
                                        TVector3 recopos = 
                                          fCrtutils->WorldToModuleCoords(HitToPosVec(simhit),fCrtutils->MacToAuxDetID(hit->feb_id[0],0));
                                        fTrueHitPoints[fNViews]->SetPoint(0,truepos.X(),truepos.Z());
                                        fSimHitPoints[fNViews]->SetPoint(0,recopos.X(),recopos.Z());
                                        float pesx[16], pesy[16];
                                        for(size_t i=0; i<16; i++){pesx[i]=0.; pesy[i]=0.;}
                                        if(simhit->pesmap.find(simhit->feb_id[0])==simhit->pesmap.end())
                                            std::cout << "FEB not found in PEs map!" << std::endl;
                                        for(auto const& chanpe : (*((simhit->pesmap).find(simhit->feb_id[0]))).second) {
                                            if(chanpe.first<16) pesx[chanpe.first] = chanpe.second;
                                            else pesy[chanpe.first-16] = chanpe.second;
                                        }
                                        for(size_t i=0; i<16; i++){
                                            for(size_t j=0; j<16; j++){
                                                fChargeViews[fNViews]->SetBinContent(i,j,pesx[i]+pesy[j]);
                                            }
                                        }
                                        fNViews++;
                                    }*/
                                    nmatch_top++;     
                                }
                                if(simreg.find("Rim")!=string::npos){
                                    fRres_rim->Fill(GetDeltaR(hit,simhit));
                                    fXres_rim->Fill(hit->x_pos - simhit->x_pos);
                                    fYres_rim->Fill(hit->y_pos - simhit->y_pos);
                                    fZres_rim->Fill(hit->z_pos - simhit->z_pos);
                                    fTres_rim->Fill(GetDeltaT(hit,simhit));
                                    nmatch++;
                                }
                                if(simreg=="Bottom") {
                                    fRres_bot->Fill(GetDeltaR(hit,simhit));
                                    fXres_bot->Fill(hit->x_pos - simhit->x_pos);
                                    fYres_bot->Fill(hit->y_pos - simhit->y_pos);
                                    fZres_bot->Fill(hit->z_pos - simhit->z_pos);
                                    fTres_bot->Fill(GetDeltaT(hit,simhit));
                                    nmatch++;
                                }
                            }
                            else {
                                nmiss++;
                                nmiss_feb++;
                            }//no feb match
                        }//c or d type
                        else {//m type
                            bool match = false;
                            for(auto const& truefeb : hit->feb_id){
                                for(auto const& simfeb : simhit->feb_id){
                                    if(truefeb==simfeb) {
                                        match = true;
                                        break;
                                    }
                                }
                                if(match) break;
                            }
                            if(match){
                                if(simreg=="North") {
                                    fRres_nor->Fill(GetDeltaR(hit,simhit));
                                    fXres_nor->Fill(hit->x_pos - simhit->x_pos);
                                    fYres_nor->Fill(hit->y_pos - simhit->y_pos);
                                    fZres_nor->Fill(hit->z_pos - simhit->z_pos);
                                    fTres_nor->Fill(GetDeltaT(hit,simhit));
                                }
                                else if(simreg=="South") {
                                    fRres_sth->Fill(GetDeltaR(hit,simhit));
                                    fXres_sth->Fill(hit->x_pos - simhit->x_pos);
                                    fYres_sth->Fill(hit->y_pos - simhit->y_pos);
                                    fZres_sth->Fill(hit->z_pos - simhit->z_pos);
                                    fTres_sth->Fill(GetDeltaT(hit,simhit));
                                }
                                else {
                                    fRres_lat->Fill(GetDeltaR(hit,simhit));
                                    fXres_lat->Fill(hit->x_pos - simhit->x_pos);
                                    fYres_lat->Fill(hit->y_pos - simhit->y_pos);
                                    fZres_lat->Fill(hit->z_pos - simhit->z_pos);
                                    fTres_lat->Fill(GetDeltaT(hit,simhit));
                                }
                            }
                            else nmiss++;
                        }//else m type
                    }//matching sim hits
                } //if muID trk found
                else if(simhasmu){
                    simhasmu = false;
                    nmiss_track++;
                    nmiss++;
                }
            }//loop over true hits
        }//if mu trk found in true map
        else
            nmiss_true++;

        //Sim Hits
        if(muToSimHits.find(trk)!=muToSimHits.end()) {
            nsimhit = muToSimHits[trk].size();
            for(auto const& hit : muToSimHits[trk]) {
                for(auto const& id : hit->feb_id)
                    fFebsSim->Fill((int)id);

                if(nregsim.find(hit->tagger)==nregsim.end())
                    nregsim[hit->tagger] = 1;
                else
                    nregsim[hit->tagger]++;

                char type = fCrtutils->GetRegTypeFromRegName(hit->tagger);
                if(type=='c') nsim_top++;
                if(type=='m') nsim_side++;
                if(type=='d') nsim_bot++;
            }
        }
        else
            nmiss_sim++;

        //Fill histos
        fNTrueHitMu->Fill(ntruehit);
        fNSimHitMu->Fill(nsimhit);
        fNTrueSimDiffMu->Fill(ntruehit-nsimhit);

        if(ntrue_top>0)  fNTrueHitMu_top->Fill(ntrue_top);
        if(ntrue_side>0) fNTrueHitMu_side->Fill(ntrue_side);
        if(ntrue_bot>0)  fNTrueHitMu_bot->Fill(ntrue_bot);
        if(nsim_top>0)   fNSimHitMu_top->Fill(nsim_top);
        if(nsim_side>0)  fNSimHitMu_side->Fill(nsim_side);
        if(nsim_bot>0)   fNSimHitMu_bot->Fill(nsim_bot);

        for(auto const& nreg : nregtrue) fNRegTrue->Fill(nreg.second);
        for(auto const& nreg : nregsim) fNRegSim->Fill(nreg.second);

    }//for muon IDs

    if(fVerbose) {
        std::cout << "found " << '\n'
                  << " * nmu: " << nmu << '\n'
                  << " * ntrue hit: " << trueHitList.size() << " with " << ntrue << " ass'd with muons" << '\n'
                  << " * nsim hit: " << simHitList.size() << " with "  << nsim << " ass'd with muons" << '\n' 
                  << std::endl;

        std::cout << "matched true/sim hits:" << '\n'
                  << " * matched: " << nmatch << '\n'
                  << " * missed: "  << nmiss  << '\n'
                  << "   - SimHits missing muon track: " << nmiss_track << '\n'
                  << "   - SimHits, c or d missing FEB match: " << nmiss_feb << '\n' 
                  << std::endl;

        std::cout << "missed muons (eff.):" << '\n'
                  << " * trueHit: " << nmiss_true << " (" << 1.0*(nmu-nmiss_true)/nmu << ")" << '\n'
                  << " * data: " << nmiss_data << " (" << 1.0*(nmu-nmiss_data)/nmu << ")" << '\n'
                  << " * simHit: " << nmiss_sim << " (" << 1.0*(nmu-nmiss_sim)/nmu << ")" << '\n' 
                  << std::endl;
    }

}

float CRTTruthMatchAnalysis::GetDeltaR(art::Ptr<CRTHit> h1, art::Ptr<CRTHit> h2){
    float dr;
    double dx = h1->x_pos - h2->x_pos;
    double dy = h1->y_pos - h2->y_pos;
    double dz = h1->z_pos - h2->z_pos;

    dr = sqrt(dx*dx + dy*dy + dz*dz);
    return dr;
}


int CRTTruthMatchAnalysis::GetDeltaT(art::Ptr<CRTHit> h1, art::Ptr<CRTHit> h2){
    int t1 = h1->ts0_ns;
    int t2 = h2->ts0_ns;
    return t1-t2; 
}

TVector3 CRTTruthMatchAnalysis::HitToPosVec(art::Ptr<CRTHit> hit){

    TVector3 vec(hit->x_pos,hit->y_pos,hit->z_pos);
    return vec;
}

DEFINE_ART_MODULE(CRTTruthMatchAnalysis)
