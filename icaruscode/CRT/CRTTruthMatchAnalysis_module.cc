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
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

//C++ includes
#include <vector>
#include <map>

//ROOT includes
#include <TH1F.h>

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

    art::InputTag fSimulationLabel;
    art::InputTag fCRTTrueHitLabel;
    art::InputTag fCRTDataLabel;
    art::InputTag fCRTSimHitLabel;
    CRTBackTracker bt;
    CRTCommonUtils* fCrtutils;

    TH1F* fNTrueHitMu; //number of CRTTrueHits per muon track
    TH1F* fNDataMu;    //number of CRTData per muon track
    TH1F* fNSimHitMu;   //number of CRTSimHits per muon track
    TH1F* fNTrueSimDiffMu; //difference in number of true - sim hits

    TH1F* fNTrueHitMu_top;
    TH1F* fNTrueHitMu_side;
    TH1F* fNTrueHitMu_bot;
    TH1F* fNDataMu_top;
    TH1F* fNDataMu_side;
    TH1F* fNDataMu_bot;
    TH1F* fNSimHitMu_top;
    TH1F* fNSimHitMu_side;
    TH1F* fNSimHitMu_bot;

    TH1F* fNRegTrue;
    TH1F* fNRegSim;

    TH1F* fFebsTrue;
    TH1F* fFebsSim;

    TH1F* fRres_top;
    TH1F* fRres_rim;
    TH1F* fRres_lat;
    TH1F* fRres_nor;
    TH1F* fRres_sth;
    TH1F* fRres_bot;

};//class definition

//constructor
CRTTruthMatchAnalysis::CRTTruthMatchAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSimulationLabel(p.get<art::InputTag>("SimulationLabel","largeant")),
    fCRTTrueHitLabel(p.get<art::InputTag>("CRTTrueHitLabel","crttruehit")),
    fCRTDataLabel(p.get<art::InputTag>("CRTDataLabel","crtdaq")),
    fCRTSimHitLabel(p.get<art::InputTag>("CRTSimHitLabel","crthit")),
    bt(p.get<fhicl::ParameterSet>("CRTBackTrack")),
    fCrtutils(new CRTCommonUtils())
{}

void CRTTruthMatchAnalysis::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    fNTrueHitMu      = tfs->make<TH1F>("hNTrueHitMu",     "Number of True Hits per #mu Track",    10,0,10);
    fNDataMu         = tfs->make<TH1F>("hNDataMu",        "Number of FEB Triggers per #mu Track", 10,0,10);
    fNSimHitMu       = tfs->make<TH1F>("hNSimHitMu",      "Number of Sim Hits per #mu Track",     10,0,10);
    fNTrueSimDiffMu  = tfs->make<TH1F>("hNTrueSimDiffMu", "Number of True-Sim Hits per #mu Track",14,-6,8); 

    fNTrueHitMu_top  = tfs->make<TH1F>("hNTrueHitMu_top", "Number of True Hits per #mu Track",  10,0,10);
    fNTrueHitMu_side = tfs->make<TH1F>("hNTrueHitMu_side","Number of True Hits per #mu Track",  10,0,10);
    fNTrueHitMu_bot  = tfs->make<TH1F>("hNTrueHitMu_bot", "Number of True Hits per #mu Track",  10,0,10);
    fNSimHitMu_top   = tfs->make<TH1F>("hNSimHitMu_top",  "Number of Sim Hits per #mu Track",   10,0,10);
    fNSimHitMu_side  = tfs->make<TH1F>("hNSimHitMu_side", "Number of Sim Hits per #mu Track",   10,0,10);
    fNSimHitMu_bot   = tfs->make<TH1F>("hNSimHitMu_bot",  "Number of Sim Hits per #mu Track",   10,0,10);

    fNRegTrue        = tfs->make<TH1F>("hNRegTrue",       "Number of True Hits per Region",     10,0,10);
    fNRegSim         = tfs->make<TH1F>("hNRegSim",        "Number of Sim Hits per Region",      10,0,10); 

    fFebsTrue        = tfs->make<TH1F>("hFebsTrue",       "mac5s in true hits",               256,0,256);
    fFebsSim         = tfs->make<TH1F>("hFebsSim",        "mac5s in sim hits",                256,0,256);

    fRres_top        = tfs->make<TH1F>("hRres_top",       "CRTHit Spatial Resolution: Top",     60,0,300);
    fRres_rim        = tfs->make<TH1F>("hRres_rim",       "CRTHit Spatial Resolution: Rim",     60,0,300);
    fRres_bot        = tfs->make<TH1F>("hRres_bot",       "CRTHit Spatial Resolution: Bottom",  60,0,300);
    fRres_lat        = tfs->make<TH1F>("hRres_lat",       "CRTHit Spatial Resolution: East/West",  60,0,500);
    fRres_nor        = tfs->make<TH1F>("hRres_nor",       "CRTHit Spatial Resolution: North",  60,0,500);
    fRres_sth        = tfs->make<TH1F>("hRres_sth",       "CRTHit Spatial Resolution: South",  60,0,500);

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

}

//main body
void CRTTruthMatchAnalysis::analyze(art::Event const& e)
{

    size_t nmu=0, nmu_dep=0, ntrue=0, ndata=0, nsim=0;
    size_t nmiss_true = 0, nmiss_data = 0, nmiss_sim = 0;

    //Intitialize CRTBackTracker
    bt.Initialize(e);

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
        //mf::LogPrint("CRTTruthMatchAnalysis") 
        //  << "found " << trueHitList.size() << " CRTTrueHits";

        for(auto const& hit : trueHitList) {
            vector<int> btIds = bt.AllTrueIds(e,*hit);
            for(const int id: btIds) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToTrueHits[id].push_back(hit);
            }//for trackIDs
        }//for CRTTrueHits

        ntrue = muToTrueHits.size();
        mf::LogInfo("CRTTruthMatchAnalysis")
          << "  found " << ntrue << " muon-associated CRTTrueHits ("
          << 100.0*ntrue/trueHitList.size() << " %)";
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
            vector<int> btIds = bt.AllTrueIds(e,*data);
            for(const int id: btIds) {
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

        //mf::LogInfo("CRTTruthMatchAnalysis")
        //  << "  found " << ndata << " muon-associated CRTData ("
        //  << 100.0*ndata/dataList.size() << " %)";
    }//if data
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTData found";

    //CRTSimHits
    art::Handle< vector<CRTHit> > simHitHandle;
    vector< art::Ptr<CRTHit> > simHitList;
    map<int,vector<art::Ptr<CRTHit>>> muToSimHits;
    if( e.getByLabel(fCRTSimHitLabel,simHitHandle) ) {
        art::fill_ptr_vector(simHitList,simHitHandle);
        //mf::LogInfo("CRTTruthMatchAnalysis")
        //  << "found " << simHitList.size() << " CRTSimHits";

        for(auto const& hit : simHitList) {
            vector<int> btIds = bt.AllTrueIds(e,*hit);
            for(const int id: btIds) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToSimHits[id].push_back(hit);
                nsim++;
            }//for trackIDs
        }//for CRTSimHits

        mf::LogPrint("CRTTruthMatchAnalysis")
          << "  found " << nsim << " muon-associated CRTSimHits ("
          << 100.0*nsim/simHitList.size() << " %)";
    }//if CRTSimHits
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTSimHits found";

    //====tallying matches==========
    // loop over muon IDs
    size_t nmatch=0, nmiss=0, nmiss_track=0;
    for(const int trk : muIds) {
        //counters
        int ntruehit=0, nsimhit=0;
        int ntrue_top=0, ntrue_side=0, ntrue_bot=0;
        int nsim_top=0, nsim_side=0, nsim_bot=0;
        map<string,size_t> nregtrue, nregsim;

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
                        string simreg = simhit->tagger; //region
                        if(simreg != truereg) continue; //make same region
                        if(truetype == 'c' || truetype=='d'){  
                            if(1) {//simhit->feb_id[0] == hit->feb_id[0]){ //same mac5

                                if(simreg=="Top") {
                                    fRres_top->Fill(GetDeltaR(hit,simhit));
                                    nmatch++;
                                }
                                if(simreg.find("Rim")!=string::npos){
                                    fRres_rim->Fill(GetDeltaR(hit,simhit));
                                    nmatch++;
                                }
                                if(simreg=="Bottom") {
                                    fRres_bot->Fill(GetDeltaR(hit,simhit));
                                    nmatch++;
                                }
                            }
                            else nmiss++;
                        }
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
                                if(simreg=="North") fRres_nor->Fill(GetDeltaR(hit,simhit));
                                else if(simreg=="South") fRres_sth->Fill(GetDeltaR(hit,simhit));
                                else fRres_lat->Fill(GetDeltaR(hit,simhit));
                            }
                            else nmiss++;
                        }//else m type
                    }//matching sim hits
                } //if muID trk found
                else {
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

    std::cout << "matched true/sim hits:" << std::endl;
    std::cout << "   matched: " << nmatch << std::endl;
    std::cout << "   missed: "  << nmiss  << ", missing SimHits for muon track: " << nmiss_track << std::endl;

    std::cout << "missed muons (eff.) out of " << nmu <<":" << std::endl;
    std::cout << "  trueHit: " << nmiss_true << " (" << 1.0*(nmu-nmiss_true)/nmu << "), data: "
              << nmiss_data << " (" << 1.0*(nmu-nmiss_data)/nmu << "), simHit: " 
              << nmiss_sim << " (" << 1.0*(nmu-nmiss_sim)/nmu << ")" << std::endl;


}

float CRTTruthMatchAnalysis::GetDeltaR(art::Ptr<CRTHit> h1, art::Ptr<CRTHit> h2){
    float dr;
    double dx = h1->x_pos - h2->x_pos;
    double dy = h1->y_pos - h2->y_pos;
    double dz = h1->z_pos - h2->z_pos;

    dr = sqrt(dx*dx + dy*dy + dz*dz);
    return dr;
}

DEFINE_ART_MODULE(CRTTruthMatchAnalysis)
