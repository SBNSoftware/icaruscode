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

    art::InputTag fSimulationLabel;
    art::InputTag fCRTTrueHitLabel;
    art::InputTag fCRTDataLabel;
    art::InputTag fCRTSimHitLabel;
    CRTBackTracker bt;

    TH1F* fNTrueHitMu; //number of CRTTrueHits per muon track
    TH1F* fNDataMu;    //number of CRTData per muon track
    TH1F* fNSimHitMu;   //number of CRTSimHits per muon track

};//class definition

//constructor
CRTTruthMatchAnalysis::CRTTruthMatchAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSimulationLabel(p.get<art::InputTag>("SimulationLabel","largeant")),
    fCRTTrueHitLabel(p.get<art::InputTag>("CRTTrueHitLabel","crttruehit")),
    fCRTDataLabel(p.get<art::InputTag>("CRTDataLabel","crtdaq")),
    fCRTSimHitLabel(p.get<art::InputTag>("CRTSimHitLabel","crthit")),
    bt(p.get<fhicl::ParameterSet>("CRTBackTrack"))
{}

void CRTTruthMatchAnalysis::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    fNTrueHitMu = tfs->make<TH1F>("NTrueHitMu","Number of True Hits per #mu Track",5,0,5);
    fNDataMu    = tfs->make<TH1F>("NDataMu", "Number of FEB Triggers per #mu Track",10,0,10);
    fNSimHitMu  = tfs->make<TH1F>("NSimHitMu","Number of Sim Hits per #mu Track",5,0,5);
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
    mf::LogPrint("CRTTruthMatchAnalysis")
      << "found " << particleMap.size() << " MCParticles";

    std::cout << "loop over AuxDetChannels" << std::endl;
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
    map<int,bool> trueRecoedTracks;
    for(auto const& adsc : adscList) {
        for(auto const& ide : adsc->AuxDetIDEs()) {
            if(particleMap.find(ide.trackID)==particleMap.end())
                continue;
            if(abs(particleMap[ide.trackID]->PdgCode())!=13) //muons only
                continue;
            muIds.push_back(ide.trackID);
            trueRecoedTracks[ide.trackID] = false;
            nmu_dep++;
        }
    }
    //remove repeat IDs
    std::sort(muIds.begin(),muIds.end());
    muIds.erase(std::unique(muIds.begin(),muIds.end()),muIds.end());
    nmu = trueRecoedTracks.size();
    mf::LogPrint("CRTTruthMatchAnalysis")
       << "found " << nmu << " muons passing through CRT"
       << " with " << nmu_dep << " associated energy deposits ("
       << 1.0*nmu_dep/nmu << ")";

    //CRTTrueHits
    art::Handle< vector<CRTHit> > trueHitHandle;
    vector< art::Ptr<CRTHit> > trueHitList;
    map<int,bool> simRecoedTracks;
    map<int,bool> dataTrueTracks;
    map<int,vector<art::Ptr<CRTHit>>> muToTrueHits; 

    if( e.getByLabel(fCRTTrueHitLabel,trueHitHandle) ) {
        art::fill_ptr_vector(trueHitList,trueHitHandle);
        mf::LogPrint("CRTTruthMatchAnalysis") 
          << "found " << trueHitList.size() << " CRTTrueHits";

        for(auto const& hit : trueHitList) {
            vector<int> btIds = bt.AllTrueIds(e,*hit);
            for(const int id: btIds) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToTrueHits[id].push_back(hit);
                trueRecoedTracks[id] = true;
                dataTrueTracks[id] = false;
                simRecoedTracks[id] = false;
                //ntrue++;
            }//for trackIDs
        }//for CRTTrueHits
        for(const int trk : muIds){
            if(muToTrueHits.find(trk)==muToTrueHits.end()){
                nmiss_true++;
                fNTrueHitMu->Fill(0);
            }
            else
                fNTrueHitMu->Fill(muToTrueHits[trk].size());
        }

        ntrue = simRecoedTracks.size();
        mf::LogPrint("CRTTruthMatchAnalysis")
          << "  found " << ntrue << " muon-associated CRTTrueHits ("
          << 100.0*ntrue/trueHitList.size() << " %)";
    }//if CRTTrueHits
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTTrueHits found";

    //CRTData
    art::Handle< vector<CRTData> > dataHandle;
    vector< art::Ptr<CRTData> > dataList;
    map<int,vector<art::Ptr<CRTData>>> muToData;
    map<int,bool> dataSimTracks;
    if( e.getByLabel(fCRTDataLabel,dataHandle)) {
        art::fill_ptr_vector(dataList,dataHandle);
        mf::LogPrint("CRTTruthMatchAnalysis")
          << "found " << dataList.size() << " CRTData";

        for(auto const& data : dataList) {
            vector<int> btIds = bt.AllTrueIds(e,*data);
            for(const int id: btIds) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToData[id].push_back(data);
                dataTrueTracks[id] = true;
                dataSimTracks[id] = false;
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

        mf::LogPrint("CRTTruthMatchAnalysis")
          << "  found " << ndata << " muon-associated CRTData ("
          << 100.0*ndata/dataList.size() << " %)";
    }//if data
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTData found";

    //CRTSimHits
    art::Handle< vector<CRTHit> > simHitHandle;
    vector< art::Ptr<CRTHit> > simHitList;
    map<int,vector<art::Ptr<CRTHit>>> muToSimHit;
    if( e.getByLabel(fCRTSimHitLabel,simHitHandle) ) {
        art::fill_ptr_vector(simHitList,simHitHandle);
        mf::LogPrint("CRTTruthMatchAnalysis")
          << "found " << simHitList.size() << " CRTSimHits";

        for(auto const& hit : simHitList) {
            vector<int> btIds = bt.AllTrueIds(e,*hit);
            for(const int id: btIds) {
                if(particleMap.find(id)==particleMap.end())
                    continue;
                if(abs(particleMap[id]->PdgCode())!=13) //muons only
                    continue;
                muToSimHit[id].push_back(hit);
                simRecoedTracks[id] = true;
                dataSimTracks[id] = true;
                nsim++;
            }//for trackIDs
        }//for CRTSimHits

        for(const int trk : muIds ){
            if(muToSimHit.find(trk)==muToSimHit.end()){
                nmiss_sim++;
                fNSimHitMu->Fill(0);
            }
            else
                fNSimHitMu->Fill(muToSimHit[trk].size());
        }

        mf::LogPrint("CRTTruthMatchAnalysis")
          << "  found " << nsim << " muon-associated CRTSimHits ("
          << 100.0*nsim/simHitList.size() << " %)";
    }//if CRTSimHits
    else
        mf::LogWarning("CRTTruthMatchAnalysis") << "no CRTSimHits found";

    //tallying matches
    std::cout << "missed muons (eff.):" << std::endl;
    std::cout << "  trueHit: " << nmiss_true << " (" << (1.0-nmiss_true)/nmu << "), data: "
              << nmiss_data << " (" << (1.0-nmiss_data)/nmu << "), simHit: " 
              << nmiss_sim << " (" << (1.0-nmiss_sim)/nmu << ")" << std::endl;

    std::cout << "matching trackIDs" << std::endl;
    size_t nmiss_truetrue=0, nmiss_truedat=0, nmiss_truesim=0, nmiss_datsim=0;
    for(auto const& trk : trueRecoedTracks) {
        if(trk.second==false) nmiss_truetrue++;
    }
    for(auto const& trk : dataTrueTracks) {
        if(trk.second==false) nmiss_truedat++;
    }
    for(auto const& trk : dataSimTracks) {
        if(trk.second==false) nmiss_datsim++;
    }
    for(auto const& trk : simRecoedTracks) {
        if(trk.second==false) nmiss_truesim++;
    }


    std::cout << ntrue-nmiss_truetrue << " (" << 100.0*(ntrue-nmiss_truetrue)/ntrue
              << " %) muon tracks matched to true hits" << std::endl;
    std::cout << ndata-nmiss_truedat << " (" << 100.0*(ndata-nmiss_truedat)/ndata
              << " %) FEB triggers matched to true hits" << std::endl;
    std::cout << nsim-nmiss_datsim << " (" << 100.0*(nsim-nmiss_datsim)/nsim
              << " %) sim hits matched to FEB triggers" << std::endl;
    std::cout << nsim-nmiss_truesim << " (" << 100.0*(nsim-nmiss_truesim)/nsim
              << " %) sim hits matched to true hits" << std::endl;

}

DEFINE_ART_MODULE(CRTTruthMatchAnalysis)
