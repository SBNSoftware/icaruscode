////////////////////////////////////////////////////////////////////////
// Class:       PhotBackground
// Plugin Type: analyzer (art v3_05_01)
// File:        PhotBackground_module.cc
//
// Generated at Tue Jul 14 12:07:00 2020 by Christopher Hilgenberg using cetskelgen
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

//LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

//local includes
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

//C++ includes
#include <vector>
#include <string>
#include <map>

//ROOT includes
#include <TTree.h>

namespace icarus {
  class PhotBackground;
}

using namespace icarus;
using std::vector;
using std::string;
using std::map;

class icarus::PhotBackground : public art::EDAnalyzer {
public:
  explicit PhotBackground(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotBackground(PhotBackground const&) = delete;
  PhotBackground(PhotBackground&&) = delete;
  PhotBackground& operator=(PhotBackground const&) = delete;
  PhotBackground& operator=(PhotBackground&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

    int EnteredLar(const TLorentzVector& v, bool& iv, bool& fv);

    crt::CRTBackTracker bt;
    crt::CRTCommonUtils* fCrtutils;

    geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();   ///< pointer to Geometry provider
    geo::CryostatGeo const& cryo0 = fGeometryService->Cryostat(0);
    geo::CryostatGeo const& cryo1 = fGeometryService->Cryostat(1);

    geo::TPCGeo const& tpc00 = cryo0.TPC(0);
    geo::TPCGeo const& tpc01 = cryo0.TPC(1);
    geo::TPCGeo const& tpc10 = cryo1.TPC(0);
    geo::TPCGeo const& tpc11 = cryo1.TPC(1);


    TTree*      fTree;
    int         fIcode;         ///< end process code [-1->do not use; 0-> pair prod; 1-> compton]
    double      fStartPos[4];   ///< photon start x,y,z,t [cm/ns]
    double      fEndPos[4];     ///< photon end x,y,z,t [cm/ns]
    double      fStartE;        ///< photon start energy [GeV]
    double      fEndE;          ///< photon end energy [GeV]
    int         fStartReg;      ///< region where photon was produced
    int         fEndReg;        ///< region where photon interacted
    bool        fPhotIV;
    bool        fPhotAV;
    vector<int> fMuRegs;
    bool        fMuIV;
    bool        fMuAV;
    bool        fMuTag;
    int         fMuNHit;
    vector<vector<double>> fMuHitPos;
    vector<double> fMuHitDx;
    vector<double> fMuHitDt;
};


PhotBackground::PhotBackground(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    bt(p.get<fhicl::ParameterSet>("CRTBackTrack")),
    fCrtutils(new crt::CRTCommonUtils())
{
    art::ServiceHandle<art::TFileService> tfs;
  
    fTree = tfs->make<TTree>("anatree","analysis tree for photon backgrounds");
    fTree->Branch("icode",    &fIcode,   "icode/I");
    fTree->Branch("startpos", fStartPos, "startpos[4]/D");
    fTree->Branch("endpos",   fEndPos,   "endpos[4]/D");
    fTree->Branch("starte",   &fStartE,  "starte/D");
    fTree->Branch("ende",     &fEndE,    "ende/D");
    fTree->Branch("startreg", &fStartReg,"startreg/I");
    fTree->Branch("endreg",   &fEndReg,  "endreg/I");
    fTree->Branch("photav",   &fPhotAV,    "photav/O");
    fTree->Branch("photiv",   &fPhotIV,    "photiv/O");
    fTree->Branch("regs",     &fMuRegs);
    fTree->Branch("muav",     &fMuAV,    "muav/O");
    fTree->Branch("muiv",     &fMuIV,    "muiv/O");
    fTree->Branch("mutag",    &fMuTag,   "mutag/O");
    fTree->Branch("nhit",     &fMuNHit,  "nhit/I");
    fTree->Branch("hitpos",   &fMuHitPos);
    fTree->Branch("dx",       &fMuHitDx);
    fTree->Branch("dt",       &fMuHitDt);
}

void PhotBackground::analyze(art::Event const& e)
{
    // Define a "handle" to point to a vector of MCParticle objects.
    art::Handle< vector<simb::MCParticle> > particleHandle;
    e.getByLabel("largeant", particleHandle);
    map<int,const simb::MCParticle*> idToMu;
    vector<const simb::MCParticle*> photList;

    //true CRT hits
    art::Handle< vector<sbn::crt::CRTHit> > trueHitHandle;
    vector< art::Ptr<sbn::crt::CRTHit> > trueHitList;
    map<int,vector<art::Ptr<sbn::crt::CRTHit>>> muToTrueHits;

    if( e.getByLabel("crttruehit",trueHitHandle) ) 
        art::fill_ptr_vector(trueHitList,trueHitHandle);

    for(auto const& particle : *particleHandle) {

        if(abs(particle.PdgCode())==13){
            idToMu[particle.TrackId()] = &particle;
        }
        
        else if(particle.PdgCode()==22 && particle.Process()=="muBrems")
            photList.push_back(&particle);

        else 
            continue;
    }//for MCParticles

    for(auto const& hit : trueHitList){
        for(const int id: bt.AllTrueIds(e,*hit)) {
            if(idToMu.find(id)==idToMu.end())
                continue;
            muToTrueHits[id].push_back(hit);
        }//for trackIDs
    }//for CRTTrueHits

    for(auto const& phot : photList) {

        const TLorentzVector& positionStart = phot->Position();
        const TLorentzVector& positionEnd   = phot->EndPosition();

        fPhotIV = false;
        fPhotAV = false;
        fEndReg = EnteredLar(positionEnd, fPhotIV, fPhotAV);
        if(fEndReg==-1)
            continue;
        bool iv = false, av = false;
        fStartReg = EnteredLar(positionStart, iv, av);

        for(size_t i=0; i<4; i++) {
            fStartPos[i] = positionStart[i];
            fEndPos[i]   = positionEnd[i];
        }
        fStartE = phot->E();
        fEndE   = phot->E(phot->NumberTrajectoryPoints()-2); 
        
        fIcode = -1;
        string endprocess = phot->EndProcess();
        if(endprocess=="conv") fIcode = 0; //pair prod
        else if(endprocess=="compt") fIcode = 1; //compton scatter
        else if(endprocess=="phot") fIcode = 2;  //photoelectric
        else std::cout << "uknown end process: " << endprocess << std::endl;    



        if(idToMu.find(phot->Mother())==idToMu.end())
            throw cet::exception("PhotBackground") <<
                  "no match muon found for mubrems photons" << std::endl;

        fMuIV = false;
        fMuAV = false;
        fMuTag = false;
        fMuNHit = 0;
        fMuRegs.clear();
        fMuHitPos.clear();
        fMuHitDx.clear();
        fMuHitDt.clear();

        if(muToTrueHits.find(phot->Mother())!=muToTrueHits.end()){
            fMuTag = true;
            fMuNHit=muToTrueHits[phot->Mother()].size();
            for(auto const& hit : muToTrueHits[phot->Mother()]){
                vector<double> xyzt = {hit->x_pos,hit->y_pos,hit->z_pos,(double)hit->ts0_ns};
                fMuHitPos.push_back(xyzt);
                double l = 0.;
                for(int j=0; j<3; j++) l+=pow(xyzt[j]-fEndPos[j],2);
                fMuHitDx.push_back(sqrt(l));
                fMuHitDt.push_back((double)hit->ts0_ns-1.6e6-fEndPos[3]);
                string reg = hit->tagger;
                //fMuRegType.push_back(fCrtutils->GetRegTypeFromRegName(reg));
                fMuRegs.push_back(fCrtutils->AuxDetRegionNameToNum(reg));
            }
        }

        for(size_t i=0; i<idToMu[phot->Mother()]->NumberTrajectoryPoints(); i++){

            const TLorentzVector& position = idToMu[phot->Mother()]->Position(i);
            fMuRegs.push_back(EnteredLar(position, fMuIV, fMuAV));

        }//for traj pts

        std::sort(fMuRegs.begin(),fMuRegs.end());
        fMuRegs.resize(std::distance(fMuRegs.begin(),std::unique(fMuRegs.begin(),fMuRegs.end())));

        fTree->Fill();

    }//for photons

}//end analyze

int PhotBackground::EnteredLar(const TLorentzVector& v, bool& iv, bool& av){

    double pos[3] = {v.X(),v.Y(),v.Z()};
    int reg = -1;

    if(cryo0.ContainsPosition(pos)){

        if(tpc00.ContainsPosition(pos)) {
            reg = 5;
            av = true;
            iv = false;
        }
        else if(tpc01.ContainsPosition(pos)) {
            reg = 6;
            av = true;
            iv = false;
        }
        else {
            reg = 10;
            if(!av)iv = true;
        }
    }
    else if(cryo1.ContainsPosition(pos)){

        if(tpc10.ContainsPosition(pos)){
            reg = 7;
            av = true;
            iv = false;
        }
        else if(tpc11.ContainsPosition(pos)) {
            reg = 8;
            av = true;
            iv = false;
        }
        else {
            reg = 12;
            if(!av)iv = true;
        }
    }

    return reg;

}

DEFINE_ART_MODULE(PhotBackground)
