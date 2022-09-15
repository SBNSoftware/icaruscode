////////////////////////////////////////////////////////////////////////
// Class:       CRTTPCMatchingAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTTPCMatchingAna_module.cc
//
// Generated at Wed Feb 16 22:07:45 2022 by Biswaranjan Behera using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"


// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// icaruscode includes
#include "icaruscode/CRT/CRTUtils/RecoUtils.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"
//#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
//#include "icaruscode/CRT/CRTUtils/CRTEventDisplay.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"


#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace icarus {
  class CRTTPCMatchingAna;
}


class icarus::CRTTPCMatchingAna : public art::EDAnalyzer {
public:
  explicit CRTTPCMatchingAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTPCMatchingAna(CRTTPCMatchingAna const&) = delete;
  CRTTPCMatchingAna(CRTTPCMatchingAna&&) = delete;
  CRTTPCMatchingAna& operator=(CRTTPCMatchingAna const&) = delete;
  CRTTPCMatchingAna& operator=(CRTTPCMatchingAna&&) = delete;

  // Called once, at start of the job
  virtual void beginJob() override;

  // Called once per event
  virtual void analyze(const art::Event& event) override;

  // Called once, at end of the job
  virtual void endJob() override;

  // Calculate the distance from the track crossing point to CRT overlap coordinates
  double DistToCrtHit(TVector3 trackPos, sbn::crt::CRTHit crtHit);

  void reconfigure(fhicl::ParameterSet const & p);
  // Required functions.
  //  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag              fCRTHitLabel;   ///< name of CRT producer
  std::vector<art::InputTag> fTPCTrackLabel; ///< labels for source of tracks
  std::vector<art::InputTag> fPFParticleLabel; ///< labels for source of PFParticle
  art::InputTag              fTriggerLabel;  ///< labels for trigger 
  bool                       fVerbose;       ///< print information about what's going on

  CRTT0MatchAlg t0Alg;
  //CRTEventDisplay evd;
  geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
  icarus::crt::CRTCommonUtils* fCrtutils;
  //  CRTCommonUtils* fCrtutils;

  TTree*                 fTree;
  int                    fEvent;        ///< number of the event being processed
  int                    fRun;          ///< number of the run being processed
  int                    fSubRun;       ///< number of the sub-run being processed
  vector<int>            fCrtRegion;    ///< CRT hit region code
  vector<double>         fDCA;          ///< Distance of closest approach between CRT hit and extended TPC track
  vector<double>         fDOL;          ///<
  vector<double>         fT0;           ///< Track T0 based on CRT hit and extended TPC Track match 
  vector<double>         fpandorat0;    ///< Track T0 based on Pandora (Cathode Crossing Track)
  vector<int>            fcryo;         ///< cryo number

  //add trigger data product vars
  unsigned int m_gate_type;
  std::string  m_gate_name;
  uint64_t     m_trigger_timestamp;
  uint64_t     m_gate_start_timestamp;
  uint64_t     m_trigger_gate_diff;
  uint64_t     m_gate_crt_diff;

  // Histograms
  std::map<std::string, TH1F*> hDCA;
  std::map<std::string, TH1F*> hMatchDCA;
  std::map<std::string, TH1F*> hNoMatchDCA;

  std::map<std::string, TH1F*> hDoL;
  std::map<std::string, TH1F*> hMatchDoL;
  std::map<std::string, TH1F*> hNoMatchDoL;

  std::map<std::string, TH1F*> hT0;
  std::map<std::string, TH1F*> hMatchT0;
  std::map<std::string, TH1F*> hNoMatchT0;
  
};


icarus::CRTTPCMatchingAna::CRTTPCMatchingAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  , t0Alg(p.get<fhicl::ParameterSet>("t0Alg"))
  //  , evd(p.get<fhicl::ParameterSet>("evd"))
  , fCrtutils(new icarus::crt::CRTCommonUtils())
  
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  reconfigure(p);
}

void icarus::CRTTPCMatchingAna::reconfigure(fhicl::ParameterSet const & p)
{
  fCRTHitLabel        =  p.get<art::InputTag> ("CRTHitLabel", "crthit");
  fTPCTrackLabel      =  p.get< std::vector<art::InputTag> >("TPCTrackLabel",             {""});
  fPFParticleLabel    =  p.get< std::vector<art::InputTag> >("PFParticleLabel",             {""});
  fTriggerLabel       =  p.get<art::InputTag>("TriggerLabel","daqTrigger");
  //fTPCTrackLabel      =  p.get< std::vector<art::InputTag> >("TPCTrackLabel", std::vector<art::InputTag>() = {""});
  fVerbose            =  p.get<bool>("Verbose");
} // CRTT0Matching::reconfigure()

void icarus::CRTTPCMatchingAna::beginJob()
{

  // Access tfileservice to handle creating and writing histograms
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("matchTree","CRTHit - TPC track matching analysis");

  fTree->Branch("Event",           &fEvent,             "Event/I");
  fTree->Branch("SubRun",          &fSubRun,           "SubRun/I");
  fTree->Branch("Run",             &fRun,                 "Run/I");
  fTree->Branch("cryo",            "std::vector<int>",      &fcryo);
  fTree->Branch("pandorat0",       "std::vector<double>", &fpandorat0);
  fTree->Branch("crtRegion",       "std::vector<int>",&fCrtRegion);
  fTree->Branch("DCA",             "std::vector<double>",   &fDCA);
  fTree->Branch("DOL",             "std::vector<double>",   &fDOL);
  fTree->Branch("t0",                   "std::vector<double>",    &fT0);
  fTree->Branch("gate_type",            &m_gate_type,    "gate_type/b");
  fTree->Branch("gate_name",            &m_gate_name);
  fTree->Branch("trigger_timestamp",    &m_trigger_timestamp,       "trigger_timestamp/l");
  fTree->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");
  fTree->Branch("trigger_gate_diff",    &m_trigger_gate_diff,       "trigger_gate_diff/l");
  fTree->Branch("gate_crt_diff",        &m_gate_crt_diff,               "gate_crt_diff/l");

  
  for(int i = 30; i < 50 + 1; i++){
    std::string tagger = "All";
    if (i >= 35 && i < 40) continue;
    if (i==48 || i==49) continue;
    
    tagger = fCrtutils->GetRegionNameFromNum(i);//fCrtGeo.GetTagger(i).name;
    //    std::cout << "tagger: " << tagger.c_str() << std::endl;
    hDCA[tagger]     = tfs->make<TH1F>(Form("DCA_%s", tagger.c_str()),        "", 50, 0, 100);
    hDoL[tagger]     = tfs->make<TH1F>(Form("DoL_%s", tagger.c_str()),        "", 100, 0, 0.25);
    hT0[tagger]      = tfs->make<TH1F>(Form("T0_%s", tagger.c_str()),        "", 600, -3000, 3000);
  }
  
}

void icarus::CRTTPCMatchingAna::analyze(const art::Event& event)
{

  fDCA.clear();
  fDOL.clear();
  fT0.clear();
  fCrtRegion.clear();
  fpandorat0.clear();
  fcryo.clear();
  // Fetch basic event info
  if(fVerbose){
    std::cout<<"============================================"<<std::endl
	     <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
	     <<"============================================"<<std::endl;
  }

  fEvent  = event.id().event();
  fRun    = event.run();
  fSubRun = event.subRun();

  //----------------------------------------------------------------------------------------------------------
  //                                          GETTING PRODUCTS
  //----------------------------------------------------------------------------------------------------------
  //add trigger info
  if( !fTriggerLabel.empty() ) {

    art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
    event.getByLabel( fTriggerLabel, trigger_handle );
    if( trigger_handle.isValid() ) {
      sbn::triggerSource bit = trigger_handle->sourceType;
      m_gate_type = (unsigned int)bit;
      m_gate_name = bitName(bit);
      m_trigger_timestamp = trigger_handle->triggerTimestamp;
      m_gate_start_timestamp =  trigger_handle->beamGateTimestamp;
      m_trigger_gate_diff = trigger_handle->triggerTimestamp - trigger_handle->beamGateTimestamp;
    }
    else{
      mf::LogError("CRTTPCMatchingAna:") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
    }
  }
  else {
    mf::LogError("CRTTPCMatchingAna:") << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
  }


  // Get CRT hits from the event
  art::Handle< std::vector<sbn::crt::CRTHit>> crtHitHandle;
  std::vector<art::Ptr<sbn::crt::CRTHit> > crtHitList;
  if (event.getByLabel(fCRTHitLabel, crtHitHandle))
    art::fill_ptr_vector(crtHitList, crtHitHandle);

  
  //  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  //art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

  std::vector<sbn::crt::CRTHit> crtHits;
  int hit_i = 0;
  double minHitTime = 99999;
  double maxHitTime = -99999;

  for(auto const& hit : (*crtHitHandle)){
    double hitTime = double(hit.ts0_ns - (m_gate_start_timestamp%1'000'000'000))/1e3;
//'
    if(hitTime<-0.5e6)      hitTime+=1e6;
    else if(hitTime>0.5e6)  hitTime-=1e6;    
  //  double hitTime = double(m_gate_start_timestamp - hit.ts0_ns)/1e3;
//    hitTime = -hitTime+1e6;
    //double hitTime = (double)(int)hit.ts0_ns * 1e-3;
    if(hitTime < minHitTime) minHitTime = hitTime;
    if(hitTime > maxHitTime) maxHitTime = hitTime;

    crtHits.push_back(hit);
    hit_i++;
  }
			    
  mf::LogError("CRTTPCMatchingAna:") << "# of TPC tracks: " << fTPCTrackLabel.size()
				     << " \t # of CRT Hits: " << crtHits.size() << "\n" ;

  //----------------------------------------------------------------------------------------------------------
  //                                DISTANCE OF CLOSEST APPROACH ANALYSIS
  //----------------------------------------------------------------------------------------------------------

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);


  // if(fVerbose) std::cout<<"----------------- DCA Analysis -------------------"<<std::endl;

  for(const auto& trackLabel : fTPCTrackLabel)
    {
      auto it = &trackLabel - fTPCTrackLabel.data();

      //std::cout << "size of the track label: " << trackLabel << std::endl;
      // Get reconstructed tracks from the event    
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;

      //Get PFParticles
      auto pfpListHandle = event.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel[it]);
      if (!pfpListHandle.isValid()) continue;

      //      art::Handle< std::vector<recob::PFParticle> > pfpListHandle; 
      //event.getByLabel({"pandoraGausCryoE", "pandoraGausCryoW"}, pfpListHandle);

      //Get PFParticle-Track association
      art::FindManyP<recob::PFParticle> fmpfp(tpcTrackHandle, event, trackLabel); 

      //Get T0-PFParticle association
      art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, event, fPFParticleLabel[it]);


      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);
      //std::cout << "# track: " << (*tpcTrackHandle).size() << std::endl;
      // Loop over reconstructed tracks
      for (auto const& tpcTrack : (*tpcTrackHandle)){
	auto idx = &tpcTrack - (*tpcTrackHandle).data();

	double t0 = -9999999;
	//Find PFParticle for track i
	//art::Ptr::key() gives the index in the vector 
	auto pfps = fmpfp.at(tpcTrack.ID());
	
	if (!pfps.empty()){
	  //Find T0 for PFParticle
	  auto t0s = fmt0pandora.at(pfps[0].key()); 
	  if (!t0s.empty()){
	  
	    t0 = t0s[0]->Time();   //Get T0
	  }
	}

	//std::cout<< "T0: " <<  t0 << std::endl;


	//if (idx == 1) break;
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
	if (hits.size() == 0) continue;

	//std::cout<< "#hits in a tpc track: " << hits.size() << std::endl;

	int const cryoNumber = hits[0]->WireID().Cryostat;
        //      matchCand closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, crtHits, event);
        //      std::vector <matchCand> closestvec = t0Alg.GetClosestCRTHit(detProp, tpcTrack, crtHits, event);
        //matchCand closest = closestvec[idx];// closestvec.back();
	//	std::cout<< "track index: " << idx << " [ 1st, last,  cryoNumber, PandoraT0 ] = [ " 
	//	 << hits[0]->WireID().TPC << " , " << hits[hits.size()-1]->WireID().TPC
	//       << " , " << cryoNumber << " , " << t0 << " ] "<< std::endl;

	matchCand closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, m_trigger_timestamp);
	if(closest.dca >=0 )
          mf::LogInfo("CRTTPCMatchingAna")
	    << "Track # " << idx  <<" Matched time = "<<closest.t0<<" [us] to track "<< tpcTrack.ID()<<" with DCA = "<<closest.dca 
	    << " for crt region " << closest.thishit.tagger;
	/*for (long unsigned int i = 0; i < closestvec.size(); i++){

	  std::cout << "[ #closestvec , closest.dca, closest.extrapLen, closest.thishit.tagger ]=  [ " 
		    << i <<" , "<< closestvec[i].dca <<" , "<< closestvec[i].extrapLen <<" , "<< closestvec[i].thishit.tagger <<"  ]" <<std::endl;
	}*/
	double sin_angle = -99999;
	if(closest.dca != -99999){
	  hDCA[closest.thishit.tagger]->Fill(closest.dca);
	  //hDCA["All"]->Fill(closest.dca);
	  sin_angle = closest.dca/closest.extrapLen;
	  hDoL[closest.thishit.tagger]->Fill(sin_angle);
	  // hDoL["All"]->Fill(sin_angle);
	  hT0[closest.thishit.tagger]->Fill(closest.t0);
	  fDCA.push_back(closest.dca);
	  fDOL.push_back(sin_angle);
	  fT0.push_back(closest.t0);
	  fcryo.push_back(cryoNumber);
	  fpandorat0.push_back(t0);
	  // fCrtRegion.push_back(closest.thishit.tagger);
	  fCrtRegion.push_back(fCrtutils->AuxDetRegionNameToNum(closest.thishit.tagger));
	} // if closest dca is physical
      } // track loops in each tracklebel
    } // trackLabel vector

//  evd.SetDrawTpc(true);
 // evd.Draw(event);
  fTree->Fill();
}
void icarus::CRTTPCMatchingAna::endJob()
{

} // CRTT0Matching::endJob()


DEFINE_ART_MODULE(icarus::CRTTPCMatchingAna)
