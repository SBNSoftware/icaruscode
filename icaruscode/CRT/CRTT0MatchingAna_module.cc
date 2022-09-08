////////////////////////////////////////////////////////////////////////
// Class:       CRTT0MatchingAna
// Module Type: analyzer
// File:        CRTT0MatchingAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// icaruscode includes
#include "icaruscode/CRT/CRTUtils/RecoUtils.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
//#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
//#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace icarus {
 
  class CRTT0MatchingAna : public art::EDAnalyzer {
  public:

    explicit CRTT0MatchingAna(fhicl::ParameterSet const& p);
    
    // Called once, at start of the job
    virtual void beginJob() override;
    
    // Called once per event
    virtual void analyze(const art::Event& event) override;
    
    // Called once, at end of the job
    virtual void endJob() override;
    
    // Calculate the distance from the track crossing point to CRT overlap coordinates
    double DistToCrtHit(TVector3 trackPos, sbn::crt::CRTHit crtHit);
    
    void reconfigure(fhicl::ParameterSet const & p);

  private:
    
    // fcl file parameters
    art::InputTag              fSimModuleLabel;///< name of detsim producer
    art::InputTag              fCRTHitLabel;   ///< name of CRT producer
    std::vector<art::InputTag> fTPCTrackLabel; ///< name of CRT producer
    art::InputTag              fTriggerLabel;  ///< labels for trigger
    bool                       fVerbose;       ///< print information about what's going on
    
    CRTT0MatchAlg t0Alg;
    
    //  CRTGeoAlg fCrtGeo;
    //TPCGeoAlg fTpcGeo;
    
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    icarus::crt::CRTCommonUtils* fCrtutils;
    //CRTBackTracker fCrtBackTrack;
    icarus::crt::CRTBackTracker bt;
    // Histograms
    std::map<std::string, TH1D*> hDCA;
    std::map<std::string, TH1D*> hMatchDCA;
    std::map<std::string, TH1D*> hNoMatchDCA;
    
    std::map<std::string, TH1D*> hDoL;
    std::map<std::string, TH1D*> hMatchDoL;
    std::map<std::string, TH1D*> hNoMatchDoL;
    
    std::map<std::string, TH1D*> hT0;
    std::map<std::string, TH1D*> hMatchT0;
    std::map<std::string, TH1D*> hNoMatchT0;
    
    std::map<std::string, TH1D*> hEffDCATotal;
    std::map<std::string, TH1D*> hEffDCAReco;
    std::map<std::string, TH1D*> hEffDoLTotal;
    std::map<std::string, TH1D*> hEffDoLReco;
    std::map<std::string, TH1D*> hEffLengthTotal;
    std::map<std::string, TH1D*> hEffLengthReco;
    
    std::map<std::string, TH1D*> hPurityDCATotal;
    std::map<std::string, TH1D*> hPurityDCAReco;
    std::map<std::string, TH1D*> hPurityDoLTotal;
    std::map<std::string, TH1D*> hPurityDoLReco;
    std::map<std::string, TH1D*> hPurityLengthTotal;
    std::map<std::string, TH1D*> hPurityLengthReco;

    //add trigger data product vars
    unsigned int m_gate_type;
    std::string  m_gate_name;
    uint64_t     m_trigger_timestamp;
    uint64_t     m_gate_start_timestamp;
    uint64_t     m_trigger_gate_diff;
    uint64_t     m_gate_crt_diff;    
  }; // class CRTT0MatchingAna


  CRTT0MatchingAna::CRTT0MatchingAna(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // ,
    , t0Alg(p.get<fhicl::ParameterSet>("t0Alg"))
    , fCrtutils(new icarus::crt::CRTCommonUtils())
    , bt(p.get<fhicl::ParameterSet>("CRTBackTrack"))
    // More initializers here.
    {
      // Call appropriate consumes<>() for any products to be retrieved by this module.
      fGeometryService = lar::providerFrom<geo::Geometry>();
      reconfigure(p);
    }


  // Constructor
  void CRTT0MatchingAna::reconfigure(fhicl::ParameterSet const& p)
  {
    
    fSimModuleLabel     =  p.get<art::InputTag> ("SimModuleLabel", "largeant");
    fCRTHitLabel        =  p.get<art::InputTag> ("CRTHitLabel", "crthit");
    fTPCTrackLabel      =  p.get< std::vector<art::InputTag> >("TPCTrackLabel",             {""});
    fTriggerLabel       =  p.get<art::InputTag>("TriggerLabel","daqTrigger");
    fVerbose            =  p.get<bool>("Verbose");
  } //CRTT0MatchingAna()


  void CRTT0MatchingAna::beginJob()
  {

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    for(int i = 30; i < 50 + 1; i++){
      std::string tagger = "All";
      if (i > 34 && i < 40) continue;
      if (i==48 || i==49) continue;
      // if(i < ){
      tagger = fCrtutils->GetRegionNameFromNum(i);//fCrtGeo.GetTagger(i).name;
      std::cout << "tagger: " << tagger.c_str() << std::endl;
      hDCA[tagger]        = tfs->make<TH1D>(Form("DCA_%s", tagger.c_str()),        "", 50, 0, 100);
      hMatchDCA[tagger]   = tfs->make<TH1D>(Form("MatchDCA_%s", tagger.c_str()),   "", 50, 0, 100);
      hNoMatchDCA[tagger] = tfs->make<TH1D>(Form("NoMatchDCA_%s", tagger.c_str()), "", 50, 0, 100);

      hDoL[tagger]        = tfs->make<TH1D>(Form("DoL_%s", tagger.c_str()),        "", 100, 0, 0.25);
      hMatchDoL[tagger]   = tfs->make<TH1D>(Form("MatchDoL_%s", tagger.c_str()),   "", 100, 0, 0.25);
      hNoMatchDoL[tagger] = tfs->make<TH1D>(Form("NoMatchDoL_%s", tagger.c_str()), "", 100, 0, 0.25);
      
      hT0[tagger]        = tfs->make<TH1D>(Form("T0_%s", tagger.c_str()),        "", 600, -3000, 3000);
      hMatchT0[tagger]   = tfs->make<TH1D>(Form("MatchT0_%s", tagger.c_str()),   "", 600, -3000, 3000);
      hNoMatchT0[tagger] = tfs->make<TH1D>(Form("NoMatchT0_%s", tagger.c_str()), "", 600, -3000, 3000);
      
      hEffDCATotal[tagger] = tfs->make<TH1D>(Form("EffDCATotal_%s", tagger.c_str()), "", 50, 0, 100);
      hEffDCAReco[tagger]  = tfs->make<TH1D>(Form("EffDCAReco_%s", tagger.c_str()),  "", 50, 0, 100);
      hEffDoLTotal[tagger] = tfs->make<TH1D>(Form("EffDoLTotal_%s", tagger.c_str()), "", 100, 0, 0.25);
      hEffDoLReco[tagger]  = tfs->make<TH1D>(Form("EffDoLReco_%s", tagger.c_str()),  "", 100, 0, 0.25);
      hEffLengthTotal[tagger] = tfs->make<TH1D>(Form("EffLengthTotal_%s", tagger.c_str()), "", 20, 0, 600);
      hEffLengthReco[tagger]  = tfs->make<TH1D>(Form("EffLengthReco_%s", tagger.c_str()),  "", 20, 0, 600);
    
      hPurityDCATotal[tagger] = tfs->make<TH1D>(Form("PurityDCATotal_%s", tagger.c_str()), "", 50, 0, 100);
      hPurityDCAReco[tagger]  = tfs->make<TH1D>(Form("PurityDCAReco_%s", tagger.c_str()),  "", 50, 0, 100);
      hPurityDoLTotal[tagger] = tfs->make<TH1D>(Form("PurityDoLTotal_%s", tagger.c_str()), "", 100, 0, 0.25);
      hPurityDoLReco[tagger]  = tfs->make<TH1D>(Form("PurityDoLReco_%s", tagger.c_str()),  "", 100, 0, 0.25);
      hPurityLengthTotal[tagger] = tfs->make<TH1D>(Form("PurityLengthTotal_%s", tagger.c_str()), "", 20, 0, 600);
      hPurityLengthReco[tagger]  = tfs->make<TH1D>(Form("PurityLengthReco_%s", tagger.c_str()),  "", 20, 0, 600);
    }

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT T0 Matching Ana Module -------------------"<<std::endl;

  } // CRTT0MatchingAna::beginJob()


  void CRTT0MatchingAna::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

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
	mf::LogError("CRTT0MatchingAna:") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
      }
    }
    else {
      mf::LogError("CRTT0MatchingAna:") << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
    }

    // Get g4 particles
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Get CRT hits from the event
    art::Handle< std::vector<sbn::crt::CRTHit>> crtHitHandle;
    std::vector<art::Ptr<sbn::crt::CRTHit> > crtHitList;
    if (event.getByLabel(fCRTHitLabel, crtHitHandle))
      art::fill_ptr_vector(crtHitList, crtHitHandle);

    // Get reconstructed tracks from the event
    //auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    //    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

    // fCrtBackTrack.Initialize(event);

    //----------------------------------------------------------------------------------------------------------
    //                                          TRUTH MATCHING
    //----------------------------------------------------------------------------------------------------------
    
    std::map<int, simb::MCParticle> particles;
    // Loop over the true particles
    for (auto const& particle: (*particleHandle)){
      
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;
      
    }

    if(fVerbose) std::cout<<"----------------- truth matching -------------------"<<std::endl;
    std::map<int, std::vector<std::string>> crtTaggerMap;
    std::vector<sbn::crt::CRTHit> crtHits;
    int hit_i = 0;
    double minHitTime = 99999;
    double maxHitTime = -99999;

    for(auto const& hit : (*crtHitHandle)){
      double hitTime = double(hit.ts0_ns - (m_gate_start_timestamp%1'000'000'000))/1e3;
      if(hitTime<-0.5e6)      hitTime+=1e6;
      else if(hitTime>0.5e6)  hitTime-=1e6;

      //      double hitTime = (double)(int)hit.ts1_ns * 1e-3;
      if(hitTime < minHitTime) minHitTime = hitTime;
      if(hitTime > maxHitTime) maxHitTime = hitTime;

      crtHits.push_back(hit);
      int trueId = bt.TrueIdFromHitId(event, hit_i);
      hit_i++;
      if(trueId == -99999) continue;
      if(crtTaggerMap.find(trueId) == crtTaggerMap.end()){
        crtTaggerMap[trueId].push_back(hit.tagger);
      }
      else if(std::find(crtTaggerMap[trueId].begin(), crtTaggerMap[trueId].end(), hit.tagger) == crtTaggerMap[trueId].end()){
        crtTaggerMap[trueId].push_back(hit.tagger);
      }
    }

    //    std::cout << " New event" << std::endl;
    //----------------------------------------------------------------------------------------------------------
    //                                DISTANCE OF CLOSEST APPROACH ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

    // if(fVerbose) std::cout<<"----------------- DCA Analysis -------------------"<<std::endl;
    for(const auto& trackLabel : fTPCTrackLabel)
      {
	auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
	if (!tpcTrackHandle.isValid()) continue;
	
	art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);

	// Loop over reconstructed tracks
	for (auto const& tpcTrack : (*tpcTrackHandle)){
	  //      std::cout<<"----------------- why only 4 times -------------------" <<std::endl;
	  // if(fVerbose) std::cout<<"----------------- # track -------------------" << (*tpcTrackHandle).size()<<std::endl;
	  // Get the associated hits
	  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
	  int trackTrueID = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);
	  if(particles.find(trackTrueID) == particles.end()) continue;
	  // Only consider primary muons
	  if(!(std::abs(particles[trackTrueID].PdgCode()) == 13 && particles[trackTrueID].Mother() == 0)) continue;
	  
	  // Only consider particles inside the reco hit time window
	  double trueTime = particles[trackTrueID].T() * 1e-3;
	  if(trueTime < minHitTime || trueTime > maxHitTime) continue;
	  //if(fVerbose) std::cout<<"----------------- line 315 -------------------"<<std::endl;
	  std::cout << "new track " << trueTime << std::endl;
	  // Calculate t0 from CRT Hit matching
	  matchCand closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, m_gate_start_timestamp);
	  // matchCand closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, crtHits, event);
	  //std::vector <matchCand> closestvec = t0Alg.GetClosestCRTHit(detProp, tpcTrack, crtHits, event);
          //matchCand closest = closestvec.back();
	  //std::cout << "closest match " << closest.t0 << std::endl;
	  //std::cout << "closest dca " << closest.dca << " ,tagger: "<< closest.thishit.tagger << " , sin angele: \t" << closest.dca/closest.extrapLen << std::endl;
	  double sin_angle = -99999;
	  if(closest.dca != -99999){ 
	    hDCA[closest.thishit.tagger]->Fill(closest.dca);
	    //hDCA["All"]->Fill(closest.dca);
	    sin_angle = closest.dca/closest.extrapLen;
	    hDoL[closest.thishit.tagger]->Fill(sin_angle);
	    // hDoL["All"]->Fill(sin_angle);
	    hT0[closest.thishit.tagger]->Fill(closest.t0);
	    //hT0["All"]->Fill(closest.t0);
	    //std::cout<< "tagger: "<< closest.thishit.tagger << " , sin angele: \t" << sin_angle << std::endl;
	    
	    // Is hit matched to that track
	    int hitTrueID = bt.TrueIdFromTotalEnergy(event, closest.thishit);
	    if(hitTrueID == trackTrueID && hitTrueID != -99999){
	      hMatchDCA[closest.thishit.tagger]->Fill(closest.dca);
	      //hMatchDCA["All"]->Fill(closest.dca);
	      hMatchDoL[closest.thishit.tagger]->Fill(sin_angle);
	      //hMatchDoL["All"]->Fill(sin_angle);
	      hMatchT0[closest.thishit.tagger]->Fill(closest.t0);
	      //hMatchT0["All"]->Fill(closest.t0);
	    }
	    else{
	      hNoMatchDCA[closest.thishit.tagger]->Fill(closest.dca);
	      //hNoMatchDCA["All"]->Fill(closest.dca);
	      hNoMatchDoL[closest.thishit.tagger]->Fill(sin_angle);
	      // hNoMatchDoL["All"]->Fill(sin_angle);
	      hNoMatchT0[closest.thishit.tagger]->Fill(closest.t0);
	      // hNoMatchT0["All"]->Fill(closest.t0);
	    }
	    //std::cout<< "---------> 2nd line tagger: "<< closest.thishit.tagger << " , sin angele: \t" << sin_angle << std::endl;
	  }else continue;
	  //if(fVerbose) std::cout<<"----------------- line 350 -------------------"<<std::endl;
	  int nbins = hEffDCATotal.begin()->second->GetNbinsX();
	  for(int i = 0; i < nbins; i++){
	    double DCAcut = hEffDCATotal.begin()->second->GetBinCenter(i);
	    
	    // Fill total efficiency histogram with each cut if track matches any hits
	    if(crtTaggerMap.find(trackTrueID) != crtTaggerMap.end()){
	      for(auto const& tagger : crtTaggerMap[trackTrueID]){
		hEffDCATotal[tagger]->Fill(DCAcut);
		
		// If closest hit is below limit and track matches any hits then fill efficiency
		if(closest.dca < DCAcut && closest.dca != -99999){
		  hEffDCAReco[tagger]->Fill(DCAcut);
		}
	      }
	      // Fill total efficiency histograms
	      // hEffDCATotal["All"]->Fill(DCAcut);
	      if(closest.dca < DCAcut && closest.dca != -99999){
		// hEffDCAReco["All"]->Fill(DCAcut);
	      }
	    }
	    
	    // Fill total purity histogram with each cut if closest hit is below limit
	    if(closest.dca < DCAcut && closest.dca != -99999){
	      hPurityDCATotal[closest.thishit.tagger]->Fill(DCAcut);
	      // hPurityDCATotal["All"]->Fill(DCAcut);
	      
	      // If closest hit is below limit and matched time is correct then fill purity
	      double hitTime = closest.thishit.ts1_ns * 1e-3;
	      if(particles.find(trackTrueID) != particles.end()){
		if(std::abs(hitTime - trueTime) < 2.){
		  hPurityDCAReco[closest.thishit.tagger]->Fill(DCAcut);
		  //  hPurityDCAReco["All"]->Fill(DCAcut);
		}
	      }
	      
	    }
	  }
	  
	  nbins = hEffDoLTotal.begin()->second->GetNbinsX();
	  //if(fVerbose) std::cout<<"----------------- line 390 -------------------"<<std::endl;
	  for(int i = 0; i < nbins; i++){
	    double DCAcut = hEffDoLTotal.begin()->second->GetBinCenter(i);
	    
	    // Fill total efficiency histogram with each cut if track matches any hits
	    if(crtTaggerMap.find(trackTrueID) != crtTaggerMap.end()){
	      for(auto const& tagger : crtTaggerMap[trackTrueID]){
		hEffDoLTotal[tagger]->Fill(DCAcut);
		
		// If closest hit is below limit and track matches any hits then fill efficiency
		if(sin_angle < DCAcut && closest.dca != -99999){
		  hEffDoLReco[tagger]->Fill(DCAcut);
		}
	      }
	      // Fill total efficiency histograms
	      //  hEffDoLTotal["All"]->Fill(DCAcut);
	      if(sin_angle < DCAcut && closest.dca != -99999){
		// hEffDoLReco["All"]->Fill(DCAcut);
	      }
	    }
	    
	    // Fill total purity histogram with each cut if closest hit is below limit
	    if(sin_angle < DCAcut && closest.dca != -99999){
	      hPurityDoLTotal[closest.thishit.tagger]->Fill(DCAcut);
	      // hPurityDoLTotal["All"]->Fill(DCAcut);
	      
	      // If closest hit is below limit and matched time is correct then fill purity
	      double hitTime = closest.thishit.ts1_ns * 1e-3;
	      if(particles.find(trackTrueID) != particles.end()){
		if(std::abs(hitTime - trueTime) < 2.){
		  hPurityDoLReco[closest.thishit.tagger]->Fill(DCAcut);
		  //  hPurityDoLReco["All"]->Fill(DCAcut);
		}
	      }
	      
	    }
	  }
	  
	  double fixedCut = 30.;
	  
	  // Fill total efficiency histogram with each cut if track matches any hits
	  if(crtTaggerMap.find(trackTrueID) != crtTaggerMap.end()){
	    for(auto const& tagger : crtTaggerMap[trackTrueID]){
	      hEffLengthTotal[tagger]->Fill(tpcTrack.Length());
	      
	      // If closest hit is below limit and track matches any hits then fill efficiency
	      if(closest.dca < fixedCut && closest.dca >=0 ){
		hEffLengthReco[tagger]->Fill(tpcTrack.Length());
	      }
	    }
	    // Fill total efficiency histograms
	    // hEffLengthTotal["All"]->Fill(tpcTrack.Length());
	    if(closest.dca < fixedCut && closest.dca >=0){
	      //hEffLengthReco["All"]->Fill(tpcTrack.Length());
	    }
	  }
	  
	  // Fill total purity histogram with each cut if closest hit is below limit
	  if(closest.dca < fixedCut && closest.dca >= 0){
	    hPurityLengthTotal[closest.thishit.tagger]->Fill(tpcTrack.Length());
	    // hPurityLengthTotal["All"]->Fill(tpcTrack.Length());
	    
	    // If closest hit is below limit and matched time is correct then fill purity
	    double hitTime = closest.thishit.ts1_ns * 1e-3;
	    if(particles.find(trackTrueID) != particles.end()){
	      double trueTime = particles[trackTrueID].T() * 1e-3;
	      if(std::abs(hitTime - trueTime) < 2.){
		hPurityLengthReco[closest.thishit.tagger]->Fill(tpcTrack.Length());
		//  hPurityLengthReco["All"]->Fill(tpcTrack.Length());
	      }
	    }
	    
	  }
	  // if(fVerbose) std::cout<<"----------------- line 463 -------------------"<<std::endl;
	}
      } // track lable
	// if(fVerbose) std::cout<<"----------------- line 466 -------------------"<<std::endl;
    
  } // CRTT0MatchingAna::analyze()
  
  
  void CRTT0MatchingAna::endJob(){
    
    
  } // CRTT0MatchingAna::endJob()
  
  
  DEFINE_ART_MODULE(CRTT0MatchingAna)
} // namespace icarus

