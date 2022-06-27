////////////////////////////////////////////////////////////////////////
// Class:       CRTTPCEvtDisp
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTTPCEvtDisp_module.cc
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
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

// LArSoft includes
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// icaruscode includes
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/RecoUtils.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg_evtdisp.h"
//#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/Decode/DataProducts/ExtraTriggerInfo.h"

#include "canvas/Persistency/Common/FindMany.h"
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
#include "TFile.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace icarus {
  class CRTTPCEvtDisp;
}


class icarus::CRTTPCEvtDisp : public art::EDAnalyzer {
public:
  explicit CRTTPCEvtDisp(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTPCEvtDisp(CRTTPCEvtDisp const&) = delete;
  CRTTPCEvtDisp(CRTTPCEvtDisp&&) = delete;
  CRTTPCEvtDisp& operator=(CRTTPCEvtDisp const&) = delete;
  CRTTPCEvtDisp& operator=(CRTTPCEvtDisp&&) = delete;

  // Called once, at start of the job
  virtual void beginJob() override;

  // Called once per event
  virtual void analyze(const art::Event& event) override;

  // Called once, at end of the job
  virtual void endJob() override;

  //Find where a track crosses the cathode
  void simple_getCatCrossXYZ(recob::Track trk, double &my_x, double &my_y, double &my_z);
  void getCatCrossXYZ(recob::Track trk, double &my_x, double &my_y, double &my_z);

  // Calculate the distance from the track crossing point to CRT overlap coordinates
  double DistToCrtHit(TVector3 trackPos, sbn::crt::CRTHit crtHit);

  void reconfigure(fhicl::ParameterSet const & p);

private:
  art::ServiceHandle<art::TFileService> tfs;
  TTree* trevtdisp;

  icarus::crt::CRTBackTracker bt;

  bool pfp_cathodecrosser, simple_cathodecrosser, is_best_DCA_startDir, is_best_DCA_endDir;
  int hit_id, num_t0s;
  int fEvent, fRun, fSubRun;
  long track_id, meta_track_id;
  double t0min, t0max, crtTime, simpleDCA_startDir, simpleDCA_endDir;
  double startDir_x, startDir_y, startDir_z;
  double endDir_x, endDir_y, endDir_z;
  double tpc_track_start_x, tpc_track_start_y, tpc_track_start_z;
  double tpc_track_end_x, tpc_track_end_y, tpc_track_end_z;
  double crt_hit_pos_x, crt_hit_pos_y, crt_hit_pos_z;
  double catcross_x, catcross_y, catcross_z;
  double simple_catcross_x, simple_catcross_y, simple_catcross_z;
  std::vector<double> track_t0s;

  //Branch variables for MC work
  int tpc_trueID, crt_trueID;
  std::vector<double> mcpart_start_x, mcpart_start_y, mcpart_start_z;
  std::vector<double> mcpart_end_x, mcpart_end_y, mcpart_end_z;
  std::vector<double> mcpart_px, mcpart_py, mcpart_pz;
  std::vector<int> mcpart_ids;// crt_trueIDs;
  std::vector<int> mcpart_pdg;
  // Declare member data here.
  art::InputTag              fCRTHitLabel;   ///< name of CRT producer
  std::vector<art::InputTag> fTPCTrackLabel; ///< labels for source of tracks
  art::InputTag              fTriggerLabel;  ///< labels for trigger 
  art::InputTag 	     fSimModuleLabel;      ///< name of detsim producer
  std::vector<art::InputTag> fPFParticleLabel; ///< labels for source of PFParticle
  bool                       fVerbose;       ///< print information about what's going on
  bool			     fIsData;        ///< switch for if this is data or MC

  TH1D *ht0_cc;// = new TH1D("ht0_cc","Cathode Crossing Tracks T0 - All Available CRT Hit times (ns)",1000,-4000000,4000000);

  CRTT0MatchAlg_evtdisp t0Alg;

  geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
  icarus::crt::CRTCommonUtils* fCrtutils;

  //add trigger data product vars
  unsigned int m_gate_type;
  std::string  m_gate_name;
  uint64_t     m_trigger_timestamp;
  uint64_t     m_gate_start_timestamp;
  uint64_t     m_trigger_gate_diff;
  uint64_t     m_gate_crt_diff;

  int            fCrtRegion;    //CRT hit region code


};


icarus::CRTTPCEvtDisp::CRTTPCEvtDisp(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
//, fSimModuleLabel(p.get<fhicl::ParameterSet>("SimModuleLabel"))
  , bt(p.get<fhicl::ParameterSet>("CRTBackTrack"))
  , t0Alg(p.get<fhicl::ParameterSet>("t0Alg"))
  , fCrtutils(new icarus::crt::CRTCommonUtils())
  
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  reconfigure(p);
}

void icarus::CRTTPCEvtDisp::reconfigure(fhicl::ParameterSet const & p)
{
  fCRTHitLabel        =  p.get<art::InputTag> ("CRTHitLabel", "crthit");
  fTPCTrackLabel      =  p.get< std::vector<art::InputTag> >("TPCTrackLabel",             {""});
  fTriggerLabel       =  p.get<art::InputTag>("TriggerLabel","daqTrigger");
  fVerbose            =  p.get<bool>("Verbose");
  fIsData	      =  p.get<bool>("IsData");
  fSimModuleLabel     =  p.get<art::InputTag> ("SimModuleLabel","largeant");
  fPFParticleLabel    =  p.get< std::vector<art::InputTag> >("PFParticleLabel",             {""});
} // CRTT0Matching::reconfigure()

void icarus::CRTTPCEvtDisp::beginJob()
{

  // Access tfileservice to handle creating and writing histograms
  trevtdisp= tfs->make<TTree>("trevtdisp","trevtdisp");

  trevtdisp->Branch("event", &fEvent, "event/I");
  trevtdisp->Branch("run", &fRun, "run/I");
  trevtdisp->Branch("subrun", &fSubRun, "subrun/I");
  trevtdisp->Branch("pfp_cathodecrosser",&pfp_cathodecrosser,"pfp_cathodecrosser/O");
  trevtdisp->Branch("simple_cathodecrosser",&simple_cathodecrosser,"simple_cathodecrosser/O");
  trevtdisp->Branch("is_best_DCA_startDir",&is_best_DCA_startDir,"is_best_DCA_startDir/O");
  trevtdisp->Branch("is_best_DCA_endDir",&is_best_DCA_endDir,"is_best_DCA_endDir/O");
  trevtdisp->Branch("hit_id",&hit_id,"hit_id/I");
  trevtdisp->Branch("track_id",&track_id,"track_id/L");
  trevtdisp->Branch("meta_track_id",&meta_track_id,"meta_track_id/L");
  trevtdisp->Branch("t0min",&t0min,"t0min/D");
  trevtdisp->Branch("t0max",&t0max,"t0max/D");
  trevtdisp->Branch("crtTime",&crtTime,"crtTime/D");
  trevtdisp->Branch("simpleDCA_startDir",&simpleDCA_startDir,"simpleDCA_startDir/D");
  trevtdisp->Branch("simpleDCA_endDir",&simpleDCA_endDir,"simpleDCA_endDir/D");
  trevtdisp->Branch("startDir_x",&startDir_x,"startDir_x/D");
  trevtdisp->Branch("startDir_y",&startDir_y,"startDir_y/D");
  trevtdisp->Branch("startDir_z",&startDir_z,"startDir_z/D");
  trevtdisp->Branch("endDir_x",&endDir_x,"endDir_x/D");
  trevtdisp->Branch("endDir_y",&endDir_y,"endDir_y/D");
  trevtdisp->Branch("endDir_z",&endDir_z,"endDir_z/D");
  trevtdisp->Branch("tpc_track_start_x",&tpc_track_start_x,"tpc_track_start_x/D");
  trevtdisp->Branch("tpc_track_start_y",&tpc_track_start_y,"tpc_track_start_y/D");
  trevtdisp->Branch("tpc_track_start_z",&tpc_track_start_z,"tpc_track_start_z/D");
  trevtdisp->Branch("tpc_track_end_x",&tpc_track_end_x,"tpc_track_end_x/D");
  trevtdisp->Branch("tpc_track_end_y",&tpc_track_end_y,"tpc_track_end_y/D");
  trevtdisp->Branch("tpc_track_end_z",&tpc_track_end_z,"tpc_track_end_z/D");
  trevtdisp->Branch("crt_hit_pos_x",&crt_hit_pos_x,"crt_hit_pos_x/D");
  trevtdisp->Branch("crt_hit_pos_y",&crt_hit_pos_y,"crt_hit_pos_y/D");
  trevtdisp->Branch("crt_hit_pos_z",&crt_hit_pos_z,"crt_hit_pos_z/D");
  trevtdisp->Branch("crtRegion",  &fCrtRegion,"crtRegion/I");

  trevtdisp->Branch("tpc_trueID",&tpc_trueID,"tpc_trueID/I");
  trevtdisp->Branch("crt_trueID",&crt_trueID,"crt_trueID/I");
  trevtdisp->Branch("mcpart_ids",&mcpart_ids);
  trevtdisp->Branch("track_t0s",&track_t0s);

  trevtdisp->Branch("num_t0s",&num_t0s,"num_t0s/I");

  trevtdisp->Branch("mcpart_start_x",&mcpart_start_x);
  trevtdisp->Branch("mcpart_start_y",&mcpart_start_y);
  trevtdisp->Branch("mcpart_start_z",&mcpart_start_z);
  trevtdisp->Branch("mcpart_end_x",&mcpart_end_x);
  trevtdisp->Branch("mcpart_end_y",&mcpart_end_y);
  trevtdisp->Branch("mcpart_end_z",&mcpart_end_z);
  trevtdisp->Branch("mcpart_px",&mcpart_px);
  trevtdisp->Branch("mcpart_py",&mcpart_py);
  trevtdisp->Branch("mcpart_pz",&mcpart_pz);
  trevtdisp->Branch("mcpart_pdg",&mcpart_pdg);

  trevtdisp->Branch("catcross_x",&catcross_x,"catcross_x/D");
  trevtdisp->Branch("catcross_y",&catcross_y,"catcross_y/D");
  trevtdisp->Branch("catcross_z",&catcross_z,"catcross_z/D");

  trevtdisp->Branch("simple_catcross_x",&simple_catcross_x,"simple_catcross_x/D");
  trevtdisp->Branch("simple_catcross_y",&simple_catcross_y,"simple_catcross_y/D");
  trevtdisp->Branch("simple_catcross_z",&simple_catcross_z,"simple_catcross_z/D");

  meta_track_id=0;
  for(int i = 30; i < 50 + 1; i++){
    std::string tagger = "All";
    if (i>=35 && i<40) continue;
    if (i==48 || i==49) continue;
    // if(i < ){
    tagger = fCrtutils->GetRegionNameFromNum(i);//fCrtGeo.GetTagger(i).name;
  }

  ht0_cc = tfs->make<TH1D>("h_t0cc","Track T0 - CRT Hit Time (ns)",1000,-5e6,5e6);

  //*/
}

void icarus::CRTTPCEvtDisp::analyze(const art::Event& event)
{

  fEvent  = event.id().event();
  fRun    = event.run();
  fSubRun = event.subRun();


  if(!fIsData){  
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    bt.Initialize(event);
    std::map<int, simb::MCParticle> particles;

    // Loop over the true particles
    mcpart_ids.clear();
    for (auto const& particle: (*particleHandle)){
      
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;
      std::cout << "Particle ID: " << partID << std::endl;
      mcpart_ids.push_back(partID);

      mcpart_start_x.push_back(particle.Position().X());
      mcpart_start_y.push_back(particle.Position().Y());
      mcpart_start_z.push_back(particle.Position().Z());

      mcpart_end_x.push_back(particle.EndPosition().X());
      mcpart_end_y.push_back(particle.EndPosition().Y());
      mcpart_end_z.push_back(particle.EndPosition().Z());

      mcpart_px.push_back(particle.Momentum().Px());
      mcpart_py.push_back(particle.Momentum().Py());
      mcpart_pz.push_back(particle.Momentum().Pz());

      mcpart_pdg.push_back(particle.PdgCode());
    }//end loop over particles assigning IDs
  }//end if(!fIsData)
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
    }//end if( trigger_handle.isValid() )
    else{
      mf::LogError("CRTTPCEvtDisp:") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
    }//end else
  }//end if( !fTriggerLabel.empty() )
  else {
    mf::LogError("CRTTPCEvtDisp:") << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
  }//end else


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
    double hitTime = double(m_gate_start_timestamp - hit.ts0_ns)/1e3;
    hitTime = -hitTime+1e6;
    //double hitTime = (double)(int)hit.ts0_ns * 1e-3;
    if(hitTime < minHitTime) minHitTime = hitTime;
    if(hitTime > maxHitTime) maxHitTime = hitTime;

    crtHits.push_back(hit);
    hit_i++;
  }

  //Temp for debugging MC
  std::cout << "# of TPC tracks:\t" << fTPCTrackLabel.size();
  std::cout << "# of CRT Hits:\t" << crtHits.size();

  //Begin building association between TPC tracks in the event and their T0s (if they have them)
  art::Handle<std::vector<recob::Track>> pandoratrkHandle;
  std::vector<art::Ptr<recob::Track>> pandoratrks;
  if(event.getByLabel("pandoraTrackGausCryoW",pandoratrkHandle)){ art::fill_ptr_vector(pandoratrks,pandoratrkHandle); }//end if(event.getByLabel("pandoraTrackGausCryoW",pandoratrkHandle)

  art::Handle< std::vector<recob::PFParticle>> pfpListHandle;
  event.getByLabel("pandoraGausCryoW",pfpListHandle);

  art::FindManyP<recob::PFParticle> fmpfp(pandoratrkHandle, event, "pandoraTrackGausCryoW");
  art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, event, "pandoraGausCryoW");

  for(size_t i=0; i<pandoratrks.size(); ++i){}//end loop over pandoratrks.size()

  //----------------------------------------------------------------------------------------------------------
  //                                DISTANCE OF CLOSEST APPROACH ANALYSIS
  //----------------------------------------------------------------------------------------------------------cc

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);


  // if(fVerbose) std::cout<<"----------------- DCA Analysis -------------------"<<std::endl;

  icarus::match_geometry thiscand;

  int temp_track_id =0;
  for(const auto& trackLabel : fTPCTrackLabel)
    {
      auto it = &trackLabel - fTPCTrackLabel.data();

      // Get reconstructed tracks from the event    
      auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(trackLabel);
      if (!tpcTrackHandle.isValid()) continue;

      std::cout << "New TPC track label!\n";  

      art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, trackLabel);

      //Begin building association between TPC tracks in the event and their T0s (if they have them)
      std::cout << "size of the track label: " << trackLabel << std::endl;

      // Get reconstructed tracks from the event
      //Get PFParticles
      //Get PFParticles
      auto pfpListHandle = event.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel[it]);
      if (!pfpListHandle.isValid()) continue;


//      art::Handle< std::vector<recob::PFParticle> > pfpListHandle;
//      event.getByLabel("pandoraGausCryoE", pfpListHandle);

      //Get PFParticle-Track association
      art::FindManyP<recob::PFParticle> fmpfp(tpcTrackHandle, event, trackLabel);

      //Get T0-PFParticle association
      art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, event, fPFParticleLabel[it]);

      std::cout << "# track: " << (*tpcTrackHandle).size() << std::endl;

      // Loop over reconstructed tracks
      for (auto const& tpcTrack : (*tpcTrackHandle)){

	track_t0s.clear();

	std::cout << "New TPC track!\n";

        //Find PFParticle for track i
        //art::Ptr::key() gives the index in the vector
        auto pfps = fmpfp.at(tpcTrack.ID());

	if(pfps.empty()) num_t0s = 0;

	pfp_cathodecrosser = false;
        if (!pfps.empty()){ 
           //Find T0 for PFParticle
           auto t0s = fmt0pandora.at(pfps[0].key());
	   num_t0s = t0s.size();
           if (!t0s.empty()){ //Get T0
	      for(size_t i=0; i<t0s.size(); i++){
                 track_t0s.push_back(t0s[0]->Time());
		 pfp_cathodecrosser = true;
	 	 catcross_x = DBL_MAX;
		 catcross_y = DBL_MAX;
		 catcross_z = DBL_MAX;
		 simple_catcross_x = DBL_MAX;
		 simple_catcross_y = DBL_MAX;
		 simple_catcross_z = DBL_MAX;
		 getCatCrossXYZ(tpcTrack, catcross_x, catcross_y, catcross_z);
		 simple_getCatCrossXYZ(tpcTrack, simple_catcross_x, simple_catcross_y, simple_catcross_z);
	      }//end loop over t0s
	   }//end if(!t0s.empty())
	}//end if(!pfps.empty())
	
	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

	if(!fIsData){
		//Using the vector of art Ptrs to recob::Hit, we can use RecoUtils.* to get the likely truth ID:
		int track_trueID = RecoUtils::TrueParticleIDFromTotalTrueEnergy(clockData,hits,false);
		std::cout << "TPC track true ID = " << track_trueID << "\n";
		tpc_trueID = track_trueID;
	}//end if(fIsData)
	std::vector<match_geometry> closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, m_trigger_timestamp, event, fIsData);

	if(!pfps.empty()){
           	auto t0s = fmt0pandora.at(pfps[0].key());
		if(!t0s.empty()){
			for(int i=0;i<(int)crtHits.size(); i++){

				double temp_crttime; 
				if(fIsData) { temp_crttime = double(m_trigger_timestamp - crtHits[i].ts0_ns);
					      temp_crttime = -temp_crttime+1e9; 				}//end if(fIsDataData)
				else {	      temp_crttime = 1600 - crtHits[i].ts0_ns/1e3; }

				double subtracted_time = t0s[0]->Time() - temp_crttime;
				std::cout << "SHOULD SEE ht0_cc->Fill(): " << subtracted_time <<"ns\n";
				ht0_cc->Fill(subtracted_time);
			}//end loop over all CRT Hits
		}//end if(!t0s.empty)
	}//end if(!pfps.empty() && !t0s.empty())
  	std::cout << "found " << closest.size() << " possible candidates\n";
	for(int i=0; i<(int)closest.size(); i++){

		sbn::crt::CRTHit checkhit = closest[i].thishit;

		if(!fIsData) crt_trueID = bt.TrueIdFromTotalEnergy(event,checkhit);
		thiscand = closest[i];
		simple_cathodecrosser = thiscand.simple_cathodecrosser;
		is_best_DCA_startDir = thiscand.is_best_DCA_startDir;
		is_best_DCA_endDir = thiscand.is_best_DCA_endDir;
		hit_id = thiscand.hit_id;
		track_id = temp_track_id; 
		t0min = thiscand.t0min;
		t0max = thiscand.t0max;
		simpleDCA_startDir = thiscand.simpleDCA_startDir;

		crtTime = thiscand.crtTime;// std::cout << "crtTime (us): " << crtTime << std::endl;
		simpleDCA_endDir = thiscand.simpleDCA_endDir;

		startDir_x = thiscand.startDir.X();
		startDir_y = thiscand.startDir.Y();
		startDir_z = thiscand.startDir.Z();

		endDir_x = thiscand.endDir.X();
		endDir_y = thiscand.endDir.Y();
		endDir_z = thiscand.endDir.Z();

		tpc_track_start_x = thiscand.tpc_track_start.X();
		tpc_track_start_y = thiscand.tpc_track_start.Y();
		tpc_track_start_z = thiscand.tpc_track_start.Z();

		tpc_track_end_x = thiscand.tpc_track_end.X();
		tpc_track_end_y = thiscand.tpc_track_end.Y();
		tpc_track_end_z = thiscand.tpc_track_end.Z();

		crt_hit_pos_x = thiscand.crt_hit_pos.X();
		crt_hit_pos_y = thiscand.crt_hit_pos.Y();
		crt_hit_pos_z = thiscand.crt_hit_pos.Z();

		fCrtRegion=(int)fCrtutils->AuxDetRegionNameToNum(thiscand.thishit.tagger);

		trevtdisp->Fill();

	}//end loop over vector of CRTHit-TPC track candidate pairs
		meta_track_id++;
		temp_track_id++;
      } // track loops in each tracklebel
    } // end loop over TPC tracks

}
void icarus::CRTTPCEvtDisp::endJob()
{
//	TFile *ftemp = new TFile("tempfile.root","update");
//	ht0_cc->Write();
//	ftemp->Close();
} // CRTT0Matching::endJob()

void icarus::CRTTPCEvtDisp::simple_getCatCrossXYZ(recob::Track trk, double &my_x, double &my_y, double &my_z){

	size_t ntrk = trk.NPoints();
	std::vector<double> x, y, z, x_diff;
//	double xsum=0, xavg;
//	bool isEast;
	std::pair<double,double> cathode_yz;
	double dist_min = DBL_MAX; int dist_min_pos = INT_MAX;

	//Begin by extracting coordinates into arrays
	for(size_t i=0; i<ntrk; i++){
		geo::Point_t thispt = trk.LocationAtPoint((int)i);
		x.push_back(thispt.X()); y.push_back(thispt.Y()); z.push_back(thispt.Z());
//		xsum+=x[i];
		x_diff.push_back(std::abs(x[i]) - 210.215);
		if(x_diff[i]<=dist_min) { dist_min=x_diff[i]; dist_min_pos = i; }
	}//end loop over track points
//	xavg=xsum/(int)ntrk;
//	if(xavg>0) isEast = false;
//	else if(xavg<0) isEast = true;
//	else if(xavg==0) { std::cout << "Average X of a track ==0?\n"; break; }

	my_x = x[dist_min_pos];
	my_y = y[dist_min_pos];
	my_z = z[dist_min_pos];

}//end definition of std::pair<int,int> getCatCrossYZ(std::vector<art::Ptr<recob::Hit>> trk_hits)
void icarus::CRTTPCEvtDisp::getCatCrossXYZ(recob::Track trk, double &my_x, double &my_y, double &my_z){

	size_t ntrk = trk.NPoints();
	std::vector<double> x, y, z, x_left_diff, x_right_diff;
	std::pair<double,double> cathode_yz;
	double left_dist_min = DBL_MAX; int left_dist_min_pos = INT_MAX;
	double right_dist_min = DBL_MAX; int right_dist_min_pos = INT_MAX;
	double xsum=0, xavg;

	//Begin by extracting coordinates into arrays
	for(size_t i=0; i<ntrk; i++){
		geo::Point_t thispt = trk.LocationAtPoint((int)i);
		x.push_back(thispt.X()); y.push_back(thispt.Y()); z.push_back(thispt.Z());
		if(std::abs(x.back())<210.215) {
			x_left_diff.push_back(std::abs(std::abs(x[i]) - 210.25));
			if(x_left_diff.back()<=left_dist_min) { left_dist_min=x_left_diff.back(); left_dist_min_pos = i; }
		}//end if(std::abs(x[i])<210.215)
		else if(std::abs(x.back())>=210.215) {
			x_right_diff.push_back(std::abs(std::abs(x[i]) - 210.25));
			if(x_right_diff.back()<=right_dist_min) { right_dist_min=x_right_diff.back(); right_dist_min_pos = i; }
		}//end if(std::abs(x[i])<210.215)
		xsum+=x.back();
	}//end loop over track points
	xavg = xsum/ntrk;

	TVector3 leftpt(x[left_dist_min_pos],y[left_dist_min_pos],z[left_dist_min_pos]);
	TVector3 rightpt(x[right_dist_min_pos],y[right_dist_min_pos],z[right_dist_min_pos]);

	//calculate parameters for a line that passes through the two points:
	double mxy = (leftpt.Y() - rightpt.Y())/(leftpt.X()-rightpt.X());
	double mxz = (leftpt.Z() - rightpt.Z())/(leftpt.X()-rightpt.X());

	double bxy = leftpt.Y() - mxy*leftpt.X();
	double bxz = leftpt.Z() - mxz*leftpt.X();

	if(xavg>0) my_x = 210.215; 
	else if(xavg<=0) my_x = -210.215;
	my_y = mxy*my_x + bxy;
	my_z = mxz*my_x + bxz;

}//end definition of std::pair<int,int> getCatCrossYZ(std::vector<art::Ptr<recob::Hit>> trk_hits)

DEFINE_ART_MODULE(icarus::CRTTPCEvtDisp)
