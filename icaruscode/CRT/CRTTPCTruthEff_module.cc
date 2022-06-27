////////////////////////////////////////////////////////////////////////
// Class:       CRTTPCTruthEff
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTTPCTruthEff_module.cc
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
#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"
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
  class CRTTPCTruthEff;
}

class icarus::CRTTPCTruthEff : public art::EDAnalyzer {
public:
  explicit CRTTPCTruthEff(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTPCTruthEff(CRTTPCTruthEff const&) = delete;
  CRTTPCTruthEff(CRTTPCTruthEff&&) = delete;
  CRTTPCTruthEff& operator=(CRTTPCTruthEff const&) = delete;
  CRTTPCTruthEff& operator=(CRTTPCTruthEff&&) = delete;

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
  art::ServiceHandle<art::TFileService> tfs;
  TTree* trevtdisp;

  icarus::crt::CRTBackTracker bt;

  int fEvent, fRun, fSubRun;

  int track_trueID, crt_trueID, ttl_tpctrks, ttl_crthits;
  int crt_pdg, trk_pdg, track_mother, crt_mother, track_ancestor, crt_ancestor;
  int crt_motherlayers, track_motherlayers;
  double this_dca;

  // Declare member data here.
  art::InputTag              fCRTHitLabel;   ///< name of CRT producer
  std::vector<art::InputTag> fTPCTrackLabel; ///< labels for source of tracks
  art::InputTag              fTriggerLabel;  ///< labels for trigger 
  art::InputTag 	     fSimModuleLabel;      ///< name of detsim producer
  std::vector<art::InputTag> fPFParticleLabel; ///< labels for source of PFParticle
  bool                       fVerbose;       ///< print information about what's going on
  bool			     fIsData;        ///< switch for if this is data or MC

  CRTT0MatchAlg t0Alg;

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


icarus::CRTTPCTruthEff::CRTTPCTruthEff(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  , bt(p.get<fhicl::ParameterSet>("CRTBackTrack"))
  , t0Alg(p.get<fhicl::ParameterSet>("t0Alg"))
  , fCrtutils(new icarus::crt::CRTCommonUtils())
  
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  reconfigure(p);
}

void icarus::CRTTPCTruthEff::reconfigure(fhicl::ParameterSet const & p)
{
  fCRTHitLabel        =  p.get<art::InputTag> ("CRTHitLabel", "crthit");
  fTPCTrackLabel      =  p.get< std::vector<art::InputTag> >("TPCTrackLabel",             {""});
  fTriggerLabel       =  p.get<art::InputTag>("TriggerLabel","daqTrigger");
  fVerbose            =  p.get<bool>("Verbose");
  fIsData	      =  p.get<bool>("IsData");
  fSimModuleLabel     =  p.get<art::InputTag> ("SimModuleLabel","largeant");
  fPFParticleLabel    =  p.get< std::vector<art::InputTag> >("PFParticleLabel",             {""});
} // CRTT0Matching::reconfigure()

void icarus::CRTTPCTruthEff::beginJob()
{
  for(int i = 30; i < 50 + 1; i++){
    std::string tagger = "All";
    if (i>=35 && i<40) continue;
    if (i==48 || i==49) continue;
    // if(i < ){
    tagger = fCrtutils->GetRegionNameFromNum(i);//fCrtGeo.GetTagger(i).name;
  }

  trevtdisp= tfs->make<TTree>("trevtdisp","trevtdisp");

  trevtdisp->Branch("fEvent",&fEvent,"fEvent/I");
  trevtdisp->Branch("crtRegion",&fCrtRegion,"crtRegion/I");
  trevtdisp->Branch("fRun",&fRun,"fRun/I");
  trevtdisp->Branch("fSubRun",&fSubRun,"fSubRun/I");
  trevtdisp->Branch("track_trueID",&track_trueID,"track_trueID/I");
  trevtdisp->Branch("track_mother",&track_mother,"track_mother/I");
  trevtdisp->Branch("crt_mother",&crt_mother,"crt_mother/I");
  trevtdisp->Branch("track_ancestor",&track_ancestor,"track_ancestor/I");
  trevtdisp->Branch("crt_ancestor",&crt_ancestor,"crt_ancestor/I");
  trevtdisp->Branch("track_motherlayers",&track_motherlayers,"track_motherlayers/I");
  trevtdisp->Branch("crt_motherlayers",&crt_motherlayers,"crt_motherlayers/I");
  trevtdisp->Branch("crt_trueID",&crt_trueID,"crt_trueID/I");
  trevtdisp->Branch("crt_pdg",&crt_pdg,"crt_pdg/I");
  trevtdisp->Branch("trk_pdg",&trk_pdg,"trk_pdg/I");
  trevtdisp->Branch("ttl_tpctrks",&ttl_tpctrks,"ttl_tpctrks/I");
  trevtdisp->Branch("ttl_crthits",&ttl_crthits,"ttl_crthits/I");
  trevtdisp->Branch("this_dca",&this_dca,"this_dca/D");

  //*/
}

void icarus::CRTTPCTruthEff::analyze(const art::Event& event)
{

  fEvent  = event.id().event();
  fRun    = event.run();
  fSubRun = event.subRun();


    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    bt.Initialize(event);
    std::map<int, simb::MCParticle> particles;

    // Loop over the true particles
//    mcpart_ids.clear();
    for (auto const& particle: (*particleHandle)){
      
      // Make map with ID
      int partID = particle.TrackId();
      particles[partID] = particle;
    }//end loop over particles in particleHandle
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
      mf::LogError("CRTTPCTruthEff:") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
    }//end else
  }//end if( !fTriggerLabel.empty() )
  else {
    mf::LogError("CRTTPCTruthEff:") << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
  }//end else


  // Get CRT hits from the event
  art::Handle< std::vector<sbn::crt::CRTHit>> crtHitHandle;
  std::vector<art::Ptr<sbn::crt::CRTHit> > crtHitList;
  if (event.getByLabel(fCRTHitLabel, crtHitHandle))
    art::fill_ptr_vector(crtHitList, crtHitHandle);

  std::vector<sbn::crt::CRTHit> crtHits;
  int hit_i = 0;
  double minHitTime = 99999;
  double maxHitTime = -99999;

  for(auto const& hit : (*crtHitHandle)){
    double hitTime = double((1600 - hit.ts0_ns)/1e3);
    //double hitTime = (double)(int)hit.ts0_ns * 1e-3;
    if(hitTime < minHitTime) minHitTime = hitTime;
    if(hitTime > maxHitTime) maxHitTime = hitTime;

    crtHits.push_back(hit);
    hit_i++;
  }

  int numtpctrks=(int)fTPCTrackLabel.size(), numcrthits=(int)crtHits.size();

  //Temp for debugging MC
  std::cout << "# of TPC tracks:\t" << numtpctrks;
  std::cout << "# of CRT Hits:\t" << numcrthits;

  ttl_tpctrks = numtpctrks; ttl_crthits = numcrthits;

  //Begin building association between TPC tracks in the event and their T0s (if they have them)
  art::Handle<std::vector<recob::Track>> pandoratrkHandle;
  std::vector<art::Ptr<recob::Track>> pandoratrks;
  if(event.getByLabel("pandoraTrackGausCryoW",pandoratrkHandle)){ art::fill_ptr_vector(pandoratrks,pandoratrkHandle); }//end if(event.getByLabel("pandoraTrackGausCryoW",pandoratrkHandle)

  //----------------------------------------------------------------------------------------------------------
  //                                DISTANCE OF CLOSEST APPROACH ANALYSIS
  //----------------------------------------------------------------------------------------------------------cc

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

  // if(fVerbose) std::cout<<"----------------- DCA Analysis -------------------"<<std::endl;

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

      //Get PFParticles
      auto pfpListHandle = event.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel[it]);
      if (!pfpListHandle.isValid()) continue;

      // Get reconstructed tracks from the event
      //Get PFParticle-Track association
      art::FindManyP<recob::PFParticle> fmpfp(tpcTrackHandle, event, trackLabel);

      //Get T0-PFParticle association
      art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, event, fPFParticleLabel[it]);

      // Loop over reconstructed tracks
      for (auto const& tpcTrack : (*tpcTrackHandle)){

	this_dca = DBL_MAX; crt_trueID = INT_MAX; track_trueID = INT_MAX; track_mother = INT_MAX;

        //Find PFParticle for track i
        //art::Ptr::key() gives the index in the vector
        auto pfps = fmpfp.at(tpcTrack.ID());

	std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());

	//Construct histogram for lookint at *all* hits
	//Using the vector of art Ptrs to recob::Hit, we can use RecoUtils.* to get the likely truth ID:
	track_trueID = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData,hits,false);

	track_mother = particles[track_trueID].Mother();
	track_ancestor = track_mother;


	std::cout << "New TPC track in the handle!\n";
	icarus::matchCand closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, hits, crtHits, m_trigger_timestamp, fIsData);
	sbn::crt::CRTHit checkhit = closest.thishit;
	fCrtRegion=(int)fCrtutils->AuxDetRegionNameToNum(closest.thishit.tagger);

	double tempdca =(double) closest.dca;
	this_dca = tempdca; std::cout << "closest.dca: " << this_dca << std::endl;
	crt_trueID =  bt.TrueIdFromTotalEnergy(event, checkhit);// std::cout << temp_crt_trueID << std::endl;
	crt_mother = particles[crt_trueID].Mother();
	crt_ancestor = crt_mother;

	int tmp_counter = 0;
	while(track_ancestor%10000000!=0&&tmp_counter<1000){

		track_ancestor = particles[track_ancestor].Mother();
		std::cout << "loop " << tmp_counter << std::endl;
		tmp_counter++;

	}//end while(track_mother%10000000!=0||tmp_counter<100)
	track_motherlayers = tmp_counter;
	tmp_counter=0;
	while(crt_ancestor%10000000!=0&&tmp_counter<1000){

		crt_ancestor = particles[crt_ancestor].Mother();
		std::cout << "loop " << tmp_counter << std::endl;
		tmp_counter++;

	}//end while(crt_mother%10000000!=0||tmp_counter<100)
	crt_motherlayers = tmp_counter;

	crt_pdg = particles[crt_trueID].PdgCode();
	trk_pdg = particles[track_trueID].PdgCode();

	trevtdisp->Fill();
	
      }// track loops in each tracklebel
    }// end loop over TPC tracks
}
void icarus::CRTTPCTruthEff::endJob()
{

} // CRTT0Matching::endJob()


DEFINE_ART_MODULE(icarus::CRTTPCTruthEff)
