/**
 * @file   icaruscode/CRT/CRTT0Tagging_module.cc
 * @author Francesco Poppi (poppi@bo.infn.it)
 * @date   October 2024
 */

#include "sbnobj/Common/CRT/CRTHit.hh"
//#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "icaruscode/CRT/CRTUtils/CRTMatchingUtils.h"
#include "icaruscode/IcarusObj/CRTTPCMatchingInfo.h"
#include "icaruscode/CRT/CRTUtils/RecoUtils.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>

// LArSoft
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/CoreUtils/zip.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

// ROOT
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TTree.h"

namespace icarus::crt {

class CRTT0Tagging;

}  // namespace icarus::crt

using namespace icarus::crt;

class icarus::crt::CRTT0Tagging : public art::EDProducer {
public:
  
  using CRTHit = sbn::crt::CRTHit;
  //using CRTPMTMatching = sbn::crt::CRTPMTMatching;

  explicit CRTT0Tagging(fhicl::ParameterSet const& p);
  
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTT0Tagging(CRTT0Tagging const&) = delete;
  CRTT0Tagging(CRTT0Tagging&&) = delete;
  CRTT0Tagging& operator=(CRTT0Tagging const&) = delete;
  CRTT0Tagging& operator=(CRTT0Tagging&&) = delete;

  // Required functions.
  bool hasModID(std::uint32_t modID, sbn::crt::CRTHit const& crthit);
  //void beginRun(art::Run& r) override;
  void beginJob() override;
  void produce(art::Event& e) override;
  void endJob() override;

private:

  // Declare member data here.

  // MC Truth
  art::InputTag fSimulationProducerLabel; 
  art::InputTag fAuxDetSimProducerLabel;
  //art::InputTag fCRTSimHitProducerLabel; 
  //art::InputTag fCRTTrueHitProducerLabel;
  //art::InputTag fCRTDetSimProducerLabel;
  //art::InputTag fCRTSimTrackProducerLabel;
  art::InputTag fSimChannelProducerLabel;

  art::InputTag fCrtHitModuleLabel;
  
  std::vector<art::InputTag> fTPCTrackLabel; ///< labels for source of tracks
  std::vector<art::InputTag> fPFParticleLabel; ///< labels for source of PFParticle
  std::vector<art::InputTag> fHitLabel; ///< labels for source of hits
  std::vector<art::InputTag> fTRKHMLabel; ///< labels for hit metadata
  
  art::ServiceHandle<art::TFileService> tfs;

  CRTCommonUtils fCrtUtils; 
  CRTMatchingAlg fMatchingAlg;
  geo::GeometryCore const* fGeometryService;  ///< pointer to Geometry provider
  
  double fMinimalTrackLength;
  int fMinimumGoodHits;
  double fMaximalCRTDistance;
  double fGoodCandidateDistance;
  double fMaximumDeltaT;
  bool fData;
  bool fSkipTruth;

  icarus::crt::TopCRTCentersMap fTopCRTCenterMap;
  icarus::crt::TopCRTTransformations fTopCRTTransformations;

  TTree* fTree;
  int fEvent;        ///< number of the event being processed
  int fRun;          ///< number of the run being processed
  int fSubRun;       ///< number of the sub-run being processed
  int fCrtRegion;
  int fCrtSys;
  int fCryo;
  double fTrackLength;
  double fTrackCrtDistance;
  double fTrackCrtDeltaX;
  double fTrackCrtDeltaY;
  double fTrackCrtDeltaZ;
  double fCrtX;
  double fCrtY;
  double fCrtZ;
  double fCrossPointX;
  double fCrossPointY;
  double fCrossPointZ;      
  bool fTrueMatch;

};

icarus::crt::CRTT0Tagging::CRTT0Tagging(fhicl::ParameterSet const& p)
    : EDProducer{p},
      fSimulationProducerLabel(p.get<art::InputTag>("SimulationLabel","largeant")),
      fAuxDetSimProducerLabel(p.get<art::InputTag>("AuxDetSimProducerLabel","genericcrt")),
      fSimChannelProducerLabel(p.get<art::InputTag>("SimChannelProducer", {"daq:simpleSC"})),
      fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit")),
      fTPCTrackLabel(p.get< std::vector<art::InputTag> >("TPCTrackLabel",             {""})),
      fPFParticleLabel(p.get< std::vector<art::InputTag> >("PFParticleLabel",             {""})),  
      fHitLabel(p.get< std::vector<art::InputTag> >("HitLabel",             {""})),
      fTRKHMLabel(p.get< std::vector<art::InputTag> > ("TRKHMLabel", {""})),
      fMatchingAlg(p.get<fhicl::ParameterSet> ("MatchingAlg")),
      fGeometryService(lar::providerFrom<geo::Geometry>()),
      fMinimalTrackLength(p.get<double>("MinimalTrackLength", 40.0)),
      fMinimumGoodHits(p.get<double>("MinimumGoodHits", 5)),
      fMaximalCRTDistance(p.get<double>("MaximalCRTDistance", 300.)),
      fGoodCandidateDistance(p.get<double>("GoodCandidateDistance", 100.)),
      fMaximumDeltaT(p.get<double>("MaximumDeltaT", 10000.)),
      fData(p.get<bool>("isData", true)),
      fSkipTruth(p.get<bool>("skipTruth", false))
{
  
  produces< std::vector<anab::T0>                   >();
  produces< art::Assns<recob::Track , anab::T0>     >();
  produces< art::Assns<sbn::crt::CRTHit, anab::T0>  >();  
  produces< std::vector<icarus::CRTTPCMatchingInfo> >();
  //produces< art::Assns<icarus::CRTTPCMatchingInfo, anab::T0>  >();  
  //produces< art::Assns<recob::Track, icarus::CRTTPCMatchingInfo>  >();
  //produces< art::Assns<sbn::crt::CRTHit, icarus::CRTTPCMatchingInfo>  >();

  if (fTPCTrackLabel.size() != fPFParticleLabel.size()) {
    throw art::Exception{ art::errors::Configuration }
      << fTPCTrackLabel.size() << " TPC track data products configured (`TPCTrackLabel`), should have been "
      << fPFParticleLabel.size();
  }
  if (fHitLabel.size() != fPFParticleLabel.size()) {
    throw art::Exception{ art::errors::Configuration }
      << fHitLabel.size() << " TPC hit data products configured (`HitLabel`), should have been "
      << fPFParticleLabel.size();
  }

  if (fTRKHMLabel.size() > fPFParticleLabel.size()) {
    throw art::Exception{ art::errors::Configuration }
      << fTRKHMLabel.size() << " track-hit metadata data products configured (`TRKHMproducer`), should have been "
      << fPFParticleLabel.size();
  }
  fTRKHMLabel.resize(fPFParticleLabel.size()); // extend with empty labels
  // replace empty defaults with the actual input tag value, assumed the same as tracks
  for (std::size_t i = 0; i < fTRKHMLabel.size(); ++i)
    if (fTRKHMLabel[i].empty()) fTRKHMLabel[i] = fTPCTrackLabel[i];

}

bool icarus::crt::CRTT0Tagging::hasModID(std::uint32_t modID, sbn::crt::CRTHit const& crthit){
  for(auto const& mactopes: crthit.pesmap){
    for(auto const& chanpe: mactopes.second){
      int thisModID=(int)fCrtUtils.MacToAuxDetID(mactopes.first, chanpe.first);
      if(thisModID==(int)modID) return true;
    }
  }
  return false;
}

//void icarus::crt::CRTT0Tagging::beginRun(art::Run& r)
void icarus::crt::CRTT0Tagging::beginJob()
{
  fTopCRTCenterMap=icarus::crt::LoadTopCRTCenters();
  fTopCRTTransformations=icarus::crt::LoadTopCRTTransformations();

  fTree = tfs->make<TTree>("matchTree","CRTHit - TPC track matching analysis");

  fTree->Branch("Event",              &fEvent,              "Event/I");
  fTree->Branch("SubRun",             &fSubRun,             "SubRun/I");
  fTree->Branch("Run",                &fRun,                "Run/I");
  fTree->Branch("Cryo",               &fCryo,               "Cryo/I");
  fTree->Branch("CrtSys",             &fCrtSys,             "CrtSys/I");
  fTree->Branch("CrtRegion",          &fCrtRegion,          "CrtRegion/I");
  fTree->Branch("TrackLength",        &fTrackLength        );
  fTree->Branch("TrackCrtDistance",   &fTrackCrtDistance   );
  fTree->Branch("TrackCrtDeltaX",     &fTrackCrtDeltaX     );
  fTree->Branch("TrackCrtDeltaY",     &fTrackCrtDeltaY     );
  fTree->Branch("TrackCrtDeltaZ",     &fTrackCrtDeltaZ     );
  fTree->Branch("CrtX",               &fCrtX               );
  fTree->Branch("CrtY",               &fCrtY               );
  fTree->Branch("CrtZ",               &fCrtZ               );
  fTree->Branch("CrossPointX",        &fCrossPointX        );
  fTree->Branch("CrossPointY",        &fCrossPointY        );
  fTree->Branch("CrossPointZ",        &fCrossPointZ        );
  fTree->Branch("TrueMatch",          &fTrueMatch          );
}

void icarus::crt::CRTT0Tagging::produce(art::Event& e)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);

  mf::LogDebug("CRTT0Tagging: ") << "beginning production" << '\n';

  auto t0col = std::make_unique< std::vector<anab::T0> > ();
  auto trackAssn = std::make_unique< art::Assns<recob::Track, anab::T0> >();
  auto t0CrtHitAssn = std::make_unique< art::Assns<sbn::crt::CRTHit, anab::T0> >();
  auto matchInfoCol = std::make_unique< std::vector<icarus::CRTTPCMatchingInfo> >();
  //auto t0matchInfoAssn = std::make_unique< art::Assns<icarus::CRTTPCMatchingInfo, anab::T0> >();
  //auto trackMatchInfoAssn = std::make_unique< art::Assns<recob::Track, icarus::CRTTPCMatchingInfo> >();
  //auto matchInfoCrtHitAssn = std::make_unique< art::Assns<sbn::crt::CRTHit, icarus::CRTTPCMatchingInfo> >();

  std::map< int, const simb::MCParticle*> particleMap;
  std::map<std::pair<int,double>,std::vector<int>> crtParticleMap;
  std::map<int,bool> isNuMap;
  std::vector<art::Ptr<sim::SimChannel>> simchannels;

  // CRTHits
  std::vector<art::Ptr<CRTHit>> CRTHitList;
  art::ValidHandle<std::vector<CRTHit>> crthits = e.getValidHandle<std::vector<CRTHit>>(fCrtHitModuleLabel);
  art::fill_ptr_vector(CRTHitList, crthits);

  // If it is not data is MC.
  // Retrieving MC truth information.
  // Three maps (particleMap, crtParticleMap and isNuMap) are filled.
  // -> particleMap maps truth level informations for all the MC particles. The key is the Geant4 TrackID
  // -> crtParticleMap maps the CRTHits at truth level. The key is a pair of reconstructed CRT Hit module and time.
  // The object of the map is a vector of Geant4 Track IDs of the particles that generated that hit.
  // -> isNuMap maps the reconstructed tracks. The key is the Track ID and the object is a boolean True/False if
  // the Track was a neutrino related interaction or not.
  if(!fData){
    if(fSkipTruth){
      mf::LogInfo("CRTT0Tagging") <<"This is MC, but MC truth is not considered!";
    } else{
      mf::LogInfo("CRTT0Tagging") <<"This is MC, MC truth is considered!";
    }   
  }
  if(!fData && !fSkipTruth){
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;

    art::Handle< std::vector<sim::SimChannel>> simChannelHandle;
    if (!e.getByLabel(fSimChannelProducerLabel, simChannelHandle)){
	    throw cet::exception("CRTT0Tagging") 
	        << " No sim::SimChannel objects in this event - "
	        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    art::fill_ptr_vector(simchannels, simChannelHandle);

    // Define "handle" to Generator level MCTruth objects
    art::Handle< vector<simb::MCTruth>> genHandle;
    // Define a "handle" to point to a vector of MCParticle objects.
    art::Handle< vector<simb::MCParticle> > particleHandle;

    if (!e.getByLabel("generator", genHandle)) {
        std::cout << "could not get handle to gen objects!!!" << std::endl;
    }

    if (!e.getByLabel(fSimulationProducerLabel, particleHandle)) {
	    // If we have no MCParticles at all in an event, but we are requiring
      // to have MC truth information, throw exception.
	    throw cet::exception("CRTT0Tagging") 
	      << " No simb::MCParticle objects in this event - "
	      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // Handle to AuxDetSimChannel (CRT module) objects generated by LArG4
    art::Handle<vector<sim::AuxDetSimChannel> > auxDetSimChannelHandle;
    if (!e.getByLabel(fAuxDetSimProducerLabel, auxDetSimChannelHandle)) {
      throw cet::exception("CRTT0Tagging")
        << " No sim::AuxDetSimChannel objects in this event - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    //if((*genHandle).size()>1) 
    //      throw cet::exception("CRTT0Tagging") << "gen stage MCParticle vector has more than 1 entry!" << '\n';
    for ( auto const& particle : (*particleHandle) ){
      // Add the address of the MCParticle to the map, with the
      // track ID as the key.
      particleMap[particle.TrackId()] = &particle;
      art::Ptr<simb::MCTruth> mcTruth=partInventory->ParticleToMCTruth_P(&particle);
      
      bool isNu=false;
      
      if (mcTruth->Origin() == simb::kBeamNeutrino) isNu=true;

      isNuMap[particle.TrackId()] = isNu;

      for ( auto const& channel : (*auxDetSimChannelHandle) ){
        auto const& auxDetIDEs = channel.AuxDetIDEs();
        for ( auto const& ide : auxDetIDEs ){
          if ( ide.trackID != particle.TrackId() ) continue;
          if ( ide.energyDeposited * 1.0e6 < 50 ) continue; // skip energy deposits of less then 50 keV
          size_t adid = channel.AuxDetID();
          uint32_t region=fCrtUtils.AuxDetRegionNameToNum(fCrtUtils.GetAuxDetRegion(adid));
          uint32_t modID=channel.AuxDetID();
          float aveT = (ide.entryT + ide.exitT) / 2.0;
          for(auto const& crthit : CRTHitList){
            if(crthit->plane!=(int)region) continue;
            if(abs(aveT-crthit->ts1_ns)>200) continue;
            bool modFound=hasModID(modID, *crthit);
            if(!modFound)continue;
            std::pair<int,double> thisMatch = std::make_pair((int)crthit->feb_id[0],crthit->ts1_ns);
            crtParticleMap[thisMatch].push_back(particle.TrackId());
          } // CRT Hits loop
        } // Energy deposits loop
      } // CRT sim channels loop
    } // MC particles loop
  } // End MC Only
  for(const auto& [ PFPLabel, TPCTrackLabel, HitLabel, TRKHMLabel ]: util::zip(fPFParticleLabel, fTPCTrackLabel, fHitLabel, fTRKHMLabel)){
    std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
    art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(PFPLabel);
    art::fill_ptr_vector(PFParticleList, pfparticles);

    // Pandora MetaData
    art::FindOne<anab::T0> fmt0pandora(pfparticles, e, PFPLabel);
    art::FindOne<larpandoraobj::PFParticleMetadata> PFPMetaDataAssoc(pfparticles, e, PFPLabel);

    // Tracks
    art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(TPCTrackLabel);

    // Track - associated data
    art::FindManyP<recob::Track> fmTracks(PFParticleList, e, TPCTrackLabel);

    // Collect all hits
    art::ValidHandle<std::vector<recob::Hit>> allhit_handle = e.getValidHandle<std::vector<recob::Hit>>(HitLabel);
    std::vector<art::Ptr<recob::Hit>> allHits;
    art::fill_ptr_vector(allHits, allhit_handle);

    std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map;
    std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map;
    if(!simchannels.empty() && !fData){
      art::ServiceHandle<cheat::BackTrackerService> btServ;
      // ID (TrackID) refers to the reconstructed TrackID, it is not the Geant4 ID.
      id_to_ide_map = RecoUtils::PrepSimChannels(simchannels, *fGeometryService);
      id_to_truehit_map = RecoUtils::buildTrackIDtoHitsMap(allHits, clockData, *btServ.get());
    }
    // Start looping on the particles
    for (art::Ptr<recob::PFParticle> const& p_pfp: PFParticleList) {
      const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(p_pfp.key());
      if (thisTrack.size() != 1) continue;    
      art::Ptr<recob::Track> trkPtr = thisTrack.at(0);
      const recob::Track &track = *trkPtr;
      if(track.Length()<fMinimalTrackLength) continue;
      art::InputTag thm_label = TRKHMLabel;
      art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e, thm_label);
      std::vector<art::Ptr<recob::Hit>> emptyHitVector;
      const std::vector<art::Ptr<recob::Hit>> &trkHits  = fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;
      std::vector<const recob::TrackHitMeta*> emptyTHMVector;
      const std::vector<const recob::TrackHitMeta*> &trkHitMetas = fmtrkHits.isValid() ? fmtrkHits.data(trkPtr.key()) : emptyTHMVector;
      int trueTrackId=-9;
      if(!fData) trueTrackId= abs(RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, trkHits, false));
      // T0
      double t0 = std::numeric_limits<float>::signaling_NaN();
      if( auto const& t0ref = fmt0pandora.at(p_pfp.key())) t0 = t0ref.ref().Time();

      int goodHits=0;
      int countE=0, countW=0;
      
      // These counters are used to determine if track is CC-E, EE, EW, CC-W, WE, WW
      // depending on the track type, the Top CRT uses the appropriate position corrections
      std::vector<double> hx, hy, hz, ht;
      for(auto const& [trkHit, trkHitMeta]: util::zip(trkHits, trkHitMetas)){
        bool badhit = (trkHitMeta->Index() == std::numeric_limits<unsigned int>::max()) ||
                    (!track.HasValidPoint(trkHitMeta->Index()));
        if(badhit) continue;
        geo::Point_t loc = track.LocationAtPoint(trkHitMeta->Index());
        hx.push_back(loc.X()); hy.push_back(loc.Y()); hz.push_back(loc.Z());
        ht.push_back(trkHit->PeakTime());
        goodHits++;
        if(trkHit->WireID().TPC==0 || trkHit->WireID().TPC==1) countE++;
        else countW++;
      }
      if(goodHits<fMinimumGoodHits) {
        mf::LogDebug("CRTT0Tagging:")<<"Track does not have the minimal requirements of good hits: "<<goodHits;
        continue;
      }
      int trackType=-1;
      int const cryo = trkHits[0]->WireID().Cryostat;
      if(countW!=0 && countE!=0 && cryo==0) trackType=0; //CCEast
      else if(countW!=0 && countE==0 && cryo==0) trackType=2; //East-West
      else if(countW==0 && countE!=0 && cryo==0) trackType=1; //East-East      
      else if(countW!=0 && countE!=0 && cryo==1) trackType=3; //CCWest
      else if(countW!=0 && countE==0 && cryo==1) trackType=5; //West-West
      else if(countW==0 && countE!=0 && cryo==1) trackType=4; //West-East     
      icarus::crt::TopCRTCorrectionMap TopCRTCorrection;
      if(trackType==0) TopCRTCorrection=fTopCRTTransformations.EastCC;
      else if(trackType==1) TopCRTCorrection=fTopCRTTransformations.EE;
      else if(trackType==2) TopCRTCorrection=fTopCRTTransformations.EW;
      else if(trackType==3) TopCRTCorrection=fTopCRTTransformations.WestCC;
      else if(trackType==4) TopCRTCorrection=fTopCRTTransformations.WE;
      else if(trackType==5) TopCRTCorrection=fTopCRTTransformations.WW;
      
      std::vector<icarus::crt::CandCRT> crtCands;
      for(art::Ptr<CRTHit> p_crthit: CRTHitList){
        const CRTHit &crtHit = *p_crthit;
        double crtTime=crtHit.ts1_ns/1e3;
        if(!isnan(t0)){
          if(fabs(t0-crtHit.ts1_ns)>fMaximumDeltaT) continue;
        }
        icarus::crt::DriftedTrack thisDriftedTrack = fMatchingAlg.DriftTrack(trkHits, trkHitMetas, fGeometryService, detProp, crtTime, track);    
        if(thisDriftedTrack.outbound>0) continue;
        icarus::crt::PCAResults driftedPCAResults=fMatchingAlg.PCAfit(thisDriftedTrack.spx, thisDriftedTrack.spy, thisDriftedTrack.spz);
        icarus::crt::TranslationVector translVector = {driftedPCAResults.eigenVector1, driftedPCAResults.mean};
        int crtSys=-1;
        if(crtHit.plane<=34) crtSys=0;
        else if (crtHit.plane==50) crtSys=2;
        else crtSys=1;
        if(crtSys==2) continue; // lets discard Bottom CRT Hits for the moment

        double deltaX=std::numeric_limits<float>::signaling_NaN();
        double deltaY=std::numeric_limits<float>::signaling_NaN();
        double deltaZ=std::numeric_limits<float>::signaling_NaN();
        double crtDistance=std::numeric_limits<float>::signaling_NaN();

        icarus::crt::CRTPlane thisCRTPlane = fMatchingAlg.DeterminePlane(crtHit);
        icarus::crt::CrossingPoint crossPoint = fMatchingAlg.DetermineProjection(translVector, thisCRTPlane);

        double crtX=crtHit.x_pos;
        double crtY=crtHit.y_pos;
        double crtZ=crtHit.z_pos;

        if(fData){ // Realignment only applies to Data, not MC
          if(crtSys==0){
            double centerDX, centerDY, centerDZ;
            centerDX=crtX-fTopCRTCenterMap[(int)crtHit.feb_id[0]].X;
            centerDY=crtY-fTopCRTCenterMap[(int)crtHit.feb_id[0]].Y;
            centerDZ=crtZ-fTopCRTCenterMap[(int)crtHit.feb_id[0]].Z;
            icarus::crt::AffineTrans thisAffine=TopCRTCorrection[(int)crtHit.feb_id[0]];
            std::pair<double,double> transCrt;
            if(crtHit.plane==30) {
              transCrt=icarus::crt::AffineTransformation(centerDX, centerDZ, thisAffine);
              crtX=transCrt.first;
              crtZ=transCrt.second;
            } else if(crtHit.plane==31 ||crtHit.plane==32) {
              transCrt=icarus::crt::AffineTransformation(centerDY, centerDZ, thisAffine);
              crtY=transCrt.first;
              crtZ=transCrt.second;            
            } else if(crtHit.plane==33 ||crtHit.plane==34){
              transCrt=icarus::crt::AffineTransformation(centerDX, centerDY, thisAffine);
              crtX=transCrt.first;
              crtY=transCrt.second;            
            }
          }
        }

        if(thisCRTPlane.first==0){
          deltaX=crtX-crossPoint.X();
          deltaY=0;
          deltaZ=crtZ-crossPoint.Z();
        } else if(thisCRTPlane.first==1){
          deltaX=0;
          deltaY=crtY-crossPoint.Y();
          deltaZ=crtZ-crossPoint.Z();
        } else if(thisCRTPlane.first==2){
          deltaX=crtX-crossPoint.X();
          deltaY=crtY-crossPoint.Y();          
          deltaZ=0;
        }
        crtDistance=sqrt(pow(deltaX,2)+pow(deltaZ,2)+pow(deltaY,2));
        if(crtDistance>fMaximalCRTDistance) continue;
        icarus::crt::CandCRT thisCrtCand={crtHit,p_crthit, thisCRTPlane.first, crtDistance, deltaX, deltaY, deltaZ, crossPoint};
        crtCands.push_back(thisCrtCand);
      } // End of CRT Hit loop
      if(crtCands.empty()) {
        mf::LogDebug("CRTT0Tagging:")<<"No Good CRT match candidates for this track";
        continue;
      }
      auto minElementIt = std::min_element(crtCands.begin(), crtCands.end(), [](const icarus::crt::CandCRT& a, const icarus::crt::CandCRT& b) {
        return a.distance < b.distance;
      });
      icarus::crt::CandCRT& bestCrtCand=*minElementIt;
      if(bestCrtCand.distance<=fGoodCandidateDistance){ 
        int matchedSys=-1;
        if(bestCrtCand.CRThit.plane<=34) matchedSys=0;
        else if (bestCrtCand.CRThit.plane==50) matchedSys=2;
        else matchedSys=1;
        if(matchedSys==2) continue; // lets discard Bottom CRT Hits for the moment
        bool trueMatch=false;
        if(!fData && !fSkipTruth){
          std::vector<int> crtTracks, crtPdgs;
          std::pair<int,double> thisMatch=std::make_pair((int)bestCrtCand.CRThit.feb_id[0], bestCrtCand.CRThit.ts1_ns);
          if(crtParticleMap.find(thisMatch) != crtParticleMap.end()){
            for(size_t j=0; j<crtParticleMap[thisMatch].size(); j++){
              crtTracks.push_back(crtParticleMap[thisMatch].at(j));
              crtPdgs.push_back(particleMap[crtParticleMap[thisMatch].at(j)]->PdgCode());
            }   
          }
          for(size_t thisID=0; thisID<crtTracks.size(); thisID++){
            if(crtTracks[thisID]==trueTrackId) trueMatch=true;
          }
        }

        fEvent=e.id().event();
        fRun=e.run();
        fRun=e.subRun();
        fCrtRegion=bestCrtCand.CRThit.plane;
        fCrtSys=matchedSys;
        fCryo=cryo;
        fTrackLength=track.Length();
        fTrackCrtDistance=bestCrtCand.distance;
        fTrackCrtDeltaX=bestCrtCand.deltaX;
        fTrackCrtDeltaY=bestCrtCand.deltaY;
        fTrackCrtDeltaZ=bestCrtCand.deltaZ;
        fCrtX=bestCrtCand.CRThit.x_pos;
        fCrtY=bestCrtCand.CRThit.y_pos;
        fCrtZ=bestCrtCand.CRThit.z_pos;
        fCrossPointX=bestCrtCand.crossPoint.X();
        fCrossPointY=bestCrtCand.crossPoint.Y();
        fCrossPointZ=bestCrtCand.crossPoint.Z();    
        fTrueMatch=trueMatch;
        fTree->Fill();
        mf::LogInfo("CRTT0Tagging")
	      <<"Matched CRT time = "<<bestCrtCand.CRThit.ts1_ns/1e3<<" [us] to track "<<track.ID()<<" with projection-hit distance = "<<bestCrtCand.distance<<" Track T0 "<<t0
	      <<"\nMatched CRT hit plane: "<<bestCrtCand.CRThit.plane<<" xpos "<<bestCrtCand.CRThit.x_pos<<" ypos "<<bestCrtCand.CRThit.y_pos<<" zpos "<<bestCrtCand.CRThit.z_pos;
        t0col->push_back(anab::T0(bestCrtCand.CRThit.ts1_ns, track.ID(), matchedSys, bestCrtCand.CRThit.plane,bestCrtCand.distance));
        util::CreateAssn(*this, e, *t0col, trkPtr, *trackAssn);
        util::CreateAssn(*this, e, *t0col, bestCrtCand.ptrCRThit, *t0CrtHitAssn);
        icarus::CRTTPCMatchingInfo matchInfo = {bestCrtCand.distance, matchedSys, bestCrtCand.CRThit.plane, bestCrtCand.CRThit.ts1_ns, bestCrtCand.deltaX, bestCrtCand.deltaY, bestCrtCand.deltaZ, bestCrtCand.crossPoint.X(), bestCrtCand.crossPoint.Y(), bestCrtCand.crossPoint.Z(), bestCrtCand.plane};        
        matchInfoCol->push_back(matchInfo);
        //util::CreateAssn(*this, e, *t0col, *matchInfoCol, *t0matchInfoAssn);
        // util::CreateAssn(*this, e, *matchInfoCol, trkPtr, *trackMatchInfoAssn);
        // util::CreateAssn(*this, e, *matchInfoCol, bestCrtCand.ptrCRThit, *matchInfoCrtHitAssn);
      }
	  } // End of Track Loop
	  
	} // End of Cryo Loop
  e.put(std::move(t0col));
  e.put(std::move(trackAssn));
  e.put(std::move(t0CrtHitAssn));
  e.put(std::move(matchInfoCol));
  //e.put(std::move(t0matchInfoAssn));
  //e.put(std::move(trackMatchInfoAssn));
  //e.put(std::move(matchInfoCrtHitAssn));
}

void CRTT0Tagging::endJob()
{

} // CRTT0Tagging::endJob()


DEFINE_ART_MODULE(CRTT0Tagging)
