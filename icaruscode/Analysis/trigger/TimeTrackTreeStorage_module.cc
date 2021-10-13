////////////////////////////////////////////////////////////////////////
// Class:       TimeTrackTreeStorage
// Plugin Type: analyzer (art v3_06_03)
// File:        TimeTrackTreeStorage_module.cc
//
// Generated at Tue Sep 21 10:33:10 2021 by Jacob Zettlemoyer using cetskelgen
// Authors: Jacob Zettlemoyer (FNAL, jzettle@fnal.gov)
//          Animesh Chatterjee (U. Pittsburgh, anc238@pitt.edu)
//          Gianluca Petrillo (SLAC, petrillo@slac.stanford.edu)
//
// Borrowed heavily from Gray Putnam's existing TrackCaloSkimmer
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#define MF_DEBUG

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TTree.h"
#include "TVector3.h"
#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include <vector>
#include <string>

namespace sbn {
  class TimeTrackTreeStorage;
  
  struct selTrackInfo {
    int trackID;
    double t0;
    float start_x;
    float start_y;
    float start_z;
    float end_x;
    float end_y;
    float end_z;
    float length;
    float dir_x;
    float dir_y;
    float dir_z;
    selTrackInfo():
      trackID(-1),
      t0(-1),
      start_x(-1),
      start_y(-1),
      start_z(-1),
      end_x(-1),
      end_y(-1),
      end_z(-1),
      length(-1),
      dir_x(-1),
      dir_y(-1),
      dir_z(-1)
    {}
  };
}

class sbn::TimeTrackTreeStorage : public art::EDAnalyzer {
public:
  explicit TimeTrackTreeStorage(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TimeTrackTreeStorage(TimeTrackTreeStorage const&) = delete;
  TimeTrackTreeStorage(TimeTrackTreeStorage&&) = delete;
  TimeTrackTreeStorage& operator=(TimeTrackTreeStorage const&) = delete;
  TimeTrackTreeStorage& operator=(TimeTrackTreeStorage&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag fPFPproducer;
  art::InputTag fT0Producer;
  art::InputTag fTrackProducer;
  art::InputTag fT0selProducer;

  std::vector<sbn::selTrackInfo> vTrackInfo;
  sbn::selTrackInfo fTrackInfo;
  
  std::vector<int> fTrackID;
  
  //std::string const fLogCategory;

  TTree *fStoreTree;

  //variables, maybe put in struct later after getting art parts correct
  unsigned int fEvent;
  unsigned int fRun;
  unsigned int fSubRun;

  int fTotalProcessed = 0;
  
};


sbn::TimeTrackTreeStorage::TimeTrackTreeStorage(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fPFPproducer = p.get< art::InputTag > ("PFPproducer", "pandoraGausCryoW");
  fT0Producer = p.get< art::InputTag > ("T0Producer", "pandoraGausCryoW");
  fT0selProducer = p.get< art::InputTag > ("T0selProducer", "pandoraGausCryoW");
  fTrackProducer = p.get< art::InputTag > ("TrackProducer", "pandoraTrackGausCryoW");
  //fLogCategory = p.get< std::string const > ("LogCategory", "");
  art::ServiceHandle<art::TFileService> tfs;
  fStoreTree = tfs->make<TTree>("TimedTrackStorage", "Timed Track Tree");
  fStoreTree->Branch("run", &fRun);
  fStoreTree->Branch("subrun", &fSubRun);
  fStoreTree->Branch("event", &fEvent);
  fStoreTree->Branch("trackID", &fTrackID);
  fStoreTree->Branch("selTracks", &vTrackInfo);
}

void sbn::TimeTrackTreeStorage::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  unsigned int run = e.run();
  unsigned int subrun = e.subRun();
  unsigned int event = e.event();
  fTrackID.clear();
  //std::cout << "Here!" << std::endl;
  if(vTrackInfo.size() > 0)
    vTrackInfo.clear();
  //std::cout << "Here!" << std::endl;
  fEvent = event;
  fSubRun = subrun;
  fRun = run;

  std::vector<art::Ptr<recob::PFParticle>> const& pfparticles = e.getByLabel<std::vector<art::Ptr<recob::PFParticle>>> (fT0selProducer);
  if(pfparticles.size() == 0)
    return;

  //std::vector<art::Ptr<recob::PFParticle>> PFParticles;
  //try {
  //art::ValidHandle<std::vector<art::Ptr<recob::PFParticle>>> pfparticles = e.getValidHandle<std::vector<art::Ptr<recob::PFParticle>>>(fT0selProducer);
    //art::fill_ptr_vector(PFParticles, pfparticles);
    //}
    //catch(...) {
    //std::cout << "PFP's with tag: " << fT0selProducer << " not present.\n";
    //throw;
    //}
  //PFParticles = pfparticles;

  /*
      for (art::InputTag const& inputTag: fTrackTimeTags) {

    auto const& T0toTrack
    = event.getByLabel<art::Assns<anab::T0, recob::PFParticle>>(inputTag);
    unsigned int const newTracks = selectTracks(T0toTrack, selectedTracks);

    mf::LogTrace(fLogCategory)
    << "From '" << inputTag.encode() << "': "
    << newTracks << " tracks selected"
      ;

  } // for 
  */

  //art::FindManyP<anab::T0> T0List(pfparticles, e, fT0selProducer);
  //art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTrackProducer);
  //auto assns = std::make_unique<art::Assns<art::Ptr<recob::PFParticle>, recob::Track>>;
  //auto const& PFPtoTrack = e.getByLabel<art::Assns<art::Ptr<recob::PFParticle>, recob::Track>> (fT0selProducer);
  //art::FindManyP<art::Ptr<recob::Track>> particleTracks(pfparticles, e, fTrackProducer);
  std::cout << "HERE!" << std::endl;
  art::FindOneP<recob::Track> particleTracks (pfparticles,e,fTrackProducer);
  std::cout << "PFParticles size: " << pfparticles.size() << " art::FindOneP Tracks Size: " << particleTracks.size() << std::endl;
  int processed = 0;
  for(unsigned int iPart = 0; iPart < pfparticles.size(); ++iPart)
  {
    art::Ptr<recob::PFParticle> particlePtr = pfparticles[iPart];
    //std::cout << particlePtr.key() << std::endl;
    art::Ptr<recob::Track> trackPtr = particleTracks.at(iPart); 
    if(!(trackPtr.isNull())) 
    {
      //std::cout << "Track Pointer: " << trackPtr << std::endl;
      std::cout << "PFP Pointer: " << particlePtr << std::endl;
      fTrackID.push_back(trackPtr->ID());
      
      fTrackInfo.trackID = trackPtr->ID();
      fTrackInfo.start_x = trackPtr->Start().X();
      fTrackInfo.start_y = trackPtr->Start().Y();
      fTrackInfo.start_z = trackPtr->Start().Z();
      fTrackInfo.end_x = trackPtr->End().X();
      fTrackInfo.end_y = trackPtr->End().Y();
      fTrackInfo.end_z = trackPtr->End().Z();
      fTrackInfo.dir_x = trackPtr->StartDirection().X();
      fTrackInfo.dir_y = trackPtr->StartDirection().Y();
      fTrackInfo.dir_z = trackPtr->StartDirection().Z();
      fTrackInfo.length = trackPtr->Length();
      vTrackInfo.push_back(fTrackInfo);
      
      ++processed;
      ++fTotalProcessed;
    }
  }
  std::cout << "Particles Processed: " << processed << std::endl;
  std::cout << "Total Particles Processed: " << fTotalProcessed << std::endl; 
  /*
  for(art::Ptr<recob::PFParticle> this_pfp: pfparticles) 
  {
    const std::vector<art::Ptr<recob::Track>> track = sel_tracks.at(this_pfp.key());
    if(track.size() != 1)
      continue;

    art::Ptr<recob::Track> trk = track.at(0);
    fTrackID.push_back(trk->ID());
    std::cout << "This track ID is: " << trk->ID() << std::endl; 
  }
  */
  fStoreTree->Fill();

}

DEFINE_ART_MODULE(sbn::TimeTrackTreeStorage)
