/**
 * @file    TimeTrackTreeStorage_module.cc
 * @authors Jacob Zettlemoyer (FNAL, jzettle@fnal.gov),
 *          Animesh Chatterjee (U. Pittsburgh, anc238@pitt.edu),
 *          Gianluca Petrillo (SLAC, petrillo@slac.stanford.edu)
 * @date    Tue Sep 21 10:33:10 2021
 * 
 * Borrowed heavily from Gray Putnam's existing TrackCaloSkimmer
 */

#define MF_DEBUG

// ICARUS libraries
#include "Objects/TrackTreeStoreObj.h"
#include "icaruscode/Decode/DataProducts/ExtraTriggerInfo.h"

// LArSoft libraries
// #include "lardata/DetectorInfoServices/DetectorClocksService.h"
// #include "larcore/Geometry/Geometry.h"
// #include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT libraries
#include "TTree.h"
#include "TVector3.h"

// C/C++ libraries
#include <iostream>
#include <vector>
#include <string>


namespace sbn {
  class TimeTrackTreeStorage;
}

class sbn::TimeTrackTreeStorage : public art::EDAnalyzer {
public:
  explicit TimeTrackTreeStorage(fhicl::ParameterSet const& p);

  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag fPFPproducer;
  art::InputTag fT0Producer;
  art::InputTag fTrackProducer;
  art::InputTag fT0selProducer;
  art::InputTag fBeamGateProducer;
  art::InputTag fTriggerProducer;

  std::vector<sbn::selTrackInfo> vTrackInfo;
  sbn::selTrackInfo fTrackInfo;
  sbn::selBeamInfo fBeamInfo;
  sbn::selTriggerInfo fTriggerInfo;
  
  //std::string const fLogCategory;

  TTree *fStoreTree;

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
  fBeamGateProducer = p.get< art::InputTag > ("BeamGateProducer", "daqTrigger"); 
  fTriggerProducer = p.get< art::InputTag > ("TriggerProducer", "daqTrigger");
  //fLogCategory = p.get< std::string const > ("LogCategory", "");
  art::ServiceHandle<art::TFileService> tfs;
  fStoreTree = tfs->make<TTree>("TimedTrackStorage", "Timed Track Tree");
  fStoreTree->Branch("run", &fRun);
  fStoreTree->Branch("subrun", &fSubRun);
  fStoreTree->Branch("event", &fEvent);
  fStoreTree->Branch("beamInfo", &fBeamInfo);
  fStoreTree->Branch("triggerInfo", &fTriggerInfo);
  fStoreTree->Branch("selTracks", &vTrackInfo);
}

void sbn::TimeTrackTreeStorage::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  unsigned int run = e.run();
  unsigned int subrun = e.subRun();
  unsigned int event = e.event(); 
  if(vTrackInfo.size() > 0)
    vTrackInfo.clear();
  fEvent = event;
  fSubRun = subrun;
  fRun = run;
  fBeamInfo = {};
  fTriggerInfo = {};
  
  std::vector<art::Ptr<recob::PFParticle>> const& pfparticles = e.getProduct<std::vector<art::Ptr<recob::PFParticle>>> (fT0selProducer);
  if(pfparticles.size() == 0)
    return;

  std::vector<sim::BeamGateInfo> const& beamgate = e.getProduct<std::vector<sim::BeamGateInfo>> (fBeamGateProducer);
  if(beamgate.size() == 0)
    std::cout << "No Beam Gate Information!" << std::endl;
  if(beamgate.size() > 1)
    std::cout << "Event has multiple beam gate info labels! (maybe this changes later to be standard)" << std::endl;
  fBeamInfo.beamGateSimStart = beamgate[0].Start();
  fBeamInfo.beamGateDuration = beamgate[0].Width();
  fBeamInfo.beamGateType = beamgate[0].BeamType();

  sbn::ExtraTriggerInfo const &triggerinfo = e.getProduct<sbn::ExtraTriggerInfo> (fTriggerProducer);
  //fTriggerInfo.beamType = triggerinfo.sourceType;
  fTriggerInfo.triggerTime = triggerinfo.triggerTimestamp;
  fTriggerInfo.beamGateTime = triggerinfo.beamGateTimestamp;
  fTriggerInfo.triggerID = triggerinfo.triggerID;
  fTriggerInfo.gateID = triggerinfo.gateID;
  //std::cout << "HERE!" << std::endl;
  art::FindOneP<recob::Track> particleTracks (pfparticles,e,fTrackProducer);
  art::FindOneP<anab::T0> t0Tracks(pfparticles,e,fT0Producer);
  //art::FindOneP<recob::SpacePoint> particleSPs(pfparticles, e, fT0selProducer);
  //std::cout << "PFParticles size: " << pfparticles.size() << " art::FindOneP Tracks Size: " << particleTracks.size() << std::endl;
  int processed = 0;
  for(unsigned int iPart = 0; iPart < pfparticles.size(); ++iPart)
  {
    //art::Ptr<recob::PFParticle> particlePtr = pfparticles[iPart];
    //std::cout << particlePtr.key() << std::endl;
    art::Ptr<recob::Track> trackPtr = particleTracks.at(iPart);
    art::Ptr<anab::T0> t0Ptr = t0Tracks.at(iPart);
    float track_t0 = -999.0;
    if(!(trackPtr.isNull())) 
    {
      track_t0 = t0Ptr->Time();
      //std::cout << "PFP Pointer: " << particlePtr << std::endl;
      fTrackInfo.trackID = trackPtr->ID();
      fTrackInfo.t0 = track_t0/1e3; //is this in nanoseconds? Will convert to seconds so I can understand better
      //if(track_t0/1e3 < 10 && track_t0/1e3 > -10)
      //std::cout << track_t0/1e3 << " Run is: " << fRun << " SubRun is: " << fSubRun << " Event is: " << fEvent << " Track ID is: " << trackPtr->ID() << std::endl;
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
      /*
      for(size_t trajp = 0; trajp < trackPtr->NumberTrajectoryPoints()-1; ++trajp)
      {
	TVector3 cur_point(trackPtr->TrajectoryPoint(traj_p).position.X(), trackPtr->TrajectoryPoint(traj_p).position.Y(), trackPtr->TrajectoryPoint(traj_p).position.Z());
	TVector3 next_point(trackPtr->TrajectoryPoint(traj_p+1).position.X(), trackPtr->TrajectoryPoint(traj_p+1).position.Y(), trackPtr->TrajectoryPoint(traj_p+1).position.Z());
	if(abs(cur_point.X()) < 170 && abs(next_point.X()) > 170)
	  //interpolate to get cathode crossing point
	  
      }	  
      */
      vTrackInfo.push_back(fTrackInfo);
      
      ++processed;
      ++fTotalProcessed;
    }
  }
  //std::cout << "Particles Processed: " << processed << std::endl;
  //std::cout << "Total Particles Processed: " << fTotalProcessed << std::endl;
  fStoreTree->Fill();

}

DEFINE_ART_MODULE(sbn::TimeTrackTreeStorage)
