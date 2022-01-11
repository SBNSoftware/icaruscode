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
#include "icaruscode/Analysis/trigger/details/TriggerResponseManager.h"
#include "icaruscode/Decode/DataProducts/ExtraTriggerInfo.h"

// LArSoft libraries
// #include "lardata/DetectorInfoServices/DetectorClocksService.h"
// #include "larcore/Geometry/Geometry.h"
// #include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
// #include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Persistency/Common/FindOneP.h"
// #include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT libraries
#include "TTree.h"

// C/C++ libraries
#include <vector>
#include <string>



// -----------------------------------------------------------------------------
namespace sbn { class TimeTrackTreeStorage; }
/**
 * @brief Fills a ROOT tree with track-based triggering information.
 * 
 * 
 * Trigger information
 * ====================
 * 
 * There are two distinct groups of information about triggers in the tree.
 * The first one is the hardware trigger, that is the trigger that was actually
 * used to decide to record the event. That is the same information for all
 * tracks within the same event (and the same in both cryostats, too).
 * The other group of information is from trigger logic simulation on the
 * collected data (optical detector waveforms): there may be multiple branches
 * in this group, one for each simulated trigger logic, and each tree entry,
 * corresponding to a selected track, may have its own value for each trigger
 * logic response.
 * 
 * 
 * Simulated trigger information
 * ------------------------------
 * 
 * A branch is added for each configured trigger logic.
 * This module actually ignores the details of the logic yielding the results:
 * a trigger result is a group of data products under the same tag.
 * For example, a tag `M1` (presumably, a very loose trigger logic requiring
 * just one PMT pair above threshold anywhere in the detector, with no hint of
 * which that threshold is) will produce a branch with name `"M1"` including
 * information from the `M1` data products (so far, only a collection of
 * `raw::Trigger` is read).
 * It is assumed that there is one trigger result for each selected track
 * (as found in the data product from `T0selProducer` configuration parameter).
 * Note that a trigger object is expected to be present _regardless whether
 * the trigger fired or not_. It is assumed that the trigger fired if at least
 * one of the `raw::Trigger::TriggerBits()` is set, not fired otherwise.
 * 
 * Each branch mirrors a data record, `TriggerInfo_t`, reporting the response
 * for one simulated trigger logic.
 * The structure and interpretation of the trigger response branch is described
 * in ROOT TTree code by the string in
 * `TriggerInfo_t::TriggerResponseBranchStructure`, and its meaning is:
 * * `fired` (boolean): whether this trigger has triggered at all during the
 *   considered trigger gate.
 * * `time` (double): time of the trigger, taken directly from
 *   `raw::Trigger::TriggerTime()`; it is expected to be measurent in
 *   microseconds and relative to the nominal time reference of the event,
 *   which for data is the hardware trigger time and for simulation is the
 *   "hardware" beam gate time (however that is defined).
 * * `gate` (double): start time of the gate where the trigger logic was
 *   simulated. It is in the same time scale as the trigger time, and
 *   directly taken from `raw::Trigger::BeamGateTime()`.
 * 
 * 
 */
class sbn::TimeTrackTreeStorage : public art::EDAnalyzer {
  
  using TriggerInputSpec_t
    = sbn::details::TriggerResponseManager::TriggerInputSpec_t;
  
public:
  
  /// Data record the trigger response branch is based on.
  using TriggerInfo_t = sbn::details::TriggerResponseManager::TriggerInfo_t;
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    /// Information on a single trigger source for the tree.
    struct TriggerSpecConfig {
      
      fhicl::Atom<std::string> Name {
        fhicl::Name("Name"),
        Comment("name of the trigger (e.g. `\"M5O3\"`)")
        };
      fhicl::Atom<art::InputTag> TriggerTag {
        fhicl::Name("TriggerTag"),
        Comment("tag of the input trigger info data product")
        };
      
    }; // TriggerSpecConfig
    using TriggerSpecConfigTable
      = fhicl::TableAs<TriggerInputSpec_t, TriggerSpecConfig>;
    
    
    fhicl::Atom<art::InputTag> PFPproducer {
      Name("PFPproducer"),
      Comment("tag of the input particle flow objects to process")
      // mandatory
      };
    
    fhicl::OptionalAtom<art::InputTag> T0Producer {
      Name("T0Producer"),
      Comment("tag of the input track time (t0) information [as PFPproducer]")
      // default: as PFPproducer
      };
    
    fhicl::OptionalAtom<art::InputTag> T0selProducer {
      Name("T0selProducer"),
      Comment
        ("tag of the selected tracks (as a collection of art::Ptr) [as PFPproducer]")
      // default: as PFPproducer
      };
    
    fhicl::Atom<art::InputTag> TrackProducer {
      Name("TrackProducer"),
      Comment("tag of the association of particles to tracks")
      // mandatory
      };
    
    fhicl::Atom<art::InputTag> BeamGateProducer {
      Name("BeamGateProducer"),
      Comment("tag of beam gate information")
      // mandatory
      };
    
    fhicl::Atom<art::InputTag> TriggerProducer {
      Name("TriggerProducer"),
      Comment("tag of hardware trigger information")
      // mandatory
      };
    
    fhicl::Atom<art::InputTag> FlashProducer {
      Name("FlashProducer"),
      Comment("tag of flash information")
      //mandatory
      };

    fhicl::Sequence<TriggerSpecConfigTable> EmulatedTriggers {
      Name("EmulatedTriggers"),
      Comment("the emulated triggers to include in the tree")
    };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("label for output messages of this module instance"),
      "TimeTrackTreeStorage" // default
      };
    
  }; // Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  explicit TimeTrackTreeStorage(Parameters const& p);

  void analyze(art::Event const& e) override;
  
  void endJob() override;

private:
  
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  art::InputTag const fPFPproducer;
  art::InputTag const fT0Producer;
  art::InputTag const fT0selProducer;
  art::InputTag const fTrackProducer;
  art::InputTag const fBeamGateProducer;
  art::InputTag const fTriggerProducer;
  art::InputTag const fFlashProducer;
  std::string const fLogCategory;
  
  // --- END ---- Configuration parameters -------------------------------------

  // --- BEGIN -- Tree buffers -------------------------------------------------
  
  unsigned int fEvent;
  unsigned int fRun;
  unsigned int fSubRun;
  
  //std::vector<sbn::selTrackInfo> fTrackInfo;
  sbn::selTrackInfo fTrackInfo; //change to one entry per track instead of per event 
  sbn::selBeamInfo fBeamInfo;
  sbn::selTriggerInfo fTriggerInfo;
  sbn::selLightInfo fFlashInfo;
  
  // --- END ---- Tree buffers -------------------------------------------------
  
  TTree *fStoreTree = nullptr;

  unsigned int fTotalProcessed = 0;
  
  
  // `convert()` needs to be a free function
  friend TriggerInputSpec_t convert(Config::TriggerSpecConfig const& config);
  
  
  // --- BEGIN -- Trigger response branches ------------------------------------
  
  ///< Manages filling of trigger result branches.
  sbn::details::TriggerResponseManager fTriggerResponses;
  
  // --- END ---- Trigger response branches ------------------------------------
  
}; // sbn::TimeTrackTreeStorage


// -----------------------------------------------------------------------------
namespace sbn {
  
  TimeTrackTreeStorage::TriggerInputSpec_t convert
    (TimeTrackTreeStorage::Config::TriggerSpecConfig const& config)
  {
    return {
        config.Name()       // name
      , config.TriggerTag() // inputTag
      };
  } // convert(sbn::TriggerSpecConfig)
  
} // namespace sbn


// -----------------------------------------------------------------------------
sbn::TimeTrackTreeStorage::TimeTrackTreeStorage(Parameters const& p)
  : EDAnalyzer{p}
  // configuration
  , fPFPproducer      { p().PFPproducer() }
  , fT0Producer       { p().T0Producer().value_or(fPFPproducer) }
  , fT0selProducer    { p().T0selProducer().value_or(fPFPproducer) }
  , fTrackProducer    { p().TrackProducer() }
  , fBeamGateProducer { p().BeamGateProducer() }
  , fTriggerProducer  { p().TriggerProducer() }
  , fFlashProducer    { p().FlashProducer() }
  , fLogCategory      { p().LogCategory() }
  // algorithms
  , fStoreTree {
      art::ServiceHandle<art::TFileService>()->make<TTree>
        ("TimedTrackStorage", "Timed Track Tree")
      }
  , fTriggerResponses
      { p().EmulatedTriggers(), consumesCollector(), *fStoreTree }
{
  
  //
  // declaration of input
  //
  
  // consumes<std::vector<recob::PFParticle>>(fPFPproducer); // not yet?
  consumes<std::vector<art::Ptr<recob::PFParticle>>>(fT0selProducer);
  consumes<sbn::ExtraTriggerInfo>(fTriggerProducer);
  consumes<std::vector<sim::BeamGateInfo>>(fBeamGateProducer);
  consumes<art::Assns<recob::PFParticle, recob::Track>>(fTrackProducer);
  consumes<art::Assns<recob::PFParticle, anab::T0>>(fT0Producer);
  consumes<recob::OpFlash>(fFlashProducer);

  //
  // tree population
  //
  fStoreTree->Branch("run", &fRun);
  fStoreTree->Branch("subrun", &fSubRun);
  fStoreTree->Branch("event", &fEvent);
  fStoreTree->Branch("beamInfo", &fBeamInfo);
  fStoreTree->Branch("triggerInfo", &fTriggerInfo);
  fStoreTree->Branch("selTracks", &fTrackInfo);
  fStoreTree->Branch("selFlashes", &fFlashInfo);
  
} // sbn::TimeTrackTreeStorage::TimeTrackTreeStorage()


void sbn::TimeTrackTreeStorage::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  unsigned int const run = e.run();
  unsigned int const subrun = e.subRun();
  unsigned int const event = e.event(); 
  
  fEvent = event;
  fSubRun = subrun;
  fRun = run;
  fBeamInfo = {};
  
  std::vector<art::Ptr<recob::PFParticle>> const& pfparticles = e.getProduct<std::vector<art::Ptr<recob::PFParticle>>> (fT0selProducer);
  if(pfparticles.empty()) {
    mf::LogDebug(fLogCategory) << "No particles in '" << fT0selProducer.encode() << "'.";
    return;
  }

  std::vector<sim::BeamGateInfo> const& beamgate = e.getProduct<std::vector<sim::BeamGateInfo>> (fBeamGateProducer);
  if(beamgate.empty())
    mf::LogWarning(fLogCategory) << "No Beam Gate Information!";
  if(beamgate.size() > 1)
    mf::LogWarning(fLogCategory) << "Event has multiple beam gate info labels! (maybe this changes later to be standard)";
  fBeamInfo.beamGateSimStart = beamgate[0].Start();
  fBeamInfo.beamGateDuration = beamgate[0].Width();
  fBeamInfo.beamGateType = beamgate[0].BeamType();

  sbn::ExtraTriggerInfo const &triggerinfo = e.getProduct<sbn::ExtraTriggerInfo> (fTriggerProducer);
  //fTriggerInfo.beamType = triggerinfo.sourceType;
  fTriggerInfo.triggerTime = triggerinfo.triggerTimestamp;
  fTriggerInfo.beamGateTime = triggerinfo.beamGateTimestamp;
  fTriggerInfo.triggerID = triggerinfo.triggerID;
  fTriggerInfo.gateID = triggerinfo.gateID;
  
  //mf::LogTrace(fLogCategory) << "HERE!";
  art::FindOneP<recob::Track> particleTracks(pfparticles,e,fTrackProducer);
  art::FindOneP<anab::T0> t0Tracks(pfparticles,e,fT0Producer);  
  std::vector<recob::OpFlash> const &particleFlashes = e.getProduct<std::vector<recob::OpFlash>>(fFlashProducer);
  //art::FindOneP<recob::SpacePoint> particleSPs(pfparticles, e, fT0selProducer);
  //mf::LogTrace(fLogCategory) << "PFParticles size: " << pfparticles.size() << " art::FindOneP Tracks Size: " << particleTracks.size();
  
  
  // get an extractor bound to this event
  sbn::details::TriggerResponseManager::Extractors triggerResponseExtractors
    = fTriggerResponses.extractorsFor(e);
  
  unsigned int processed = 0;
  for(unsigned int iPart = 0; iPart < pfparticles.size(); ++iPart)
  {
    //art::Ptr<recob::PFParticle> particlePtr = pfparticles[iPart];
    //mf::LogTrace(fLogCategory) << particlePtr.key();
    fFlashInfo = {};
    art::Ptr<recob::Track> const trackPtr = particleTracks.at(iPart);
    if(trackPtr.isNull()) continue;
    
    art::Ptr<anab::T0> const t0Ptr = t0Tracks.at(iPart);
    float const track_t0 = t0Ptr->Time();
    if(!particleFlashes.empty())
    {
      float min_flash_t0_diff = 999999.0; 
      for(unsigned int iFlash = 0; iFlash < particleFlashes.size(); ++iFlash)
      {
	recob::OpFlash const flashPtr = particleFlashes.at(iFlash);
	float const flash_pe = flashPtr.TotalPE();
	float const flash_time = flashPtr.Time();
	float flash_t0_diff = flash_time - track_t0/1e3;
	if(std::abs(flash_t0_diff) < min_flash_t0_diff)
	{ 
	  fFlashInfo.flash_id = iFlash;
	  fFlashInfo.sum_pe = flash_pe;
	  fFlashInfo.flash_time = flash_time;
	  fFlashInfo.diff_flash_t0 = flash_t0_diff;
	  min_flash_t0_diff = std::abs(flash_t0_diff);
	}
      }
    }
    sbn::selTrackInfo trackInfo;
    trackInfo.trackID = trackPtr->ID();
    trackInfo.t0 = track_t0/1e3; //is this in nanoseconds? Will convert to seconds so I can understand better
    //if(track_t0/1e3 < 10 && track_t0/1e3 > -10)
    //mf::LogTrace(fLogCategory) << track_t0/1e3 << " Run is: " << fRun << " SubRun is: " << fSubRun << " Event is: " << fEvent << " Track ID is: " << trackPtr->ID();
    trackInfo.start_x = trackPtr->Start().X();
    trackInfo.start_y = trackPtr->Start().Y();
    trackInfo.start_z = trackPtr->Start().Z();
    trackInfo.end_x = trackPtr->End().X();
    trackInfo.end_y = trackPtr->End().Y();
    trackInfo.end_z = trackPtr->End().Z();
    trackInfo.dir_x = trackPtr->StartDirection().X();
    trackInfo.dir_y = trackPtr->StartDirection().Y();
    trackInfo.dir_z = trackPtr->StartDirection().Z();
    trackInfo.length = trackPtr->Length();
    fTrackInfo = trackInfo;
    /*
    for(size_t trajp = 0; trajp < trackPtr->NumberTrajectoryPoints()-1; ++trajp)
    {
      TVector3 cur_point(trackPtr->TrajectoryPoint(traj_p).position.X(), trackPtr->TrajectoryPoint(traj_p).position.Y(), trackPtr->TrajectoryPoint(traj_p).position.Z());
      TVector3 next_point(trackPtr->TrajectoryPoint(traj_p+1).position.X(), trackPtr->TrajectoryPoint(traj_p+1).position.Y(), trackPtr->TrajectoryPoint(traj_p+1).position.Z());
      if(abs(cur_point.X()) < 170 && abs(next_point.X()) > 170)
        //interpolate to get cathode crossing point
        
    }
    */
    //fTrackInfo.push_back(trackInfo);
    
    triggerResponseExtractors.fetch(iPart); // TODO no check performed; need to find a way
  
    ++processed;
    ++fTotalProcessed;
    fStoreTree->Fill();
  } // for particle
  //mf::LogInfo(fLogCategory) << "Particles Processed: " << processed
  //mf::LogTrace(fLogCategory) << "Total Particles Processed: " << fTotalProcessed;
  //fStoreTree->Fill();
  //fTrackInfo.clear();

} // sbn::TimeTrackTreeStorage::analyze()


// -----------------------------------------------------------------------------
void sbn::TimeTrackTreeStorage::endJob() {
  
  mf::LogInfo(fLogCategory) << "Processed " << fTotalProcessed << " tracks.";
  
} // sbn::TimeTrackTreeStorage::endJob()


// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(sbn::TimeTrackTreeStorage)
