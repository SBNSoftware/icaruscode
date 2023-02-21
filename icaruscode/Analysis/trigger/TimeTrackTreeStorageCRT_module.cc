/**
 * @file    TimeTrackTreeStorage_module.cc
 * @authors Jacob Zettlemoyer (FNAL, jzettle@fnal.gov),
 *          Animesh Chatterjee (U. Pittsburgh, anc238@pitt.edu),
 *          Gianluca Petrillo (SLAC, petrillo@slac.stanford.edu)
 * @date    Tue Sep 21 10:33:10 2021
 * 
 * Borrowed heavily from Gray Putnam's existing TrackCaloSkimmer
 */

// ICARUS libraries
#include "icaruscode/Analysis/trigger/details/TriggerResponseManager.h"
#include "icaruscode/Analysis/trigger/details/CathodeCrossingUtils.h"
#include "icaruscode/Analysis/trigger/Objects/TrackTreeStoreObj.h"
#include "icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h"
#include "icaruscode/Utilities/TrajectoryUtils.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// LArSoft libraries
// #include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"
// #include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
// #include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
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
#include <utility> // std::move(), std::pair
#include <cmath>


namespace {
  
  // ---------------------------------------------------------------------------
  /// Returns the sequence of `track` valid points (as geometry points).
  std::vector<geo::Point_t> extractTrajectory
    (recob::Track const& track, bool reverse = false)
  {
    std::vector<geo::Point_t> trackPath;
    std::size_t index = track.FirstValidPoint();
    while (index != recob::TrackTrajectory::InvalidIndex) {
      trackPath.push_back(track.LocationAtPoint(index));
      if (++index >= track.NPoints()) break;
      index = track.NextValidPoint(index);
    }
    if (reverse) std::reverse(trackPath.begin(), trackPath.end());
    return trackPath;
  } // extractTrajectory()
  
} // local namespace


// -----------------------------------------------------------------------------
namespace sbn { class TimeTrackTreeStorage; }
/**
 * @brief Fills a ROOT tree with track-based triggering information.
 * 
 * 
 * Track selection criteria
 * =========================
 * 
 * This module does not perform almost any selection: it processes all the
 * reconstructed particle flow objects (PFO) (typically from Pandora) selected
 * by the module in `T0selProducer`, which may perform the desired selection.
 * For every PFO which is associated to a track (`recob::Track`), that
 * associated track and the associated time (`anab::T0`) are saved in a tree
 * entry.
 * 
 * 
 * Optical detector information
 * -----------------------------
 * 
 * Currently, _all_ reconstructed flashes in the cryostat are saved together
 * with _each_ tree entry, i.e. to each selected track. There is no attempt to
 * match a single flash to a track.
 * 
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

    fhicl::Atom<art::InputTag> TrackFitterProducer {
      Name("TrackFitterProducer"),
        Comment("tag of the association of the tracks with the hits")
        // mandatory                                                                                    
        };

    fhicl::Atom<art::InputTag> CaloProducer {
      Name("CaloProducer"),
        Comment("tag of the association of the tracks with the calorimetry module")
        //mandatory
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

    fhicl::Atom<float> MODA {
      Name("MODA"),
      Comment("first recombination parameter for dE/dx calculations"),
      0.930 //default
      };
    
    fhicl::Atom<float> MODB {
      Name("MODB"),
      Comment("second recombination parameter for dE/dx calculations"),
      0.212
      };
    
    fhicl::Atom<float> Wion {
      Name("Wion"),
      Comment("work function for recombination"),
      0.0000236016
      };
    
    fhicl::Atom<float> Efield {
      Name("Efield"),
      Comment("Electric field in kV/cm"),
      0.5
      };
    
    fhicl::Atom<bool> ForceDowngoing {
      Name("ForceDowngoing"),
      Comment("force all tracks to be downgoing, flipping them when necessary"),
      false
      };
    
  }; // Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  explicit TimeTrackTreeStorage(Parameters const& p);

  sbn::selHitInfo makeHit(const recob::Hit &hit,
                          unsigned hkey,
                          const recob::Track &trk,
                          const recob::TrackHitMeta &thm,
                          const std::vector<art::Ptr<anab::Calorimetry>> &calo,
                          const geo::GeometryCore *geo);

  float dEdx_calc(float dQdx, 
                  float A,
                  float B,
                  float Wion,
                  float E);

  void analyze(art::Event const& e) override;
  
  void endJob() override;

private:
  
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  art::InputTag const fPFPproducer;
  art::InputTag const fT0Producer;
  art::InputTag const fT0selProducer;
  art::InputTag const fTrackProducer;
  art::InputTag const fTrackFitterProducer;
  art::InputTag const fCaloProducer;
  art::InputTag const fBeamGateProducer;
  art::InputTag const fTriggerProducer;
  art::InputTag const fFlashProducer;
  std::string const fLogCategory;
  float const fMODA;
  float const fMODB;
  float const fWion;
  float const fEfield;
  bool const fForceDowngoing; ///< Whether to force all tracks to be downgoing.
  
  // --- END ---- Configuration parameters -------------------------------------

  // --- BEGIN -- Tree buffers -------------------------------------------------
  
  unsigned int fEvent;
  unsigned int fRun;
  unsigned int fSubRun;
  
  sbn::selTrackInfo fTrackInfo; //change to one entry per track instead of per event 
  sbn::selBeamInfo fBeamInfo;
  sbn::selTriggerInfo fTriggerInfo;
  sbn::selLightInfo fFlashInfo;
  std::vector<sbn::selLightInfo> fFlashStore;
  sbn::selHitInfo fHitInfo;
  std::vector<sbn::selHitInfo> fHitStore;
  // --- END ---- Tree buffers -------------------------------------------------
  
  // --- BEGIN -- Cached information -------------------------------------------
  
  /// PMT geometry objects, grouped by wall (drift) coordinate.
  std::vector<std::pair<double, std::vector<raw::Channel_t>>> fPMTwalls;
  
  // --- BEGIN -- Cached information -------------------------------------------
  
  
  TTree *fStoreTree = nullptr;

  unsigned int fTotalProcessed = 0;
  
  
  // `convert()` needs to be a free function
  friend TriggerInputSpec_t convert(Config::TriggerSpecConfig const& config);
  
  
  // --- BEGIN -- Trigger response branches ------------------------------------
  
  /// Manages filling of trigger result branches.
  sbn::details::TriggerResponseManager fTriggerResponses;
  
  // --- END ---- Trigger response branches ------------------------------------
  
  /// Accesses detector geometry to return all PMT split by wall.
  std::vector<std::pair<double, std::vector<raw::Channel_t>>>
    computePMTwalls() const;

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
  , fTrackFitterProducer { p().TrackFitterProducer() }
  , fCaloProducer     { p().CaloProducer() }
  , fBeamGateProducer { p().BeamGateProducer() }
  , fTriggerProducer  { p().TriggerProducer() }
  , fFlashProducer    { p().FlashProducer() }
  , fLogCategory      { p().LogCategory() }
  , fMODA             { p().MODA() }
  , fMODB             { p().MODB() }
  , fWion             { p().Wion() }
  , fEfield           { p().Efield() }
  , fForceDowngoing    { p().ForceDowngoing() }
  // algorithms
  , fPMTwalls         { computePMTwalls() }
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
  consumes<std::vector<art::Ptr<recob::Track>>>(fT0selProducer);
  consumes<sbn::ExtraTriggerInfo>(fTriggerProducer);
  consumes<std::vector<sim::BeamGateInfo>>(fBeamGateProducer);
  consumes<art::Assns<recob::PFParticle, recob::Track>>(fTrackProducer);
  consumes<art::Assns<recob::Track, anab::T0>>(fT0Producer);
  consumes<art::Assns<recob::Hit, recob::Track, recob::TrackHitMeta>>(fTrackProducer);
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
  fStoreTree->Branch("selFlashes", &fFlashStore); //store all flashes in an event for all tracks
  fStoreTree->Branch("selHits", &fHitStore);
  //fStoreTree->Branch("selFlashes", &fFlashInfo);
  
} // sbn::TimeTrackTreeStorage::TimeTrackTreeStorage()


void sbn::TimeTrackTreeStorage::analyze(art::Event const& e)
{
  const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();
  // Implementation of required member function here.
  fEvent = e.event();
  fSubRun = e.subRun();
  fRun = e.run();
  fBeamInfo = {};
  
  std::vector<art::Ptr<recob::Track>> const& tracks = e.getProduct<std::vector<art::Ptr<recob::Track>>> (fT0selProducer);
  if(tracks.empty()) {
    mf::LogDebug(fLogCategory) << "No tracks in '" << fT0selProducer.encode() << "'.";
    return;
  }

  std::vector<sim::BeamGateInfo> const& beamgate = e.getProduct<std::vector<sim::BeamGateInfo>> (fBeamGateProducer);
  if(beamgate.empty())
    mf::LogWarning(fLogCategory) << "No Beam Gate Information!";
  else {
    if(beamgate.size() > 1)
      mf::LogWarning(fLogCategory) << "Event has multiple beam gate info labels! (maybe this changes later to be standard)";
    sim::BeamGateInfo const& bg = beamgate[0];
    fBeamInfo.beamGateSimStart = bg.Start();
    fBeamInfo.beamGateDuration = bg.Width();
    fBeamInfo.beamGateType = bg.BeamType();
  }

  sbn::ExtraTriggerInfo const &triggerinfo = e.getProduct<sbn::ExtraTriggerInfo> (fTriggerProducer);
  fTriggerInfo.beamType = value(triggerinfo.sourceType);
  fTriggerInfo.triggerTime = triggerinfo.triggerTimestamp;
  fTriggerInfo.beamGateTime = triggerinfo.beamGateTimestamp;
  fTriggerInfo.triggerID = triggerinfo.triggerID;
  fTriggerInfo.gateID = triggerinfo.gateID;
  
  //mf::LogTrace(fLogCategory) << "HERE!";
  //art::FindOneP<recob::Track> particleTracks(pfparticles,e,fTrackProducer);
  art::FindOneP<anab::T0> t0Tracks(tracks,e,fT0Producer);  
  std::vector<recob::OpFlash> const &particleFlashes = e.getProduct<std::vector<recob::OpFlash>>(fFlashProducer);
  //art::FindOneP<recob::SpacePoint> particleSPs(pfparticles, e, fT0selProducer);
  //mf::LogTrace(fLogCategory) << "PFParticles size: " << pfparticles.size() << " art::FindOneP Tracks Size: " << particleTracks.size();
  art::ValidHandle<std::vector<recob::Track>> allTracks = e.getValidHandle<std::vector<recob::Track>>(fTrackProducer);
  //auto particles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  //art::FindManyP<anab::T0> pandora_t0Tracks(particles, e, fPFPproducer);
  art::FindManyP<recob::Hit,recob::TrackHitMeta> trkht(allTracks,e,fTrackProducer); //for track hits
  art::FindManyP<anab::Calorimetry> calorim(allTracks, e, fCaloProducer);  
  
  // get an extractor bound to this event
  sbn::details::TriggerResponseManager::Extractors triggerResponseExtractors
    = fTriggerResponses.extractorsFor(e);
  unsigned int processed = 0;
  
  for(unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack)
  {
    art::Ptr<recob::Track> const& trackPtr = tracks.at(iTrack);
    if(trackPtr.isNull()) continue;
    
    fFlashInfo = {};
    fFlashStore.clear();
    fHitStore.clear();
    
    art::Ptr<anab::T0> const& t0Ptr = t0Tracks.at(iTrack);
    //auto pandora_t0Ptr = pandora_t0Tracks.at(trackPtr.key());
    float const track_t0 = t0Ptr->Time();
    //float pandora_t0 = -99999.0;
    //if(!pandora_t0Ptr.empty())
    //{
    //float const pandora_t0 = pandora_t0Ptr->Time();
    //}
    //std::cout << "CRT-TPC matching algo T0: " << track_t0/1e3 << "Pandora t0 value: " << pandora_t0/1e3 << std::endl; 
    if(!particleFlashes.empty())
    {
      //float min_flash_t0_diff = 999999.0; 
      for(unsigned int iFlash = 0; iFlash < particleFlashes.size(); ++iFlash)
      {
        fFlashInfo = {};
        recob::OpFlash const flashPtr = particleFlashes.at(iFlash);
        float const flash_pe = flashPtr.TotalPE();
        float const flash_time = flashPtr.Time();
        float const flash_x = flashPtr.XCenter();
        float const flash_y = flashPtr.YCenter();
        float const flash_z = flashPtr.ZCenter();
        float flash_t0_diff = flash_time - track_t0/1e3;
        //if(std::abs(flash_t0_diff) < min_flash_t0_diff)
          //{ 
        fFlashInfo.flash_id = iFlash;
        fFlashInfo.sum_pe = flash_pe;
        fFlashInfo.flash_time = flash_time;
        fFlashInfo.flash_x = flash_x;
        fFlashInfo.flash_y = flash_y;
        fFlashInfo.flash_z = flash_z;
        fFlashInfo.diff_flash_t0 = flash_t0_diff;
        //min_flash_t0_diff = std::abs(flash_t0_diff);
        fFlashStore.push_back(fFlashInfo);
        //}
      }
    }
    sbn::selTrackInfo trackInfo;
    trackInfo.trackID = trackPtr->ID();
    trackInfo.t0 = track_t0/1e3; //is this in nanoseconds? Will convert to seconds so I can understand better
    //if(track_t0/1e3 < 10 && track_t0/1e3 > -10)
    //mf::LogTrace(fLogCategory) << track_t0/1e3 << " Run is: " << fRun << " SubRun is: " << fSubRun << " Event is: " << fEvent << " Track ID is: " << trackPtr->ID();
    
    recob::tracking::Point_t startPoint = trackPtr->Start();
    recob::tracking::Point_t endPoint = trackPtr->End();
    recob::tracking::Vector_t startDir = trackPtr->StartDirection();
    bool const flipTrack = fForceDowngoing && (startDir.Y() > 0.0);
    if (flipTrack) {
      std::swap(startPoint, endPoint);
      startDir = -trackPtr->EndDirection();
    }
    
    trackInfo.start_x = startPoint.X();
    trackInfo.start_y = startPoint.Y();
    trackInfo.start_z = startPoint.Z();
    trackInfo.end_x = endPoint.X();
    trackInfo.end_y = endPoint.Y();
    trackInfo.end_z = endPoint.Z();
    trackInfo.dir_x = startDir.X();
    trackInfo.dir_y = startDir.Y();
    trackInfo.dir_z = startDir.Z();
    trackInfo.length = trackPtr->Length();
    
    std::vector<geo::Point_t> const trackPath
      = extractTrajectory(*trackPtr, flipTrack);
    geo::Point_t const middlePoint
      = util::pathMiddlePoint(trackPath.begin(), std::prev(trackPath.end()));
    trackInfo.middle_x = middlePoint.X();
    trackInfo.middle_y = middlePoint.Y();
    trackInfo.middle_z = middlePoint.Z();
    
    //
    // determination of cathode crossing
    //
    // this determination should be independent of track direction;
    icarus::CathodeDesc_t const cathode
      = icarus::findTPCcathode(middlePoint, *geom);
    
    icarus::CathodeCrossing_t crossInfo
      = icarus::detectCrossing(trackPath.begin(), trackPath.end(), cathode);
    
    if (crossInfo) {
      
      auto itBeforeCathode = trackPath.begin() + crossInfo.indexBefore;
      auto itAfterCathode = trackPath.begin() + crossInfo.indexAfter;
      
      geo::Point_t middleBeforeCathode
        = util::pathMiddlePoint(trackPath.begin(), itBeforeCathode);
      geo::Point_t middleAfterCathode
        = util::pathMiddlePoint(itAfterCathode, std::prev(trackPath.end()));
      
        // "before" is defined as "smaller x", so:
      if (distance(middleAfterCathode, cathode) < 0.0) {
        assert(distance(middleBeforeCathode, cathode) >= 0.0);
        std::swap(crossInfo.indexBefore, crossInfo.indexAfter);
        std::swap(itBeforeCathode, itAfterCathode);
        std::swap(crossInfo.before, crossInfo.after);
        std::swap(middleBeforeCathode, middleAfterCathode);
      }
      
      trackInfo.beforecathode = crossInfo.before;
      trackInfo.aftercathode = crossInfo.after;
      
      geo::Point_t const& atCathodePoint = crossInfo.crossingPoint;
      trackInfo.atcathode_x = atCathodePoint.X();
      trackInfo.atcathode_y = atCathodePoint.Y();
      trackInfo.atcathode_z = atCathodePoint.Z();
      
      trackInfo.midbeforecathode_x = middleBeforeCathode.X();
      trackInfo.midbeforecathode_y = middleBeforeCathode.Y();
      trackInfo.midbeforecathode_z = middleBeforeCathode.Z();
      
      trackInfo.midaftercathode_x = middleAfterCathode.X();
      trackInfo.midaftercathode_y = middleAfterCathode.Y();
      trackInfo.midaftercathode_z = middleAfterCathode.Z();
      
    } // if crossing
    
    //Animesh added hit information - 2/8/2022

    unsigned int plane = 0; //hit plane number

    std::vector<art::Ptr<recob::Hit>> const& allHits = trkht.at(trackPtr.key());
    std::vector<const recob::TrackHitMeta*> const& trkmetas = trkht.data(trackPtr.key());
    std::vector<art::Ptr<anab::Calorimetry>> const& calorimetrycol = calorim.at(trackPtr.key());
    std::vector<std::vector<unsigned int>> hits(plane);

    // art::FindManyP<recob::SpacePoint> fmspts(allHits, e, fT0selProducer);
    // or art::FindManyP<recob::SpacePoint> fmspts(allHits, e, fSpacePointModuleLabel); // Not sure how to define it
    for (size_t ih = 0; ih < allHits.size(); ++ih)
    {
      //hits[allHits[ih]->WireID().Plane].push_back(ih);
      sbn::selHitInfo hinfo = makeHit(*allHits[ih], allHits[ih].key(), *trackPtr, *trkmetas[ih], calorimetrycol, geom);
      if(hinfo.plane == 2)
        fHitStore.push_back(hinfo);

    }
    float totE = 0;
    float totq_int = 0;
    float totq_dqdx = 0;
    for (size_t i = 0; i < fHitStore.size(); ++i)
    {
      if(fHitStore[i].dEdx > -1)
      {
        float E_hit = fHitStore[i].dEdx*fHitStore[i].pitch; //energy of hit, in MeV?
        totE += E_hit;
      }

      if(fHitStore[i].dqdx > -1)
      {
        float q_hit = fHitStore[i].integral;
        totq_int += q_hit;
        float q_hit_dqdx = fHitStore[i].dqdx*fHitStore[i].pitch;
        totq_dqdx += q_hit_dqdx;
      }
      /*
      if(fHitStore[i].pitch > 1 && fHitStore[i].dqdx < 100)
      {
        std::cout << "In strange peak! Event: " << fEvent << " Track ID: " << trackInfo.trackID << " Hit ID: " << i << " Total Number of hits: " << fHitStore.size() << std::endl;
      }
      */
    }
    trackInfo.energy = totE;
    trackInfo.charge_int = totq_int;
    trackInfo.charge_dqdx = totq_dqdx;
    fTrackInfo = std::move(trackInfo);
    /*
    for(size_t trajp = 0; trajp < trackPtr->NumberTrajectoryPoints()-1; ++trajp)
    {
      TVector3 cur_point(trackPtr->TrajectoryPoint(traj_p).position.X(), trackPtr->TrajectoryPoint(traj_p).position.Y(), trackPtr->TrajectoryPoint(traj_p).position.Z());
      TVector3 next_point(trackPtr->TrajectoryPoint(traj_p+1).position.X(), trackPtr->TrajectoryPoint(traj_p+1).position.Y(), trackPtr->TrajectoryPoint(traj_p+1).position.Z());
      if(abs(cur_point.X()) < 170 && abs(next_point.X()) > 170)
        //interpolate to get cathode crossing point
        
    }
    */
    
    triggerResponseExtractors.fetch(iTrack); // TODO no check performed; need to find a way
  
    ++processed;
    ++fTotalProcessed;
    fStoreTree->Fill();
  } // for particle
  
  mf::LogInfo(fLogCategory) << "Particles Processed: " << processed;

} // sbn::TimeTrackTreeStorage::analyze()


sbn::selHitInfo sbn::TimeTrackTreeStorage::makeHit(const recob::Hit &hit,
                                                   unsigned hkey,
                                                   const recob::Track &trk,
                                                   const recob::TrackHitMeta &thm,
                                                   const std::vector<art::Ptr<anab::Calorimetry>> &calo,
                                                   const geo::GeometryCore *geo)
{

  // TrackHitInfo to save
  sbn::selHitInfo hinfo;

  // information from the hit object
  hinfo.integral = hit.Integral();
  hinfo.sumadc = hit.SummedADC();
  hinfo.width = hit.RMS();
  hinfo.pk_time = hit.PeakTime();
  hinfo.mult = hit.Multiplicity();
  hinfo.wire = hit.WireID().Wire;
  hinfo.plane = hit.WireID().Plane;
  hinfo.channel = geo->PlaneWireToChannel(hit.WireID());
  hinfo.tpc = hit.WireID().TPC;
  hinfo.end = hit.EndTick();
  hinfo.start = hit.StartTick();
  hinfo.id = (int)hkey;

  bool const badhit = (thm.Index() == std::numeric_limits<unsigned int>::max()) ||
    (!trk.HasValidPoint(thm.Index()));

  //hinfo.ontraj = !badhit;
  // Save trajectory information if we can
  if(!badhit)
  {
    geo::Point_t const& loc = trk.LocationAtPoint(thm.Index());
    hinfo.px = loc.X();
    hinfo.py = loc.Y();
    hinfo.pz = loc.Z();
  
    geo::Vector_t const& dir = trk.DirectionAtPoint(thm.Index());
    hinfo.dirx = dir.X();
    hinfo.diry = dir.Y();
    hinfo.dirz = dir.Z();
    
    // And determine if the Hit is on a Calorimetry object
    for (const art::Ptr<anab::Calorimetry> &c: calo) {
      if (c->PlaneID().Plane != hinfo.plane) continue;
      
      // Found the plane! Now find the hit:
      for (unsigned i_calo = 0; i_calo < c->dQdx().size(); i_calo++) {
        // "TpIndices" match to the hit key
        if (c->TpIndices()[i_calo] != hkey) continue;
        // Fill the calo information associated with the hit 
        hinfo.oncalo = true;
        hinfo.pitch = c->TrkPitchVec()[i_calo];
        hinfo.dqdx = c->dQdx()[i_calo];
        hinfo.dEdx = dEdx_calc(hinfo.dqdx, fMODA, fMODB, fWion, fEfield);
        hinfo.rr = c->ResidualRange()[i_calo];
        break;
      } // for i_calo
      break;
    } // for c
  }
  
  return hinfo;
}

float sbn::TimeTrackTreeStorage::dEdx_calc(float dQdx,
                                           float A,
                                           float B,
                                           float Wion,
                                           float E) 
{
  float LAr_density_gmL = 1.389875; //LAr density in g/cm^3
  float alpha = A;
  float beta = B/(LAr_density_gmL*E);
  float dEdx = ((std::exp(dQdx*Wion*beta) - alpha)/beta)*3.278;
 
  return dEdx;
  
}


// -----------------------------------------------------------------------------
void sbn::TimeTrackTreeStorage::endJob() {
  
  mf::LogInfo(fLogCategory) << "Processed " << fTotalProcessed << " tracks.";
  
} // sbn::TimeTrackTreeStorage::endJob()


// -----------------------------------------------------------------------------
std::vector<std::pair<double, std::vector<raw::Channel_t>>>
sbn::TimeTrackTreeStorage::computePMTwalls() const {
  
  geo::GeometryCore const& geom { *lar::providerFrom<geo::Geometry>() };
  
  // run the algorithm to identify the PMT walls (as groups of geo::OpDetGeo)
  std::vector<std::pair<double, std::vector<geo::OpDetGeo const*>>> opDetWalls
    = icarus::trigger::PMTverticalSlicingAlg{}.PMTwalls(geom);
  
  // and weirdly, the only portable way to go from a OpDetGeo to its channel
  // is to build a map (maybe because it's not guaranteed to be 1-to-1?)
  std::map<geo::OpDetGeo const*, raw::Channel_t> opDetToChannel;
  for (auto const channel: util::counter<raw::Channel_t>(geom.MaxOpChannel()))
    opDetToChannel[&geom.OpDetGeoFromOpChannel(channel)] = channel;
  
  // rewrite the data structure replacing each detector with its readout channel
  std::vector<std::pair<double, std::vector<raw::Channel_t>>> channelWalls;
  for (auto const& [ coord, PMTs ]: opDetWalls) {
    std::vector<raw::Channel_t> channels;
    channels.reserve(PMTs.size());
    for (geo::OpDetGeo const* PMT: PMTs)
      channels.push_back(opDetToChannel.at(PMT));
    
    channelWalls.emplace_back(coord, std::move(channels));
  } // for walls
  
  return channelWalls;
} // sbn::TimeTrackTreeStorage::computePMTwalls()


// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(sbn::TimeTrackTreeStorage)
