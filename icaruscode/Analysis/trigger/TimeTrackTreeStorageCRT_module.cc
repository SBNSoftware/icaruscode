/**
 * @file    TimeTrackTreeStorageCRT_module.cc
 * @authors Jacob Zettlemoyer (FNAL, jzettle@fnal.gov),
 *          Animesh Chatterjee (U. Pittsburgh, anc238@pitt.edu),
 *          Gianluca Petrillo (SLAC, petrillo@slac.stanford.edu)
 * @date    February 2023
 * 
 * Originally borrowed heavily from Gray Putnam's existing `TrackCaloSkimmer`.
 */

// ICARUS libraries
#include "icaruscode/Analysis/trigger/details/TriggerResponseManager.h"
#include "icaruscode/Analysis/trigger/details/CathodeCrossingUtils.h"
#include "icaruscode/Analysis/trigger/Objects/TrackTreeStoreObj.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h"
#include "icaruscode/Utilities/TrajectoryUtils.h"
#include "icaruscode/Utilities/ArtAssociationCaches.h" // OneToOneAssociationCache
#include "icarusalg/Utilities/TrackTimeInterval.h"
#include "icarusalg/Utilities/AssnsCrosser.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "lardataalg/Utilities/StatCollector.h" // lar::util::MinMaxCollector
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
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
  std::pair<std::vector<geo::Point_t>, std::vector<geo::Vector_t>>
  extractTrajectory
    (recob::Track const& track, bool reverse = false, geo::Vector_t shift = {})
  {
    std::vector<geo::Point_t> trackPath;
    std::vector<geo::Vector_t> trackMom;
    std::size_t index = track.FirstValidPoint();
    while (index != recob::TrackTrajectory::InvalidIndex) {
      trackPath.push_back(track.LocationAtPoint(index) + shift);
      trackMom.push_back(track.MomentumVectorAtPoint(index));
      if (++index >= track.NPoints()) break;
      index = track.NextValidPoint(index);
    }
    if (reverse) {
      std::reverse(trackPath.begin(), trackPath.end());
      std::reverse(trackMom.begin(), trackMom.end());
      for (geo::Vector_t& mom: trackMom) mom *= -1.0;
    }
    return { std::move(trackPath), std::move(trackMom) };
  } // extractTrajectory()
  
  
  /// Returns the input tag of data product the `ptr` points to.
  template <typename T>
  art::InputTag inputTagOf(art::Ptr<T> ptr, art::Event const& event)
    { return event.getProductDescription(ptr.id())->inputTag(); }
  
} // local namespace



// -----------------------------------------------------------------------------
namespace sbn {
  class TimeTrackTreeStorage; 
  
}
/**
 * @brief Fills a ROOT tree with track-based triggering information.
 * 
 * 
 * Track selection criteria
 * =========================
 * 
 * This module does not perform almost any selection: it processes all the
 * selected track, output of the module in `T0selProducer`, which include a
 * track (`recob::Track`) and its time (`anab::T0`).
 * That module may perform the desired selection; this one will use all tracks
 * that have a time, even if that time is apparently invalid.
 * 
 * For every track which is associated to a PFO (`recob::PFParticle`),
 * the time (`anab::T0`) associated to that PFO is saved in a tree entry.
 * 
 * For every track which is matched to a CRT hit (from `CRTMatchingProducer`
 * module), the time of that match is saved, together with some information
 * about the hit itself.
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
 * Data stored in the tree
 * ========================
 * 
 * Tracks
 * -------
 * 
 * ### Position in the drift direction
 * 
 * The input tracks are expected to be processed in a Pandora-like way.
 * That is, if a TPC time is associated with them, then they have been already
 * translated along the drift direction (_x_) to match the assigned time.
 * 
 * In case the tracks have no TPC time associated with them, instead, the time
 * associated from the CRT (not the one from the `T0selProducer` source),
 * is considered. If that time is valid, then it is used to translate the track
 * in the position that pertains a track of that time. The track in input is
 * assumed to have been placed on the drift direction in the position that would
 * pertain a track arrived at trigger time. No correction for the time of flight
 * between the CRT hit and the track is performed (this correction, of the order
 * of nanoseconds, is irrelevant to the purpose of trigger efficiency
 * evaluation).
 * 
 * 
 * 
 * Calorimetry
 * ------------
 * 
 * Track energy and energy profile are directly taken from the results of
 * "standard" ICARUS calorimetry application.
 * Calibration hit by hit is currently still custom, using recombination
 * correction.
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
 * Service dependencies
 * =====================
 * 
 * * `Geometry` for determinations related to the position of the cathode and
 *   anodes, and for association of track and hits.
 * * `LArProperties` as dependency for `DetectorPropertiesService`
 * * `DetectorClocksService` as dependency for `DetectorPropertiesService`
 * * `DetectorPropertiesService` for the drift velocity to convert a time shift
 *   into a drift direction shift.
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
    
    
    fhicl::Atom<art::InputTag> T0selProducer {
      Name("T0selProducer"),
      Comment("tag of the tracks to process (as a collection of art::Ptr)")
      };
    
    fhicl::OptionalAtom<art::InputTag> PFPproducer {
      Name("PFPproducer"),
      Comment("tag of associations of particle flow objects and input tracks")
      };
    
    fhicl::OptionalAtom<art::InputTag> T0Producer {
      Name("T0Producer"),
      Comment("tag of the time (t0) of particle flow objects [as PFPproducer]")
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

    fhicl::Atom<art::InputTag> CRTMatchingProducer {
      Name("CRTMatchingProducer"),
      Comment("tag of module used to associate TPC tracks and CRT hits"),
      ""
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
    
    fhicl::Atom<int> CalorimetryPlane {
      Name("CalorimetryPlane"),
      Comment(
        "number of the plane to pick calorimetry from"
        " (0 is closest to cathode; -1: the highest available)"
        ),
      -1
      };
    
  }; // Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  explicit TimeTrackTreeStorage(Parameters const& p);

  sbn::selHitInfo makeHit(const recob::Hit &hit,
                          unsigned hkey,
                          const recob::Track &trk,
                          const recob::TrackHitMeta &thm,
                          bool flippedTrack,
                          const std::vector<anab::Calorimetry const*> &calo,
                          const geo::GeometryCore *geo);

  float dEdx_calc(float dQdx, 
                  float A,
                  float B,
                  float Wion,
                  float E);

  void analyze(art::Event const& e) override;
  
  void endJob() override;

private:
  
  using electronics_time = detinfo::timescales::electronics_time;
  using microseconds = util::quantities::intervals::microseconds;
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  art::InputTag const fT0selProducer;
  art::InputTag const fTrackProducer;
  art::InputTag const fPFPproducer; // producer of PFP/track associations
  art::InputTag const fT0Producer;
  art::InputTag const fTrackFitterProducer;
  art::InputTag const fCaloProducer;
  art::InputTag const fBeamGateProducer;
  art::InputTag const fTriggerProducer;
  art::InputTag const fFlashProducer;
  art::InputTag const fCRTMatchProducer;
  std::string const fLogCategory;
  float const fMODA;
  float const fMODB;
  float const fWion;
  float const fEfield;
  bool const fForceDowngoing; ///< Whether to force all tracks to be downgoing.
  int const fCalorimetryPlaneNumber; ///< Which plane to use for calorimetry.
  
  // --- END ---- Configuration parameters -------------------------------------

  // --- BEGIN -- Tree buffers -------------------------------------------------
  
  unsigned int fEvent;
  unsigned int fRun;
  unsigned int fSubRun;
  
  sbn::selTrackInfo fTrackInfo; //change to one entry per track instead of per event 
  sbn::selBeamInfo fBeamInfo;
  sbn::selTriggerInfo fTriggerInfo;
  sbn::selLightInfo fClosestFlash;
  sbn::selLightInfo fNearestFlash;
  std::vector<sbn::selLightInfo> fFlashStore;
  sbn::selHitInfo fHitInfo;
  std::vector<sbn::selHitInfo> fHitStore;
  sbn::selCRTInfo fCRTInfo;
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
  
  /// Extracts CRT information pertaining the specified `track`.
  sbn::selCRTInfo extractCRTinfoFor(
    art::Ptr<recob::Track> const& track, art::Event const& event,
    util::OneToOneAssociationCache<recob::Track, anab::T0> const& trackToCRTt0,
    util::OneToOneAssociationCache<anab::T0, sbn::crt::CRTHit> const& t0ToCRThits
    ) const;
  
  /// Accesses detector geometry to return all PMT split by wall.
  std::vector<std::pair<double, std::vector<raw::Channel_t>>>
    computePMTwalls() const;
  
  /// Returns the position correction of a `track` based on a time shift.
  /// ("actual" track time minus the time assumed during reconstruction).
  /// A non-empty list of the TPCs crossed by the track is required.
  geo::Vector_t posShiftFromCRTtime(
    recob::Track const& track,
    double timeShift,
    std::vector<geo::TPCID> const& TPCs,
    geo::GeometryCore const& geom,
    detinfo::DetectorPropertiesData const& detProp
    ) const;
  
  /**
   * @brief Returns the most significative distance of a `time` from a `range`.
   * @param time time to be checked
   * @param range time range `time` is going to be compared to
   * @returns the most significative distance of a `time` from a `range`
   * 
   * If `time` is contained in the range, the return value is negative and its
   * modulus is the distance from the closest `range` boundary; no information
   * is encoded about which of the boundaries that is.
   * 
   * If `time` is not contained in the range, the return value is positive and
   * its modulus is again the distance from the closest `range` boundary, which
   * is also the one crossed by `time`.
   * 
   * The rule also holds for ranges where the start time is larger than the stop
   * time. If `range` is not valid, it is interpreted as a range containing all
   * time, and the return value is `0`.
   */
  static microseconds distanceFromTimeRange(
    electronics_time time, lar::util::TrackTimeInterval::TimeRange const& range
    );
    
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
namespace {
  
  /// A vector-like collection which inserts elements only if not there yet.
  /// It's designed for small collections.
  template <typename... Args>
  class UniqueVector: private std::vector<Args...> {
    using Vector_t = std::vector<Args...>;
      public:
    using value_type = typename Vector_t::value_type;
    using const_iterator = typename Vector_t::const_iterator;
    
    /// Adds a copy of `value` if it's not present yet.
    void push_back(value_type value)
      {
        auto const end = Vector_t::crend();
        if (std::find(Vector_t::crbegin(), end, value) != end) return;
        Vector_t::push_back(std::move(value));
      }
    
    using Vector_t::empty;
    using Vector_t::size;
    using Vector_t::back;
    using Vector_t::front;
    using Vector_t::cbegin;
    using Vector_t::cend;
    auto begin() const { return cbegin(); }
    auto end() const { return cend(); }
    
    Vector_t const& vector() const { return *this; }
    
  }; // class UniqueVector
  
} // local namespace


// -----------------------------------------------------------------------------
sbn::TimeTrackTreeStorage::TimeTrackTreeStorage(Parameters const& p)
  : EDAnalyzer{p}
  // configuration
  , fT0selProducer    { p().T0selProducer() }
  , fTrackProducer    { p().TrackProducer() }
  , fPFPproducer      { p().PFPproducer().value_or(fTrackProducer) }
  , fT0Producer       { p().T0Producer().value_or(fPFPproducer) }
  , fTrackFitterProducer { p().TrackFitterProducer() }
  , fCaloProducer     { p().CaloProducer() }
  , fBeamGateProducer { p().BeamGateProducer() }
  , fTriggerProducer  { p().TriggerProducer() }
  , fFlashProducer    { p().FlashProducer() }
  , fCRTMatchProducer { p().CRTMatchingProducer() }
  , fLogCategory      { p().LogCategory() }
  , fMODA             { p().MODA() }
  , fMODB             { p().MODB() }
  , fWion             { p().Wion() }
  , fEfield           { p().Efield() }
  , fForceDowngoing    { p().ForceDowngoing() }
  , fCalorimetryPlaneNumber { p().CalorimetryPlane() }
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
  
  consumes<std::vector<art::Ptr<recob::Track>>>(fT0selProducer);
  consumes<art::Assns<recob::Track, anab::T0>>(fT0selProducer);
  consumes<std::vector<sim::BeamGateInfo>>(fBeamGateProducer);
  consumes<sbn::ExtraTriggerInfo>(fTriggerProducer);
  consumes<art::Assns<recob::Track, recob::PFParticle>>(fTrackProducer);
  consumes<art::Assns<recob::Track, anab::Calorimetry>>(fCaloProducer);
  consumes<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>(fTrackProducer);
  consumes<recob::OpFlash>(fFlashProducer);
  if (!fCRTMatchProducer.empty()) {
    consumes<art::Assns<recob::Track, anab::T0>>(fCRTMatchProducer);
    consumes<art::Assns<anab::T0, sbn::crt::CRTHit>>(fCRTMatchProducer);
  }

  //
  // tree population
  //
  fStoreTree->Branch("run", &fRun);
  fStoreTree->Branch("subrun", &fSubRun);
  fStoreTree->Branch("event", &fEvent);
  fStoreTree->Branch("beamInfo", &fBeamInfo);
  fStoreTree->Branch("triggerInfo", &fTriggerInfo);
  fStoreTree->Branch("selTracks", &fTrackInfo);
  fStoreTree->Branch("closestFlash", &fClosestFlash);
  fStoreTree->Branch("nearestFlash", &fNearestFlash);
  fStoreTree->Branch("selFlashes", &fFlashStore); //store all flashes in an event for all tracks
  fStoreTree->Branch("selHits", &fHitStore);
  if (!fCRTMatchProducer.empty())
    fStoreTree->Branch("selCRTHits", &fCRTInfo); // one entry per track
  
} // sbn::TimeTrackTreeStorage::TimeTrackTreeStorage()


void sbn::TimeTrackTreeStorage::analyze(art::Event const& e)
{
  /*
   * The flow is the following:
   *  * driving the process it is the T0 selector output:
   *      * its list of tracks is processed
   *      * if there is an association to a T0, that becomes `t0` of the track
   *  * time from TPC reconstruction (`t0_TPC`) is taken from PFParticle reco
   *  * time from CRT (`t0_CRT`) is taken from the CRT/TPC matching
   * Usually `t0` is equal to either one of the other two times, and one of them
   * can be missing (or better, with a special invalid value).
   * 
   */
  
  using namespace util::quantities::time_literals; // ""_us
  using trigger_time = detinfo::timescales::trigger_time;
  
  //
  // reach for all the needed services (geometry, timing and detector state)
  //
  const geo::GeometryCore *geom = lar::providerFrom<geo::Geometry>();
  
  detinfo::DetectorClocksData const detClocks
    = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  
  detinfo::DetectorPropertiesData const detProp
    = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor
      (e, detClocks)
    ;
  
  detinfo::DetectorTimings const detTimings{ detClocks };
  
  lar::util::TrackTimeInterval const allowedTrackTime
    { *geom, detProp, detTimings };
  
  //
  // event identification
  //
  fEvent = e.event();
  fSubRun = e.subRun();
  fRun = e.run();
  fBeamInfo = {};
  
  //
  // hardware beam gate information
  //
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

  //
  // hardware trigger information
  //
  sbn::ExtraTriggerInfo const &triggerinfo = e.getProduct<sbn::ExtraTriggerInfo> (fTriggerProducer);
  fTriggerInfo.beamType = value(triggerinfo.sourceType);
  fTriggerInfo.triggerTime = triggerinfo.triggerTimestamp;
  fTriggerInfo.beamGateTime = triggerinfo.beamGateTimestamp;
  fTriggerInfo.triggerID = triggerinfo.triggerID;
  fTriggerInfo.gateID = triggerinfo.gateID;
  
  //
  // trigger emulation support
  //
  // get an extractor bound to this event
  sbn::details::TriggerResponseManager::Extractors triggerResponseExtractors
    = fTriggerResponses.extractorsFor(e);
  
  //
  // tracks
  //
  
  std::vector<art::Ptr<recob::Track>> const& tracks = e.getProduct<std::vector<art::Ptr<recob::Track>>> (fT0selProducer);
  if(tracks.empty()) {
    mf::LogDebug(fLogCategory) << "No tracks in '" << fT0selProducer.encode() << "'.";
    return;
  }

  
  //
  // additional associations to navigate from tracks to other information
  //
  
  // t0 from selection
  art::FindOneP<anab::T0> TrackT0s(tracks,e,fT0selProducer);
  
  // t0 from TPC (cathode crossing tracks)
  icarus::ns::util::AssnsCrosser<recob::Track, recob::PFParticle, anab::T0> const
  TrackTPCt0s{ e, fPFPproducer, fT0Producer };
  
  // optical flashes
  std::vector<recob::OpFlash> const &particleFlashes = e.getProduct<std::vector<recob::OpFlash>>(fFlashProducer);
  
  // track calorimetry
  art::FindMany<anab::Calorimetry> trackCalorimetry(tracks, e, fCaloProducer);
  
  // track TPC hits
  art::FindManyP<recob::Hit,recob::TrackHitMeta> trkht(tracks,e,fTrackProducer);
  
  // track t0 from CRT
  std::optional<util::OneToOneAssociationCache<recob::Track, anab::T0> const>
  TrackCRTt0s;
  if (!fCRTMatchProducer.empty()) {
    TrackCRTt0s.emplace
      (e.getProduct<art::Assns<recob::Track, anab::T0>>(fCRTMatchProducer));
  }
  // track t0 from CRT
  std::optional<util::OneToOneAssociationCache<anab::T0, sbn::crt::CRTHit> const>
  T0CRThits;
  if (!fCRTMatchProducer.empty()) {
    T0CRThits.emplace
      (e.getProduct<art::Assns<anab::T0, sbn::crt::CRTHit>>(fCRTMatchProducer));
  }
  
  unsigned int processed = 0;
  for(unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack)
  {
    art::Ptr<recob::Track> const& trackPtr = tracks.at(iTrack);
    if(trackPtr.isNull()) continue;
    
    // decide immediately if the track needs to be flipped
    bool const flipTrack
      = fForceDowngoing && (trackPtr->StartDirection().Y() > 0.0);
    
    //
    // matched CRT information
    //
    if (!fCRTMatchProducer.empty())
      fCRTInfo = extractCRTinfoFor(trackPtr, e, *TrackCRTt0s, *T0CRThits);
    bool const hasCRTT0 = (fCRTInfo.time != sbn::selCRTInfo::NoTime);
    
    //
    // hit information
    //
    
    //Animesh added hit information - 2/8/2022
    std::vector<anab::Calorimetry const*> const& calorimetry
      = trackCalorimetry.at(iTrack);
    std::vector<art::Ptr<recob::Hit>> const& allHits = trkht.at(iTrack);
    std::vector<const recob::TrackHitMeta*> const& trkmetas = trkht.data(iTrack);
    
    lar::util::MinMaxCollector<float> hitTimeRange;

    fHitStore.clear();
    UniqueVector<geo::TPCID> TPCs; // collect the list of TPC for later
    for (auto const& [ hitPtr, pTrkMeta ]: util::zip(allHits, trkmetas))
    {
      recob::Hit const& hit = *hitPtr;
      geo::WireID const wire = hit.WireID();
      if (!wire.isValid) continue;
      TPCs.push_back(wire);
      hitTimeRange.add(hit.PeakTime());
      if (wire.Plane != 2) continue;
      fHitStore.push_back(makeHit
        (hit, hitPtr.key(), *trackPtr, *pTrkMeta, flipTrack, calorimetry, geom)
        );
    }
    if (TPCs.empty()) {
      // this is needed for corrections; if needed, this check may be made optional,
      // but then a fix for the position correction needs to be devised.
      throw art::Exception{ art::errors::NotFound }
        << "Track ID=" << trackPtr->ID() << " is not associated with any hit.\n";
    }
    
    //
    // track information (at last!!)
    //
    sbn::selTrackInfo trackInfo;
    trackInfo.trackID = trackPtr->ID();
    trackInfo.flipped = flipTrack;
    
    art::Ptr<anab::T0> const& t0Ptr = TrackT0s.at(iTrack);
    if (!t0Ptr) {
      // this is the t0 that must not be missing,
      // since it's the one the trigger emulation is based on
      throw art::Exception{ art::errors::NotFound }
        << "Track #" << trackPtr.key() << " from '" << inputTagOf(trackPtr, e).encode()
        << "' (ID=" << trackPtr->ID() << ") is not associated to any anab::T0 from '"
        << fT0selProducer.encode() << "'!\n";
    }
    // it may happen that a T0 object is present, but with an invalid time;
    // in that case, the ID is expected not to be valid either;
    // here we assume that since the tracks were selected they are desired;
    // on the other end, all the trigger efficiency figures will be bogus
    bool const hasT0 = (t0Ptr->ID() >= 0);
    if (hasT0) trackInfo.t0 = t0Ptr->Time() / 1000.0; // otherwise stays NoTime
    
    art::Ptr<anab::T0> const& t0TPCPtr = TrackTPCt0s.assPtr(trackPtr);
    // differently from CRT, we assume that all available T0 from TPC are valid
    bool const hasTPCT0 = t0TPCPtr.isNonnull();
    if (hasTPCT0)
      trackInfo.t0_TPC = t0TPCPtr->Time() / 1000.0; // nanoseconds -> microseconds
    
    trackInfo.t0_CRT = fCRTInfo.time; // replica
    
    lar::util::TrackTimeInterval::TimeRange const trackTimeRange
      = allowedTrackTime.timeRangeOfHits(allHits);
    
    trackInfo.t0_TPC_min
      = detTimings.toTriggerTime(trackTimeRange.start).value();
    trackInfo.t0_TPC_max
      = detTimings.toTriggerTime(trackTimeRange.stop).value();
    
    electronics_time const t0_elec
      = detTimings.toElectronicsTime(trigger_time{ trackInfo.t0 });
    trackInfo.t0_diff = distanceFromTimeRange(t0_elec, trackTimeRange).value();
    
    electronics_time const t0_CRT_elec
      = detTimings.toElectronicsTime(trigger_time{ trackInfo.t0_CRT });
    trackInfo.t0_CRT_diff
      = distanceFromTimeRange(t0_CRT_elec, trackTimeRange).value();
    
    trackInfo.hitTick_min = hitTimeRange.min();
    trackInfo.hitTick_max = hitTimeRange.max();
    
    //
    // correction on position from time:
    //   if we have a t0 from TPC (cathode-crossing track) we assume it to be
    //   good enough and the track to be already stitched and correct.
    //   If instead we have at least the CRT time, correction is upon us;
    //   reference time for t0_CRT is trigger time,
    //   which is also the one assumed by TPC track reconstruction.
    //
    geo::Vector_t const positionShift = (!hasTPCT0 && hasCRTT0)
      ? posShiftFromCRTtime
        (*trackPtr, trackInfo.t0_CRT, TPCs.vector(), *geom, detProp)
      : geo::Vector_t{}
      ;
    trackInfo.driftCorrX = positionShift.X();

    // the trajectory in space is thoroughly flipped if required;
    // directions/momenta are not, and explicit treatment is needed.
    auto const& [ trackPath, trackMom ]
      = extractTrajectory(*trackPtr, flipTrack, positionShift);
    
    recob::tracking::Point_t const& startPoint = trackPath.front();
    recob::tracking::Point_t const& endPoint = trackPath.back();
    recob::tracking::Vector_t const& startDir = trackMom.front();
    
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
      if (positionShift != geo::Vector_t{}) {
        // in general, we correct only tracks that don't cross the cathode,
        // and for the others the correction should be 0;
        // but a contained track may have been pushed to cross cathode by the
        // correction (or maybe because of reconstruction imprecisions) by a bit
        constexpr double aBit = 5.0; // say, no more than this mucj
        if (std::min(
          std::abs(distance(trackPath.front(), cathode)),
          std::abs(distance(trackPath.back(), cathode))
          ) > aBit
        ) {
          mf::LogPrint(fLogCategory) 
            << "Track ID=" << trackPtr->ID()
            << " (" << startPoint << " to " << endPoint
            << ", crossing at " << crossInfo.crossingPoint
            << ") is scheduled for position correction by " << positionShift
            << " cm! (times: TPC=" << trackInfo.t0_TPC << ", CRT="
            << trackInfo.t0_CRT << " us).";
        }
      } // position shift check
      
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
    
    //
    // calorimetry information
    //
    
    // pick the plane of choice:
    anab::Calorimetry const* trackCalo = nullptr;
    if (fCalorimetryPlaneNumber < 0) {
      // pick the plane with the highest number (i.e. farther from cathode):
      for (anab::Calorimetry const* candCalo: calorimetry) {
        if (trackCalo && (candCalo->PlaneID().Plane < trackCalo->PlaneID().Plane))
          continue;
        trackCalo = candCalo;
      } // for
    }
    else {
      // find the calorimetry object for the requested plane:
      for (anab::Calorimetry const* candCalo: calorimetry) {
        if ((int) candCalo->PlaneID().Plane != fCalorimetryPlaneNumber) continue;
        trackCalo = candCalo;
        break;
      } // for
    }
    
    bool const hasCaloInfo = trackCalo && !trackCalo->TrkPitchVec().empty();
    if (hasCaloInfo) {
      trackInfo.energy_range = trackCalo->Range();
      trackInfo.energy = trackCalo->KineticEnergy();
      trackInfo.energy_int = 0.0;
      trackInfo.charge_int = 0.0;
      for (auto const [ dx, dEdx, dQdx ]
        : util::zip(trackCalo->TrkPitchVec(), trackCalo->dEdx(), trackCalo->dQdx())
      ) {
        if (dx <= 0.0) continue;
        trackInfo.energy_int += dEdx * dx;
        trackInfo.charge_int += dQdx * dx;
      } // for
    }
    else {
      mf::LogTrace{ "TimeTrackTreeStorage" }
        << "Track '" << inputTagOf(trackPtr, e).encode() << "' ID="
        << trackInfo.trackID << " has no suitable calorimetry information.";
    }
    
    fTrackInfo = std::move(trackInfo);
    
    //
    // optical flash information
    //
    fFlashStore.clear();
    // to use pointers below, need to avoid reallocation of vector data
    fFlashStore.reserve(particleFlashes.size());
    sbn::selLightInfo* closestFlash = nullptr;
    sbn::selLightInfo* nearestFlash = nullptr;
    // we load all of them in the tree so far
    for (auto const& [ iFlash, flash ]: util::enumerate(particleFlashes))
    {
      geo::Point_t const center{ flash.XCenter(), flash.YCenter(), flash.ZCenter() };
      float const flash_pe = flash.TotalPE();
      float const flash_time = flash.Time();
      float const flash_x = center.X();
      float const flash_y = center.Y();
      float const flash_z = center.Z();
      
      geo::Vector_t const diff_flash_pos = center - middlePoint;
      
      sbn::selLightInfo flashInfo;
      flashInfo.flash_id = iFlash;
      flashInfo.sum_pe = flash_pe;
      flashInfo.flash_time = flash_time;
      flashInfo.flash_x = flash_x;
      flashInfo.flash_y = flash_y;
      flashInfo.flash_z = flash_z;
      if (hasT0)    flashInfo.diff_flash_t0    = flash_time - fTrackInfo.t0;
      if (hasTPCT0) flashInfo.diff_flash_TPCt0 = flash_time - fTrackInfo.t0_TPC;
      if (hasCRTT0) flashInfo.diff_flash_CRTt0 = flash_time - fTrackInfo.t0_CRT;
      flashInfo.diff_flash_pos = std::hypot(diff_flash_pos.Y(), diff_flash_pos.Z());
      
      fFlashStore.push_back(flashInfo);
      
      if (hasT0 && (!closestFlash
        || std::abs(flashInfo.diff_flash_t0) < std::abs(closestFlash->diff_flash_t0))
      ) {
        closestFlash = &fFlashStore.back();
      }
      if (!nearestFlash
        || std::abs(flashInfo.diff_flash_pos) < std::abs(nearestFlash->diff_flash_pos)
      ) {
        nearestFlash = &fFlashStore.back();
      }
      
    } // for flashes
    
    if (closestFlash) closestFlash->flash_closest_to_track = true;
    if (nearestFlash) nearestFlash->flash_nearest_to_track = true;
    
    fClosestFlash = closestFlash? *closestFlash: sbn::selLightInfo{};
    fNearestFlash = nearestFlash? *nearestFlash: sbn::selLightInfo{};
    
    
    //
    // trigger emulation information
    //
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
                                                   bool flippedTrack,
                                                   const std::vector<anab::Calorimetry const*> &calo,
                                                   const geo::GeometryCore *geom)
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
  hinfo.channel = geom->PlaneWireToChannel(hit.WireID());
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
    // note that the track `trk` never comes flipped:
    // indices are still consistent between trajectory point, hit and metadata;
    // but all flipping needs to be done by hand here
    geo::Point_t const& loc = trk.LocationAtPoint(thm.Index());
    hinfo.px = loc.X();
    hinfo.py = loc.Y();
    hinfo.pz = loc.Z();
  
    geo::Vector_t const dir
      = (flippedTrack? -1: +1) * trk.DirectionAtPoint(thm.Index());
    hinfo.dirx = dir.X();
    hinfo.diry = dir.Y();
    hinfo.dirz = dir.Z();
    
    // And determine if the Hit is on a Calorimetry object
    for (anab::Calorimetry const* c: calo) {
      if (c->PlaneID().Plane != hinfo.plane) continue;
      
      auto const sortRange = [fullRange=c->Range(), flippedTrack]
        (float rr){ return flippedTrack? (fullRange - rr): rr; };
      
      // Found the plane! Now find the hit:
      for (unsigned i_calo = 0; i_calo < c->dQdx().size(); i_calo++) {
        // "TpIndices" match to the hit key
        if (c->TpIndices()[i_calo] != hkey) continue;
        // Fill the calo information associated with the hit 
        hinfo.oncalo = true;
        hinfo.pitch = c->TrkPitchVec()[i_calo];
        hinfo.dqdx = c->dQdx()[i_calo];
        hinfo.dEdx = c->dEdx()[i_calo];
        hinfo.rr = sortRange(c->ResidualRange()[i_calo]);
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


sbn::selCRTInfo sbn::TimeTrackTreeStorage::extractCRTinfoFor(
  art::Ptr<recob::Track> const& trackPtr, art::Event const& event,
  util::OneToOneAssociationCache<recob::Track, anab::T0> const& trackToCRTt0,
  util::OneToOneAssociationCache<anab::T0, sbn::crt::CRTHit> const& t0ToCRThits
) const {
  assert(!fCRTMatchProducer.empty());
  assert(trackPtr);
  
  sbn::selCRTInfo CRTinfo;
  
  //
  // fetch the CRT time associated to the track
  //
  art::Ptr<anab::T0> const t0Ptr = trackToCRTt0(trackPtr);
  if (!t0Ptr) return {}; // no CRT information associated with this track
  
  if (t0Ptr->ID() < 0) return {}; // supposedly invalid t0 content
  
  CRTinfo.time = t0Ptr->Time() / 1000.0; // nanoseconds -> microseconds
  
  //
  // fetch the CRT hit(s) information
  //
  
  // so far, we expect only one CRT hit;
  // when this changes, the cache one-to-one won't do any more...
  //   (then ask petrillo@slac.stanford.edu to write the one-to-many version)
  std::vector<sbn::crt::CRTHit const*> hits { t0ToCRThits(t0Ptr).get() };
  
  if (hits.empty()) {
    // we expect that the matching module between CRT hits and TPC tracks
    // produced a T0 and its associations to both track and hits;
    // if we find a T0 but we can't find the associated track or hit are missing,
    // we are doing something wrong
    throw art::Exception{ art::errors::LogicError }
      << "No CRT hit from '" << fCRTMatchProducer.encode() << "' associated to track #"
      << trackPtr.key() << " of '" << inputTagOf(trackPtr, event).encode() << "' (ID="
      << trackPtr->ID() << ") via anab::T0 #" << t0Ptr.key() << " from '"
      << inputTagOf(t0Ptr, event).encode() << "' (ID=" << t0Ptr->ID()
      << " time=" << t0Ptr->Time() << " ns)!\n";
  }
  
  // currently (v09_67_00) only one CRT hit per track is found; if we found more than one,
  // we should go back to the documentation or the authors to add the new information
  sbn::crt::CRTHit const& entryHit = *(hits.front());
  CRTinfo.entry = {
      entryHit.ts1_ns / 1000.0                                    // time [us]
    , entryHit.x_pos                                              // x
    , entryHit.y_pos                                              // y
    , entryHit.z_pos                                              // z
    , std::hypot(entryHit.x_err, entryHit.y_err, entryHit.z_err)  // resolution
    , static_cast<float>(t0Ptr->TriggerConfidence())              // DCA FIXME DCA is shoved into anab::T0 "trigger confidence" so far
    , icarus::crt::CRTCommonUtils{}.AuxDetRegionNameToNum(entryHit.tagger)
                                                                  // region
    };
  
  if (hits.size() > 1) {
    // this means that something has changed in the algorithm, and we should know what
    throw art::Exception{ art::errors::LogicError }
      << hits.size() << " CRT hits from '" << fCRTMatchProducer.encode()
      << "' associated to track #" << trackPtr.key() << " of '"
      << inputTagOf(trackPtr, event).encode() << "' (ID=" << trackPtr->ID()
      << ") via anab::T0 #" << t0Ptr.key() << " from '"
      << inputTagOf(t0Ptr, event).encode() << "' (ID=" << t0Ptr->ID()
      << " time=" << t0Ptr->Time() << "); expecting no more than 1.\n";
  }
  
  return CRTinfo;
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
geo::Vector_t sbn::TimeTrackTreeStorage::posShiftFromCRTtime(
  recob::Track const& track,
  double timeShift,
  std::vector<geo::TPCID> const& TPCs,
  geo::GeometryCore const& geom,
  detinfo::DetectorPropertiesData const& detProp
) const {
  if (TPCs.empty()) return {};
  
  auto iTPC = TPCs.begin();
  geo::Vector_t const driftDir = geom.TPC(*iTPC).DriftDir(); // toward anode
  
  // check for consistency of the drift direction
  while (++iTPC != TPCs.end()) {
    if (driftDir == geom.TPC(*iTPC).DriftDir()) continue;
    mf::LogPrint log(fLogCategory);
    log << "Warning: track ID=" << track.ID()
      << " (from " << track.Start() << " to " << track.End()
      << " cm) with no TPC T0 crosses TPCs with different drifts:";
    for (geo::TPCID const& tpcid: TPCs)
      log << " " << tpcid << " " << geom.TPC(tpcid).DriftDir();
    return {};
  } // while
  
  double const shiftTowardAnode = timeShift * detProp.DriftVelocity();
  
  return driftDir * shiftTowardAnode;
} // sbn::TimeTrackTreeStorage::posShiftFromCRTtime()


// -----------------------------------------------------------------------------
auto sbn::TimeTrackTreeStorage::distanceFromTimeRange
  (electronics_time time, lar::util::TrackTimeInterval::TimeRange const& range)
  -> microseconds
{
  using namespace util::quantities::time_literals;
  
  if (!range.isValid()) return 0_us; // all-inclusive range, perfect match
  
  // how much beyond the boundary the time is (negative = within boundary):
  microseconds const startDiff = range.start - time;
  microseconds const stopDiff = time - range.stop;
  
  if (range.duration() >= 0_us) {
    assert(startDiff <= 0_us || stopDiff <= 0_us);
    // if range contains time, both diffs are negative;
    //   we want to return the one with the smallest modulus, i.e. the largest
    // otherwise, if time is earlier than the range, startDiff is positive and
    //   stopDiff is negative, we want to return the first one, which is again
    //   the largest; and if time is after the range, the roles are inverted
    //   but we still want the largest (and only positive) one to be returned
    return std::max(startDiff, stopDiff);
  }
  else {
    assert(startDiff > 0_us || stopDiff > 0_us);
    /* It is impossible for this result to be negative, since the range is more
     * than empty.
     * If the time is earlier than the stop bound, then the exceeded boundary is
     * the start one, and we want to return the start difference; likewise if
     * the time is later than the start bound, then we want to return the stop
     * difference. In both case, the one we want to return is positive and the
     * other is negative, so we can return the highest value among them.
     * If the time is in between, both differences are positive and we want to
     * return the one from the closest boundary, that is the smallest one.
     */
    if ((startDiff > 0_us) && (stopDiff > 0_us))
      return std::min(startDiff, stopDiff);
    else
      return std::max(startDiff, stopDiff);
  }
  
} // distanceFromTimeRange()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::TimeTrackTreeStorage)


