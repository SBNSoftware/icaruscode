/**
 * @file    icaruscode/Analysis/trigger/TimedTrackSelector_module.cc
 * @date    September 17, 2021
 * @authors Animesh Chatterjee (ANC238@pitt.edu),
 *          Gianluca Petrillo (petrillo@slac.stanford.edu),
 *          Jacob Zettlemoyer (jzettle@fnal.gov)
 * 
 */

#define MF_DEBUG

// SBN and ICARUS libraries
#ifdef USE_ATOMICPASSCOUNTER
#include "icarusalg/Utilities/AtomicPassCounter.h"
#else // !USE_ATOMICPASSCOUNTER
#include "icarusalg/Utilities/PassCounter.h"
#endif // USE_ATOMICPASSCOUNTER

// LArSoft libraries
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"


// framework libraries
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard library
#include <vector>
#include <string>
#include <atomic>
#include <memory>
#include <limits>

// -----------------------------------------------------------------------------
namespace sbn { class TimedTrackSelector; }
/**
 * @brief Selects tracks with time information.
 * 
 * This module produces a list of "tracks" that are associated to a time.
 * Optionally, it can select only the tracks with that time in a specified
 * range.
 * 
 * This module is a filter that will return `true` for an event if at least one
 * of its tracks is selected. This threshold can be adjusted by configuration.
 * 
 * 
 * Input
 * ------
 * 
 * The input file is expected to contain the time information and association
 * to "tracks", which are actually represented by `recob::PFParticle` rather
 * than `recob::Track`.
 * Note that the track data products are not explicitly required.
 * 
 * 
 * Output
 * -------
 * 
 * The filter _produces_ a list of pointers to the selected tracks; the list
 * will mix tracks from different data products if multiple input collections
 * are specified. A track is selected if all the following apply:
 * 
 * 1. the track is associated to a time;
 * 2. only if a time range is specified, that track time must fall within.
 * 
 * The filter passes the event if:
 *  * the number of selected tracks is within the configured range of requested
 *    tracks per event.
 * 
 * 
 * Configuration options
 * ----------------------
 * 
 * * `TrackTimeTags` (list of data product tags, required): data product of the
 *    times and their associations tracks.
 * * `MinT0` (real, optional): if specified, tracks are selected only if their
 *     associated time is not earlier than this value. Time is in the same time
 *     scale as the associated track time, which is expected to be the
 *     @ref DetectorClocksElectronicsStartTime "electronics time scale [us]".
 * * `MaxT0` (real, optional): if specified, tracks are selected only if their
 *     associated time is earlier than this value. Time is in the same time
 *     scale as the associated track time, which is expected to be the
 *     @ref DetectorClocksElectronicsStartTime "electronics time scale [us]".
 * * `MinTracks` (integer, default: `1`): the filter "passes" the event only if
 *     at least these many tracks are selected; disable this by setting it to
 *     `0`.
 * * `MaxTracks` (integer, default: a large number): if specified, filter
 *     "passes" the event only if at most these many tracks are selected.
 * 
 */
class sbn::TimedTrackSelector: public art::SharedFilter {
  
    public:
  
  static constexpr double NoMinTime { std::numeric_limits<double>::lowest() };
  static constexpr double NoMaxTime { std::numeric_limits<double>::max() };
  static constexpr unsigned int NoMinTracks { 0U };
  static constexpr unsigned int NoMaxTracks
    { std::numeric_limits<unsigned int>::max() };
  
  
  /// Module configuration parameters.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
    fhicl::Sequence<art::InputTag> TrackTimeTags {
      Name{ "TrackTimeTags" },
      Comment{ "data products with trach time information (and associations)" }
      };
    
    fhicl::Atom<double> MinT0 {
      Name{ "MinT0" },
      Comment{ "Select only tracks not earlier than this time [us]" },
      NoMinTime
      };
    
    fhicl::Atom<double> MaxT0 {
      Name{ "MaxT0" },
      Comment{ "Select only tracks earlier than this time [us]" },
      NoMaxTime
      };
    
    fhicl::Atom<unsigned int> MinTracks {
      Name{ "MinTracks" },
      Comment{ "Pass only events with at least these many selected tracks" },
      1U
      };
    
    fhicl::Atom<unsigned int> MaxTracks {
      Name{ "MaxTracks" },
      Comment{ "Pass only events with at most these many selected tracks" },
      NoMaxTracks
      };
    
    fhicl::Atom<bool> SaveTracks {
      Name{ "SaveTracks" },
      Comment{ "Whether to write the list of selected tracks" },
      true
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of a message facility stream for this module" },
      "TimedTrackSelector"
      };
    
  }; // Config
  
  
  using Parameters = art::SharedFilter::Table<Config>;
  
  explicit TimedTrackSelector
    (Parameters const& params, art::ProcessingFrame const&);
  
  bool filter(art::Event& event, art::ProcessingFrame const&) override;
  
  /// Prints end-job summary.
  void endJob(art::ProcessingFrame const&) override;
  

    private:

  // --- BEGIN -- Configuration parameters -------------------------------------
  
  /// List of track-time association input tags.
  std::vector<art::InputTag> const fTrackTimeTags;
  
  double const fMinT0;  ///< Minimum track time for track selection
  double const fMaxT0;  ///< Maximum track time for track selection
  
  /// Minimum selected tracks for event selection.
  unsigned int const fMinTracks;
  /// Maximum selected tracks for event selection.
  unsigned int const fMaxTracks;
  
  bool const fSaveTracks;  ///< Whether to save selected tracks into the event.
  
  std::string const fLogCategory; ///< Message facility stream name.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  /// Counter of passed events (not thread-safe).
#ifdef USE_ATOMICPASSCOUNTER
  icarus::ns::util::AtomicPassCounter<> fPassRate;
#else // !USE_ATOMICPASSCOUNTER
  icarus::ns::util::PassCounter<> fPassRate;
#endif // USE_ATOMICPASSCOUNTER
  
  
  /**
   * @brief Adds to `selectedTracks` qualifying tracks from `timeTracks`.
   * @param timeTracks time/track associations
   * @param[out] selectedTracks collection to expand with the qualifying tracks
   * @return the number of qualifying tracks found in `timeTracks` and added
   */
  unsigned int selectTracks(
    art::Assns<anab::T0, recob::PFParticle> const& timeTracks, std::vector<art::Ptr<recob::PFParticle>>& selectedTracks
    ) const;
  
  /// Returns whether the specified track (with specified time) qualifies.
  bool isTrackSelected
    (recob::PFParticle const& track, anab::T0 const& time) const;

  /// Returns if the number of tracks `nTracks` satisfies filter requirements.
  bool selectedTracksRequirement(unsigned int nTracks) const;
  

}; // sbn::TimedTrackSelector


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
sbn::TimedTrackSelector::TimedTrackSelector
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedFilter{ params }
  , fTrackTimeTags{ params().TrackTimeTags() }
  , fMinT0{ params().MinT0() }
  , fMaxT0{ params().MaxT0() }
  , fMinTracks{ params().MinTracks() }
  , fMaxTracks{ params().MaxTracks() }
  , fSaveTracks{ params().SaveTracks() }
  , fLogCategory{ params().LogCategory() }
{
  async<art::InEvent>();
  
  if (fSaveTracks)
    produces<std::vector<art::Ptr<recob::PFParticle>>>();
  
  //
  // configuration dump
  //
  {
    mf::LogInfo log{ fLogCategory };
    log << "Configuration:"
      << "\n  * tracks required to be associated to a time (cathode-crossers)"
      ;
    log << "\n  * track time:";
    if (fMinT0 == NoMinTime) {
      if (fMaxT0 == NoMaxTime) log << " any";
      else                     log << " before " << fMaxT0;
    }
    else {
      log << " from " << fMinT0;
      if (fMaxT0 == NoMaxTime) log << " on";
      else                     log << " to " << fMaxT0;
    }
    
    log << "\n  * selected tracks per event: ";
    if (fMinTracks == NoMinTracks) {
      if (fMaxTracks == NoMaxTracks) log << "any";
      else log << fMaxTracks << " or less";
    }
    else {
      if (fMaxTracks == NoMaxTracks) log << fMinTracks << " or more";
      else log << "between " << fMinTracks << " and " << fMaxTracks;
    }
    log << "\n  * selected tracks will" << (fSaveTracks? "": " not")
      << " be saved";
  } // end local scope
  
} // sbn::TimedTrackSelector::TimedTrackSelector()


// -----------------------------------------------------------------------------
bool sbn::TimedTrackSelector::filter
  (art::Event& event, art::ProcessingFrame const&)
{
  
  mf::LogDebug(fLogCategory) << "Processing " << event.id();
  
  //
  // select tracks from each of the input tags
  //
  std::vector<art::Ptr<recob::PFParticle>> selectedTracks;
  for (art::InputTag const& inputTag: fTrackTimeTags) {
    
    auto const& T0toTrack
      = event.getByLabel<art::Assns<anab::T0, recob::PFParticle>>(inputTag);
    
    unsigned int const newTracks = selectTracks(T0toTrack, selectedTracks);
    
    mf::LogTrace(fLogCategory)
      << "From '" << inputTag.encode() << "': "
      << newTracks << " tracks selected"
      ;
    
  } // for
  unsigned int const nSelectedTracks = selectedTracks.size();
  
  
  //
  // save track list in the event
  //
  
  // after this, selectedTracks may be empty
  if (fSaveTracks) {
    event.put(
      std::make_unique<std::vector<art::Ptr<recob::PFParticle>>>(selectedTracks)
      );
  }
  
  
  //
  // filter logic
  //
  
  bool const passed = selectedTracksRequirement(nSelectedTracks);
  fPassRate.add(passed);
  mf::LogTrace(fLogCategory) << event.id()
    << ' ' << (passed? "passed": "rejected") << " (" << nSelectedTracks
    << " selected tracks)."; // funny fact: we don't know the total track count
  
  
  mf::LogDebug(fLogCategory) << "Completed " << event.id();
  return passed;
  
} // sbn::TimedTrackSelector::filter()


// -----------------------------------------------------------------------------
void sbn::TimedTrackSelector::endJob(art::ProcessingFrame const&) {
  
  mf::LogInfo(fLogCategory) << "Selected " << fPassRate.passed()
    << '/' << fPassRate.total() << " events with qualifying tracks.";
  
} // sbn::TimedTrackSelector::endJob()


// -----------------------------------------------------------------------------
unsigned int sbn::TimedTrackSelector::selectTracks(
  art::Assns<anab::T0, recob::PFParticle> const& timeTracks, std::vector<art::Ptr<recob::PFParticle>>& selectedTracks
) const {
  
  unsigned int nSelectedTracks { 0U };
  for (auto const& [ t0Ptr, trackPtr ]: timeTracks) {
    
    // Q: why am I doing this?
    // A: currently, we don't need to read the track; but I leave it in the
    //    interface for the future; unless this is fixes, it will cause
    //    a segmentation violation at the first access.
    static recob::PFParticle const* dummyTrack { nullptr };
    
    if (!isTrackSelected(/* *trackPtr */ *dummyTrack, *t0Ptr)) continue;
    
    MF_LOG_TRACE(fLogCategory) << "Track #" << trackPtr.key() << " selected.";
    
    selectedTracks.push_back(trackPtr);
    ++nSelectedTracks;
    
  } // for
  
  return nSelectedTracks;
  
} // sbn::TimedTrackSelector::selectTracks()


// -----------------------------------------------------------------------------
bool sbn::TimedTrackSelector::isTrackSelected
  (recob::PFParticle const& track, anab::T0 const& time) const
{
  
  double const T0 = time.Time();
  MF_LOG_TRACE(fLogCategory) << "Track time: " << T0;
  
  if ((T0 < fMinT0) || (T0 >= fMaxT0)) {
    MF_LOG_TRACE(fLogCategory) << "Time out of range [ " << fMinT0 << "; "
      << fMaxT0 << " ] => discarded!";
    return false;
  }
  
  return true;
} // sbn::TimedTrackSelector::isTrackSelected()


// -----------------------------------------------------------------------------
bool sbn::TimedTrackSelector::selectedTracksRequirement
  (unsigned int nTracks) const
{
  return (nTracks >= fMinTracks) && (nTracks <= fMaxTracks);
} // sbn::TimedTrackSelector::selectedTracks()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::TimedTrackSelector)


// -----------------------------------------------------------------------------
