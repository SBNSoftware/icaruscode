/**
 * @file   CopyBeamTimePMTwaveforms_module.cc
 * @brief  Selects `raw::OpDetWaveform` with a specific time.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 1, 2022
 */


// ICARUS libraries
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <memory> // std::make_unique()
#include <utility> // std::move()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus { class CopyBeamTimePMTwaveforms; }

/**
 * @brief Creates a new list of PMT waveforms covering the specified time.
 * 
 * This module creates a new data product with all and only the PMT waveforms
 * from an existing data product that cover the configured time point.
 * Some associations are replicated: a new association is produced from the
 * elements of the new data product to the old elements associated to the
 * original PMT waveform data product; the associated data is not replicated.
 * 
 * The replicated associations currently include:
 * 
 * * `sbn::OpDetWaveformMeta` (configuration: `OpDetWaveformMetaAssns`)
 * 
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::OpDetWaveform>` (tag from `Waveforms`):
 *   all optical detector waveforms to apply the selection on
 * * `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>`
 *   (configuration: `OpDetWaveformMetaAssns`): association between the original
 *   waveforms and their metadata
 * 
 * 
 * Output data products
 * =====================
 *
 * * `std::vector<raw::OpDetWaveform>`: a copy of the selected waveforms from
 *   the original collection (`Waveforms` configuration parameter).
 * * `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>` (optional):
 *   association of the copied waveforms to their existing metadata
 *   (`OpDetWaveformMetaAssns` configuration parameter)
 * 
 * 
 * Service dependences
 * ====================
 * 
 * This module requires:
 * 
 * * `DetectorClocksService`: we learn the trigger and beam time from here, and
 *   use it (and its associate, `detinfo::DetectorTimings`) to convert times.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse online description of the parameters is printed by running
 * `lar --print-description CopyBeamTimePMTwaveforms`.
 * 
 * * `Waveforms` (input tag, mandatory): the list of optical detector waveforms
 *   to be processed.
 * * `SelectTime` (time point, default: `0 ns`): the time that the waveforms
 *   need to cover for them to be selected. The time is by default referring to
 *   the beam time (from `TriggerTag`), but a different reference time can be
 *   specified (`TimeReference` configuration parameter).
 *   The time requires a unit.
 * * `TimeReference` (one of: `BeamGate` (default), `electronics`) time scale
 *   the selected time lies on:
 *    * `electronics`:
 *      @ref DetectorClocksElectronicsStartTime "electronics time scale".
 *    * `BeamGate`: time of the opening of the beam gate; beam time is learnt
 *      from `DetectorClocksService`.
 * * association parameters:
 *     * `OpDetWaveformMetaAssns` (input tag, optional): the tag of the
 *       associations between the original data product above and its metadata.
 *       If specified empty, no metadata association is created. If the
 *       parameter is omitted, instead, the associated metadata is found with
 *       the same tag as `Waveforms`.
 * * `LogCategory` (string, default: `CopyBeamTimePMTwaveforms`): name of the
 *     output stream category for console messages (managed by MessageFacility
 *     library).
 * 
 */
class icarus::CopyBeamTimePMTwaveforms: public art::SharedProducer {
  
  enum class TimeReference_t {
      kElectronicsTime
    , kBeamGate
    , kTrigger
    , kDefault = kBeamGate
  }; // TimeReference_t
  
    public:
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using microseconds = util::quantities::intervals::microseconds; // alias
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> Waveforms {
      Name{ "Waveforms" },
      Comment{ "tag of input optical detector waveforms" }
      // mandatory
      };
    
    fhicl::Atom<microseconds> SelectTime {
      Name{ "SelectTime" },
      Comment{ "time the selected waveforms need to include" },
      microseconds{ 0 }
      };
    
    fhicl::Atom<std::string> TimeReference {
      Name{ "TimeReference" },
      Comment{ "time scale the selection time refers to: "
        + TimeReferenceSelector.optionListString() },
      TimeReferenceSelector.get(TimeReference_t::kDefault).name()
      };
    
    fhicl::OptionalAtom<art::InputTag> OpDetWaveformMetaAssns {
      Name{ "OpDetWaveformMetaAssns" },
      Comment{ "tag of waveform metadata association"
        " (default as Waveforms, empty to skip)" 
        }
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of the category used for the output" },
      "CopyBeamTimePMTwaveforms" // default
      };
    
    
    /// Selector for `TimeReference` parameter.
    static util::MultipleChoiceSelection<TimeReference_t> const
      TimeReferenceSelector;
    
    TimeReference_t getTimeReference() const; ///< Returns `TimeReference`.
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit CopyBeamTimePMTwaveforms
    (Parameters const& config, art::ProcessingFrame const&);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // aliases
  using nanoseconds = util::quantities::intervals::nanoseconds;
  using electronics_time = detinfo::timescales::electronics_time;
  
  
  // --- BEGIN -- Configuration variables --------------------------------------
  
  art::InputTag const fWaveformTag; ///< Input waveforms.
  art::InputTag const fWaveMetaTag; ///< Input waveform metadata associations.
  
  nanoseconds const fTargetTime; ///< Time the selected waveforms contain.
  TimeReference_t const fTimeReference; ///< Reference for the target time.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END ---- Configuration variables --------------------------------------
  
  
  // --- BEGIN -- Cached quantities --------------------------------------------
  
  nanoseconds const fOpticalTick; ///< PMT digitiser sampling period.
  
  // --- END ---- Cached quantities --------------------------------------------
  
  
  /// Returns the configured reference value of the selection time.
  electronics_time getReferenceTime
    (detinfo::DetectorTimings const& detTimings) const;
  
  /// Returns whether `waveform` includes `time`.
  bool contains
    (electronics_time time, raw::OpDetWaveform const& waveform) const;
  
  
  /// Returns whether the configuration requested copying metadata associations.
  bool doWaveMetaAssns() const { return !fWaveMetaTag.empty(); }
  
}; // class icarus::CopyBeamTimePMTwaveforms


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& data)
    { return std::make_unique<T>(std::move(data)); }

} // local namespace


//------------------------------------------------------------------------------
//--- icarus::CopyBeamTimePMTwaveforms
//------------------------------------------------------------------------------
util::MultipleChoiceSelection<icarus::CopyBeamTimePMTwaveforms::TimeReference_t>
  const
icarus::CopyBeamTimePMTwaveforms::Config::TimeReferenceSelector{
    { TimeReference_t::kElectronicsTime, "Electronics", "ElectronicsTime" }
  , { TimeReference_t::kBeamGate,        "BeamGate", "Beam", "BeamGateTime" }
  , { TimeReference_t::kTrigger,         "Trigger", "TriggerTime" }
  };


//------------------------------------------------------------------------------
auto icarus::CopyBeamTimePMTwaveforms::Config::getTimeReference() const
  -> TimeReference_t
{
  try {
    return static_cast<TimeReference_t>
      (TimeReferenceSelector.parse(TimeReference()).value());
  }
  catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e)
  {
    throw art::Exception(art::errors::Configuration)
      << "Invalid value for '" << TimeReference.name()
      << "' parameter: '" << e.label() << "'; valid options: "
      << TimeReferenceSelector.optionListString() << ".\n";
  }
} // icarus::CopyBeamTimePMTwaveforms::Config::getTimeReference()


//------------------------------------------------------------------------------
icarus::CopyBeamTimePMTwaveforms::CopyBeamTimePMTwaveforms
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedProducer(config)
  // configuration
  , fWaveformTag  { config().Waveforms() }
  , fWaveMetaTag  { config().OpDetWaveformMetaAssns().value_or(fWaveformTag) }
  , fTargetTime   { config().SelectTime() }
  , fTimeReference{ config().getTimeReference() }
  , fLogCategory  { config().LogCategory() }
  // cached
  , fOpticalTick{
    detinfo::makeDetectorTimings(
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob()
      ).OpticalClockPeriod()
    }
{
  
  async<art::InEvent>();
  
  //
  // input data declaration
  //
  consumes<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  if (doWaveMetaAssns()) {
    consumes<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
      (fWaveMetaTag);
  }
  
  
  //
  // output data declaration
  //
  produces<std::vector<raw::OpDetWaveform>>();
  if (doWaveMetaAssns()) {
    produces<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>();
  }
  
  
  //
  // configuration report (short)
  //
  
  mf::LogInfo log{ fLogCategory };
  log
    << "Configuration:"
    << "\n - input waveforms: '" << fWaveformTag.encode() << '\''
    ;
  if (doWaveMetaAssns()) {
    log << "\n - waveform metadata: associated in '" << fWaveMetaTag.encode()
      << "'";
  }
  log << "\n - selection time:";
  if (fTargetTime != nanoseconds{ 0.0 }) log << " " << fTargetTime << " from";
  log << Config::TimeReferenceSelector.get(fTimeReference).name();
  
} // icarus::CopyBeamTimePMTwaveforms::CopyBeamTimePMTwaveforms()


//------------------------------------------------------------------------------
void icarus::CopyBeamTimePMTwaveforms::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // get the timing information for this event
  //
  detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings
    (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event))
    ;
  assert(detTimings.OpticalClockPeriod() == fOpticalTick);
  
  //
  // fetch input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fWaveformTag);

  art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> const* waveMetaAssns
    = doWaveMetaAssns()
    ? &(event.getProduct<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
      (fWaveMetaTag))
    : nullptr
    ;
  
  if (waveMetaAssns && (waveMetaAssns->size() != waveformHandle->size())) {
    // the assumption: the associations are in the same order as the
    // original waveforms, and none is missing; this allows us to skip
    // art::FindOneP calls. If not true... art::FindOneP is a way.
    throw art::Exception{ art::errors::LogicError }
      << "CopyBeamTimePMTwaveforms association logic assumption is broken (I)."
      "\nPlease contact the author for a fix.\n";
  }
  
  //
  // determine the target time
  //
  
  {
    electronics_time const triggerTime = detTimings.TriggerTime();
    electronics_time const beamGateTime = detTimings.BeamGateTime();
    mf::LogDebug(fLogCategory)
      << "Event " << event.id() << " has beam gate starting at " << beamGateTime
      << " and trigger at " << triggerTime << "."
      << "\nNow extracting information from " << waveformHandle->size()
        << " waveforms."
      ;
  }
  
  electronics_time const targetTime
    = getReferenceTime(detTimings) + fTargetTime;
  mf::LogDebug{ fLogCategory } << "Target time: " << targetTime;
  
  //
  // selection and copies
  //
  std::vector<raw::OpDetWaveform> selWaveforms;
  std::unique_ptr<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
    selWaveMetaAssns
    = doWaveMetaAssns()
    ? std::make_unique<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>()
    : nullptr
    ;
  
  art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr
    { event, waveformHandle.id() };
  
  for (auto const& [ iWaveform, waveform ]: util::enumerate(*waveformHandle)) {
    
    if (!contains(targetTime, waveform)) continue;
    
    mf::LogTrace{ fLogCategory } << "Waveform #" << iWaveform
      << " selected (channel " << waveform.ChannelNumber() << ", time "
      << waveform.TimeStamp() << " us, duration: "
      << (waveform.size() * fOpticalTick) << ")";
    
    selWaveforms.push_back(waveform);
    
    //
    // associations
    //
    if (doWaveMetaAssns()) {
      assert(selWaveMetaAssns);
      
      auto waveAssn = waveMetaAssns->at(iWaveform);
      if (waveAssn.second.key() != iWaveform) {
        // the assumption: the associations are in the same order as the
        // original waveforms, and none is missing; this allows us to skip
        // art::FindOneP calls. If not true... art::FindOneP is a way.
        throw art::Exception{ art::errors::LogicError }
          << "CopyBeamTimePMTwaveforms association logic assumption is broken"
            " (II)."
            "\nPlease contact the author for a fix.\n";
      }
      
      art::Ptr<raw::OpDetWaveform> const wavePtr = makeWaveformPtr(iWaveform);
      
      selWaveMetaAssns->addSingle(wavePtr, waveAssn.second);
      
    } // if waveform metadata association
    
  } // for waveforms
  
  
  //
  // store output
  //
  event.put(moveToUniquePtr(selWaveforms));
  if (doWaveMetaAssns()) event.put(std::move(selWaveMetaAssns));
  
} // icarus::CopyBeamTimePMTwaveforms::produce()


//------------------------------------------------------------------------------
auto icarus::CopyBeamTimePMTwaveforms::getReferenceTime
  (detinfo::DetectorTimings const& detTimings) const -> electronics_time
{
  
  switch (fTimeReference) {
    
    case TimeReference_t::kElectronicsTime: return electronics_time{ 0.0 };
    case TimeReference_t::kBeamGate:        return detTimings.BeamGateTime();
    case TimeReference_t::kTrigger:         return detTimings.TriggerTime();
    
  } // switch(fTimeReference)
  
  // compiler should issue an error if a case of a enum class is not handled...
  throw art::Exception{ art::errors::LogicError }
    << "icarus::CopyBeamTimePMTwaveforms::getReferenceTime(): case '"
    << Config::TimeReferenceSelector.get(fTimeReference).name()
    << "' did not get handled, and the compiler missed it.\n";
  
} // icarus::CopyBeamTimePMTwaveforms::getReferenceTime()


//------------------------------------------------------------------------------
bool icarus::CopyBeamTimePMTwaveforms::contains
  (electronics_time time, raw::OpDetWaveform const& waveform) const
{
  
  using util::quantities::points::microsecond;
  
  electronics_time const startWaveformTime
    { microsecond{ waveform.TimeStamp() } };
  electronics_time const endWaveformTime
    = startWaveformTime + fOpticalTick * waveform.size();
  
  return (time >= startWaveformTime) && (time < endWaveformTime);
  
} // icarus::CopyBeamTimePMTwaveforms::contains()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::CopyBeamTimePMTwaveforms)


//------------------------------------------------------------------------------
