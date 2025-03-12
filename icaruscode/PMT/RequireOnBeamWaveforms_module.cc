/**
 * @file   RequireOnBeamWaveforms_module.cc
 * @brief  Module requiring the presence of waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 31, 2023
 *
 */

// ICARUS/LArSoft libraries
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icarusalg/Utilities/AtomicPassCounter.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h" // DetectorClocksWithUnits
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds, ...
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"
// 
// C/C++ standard libraries
#include <algorithm> // std::count()
#include <optional>
#include <vector>
#include <string>
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
/**
 * @brief Requires the presence of one waveform per PMT covering beam gate time.
 *
 * This module passes only the events with one optical detector waveform
 * per channel covering the beam gate opening time.
 * 
 * It attempts to read optical detector metadata and, failing that, falls back
 * to the waveforms themselves.
 * 
 * Optionally, the availability of the waveforms is also verified
 * (`ForceWaveformPresence`). With this option, if a waveform is not available
 * it is considered like it were missing completely.
 * 
 * Note that the timestamp in the waveforms is assumed to be on the same scale
 * as the beam gate opening time, and, where it matters, expressed in
 * microseconds.
 * 
 *
 * Configuration parameters
 * =========================
 *
 * * `WaveformTag` (input tag, optional): tag of the optical detector waveforms;
 *   if not specified, optical detector waveforms are not read at all
 *   (unless `ForceWaveformPresence` is set).
 * * `WaveformMetaTag` (input tag, optional): tag of the optical detector
 *   waveform metadata; if empty, metadata is not used, and if not specified,
 *   the value of `WaveformTag` is used instead; in both special cases,
 *   `WaveformTag` must be specified.
 * * `TriggerTag` (input tag, optional): data product to use to learn when the
 *   beam gate was opened; if omitted, `DetectorClocksService` is used instead.
 * * `WaveformMetaAssnsTag` (input tag, default: like `WaveformMetaTag`): tag of
 *   the association between waveforms and their metadata. Normally, this is the
 *   same as the tag of the metadata itself, but in some old data that was not
 *   the case (`ReassociatePMTbaselines` would reassociate the existing one).
 * * `ForceWaveformPresence` (flag, default: `false`): if set, waveforms from
 *   `WaveformTag` are always read even when the check is performed with
 *   metadata; this detects the case where optical detector waveforms _were_
 *   present but have since been dropped.
 * * `ThrowOnFailure` (string, default: empty): if non-empty, instead of just
 *   returning `false`, when an event does not qualify, an exception is thrown
 *   with the category specified here, which can then be handled by _art_'s
 *   scheduler "service".
 * * `MissingChannels` (sequence of numbers, default: empty): the channels in
 *   this list are allowed not to be present.
 * * `LogCategory` (string, default: `RequireOnBeamWaveforms`): message facility
 *   stream name where to write console messages to.
 *
 * 
 * Service dependencies
 * =====================
 * 
 * * `Geometry`: to learn how many channels we are looking for.
 * * `DetectorClocksService`: depending on configuration, to get the beam
 *    opening time and for optical tick duration.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::OpDetWaveform>` (`WaveformTag`, if specified): the
 *   waveforms which must have been present (and, if required, available).
 * * `std::vector<sbn::OpDetWaveformMeta>` (`WaveformTagMeta`, if specified):
 *   the metadata of the waveforms which must be present.
 * * `std::vector<raw::Trigger>` (`TriggerTag`, if specified) to read the beam
 *   gate opening time.
 * * `art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>` (`WaveformTagMeta`
 *   if `ForceWaveformPresence` is set and metadata is used) to verify the
 *   availability of waveforms.
 * 
 * 
 * Output data products
 * =====================
 * 
 * None.
 * 
 */
class RequireOnBeamWaveforms: public art::SharedFilter {
  
    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::OptionalAtom<art::InputTag> WaveformTag{
      Name{ "WaveformTag" },
      Comment{ "optical detector waveform data product to be checked" }
      };
    
    fhicl::OptionalAtom<art::InputTag> WaveformMetaTag{
      Name{ "WaveformMetaTag" },
      Comment{ "optical detector waveform metadata product to be checked" }
      };
    
    fhicl::Atom<art::InputTag> TriggerTag{
      Name{ "TriggerTag" },
      Comment{ "trigger data product to use to discover the beam gate time" },
      art::InputTag{}
      };
    
    fhicl::OptionalAtom<art::InputTag> WaveformMetaAssnsTag{
      Name{ "WaveformMetaAssnsTag" },
      Comment{ "optical detector waveform and metadata associations" }
      };
    
    fhicl::Atom<bool> ForceWaveformPresence{
      Name{ "ForceWaveformPresence" },
      Comment{ "ensures waveforms are available in the input file" },
      false
      };
    
    fhicl::Sequence<raw::Channel_t> MissingChannels{
      Name{ "MissingChannels" },
      Comment{ "" },
      std::vector<raw::Channel_t>{}
      };
    
    fhicl::Atom<std::string> ThrowOnFailure{
      Name{ "ThrowOnFailure" },
      Comment{ "if criteria are not met throw exception with this category" },
      ""
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "name of the message facility stream used by the module" },
      "RequireOnBeamWaveforms"
      };
    
  }; // Config
  
  using Parameters = art::SharedFilter::Table<Config>;
  
  
  RequireOnBeamWaveforms
    (Parameters const& params, art::ProcessingFrame const&);
  
  /// Evaluate the filtering logic.
  virtual bool filter(art::Event& event, art::ProcessingFrame const&) override;
  
  /// Prints a summary.
  virtual void endJob(art::ProcessingFrame const&) override;
  
    private:
  
  using microsecond = util::quantities::points::microsecond;
  using nanoseconds = util::quantities::intervals::nanoseconds;
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  art::InputTag const fWaveformTag; ///< Input waveform data product.
  
  art::InputTag const fWaveformMetaTag; ///< Input waveform metadata product.
  
  art::InputTag const fTriggerTag; ///< Data product for beam gate.
  
  art::InputTag const fWaveformMetaAssnsTag; ///< Input waveform/metadata assns.
  
  bool const fForceWaveformPresence; ///< Whether to always read waveforms.
  
  /// Category of the exception to throw on failure (empty: does not throw).
  std::string const fThrowOnFailure;
  
  std::string const fLogCategory; ///< Name of message facility stream.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  // --- BEGIN -- Caches -------------------------------------------------------
  
  /// Bits of required channels (index is channel number).
  std::vector<bool> const fRequiredChannelMask;
  
  std::ptrdiff_t const fNRequired; ///< Number of channels needed on beam.
  
  // --- END ---- Caches -------------------------------------------------------
  
  /// Count of passed and seen events.
  icarus::ns::util::AtomicPassCounter<> fPassed;
  
  /// Reads the beam gate opening time from `DetectorClocksService`.
  microsecond beamGateTimeFromService(art::Event const& event) const;
  
  /// Reads the beam gate opening time from the configured trigger data product.
  microsecond beamGateTimeFromData(art::Event const& event) const;

  /// Returns a mask with positions of the missing channels set to `true`.
  std::vector<bool> missingFromMetadata(
    microsecond beamGateTime,
    std::vector<sbn::OpDetWaveformMeta> const& metadata,
    art::FindOneP<raw::OpDetWaveform> const* toWaveforms = nullptr
    ) const;


  /// Returns a mask with positions of the missing channels set to `true`.
  std::vector<bool> missingFromWaveforms(
    microsecond beamGateTime,
    std::vector<raw::OpDetWaveform> const& waveforms,
    nanoseconds clockTick
    ) const;
  
  /// Composes and returns the list of channels to require.
  static std::vector<bool> buildRequiredChannelList
    (std::vector<raw::Channel_t> const& skipChannels = {});
  
  
}; // class RequireOnBeamWaveforms


//------------------------------------------------------------------------------

RequireOnBeamWaveforms::RequireOnBeamWaveforms
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedFilter{ params }
  // configuration
  , fWaveformTag          { params().WaveformTag().value_or(art::InputTag{}) }
  , fWaveformMetaTag      { params().WaveformMetaTag().value_or(fWaveformTag) }
  , fTriggerTag           { params().TriggerTag() }
  , fWaveformMetaAssnsTag
    { params().WaveformMetaAssnsTag().value_or(fWaveformMetaTag) }
  , fForceWaveformPresence{ params().ForceWaveformPresence() }
  , fThrowOnFailure       { params().ThrowOnFailure() }
  , fLogCategory          { params().LogCategory() }
  // caches
  , fRequiredChannelMask{ buildRequiredChannelList(params().MissingChannels()) }
  , fNRequired{
      std::count
        (fRequiredChannelMask.begin(), fRequiredChannelMask.end(), true)
      }
{
  
  async<art::InEvent>();
  
  //
  // configuration check
  //
  if (fWaveformTag.empty() && fWaveformMetaTag.empty()) {
    throw art::Exception{ art::errors::Configuration }
      << "At least one of '" << params().WaveformTag.name() << "' and '"
      << params().WaveformMetaTag.name() << "' must be specified non-empty.\n";
  }
  
  
  //
  // consume declaration
  //
  if (!fWaveformMetaTag.empty()) {
    consumes<std::vector<sbn::OpDetWaveformMeta>>(fWaveformMetaTag);
    if (!fWaveformTag.empty())
      mayConsume<std::vector<raw::OpDetWaveform>>(fWaveformTag);
    if (fForceWaveformPresence) {
      consumes<art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>>
        (fWaveformMetaTag);
    }
  }
  else {
    assert(!fWaveformTag.empty());
    consumes<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  }
  
  if (!fTriggerTag.empty()) consumes<std::vector<raw::Trigger>>(fTriggerTag);
  
  
  //
  // dump configuration
  //
  mf::LogInfo log{ fLogCategory };
  log << "Configuration:";
  if (!fWaveformMetaTag.empty()) {
    log << "\n * check performed on waveform metadata '"
      << fWaveformMetaTag.encode() << "'";
    if (!fWaveformTag.empty()) {
      log << "\n   * if not available, check is performed on waveforms '"
        << fWaveformTag.encode() << "'";
    }
  }
  else {
    log << "\n * check performed on waveforms '"
      << fWaveformTag.encode() << "'";
  }
  if (fForceWaveformPresence) {
    log
      << "\n * presence of target waveforms checked (via metadata association '"
      << fWaveformMetaAssnsTag.encode() << "')";
  }
  if (fTriggerTag.empty()) {
    log << "\n * beam gate time from DetectorClocksService";
  }
  else {
    log << "\n * beam gate time from trigger data product '"
      << fTriggerTag.encode() << "'";
  }
  log << "\n * requiring the presence of " << fNRequired << " channels";
  if (!fThrowOnFailure.empty())
    log << "\n * if an event fails the criteria, throw an exception of category '" << fThrowOnFailure << "'";
  
  
} // RequireOnBeamWaveforms::RequireOnBeamWaveforms()


//------------------------------------------------------------------------------
std::vector<bool> RequireOnBeamWaveforms::buildRequiredChannelList
  (std::vector<raw::Channel_t> const& skipChannels /* = {} */)
{
  std::vector<bool> mask
    (lar::providerFrom<geo::Geometry>()->MaxOpChannel(), true);
  
  // they will be never missed
  for (raw::Channel_t const channel: skipChannels) mask[channel] = false;
  
  return mask;
} // RequireOnBeamWaveforms::buildRequiredChannelList()


//------------------------------------------------------------------------------
bool RequireOnBeamWaveforms::filter
  (art::Event& event, art::ProcessingFrame const&)
{
  
  std::vector<std::size_t> WaveformIndices;
  std::vector<bool> missingWaveforms;
  
  //
  // discover the time
  //
  microsecond const beamGateTime = fTriggerTag.empty()
    ? beamGateTimeFromService(event): beamGateTimeFromData(event);
  
  if (!fWaveformMetaTag.empty()) {
    auto waveformMetadata
      = event.getHandle<std::vector<sbn::OpDetWaveformMeta>>(fWaveformMetaTag);
    if (waveformMetadata) {
      
      auto const toWaveforms = fForceWaveformPresence
        ? std::make_optional<art::FindOneP<raw::OpDetWaveform>>
          (waveformMetadata, event, fWaveformMetaAssnsTag)
        : std::nullopt
        ;
      
      missingWaveforms = missingFromMetadata
        (beamGateTime, *waveformMetadata, toWaveforms? &*toWaveforms: nullptr);
      
    }
    else {
      mf::LogTrace{ fLogCategory } << "Waveform metadata '"
        << fWaveformMetaTag.encode() << "' not available.";
    }
  }
  
  if (missingWaveforms.empty()) {
    if (fWaveformTag.empty()) {
      throw art::Exception{ art::errors::ProductNotFound }
        << "Waveform metadata '" << fWaveformMetaTag.encode()
        << "' not available and no waveform tag was configured!\n";
    }
    
    detinfo::DetectorClocksWithUnits const detClocks{
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
      };
    
    
    missingWaveforms = missingFromWaveforms(
      beamGateTime,
      event.getProduct<std::vector<raw::OpDetWaveform>>(fWaveformTag),
      detClocks.OpticalClockPeriod()
      );
  }
  
  
  //
  // draw conclusions
  //
  unsigned int const nMissing
    = std::count(missingWaveforms.begin(), missingWaveforms.end(), true);
  bool const accepted = (nMissing == 0);
  
  fPassed.add(accepted);
  
  if (!accepted) {
    mf::LogInfo{ fLogCategory }
      << event.id() << ": " << nMissing << "/" << fNRequired << " '"
      << fWaveformTag.encode()
      << "' optical detector waveforms did not meet the on-beam requirements.";
    
    if (!fThrowOnFailure.empty()) {
      throw cet::exception{ fThrowOnFailure }
        << event.id() << ": " << nMissing << "/" << fNRequired << " '"
      << fWaveformTag.encode()
      << "' optical detector waveforms did not meet the on-beam requirements\n";
    }
  }
  
  return accepted;
  
} // RequireOnBeamWaveforms::filter()


//------------------------------------------------------------------------------
void RequireOnBeamWaveforms::endJob(art::ProcessingFrame const&) {
  
  if (fPassed.failed() > 0) {
    mf::LogInfo{ fLogCategory }
      << fPassed.failed() << "/" << fPassed.total()
      << " events failed the requirements.";
  }
  
} // RequireOnBeamWaveforms::endJob()


//------------------------------------------------------------------------------
auto RequireOnBeamWaveforms::beamGateTimeFromService
  (art::Event const& event) const -> microsecond
{
  detinfo::DetectorClocksWithUnits const detClocks{
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  return detClocks.BeamGateTime();
} // RequireOnBeamWaveforms::beamGateTimeFromService()


//------------------------------------------------------------------------------
auto RequireOnBeamWaveforms::beamGateTimeFromData(art::Event const& event) const
  -> microsecond
{
  assert(!fTriggerTag.empty());
  auto const& triggers
    = event.getProduct<std::vector<raw::Trigger>>(fTriggerTag);
  if (triggers.empty()) {
    throw cet::exception("RequireOnBeamWaveforms")
      << "No triggers in '" << fTriggerTag.encode() << "'!\n";
  }
  return microsecond::castFrom(triggers.front().BeamGateTime());
} // RequireOnBeamWaveforms::beamGateTimeFromData()


//------------------------------------------------------------------------------
std::vector<bool> RequireOnBeamWaveforms::missingFromMetadata(
  microsecond beamGateTime,
  std::vector<sbn::OpDetWaveformMeta> const& metadata,
  art::FindOneP<raw::OpDetWaveform> const* toWaveforms /* = nullptr */
) const {
  
  std::vector<bool> missingWaveforms{ fRequiredChannelMask };
  
  for (auto const& [ iWaveform, waveformMeta ]: util::enumerate(metadata)) {
    
    // we could also check the on-beam flag here (`withBeamGate()`),
    // but it may not be reliable for the original simulation waveforms
    // (when metadata is created, beam gate time might be not known yet)
    
    if (!waveformMeta.includes(beamGateTime.value())) continue;
    
    if (fForceWaveformPresence) {
      assert(toWaveforms);
      
      art::Ptr<raw::OpDetWaveform> const waveform = toWaveforms->at(iWaveform);
      
      if (waveform.isNull()) {
        throw art::Exception{ art::errors::LogicError }
          << "Metadata #" << iWaveform
          << " is not associated to any waveform!\n";
      }
      if (!waveform.isAvailable()) continue; // like if we hadn't found it
    }
    
    // found
    if (waveformMeta.channel < missingWaveforms.size())
      missingWaveforms[waveformMeta.channel] = false;
    
  } // for
  
  return missingWaveforms;
} // RequireOnBeamWaveforms::missingFromMetadata()


//------------------------------------------------------------------------------
std::vector<bool> RequireOnBeamWaveforms::missingFromWaveforms(
  microsecond beamGateTime, std::vector<raw::OpDetWaveform> const& waveforms,
  nanoseconds clockTick
) const {
  
  std::vector<bool> missingWaveforms{ fRequiredChannelMask };
  
  for (raw::OpDetWaveform const& waveform: waveforms) {
    
    // check: contains beam gate
    microsecond const startTime{ waveform.TimeStamp() };
    microsecond const endTime{ startTime + clockTick * waveform.size() };
    if ((beamGateTime < startTime) || (beamGateTime >= endTime)) continue;
    
    // that was enough
    if (waveform.ChannelNumber() < missingWaveforms.size())
      missingWaveforms[waveform.ChannelNumber()] = false;
    
  } // for
  
  return missingWaveforms;
} // RequireOnBeamWaveforms::missingFromWaveforms()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(RequireOnBeamWaveforms)


//------------------------------------------------------------------------------

