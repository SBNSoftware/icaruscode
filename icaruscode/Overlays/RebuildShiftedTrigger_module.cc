/**
 * @file   RebuildShiftedTrigger_module.cc
 * @brief  Makeup module to produce shifted trigger data product.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 8, 2026
 * 
 */

// LArSoft and ICARUS libraries
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/TriggerData.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr
#include <string>
#include <utility> // std::move()
#include <vector>


/**
 * @brief Recreates a trigger data product after `AdjustSimForTrigger` shifting.
 * 
 * Old versions of module `AdjustSimForTrigger`, used to shift time reference
 * for some simulated data products, did not produce a new, updated trigger
 * data product that would be useful as time reference for the following stages.
 * And the new module still accommodates the inconsiderate choice of not
 * producing it.
 * 
 * This module tries to replicate the logic in the shifting module to reproduce
 * the appropriate trigger time and data product, and then produces them.
 * 
 * In addition to the configuration that was used in `AdjustSimForTrigger`, this
 * module also supports the specification of the reference time scale of the
 * trigger data product used as new reference, that in the original module was
 * taken from `DetectorClocksService`.
 * 
 * Note that if the input trigger is invalid, the output trigger will be still
 * given a valid time (same as the reference time).
 * 
 * 
 * Input
 * ------
 * 
 * * `std::vector<raw::Trigger>` (`InputTriggerLabel`): the trigger data product
 *   that was used as new time reference.
 * * `sbn::ExtraTriggerInfo` (`InputTriggerLabel`, if `SkipExtraTriggerInfo` is
 *   not set): auxiliary trigger data product that will also be produced.
 * * `std::vector<raw::Trigger>` (`TimeReferenceLabel`): the time reference
 *   the trigger in `InputTriggerLabel` was measured with respect to. This one
 *   is in electronics time scale.
 * 
 * 
 * Output
 * -------
 * 
 * * `std::vector<raw::Trigger>`: the shifted trigger data product, with the
 *   same bits as the original and both trigger and beam gate time shifted.
 * * `sbn::ExtraTriggerInfo` (if `SkipExtraTriggerInfo` is not set): the
 *   "shifted" version of the input `sbn::ExtraTriggerInfo` data product.
 *   Currently it is copied verbatim, just like `AdjustSimForTrigger` does.
 * 
 * Note that no `sim::BeamGateInfo` data product is currently produced.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * The following configuration parameters are supported:
 * 
 * * `InputTriggerLabel` (input tag, mandatory): tag of the trigger object used
 *   as new time reference by `AdjustSimForTrigger`.
 * * `TimeReferenceLabel` (input tag, optional): tag of the time reference the
 *   trigger in `InputTriggerLabel` is measured against. If not specified, the
 *   current trigger time from `DetectorClocksService` will be used (as it was
 *   in `AdjustSimForTrigger`).
 * * `AdditionalOffset` (real number, default: `0`): additional offset that was
 *   used in `AdjustSimForTrigger` configuration, in microseconds.
 * * `SkipExtraTriggerInfo` (flag, default: unset): if set, the auxiliary data
 *   product `sbn::ExtraTriggerInfo` won't be "shifted".
 * * `LogCategory` (text, default: `"RebuildShiftedTrigger"`): name of the
 *   message facility stream used for module messages on screen.
 * 
 * The parameters `InputTriggerLabel` and `AdditionalOffset` can be copied
 * verbatim from the original `AdjustSimForTrigger` configuration.
 * 
 */
class RebuildShiftedTrigger
  : public art::SharedProducer
  , private icarus::ns::util::mfLoggingClass
{
    public:
  
  /// FHiCL configuration interface.
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> InputTriggerLabel {
      Name{ "InputTriggerLabel" },
      Comment{ "tag of trigger data product used as new time reference" }
      };
    
    fhicl::OptionalAtom<art::InputTag> TimeReferenceLabel {
      Name{ "TimeReferenceLabel" },
      Comment
        { "tag of trigger data product that was reference to `InputTriggerLabel`" }
      };
    
    fhicl::Atom<double> AdditionalOffset {
      Name{ "AdditionalOffset" },
      Comment{ "additional shift that was added [us]" },
      0.0
      };
    
    fhicl::Atom<bool> SkipExtraTriggerInfo {
      Name{ "SkipExtraTriggerInfo" },
      Comment{ "do not shift the extra trigger information in `InputTriggerLabel`" },
      false
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of the message facility stream used by this module" },
      "RebuildShiftedTrigger"
      };
    
  }; // Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  RebuildShiftedTrigger(Parameters const& params, art::ProcessingFrame const&);
  
  void produce(art::Event& event, art::ProcessingFrame const&) override;
  
    private:
  
  // --- BEGIN ---  Configuration  ---------------------------------------------
  art::InputTag const fInputTriggerLabel; ///< Trigger tag for new ref.
  art::InputTag const fTimeReferenceLabel; ///< Trigger tag for new ref ref.
  double const fAdditionalOffset; ///< Additional offset on reference [&micro;s]
  bool const fSkipExtraTriggerInfo; ///< Whether to skip extra data production.
  // ---  END  ---  Configuration  ---------------------------------------------
  
  /**
   * @brief Computes the shift to be applied.
   * @param newRef the new time reference in electronics time
   * @param oldRef the old time reference in electronics time
   * @return a time shift [&micro;s] or `nullopt` if the reference is not valid
   */
  std::optional<double> computeShift(double newRef, double oldRef) const;
  
}; // RebuildShiftedTrigger


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
namespace {
  
  /// Moves the content of `obj` into an `std::unique_ptr` object.
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& obj)
    { return std::make_unique<T>(std::move(obj)); }
  
} // local namespace 


// -----------------------------------------------------------------------------
RebuildShiftedTrigger::RebuildShiftedTrigger
  (Parameters const& params, const art::ProcessingFrame&)
  : art::SharedProducer{ params }
  , icarus::ns::util::mfLoggingClass{ params().LogCategory() }
  // configuration parameters
  , fInputTriggerLabel{ params().InputTriggerLabel() }
  , fTimeReferenceLabel{ params().TimeReferenceLabel().value_or("") }
  , fAdditionalOffset{ params().AdditionalOffset() }
  , fSkipExtraTriggerInfo{ params().SkipExtraTriggerInfo() }
{
  
  async<art::InEvent>();
  
  //
  // input declaration
  //
  consumes<std::vector<raw::Trigger>>(fInputTriggerLabel);
  if (!fTimeReferenceLabel.empty())
    consumes<std::vector<raw::Trigger>>(fTimeReferenceLabel);
  if (!fSkipExtraTriggerInfo)
    consumes<sbn::ExtraTriggerInfo>(fInputTriggerLabel);
  
  //
  // output declaration
  //
  produces<std::vector<raw::Trigger>>();
  if (!fSkipExtraTriggerInfo) produces<sbn::ExtraTriggerInfo>();
  
  //
  // configuration report
  //
  {
    auto log = mfLogInfo();
    log
      <<   "Recreate a trigger data product for shift from '"
        << fInputTriggerLabel.encode() << "':"
      << "\n - use the trigger time from " << (fTimeReferenceLabel.empty()
        ? "DetectorClocksService": "'" + fTimeReferenceLabel.encode() + "'")
        << " as reference";
    if (fAdditionalOffset != 0)
      log << "\n - add a " << fAdditionalOffset << " us offset";
    log << "\n - " << (fSkipExtraTriggerInfo? "do not ": "")
      << "produce a sbn::ExtraTriggerInfo too";
    
  } // local scope
  
} // RebuildShiftedTrigger::RebuildShiftedTrigger()


// -----------------------------------------------------------------------------
void RebuildShiftedTrigger::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // input
  //
  
  auto const& triggers
    = event.getProduct<std::vector<raw::Trigger>>(fInputTriggerLabel);
  
  if (triggers.size() != 1) {
    if (triggers.empty()) {
      throw art::Exception(art::errors::EventProcessorFailure)
        << "No trigger '" << fInputTriggerLabel.encode() << "' available";
    }
    throw art::Exception(art::errors::EventProcessorFailure)
      << triggers.size() << " trigger objects in '"
      << fInputTriggerLabel.encode() << ", only one expected.";
  }
  
  const double oldReferenceTime = fTimeReferenceLabel.empty()
    ? art::ServiceHandle<detinfo::DetectorClocksService const>()
      ->DataFor(event).TriggerTime()
    : event.getProduct<std::vector<raw::Trigger>>(fTimeReferenceLabel)
      .front().TriggerTime()
    ;
  
  //
  // recompute shift
  // (based on the logic in AdjustSimForTrigger from sbncode v10_00_06)
  //
  
  raw::Trigger const& unshiftedTrigger = triggers.front();
  
  std::optional<double> const timeShift
    = computeShift(unshiftedTrigger.TriggerTime(), oldReferenceTime);
  if (timeShift) {
    mfLogInfo()
      << "Shifted trigger rebuilt applying a shift of " << *timeShift << " us";
  }
  else {
    mfLogInfo() << "Shifted trigger copied since new reference is invalid.";
  }
  
  //
  // output data products
  //
  
  std::vector<raw::Trigger> shiftedTriggers = {
    raw::Trigger{
      unshiftedTrigger.TriggerNumber(),
      timeShift? (unshiftedTrigger.TriggerTime() + *timeShift): oldReferenceTime,
        // trigger_time
      unshiftedTrigger.BeamGateTime() + timeShift.value_or(0.0),// beamgate_time
      unshiftedTrigger.TriggerBits()
    }
  };

  event.put(moveToUniquePtr(shiftedTriggers));
  
  
  // sbn::ExtraTriggerInfo (and raw::ExternalTrigger) are not shifted (so far);
  // it's debatable if they should be, since they only hold absolute timestamps
  if (!fSkipExtraTriggerInfo) {
    auto extraTrigger
      = event.getProduct<sbn::ExtraTriggerInfo>(fInputTriggerLabel);
    // no shift performed at this time
    event.put(moveToUniquePtr(extraTrigger));
  }
  
} // RebuildShiftedTrigger::produce()


// -----------------------------------------------------------------------------
std::optional<double> RebuildShiftedTrigger::computeShift
  (double newRef, double oldRef) const
{
  // the logic here is somehow questionable but replicates the original one:
  const bool hasValidTriggerTime =
    newRef >
      (std::numeric_limits<double>::min() + std::numeric_limits<double>::epsilon()) 
    && newRef <
      (std::numeric_limits<double>::max() - std::numeric_limits<double>::epsilon())
    ;
  if (!hasValidTriggerTime) return std::nullopt;
  
  const double newReferenceTime = newRef - fAdditionalOffset;
  const double timeShift = oldRef - newReferenceTime; // us
  return timeShift;
} // RebuildShiftedTrigger::computeShift()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(RebuildShiftedTrigger)


// -----------------------------------------------------------------------------

