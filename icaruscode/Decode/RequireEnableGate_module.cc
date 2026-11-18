/**
 * @file   icaruscode/Decode/RequireEnableGate_module.cc
 * @brief  Selects only events with a reasonable enable gate time.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 17, 2026
 *
 */

// ICARUS/LArSoft libraries
#include "icarusalg/Utilities/AtomicPassCounter.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "sbnobj/Common/Trigger/BeamBits.h" // sbn::bits::name()
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/intervals_fhicl.h" // nanoseconds from FHiCL

// framework libraries
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <cmath> // std::nan()
#include <cstdint> // std::int64_t
#include <string>


//------------------------------------------------------------------------------
using nanoseconds = util::quantities::intervals::nanoseconds;
using namespace util::quantities::time_literals;


//------------------------------------------------------------------------------
/**
 * @brief Requires the presence of an enable gate with a reasonable time.
 *
 * This filter accepts only events that have an enable gate time within a
 * configured time from the beam gate (and never after it).
 * 
 * The main customer is the data from ICARUS Run3 after April 30, in some events
 * of which the enable gate had not been opened, leaving the gate time
 * information to the one of the previous enable gate.
 * 
 * 
 * Configuration parameters
 * =========================
 *
 * * `TriggerTag` (input tag, default: `daqTrigger`): data product of the
 *   decoded trigger.
 * * `MaximumEnableGateDelay` (time, default: `1.1 ms`): maximum interval
 *   between the enable gate and the beam gate; note that this delay also covers
 *   any veto time where the beam gate is already opened but triggers are
 *   vetoed, and it is therefore larger than the customary 1 ms. It should also
 *   cover the largest possible jitter between the enable and beam gates.
 * * `LogCategory` (string, default: `RequireEnableGate`):
 *    message facility stream name where to write console messages to.
 *
 * 
 * Input data products
 * ====================
 * 
 * * `sbn::ExtraTriggerInfo` (from `TriggerTag`):
 *     the decoded hardware trigger information.
 * 
 * 
 * Service dependencies
 * =====================
 * 
 * Only message facility.
 * 
 * 
 * Output data products
 * =====================
 * 
 * None.
 * 
 */
class RequireEnableGate
  : public art::SharedFilter, icarus::ns::util::mfLoggingClass
{
  
    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> TriggerTag{
      Name{ "TriggerTag" },
      Comment{ "data product of the decoded trigger" },
      art::InputTag{ "daqTrigger" }
      };
    
    fhicl::Atom<nanoseconds> MaximumEnableGateDelay{
      Name{ "MaximumEnableGateDelay" },
      Comment{ "maximum interval between the enable gate and the beam gate" },
      1.1_ms
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "name of the message facility stream used by the module" },
      "RequireEnableGate"
      };
    
  }; // Config
  
  using Parameters = art::SharedFilter::Table<Config>;
  
  
  RequireEnableGate(Parameters const& params, art::ProcessingFrame const&);
  
  /// Evaluate the filtering logic.
  virtual bool filter(art::Event& event, art::ProcessingFrame const&) override;
  
  /// Prints a summary.
  virtual void endJob(art::ProcessingFrame const&) override;
  
    private:
  
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  art::InputTag const fTriggerTag; ///< Input trigger data product.
  
  nanoseconds const fMaximumEnableGateDelay; ///< Time from enable to beam gate.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  /// Count of passed and seen events.
  icarus::ns::util::AtomicPassCounter<> fPassed;
  
  
  /// Returns the difference between enable and beam gate times
  /// and whether it passes the requirements. Prints debug info.
  std::pair<nanoseconds, bool> testTrigger
    (sbn::ExtraTriggerInfo const& info) const;
  
  /// Returns the difference between two timestamps: `to - from`.
  static std::int64_t timestampDiff(std::uint64_t from, std::uint64_t to)
    { return static_cast<std::int64_t>(to) - static_cast<std::int64_t>(from); }
  
}; // class RequireEnableGate


//------------------------------------------------------------------------------
RequireEnableGate::RequireEnableGate
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedFilter{ params }
  , icarus::ns::util::mfLoggingClass{ params().LogCategory() }
  // configuration
  , fTriggerTag            { params().TriggerTag() }
  , fMaximumEnableGateDelay{ params().MaximumEnableGateDelay() }
{
  
  async<art::InEvent>();
  
  //
  // consume declaration
  //
  consumes<sbn::ExtraTriggerInfo>(fTriggerTag);
  
  //
  // dump configuration
  //
  auto log = mfLogInfo();
  log << "Configuration:"
    << "\n * gate times from trigger " << fTriggerTag.encode()
    << "\n * maximum time from enable to beam gate: "
      << util::quantities::intervals::milliseconds{ fMaximumEnableGateDelay }
    ;
  
} // RequireEnableGate::RequireEnableGate()


//------------------------------------------------------------------------------
auto RequireEnableGate::testTrigger(sbn::ExtraTriggerInfo const& info) const
  -> std::pair<nanoseconds, bool>
{
  
  auto log = mfLogDebug();
  log << "Gate #" << info.gateID << " (type: " << name(info.sourceType)
    << ") trigger #" << info.triggerID
    << " (type " << name(info.triggerType) << ")"
    ;
  
  if (!sbn::ExtraTriggerInfo::isValidTimestamp(info.enableGateTimestamp)) {
    log << " has invalid enable gate timestamp";
    return { nanoseconds{ std::nan("") }, false };
  }
  
  if (!sbn::ExtraTriggerInfo::isValidTimestamp(info.beamGateTimestamp)) {
    log << " has invalid beam gate timestamp";
    return { nanoseconds{ std::nan("") }, false };
  }
  
  auto const enableToBeam = nanoseconds::castFrom
    (timestampDiff(info.enableGateTimestamp, info.beamGateTimestamp));
  
  bool const passed
    = (enableToBeam >= 0_ns) && (enableToBeam <= fMaximumEnableGateDelay);
  
  log << " has " << enableToBeam
    << " between enable gate at " << info.enableGateTimestamp
    << " and beam gate at " << info.beamGateTimestamp;
  
  return { enableToBeam, passed };
} // RequireEnableGate::testTrigger()


//------------------------------------------------------------------------------
bool RequireEnableGate::filter
  (art::Event& event, art::ProcessingFrame const&)
{
  
  auto const& info = event.getProduct<sbn::ExtraTriggerInfo>(fTriggerTag);
  
  auto const [ enableToBeam, accepted ] = testTrigger(info);
  fPassed.add(accepted);
  
  if (accepted) {
    mfLogTrace() << "Enable gate is valid => ACCEPTED.";
  }
  else {
    mfLogInfo()
      << "Gate " << name(info.sourceType)
      << " trigger (type " << name(info.triggerType)
      << ") has " << enableToBeam
      << " between enable gate at " << info.enableGateTimestamp
      << " and beam gate at " << info.beamGateTimestamp << " => REJECTED"
      ;
  }
  
  return accepted;
  
} // RequireEnableGate::filter()


//------------------------------------------------------------------------------
void RequireEnableGate::endJob(art::ProcessingFrame const&) {
  
  if (fPassed.total() > 0) {
    mfLogInfo()
      << fPassed.failed() << "/" << fPassed.total()
      << " events (" << (fPassed.failed() * 100.0 / fPassed.total())
      << "%) were discarded.";
  }
  
} // RequireEnableGate::endJob()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(RequireEnableGate)


//------------------------------------------------------------------------------

