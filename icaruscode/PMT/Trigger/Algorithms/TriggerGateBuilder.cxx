/**
 * @file   icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.cxx
 * @brief  Algorithm to produce trigger gates out of optical readout waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h`
 * 
 */


// class header
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h"
#include "icarusalg/Utilities/WaveformOperations.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electromagnetism.h" // volt...
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <ostream>
#include <algorithm> // std::lower_bound(), std::transform()
#include <cctype> // std::isblank()
#include <iterator> // std::back_inserter()


//------------------------------------------------------------------------------
namespace {
  
  //----------------------------------------------------------------------------
  template <typename Coll>
  Coll sorted(Coll coll) {
    using std::begin, std::end;
    std::sort(begin(coll), end(coll));
    return coll;
  }

  //----------------------------------------------------------------------------
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::trigger::TriggerGateBuilder::TriggerGates
//------------------------------------------------------------------------------
auto icarus::trigger::TriggerGateBuilder::TriggerGates::findGateFor
  (raw::Channel_t const channel) const
{
  // keeping the gates sorted by channel
  auto const gend = fGates.end();
  auto const iGate
    = std::lower_bound(fGates.begin(), gend, channel, ChannelComparison<>());
  return ((iGate != gend) && (ChannelExtractor::channelOf(*iGate) == channel))
    ? iGate: gend;
} // icarus::trigger::TriggerGateBuilder::TriggerGates::findGateFor()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerGateBuilder::TriggerGates::findGateFor
  (raw::Channel_t const channel)
{
  return fGates.begin()
    + (std::as_const(*this).findGateFor(channel) - fGates.cbegin()); // weird
} // icarus::trigger::TriggerGateBuilder::TriggerGates::findGateFor()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerGateBuilder::TriggerGates::getGateFor
  (raw::Channel_t const channel) const -> triggergate_t const*
{
  auto const it = findGateFor(channel);
  return (it == fGates.end())? nullptr: &*it;
}

//------------------------------------------------------------------------------
auto icarus::trigger::TriggerGateBuilder::TriggerGates::gateFor
  (raw::OpDetWaveform const& waveform) -> triggergate_t&
{
  // keeping the gates sorted by channel
  raw::Channel_t const channel = waveform.ChannelNumber();
  auto iGate = findGateFor(channel);
  if (iGate != fGates.end()) {
    // found, it's there already
    MF_LOG_TRACE(details::TriggerGateDebugLog)
      << "Appending waveform to trigger gate (thr=" << threshold()
      << ") of channel " << channel;
    iGate->tracking().add(&waveform);
    return *iGate; 
  }
  // add the new gate in channel order
  MF_LOG_TRACE(details::TriggerGateDebugLog)
    << "Creating a new trigger gate (thr=" << threshold()
    << ") for channel " << channel;
  // wrap a new trigger gate on `channel` in a also new tracking gate
  iGate = fGates.emplace(iGate, triggergate_t::TriggerGate_t{ channel });
//   iGate->addChannel(channel);
  iGate->tracking().add(&waveform);
  return *iGate;
} // icarus::trigger::TriggerGateBuilder::TriggerGates::gateFor()


//------------------------------------------------------------------------------
//--- icarus::trigger::TriggerGateBuilder
//------------------------------------------------------------------------------
icarus::trigger::TriggerGateBuilder::TriggerGateBuilder(Config const& config)
  : fChannelThresholds(sorted(parseThresholds(config.ChannelThresholds())))
  , fPolarity{ config.Polarity() }
{
  
  //
  // configuration check
  //
  if (!fChannelThresholds.empty()
    && (fChannelThresholds.front() < ADCCounts_t{0})
  ) {
    throw cet::exception("TriggerGateBuilder")
     << "icarus::trigger::TriggerGateBuilder does not support"
        " negative thresholds (like "
     << fChannelThresholds.front() << ").\n";
  }
  
} // icarus::trigger::TriggerGateBuilder::TriggerGateBuilder()
  
  
//------------------------------------------------------------------------------
void icarus::trigger::TriggerGateBuilder::setup
  (detinfo::DetectorTimings const& timings)
{
  fDetTimings = timings;
} // icarus::trigger::TriggerGateBuilder::setup()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerGateBuilder::timeToOpticalTick
  (microsecond time) const -> optical_tick
{
  using optical_time = detinfo::timescales::optical_time;
  return detTimings().toTick<optical_tick>(optical_time{ time }); 
}


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerGateBuilder::timeStampToOpticalTick
  (raw::TimeStamp_t time) const -> optical_tick
{
  return timeToOpticalTick(microsecond{ time });
} // icarus::trigger::TriggerGateBuilder::timeStampToOpticalTick()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerGateBuilder::prepareAllGates() const
  -> std::vector<TriggerGates>
{

  // create an empty TriggerGates object for each threshold;
  // thresholds are kept relative
  std::vector<TriggerGates> allGates;
  allGates.reserve(nChannelThresholds());
  for (auto threshold: channelThresholds()) allGates.emplace_back(threshold);
  return allGates;
  
} // icarus::trigger::TriggerGateBuilder::prepareAllGates()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerGateBuilder::doDumpConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  dumpLocalConfiguration(out, indent, firstIndent);
  
} // icarus::trigger::TriggerGateBuilder::doDumpConfiguration()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerGateBuilder::dumpLocalConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  out << firstIndent << " * configured " << fChannelThresholds.size()
    << " discrimination thresholds:";
  for (ADCCounts_t const thr: fChannelThresholds)
    out << " " << thr;
  out << indent << "\n * signal polarity: "
    << util::StandardSelectorFor<util::SignalPolarity>{}.get(fPolarity).name();
  
} // icarus::trigger::TriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------
