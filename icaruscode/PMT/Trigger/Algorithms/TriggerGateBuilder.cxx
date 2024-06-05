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
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icarusalg/Utilities/WaveformOperations.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <algorithm> // std::lower_bound(), std::transform()
#include <iterator> // std::back_inserter()


//------------------------------------------------------------------------------
namespace {
  
  //----------------------------------------------------------------------------
  // comparison using the channel number (special comparison operator)
  template <typename Comp = std::less<raw::Channel_t>>
  struct ChannelComparison {
    
    using Comparer_t = Comp;
    
    Comparer_t comp;
    
    constexpr ChannelComparison(): comp() {}
    constexpr ChannelComparison(Comparer_t comp): comp(comp) {}
    
    template <typename A, typename B>
    constexpr bool operator() (A const& a, B const& b) const
      { return comp(channelOf(a), channelOf(b)); }
    
    static constexpr raw::Channel_t channelOf(raw::Channel_t channel)
      { return channel; }
    static raw::Channel_t channelOf(raw::OpDetWaveform const& waveform)
      { return waveform.ChannelNumber(); }
    static raw::Channel_t channelOf(sbn::OpDetWaveformMeta const& waveformMeta)
      { return waveformMeta.ChannelNumber(); }
    static raw::Channel_t channelOf
      (icarus::trigger::OpticalTriggerGateData_t const& gate)
      { return gate.channel(); }
    template <typename Gate, typename OpDetInfo>
    static raw::Channel_t channelOf
      (icarus::trigger::TrackedTriggerGate<Gate, OpDetInfo> const& gate)
      { return channelOf(gate.gate()); }
    
  }; // struct ChannelComparison
  
  
  //----------------------------------------------------------------------------
  /*
   * Currently util::Quantity objects have serious issues with FHiCL validation,
   * so we are not reading them directly as such. This pulls in a number of
   * workarounds for:
   *  * every `Quantity` parameter: its FHiCL parameter must be declared as
   *    the base type `Quantity::value_t`
   *  * every optional parameter must be then read indirectly since a reference
   *    to the exact type of the parameter is required in such reading
   *  * for sequences, direct vector assignment
   *    (`vector<Quantity> = vector<Quantity::value_t>`) won't work, and since
   *    `Sequence::operator()` returns a temporary, this needs to be wraped too
   */
  template <typename SeqValueType>
  struct FHiCLsequenceWrapper {
    SeqValueType seqValue;
    
    FHiCLsequenceWrapper(SeqValueType seqValue): seqValue(seqValue) {}
    
    template <typename T>
    operator std::vector<T>() const
      { return { seqValue.cbegin(), seqValue.cend() }; }
    
  }; // FHiCLsequenceWrapper
  
  
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
    = std::lower_bound(fGates.begin(), gend, channel, ::ChannelComparison<>());
  return (
      (iGate != gend)
      && (::ChannelComparison<>::channelOf(*iGate) == channel)
    )
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
  : fChannelThresholds(FHiCLsequenceWrapper(config.ChannelThresholds()))
{
  std::sort(fChannelThresholds.begin(), fChannelThresholds.end());
  
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
