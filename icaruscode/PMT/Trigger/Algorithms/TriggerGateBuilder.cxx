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
#include "icarusalg/Utilities/WaveformOperations.h"

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
    static raw::Channel_t channelOf
      (icarus::trigger::SingleChannelOpticalTriggerGate const& gate)
      { return gate.channel(); }
    
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
icarus::trigger::SingleChannelOpticalTriggerGate&
icarus::trigger::TriggerGateBuilder::TriggerGates::gateFor
  (raw::OpDetWaveform const& waveform)
{
  // keeping the gates sorted by channel
  auto iGate = std::lower_bound
    (fGates.begin(), fGates.end(), waveform, ::ChannelComparison<>());
  if (iGate != fGates.end()) { // found, it's there already
    MF_LOG_TRACE(details::TriggerGateDebugLog)
      << "Appending waveform to trigger gate (thr=" << threshold()
      << ") of channel " << waveform.ChannelNumber();
    iGate->add(waveform);
    return *iGate; 
  }
  // add the new gate in channel order
  MF_LOG_TRACE(details::TriggerGateDebugLog)
    << "Creating a new trigger gate (thr=" << threshold()
    << ") for channel " << waveform.ChannelNumber();
  iGate = fGates.emplace(iGate, waveform);
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
