/**
 * @file   icaruscode/PMT/Trigger/Data/OpticalTriggerGate.cxx
 * @brief  A trigger gate data object for optical detector electronics.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see   icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h
 * 
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"

// C/C++ standard libraries
#include <ostream>
#include <algorithm> // std::inplace_merge(), std::unique(), ...


//------------------------------------------------------------------------------
namespace {
  /// Comparison for `raw::OpDetWaveform`.
  struct OpDetWaveformComp {
    
    /// Returns whether `a < b`.
    static bool less(raw::OpDetWaveform const& a, raw::OpDetWaveform const& b)
      {
        if (a.ChannelNumber() < b.ChannelNumber()) return true;
        if (a.ChannelNumber() > b.ChannelNumber()) return false;
        
        if (a.TimeStamp() < b.TimeStamp()) return true;
        if (a.TimeStamp() > b.TimeStamp()) return false;
        
        return false; // they're equivalent
      } // less()
    
    bool operator()
      (raw::OpDetWaveform const& a, raw::OpDetWaveform const& b) const
      { return less(a, b); }
    
    bool operator()
      (raw::OpDetWaveform const* a, raw::OpDetWaveform const* b) const
      { return less(*a, *b); }
    
  }; // OpDetWaveformComp
  
  
  //----------------------------------------------------------------------------
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::trigger::OpticalTriggerGate
//------------------------------------------------------------------------------
void icarus::trigger::OpticalTriggerGate::registerWaveforms
  (Waveforms_t const& moreWaveforms)
{
  //
  // add channels
  //
  associateChannels(extractChannels(moreWaveforms));

  //
  // merge the two lists, then sort them
  //
  auto const middle = fWaveforms.insert
    (fWaveforms.end(), moreWaveforms.begin(), moreWaveforms.end());
  std::inplace_merge
   (fWaveforms.begin(), middle, fWaveforms.end(), ::OpDetWaveformComp());

  // finally remove any duplicate (there should be none, should it?)
  auto const actualEnd
    = std::unique(fWaveforms.begin(), fWaveforms.end()); // pointer comparison
  fWaveforms.erase(actualEnd, fWaveforms.end());


} // OpticalTriggerGate::registerWaveforms()


//------------------------------------------------------------------------------
auto icarus::trigger::OpticalTriggerGate::mergeWaveforms
  (Waveforms_t const& a, Waveforms_t const& b) -> Waveforms_t
{
  Waveforms_t merged;
  merged.reserve(a.size() + b.size());
  std::merge(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(merged),
    ::OpDetWaveformComp()
    );
  
  // finally remove any duplicate (there should be non, should it?)
  auto const actualEnd
    = std::unique(merged.begin(), merged.end()); // pointer comparison
  merged.erase(actualEnd, merged.end());
  
  return merged;
} // icarus::trigger::OpticalTriggerGate::mergeWaveforms()


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate& icarus::trigger::OpticalTriggerGate::Min
  (OpticalTriggerGate const& other)
{
  GateData_t::Min(other);
  mergeWaveformsFromGate(other);
  return *this;
} // icarus::trigger::OpticalTriggerGate::Min()


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate& icarus::trigger::OpticalTriggerGate::Max
  (OpticalTriggerGate const& other)
{
  GateData_t::Max(other);
  mergeWaveformsFromGate(other);
  return *this;
} // icarus::trigger::OpticalTriggerGate::Max()


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate& icarus::trigger::OpticalTriggerGate::Sum
  (OpticalTriggerGate const& other)
{
  GateData_t::Sum(other);
  mergeWaveformsFromGate(other);
  return *this;
} // icarus::trigger::OpticalTriggerGate::Sum()


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate& icarus::trigger::OpticalTriggerGate::Mul
  (OpticalTriggerGate const& other)
{
  GateData_t::Mul(other);
  mergeWaveformsFromGate(other);
  return *this;
} // icarus::trigger::OpticalTriggerGate::Mul()


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate icarus::trigger::OpticalTriggerGate::Min
  (OpticalTriggerGate const& a, OpticalTriggerGate const& b)
  { auto combination { a }; combination.Min(b); return combination; }


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate icarus::trigger::OpticalTriggerGate::Max
  (OpticalTriggerGate const& a, OpticalTriggerGate const& b)
  { auto combination { a }; combination.Max(b); return combination; }


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate icarus::trigger::OpticalTriggerGate::Sum
  (OpticalTriggerGate const& a, OpticalTriggerGate const& b)
  { auto combination { a }; combination.Sum(b); return combination; }


//------------------------------------------------------------------------------
icarus::trigger::OpticalTriggerGate icarus::trigger::OpticalTriggerGate::Mul
  (OpticalTriggerGate const& a, OpticalTriggerGate const& b)
  { auto combination { a }; combination.Mul(b); return combination; }


//------------------------------------------------------------------------------
template <typename Op>
icarus::trigger::OpticalTriggerGate
icarus::trigger::OpticalTriggerGate::SymmetricCombination(
  Op&& op, OpticalTriggerGate const& a, OpticalTriggerGate const& b,
  TriggerGateTicks_t aDelay /* = TriggerGateTicks_t{ 0 } */,
  TriggerGateTicks_t bDelay /* = TriggerGateTicks_t{ 0 } */
  )
{
  return { 
    GateData_t::SymmetricCombination
      (std::forward<Op>(op), a, b, aDelay, bDelay),
    OpticalTriggerGate::mergeWaveforms(a.waveforms(), b.waveforms())
  };
  
} // icarus::trigger::OpticalTriggerGate::SymmetricCombination()


//------------------------------------------------------------------------------
bool icarus::trigger::OpticalTriggerGate::add
  (raw::OpDetWaveform const& waveform)
{
  // insertion keeps the list ordered and the elements unique
  auto const& insertionPoint = std::lower_bound
    (fWaveforms.begin(), fWaveforms.end(), &waveform, ::OpDetWaveformComp());
  if ((insertionPoint != fWaveforms.end()) && (*insertionPoint == &waveform))
    return false;
  fWaveforms.insert(insertionPoint, &waveform);
  addChannel(waveform.ChannelNumber());
  return true;
} // icarus::trigger::OpticalTriggerGate::add()


//------------------------------------------------------------------------------
auto icarus::trigger::OpticalTriggerGate::extractChannels
  (Waveforms_t const& waveforms) -> GateData_t::ChannelList_t
{
  GateData_t::ChannelList_t channels;
  channels.reserve(waveforms.size());
  std::transform(
    waveforms.begin(), waveforms.end(), std::back_inserter(channels),
    std::mem_fn(&raw::OpDetWaveform::ChannelNumber)
    );
  return channels;
} // icarus::trigger::OpticalTriggerGate::extractChannels()


//------------------------------------------------------------------------------
auto icarus::trigger::OpticalTriggerGate::waveformChannels
  (Waveforms_t const& waveforms) -> GateData_t::ChannelList_t
{
  return GateData_t::normalizeChannels(extractChannels(waveforms));
} // icarus::trigger::OpticalTriggerGate::waveformChannels()


//------------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, icarus::trigger::OpticalTriggerGate const& gate)
{
  out << gate.gateLevels();
  return out;
} // icarus::trigger::operator<< (icarus::trigger::OpticalTriggerGate)


//------------------------------------------------------------------------------

