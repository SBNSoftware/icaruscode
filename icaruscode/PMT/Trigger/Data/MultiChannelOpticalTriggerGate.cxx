/**
 * @file   icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.cxx
 * @brief  Logical multi-level gate associated to one or more channels.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h`
 * 
 */


// class header
#include "icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h"

// C/C++ standard libraries
#include <ostream>


//------------------------------------------------------------------------------
//--- icarus::trigger::MultiChannelOpticalTriggerGate
//------------------------------------------------------------------------------
std::vector<raw::Channel_t>
icarus::trigger::MultiChannelOpticalTriggerGate::channels() const {
  
  std::vector<raw::Channel_t> channels;
  
  // we skip invalid channels
  for (auto const* waveform: waveforms()) {
    raw::Channel_t const channel = waveform->ChannelNumber();
    
    if (!raw::isValidChannel(channel)) continue;
    
    if (!channels.empty() && (channel == channels.back())) continue;
    
    channels.push_back(channel);
  } // for
  
  return channels;
  
} // icarus::trigger::MultiChannelOpticalTriggerGate::channels()


//------------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<< (
  std::ostream& out,
  icarus::trigger::MultiChannelOpticalTriggerGate const& gate
) {
  auto const& channels = gate.channels();
  if (channels.empty())
    out << "(no channels)";
  else {
    out << "on " << channels.size() << " channels {";
    for (auto channel: channels) out << " " << channel;
    out << " },";
  }
  out << " " << static_cast<icarus::trigger::OpticalTriggerGate const&>(gate);
  return out;
} // icarus::trigger::operator<< (MultiChannelOpticalTriggerGate)


//------------------------------------------------------------------------------
icarus::trigger::MultiChannelOpticalTriggerGate
icarus::trigger::sumTriggerGates
  (std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const& gates)
{
  icarus::trigger::MultiChannelOpticalTriggerGate sum;
  for (auto const& gate: gates) {
    sum.Sum(gate);
  } // for
  return sum;
} // icarus::trigger::sumTriggerGates()


//------------------------------------------------------------------------------
