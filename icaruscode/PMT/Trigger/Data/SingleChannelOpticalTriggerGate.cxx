/**
 * @file   icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.cxx
 * @brief  Logical multi-level gate associated to a optical detector channel.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h`
 * 
 */

// class header
#include "icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h"

// framework libraries
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <ostream>


//------------------------------------------------------------------------------
//--- icarus::trigger::SingleChannelOpticalTriggerGate
//------------------------------------------------------------------------------
bool icarus::trigger::SingleChannelOpticalTriggerGate::add
  (raw::OpDetWaveform const& waveform)
{
  if (raw::isValidChannel(channel()) && (channel() != waveform.ChannelNumber()))
  {
    // currently we require each gate to be on the same channel
    throw cet::exception("TriggerGateBuilder")
      << "icarus::trigger::SingleChannelOpticalTriggerGate::add(): "
      << "can't add a waveform on channel " << waveform.ChannelNumber()
      << " to a gate on channel " << channel() << "\n";
  }
  return OpticalTriggerGate::add(waveform);
} // icarus::trigger::SingleChannelOpticalTriggerGate::add()


//------------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<< (
  std::ostream& out,
  icarus::trigger::SingleChannelOpticalTriggerGate const& gate
) {
  out << "on channel " << gate.channel() << " " <<
    static_cast<icarus::trigger::OpticalTriggerGate const&>(gate)
    ;
  return out;
} // icarus::trigger::operator<< (SingleChannelOpticalTriggerGate)


//------------------------------------------------------------------------------
