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
