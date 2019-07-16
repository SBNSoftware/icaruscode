/**
 * @file   icaruscode/PMT/Trigger/Algorithms/FixedTriggerGateBuilder.cxx
 * @brief  Fixed-length gate builder.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/FixedTriggerGateBuilder.h`
 * 
 */


// class header
#include "icaruscode/PMT/Trigger/Algorithms/FixedTriggerGateBuilder.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"


//------------------------------------------------------------------------------
//--- icarus::trigger::FixedTriggerGateBuilder
//------------------------------------------------------------------------------
icarus::trigger::FixedTriggerGateBuilder::FixedTriggerGateBuilder
  (Config const& config)
  : Base_t(config)
  , fGateDuration(config.GateDuration())
{}
  
  
//------------------------------------------------------------------------------
void icarus::trigger::FixedTriggerGateBuilder::setup
  (detinfo::DetectorTimings const& timings)
{
  Base_t::setup(timings);
  
  fGateTicks = timings.toTicks<optical_time_ticks>(fGateDuration);
  
} // icarus::trigger::FixedTriggerGateBuilder::setup()


//------------------------------------------------------------------------------
