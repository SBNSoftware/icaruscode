/**
 * @file   icaruscode/PMT/Trigger/Algorithms/TrailingFixedTriggerGateBuilder.cxx
 * @brief  Fixed-length gate builder.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/TrailingFixedTriggerGateBuilder.h`
 * 
 */


// class header
#include "icaruscode/PMT/Trigger/Algorithms/TrailingFixedTriggerGateBuilder.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"


//------------------------------------------------------------------------------
//--- icarus::trigger::TrailingFixedTriggerGateBuilder
//------------------------------------------------------------------------------
icarus::trigger::TrailingFixedTriggerGateBuilder::TrailingFixedTriggerGateBuilder
  (Config const& config)
  : Base_t(config)
  , fGateDuration(config.GateDuration())
  , fExtendGate(config.ExtendGate())
{}
  
  
//------------------------------------------------------------------------------
void icarus::trigger::TrailingFixedTriggerGateBuilder::setup
  (detinfo::DetectorTimings const& timings)
{
  Base_t::setup(timings);
  
  fGateTicks = timings.toTicks<optical_time_ticks>(fGateDuration);
  
} // icarus::trigger::TrailingFixedTriggerGateBuilder::setup()


//------------------------------------------------------------------------------
void icarus::trigger::TrailingFixedTriggerGateBuilder::doDumpConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  dumpLocalConfiguration(out, indent, firstIndent);
  Base_t::doDumpConfiguration(out, indent, "\n" + indent);
  
} // icarus::trigger::TrailingFixedTriggerGateBuilder::doDumpConfiguration()


//------------------------------------------------------------------------------
void icarus::trigger::TrailingFixedTriggerGateBuilder::dumpLocalConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  out << firstIndent << " * algorithm: TrailingFixedTriggerGateBuilder"
    << "\n" << indent << " * open a " << fGateDuration
      << " gate when going under threshold"
    << "\n" << indent << " * if a crossing happens in while the gate is open, "
      << (fExtendGate? "extend the gate": "ignore it")
    ;
  
} // icarus::trigger::TrailingFixedTriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------
