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
  , fExtendGate(config.ExtendGate())
{}
  
  
//------------------------------------------------------------------------------
void icarus::trigger::FixedTriggerGateBuilder::setup
  (detinfo::DetectorTimings const& timings)
{
  Base_t::setup(timings);
  
  fGateTicks = timings.toTicks<optical_time_ticks>(fGateDuration);
  
} // icarus::trigger::FixedTriggerGateBuilder::setup()


//------------------------------------------------------------------------------
void icarus::trigger::FixedTriggerGateBuilder::doDumpConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  dumpLocalConfiguration(out, indent, firstIndent);
  Base_t::doDumpConfiguration(out, indent, "\n" + indent);
  
} // icarus::trigger::FixedTriggerGateBuilder::doDumpConfiguration()


//------------------------------------------------------------------------------
void icarus::trigger::FixedTriggerGateBuilder::dumpLocalConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  out << firstIndent << " * algorithm: FixedTriggerGateBuilder"
    << "\n" << indent << " * open a " << fGateDuration
      << " gate at threshold crossing"
    << "\n" << indent << " * if a crossing happens in while the gate is open, "
      << (fExtendGate? "extend the gate": "ignore it")
    ;
  
} // icarus::trigger::FixedTriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------
