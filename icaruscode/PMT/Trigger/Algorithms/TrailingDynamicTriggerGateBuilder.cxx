/**
 * @file   icaruscode/PMT/Trigger/Algorithms/TrailingDynamicTriggerGateBuilder.cxx
 * @brief  Dynamic gate builder.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/TrailingDynamicTriggerGateBuilder.h`
 * 
 */

// class header
#include "icaruscode/PMT/Trigger/Algorithms/TrailingDynamicTriggerGateBuilder.h"

//------------------------------------------------------------------------------
void icarus::trigger::TrailingDynamicTriggerGateBuilder::doDumpConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  dumpLocalConfiguration(out, indent, firstIndent);
  Base_t::doDumpConfiguration(out, indent, "\n" + indent);
  
} // icarus::trigger::TrailingDynamicTriggerGateBuilder::doDumpConfiguration()


//------------------------------------------------------------------------------
void icarus::trigger::TrailingDynamicTriggerGateBuilder::dumpLocalConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  out << firstIndent << " * algorithm: TrailingDynamicTriggerGateBuilder"
    << "\n" << indent
      << " * open a gate when signal goes under threshold,"
      << " close it when signal crosses the threshold"
    ;
  
} // icarus::trigger::TrailingDynamicTriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------

