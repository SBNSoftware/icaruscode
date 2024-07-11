/**
 * @file   icaruscode/PMT/Trigger/Algorithms/DynamicTriggerGateBuilder.cxx
 * @brief  Dynamic gate builder.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/DynamicTriggerGateBuilder.h`
 * 
 */

// class header
#include "icaruscode/PMT/Trigger/Algorithms/DynamicTriggerGateBuilder.h"

//------------------------------------------------------------------------------
void icarus::trigger::DynamicTriggerGateBuilder::doDumpConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  dumpLocalConfiguration(out, indent, firstIndent);
  Base_t::doDumpConfiguration(out, indent, "\n" + indent);
  
} // icarus::trigger::DynamicTriggerGateBuilder::doDumpConfiguration()


//------------------------------------------------------------------------------
void icarus::trigger::DynamicTriggerGateBuilder::dumpLocalConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  out << firstIndent << " * algorithm: DynamicTriggerGateBuilder"
    << "\n" << indent
      << " * open a gate when signal crosses a threshold,"
      << " close it when signal is back under threshold"
    ;
  
} // icarus::trigger::DynamicTriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------

