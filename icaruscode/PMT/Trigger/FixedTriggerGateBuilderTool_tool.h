/**
 * @file   icaruscode/PMT/Trigger/FixedTriggerGateBuilderTool_tool.h
 * @brief  Toolization of `icarus::trigger::FixedTriggerGateBuilder`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algoritmhs/FixedTriggerGateBuilder.h`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_FIXEDTRIGGERGATEBUILDERTOOL_TOOL_H
#define ICARUSCODE_PMT_TRIGGER_FIXEDTRIGGERGATEBUILDERTOOL_TOOL_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/FixedTriggerGateBuilder.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"


namespace icarus::trigger {
  
  //----------------------------------------------------------------------------
  /**
   * @brief Fixed gate builder tool.
   * @see `icarus::trigger::FixedTriggerGateBuilder`
   */
  struct FixedTriggerGateBuilderTool
    : icarus::trigger::FixedTriggerGateBuilder
  {
    
    using Parameters =
      art::ToolConfigTable<icarus::trigger::FixedTriggerGateBuilder::Config>;
    
    /// Constructor: sets the configuration.
    FixedTriggerGateBuilderTool(Parameters const& config)
      : FixedTriggerGateBuilder(config()) {}
    
    
  }; // class icarus::trigger::FixedTriggerGateBuilderTool
  
  
  //----------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_FIXEDTRIGGERGATEBUILDERTOOL_TOOL_H
