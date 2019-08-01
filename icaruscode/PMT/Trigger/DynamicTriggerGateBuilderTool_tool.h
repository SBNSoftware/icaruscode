/**
 * @file   icaruscode/PMT/Trigger/DynamicTriggerGateBuilderTool_tool.h
 * @brief  Toolization of `icarus::trigger::DynamicTriggerGateBuilder`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algoritmhs/DynamicTriggerGateBuilder.h`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DYNAMICTRIGGERGATEBUILDERTOOL_TOOL_H
#define ICARUSCODE_PMT_TRIGGER_DYNAMICTRIGGERGATEBUILDERTOOL_TOOL_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/DynamicTriggerGateBuilder.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"


namespace icarus::trigger {
  
  //----------------------------------------------------------------------------
  /**
   * @brief Dynamic gate builder tool.
   * @see `icarus::trigger::DynamicTriggerGateBuilder`
   */
  struct DynamicTriggerGateBuilderTool
    : icarus::trigger::DynamicTriggerGateBuilder
  {
    
    using Parameters =
      art::ToolConfigTable<icarus::trigger::DynamicTriggerGateBuilder::Config>;
    
    /// Constructor: sets the configuration.
    DynamicTriggerGateBuilderTool(Parameters const& config)
      : DynamicTriggerGateBuilder(config()) {}
    
    
  }; // class icarus::trigger::DynamicTriggerGateBuilderTool
  
  
  //----------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_DYNAMICTRIGGERGATEBUILDERTOOL_TOOL_H
