/**
 * @file   icaruscode/PMT/Trigger/TrailingDynamicTriggerGateBuilderTool_tool.h
 * @brief  Toolization of `icarus::trigger::TrailingDynamicTriggerGateBuilder`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algoritmhs/TrailingDynamicTriggerGateBuilder.h`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_TRAILINGDYNAMICTRIGGERGATEBUILDERTOOL_TOOL_H
#define ICARUSCODE_PMT_TRIGGER_TRAILINGDYNAMICTRIGGERGATEBUILDERTOOL_TOOL_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TrailingDynamicTriggerGateBuilder.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"


namespace icarus::trigger {
  
  //----------------------------------------------------------------------------
  /**
   * @brief Dynamic gate builder tool.
   * @see `icarus::trigger::TrailingDynamicTriggerGateBuilder`
   */
  struct TrailingDynamicTriggerGateBuilderTool
    : icarus::trigger::TrailingDynamicTriggerGateBuilder
  {
    
    using Parameters =
      art::ToolConfigTable<icarus::trigger::TrailingDynamicTriggerGateBuilder::Config>;
    
    /// Constructor: sets the configuration.
    TrailingDynamicTriggerGateBuilderTool(Parameters const& config)
      : TrailingDynamicTriggerGateBuilder(config()) {}
    
    
  }; // class icarus::trigger::TrailingDynamicTriggerGateBuilderTool
  
  
  //----------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_TRAILINGDYNAMICTRIGGERGATEBUILDERTOOL_TOOL_H
