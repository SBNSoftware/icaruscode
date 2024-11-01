/**
 * @file   icaruscode/PMT/Trigger/TrailingFixedTriggerGateBuilderTool_tool.h
 * @brief  Toolization of `icarus::trigger::TrailingFixedTriggerGateBuilder`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algoritmhs/TrailingFixedTriggerGateBuilder.h`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_TRAILINGFIXEDTRIGGERGATEBUILDERTOOL_TOOL_H
#define ICARUSCODE_PMT_TRIGGER_TRAILINGFIXEDTRIGGERGATEBUILDERTOOL_TOOL_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TrailingFixedTriggerGateBuilder.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"


namespace icarus::trigger {
  
  //----------------------------------------------------------------------------
  /**
   * @brief Fixed gate builder tool.
   * @see `icarus::trigger::TrailingFixedTriggerGateBuilder`
   */
  struct TrailingFixedTriggerGateBuilderTool
    : icarus::trigger::TrailingFixedTriggerGateBuilder
  {
    
    using Parameters =
      art::ToolConfigTable<icarus::trigger::TrailingFixedTriggerGateBuilder::Config>;
    
    /// Constructor: sets the configuration.
    TrailingFixedTriggerGateBuilderTool(Parameters const& config)
      : TrailingFixedTriggerGateBuilder(config()) {}
    
    
  }; // class icarus::trigger::TrailingFixedTriggerGateBuilderTool
  
  
  //----------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_TRAILINGFIXEDTRIGGERGATEBUILDERTOOL_TOOL_H
