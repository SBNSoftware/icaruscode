/**
 * @file   icaruscode/PMT/Trigger/Algorithms/DynamicTriggerGateBuilder.h
 * @brief  Dynamic gate builder.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/DynamicTriggerGateBuilder.cxx`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DYNAMICTRIGGERGATEBUILDER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DYNAMICTRIGGERGATEBUILDER_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/Utilities/quantities/electronics.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>


namespace icarus::trigger {
  
  // ---------------------------------------------------------------------------
  //
  // declarations
  //
  
  class DynamicTriggerGateBuilder;
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Dynamic gate builder.
 * @see `icarus::trigger::FixedTriggerGateBuilder`
 * 
 * Gates are opened when the waveform goes above threshold, and closed when it
 * goes below.
 * 
 * By construction a gate with higher threshold will be fully included in one
 * with lower threshold.
 * 
 */
class icarus::trigger::DynamicTriggerGateBuilder
  : public icarus::trigger::ManagedTriggerGateBuilder
{
  using Base_t = icarus::trigger::ManagedTriggerGateBuilder;
  
  class DynamicGateManager: private GateManager {
    
    struct DynamicGateInfo: public GateInfoBase {
      DynamicGateInfo(TriggerGate_t& gate): GateInfoBase(gate) {}
      
      void belowThresholdAt(optical_tick tick) { gate().closeAt(tick); }
      void aboveThresholdAt(optical_tick tick) { gate().openAt(tick); }
      
    }; // struct DynamicGateInfo
    
      public:
    using GateInfo_t = DynamicGateInfo;
    
    DynamicGateManager() = default;
    
    GateInfo_t create(GateInfo_t::TriggerGate_t& gate) const
      { return { gate }; }
    
  }; // struct DynamicGateManager
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config: public Base_t::Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
  }; // struct Config
  // --- END Configuration -----------------------------------------------------
  
  
  /// Constructor: sets the configuration.
  DynamicTriggerGateBuilder(Config const& config)
    : Base_t(config) {}
  
  
  /// Returns a collection of `TriggerGates` objects sorted by threshold.
  virtual std::vector<TriggerGates> build
    (std::vector<raw::OpDetWaveform> const& waveforms) const override
    { return unifiedBuild(DynamicGateManager(), waveforms); }
  
  
}; // class icarus::trigger::DynamicTriggerGateBuilder


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DYNAMICTRIGGERGATEBUILDER_H
