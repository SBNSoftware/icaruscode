/**
 * @file   icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h
 * @brief  A simple alias for a most commonly used `TrackedTriggerGate` type.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 7, 2022
 * @see    icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h
 * 
 * This library is header-only.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_TRACKEDOPTICALTRIGGERGATE_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_TRACKEDOPTICALTRIGGERGATE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"


// -----------------------------------------------------------------------------

namespace icarus::trigger {
  
  /**
   * @brief `TrackedTriggerGate` wrapper to trigger gate type commonly stored.
   * @tparam TrackedType type of tracked objects
   * 
   * The data in the gate, of type `icarus::trigger::OpticalTriggerGateData_t`,
   * is accessed via `gate()` and it is owned (i.e. if its source is a data
   * product, this is a copy of that data and not a reference to it).
   * 
   * Tracking information stores pointers to existing `TrackedType` objects
   * (typically `raw::OpDetWaveform` or `sbn::OpDetWaveformMeta`).
   * 
   */
  template <typename TrackedType>
  using TrackedOpticalTriggerGate = icarus::trigger::TrackedTriggerGate
    <icarus::trigger::OpticalTriggerGateData_t, TrackedType const*>;
  
    
} // namespace icarus::trigger

// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRACKEDOPTICALTRIGGERGATE_H
