/**
 * @file   icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h
 * @brief  Logical multi-level gate associated to one or more channels.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.cxx`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DATA_MULTICHANNELOPTICALTRIGGERGATE_H
#define ICARUSCODE_PMT_TRIGGER_DATA_MULTICHANNELOPTICALTRIGGERGATE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h" // raw::Channel_t TODO should be moved to OpDetWaveforms

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <iosfwd> // std::ostream


namespace icarus::trigger {
  
  class SingleChannelOpticalTriggerGate; // external
  
  // ---------------------------------------------------------------------------
  //
  // declarations
  //
  
  class MultiChannelOpticalTriggerGate;
  
  std::ostream& operator<<
    (std::ostream&, MultiChannelOpticalTriggerGate const&);

  MultiChannelOpticalTriggerGate sumTriggerGates
    (std::vector<SingleChannelOpticalTriggerGate> const& gates);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Logical multi-level gate associated to one or more channels.
 * 
 * This object is a trigger gate associated with one or more readout channels.
 * The object contains a single trigger gate, representing the combination of
 * all involved channels.
 */
class icarus::trigger::MultiChannelOpticalTriggerGate
  : public icarus::trigger::OpticalTriggerGate
{
  
    public:
  
  // --- BEGIN Query -----------------------------------------------------------
  /// @name Query
  /// @{
  
  /// Returns a list of channels contributing to this gate.
  std::vector<raw::Channel_t> channels() const;
  
  // --- END Query -------------------------------------------------------------
  
}; // class icarus::trigger::MultiChannelOpticalTriggerGate


#endif // ICARUSCODE_PMT_TRIGGER_DATA_MULTICHANNELOPTICALTRIGGERGATE_H
