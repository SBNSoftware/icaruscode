/**
 * @file   icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h
 * @brief  Logical multi-level gate associated to a optical detector channel.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.cxx`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DATA_SINGLECHANNELOPTICALTRIGGERGATE_H
#define ICARUSCODE_PMT_TRIGGER_DATA_SINGLECHANNELOPTICALTRIGGERGATE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <iosfwd> // std::ostream


namespace icarus::trigger {
  
  // ---------------------------------------------------------------------------
  //
  // declarations
  //
  
  class SingleChannelOpticalTriggerGate;
  
  std::ostream& operator<<
    (std::ostream&, SingleChannelOpticalTriggerGate const&);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
// TODO move this into `lardataobj/RawData/OpDetWaveform.h`
namespace raw {
  
  /// Mnemonics for an invalid channel number.
  constexpr Channel_t InvalidChannel = std::numeric_limits<Channel_t>::max();
  
  /// Returns whether the specified one is a valid `channel` number.
  constexpr bool isValidChannel(Channel_t channel)
    { return channel != InvalidChannel; }
  
} // namespace raw


//------------------------------------------------------------------------------
/**
 * @brief Logical multi-level gate associated to a optical detector channel.
 * 
 * This object is a trigger gate associated with one readout channel.
 */
class icarus::trigger::SingleChannelOpticalTriggerGate
  : public icarus::trigger::OpticalTriggerGate
{
  //
  // NOTE while there is some infrastructure allowing for the presence of more
  //      than one waveform, this class is designed with a single channel in
  //      mind.
  //
  
    public:
  
  /// Constructor: a closed gate for the channel in `waveform`.
  SingleChannelOpticalTriggerGate(raw::OpDetWaveform const& waveform)
    : icarus::trigger::OpticalTriggerGate(waveform)
    {}
  
  /// Adds another waveform to the channel (unless it has just been added).
  bool add(raw::OpDetWaveform const& waveform);
  
  // --- BEGIN Query -----------------------------------------------------------
  /// @name Query
  /// @{
  
  /// Returns the channel this gate is on.
  raw::Channel_t channel() const
    { 
      return waveforms().empty()
        ? raw::InvalidChannel: refWaveform().ChannelNumber();
    }
  
  // --- END Query -------------------------------------------------------------
  
  
  
  /// Comparison operator: sorts by increasing channel number.
  bool operator< (SingleChannelOpticalTriggerGate const& other) const
    { return channel() < other.channel(); }
  
    private:
  
  friend std::ostream& operator<<
    (std::ostream&, SingleChannelOpticalTriggerGate const&);
  
  
  /// Returns the "reference" waveform, used when a single waveform is needed.
  bool hasRefWaveform() const { return !waveforms().empty(); }
  
  /// Returns the "reference" waveform, used when a single waveform is needed.
  raw::OpDetWaveform const& refWaveform() const
    { return *(waveforms().front()); }
  
  /// Returns the "reference" waveform, used when a single waveform is needed.
  raw::OpDetWaveform const* refWaveformPtr() const
    { return hasRefWaveform()? nullptr: &(refWaveform()); }
  

}; // class icarus::trigger::SingleChannelOpticalTriggerGate


#endif // ICARUSCODE_PMT_TRIGGER_DATA_SINGLECHANNELOPTICALTRIGGERGATE_H
