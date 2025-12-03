/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderChannelID.h
 * @brief  Simple object representing an ID for an adder channel.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 6, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderChannelID.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELID_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELID_H

// LArSoft and framework libraries
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t
#include "cetlib_except/exception.h" // cet::exception

// C/C++ standard libraries
#include <iosfwd>
#include <string>

// -----------------------------------------------------------------------------
// forward declarations
namespace icarus::trigger {
  
  struct AdderChannelID;
  
  std::ostream& operator<< (std::ostream& out, AdderChannelID channel);
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
/**
 * @brief Adder channel ID type.
 * 
 * It's a glorified integral value which prints as hexadecimal.
 * Actually, it does that only for "special" channels (`0x1000` and above):
 * the PMT channel numbers are printed in base 10.
 * 
 */
class icarus::trigger::AdderChannelID {
  
  raw::Channel_t fChannel; ///< The actual channel value.
  
    public:
  
  /// Base value of the IDs assigned to adder channels.
  static constexpr raw::Channel_t BaseAdderChannelID = 0xA000;
  
  /// (Default) constructor: acquires an ID value (`0` by default).
  constexpr AdderChannelID(raw::Channel_t channel = 0): fChannel{ channel } {}
  
  /// Returns whether this object has a valid adder channel ID.
  constexpr bool isValid() const { return isValid(*this); }
  
  /// Returns the channel number directly.
  constexpr raw::Channel_t channel() const { return fChannel; }
  
  /// Converts to a `raw::Channel_t` type.
  constexpr explicit operator raw::Channel_t() const { return fChannel; } // FIXME remove explicit
  
  /// Converts to a string (same as free-function `operator<<`).
  explicit operator std::string() const;
  
  /// Direct access to the channel value.
  explicit operator raw::Channel_t&() { return fChannel; } // FIXME remove explicit
  
  /// @name Comparison operators
  /// @{
  
  constexpr bool operator== (icarus::trigger::AdderChannelID other) const
    { return fChannel == other.fChannel; }
  
  constexpr bool operator!= (icarus::trigger::AdderChannelID other) const
    { return fChannel != other.fChannel; }
  
  constexpr bool operator< (icarus::trigger::AdderChannelID other) const
    { return fChannel < other.fChannel; }
  
  constexpr bool operator> (icarus::trigger::AdderChannelID other) const
    { return fChannel > other.fChannel; }
  
  constexpr bool operator<= (icarus::trigger::AdderChannelID other) const
    { return fChannel <= other.fChannel; }
  
  constexpr bool operator>= (icarus::trigger::AdderChannelID other) const
    { return fChannel >= other.fChannel; }
  
  /// @}
  
  /// Returns whether the channel ID is a valid adder channel ID.
  /// It does not acknowledge the "legacy" adder channel numbering.
  static constexpr bool isValid(AdderChannelID const channel)
    { return (raw::Channel_t{ channel } & 0xF000) == BaseAdderChannelID; }
  
}; // class icarus::trigger::AdderChannelID


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELID_H
