/**
 * @file icaruscode/Analysis/trigger/Objects/PMTcoverageInfo.h
 * @brief Derivative information from `raw::OpDetWaveform` data.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date November 29, 2021
 * 
 */

#ifndef ICARUSCODE_ANALYSIS_TRIGGER_OBJECTS_PMTCOVERAGEINFO_H
#define ICARUSCODE_ANALYSIS_TRIGGER_OBJECTS_PMTCOVERAGEINFO_H


// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <limits>


// -----------------------------------------------------------------------------
namespace sbn { struct PMTcoverageInfo; }
/**
 * @brief Derivative information from `raw::OpDetWaveform` data.
 * 
 * This objects stores some of the information from `raw::OpDetWaveform`,
 * with the notable exception of the content of the waveform.
 * In particular, it reports the time range covered by one waveform, and it may
 * store whether the range includes the beam gate opening time.
 * 
 * Times are in the same scale as for `raw::OpDetWaveform`, that is the
 * @ref DetectorClocksElectronicsStartTime "electronics time scale".
 */
struct sbn::PMTcoverageInfo {
  
  using Bits_t = unsigned int; ///< Type for bits.
  
  /// Magic value denoting the absence of time information.
  static constexpr double NoTime = std::numeric_limits<double>::lowest();
  static constexpr raw::Channel_t NoChannel
    = std::numeric_limits<raw::Channel_t>::max();
  
  /// Namespace for bits in the flags.
  struct bits {
    
    /// Whether this time interval includes the nominal beam gate opening.
    static constexpr Bits_t WithBeamGate = (1U << 0);
    
    /// Whether this time interval includes the hardware trigger time.
    static constexpr Bits_t WithTrigger = (1U << 1);
    
  }; // struct bits
  
  raw::Channel_t channel = NoChannel; ///< ID of the PMT channel.
  
  /// Time of the first sample in the waveform [us]
  double startTime = NoTime;
  /// Time at the end of the last sample in the waveform [us]
  double endTime = NoTime;
  
  /// Bitmask of the flags with state set (`1`).
  Bits_t setFlags = 0;
  /// Bitmask of the flags which have been assigned a value.
  Bits_t definedFlags = 0;
  
  
  // --- BEGIN -- Short bit managing interface ---------------------------------
  
  /// @name Bit query
  /// @{
  
  /// Returns whether the time interval includes the beam opening time.
  bool withBeamGate() const { return allSet(bits::WithBeamGate); }
  /// Returns whether the time interval does not include the beam opening time.
  bool withoutBeamGate() const { return allUnset(bits::WithBeamGate); }
  
  /// Returns whether the time interval includes the trigger time.
  bool withTrigger() const { return allSet(bits::WithTrigger); }
  /// Returns whether the time interval does not include the trigger time.
  bool withoutTrigger() const { return allUnset(bits::WithTrigger); }
  
  
  /// Returns whether all specified bits are defined and sets.
  bool allDefined(Bits_t mask) const { return (definedFlags & mask) == mask; }
  
  /// Returns whether all specified bits are defined and sets.
  bool allSet(Bits_t mask) const
    { return allDefined(mask) && ((setFlags & mask) == mask); }
  
  /// Returns whether all specified bits are defined and not sets.
  bool allUnset(Bits_t mask) const
    { return allDefined(mask) && ((setFlags & mask) == 0); }
  
  /// @}
  
  /// @name Bit setting
  /// @{
  
  /// Set all the bits in the mask to the specified values.
  void setBits(Bits_t mask, Bits_t values)
    {
      definedFlags |= mask;
      setFlags |= values &  mask;
      setFlags &= values | ~mask;
    }
  
  /// Sets or unsets (according to `set`) all the bits in the mask.
  void set(Bits_t mask, bool doSet) { if (doSet) set(mask); else unset(mask); }
  
  /// Set all the bits in the mask.
  void set(Bits_t mask) { setBits(mask, mask); }
  
  /// Unset all the bits in the mask.
  void unset(Bits_t mask) { setBits(mask, ~mask); }
  
  /// Sets all bits in the mask as undefined.
  void forget(Bits_t mask) { definedFlags &= ~mask; setFlags &= ~mask; }
  
  /// @}
  
  // --- END ---- Short bit managing interface ---------------------------------
  
}; // sbn::PMTcoverageInfo


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_ANALYSIS_TRIGGER_OBJECTS_PMTCOVERAGEINFO_H
