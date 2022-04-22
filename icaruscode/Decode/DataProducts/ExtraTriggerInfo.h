/**
 * @file sbnobj/Common/Trigger/ExtraTriggerInfo.h
 * @brief Data product holding additional trigger information.
 * @authors Gianluca Petrillo (petrillo@slac.stanford.edu),
 *          Jacob Zettlemoyer (jzettle@fnal.gov>)
 * @date June 18, 2021
 * @see sbnobj/Common/Trigger/ExtraTriggerInfo.cxx
 */

#ifndef SBNOBJ_COMMON_TRIGGER_EXTRATRIGGERINFO_H
#define SBNOBJ_COMMON_TRIGGER_EXTRATRIGGERINFO_H


// SBN libraries
#include "icaruscode/Decode/BeamBits.h"

// C/C++ standard libraries
#include <array>
#include <iosfwd> // std::ostream
#include <limits> // std::numeric_limits<>
#include <cstdint> // std::uint64_t


// -----------------------------------------------------------------------------
namespace sbn { struct ExtraTriggerInfo; }
/**
 * @brief Additional information on trigger.
 * 
 * This data structure holds information from the trigger that has no place in
 * the standard LArSoft data products (`raw::Trigger` and
 * `raw::ExternalTrigger`).
 * 
 * Absolute times are counted since The Epoch and are in UTC time.
 * 
 * Some "historical" information is available:
 * *  count of gates and triggers from the start of the run;
 * *  count of gates and triggers from the previous trigger;
 * *  timestamps of previous gates and triggers.
 * 
 * Some historical information is available for different gate selections:
 * * including only triggers and gates the same source as this trigger source;
 * * including triggers and gates from any source (tagged as "any").
 * 
 */
struct sbn::ExtraTriggerInfo {
  
  /// Special ID value indicating the absence of an ID.
  static constexpr unsigned int NoID = std::numeric_limits<unsigned int>::max();
  
  /// Special timestamp value indicating the absence of timestamp.
  static constexpr std::uint64_t NoTimestamp
    = std::numeric_limits<std::uint64_t>::max();
  
  /// Special timestamp correction value indicating an unknown correction.
  static constexpr std::int64_t UnknownCorrection
    = std::numeric_limits<std::int64_t>::max();
  
  
  /// Type of this gate (`sbn::triggerSource::NBits` marks this object invalid).
  sbn::triggerSource sourceType { sbn::triggerSource::NBits };
  
  
  // --- BEGIN -- Since the beginning of the run -------------------------------
  /// @name Counters and times since the beginning of the run.
  /// @{
  
  /// Absolute timestamp of this trigger [ns]
  std::uint64_t triggerTimestamp { NoTimestamp };
  
  /// Absolute timestamp of the opening of this beam gate [ns]
  std::uint64_t beamGateTimestamp { NoTimestamp };
  
  /// The identifier of this trigger (usually matching the event number).
  unsigned int triggerID { NoID };
  
  /// Incremental counter of gates from any source opened from start of the run.
  unsigned int gateID { 0 };
  
  /// Incremental counter of triggers from this source from start of the run.
  unsigned int triggerCount { 0 };
  
  /// Incremental counter of gates from this source opened from start of the run.
  unsigned int gateCount { 0 };
  
  /// @}
  // --- END ---- Since the beginning of the run -------------------------------
  
  
  
  // --- BEGIN -- Since the previous trigger from this same source -------------
  /// @name Counters since the previous trigger from this same source
  /// @{
  
  /// Gates from this source since the previous trigger also of this source
  /// (minimum is `1`: the current gate).
  unsigned int gateCountFromPreviousTrigger { 0 };
  
  /// Triggers from any source since the previous trigger also of this source
  /// (minimum is `1`: the current trigger).
  unsigned int anyTriggerCountFromPreviousTrigger { 0 };
  
  /// Gates from any source since the previous trigger also of this source
  /// (minimum is `1`: the current trigger).
  unsigned int anyGateCountFromPreviousTrigger { 0 };
  
  /// @}
  // --- END ---- Since the previous trigger from this same source -------------
  
  
  // --- BEGIN -- Since the previous trigger from any source -------------------
  /// @name Counters since the previous trigger from any source
  /// @{
  
  /// Type of the previous trigger.
  sbn::triggerSource anyPreviousTriggerSourceType
    { sbn::triggerSource::Unknown };
  
  /// Gates from any source since the previous trigger also from any source
  /// (minimum is `1`: the current gate).
  unsigned int anyGateCountFromAnyPreviousTrigger { 0 };
  
  /// @}
  // --- END ---- Since the previous trigger from any source -------------------
  
  
  // --- BEGIN -- Additional timestamps ----------------------------------------
  /// @name Additional timestamps
  /// @{
  
  // we do not have timestamps of gates not associated to triggers
  
  /// Absolute timestamp of the previous trigger from this same source [ns]
  std::uint64_t previousTriggerTimestamp { NoTimestamp };
  
  /// Absolute timestamp of the previous trigger from any source [ns]
  std::uint64_t anyPreviousTriggerTimestamp { NoTimestamp };
  
  /// @}
  // --- END ---- Additional timestamps ----------------------------------------
  
  
  // --- BEGIN -- Decoding information -----------------------------------------
  /// @name Decoding information
  /// @{
  
  /// Correction added to the White Rabbit time to get the trigger time.
  std::int64_t WRtimeToTriggerTime { UnknownCorrection };
  
  /// @}
  // --- END ---- Decoding information -----------------------------------------

  // --- BEGIN -- Trigger topology ---------------------------------------------
  /**
   * @name Trigger topology
   * 
   * The information is represented in groups of cryostats and, within each of
   * them, of PMT "walls", that is groups of PMT lying on the same geometric
   * plane (in SBN detectors that is behind each anode).
   * 
   * Currently the information is represented by fixed size arrays, because the
   * overhead of a variable length container (`std::vector`) is comparable to
   * the data itself.
   */
  /// @{
  
  /// Maximum number of cryostats in the detector.
  static constexpr std::size_t MaxCryostats { 2 };
  
  /// Maximum number of PMT "walls" in a cryostat.
  static constexpr std::size_t MaxWalls { 2 };

  /// Mnemonic index for the east cryostat.
  static constexpr std::size_t EastCryostat { 0 };
  
  /// Mnemonic index for the west cryostat.
  static constexpr std::size_t WestCryostat { 1 };
  
  /// Mnemonic index for the east PMT wall within a cryostat.
  static constexpr std::size_t EastPMTwall { 0 };
  
  /// Mnemonic index for the west PMT wall within a cryostat.
  static constexpr std::size_t WestPMTwall { 1 };
  
  
  /// Trigger data pertaining a single cryostat.
  struct CryostatInfo {
    /// Count of triggers in this cryostat.
    unsigned long int triggerCount { 0 };
    
    /**
     * @brief Status of LVDS signals from PMT pairs at the time of the trigger.
     * 
     * There is one status per PMT wall (index mnemonic constants: `EastPMTwall`
     * and `WestPMTwall`).
     * 
     * 
     * ICARUS
     * -------
     * 
     * Bits are 48 per wall, 8 pairs from each on 6 PMT readout boards.
     * For the position of the PMT generating a LVDS channel, refer to the
     * configuration of the trigger emulation.
     * 
     * The 48 bits are distributed in the value as follow:
     * * the south part (lower _z_, lower channel number) is in the first 32-bit
     * * the north part (upper _z_, higher channel number) is in the last 32-bit
     * 
     * The 8 most significant bits of each 32-bit half-word are zeroed.
     * Then the most significant bits match the lowest channel numbers (lower
     * _z_). Splitting the detector in six regions (of 15 PMT each) so that it
     * looks, from south to north: `AA|BB|CC|DD|EE|FF`, the bit mask of the six
     * parts (each with 15 PMT, 8 LVDS channels, 8 bits) looks like:
     * `00AABBCC'00DDEEFF`.
     */
    std::array<std::uint64_t, MaxWalls> LVDSstatus { 0U, 0U };
    
    /// Returns whether there is some recorded LVDS activity.
    constexpr bool hasLVDS() const;
    
  }; // CryostatInfo
  
  
  /// Bits for the trigger location (@see `triggerLocation()`).
  unsigned int triggerLocationBits { 0U };
  
  /// Status of each LVDS channel in each PMT wall at trigger time.
  std::array<CryostatInfo, MaxCryostats> cryostats;
  
  /**
   * @brief Returns the location of the trigger.
   * 
   * The returned value is a mask of `sbn::bits::triggerLocation` bits.
   * To test the state of a bit, it needs to be converted into a mask, e.g.:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * bool onTPCEE
   *   = extraTriggerInfo.triggerLocation() & mask(sbn::triggerLocation::TPCEE);
   * bool onTPCxE = extraTriggerInfo.triggerLocation()
   *   & mask(sbn::triggerLocation::TPCEE, sbn::triggerLocation::TPCWE);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  sbn::bits::triggerLocationMask triggerLocation() const
    { return { triggerLocationBits }; }
  
  
  /// @}
  // --- END ---- Trigger topology ---------------------------------------------
  
  
  
  /// Returns whether this object contains any valid information.
  constexpr bool isValid() const noexcept
    { return sourceType != sbn::triggerSource::NBits; }
  
  
  /// Returns whether the specified `ID` is valid (i.e. is not `NoID`).
  static constexpr bool isValidID(unsigned int ID) noexcept
    { return ID != NoID; }
  
  /// Returns whether the timestamp `ts` is valid (i.e. is not `NoTimestamp`).
  static constexpr bool isValidTimestamp(std::uint64_t ts) noexcept
    { return ts != NoTimestamp; }
  
  /// Returns whether the `count` is valid (i.e. is not `0`).
  static constexpr bool isValidCount(unsigned int count) noexcept
    { return count != 0U; }
  
  
}; // sbn::ExtraTriggerInfo


namespace sbn {
  std::ostream& operator<< (std::ostream& out, ExtraTriggerInfo const& info);
}

// -----------------------------------------------------------------------------

constexpr bool sbn::ExtraTriggerInfo::CryostatInfo::hasLVDS() const {
  // C++20:
  //  return std::ranges::any_of(LVDSstatus, std::identity{});   
  for (std::uint64_t bits: LVDSstatus) if (bits) return true;
  return false;
} // sbn::ExtraTriggerInfo::CryostatInfo::empty()




#endif // SBNOBJ_COMMON_TRIGGER_EXTRATRIGGERINFO_H
