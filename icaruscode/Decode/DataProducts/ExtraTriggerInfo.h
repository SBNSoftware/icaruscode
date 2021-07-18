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


#endif // SBNOBJ_COMMON_TRIGGER_EXTRATRIGGERINFO_H
