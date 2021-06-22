/**
 * @file sbnobj/Common/Trigger/ExtraTriggerInfo.cxx
 * @brief Data product holding additional trigger information.
 * @authors Gianluca Petrillo (petrillo@slac.stanford.edu),
 *          Jacob Zettlemoyer (jzettle@fnal.gov>)
 * @date June 18, 2021
 * @see sbnobj/Common/Trigger/ExtraTriggerInfo.h
 */

// #include "sbnobj/Common/Trigger/ExtraTriggerInfo.h" // future location of:
#include "icaruscode/Decode/DataProducts/ExtraTriggerInfo.h"

// C/C++ standard library
#include <ostream>
#include <iomanip> // std::setfill(), std::setw()


// -----------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  struct TimestampDumper { std::uint64_t timestamp; };
  
  TimestampDumper dumpTimestamp(std::uint64_t timestamp)
    { return { timestamp }; }
  
  std::ostream& operator<< (std::ostream& out, TimestampDumper wrapper) {
    std::uint64_t const timestamp = wrapper.timestamp;
    if (sbn::ExtraTriggerInfo::isValidTimestamp(timestamp)) {
      out << (timestamp / 1'000'000'000) << "."
        << std::setfill('0') << std::setw(9) << (timestamp % 1'000'000'000)
        << " s";
    }
    else out << "<n/a>";
    return out;
  } // operator<< (TimestampDumper)
  
  
  // ---------------------------------------------------------------------------
  struct TriggerIDdumper { unsigned int ID; };
  
  TriggerIDdumper dumpTriggerID(unsigned int ID) { return { ID }; }
  
  std::ostream& operator<< (std::ostream& out, TriggerIDdumper wrapper) {
    unsigned int const ID = wrapper.ID;
    if (sbn::ExtraTriggerInfo::isValidID(ID)) out << ID;
    else                                      out << "<n/a>";
    return out;
  } // operator<< (TriggerIDdumper)
  
  
  // ---------------------------------------------------------------------------
  struct TriggerCountDumper { unsigned int count; };
  
  TriggerCountDumper dumpTriggerCount(unsigned int count) { return { count }; }
  
  std::ostream& operator<< (std::ostream& out, TriggerCountDumper wrapper) {
    unsigned int const count = wrapper.count;
    if (sbn::ExtraTriggerInfo::isValidCount(count)) out << count;
    else                                            out << "<n/a>";
    return out;
  } // operator<< (TriggerCountDumper)
  
  
  // ---------------------------------------------------------------------------
  long long int timestampDiff(std::uint64_t timestamp, std::uint64_t ref) {
    return static_cast<long long int>
      ((timestamp > ref)? (timestamp - ref): (ref - timestamp));
  }
  
  
  // ---------------------------------------------------------------------------
  
} // local namespace


// -----------------------------------------------------------------------------
std::ostream& sbn::operator<< (std::ostream& out, ExtraTriggerInfo const& info)
{
  if (!info.isValid()) {
    out << "<invalid>";
    return out;
  }
  
  // quite a load:
  out
    <<   "trigger ID=" << dumpTriggerID(info.triggerID) << " from source "
      << name(info.sourceType)
      << " at " << dumpTimestamp(info.triggerTimestamp)
      << " on beam gate ID=" << dumpTriggerID(info.gateID)
      << " at " << dumpTimestamp(info.beamGateTimestamp)
      << " (diff: "
      << timestampDiff(info.beamGateTimestamp, info.triggerTimestamp) << " ns)"
    << "\n"
      << "counts from this source: trigger="
        << dumpTriggerCount(info.triggerCount)
      << " beam=" << dumpTriggerCount(info.gateCount)
    << "\n"
      << "previous trigger from this source at "
      << dumpTimestamp(info.previousTriggerTimestamp)
      << ", triggers since: "
      << dumpTriggerCount(info.anyTriggerCountFromPreviousTrigger)
      << ", gates since: "
      << dumpTriggerCount(info.anyGateCountFromPreviousTrigger)
      << " ("
      << dumpTriggerCount(info.gateCountFromPreviousTrigger)
      << " from this same source)"
    << "\n"
      << "most recent trigger was from source "
      << name(info.anyPreviousTriggerSourceType) << " at "
      << dumpTimestamp(info.anyPreviousTriggerTimestamp)
      << " ("
      << dumpTimestamp(info.triggerTimestamp - info.anyPreviousTriggerTimestamp)
      << " earlier), and "
      << dumpTriggerCount(info.anyGateCountFromAnyPreviousTrigger)
      << " gates from any source have opened since"
    ;
  
  return out;
} // sbn::operator<< (ExtraTriggerInfo)


// -----------------------------------------------------------------------------
