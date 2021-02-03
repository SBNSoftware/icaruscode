/**
 * @file   sbncode/Utilities/DuplicateEventTracker_service.cc
 * @brief  Service keeping track of _art_ event IDs.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 22, 2021
 * @see sbncode/Utilities/DuplicateEventTracker.h
 */


// library header
#include "icaruscode/Utilities/DuplicateEventTracker.h"

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Provenance/EventID.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// C/C++ standard libraries
#include <set>
#include <cassert>


// -----------------------------------------------------------------------------
sbn::DuplicateEventTracker::DuplicateEventTracker
  (Parameters const& config, art::ActivityRegistry& reg)
  // configuration parameters
  : fWarningOnly(config().WarningOnly())
  , fSkipSummary(config().SkipSummary())
  , fExceptionAtEnd(config().ExceptionAtEnd())
  , fLogCategory(config().LogCategory())
{
  
  reg.sPostOpenFile.watch(this, &DuplicateEventTracker::postOpenFile);
  
  reg.sPostSourceEvent.watch(this, &DuplicateEventTracker::postEventReading);
  
  reg.sPostEndJob.watch(this, &DuplicateEventTracker::postEndJob);

} // sbn::DuplicateEventTracker::DuplicateEventTracker()


// -----------------------------------------------------------------------------
void sbn::DuplicateEventTracker::printSummary() const {
  
  mf::LogInfo log(fLogCategory);
  log << "Summary of duplicate events encountered:";
  
  std::set<sbn::EventRegistry::FileID_t> duplicateFileIDs;
  unsigned int nDuplicateEvents = 0U;
  auto eventRecords = fEventRegistry.records();
  std::sort(eventRecords.begin(), eventRecords.end(),
    [](auto const& A, auto const& B){ return A.first < B.first; });
  for (auto const& [ eventID, record ]: eventRecords) {
    if (record.sourceFiles.size() <= 1U) continue;
    ++nDuplicateEvents;
    
    log << "\n  " << eventID << " from " << record.sourceFiles.size()
      << " sources:";
    for (auto const fileID: record.sourceFiles) {
      log << " [" << fileID << "]";
      duplicateFileIDs.insert(fileID);
    }
  } // for
  
  if (!duplicateFileIDs.empty()) {
    log << "\nDuplicate events from " << duplicateFileIDs.size() << " files:";
    for (auto const fileID: duplicateFileIDs) {
      // a bit cheating multithreading: we assume that between the copy of the
      // event records and the source file list no source has been moved or
      // removed; with the current interface, this is assured
      auto const sourceName = fEventRegistry.sourceName(fileID);
      log << "\n  [" << fileID << "]";
      if (sourceName) log << " '" << *sourceName << "'";
      else log << " N/A"; // should not happen though
    } // for sources
  } // if duplicates
  
  if (nDuplicateEvents > 0U) {
    log << "\nCounted " << nDuplicateEvents << " duplicate events from "
      << duplicateFileIDs.size() << " source files.";
  }
  else log << "\nNo duplicate events found.";
  
} // sbn::DuplicateEventTracker::printSummary()


// -----------------------------------------------------------------------------
void sbn::DuplicateEventTracker::postOpenFile(std::string const& fileName) {
  
  // register the current input file; if it exists already,
  // we just get its ID again (but that's troublesome anyway)
  fCurrentInputFileID = fEventRegistry.recordSource(fileName);
  
} // sbn::DuplicateEventTracker::postOpenFile()


// -----------------------------------------------------------------------------
void sbn::DuplicateEventTracker::postEventReading
  (art::Event const& event, art::ScheduleContext)
{
  
  sbn::EventRegistry::EventRecord_t const eventInfo
    = fEventRegistry.recordEvent(event.id(), fCurrentInputFileID);
  assert(eventInfo.sourceFiles.size() > 0U);
  if (eventInfo.sourceFiles.size() == 1U) return; // all done, no new duplicates
  
  ++fNDuplicateEvents;
  
  if (fWarningOnly) {
    mf::LogWarning(fLogCategory) 
      << "WARNING: event " << event.id() << " encountered "
      << eventInfo.sourceFiles.size() << " times,"
      << "\n  now from '"
      << fEventRegistry.sourceNameOr(eventInfo.sourceFiles.back(), "<?>")
      << "',"
      << "\n  first time from '"
      << fEventRegistry.sourceNameOr(eventInfo.sourceFiles.front(), "<?>")
      << "'"
      ;
    return;
  }
  
  if (fExceptionAtEnd) return;
  
  assert(eventInfo.sourceFiles.size() == 2U);
  throw cet::exception(fLogCategory)
    << "Duplicate event " << event.id() << " encountered"
    << "\n  first in '"
    << fEventRegistry.sourceNameOr(eventInfo.sourceFiles.front(), "<?>") << "'"
    << "\n  and now in '"
    << fEventRegistry.sourceNameOr(eventInfo.sourceFiles.back(), "<?>") << "'\n"
    ;
  
} // sbn::DuplicateEventTracker::postEventReading()


// -----------------------------------------------------------------------------
void sbn::DuplicateEventTracker::postEndJob() {
  
  if (!fSkipSummary) printSummary();
  
  if (fExceptionAtEnd && (fNDuplicateEvents > 0U)) {
    throw cet::exception(fLogCategory)
      << "Found " << fNDuplicateEvents << " duplicate events in the job.\n";
  }
  
} // sbn::DuplicateEventTracker::postEndJob()


// -----------------------------------------------------------------------------
DEFINE_ART_SERVICE(sbn::DuplicateEventTracker)

// -----------------------------------------------------------------------------

