/**
 * @file   icaruscode/Utilities/EventRegistry.cxx
 * @brief  Class keeping track of _art_ event IDs (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 22, 2021
 * @see    icaruscode/Utilities/EventRegistry.h
 * 
 * 
 * 
 * 
 */

// library header
#include "icaruscode/Utilities/EventRegistry.h"


// C/C++ standard libraries
#include <algorithm> // std::copy()
#include <iterator> // std::distance()


// -----------------------------------------------------------------------------
auto sbn::EventRegistry::recordSource(std::string const& fileName) -> FileID_t {
  
  auto const fileID = sourceID(fileName);
  if (fileID) return fileID.value();
  
  fSourceFiles.push_back(fileName);
  return indexToFileID(fSourceFiles.size() - 1U);
  
} // sbn::EventRegistry::recordSource()


// -----------------------------------------------------------------------------
bool sbn::EventRegistry::hasSource(FileID_t const& fileID) const
  { return fileIDtoIndex(fileID) < fSourceFiles.size(); }


// -----------------------------------------------------------------------------
bool sbn::EventRegistry::hasSource(std::string const& fileName) const
  { return findSource(fileName) != fSourceFiles.end(); }


// -----------------------------------------------------------------------------
auto sbn::EventRegistry::sourceID(std::string const& fileName) const
  -> std::optional<FileID_t>
{
  auto const iSource = findSource(fileName);
  return (iSource != fSourceFiles.end())
    ? std::optional{indexToFileID(std::distance(fSourceFiles.begin(), iSource))}
    : std::nullopt
    ;
} // sbn::EventRegistry::sourceID()


// -----------------------------------------------------------------------------
std::optional<std::string> sbn::EventRegistry::sourceName
  (FileID_t const& fileID) const
{
  return hasSource(fileID)
    ? std::optional{ fSourceFiles[fileIDtoIndex(fileID)] }: std::nullopt;
} // sbn::EventRegistry::sourceName()


// -----------------------------------------------------------------------------
std::string sbn::EventRegistry::sourceNameOr
  (FileID_t const& fileID, std::string const& defName) const
  { return sourceName(fileID).value_or(defName); }


// -----------------------------------------------------------------------------
auto sbn::EventRegistry::records() const -> std::vector<EventIDandRecord_t> {
  
  std::vector<EventIDandRecord_t> recordCopy;
  copyEventRecordsInto(recordCopy);
  return recordCopy;
  
} // sbn::EventRegistry::records()


// -----------------------------------------------------------------------------
auto sbn::EventRegistry::eventRecord(art::EventID const& event) const
  -> std::optional<EventRecord_t>
{
  
  auto const lg = lockEventRegistry();
  // BEGIN needs lock
  auto const iRecord = fEventRegistry.find(event);
  return (iRecord != fEventRegistry.end())
    ? std::optional{ iRecord->second }: std::nullopt;
  // END needs lock
  
} // sbn::EventRegistry::eventRecord()


// -----------------------------------------------------------------------------
auto sbn::EventRegistry::recordEvent
  (EventID_t const& event, FileID_t sourceFileID) -> EventRecord_t
{
  auto const lg = lockEventRegistry();
  // BEGIN needs lock
  auto& record = fEventRegistry[event];
  record.sourceFiles.push_back(sourceFileID);
  return record;
  // END needs lock
  
} // sbn::EventRegistry::recordEvent()


// -----------------------------------------------------------------------------
auto sbn::EventRegistry::findSource(std::string const& fileName)
  -> FileRegistry_t::iterator
  { return std::find(fSourceFiles.begin(), fSourceFiles.end(), fileName); }

auto sbn::EventRegistry::findSource(std::string const& fileName) const
  -> FileRegistry_t::const_iterator
  { return std::find(fSourceFiles.begin(), fSourceFiles.end(), fileName); }


// -----------------------------------------------------------------------------
void sbn::EventRegistry::copyEventRecordsInto
  (std::vector<EventIDandRecord_t>& recordCopy) const
{
  recordCopy.reserve(fEventRegistry.size()); // bet size does not change
  auto const lg = lockEventRegistry();
  // BEGIN needs lock
  recordCopy.reserve(fEventRegistry.size()); // has size() changed already?
  std::copy(fEventRegistry.cbegin(), fEventRegistry.cend(),
    std::back_inserter(recordCopy));
  // END needs lock
} // sbn::EventRegistry::records()


// -----------------------------------------------------------------------------
std::lock_guard<std::mutex> sbn::EventRegistry::lockEventRegistry() const
  { return std::lock_guard{ fEventRegistryLock }; }


// -----------------------------------------------------------------------------
auto sbn::EventRegistry::indexToFileID(std::size_t index) -> FileID_t
  { return static_cast<FileID_t>(index); }


// -----------------------------------------------------------------------------
std::size_t sbn::EventRegistry::fileIDtoIndex(FileID_t fileID)
  { return static_cast<std::size_t>(fileID); }


// -----------------------------------------------------------------------------

