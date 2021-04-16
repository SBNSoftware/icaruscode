/**
 * @file   icaruscode/Utilities/EventRegistry.h
 * @brief  Class keeping track of _art_ event IDs.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 22, 2021
 * @see    icaruscode/Utilities/EventRegistry.cxx
 * 
 * 
 * 
 * 
 */

#ifndef ICARUSCODE_UTILITIES_EVENTREGISTRY_H
#define ICARUSCODE_UTILITIES_EVENTREGISTRY_H


// framework libraries
#include "canvas/Persistency/Provenance/EventID.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <unordered_map>
#include <vector>
#include <mutex>
#include <utility> // std::pair<>
#include <string>
#include <functional> // std::hash<>
#include <optional>
#include <limits> // std::numeric_limits<>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace std {
  
  template <>
  struct hash<art::EventID> {
    
    std::size_t operator() (art::EventID const& ID) const noexcept
      {
        return std::hash<std::uint64_t>{}(
            (std::uint64_t{ ID.run() } << 40U)     // run:    24 bits (16M)
          + (std::uint64_t{ ID.subRun() } << 22U)  // subrun: 18 bits (256k)
          + (std::uint64_t{ ID.event() })          // event:  22 bits (4M)
          );
      }
    
  }; // hash<art::EventID>
  
} // namespace std


// -----------------------------------------------------------------------------
namespace sbn { class EventRegistry; }
/**
 * @brief Keeps a record of all registered events and their source.
 * 
 * This registry object will keep track of every time an event is "registered".
 * An event is represented by its ID (`art::EventID`), and it is associated to
 * the input file(s) it is stored in.
 * 
 * Input files can be registered separately by file name, or at the same time
 * with the event.
 * 
 * 
 * Thread safety
 * --------------
 * 
 * Registration of events can happen concurrently (thread-safe).
 * Registration of source files, instead, is not protected.
 * 
 * It is guaranteed that the source file records are never modified once
 * registered (i.e. neither removed from the registry, nor the path of a source
 * record changed).
 * 
 */
class sbn::EventRegistry {
  
    public:
  
  using EventID_t = art::EventID; ///< Type used to identify an event.
  
  using FileID_t = std::size_t; ///< Type used to identify a source file.
  
  /// Element of the registry for an event.
  struct EventRecord_t {
    
    std::vector<FileID_t> sourceFiles; ///< List of ID of source files.
    
  }; // EventRecord_t
  
  /// Type with event ID (`first`) and event record information (`second`).
  using EventIDandRecord_t = std::pair<EventID_t, EventRecord_t>;
  
  
  /// Mnemonic for no file ID.
  static constexpr FileID_t NoFileID = std::numeric_limits<FileID_t>::max();
  
  
  // -- BEGIN -- Source interface ----------------------------------------------
  /// @name Source interface
  /// @{
  
  /**
   * @brief Registers a source file and returns its ID in the registry.
   * @param fileName the name of the source file to be registered
   * @return the internal ID this registry will refer to `fileName` with
   * 
   * If `fileName` has already been registered, the existing ID is returned and
   * no other action is performed.
   */
  FileID_t recordSource(std::string const& fileName);
  
  /// Returns whether the specified file ID is registered as a source.
  bool hasSource(FileID_t const& fileID) const;
  
  /// Returns whether the specified file name is registered as a source.
  bool hasSource(std::string const& fileName) const;
  
  /// Returns the ID of the source with the specified file name (slow!).
  /// @return the ID of the source, or unset if `fileID` is not registered
  std::optional<FileID_t> sourceID(std::string const& fileName) const;
  
  /// Returns the name of the source associated to the specified file ID.
  /// @return the name of the source, or unset if `fileID` is not registered
  std::optional<std::string> sourceName(FileID_t const& fileID) const;
  
  /// Returns the name of the source associated to the specified file ID.
  /// @return the name of the source, or `defName` if `fileID` is not registered
  std::string sourceNameOr
    (FileID_t const& fileID, std::string const& defName) const;
  
  /// @}
  // -- END -- Source interface ------------------------------------------------
  
  
  // -- BEGIN -- Event interface -----------------------------------------------
  /// @name Event interface
  /// @{
  
  /// Returns a copy of all event records.
  std::vector<EventIDandRecord_t> records() const;
  
  /// Returns a copy of the specified event record.
  /// @return the record for `event`, or unset if `event` is not registered
  std::optional<EventRecord_t> eventRecord(art::EventID const& event) const;
  
  /**
   * @brief Registers an event and returns a copy of its record.
   * @param event ID of the event to be registered
   * @param sourceFileID ID of the source file to be associated with the `event`
   * @return the current copy, just updated, of the record of the `event`
   * @throw cet::exception (category: `"sbn::EventRegistry"`) if `sourceFileID`
   *        dose not match a registered source file
   * @see `recordEvent(EventID_t const&, std::string const&)`
   * 
   * The record of the specified `event` is updated adding `sourceFileID` among
   * its sources. A single `sourceFileID` can appear multiple times for the same
   * event, indicating a duplicate event in the same source file.
   * 
   * A source with `sourceFileID` must have been registered already
   * (`recordSource()`).
   * 
   */
  EventRecord_t recordEvent(EventID_t const& event, FileID_t sourceFileID);
  
  /// @}
  // -- END -- Event interface -------------------------------------------------
  
  
    private:
  
  /// Type for source file registry.
  using FileRegistry_t = std::vector<std::string>;
  
  /// Registered source file, by file ID key.
  std::vector<std::string> fSourceFiles;
  
  /// Registry of all events.
  std::unordered_map<EventID_t, EventRecord_t> fEventRegistry;
  
  mutable std::mutex fEventRegistryLock; ///< Lock for `fEventRegistry`.
  
  //@{
  /// Returns an iterator pointing to the specified file registry entry.
  FileRegistry_t::iterator findSource(std::string const& fileName);
  FileRegistry_t::const_iterator findSource(std::string const& fileName) const;
  //@}
  
  /// Copies all event records into `recordCopy`.
  void copyEventRecordsInto(std::vector<EventIDandRecord_t>& recordCopy) const;

  /// Returns a lock guard around `fEventRegistry`.
  std::lock_guard<std::mutex> lockEventRegistry() const;
  
  
  /// Converts an internal index in file source registry into a `FileID_t`.
  static FileID_t indexToFileID(std::size_t index);

  /// Converts a `FileID_t` into an internal index in file source registry.
  static std::size_t fileIDtoIndex(FileID_t fileID);

}; // sbn::EventRegistry


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_EVENTREGISTRY_H
