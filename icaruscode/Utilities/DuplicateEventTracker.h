/**
 * @file   icaruscode/Utilities/DuplicateEventTracker.h
 * @brief  Service keeping track of _art_ event IDs.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 22, 2021
 * @see    icaruscode/Utilities/DuplicateEventTracker_service.cc
 * 
 * Usage of the basic service functionality does not require this header:
 * it suffices to configure the service for execution in _art_ configuration
 * file.
 */

#ifndef ICARUSCODE_UTILITIES_DUPLICATEEVENTTRACKER_H
#define ICARUSCODE_UTILITIES_DUPLICATEEVENTTRACKER_H


// library header
#include "icaruscode/Utilities/EventRegistry.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h" // courtesy

// C/C++ standard libraries
#include <string>
#include <atomic>


// -----------------------------------------------------------------------------
// forward declarations
namespace art {
  class ActivityRegistry;
  class ScheduleContext; // actually unused
  class Event;
} // namespace art

// -----------------------------------------------------------------------------
namespace sbn { class DuplicateEventTracker; }
/**
 * @brief Keeps track of duplicate events in the job.
 * 
 * This _art_ service maintains a record of all the processed events, and it
 * takes action when one is processed more than once.
 * The main purpose of the service is to detect the duplication as soon as
 * possible and in such case to interrupt the job with an error, so that the
 * input can be amended.
 * 
 * The default action is to immediately throw an exception.
 * Alternatively, a warning message can be emitted each time an event is
 * reprocessed, and a summary of all duplication can be emitted at the end of
 * the job.
 * 
 * The service also offers an interface, to receive the current list of events
 * that have been processed so far.
 * 
 * See `sbn::EventRegistry` for a description of the (large) memory usage of
 * this service.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * * **WarningOnly** (flag, default: `false`): this flag has effect on warning
 *     messages and on exception throwing:
 *     * warning messages: if this flag is set to `true`, the service will emit
 *       a warning on each duplicate event encountered, otherwise it will not
 *       emit any warning (summary will still describe the duplicates if not
 *       skipped);
 *     * exception in case of duplicate events: if set to `false`, an exception
 *       will be thrown at the first duplicate event, unless `ExceptionAtEnd` is
 *       set, in which case the exception will be thrown at the end; if the flag
 *       is set to `true`, an exception is thrown at the end if `ExceptionAtEnd`
 *       is set, otherwise no exception is thrown at all
 * * **SkipSummary** (flag, default: `false`): when set to `true`, it disables
 *     the printout of the duplicate events at the end of the job; unless
 *     `WarningOnly` is also set, this summary is going to be one-entry only
 * * **ExceptionAtEnd** (flag, default: `false`): if duplicate events are
 *     detected, an exception is thrown _at the end of the job_ instead of
 *     at the first duplicate event; this allows to assess all the duplicate
 *     events at once; if set, it is recommended that `SkipSummary` be disabled;
 *     note that with this option enabled, if duplicate events are detected an
 *     exception is _always_ thrown (at the end of the job), regardless the
 *     setting of `WarningOnly`
 *     end of the job) 
 * * **LogCategory** (text, default: `"DuplicateEventTracker"`): tag of the
 *     output category used by this service to the message facility
 * 
 * 
 * Multi-threading notes
 * ======================
 * 
 * This service is designed to work with the concurrent processing of events
 * _from the same file_. Reading multiple files is not supported (but _art_ 3
 * doesn't do that anyway).
 * 
 */
class sbn::DuplicateEventTracker {
  
    public:
  
  /// Service configuration.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
    fhicl::Atom<bool> WarningOnly {
      Name("WarningOnly"),
      Comment("just emit a warning on screen instead of throwing an exception"),
      false // default
      };
    
    fhicl::Atom<bool> SkipSummary {
      Name("SkipSummary"),
      Comment("do not print a summary of the duplicate events at end of job"),
      false // default
      };
    
    fhicl::Atom<bool> ExceptionAtEnd {
      Name("ExceptionAtEnd"),
      Comment
        ("in case of duplicate, wait until the end of job to throw an exception"),
      false // default
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment
        ("tag of output category used by this service to the message facility"),
      "DuplicateEventTracker" // default
      };
    
    
  }; // Config
  
  using Parameters = art::ServiceTable<Config>;
  
  /// Constructor: reads the configuration.
  DuplicateEventTracker(Parameters const& config, art::ActivityRegistry& reg);

  
  /// Prints a summary of the current duplicates.
  void printSummary() const;


    protected:
  
  // --- BEGIN -- Configuration variables --------------------------------------
      
  bool const fWarningOnly; ///< Only warning, no exception.
  bool const fSkipSummary; ///< Do not print the summary of duplicate events.
  bool const fExceptionAtEnd; ///< Throw exception as late as possible.
  
  std::string const fLogCategory;  ///< Message service category tag.
  
  // --- END -- Configuration variables ----------------------------------------
  
  /// ID in the event registry of the current input file.
  sbn::EventRegistry::FileID_t fCurrentInputFileID
    = sbn::EventRegistry::NoFileID;
  
  sbn::EventRegistry fEventRegistry; ///< Record of all events and their source.
  
  std::atomic<unsigned int> fNDuplicateEvents{ 0U }; ///< Duplicate event count.
  
  
  // --- BEGIN -- Callback functions -------------------------------------------
  /// @name Callback functions
  /// @{
  
  /// Records the event and throws an exception depending on the configuration.
  void postEventReading(art::Event const& event, art::ScheduleContext);
  
  /// Records the current file name.
  void postOpenFile(std::string const& fileName);
  
  /// Prints the summary and throws an exception depending on the configuration.
  void postEndJob();
  
  /// @}
  // --- END -- Callback functions ---------------------------------------------
  
  
}; // sbn::DuplicateEventTracker


// -----------------------------------------------------------------------------
DECLARE_ART_SERVICE(sbn::DuplicateEventTracker, SHARED)


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_DUPLICATEEVENTTRACKER_H
