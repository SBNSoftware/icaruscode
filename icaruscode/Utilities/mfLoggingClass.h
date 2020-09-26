/**
 * @file   icaruscode/Utilities/mfLoggingClass.h
 * @brief  Base class facilitating logging to message facility.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 16, 2020
 */

#ifndef ICARUSCODE_UTILITIES_MFLOGGINGCLASS_H
#define ICARUSCODE_UTILITIES_MFLOGGINGCLASS_H

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>
#ifdef __cpp_lib_source_location
#  include <source_location>
#endif //  __cpp_lib_source_location


//------------------------------------------------------------------------------
namespace icarus::ns::util { class mfLoggingClass; }

/**
 * @brief Helper for logging classes.
 * 
 * A derived class can utilize the member functions of this class for easier
 * tracking of the boilerplate category:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * struct Algorithm: private icarus::mfLoggingClass {
 * 
 *   Algorithm(): icarus::ns::util::mfLoggingClass("Algorithm") {}
 * 
 *   double compute() const
 *     {
 *        mfLogInfo() << "Starting computation()";
 *        return 0.0;
 *     }
 *   
 * }; // class Algorithm
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */
class icarus::ns::util::mfLoggingClass {
  
  std::string fLogCategory; ///< Logging category string used for the messages.
  
    public:
  
  /// Constructor: initializes with the specified log category.
  mfLoggingClass(std::string const& logCategory): fLogCategory(logCategory) {}
  
  /// Returns the logging category string for this object.
  std::string logCategory() const { return fLogCategory; }
  
  /// Returns this object (as a logging class object).
  mfLoggingClass const& loggingClass() const { return *this; }
  
  // --- BEGIN -- Access to temporary loggers ----------------------------------
  /**
   * @name Access to temporary loggers
   * 
   * These methods return a temporary logger for fast logging:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * mfLogError() << "That was not a smart thing to do!";
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * The returned log can also be made a bit less temporary if some more complex
   * output is required:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * if (!reasons.empty()) {
   *   auto log = mfLogError();
   *   log << "That was not a smart thing to do, for "
   *     << size(reasons) << " reasons:";
   *   for (auto const& [ iReason, reason ]: util::enumerate(reasons))
   *     log << "\n " << (iReason+1) << ": " << reason;
   * } // if
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * The `file` and `lineNumber` argument are optional and passed directly to
   * the logger on construction. If specified, they provide information about
   * the location of the message source:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * mfLogError(__FILE__, __LINE__) << "That was not a smart thing to do!";
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * @note With C++20 it will be possible to extend the implementation to
   *       allow for automatic detection of the call location.
   */
  /// @{
  
  /// Returns a mf::LogError() stream for logging.
  mf::LogError mfLogError
    (std::string const& file = {}, int const lineNumber = 0) const
    { return { logCategory(), file, lineNumber }; }
  
  /// Returns a mf::LogWarning() stream for logging.
  mf::LogWarning mfLogWarning
    (std::string const& file = {}, int const lineNumber = 0) const
    { return { logCategory(), file, lineNumber }; }
  
  /// Returns a mf::LogProblem() stream for logging.
  mf::LogProblem mfLogProblem
    (std::string const& file = {}, int const lineNumber = 0) const
    { return { logCategory(), file, lineNumber }; }
  
  /// Returns a mf::LogInfo() stream for logging.
  mf::LogInfo mfLogInfo
    (std::string const& file = {}, int const lineNumber = 0) const
    { return { logCategory(), file, lineNumber }; }
  
  /// Returns a mf::LogVerbatim() stream for logging.
  mf::LogVerbatim mfLogVerbatim
    (std::string const& file = {}, int const lineNumber = 0) const
    { return { logCategory(), file, lineNumber }; }
  
  /// Returns a mf::LogDebug() stream for logging.
  mf::LogDebug mfLogDebug
    (std::string const& file = {}, int const lineNumber = 0) const
    { return { logCategory(), file, lineNumber }; }
  
  /// Returns a mf::LogTrace() stream for logging.
  mf::LogTrace mfLogTrace
    (std::string const& file = {}, int const lineNumber = 0) const
    { return { logCategory(), file, lineNumber }; }
  
#ifdef __cpp_lib_source_location
  
  /// Returns a mf::LogDebug() with information about the calling location.
  mf::LogDebug mfLogDebugLine
    (std::source_location const loc = std::source_location::current()) const
    { return { logCategory(), loc.file_name(), loc.line_number() }; }
  
#endif // __cpp_lib_source_location
  
  /// @}
  // --- END -- Access to temporary loggers ------------------------------------
  
}; // class icarus::ns::util::mfLoggingClass


#endif // ICARUSCODE_UTILITIES_MFLOGGINGCLASS_H
