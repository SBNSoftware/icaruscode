/**
 * @file   icaruscode/TPC/Utilities/TPCSIOVchannelStatusDBdumper.h
 * @brief  Defines the `lariov::TPCSIOVchannelStatusDBdumper` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanfird.edu)
 * @date   March 11, 2025
 * @see    icaruscode/TPC/Utilities/TPCSIOVchannelStatusDBdumper.cxx
 * 
 */

#ifndef ICARUSCODE_TPC_UTILITIES_TPCSIOVCHANNELSTATUSDBDUMPER_H
#define ICARUSCODE_TPC_UTILITIES_TPCSIOVCHANNELSTATUSDBDUMPER_H

// SBN/ICARUS libraries
#include "icaruscode/TPC/Utilities/TPCchannelStatusDBdumper.h"

// LArSoft libraries
#include "larevt/CalibrationDBI/Providers/SIOVChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/CalibrationDBIFwd.h" // lariov::DBTimeStamp_t

// framework libraries
#include "fhiclcpp/types/TableFragment.h"

// C/C++ standard libraries
#include <ostream>


// -----------------------------------------------------------------------------
namespace lariov { class TPCSIOVchannelStatusDBdumper; }
/**
 * @brief Dumps on screen the state of all channels at specified times.
 * 
 * Dumps the content of the channel status database interfaced with 
 * `SIOVChannelStatusProvider`, for a series of timestamps specified in the
 * configuration.
 * 
 * The dump is format is a header in the form:
 * ```
 * Status of XXXX channels for N timestamps:
 * ```
 * (`XXX` the total number of channels in the geometry, N the number of
 * timestamps in the configuration) and then, for each timestamp:
 * ```
 * === BEGIN TIMESTAMP: xxxxxxxxxxxxxxxxxxx =======================================
 * CH=0: status
 * CH=1: status
 * [...]
 * 
 * Counting B BAD channels: BadCH BadCH [...]
 * Counting N NOISY channels: NoisyCH NoisyCh [...]
 * Counting G good channels.
 * === END   TIMESTAMP: xxxxxxxxxxxxxxxxxxx =======================================
 * ```
 * 
 * 
 * Timestamp format
 * -----------------
 * 
 * The timestamp format is really defined by the database conventions.
 * In ICARUS database, for example, that is in nanoseconds and in UTC.
 * 
 * 
 * Configuration
 * --------------
 * 
 * 
 * * `PrintIndividualChannels` (flag, default: `false`): if set, one line will
 *    be printed for each channel, showing its status.
 * * `PrintSummary` (flag, default: `true`): if set, the list of noisy and of
 *    bad channels will be printed, together with a count of good channels.
 * 
 */
class lariov::TPCSIOVchannelStatusDBdumper {
  
    public:
  
  using DBTimeStamp_t = lariov::DBTimeStamp_t;
  
  struct Config {
    
    // using Name = fhicl::Name;
    // using Comment = fhicl::Comment;
    
    fhicl::TableFragment<lariov::TPCchannelStatusDBdumper::Config>
      TPCchannelStatusDBdumperConfig;
    
    // no further configuration so far
    
  }; // Config
  
  
  /// Configures and sets up the algorithm.
  TPCSIOVchannelStatusDBdumper(
    Config const& config,
    lariov::SIOVChannelStatusProvider& channelStatus, unsigned int nChannels
    );
  
  /// Returns the configured number of channels.
  unsigned int nChannels() const { return fDumper.nChannels(); }
  
  /**
   * @brief Dumps all the information on channels at the given `timestamp`.
   * @param out the output stream to write information into
   * @param timestamp the timestamp to snapshot the channel status at
   * 
   * @note While this method is `const` and does not change this object, it may
   *       (and usually does) change the state of the service provider that it
   *       was passed at construction.
   * 
   */
  void dumpTimestamp(std::ostream& out, DBTimeStamp_t timestamp) const;
  
  /// Dumps all the information on channels at the timestamp set in the provider.
  void dumpCurrentTimestamp(std::ostream& out) const;
  
  
  /**
   * @brief Enables dumps into `std::ostream`.
   * @param timestamp the timestamp to dump the channel status at
   * @return an opaque object that will trigger the dump
   * 
   * Example (with `dumper` an instance of `TPCSIOVchannelStatusDBdumper`):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << dumper.timestampToStream(1'650'000'000'000'000'000);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * It is equivalent to `dumpTimestamp()`:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * dumper.dumpTimestamp(std::cout, 1'650'000'000'000'000'000);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  auto timestampToStream(DBTimeStamp_t timestamp) const;
  
  
    private:
  
  // --- BEGIN ---  Configuration parameters  ----------------------------------
  
  // --- END -----  Configuration parameters  ----------------------------------
  
  // --- BEGIN ---  Cached service information  --------------------------------
  
  /// Channel status service provider.
  lariov::SIOVChannelStatusProvider& fChannelStatus;
  
  // --- END   ---  Cached service information  --------------------------------
  
  lariov::TPCchannelStatusDBdumper const fDumper; ///< Inner dumper algorithm.
  
  
  /// Helper structure for insertion in `std::ostream`.
  struct timestampDump {
    TPCSIOVchannelStatusDBdumper const* dumper = nullptr;
    DBTimeStamp_t timestamp = 0;
    
    void dump(std::ostream& out) const { dumper->dumpTimestamp(out, timestamp); }
  };
  
  friend std::ostream& operator<<
    (std::ostream&, TPCSIOVchannelStatusDBdumper::timestampDump);
  
}; // class lariov::TPCSIOVchannelStatusDBdumper


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace lariov {
  
  /// Helper function for insertion.
  /// @see `lariov::TPCSIOVchannelStatusDBdumper::timestampToStream()`
  std::ostream& operator<<
    (std::ostream&, TPCSIOVchannelStatusDBdumper::timestampDump);
  
} // namespace lariov


// -----------------------------------------------------------------------------
// ---  Inline definitions
// -----------------------------------------------------------------------------
inline auto lariov::TPCSIOVchannelStatusDBdumper::timestampToStream
  (DBTimeStamp_t timestamp) const
{
  return timestampDump{ this, timestamp };
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_TPC_UTILITIES_TPCSIOVCHANNELSTATUSDBDUMPER_H
