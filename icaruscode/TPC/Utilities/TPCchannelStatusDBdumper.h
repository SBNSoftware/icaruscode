/**
 * @file   icaruscode/TPC/Utilities/TPCchannelStatusDBdumper.h
 * @brief  Defines the `lariov::TPCchannelStatusDBdumper` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanfird.edu)
 * @date   March 11, 2025
 * @see    icaruscode/TPC/Utilities/TPCchannelStatusDBdumper.cxx
 * 
 */

#ifndef ICARUSCODE_TPC_UTILITIES_TPCCHANNELSTATUSDBDUMPER_H
#define ICARUSCODE_TPC_UTILITIES_TPCCHANNELSTATUSDBDUMPER_H

// LArSoft libraries
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <ostream>


// -----------------------------------------------------------------------------
namespace lariov { class TPCchannelStatusDBdumper; }
/**
 * @brief Dumps on screen the state of all channels at specified times.
 * 
 * Dumps the content of the channel status database interfaced with 
 * `ChannelStatusProvider`.
 * 
 * This algorithm dumps the status at the time previously selected in the
 * service provider. It may be used with `ChannelStatusService` service managed
 * by _art_.
 * If control on which timestamp to dump is needed, and the service provider in
 * use is derived from `lariov::SIOVChannelStatusProvider`, then the algorithm
 * `lariov::TPCSIOVchannelStatusDBdumper` can be used instead.
 * 
 * The dump format is in the form:
 * ```
 * CH=0: status
 * CH=1: status
 * [...]
 * 
 * Counting B BAD channels: BadCH BadCH [...]
 * Counting N NOISY channels: NoisyCH NoisyCh [...]
 * Counting G good channels.
 * ```
 * where the first part is printed if `PrintIndividualChannels` is set, while
 * the second part is printed if `PrintSummary` is set.
 * 
 * 
 * Configuration
 * --------------
 * 
 * * `PrintIndividualChannels` (flag, default: `false`): if set, one line will
 *    be printed for each channel, showing its status.
 * * `PrintSummary` (flag, default: `true`): if set, the list of noisy and of
 *    bad channels will be printed, together with a count of good channels.
 * 
 */
class lariov::TPCchannelStatusDBdumper {
  
    public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<bool> PrintIndividualChannels {
      Name{ "PrintIndividualChannels" },
      Comment{ "print one status line per channel" },
      false
      };
    
    fhicl::Atom<bool> PrintSummary {
      Name{ "PrintSummary" },
      Comment{ "print a summary with the list of non-good channels" },
      true
      };
    
  }; // Config
  
  
  /// Configures and sets up the algorithm.
  TPCchannelStatusDBdumper(
    Config const& config,
    lariov::ChannelStatusProvider& channelStatus, unsigned int nChannels
    );
  
  /// Returns the configured number of channels.
  unsigned int nChannels() const { return fNChannels; }
  
  /// Dumps all the information on channels at the timestamp set in the provider.
  void dump(std::ostream& out) const;
  
  
  /**
   * @brief Enables dumps into `std::ostream`.
   * @return an opaque object that will trigger the dump
   * 
   * Example (with `dumper` an instance of `TPCchannelStatusDBdumper`):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << dumper.toStream();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * It is equivalent to `dumpTimestamp()`:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * dumper.dumpTimestamp(std::cout);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (but the latter will not work with other loggers like `mf::LogInfo`).
   */
  auto toStream() const;
  
  
    private:
  
  // --- BEGIN ---  Configuration parameters  ----------------------------------
  
  bool const fPrintChannels; ///< Print individual channels.
  
  bool const fPrintSummary; ///< Print non-good channel summary.
  
  // --- END -----  Configuration parameters  ----------------------------------
  
  // --- BEGIN ---  Cached service information  --------------------------------
  unsigned int const fNChannels; ///< Number of TPC channels in the detector.
  
  /// Channel status service provider.
  lariov::ChannelStatusProvider const& fChannelStatus;
  
  // --- END   ---  Cached service information  --------------------------------
  
  /// Helper structure for insertion in `std::ostream`.
  struct dumpStruct {
    TPCchannelStatusDBdumper const* dumper = nullptr;
    
    void dump(std::ostream& out) const { dumper->dump(out); }
  };
  
  friend std::ostream& operator<<
    (std::ostream&, TPCchannelStatusDBdumper::dumpStruct);
  
}; // class lariov::TPCchannelStatusDBdumper


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace lariov {
  
  /// Helper function for insertion.
  /// @see `lariov::TPCchannelStatusDBdumper::toStream()`
  std::ostream& operator<<
    (std::ostream&, TPCchannelStatusDBdumper::dumpStruct);
  
} // namespace lariov


// -----------------------------------------------------------------------------
// ---  Inline definitions
// -----------------------------------------------------------------------------
inline auto lariov::TPCchannelStatusDBdumper::toStream() const
{
  return dumpStruct{ this };
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_TPC_UTILITIES_TPCCHANNELSTATUSDBDUMPER_H
