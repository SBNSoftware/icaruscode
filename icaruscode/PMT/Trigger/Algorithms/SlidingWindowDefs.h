/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.h
 * @brief  Definition for PMT sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 26, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFS_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFS_H


// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus::trigger {
  
  // --- BEGIN -- Optical detector windows -------------------------------------
  /**
   * @name Optical detector windows
   * 
   * An optical detector window is just a group of optical detectors.
   * For trigger, these windows comprise contiguous optical detectors, and may
   * overlap for better coverage.
   * 
   * The algorithm `SlidingWindowDefinitionAlg` allows the creation of "sliding"
   * windows. The information on a single window is encoded in a standard
   * container (`TriggerWindowChannels_t`), with a unique channel number
   * representing each detector in the window in no particular order.
   * 
   * Definitions ("aliases") are here provided for convenience, together with
   * a couple of functions to print the content of a window or a set of windows.
   */
  /// @{
  
  /// Type of optical detector channel list in a window.
  using TriggerWindowChannels_t = std::vector<raw::Channel_t>;

  /// Definition of all windows.
  using TriggerWindowDefs_t = std::vector<TriggerWindowChannels_t>;
  
  // --- BEGIN -- Optical detector window dumping on stream --------------------
  
  /// Prints the composition of the optical detector `window` inline.
  void printTriggerWindowChannels
    (std::ostream& out, TriggerWindowChannels_t const& window);
  
  /// Prints the composition of all `windows` in long format.
  void printTriggerWindowDefs
    (std::ostream& out, TriggerWindowDefs_t const& windows);
  
  
  // ---------------------------------------------------------------------------
  namespace details {
    
    struct DumpTriggerWindowChannelWrapper
      { TriggerWindowChannels_t const* window; };
    struct DumpTriggerWindowDefWrapper
      { TriggerWindowDefs_t const* windows; };
    
    std::ostream& operator<<
      (std::ostream& out, DumpTriggerWindowChannelWrapper window);
    
    std::ostream& operator<<
      (std::ostream& out, DumpTriggerWindowDefWrapper windows);
    
  } // namespace details
  // ---------------------------------------------------------------------------
  
  
  /**
   * Helper for printing a `TriggerWindowChannels_t` into a stream.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * TriggerWindowChannels_t window; // ... filled with some channels
   * std::cout << "Window: " << icarus::trigger::dumpTriggerWindowChannels(window);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  auto dumpTriggerWindowChannels(TriggerWindowChannels_t const& window)
    -> details::DumpTriggerWindowChannelWrapper;
  
  /**
   * Helper for printing a TriggerWindowDefs_t into a stream.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * TriggerWindowDefs_t windows; // ... filled with some content
   * std::cout << "Windows: " << icarus::trigger::dumpTriggerWindowDefs(windows);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  auto dumpTriggerWindowDefs(TriggerWindowDefs_t const& windows)
    -> details::DumpTriggerWindowDefWrapper;
  
  /// @}
  // --- END ---- Optical detector windows -------------------------------------
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
// --- inline implementation
// -----------------------------------------------------------------------------
inline std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, DumpTriggerWindowChannelWrapper window)
  { printTriggerWindowChannels(out, *(window.window)); return out; }


// -----------------------------------------------------------------------------
inline auto icarus::trigger::dumpTriggerWindowChannels
  (TriggerWindowChannels_t const& window)
  -> details::DumpTriggerWindowChannelWrapper
  { return { &window }; }


// -----------------------------------------------------------------------------
inline std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, DumpTriggerWindowDefWrapper windows)
  { printTriggerWindowDefs(out, *(windows.windows)); return out; }


// -----------------------------------------------------------------------------
inline auto icarus::trigger::dumpTriggerWindowDefs
  (TriggerWindowDefs_t const& windows) -> details::DumpTriggerWindowDefWrapper
  { return { &windows }; }


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFS_H
