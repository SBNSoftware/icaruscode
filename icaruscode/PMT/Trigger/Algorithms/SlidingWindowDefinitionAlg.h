/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.h
 * @brief  Algorithm composing PMT sliding windows from geometry information.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 26, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFINITIONALG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFINITIONALG_H


// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t

// framework libraries
#include "cetlib_except/exception.h" // convenience

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>
#include <string>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
//
// forward declarations
// 
namespace geo { class GeometryCore; }

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
  
  
  class SlidingWindowDefinitionAlg;
  
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
class icarus::trigger::SlidingWindowDefinitionAlg {
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN -- Setup --------------------------------------------------------
  
  geo::GeometryCore const& fGeom; ///< Geometry service provider.
  
  // --- END ---- Setup --------------------------------------------------------
  
    public:
  
  /// Type of optical detector channel list in a window.
  using WindowChannels_t = icarus::trigger::TriggerWindowChannels_t;

  /// Definition of all windows.
  using WindowDefs_t = icarus::trigger::TriggerWindowDefs_t;
  
  
  /**
   * @brief Constructor: complete algorithm setup.
   * @param geom LArSoft geometry service provider
   * @param logCategory tag to use for messages to message facility service
   */
  SlidingWindowDefinitionAlg(
    geo::GeometryCore const& geom,
    std::string logCategory = "SlidingWindowDefinitionAlg"
    );
  
  
  // @{
  /**
   * @brief Performs the calculation and returns the sliding window composition.
   * @param windowSize number of optical detectors in each window
   * @param windowStride new window every `windowStride` optical detectors
   */
  WindowDefs_t makeWindows
    (unsigned int windowSize, unsigned int windowStride) const;
  
  WindowDefs_t operator()
    (unsigned int windowSize, unsigned int windowStride) const
    { return makeWindows(windowSize, windowStride); }
  
  //@}
  
  // @{
  /**
   * @brief Performs the calculation and returns the sliding window composition.
   * @param windowSize number of optical detectors in each window
   * 
   * Windows are contiguous and not overlapping.
   */
  WindowDefs_t makeWindows(unsigned int windowSize) const
    { return makeWindows(windowSize, windowSize); }
  
  WindowDefs_t operator() (unsigned int windowSize) const
    { return makeWindows(windowSize); }
  
  //@}
  
}; // icarus::trigger::SlidingWindowDefinitionAlg



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
// ---  icarus::trigger::SlidingWindowDefinitionAlg
// -----------------------------------------------------------------------------
inline icarus::trigger::SlidingWindowDefinitionAlg::SlidingWindowDefinitionAlg(
  geo::GeometryCore const& geom,
  std::string logCategory /* = "SlidingWindowDefinitionAlg" */
  )
  : fLogCategory(std::move(logCategory))
  , fGeom(geom)
  {}


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFINITIONALG_H
