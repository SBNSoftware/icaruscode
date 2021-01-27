/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.h
 * @brief  Algorithm composing PMT sliding windows from geometry information.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 26, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFINITIONALG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFINITIONALG_H

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.h"

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
namespace icarus::trigger { class SlidingWindowDefinitionAlg; }
/**
 * @brief Groups optical detector channels into windows based on position.
 *
 * This algorithm groups all the optical detectors in the specified detector
 * geometry description into "sliding windows".
 *
 * Each optical detector "wall" (detectors at the same drift coordinate, that is
 * on the same plane) is sliced in windows of a given size (`windowSize`, the
 * number of optical detector channels within) starting one after the other at
 * fixed intervals (`windowStride`).
 *
 * For example, a partition with window size 30 channels and stride also 30
 * channels on a wall of 90 optical detector channels will create 3 windows
 * with 30 channels each. A splitting with size 30 channels but stride only 15
 * channels will create 5 windows of 30 channels each, which overlap (like in
 * 0-29, 15-44, 30-59, 45-74 and 60-89).
 *
 * The windows are returned in the formats defined in
 * `icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.h`, that is a list
 * of windows, each being a list of optical detector channels.
 *
 *
 * Algorithm details
 * ------------------
 *
 * Optical detectors are split into "walls" and within each wall into "towers"
 * (a wall being a set of optical detectors at the same drift coordinate, i.e.
 * on a plane, and a tower being a set of detector in a wall which share the
 * horizontal position but are piled in the vertical one, _y_).
 * The task of separating the detectors in walls and towers is delegated to
 * the algorithm `icarus::trigger::PMTverticalSlicingAlg`.
 * 
 * Each wall is processed separately. The towers in the wall are sorted in the
 * non-vertical direction (in ICARUS they will have a pattern of 2 channels,
 * 3 channels, 3 channels, 2 channels, repeated 9 times), and the windows are
 * created starting from one end. Towers are progressively added to the window
 * until the window reaches the desided size (if it overshoots it, the window
 * size is not compatible with the geometry and an exception is thrown).
 * The process is repeated for each window, starting with the first tower,
 * then the tower starting after `windowStride` optical detectors, then the
 * tower starting after twice `windowStride` detectors, and so on, until a
 * a window is reached that can't be completed because we ran out of towers.
 * If there is no tower starting at the exact multiple of `windowStride`,
 * the stride parameter is not compatible with the detector geometry, and an
 * exception is thrown.
 *
 */
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
inline icarus::trigger::SlidingWindowDefinitionAlg::SlidingWindowDefinitionAlg(
  geo::GeometryCore const& geom,
  std::string logCategory /* = "SlidingWindowDefinitionAlg" */
  )
  : fLogCategory(std::move(logCategory))
  , fGeom(geom)
  {}


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWDEFINITIONALG_H
