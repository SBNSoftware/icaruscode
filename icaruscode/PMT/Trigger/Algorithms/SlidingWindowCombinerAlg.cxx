/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowCombinerAlg.cxx
 * @brief  Algorihtm to combine trigger channels into sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 6, 2020
 * @see    icaruscode/PMT/Algorithms/Trigger/SlidingWindowCombinerAlg.h
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowCombinerAlg.h"

// C/C++ standard libraries
#include <algorithm> // std::sort()
#include <utility> // std::move()


// -----------------------------------------------------------------------------
icarus::trigger::SlidingWindowCombinerAlg::SlidingWindowCombinerAlg(
  Windows_t const& windows,
  std::string logCategory /* = "SlidingWindowCombinerAlg" */
  )
  : fWindowChannels(sortedWindowChannels(windows))
  , fLogCategory(std::move(logCategory))
  {}


// -----------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowCombinerAlg::sortedWindowChannels
  (Windows_t const& windows) -> Windows_t
{
  Windows_t newWindows { windows };
  for (auto& window: newWindows) std::sort(window.begin(), window.end());
  return newWindows;
} // icarus::trigger::SlidingWindowCombinerAlg::sortWindowChannels()


// -----------------------------------------------------------------------------
