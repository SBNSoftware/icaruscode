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
  std::vector<raw::Channel_t> missingChannels /* = {} */,
  bool requireFullCoverage /* = true */,
  std::string logCategory /* = "SlidingWindowCombinerAlg" */
  )
  : fWindowChannels(sortedWindowChannels(windows))
  , fMissingChannels(sortChannels(std::move(missingChannels)))
  , fRequireFullCoverage(requireFullCoverage)
  , fLogCategory(std::move(logCategory))
  {}


// -----------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowCombinerAlg::firstChannelPresent
  (WindowChannels_t const& channels) const -> WindowChannels_t::const_iterator
{
  auto iChannel = channels.begin();
  auto const cend = channels.end();
  while (iChannel != cend) if (!isMissingChannel(*iChannel)) return iChannel;
  return cend;
} // icarus::trigger::SlidingWindowCombinerAlg::firstChannelPresent()


// -----------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowCombinerAlg::sortedWindowChannels
  (Windows_t const& windows) -> Windows_t
{
  Windows_t newWindows { windows };
  for (auto& window: newWindows) std::sort(window.begin(), window.end());
  return newWindows;
} // icarus::trigger::SlidingWindowCombinerAlg::sortWindowChannels()


// -----------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowCombinerAlg::sortChannels
  (std::vector<raw::Channel_t> channels) -> std::vector<raw::Channel_t>
{
  std::sort(channels.begin(), channels.end());
  return std::move(channels); // is move() needed?
} // icarus::trigger::SlidingWindowCombinerAlg::sortChannels


// -----------------------------------------------------------------------------
