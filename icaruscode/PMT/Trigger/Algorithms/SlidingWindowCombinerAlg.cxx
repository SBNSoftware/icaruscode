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
auto icarus::trigger::SlidingWindowCombinerAlg::combine
  (std::vector<TrackedTriggerGate_t> const& gates) const
  -> std::vector<TrackedTriggerGate_t>
{

  TriggerGateIndex_t gateIndex { gates };

  checkInput(gateIndex);

  std::vector<TrackedTriggerGate_t> combinedGates;
  for (std::vector<raw::Channel_t> const& channels: fWindowChannels) {
    assert(!channels.empty());

    combinedGates.push_back(combineChannels(gateIndex, channels));

  } // for

  return combinedGates;
} // icarus::trigger::SlidingWindowCombinerAlg::combine()


// -----------------------------------------------------------------------------
void icarus::trigger::SlidingWindowCombinerAlg::checkInput
  (TriggerGateIndex_t const& gateIndex) const
{
  /*
   * Checks:
   *
   * * no gate is only partially contained in a window
   * * no channel in the gates is skipped (optional)
   *
   */

  auto const isChannelInWindow
    = [](raw::Channel_t channel, std::vector<raw::Channel_t> const& window)
      { return std::binary_search(window.begin(), window.end(), channel); };

  std::set<TrackedTriggerGate_t const*> usedGates;
  for (std::vector<raw::Channel_t> const& window: fWindowChannels) {

    for (raw::Channel_t const channel: window) {

      if (isMissingChannel(channel)) continue;

      //
      // check that we have all the channels we need
      //
      TrackedTriggerGate_t const* gate = gateIndex.find(channel);
      if (!gate) {
        cet::exception e("SlidingWindowCombinerAlg");
        e << "SlidingWindowCombinerAlg::checkInput(): No gate has channel "
          << channel << ", required in the window including channels: ";
        for (raw::Channel_t const channel: window) e << " " << channel;
        e << ".\n";
        throw e;
      }

      //
      // check that all channels of this gate are in this window
      //
      for (raw::Channel_t const gateChannel: gate->channels()) {
        if (isChannelInWindow(gateChannel, window)) continue;

        cet::exception e("SlidingWindowCombinerAlg");
        e << "SlidingWindowCombinerAlg::checkInput(): gate with channels {";
        for (raw::Channel_t const channel: gate->channels())
          e << " " << channel;
        e << " } is not fully included in window including channels: ";
        for (raw::Channel_t const channel: window) e << " " << channel;
        e << " (" << gateChannel << " is missing).\n";
        throw e;
      } // for gate channels

      usedGates.insert(gate);
    } // for channel in window

  } // for windows

  //
  // check for channels that should be missing but are here
  //
  std::vector<raw::Channel_t> spuriousChannels;
  for (raw::Channel_t const channel: fMissingChannels) {
    if (gateIndex.find(channel)) spuriousChannels.push_back(channel);
  } // for
  if (!spuriousChannels.empty()) {
    cet::exception e("SlidingWindowCombinerAlg");
    e << "SlidingWindowCombinerAlg::checkInput(): "
      << spuriousChannels.size() << " of the " << fMissingChannels.size()
      << " supposedly missing channels are actually included:";
    for (raw::Channel_t const channel: spuriousChannels) e << " " << channel;
    e << "\n";
    throw e;
  }
  
  //
  // check that all channels are here (except known missing ones)
  //
  if (fRequireFullCoverage) {
    std::set<TrackedTriggerGate_t const*> allGates;
    for (auto const channel: util::counter<raw::Channel_t >(gateIndex.nChannels())) {
      if (isMissingChannel(channel)) continue;
      allGates.insert(&(gateIndex[channel]));
    }

    std::vector<TrackedTriggerGate_t const*> unusedGates;
    std::set_difference(
      allGates.cbegin(), allGates.cend(), usedGates.cbegin(), usedGates.cend(),
      std::back_inserter(unusedGates)
      );
    if (!unusedGates.empty()) {
      // we don't have much information to identify "which" gates...
      cet::exception e("SlidingWindowCombinerAlg");
      e << "SlidingWindowCombinerAlg::checkInput(): "
        << "used only " << usedGates.size() << "/" << allGates.size()
        << " input trigger gates, " << unusedGates.size() << " were not used:";
      for (TrackedTriggerGate_t const* gate: unusedGates) {
        e << "{";
        for (raw::Channel_t const channel: gate->channels())
          e << " ch:" << std::to_string(channel);
        e << " }";
      } // for
      e << "\n";
      throw e;
    } // if unused
  } // if full coverage

} // icarus::trigger::SlidingWindowCombinerAlg::checkInput()


// -----------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowCombinerAlg::combineChannels
  (TriggerGateIndex_t const& gates, WindowChannels_t const& channels) const
  -> TrackedTriggerGate_t
{
  assert(!channels.empty());

  auto iChannel = firstChannelPresent(channels);
  auto const cend = channels.end();
  if (iChannel == cend) return {}; // empty gate, no channels inside

  TrackedTriggerGate_t const& firstGate = gates[*iChannel];

  // add the first gate
  TrackedTriggerGate_t gate { firstGate };

  mf::LogTrace(fLogCategory) << "Input:  " << firstGate;
  while (++iChannel != cend) {
    if (isMissingChannel(*iChannel)) continue;

    TrackedTriggerGate_t const& inputGate = gates[*iChannel];

    mf::LogTrace(fLogCategory) << "Input:  " << inputGate;
    mergeGateInto(gate, inputGate);

  } // while
  mf::LogTrace(fLogCategory) << "Output: " << gate;

  return gate;

} // icarus::trigger::SlidingWindowCombinerAlg::combineChannel()


// -----------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowCombinerAlg::firstChannelPresent
  (WindowChannels_t const& channels) const -> WindowChannels_t::const_iterator
{
  auto iChannel = channels.begin();
  auto const cend = channels.end();
  while (iChannel != cend) if (!isMissingChannel(*iChannel)) return iChannel;
  return cend;
} // icarus::trigger::SlidingWindowCombinerAlg::firstChannelPresent()


//------------------------------------------------------------------------------
bool icarus::trigger::SlidingWindowCombinerAlg::mergeGateInto
  (TrackedTriggerGate_t& dest, TrackedTriggerGate_t const& input)
{
  // loose check: assumes channel sets of all possible inputs do not overlap
  // and checks whether the destination has already any one of input channels
  auto const& channels = dest.gate().channels();
  if (inList(channels, *(input.gate().channels().begin()))) return false;

  dest = icarus::trigger::sumGates(std::move(dest), input);

  return true;

} // icarus::trigger::SlidingWindowCombinerAlg::mergeGateInto()


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
  return channels;
} // icarus::trigger::SlidingWindowCombinerAlg::sortChannels


// -----------------------------------------------------------------------------
