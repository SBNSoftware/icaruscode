/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowCombinerAlg.h
 * @brief  Algorithm to combine trigger channels into sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 6, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/SlidingWindowCombinerAlg.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWCOMBINERALG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWCOMBINERALG_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // TriggerGateIndex

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"


// C/C++ standard libraries
#include <set>
#include <vector>
#include <string>
#include <algorithm> // std::find(), std::binary_search()


// -----------------------------------------------------------------------------
// --- forward declarations
// ---
namespace icarus::trigger {
  template <typename GateObject> class TriggerGateIndex;
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
namespace icarus::trigger { class SlidingWindowCombinerAlg; }
class icarus::trigger::SlidingWindowCombinerAlg {

    public:
  /// Type of optical detector channel list in a window.
  using WindowChannels_t = std::vector<raw::Channel_t>;

  /// Type of content of all windows.
  using Windows_t = std::vector<WindowChannels_t>;

  /// Constructor: learns about the window pattern (keeps a reference).
  SlidingWindowCombinerAlg(
    Windows_t const& windows,
    std::string logCategory = "SlidingWindowCombinerAlg"
    );

  /// Combines the `gates` according to the configured grouping.
  template <typename GateObject>
  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> combine
    (std::vector<GateObject> const& gates) const;

    private:

  /// Content of channels of each window.
  Windows_t const fWindowChannels;

  std::string fLogCategory; ///< Category for messages to MessageFacility.

  /// Throws an exception if the gates are not suitable for input.
  template <typename GateObject>
  void checkInput
    (icarus::trigger::TriggerGateIndex<GateObject> const& gates) const;

  /// Returns the combination of the `channels` selected from the `gates`.
  template <typename GateObject>
  icarus::trigger::MultiChannelOpticalTriggerGate combineChannels(
    icarus::trigger::TriggerGateIndex<GateObject> const& gates,
    WindowChannels_t const& channels
    ) const;

  /// Adds the gate data of `input` to `dest`, unless it's already included.
  /// @return whether the addition happened
  template <typename GateObject>
  static bool mergeGateInto(
    icarus::trigger::MultiChannelOpticalTriggerGate& dest,
    GateObject const& input
    );

  /// Returns windows with numerically sorted channel numbers.
  static Windows_t sortedWindowChannels(Windows_t const& windows);

  /// Returns whether the container `c` has `value`.
  template <typename Cont, typename T>
  static bool inList(Cont const& c, T const& value)
    { auto cend = c.end(); return std::find(c.begin(), cend, value) != cend; }

}; // class icarus::trigger::SlidingWindowCombinerAlg



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename GateObject>
std::vector<icarus::trigger::MultiChannelOpticalTriggerGate>
icarus::trigger::SlidingWindowCombinerAlg::combine
  (std::vector<GateObject> const& gates) const
{

  icarus::trigger::TriggerGateIndex gateIndex(gates);

  checkInput(gateIndex);

  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> combinedGates;
  for (std::vector<raw::Channel_t> const& channels: fWindowChannels) {
    assert(!channels.empty());

    combinedGates.push_back(combineChannels(gateIndex, channels));

  } // for

  return combinedGates;
} // icarus::trigger::SlidingWindowCombinerAlg::combine()


// -----------------------------------------------------------------------------
template <typename GateObject>
void icarus::trigger::SlidingWindowCombinerAlg::checkInput
  (icarus::trigger::TriggerGateIndex<GateObject> const& gateIndex) const
{
  /*
   * Checks:
   *
   * * no gate is only partially contained in a window
   * * no channel in the gates is skipped
   *
   */

  auto const isChannelInWindow
    = [](raw::Channel_t channel, std::vector<raw::Channel_t> const& window)
      { return std::binary_search(window.begin(), window.end(), channel); };

  std::set<GateObject const*> usedGates;
  for (std::vector<raw::Channel_t> const& window: fWindowChannels) {

    for (raw::Channel_t const channel: window) {

      //
      // check that we have all the channels we need
      //
      GateObject const* gate = gateIndex.find(channel);
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

  std::set<GateObject const*> allGates;
  for (auto const channel: util::counter<raw::Channel_t >(gateIndex.nChannels()))
    allGates.insert(&(gateIndex[channel]));

  std::vector<GateObject const*> unusedGates;
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
    for (GateObject const* gate: unusedGates) {
      e << "{";
      for (raw::Channel_t const channel: gate->channels())
        e << " ch:" << std::to_string(channel);
      e << " }";
    } // for
    e << "\n";
    throw e;
  } // if unused

} // icarus::trigger::SlidingWindowCombinerAlg::checkInput()


// -----------------------------------------------------------------------------
template <typename GateObject>
icarus::trigger::MultiChannelOpticalTriggerGate
icarus::trigger::SlidingWindowCombinerAlg::combineChannels(
  icarus::trigger::TriggerGateIndex<GateObject> const& gates,
  WindowChannels_t const& channels
  ) const
{
  using TriggerGate_t = GateObject;

  assert(!channels.empty());

  auto iChannel = channels.begin();
  auto cend = channels.end();

  TriggerGate_t const& firstGate = gates[*iChannel];

  // add the first gate
  icarus::trigger::MultiChannelOpticalTriggerGate gate
    { static_cast<icarus::trigger::OpticalTriggerGate const&>(firstGate) };

  mf::LogTrace(fLogCategory) << "Input:  " << firstGate;
  while (++iChannel != cend) {

    TriggerGate_t const& inputGate = gates[*iChannel];

    mf::LogTrace(fLogCategory) << "Input:  " << inputGate;
    mergeGateInto(gate, inputGate);

  } // while
  mf::LogTrace(fLogCategory) << "Output: " << gate;

  return gate;

} // icarus::trigger::SlidingWindowCombinerAlg::combineChannel()


//------------------------------------------------------------------------------
template <typename GateObject>
bool icarus::trigger::SlidingWindowCombinerAlg::mergeGateInto(
  icarus::trigger::MultiChannelOpticalTriggerGate& dest,
  GateObject const& input
  )
{
  // loose check: assumes channel sets of all possible inputs do not overlap
  // and checks whether the destination has already any one of input channels
  auto const& channels = dest.channels();
  if (inList(channels, *(input.channels().begin()))) return false;

  dest.Sum(input);

  return true;

} // icarus::trigger::SlidingWindowCombinerAlg::mergeGateInto()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWCOMBINERALG_H
