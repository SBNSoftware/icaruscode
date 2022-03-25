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
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // TriggerGateIndex
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h" // sumGates()
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t

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
  using WindowChannels_t = icarus::trigger::TriggerWindowChannels_t;

  /// Type of content of all windows.
  using Windows_t = icarus::trigger::TriggerWindowDefs_t;

  /// Trigger gate data type (tracks sources via `sbn::OpDetWaveformMeta`).
  using TrackedTriggerGate_t
    = icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>;
  
  /// Constructor: learns about the window pattern (keeps a reference).
  SlidingWindowCombinerAlg(
    Windows_t const& windows,
    std::vector<raw::Channel_t> missingChannels = {},
    bool requireFullCoverage = true,
    std::string logCategory = "SlidingWindowCombinerAlg"
    );

  /// Combines the `gates` according to the configured grouping.
  std::vector<TrackedTriggerGate_t> combine
    (std::vector<TrackedTriggerGate_t> const& gates) const;

  /// Returns if `channel` is configured to be missing.
  bool isMissingChannel(raw::Channel_t channel) const;


    private:
  
  using TriggerGateIndex_t
    = icarus::trigger::TriggerGateIndex<TrackedTriggerGate_t>;

  /// Content of channels of each window.
  Windows_t const fWindowChannels;
  
  /// Channels known (and required) to be missing (sorted).
  std::vector<raw::Channel_t> fMissingChannels;

  /// Whether to require all channels to be used.
  bool const fRequireFullCoverage;
  
  std::string fLogCategory; ///< Category for messages to MessageFacility.

  /// Throws an exception if the gates are not suitable for input.
  void checkInput(TriggerGateIndex_t const& gates) const;

  /// Returns the combination of the `channels` selected from the `gates`.
  TrackedTriggerGate_t combineChannels
    (TriggerGateIndex_t const& gates, WindowChannels_t const& channels) const;

  /// Returns an iterator to the first of the `channels` which is not missing.
  WindowChannels_t::const_iterator firstChannelPresent
    (WindowChannels_t const& channels) const;


  /// Adds the gate data of `input` to `dest`, unless it's already included.
  /// @return whether the addition happened
  static bool mergeGateInto
    (TrackedTriggerGate_t& dest, TrackedTriggerGate_t const& input);

  /// Returns windows with numerically sorted channel numbers.
  static Windows_t sortedWindowChannels(Windows_t const& windows);

  /// Returns a sorted copy of `channels`.
  static std::vector<raw::Channel_t> sortChannels
    (std::vector<raw::Channel_t> channels);

  /// Returns whether the container `c` has `value`.
  template <typename Cont, typename T>
  static bool inList(Cont const& c, T const& value)
    { auto cend = c.end(); return std::find(c.begin(), cend, value) != cend; }

}; // class icarus::trigger::SlidingWindowCombinerAlg



// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------
inline bool icarus::trigger::SlidingWindowCombinerAlg::isMissingChannel
  (raw::Channel_t channel) const
{
  return std::find(fMissingChannels.begin(), fMissingChannels.end(), channel)
    != fMissingChannels.end();
} // icarus::trigger::SlidingWindowCombinerAlg::isMissingChannel()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_SLIDINGWINDOWCOMBINERALG_H
