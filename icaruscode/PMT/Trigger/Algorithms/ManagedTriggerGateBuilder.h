/**
 * @file   icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h
 * @brief  Algorithm to produce trigger gates out of optical readout waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.tcc`
 *         `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.cxx`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icarusalg/Utilities/CommonChoiceSelectors.h" // util::SignalPolarity

// LArSoft/framework libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "fhiclcpp/types/TableFragment.h"


// framework libraries
#include "fhiclcpp/types/Atom.h"

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus::trigger { class ManagedTriggerGateBuilder; }
/**
 * @brief Base interface for gate builders.
 * @see `icarus::trigger::DynamicTriggerGateBuilder`,
 *      `icarus::trigger::FixedTriggerGateBuilder`
 * 
 * This base class provides a skeleton building algorithm that can be customised
 * by templates.
 * The allowed customization includes what to do when a threshold is crossed
 * in a gate.
 * 
 * The algorithm considers samples grouped in blocks of fixed (but configurable)
 * size. Each block is treated as a "macro-sample" and assigned a level based on
 * the value of its samples and a time based on a fixed, configurable offset
 * from its first sample (`BlockTimeReference`).
 * Currently, the level of the block is always taken as the _largest_ of the
 * levels of the enabled samples in the block (`SamplingPattern`), where the
 * level of the sample is how much above the baseline the sample is.
 * 
 * Note that actions are performed only when a block crosses a threshold.
 * The algorithm keeps track at each time of which are the thresholds enclosing
 * the signal level, and if the level crosses one of them, the gates associated
 * to those thresholds, and only them, are offered a chance to react.
 * 
 * 
 * 
 * ### Time alignment of the waveforms
 * 
 * The algorithm uses "blocks" of samples to determine the state of the gate.
 * The start of the first block is always aligned with the start of the waveform
 * (effectively, the algorithm uses the index of the sample in the waveform).
 * While this is arbitrary, it happens to be appropriate for the CAEN V1730
 * boards, which can "trigger" only in multiples of 16 nanoseconds from an
 * absolute reference, provided that the required alignment is a divisor
 * 16 ns.
 * 
 * Note that the timestamp of the waveforms may have an arbitrary shift,
 * which is common to all waveforms on the same channel in the entire event.
 * This timestamp is eventually used to turn the block tick into a time on the
 * output gate, so that the gates from different channels are correctly
 * synchronised (and can be used to simulate real time coincidence).
 * The blocks themselves may appear not to be synchronized between channels due
 * to the channel-dependent time offsets and corrections.
 * 
 * For simulation, the input waveforms need to preserve this discretization
 * (16 ns) in order for this mechanism to exactly reproduce the hardware.
 * If that's not the case, simulation will not yield the exact result that
 * hardware would, however the results should be statistically equivalent.
 * 
 * 
 * Specific configuration
 * -----------------------
 * 
 * All classes derived by this base algorithm should support the following
 * configuration parameters:
 *  * `SamplingPattern` (sequence of switches, default: `[ true ]`): the sample
 *    pattern used to detect a single gate change. The size of this sequence
 *    determines the size of the block of samples used to detect a change.
 *    The state is open if any of the samples in the block marked with `true`
 *    passes the threshold, and closed if all the samples in the block marked
 *    with `true` do not pass the threshold; the values of samples marked with
 *    `false` is always ignored. The first block starts with the sample `0` of
 *    the input, unless `SampleOffset` is specified, in which case it starts
 *    from it instead (see below).
 *    For example, a pattern of `[ true, true, true, true ]` will determine the
 *    state based on the values of four samples, declaring an open state if any
 *    of the four samples passes the threshold, and a closed state if no sample
 *    out of four passes that threshold. A pattern of
 *    `[ true, false, false, false ]` will instead only use one sample every
 *    four to determine the state of the gate for all four samples, ignoring the
 *    other three. The default, `[ true ]`, treats all the samples
 *    independently.
 *  * `BlockTimeReference` (integer, default: `0`): for each pattern, which
 *    is the sample providing its time; for example, if `BlockTimeReference`
 *    is `2` when an opening or closing is caused by a pattern, that change is
 *    recorded at the time of its third (#2) sample.
 *  * `SampleOffset` (non-negative integer; default: `0`): skips this many
 *    samples at the beginning of each waveform. This parameter is intended to
 *    provide an offset for the block alignment.
 * 
 * The configuration option `SamplePrescale` has been discontinued. It is
 * possible to achieve the same functionality of `SamplePrescale: N` by
 * specifying a pattern of `N` elements with only the first set, like in
 * `[ 1, 0, 0, ... ]`.
 */
class icarus::trigger::ManagedTriggerGateBuilder
  : public icarus::trigger::TriggerGateBuilder
{
  using Base_t = icarus::trigger::TriggerGateBuilder;
  
    public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::TableFragment<Base_t::Config> baseConfig;
    
    fhicl::Atom<util::SignalPolarity> Polarity {
      Name{ "Polarity" },
      Comment{ "waveform polarity" },
      util::SignalPolarity::Negative
      };
    
    fhicl::Sequence<bool> SamplingPattern {
      Name("SamplingPattern"),
      Comment("which samples to test within each block (default: one sample per block)"),
      std::vector{ true } // default: all, independently
      };
    
    fhicl::Atom<std::ptrdiff_t> BlockTimeReference {
      Name("BlockTimeReference"),
      Comment("the sample within the block which represents the block time"),
      0 // default: the first sample in the block
      };
      
    fhicl::Atom<std::size_t> SampleOffset {
      Name("SampleOffset"),
      Comment("skip this many samples from the beginning of each waveform"),
      0 // default: skip none
      };
    
  }; // Config
  
  
  /// Constructor: configures from FHiCL.
  ManagedTriggerGateBuilder(Config const& config);
  
  
    protected:
  
  // This class describes the interface of a gate manager but is incomplete.
  struct GateManager {
    
      protected:
    struct GateInfoBase {
      // OpticalTriggerGateData_t:
      using TriggerGate_t = Base_t::TriggerGates::triggergate_t;
      // icarus::trigger::ReadoutTriggerGate:
      using TriggerGateData_t = TriggerGate_t::TriggerGate_t;
      
      TriggerGate_t* pGate = nullptr; ///< Pointer to the gate.
      
      GateInfoBase(TriggerGate_t& gate): pGate(&gate) {}
      TriggerGateData_t* operator->() const { return &gate(); }
      TriggerGateData_t& gate() const { return pGate->gate(); }
      void addTrackingInfo(raw::OpDetWaveform const& waveform) const
        { return pGate->tracking().add(&waveform); }
      
      void belowThresholdAt(optical_tick tick); // undefined: don't call!
      void aboveThresholdAt(optical_tick tick); // undefined: don't call!
      
    }; // struct GateInfoBase
    
      public:
    using GateInfo_t = GateInfoBase;
    
    GateInfo_t create(GateInfo_t::TriggerGate_t& gate) const;
    
  }; // struct GateManager
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  
  util::SignalPolarity const fPolarity; ///< Polarity of the input waveforms.
  
  /// Offsets of the samples included in the pattern.
  std::vector<std::ptrdiff_t> fPatternIndices;
  
  std::size_t const fBlockSize; ///< How many samples are in each block.
  
  /// Which sample to use for the time of a block.
  std::ptrdiff_t const fBlockTimeReference;
  
  std::size_t const fSampleOffset; ///< Skip this many samples at the beginning.
  
  // --- END Configuration -----------------------------------------------------
  
  
  /// Returns a collection of `TriggerGates` objects sorted by threshold.
  template <typename GateMgr>
  std::vector<TriggerGates> unifiedBuild
    (GateMgr&& gateManager, std::vector<WaveformWithBaseline> const& waveforms)
    const;
  
  /// Returns the `buildChannelGates()` function with the appropriate polarity.
  template <typename GateInfo, typename GatesByChannel>
  auto selectChannelGateFunc(GatesByChannel const&) const;
  
  /// Computes the gates for all the waveforms in one optical channel.
  template <typename Ops, typename GateInfo, typename Waveforms>
  void buildChannelGates
    (std::vector<GateInfo>& channelGates, Waveforms const& channelWaveforms)
    const;
  
  /// Prints the class configuration.
  /// 
  /// Assumes the start of a new line and does not break the last one.
  virtual void doDumpConfiguration(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const override;
  
  /// Prints the algorithm-specific configuration.
  ///
  /// Assumes the start of a new line and does not break the last one.
  void dumpLocalConfiguration(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const;
  
  
    private:
  
  /// Returns the indices of non-null elements of the `pattern`.
  static std::vector<std::ptrdiff_t> patternToIndices
    (std::vector<bool> const& pattern);
  

}; // class icarus::trigger::ManagedTriggerGateBuilder


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
#include "icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.tcc"

//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_H
