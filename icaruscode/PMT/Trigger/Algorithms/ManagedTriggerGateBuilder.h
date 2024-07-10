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
 * Note that actions are performed only when the sample crosses a threshold.
 * The algorithm keeps track at each time of which are the thresholds enclosing
 * the signal level, and if the level crosses one of them, the gates associated
 * to those thresholds, and only them, are offered a chance to react.
 * 
 * 
 * Specific configuration
 * -----------------------
 * 
 * All classes derived by this base algorithm should support the following
 * configuration parameters:
 *  * `SamplePrescale` (positive integer; default: `1`): only consider one
 *    sample every `SamplePrescale` for discrimination. For example, if
 *    `SamplePrescale` is set to `4` (and `SampleOffset` is `0`), each waveform
 *    in input will be discriminated considering only samples #0, #4, #8 and so
 *    on (with a sampling rate of 2 ns this means the discrimination is
 *    performed only every 8 nanoseconds).
 *    This implies that if the waveform passes the threshold at sample #1 and by
 *    sample #4 is back to not passing the threshold, the crossing is not
 *    detected. The same holds if at sample #0 the threshold is already passed,
 *    at sample #1 it is not any more, but by sample #4 it passed again.
 *    If the parameter is set to `1` (default value), or `0` (special case),
 *    all samples are considered.
 *  * `SampleOffset` (non-negative integer; default: `0`): skips this many
 *    samples at the beginning of each waveform. This parameter is intended to
 *    provide an offset for the prescale: for example, in case the
 *    discrimination were desired for the last sample in a group of four instead
 *    of the first one, this could be achieved by setting `SamplePrescale` to
 *    `4` and `SampleOffset` to `3`.
 *    Note however that this parameter is honoured regardless of the prescale
 *    settings (i.e. even if prescale is `1`, and even if it is smaller than
 *    the specified offset).
 * 
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
    
    fhicl::Atom<std::size_t> SamplePrescale {
      Name("SamplePrescale"),
      Comment("only consider one sample out of this many (1 = consider all)"),
      1 // default: all
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
  
  std::size_t const fSamplePrescale; ///< Use only one out of this many samples.
  
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
  
}; // class icarus::trigger::ManagedTriggerGateBuilder


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
#include "icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.tcc"

//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_MANAGEDTRIGGERGATEBUILDER_H
