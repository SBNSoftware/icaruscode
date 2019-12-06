/**
 * @file   icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h
 * @brief  Algorithm to produce trigger gates out of optical readout waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.cxx`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERGATEBUILDER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERGATEBUILDER_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // icarus::trigger::ADCCounts_t
#include "icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::count_if()
#include <functional> // std::mem_fn()
#include <optional>
#include <cstddef> // std::size_t


namespace icarus::trigger {
  
  // ---------------------------------------------------------------------------
  //
  // declarations
  //
  
  class TriggerGateBuilder;
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Algorithm to produce trigger gates out of optical readout waveforms.
 * @see `icarus::trigger::DynamicTriggerGateBuilder`,
 *      `icarus::trigger::FixedTriggerGateBuilder`
 *
 * This is an abstract class.
 * Derived algorithms need to provide a way to actually `build()` the gates.
 */
class icarus::trigger::TriggerGateBuilder {
  
    public:
  
  // --- BEGIN Data types ------------------------------------------------------
  /// Container of logical gates for all triggering channels for a threshold.
  class TriggerGates {
    
      public:
    using triggergate_t = icarus::trigger::SingleChannelOpticalTriggerGate;
    using GateData_t= std::vector<triggergate_t>;
    
    /// Constructor: acquires data of the specified gates.
    TriggerGates(ADCCounts_t threshold, GateData_t&& gates)
      : fThreshold(threshold), fGates(std::move(gates))
      {}
    
    /// Constructor: no trigger gate added.
    TriggerGates(ADCCounts_t threshold)
      : TriggerGates(threshold, {})
      {}
    
    /// Returns the threshold used for the gates.
    ADCCounts_t threshold() const { return fThreshold; }
    
    /// Returns the complete collection of trigger gates on all channels.
    GateData_t const& gates() const& { return fGates; }
    
    /**
     * @brief Yields the complete collection of trigger gates on all channels.
     * @return the gates that were contained in this object
     * 
     * After this call, the gate information is lost (except for the threshold).
     */
    GateData_t gates() && { return std::move(fGates); }
    
    /// Returns (and creates, if necessary) the gate for the specified waveform.
    icarus::trigger::SingleChannelOpticalTriggerGate& gateFor
      (raw::OpDetWaveform const& waveform);
    
    /// Dumps the content of this set of gates into the `out` stream.
    template <typename Stream>
    void dump(Stream& out) const;
    
    /// Comparison: sorts by increasing threshold.
    bool operator< (TriggerGates const& other) const
      { return threshold() < other.threshold(); }
    
    
      private:
    ADCCounts_t fThreshold; ///< The threshold for all the gates.
    GateData_t fGates; ///< All the gates, at most one per channel.
    
  }; // class TriggerGates
  // --- END Data types --------------------------------------------------------
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
    fhicl::Atom<ADCCounts_t::value_t> Baseline {
      Name("Baseline"),
      Comment("the baseline level, fixed for all channels [ADC counts]"),
      0 // default: null baseline (!)
      };
    
    fhicl::Sequence<ADCCounts_t::value_t> ChannelThresholds {
      Name("ChannelThresholds"),
      Comment("triggering thresholds [ADC counts]")
      // mandatory
      };
    
    
  }; // struct Config
  // --- END Configuration -----------------------------------------------------
  
  
  /// Constructor: sets the configuration.
  TriggerGateBuilder(Config const& config);
  
  /// Virtual destructor. Nothing special about it, except that it's virtual.
  virtual ~TriggerGateBuilder() = default;
  
  /// Algorithm setup.
  virtual void setup(detinfo::DetectorTimings const& timings);
  
  /// Returns a collection of `TriggerGates` objects sorted by threshold.
  virtual std::vector<TriggerGates> build
    (std::vector<raw::OpDetWaveform> const& waveforms) const = 0;
  
  /// Returns all the configured thresholds.
  std::vector<ADCCounts_t> const& channelThresholds() const
    { return fChannelThresholds; }
  
  /// Returns the assumed waveform baseline [ADC counts].
  ADCCounts_t baseline() const { return fBaseline; }
  
  /// Returns the number of configured thresholds.
  std::size_t nChannelThresholds() const { return channelThresholds().size(); }
  
  
  /// Converts a time [&micro;s] into optical ticks.
  optical_tick timeToOpticalTick(microsecond time) const;
  
  /// Converts a timestamp from `raw::OpDetWaveform` into optical ticks.
  optical_tick timeStampToOpticalTick(raw::TimeStamp_t time) const;

  
    protected:
  
  /// Returns a detector timings object.
  detinfo::DetectorTimings const& detTimings() const { return *fDetTimings; }
  
  /// Returns the absolute thresholds,
  /// in the same order as `channelThresholds()`.
  std::vector<ADCCounts_t> const& absoluteThresholds() const
    { return fAbsoluteThresholds; }
  
  
  /// Creates an empty TriggerGates object for each threshold;
  /// thresholds are kept relative.
  std::vector<TriggerGates> prepareAllGates() const;

  
    private:
  
  // --- BEGIN Configuration parameters ----------------------------------------
  ADCCounts_t fBaseline; ///< The baseline for all channels.
  
  /// All single channel thresholds, sorted in increasing order.
  std::vector<ADCCounts_t> fChannelThresholds;
  
  // --- END Configuration parameters ------------------------------------------
  
  
  // --- BEGIN Setup -----------------------------------------------------------
  
  /// LArSoft detector timing utility.
  std::optional<detinfo::DetectorTimings> fDetTimings;
  
  // --- END Setup -------------------------------------------------------------
  
  
  // --- BEGIN Cache -----------------------------------------------------------
  
  /// Trigger thresholds in absolute ADC counts.
  std::vector<ADCCounts_t> fAbsoluteThresholds;
  
  // --- END Cache -------------------------------------------------------------
  
  
  /// Converts the thresholds from relative to baseline to absolute.
  std::vector<ADCCounts_t> makeAbsoluteThresholds() const;
  
  
}; // class icarus::trigger::TriggerGateBuilder


//------------------------------------------------------------------------------
namespace icarus::trigger::details {
  
  // Category used for debugging information output
  constexpr auto TriggerGateDebugLog = "TriggerSimulation";
  
} // namespace icarus::trigger::details


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void icarus::trigger::TriggerGateBuilder::TriggerGates::dump(Stream& out) const
{
  
  auto const nOpenGates = fGates.size() - std::count_if(
    fGates.begin(), fGates.end(), std::mem_fn(&triggergate_t::alwaysClosed)
    );
  out << nOpenGates << "/" << fGates.size() << " trigger gates on threshold "
    << fThreshold;
  if (nOpenGates) {
    out << ":";
    for (auto const& gate: fGates) {
      if (gate.alwaysClosed()) continue;
      out << "\n  " << gate;
    }
  }
  out << "\n";
} // icarus::trigger::TriggerGateBuilder::TriggerGates::dump()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERGATEBUILDER_H
