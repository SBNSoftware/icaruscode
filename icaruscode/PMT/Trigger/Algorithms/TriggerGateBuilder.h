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
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h" // gatesIn()
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// C/C++ standard libraries
#include <vector>
#include <tuple>
#include <algorithm> // std::count_if()
#include <functional> // std::mem_fn()
#include <optional>
#include <cstddef> // std::size_t


namespace icarus::trigger {
  
  // ---------------------------------------------------------------------------
  //
  // declarations
  //
  
  struct WaveformWithBaseline;
  
  class TriggerGateBuilder;
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/// Object to carry around waveform ant its baseline.
struct icarus::trigger::WaveformWithBaseline
  : std::tuple<raw::OpDetWaveform const*, icarus::WaveformBaseline const*>
{
  
  using Waveform_t = raw::OpDetWaveform;
  using Baseline_t = icarus::WaveformBaseline;
  using Base_t = std::tuple<Waveform_t const*, Baseline_t const*>;
  
  using Base_t::Base_t; // inherit constructors
  
  /// Returns a reference to the waveform.
  Waveform_t const& waveform() const { return *waveformPtr(); }
  
  /// Returns a reference to the waveform baseline.
  Baseline_t const& baseline() const { return *baselinePtr(); }
  
  /// Returns a pointer to the waveform.
  Waveform_t const* waveformPtr() const
    { return std::get<Waveform_t const*>(*this); }
  
  /// Returns a pointer to the waveform baseline.
  Baseline_t const* baselinePtr() const
    { return std::get<Baseline_t const*>(*this); }
  
  /// Returns whether the baseline is available for this waveform.
  bool hasBaseline() const { return baselinePtr() != nullptr; }
  
  
  // questionable practices...
  operator Waveform_t const& () const { return waveform(); }
  operator Waveform_t const* () const { return waveformPtr(); }
  operator Baseline_t const& () const { return baseline(); }
  operator Baseline_t const* () const { return baselinePtr(); }
  
}; // struct icarus::trigger::WaveformWithBaseline


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
  
  using Channel_t = raw::Channel_t;
  
  /// Mnemonic for an invalid optical detector channel.
  static constexpr Channel_t InvalidChannel
    = std::numeric_limits<Channel_t>::max();
  
  
  /// Container of logical gates for all triggering channels for a threshold.
  class TriggerGates {
    
      public:
    
    // we need a gate type able to track its source, and the source here is
    // unequivocally a waveform
    using triggergate_t
      = icarus::trigger::TrackedOpticalTriggerGate<raw::OpDetWaveform>;
    using GateData_t = std::vector<triggergate_t>;
    
    
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
    
    /// Returns the gate for the specified waveform `channel`, `nullptr` if n/a.
    triggergate_t const* getGateFor(raw::Channel_t const channel) const;
    
    /// Returns (and creates, if necessary) the gate for the specified waveform.
    triggergate_t& gateFor(raw::OpDetWaveform const& waveform);
    
    /// Dumps the content of this set of gates into the `out` stream.
    template <typename Stream>
    void dump(Stream& out) const;
    
    /// Comparison: sorts by increasing threshold.
    bool operator< (TriggerGates const& other) const
      { return threshold() < other.threshold(); }
    
    
    static constexpr bool isValidChannel(Channel_t channel)
      { return channel != InvalidChannel; }
    
      private:
    ADCCounts_t fThreshold; ///< The threshold for all the gates.
    GateData_t fGates; ///< All the gates, at most one per channel.
    
    /// Returns an iterator to the gate for `channel`, `fGates.end()` if none.
    auto findGateFor(raw::Channel_t const channel);
    
    /// Returns an iterator to the gate for `channel`, `fGates.end()` if none.
    auto findGateFor(raw::Channel_t const channel) const;
    
  }; // class TriggerGates
  // --- END Data types --------------------------------------------------------
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
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
  
  /// Algorithm reset. It will require a new setup before using it again.
  virtual void reset() {}
  
  /// Resets and sets up.
  virtual void resetup(detinfo::DetectorTimings const& timings)
    { reset(); setup(timings); }
  
  /// Resets and sets up (including a new set of thresholds).
  virtual void resetup(
    detinfo::DetectorTimings const& timings,
    std::vector<ADCCounts_t> const& thresholds
    )
    { resetup(timings); doSetThresholds(thresholds); }
  
  /// Returns a collection of `TriggerGates` objects sorted by threshold.
  virtual std::vector<TriggerGates> build
    (std::vector<WaveformWithBaseline> const& waveforms) const = 0;
  
  /// Returns all the configured thresholds.
  std::vector<ADCCounts_t> const& channelThresholds() const
    { return fChannelThresholds; }
  
  /// Returns the number of configured thresholds.
  std::size_t nChannelThresholds() const { return channelThresholds().size(); }
  
  
  /// Converts a time [&micro;s] into optical ticks.
  optical_tick timeToOpticalTick(microsecond time) const;
  
  /// Converts a timestamp from `raw::OpDetWaveform` into optical ticks.
  optical_tick timeStampToOpticalTick(raw::TimeStamp_t time) const;

  
  /// Returns whether `channel` is valid.
  static constexpr bool isValidChannel(Channel_t channel)
    { return channel != InvalidChannel; }
  
    protected:
  
  /// Returns a detector timings object.
  detinfo::DetectorTimings const& detTimings() const { return *fDetTimings; }
  
  /// Creates an empty TriggerGates object for each threshold;
  /// thresholds are kept relative.
  std::vector<TriggerGates> prepareAllGates() const;

  /// Sets all thresholds anew.
  virtual void doSetThresholds(std::vector<ADCCounts_t> const& thresholds)
    { fChannelThresholds = thresholds; }
  
    private:
  
  // --- BEGIN Configuration parameters ----------------------------------------
  
  /// All single channel thresholds, sorted in increasing order.
  std::vector<ADCCounts_t> fChannelThresholds;
  
  // --- END Configuration parameters ------------------------------------------
  
  
  // --- BEGIN Setup -----------------------------------------------------------
  
  /// LArSoft detector timing utility.
  std::optional<detinfo::DetectorTimings> fDetTimings;
  
  // --- END Setup -------------------------------------------------------------
  
  
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
    fGates.begin(), fGates.end(),
    [](auto const& gate){ return gate.gate().alwaysClosed(); }
    );
  out << nOpenGates << "/" << fGates.size() << " trigger gates on threshold "
    << fThreshold;
  if (nOpenGates) {
    out << ":";
    for (auto const& gate: icarus::trigger::gatesIn(fGates)) {
      if (gate.alwaysClosed()) continue;
      out << "\n  " << gate;
    }
  }
  out << "\n";
} // icarus::trigger::TriggerGateBuilder::TriggerGates::dump()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERGATEBUILDER_H
