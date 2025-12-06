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
#include "icaruscode/PMT/Trigger/Algorithms/WaveformWithBaseline.h"
#include "icaruscode/PMT/Algorithms/ADCsettings.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"

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


//------------------------------------------------------------------------------
namespace icarus::trigger { class TriggerGateBuilder; }
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
    
    
    fhicl::Sequence<std::string> ChannelThresholds {
      Name("ChannelThresholds"),
      Comment("triggering thresholds [voltage or ADC counts]")
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
  
  /**
   * @brief Returns a collection of `TriggerGates` objects sorted by threshold.
   * @param waveforms the digitized waveforms, together with their baseline
   * @return a list of trigger gates
   * 
   * The waveforms are grouped by channel, and for each channel
   * and threshold level a on/off trigger gate is produced.
   * 
   * The waveforms must be already sorted by channel and then by time.
   */
  virtual std::vector<TriggerGates> build
    (std::vector<WaveformWithBaseline> const& waveforms) const = 0;
  
  /// Returns all the configured thresholds, sorted from smaller to higher.
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
  
  // --- BEGIN ---  Utilities  -------------------------------------------------
  /**
   * @name Utilities
   * 
   * These are generic utilities that are used here but turned out to be of use
   * also in related contexts.
   * 
   * If they become useful in _unrelated_ contexts, they should win their own
   * header.
   */
  /// @{
  
  /// Channel number extraction for several optical-related data objects.
  struct ChannelExtractor {
    
    static constexpr raw::Channel_t channelOf(raw::Channel_t channel)
      { return channel; }
    
    static raw::Channel_t channelOf(raw::OpDetWaveform const& waveform)
      { return waveform.ChannelNumber(); }
    
    static raw::Channel_t channelOf(sbn::OpDetWaveformMeta const& waveformMeta)
      { return waveformMeta.ChannelNumber(); }
    
    static raw::Channel_t channelOf
      (icarus::trigger::OpticalTriggerGateData_t const& gate)
      { return gate.channel(); }
    
    template <typename Gate, typename OpDetInfo>
    static raw::Channel_t channelOf
      (icarus::trigger::TrackedTriggerGate<Gate, OpDetInfo> const& gate)
      { return channelOf(gate.gate()); }
    
    template <typename T>
    raw::Channel_t operator() (T&& obj) const { return channelOf(obj); }
    
  }; // struct ChannelExtractor
  
  // comparison using the channel number (special comparison operator)
  template <typename Comp = std::less<raw::Channel_t>>
  struct ChannelComparison {
    
    using Comparer_t = Comp;
    Comparer_t comp;
    
    constexpr ChannelComparison(Comparer_t comp = Comparer_t{})
      : comp{ std::move(comp) } {}
    
    template <typename A, typename B>
    bool operator() (A const& a, B const& b) const
      { ChannelExtractor channelOf; return comp(channelOf(a), channelOf(b)); }
    
  }; // struct ChannelComparison
  
  
  /**
   * @brief Converts a voltage or number string into ADC counts.
   * @tparam T type of the `ADCsettings` data
   * @param s the string to be converted
   * @param ADCsettings settings used to convert voltage to ADC
   * @return the value in ADC counts
   * 
   * If the string ends with the volt symbol, it is interpreted as a voltage
   * quantity specification. Otherwise, it is interpreted as a pure (real)
   * number representing an ADC value.
   * In the former case, `ADCsettings()` (default-constructed) is used to
   * convert it into ADC counts. Then, in both cases it is converted to integer
   * via `ADCsettings.roundADC()` (which may actually truncate it) and then
   * into `ADCCounts_t` quantity type.
   */
  template <typename T>
  static ADCCounts_t numberOrVoltageToADC
    (std::string const& s, icarus::ADCsettings<T> const& settings);
  
  /// Converts all threshold specifications with `numberOrVoltageToADC()`.
  template <typename T>
  static std::vector<ADCCounts_t> parseThresholds
    (std::vector<std::string> const& specs, icarus::ADCsettings<T> const& settings);
  
  /// @}
  // ---- END ----  Utilities  -------------------------------------------------
  
    protected:
  
  /// Returns a detector timings object.
  detinfo::DetectorTimings const& detTimings() const { return *fDetTimings; }
  
  /// Returns the settings of the ADC for threshold conversions.
  icarus::ADCsettings<> const& ADCsettings() const { return fADCsettings; }
  
  /// Creates an empty TriggerGates object for each threshold;
  /// thresholds are kept relative.
  std::vector<TriggerGates> prepareAllGates() const;

  /// Sets all thresholds anew.
  virtual void doSetThresholds(std::vector<ADCCounts_t> const& thresholds)
    { fChannelThresholds = thresholds; }
  
    private:
  
  // --- BEGIN Configuration parameters ----------------------------------------
  
  /// Settings of the ADC for threshold conversions (default values).
  icarus::ADCsettings<> fADCsettings;
  
  /// All single channel thresholds, sorted in increasing order.
  std::vector<ADCCounts_t> fChannelThresholds;
  
  // --- END Configuration parameters ------------------------------------------
  
  
  // --- BEGIN Setup -----------------------------------------------------------
  
  /// LArSoft detector timing utility.
  std::optional<detinfo::DetectorTimings> fDetTimings;
  
  // --- END Setup -------------------------------------------------------------
  
  /// Converts all threshold specifications with `numberOrVoltageToADC()`.
  std::vector<ADCCounts_t> parseThresholds
    (std::vector<std::string> const& specs) const
    { return parseThresholds(specs, fADCsettings); }
  
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
template <typename T>
auto icarus::trigger::TriggerGateBuilder::numberOrVoltageToADC
  (std::string const& s, icarus::ADCsettings<T> const& settings) -> ADCCounts_t
{
  short int ADC;
  
  std::string_view sv = s;
  while (!sv.empty()) {
    if (!std::isblank(sv.back())) break;
    sv.remove_suffix(1);
  }
  // C++20: sv.ends_with(unitStr)!!
  std::string_view const& unitStr = util::quantities::units::Volt::symbol;
  if (sv.substr(sv.length() - unitStr.length()) == unitStr) { // voltage unit
    auto const V = util::quantities::makeQuantity<util::quantities::volt>(sv);
    ADC = settings.to_ADC(V);
  }
  else { // ADC plain number
    std::size_t lastChar;
    ADC = settings.roundADC(stod(std::string{ sv }, &lastChar));
    if (lastChar != sv.length()) {
      throw std::invalid_argument
        { "Failed to convert '" + std::string{ sv } + "' to an ADC value." };
    }
  }
  return icarus::trigger::ADCCounts_t{ ADC };
} // icarus::trigger::TriggerGateBuilder::numberOrVoltageToADC()


//------------------------------------------------------------------------------
template <typename T>
auto icarus::trigger::TriggerGateBuilder::parseThresholds
  (std::vector<std::string> const& specs, icarus::ADCsettings<T> const& settings)
  -> std::vector<ADCCounts_t>
{
  std::vector<icarus::trigger::ADCCounts_t> thresholds;
  for (std::string const& spec: specs)
    thresholds.push_back(numberOrVoltageToADC(spec, settings));
  return thresholds;
}


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERGATEBUILDER_H
