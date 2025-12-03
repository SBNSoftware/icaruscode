/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulationTypes.h
 * @brief  Provides defintions for adder signal simulation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALSIMULATIONTYPES_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALSIMULATIONTYPES_H

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/ADCsettings.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond
#include "lardataalg/Utilities/quantities/frequency.h" // megahertz
#include "lardataalg/Utilities/quantities/electromagnetism.h" // volt...

// C/C++ standard libraries
#include <complex>
#include <ratio> // std::milli
#include <string>
#include <valarray>


// -----------------------------------------------------------------------------
/// Data types used in adder signal simulation.
namespace icarus::trigger::adder::types {
  
  /// Base data type used for waveform math.
  using Data_t = float;
  
  /// Type representing voltage for internal computation. Not a `Quantity`.
  using Voltage_t = Data_t;
  
  
  /// Relevant specifications of the PMT analog-to-digital converter.
  /// Defaults are from CAEN V1730B "default" settings.
  using ADCsettings_t = ::icarus::ADCsettings<Data_t>;
#if 0
  struct ADCsettings_t {
    using nanosecond = util::quantities::nanosecond;
    
    /// Type of voltage, as a `util::Quantity` type.
    using volt_t = util::quantities::volt_as<Voltage_t>;
    
    unsigned int bits { 14 }; ///< ADC bits.
    volt_t range { 2.0 }; ///< ADC range.
    nanosecond samplingTime { 2.0 }; ///< Duration of a sampling tick.
    
    /// Returns the ADC sampling rate.
    constexpr util::quantities::megahertz samplingRate() const;
    
    /// Returns the ADC sampling period.
    constexpr util::quantities::nanosecond samplingPeriod() const;
    
    // --- BEGIN ---  ADC/voltage conversions  ---------------------------------
    /// @name ADC/voltage conversions
    /// @{
    
    /// Returns the number of possible ADC values.
    constexpr unsigned int ADCrange() const { return 1 << bits; }
    
    /// ADC-to-millivolt conversion factor (multiply for ADC to get millivolt).
    constexpr Voltage_t ADC2mV() const;
    
    /// Millivolt-to-ADC conversion factor (multiply for millivolt to get ADC).
    constexpr Voltage_t mV2ADC() const { return 1.0 / ADC2mV(); }
    
    /// Converts `ADC` counts into millivolt.
    constexpr Voltage_t ADC2mV(short int ADC) const { return ADC * ADC2mV(); }
    
    /// Converts `millivolt` into ADC counts.
    constexpr Voltage_t mV2ADC(Voltage_t mV) const { return mV / mV2ADC(); }
    
    /// Converts `millivolt` into ADC counts.
    constexpr short int toADC(volt_t mV) const { return static_cast<short_int>(mV / range * ADCrange()); }
    
    /// @}
    // ---- END ----  ADC/voltage conversions  ---------------------------------
    
  }; // ADCsettings_t
#endif // 0
  
  
  /// Representation of samples in time domain for internal computation.
  using WaveformSamples_t = std::valarray<Voltage_t>;
  
  /// Representation of amplitudes in frequency domain for internal computation.
  using WaveformAmplitudes_t = std::valarray<std::complex<Voltage_t>>;
  
  /// Type for representing a frequency.
  using Frequency_t = std::complex<Data_t>;
  
  /// Representation of amplitudes in frequency domain for internal computation.
  using Frequencies_t = std::valarray<Frequency_t>;
  
  
}; // namespace icarus::trigger::adder::types

#if 0
// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------
namespace icarus::trigger {
  
  // these were not meant to be part of the public interface,
  // but they do have general usefulness
  
  /// Returns a quantity with the same value as `q` with a unit based on `Prefix`.
  template <typename Prefix, typename Quantity>
  constexpr std::enable_if_t<
    !util::quantities::concepts::is_interval_or_point<Quantity>::value,
    util::quantities::concepts::rescale<Quantity, Prefix>
    >
  scaleUnitTo(Quantity const& q)
    { return q; }
  
  /// Returns an interval with the same value as `iv` with a unit based on `Prefix`.
  template <typename Prefix, typename Interval>
  constexpr std::enable_if_t<
    util::quantities::concepts::is_interval<Interval>::value,
    util::quantities::concepts::rescale_interval<Interval, Prefix>
    >
  scaleUnitTo(Interval const& iv)
    { return iv; }
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
inline constexpr util::quantities::megahertz
icarus::trigger::adder::types::ADCsettings_t::samplingRate() const {
  return 1.0 / samplingTime;
}


// -----------------------------------------------------------------------------
inline constexpr util::quantities::nanosecond
icarus::trigger::adder::types::ADCsettings_t::samplingPeriod() const {
  return samplingTime;
}


// -----------------------------------------------------------------------------
inline constexpr auto icarus::trigger::adder::types::ADCsettings_t::ADC2mV
  () const -> Voltage_t
{
  return scaleUnitTo<std::milli>(range).value() / ADCrange();
}


// -----------------------------------------------------------------------------
constexpr short int icarus::trigger::adder::types::ADCsettings_t::toADC
  (volt_t mV) const
{
  return static_cast<short int>(mV / range * ADCrange()); 
}

#endif // 0
// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALSIMULATIONTYPES_H
