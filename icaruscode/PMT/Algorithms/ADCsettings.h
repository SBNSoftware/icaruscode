/**
 * @file   icaruscode/PMT/Algorithms/ADCsettings.h
 * @brief  Provides ADC settings.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_ADCSETTINGS_H
#define ICARUSCODE_PMT_ALGORITHMS_ADCSETTINGS_H


// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond
#include "lardataalg/Utilities/quantities/frequency.h" // megahertz
#include "lardataalg/Utilities/quantities/electromagnetism.h" // volt...
#include "lardataalg/Utilities/quantities/electronics.h" // counts

// C/C++ standard libraries
#include <ratio> // std::milli
#include <type_traits> // std::enable_if_t, ...


// -----------------------------------------------------------------------------
namespace icarus::details { template <typename Q> struct is_voltage; }
// -----------------------------------------------------------------------------
namespace icarus { template <typename T> struct ADCsettings; }
/**
 * @brief Relevant specifications of the PMT analog-to-digital converter.
 * @tparam T (default: `double`) the data type used for numerical conversions
 * 
 * Defaults are from CAEN V1730B "default" settings.
 */
template <typename T = double>
struct icarus::ADCsettings {
  
  /// Base data type used for waveform (voltage, ADC) math.
  using Data_t = T;
  
  using nanosecond = util::quantities::nanosecond;
  using millivolt = util::quantities::millivolt;
  
  /// Whether `Volt` is a voltage quantity (`util::quantities::units::Volt`).
  template <typename Volt>
  static constexpr bool is_voltage_v
    = util::quantities::concepts::is_quantity_v<Volt>
    && Volt::template sameBaseUnitAs<millivolt>()
    ;
  
  
  // --- BEGIN ---  Time axis  -------------------------------------------------
  /// @name Time and sampling
  /// @{
  
  nanosecond samplingTime { 2.0 }; ///< Duration of a sampling tick.
  
  /// Returns the ADC sampling rate.
  constexpr util::quantities::megahertz samplingRate() const;
  
  /// Returns the ADC sampling period.
  constexpr util::quantities::nanosecond samplingPeriod() const;
  
  /// @}
  // ---- END ----  Time axis  -------------------------------------------------
  
  
  // --- BEGIN ---  Conversion axis  -------------------------------------------
  /// @name Conversion axis
  /// @{
  
  /// Type representing voltage for internal computation. Not a `Quantity`.
  using Voltage_t = Data_t;
  
  /// Type representing integral ADC count. Not a `Quantity`.
  using ADC_t = short int;
  
  /// Type of voltage, as a `util::Quantity` type.
  using millivolt_t = util::quantities::millivolt_as<Voltage_t>;
  
  
  unsigned int bits = 14; ///< ADC bits.
  
  /// ADC voltage conversion range (0 to `range`).
  millivolt_t range = util::quantities::volt_as<Voltage_t>{ 2.0 };
  
  /// Number of possible ADC values (0 to `ADCrange() - 1`).
  constexpr ADC_t ADCrange() const { return static_cast<ADC_t>(1 << bits); }
  
  /// Rounds the `ADC` value to an integer (currently, by truncation).
  constexpr ADC_t roundADC(Data_t ADC) const { return static_cast<ADC_t>(ADC); }
  
  /// @}
  // ---- END ----  Conversion axis  -------------------------------------------
  
  
  // --- BEGIN ---  ADC/voltage conversions  ---------------------------------
  /// @name ADC/voltage conversions
  /// @{
  
  /// ADC-to-millivolt conversion factor (multiply for ADC to get millivolt).
  constexpr Voltage_t ADC2mV() const;
  
  /// Millivolt-to-ADC conversion factor (multiply for millivolt to get ADC).
  constexpr Voltage_t mV2ADC() const { return 1.0 / ADC2mV(); }
  
  /// Converts `ADC` counts into millivolt.
  constexpr Voltage_t ADC2mV(ADC_t ADC) const { return ADC * ADC2mV(); }
  
  /// Converts `millivolt` into ADC counts.
  constexpr Voltage_t mV2ADC(Voltage_t mV) const { return mV / mV2ADC(); }
  
  /// Converts `voltage` into ADC counts.
  template <typename R, typename V>
  constexpr short int to_ADC(util::quantities::scaled_volt<R, V> voltage) const;
  
  /// Converts `ADC` counts into millivolt.
  template <typename ADCcount>
  constexpr millivolt_t to_mV(ADCcount ADC) const;
  
  /// @}
  // ---- END ----  ADC/voltage conversions  ---------------------------------
    
}; // icarus::ADCsettings


// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------
namespace icarus {
  
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
template <typename T>
constexpr util::quantities::megahertz
icarus::ADCsettings<T>::samplingRate() const {
  return 1.0 / samplingTime;
}


// -----------------------------------------------------------------------------
template <typename T>
constexpr util::quantities::nanosecond
icarus::ADCsettings<T>::samplingPeriod() const {
  return samplingTime;
}


// -----------------------------------------------------------------------------
template <typename T>
constexpr auto icarus::ADCsettings<T>::ADC2mV() const -> Voltage_t
{
  return scaleUnitTo<std::milli>(range).value() / ADCrange();
}


// -----------------------------------------------------------------------------
template <typename T>
template <typename R, typename V>
constexpr short int icarus::ADCsettings<T>::to_ADC
  (util::quantities::scaled_volt<R, V> voltage) const
{
  return roundADC(voltage / range * ADCrange());
}


// -----------------------------------------------------------------------------
template <typename T>
template <typename ADCcount>
constexpr auto icarus::ADCsettings<T>::to_mV(ADCcount ADC) const -> millivolt_t
{
  if constexpr(util::quantities::concepts::is_quantity_v<ADCcount>) {
    static_assert(
      ADCcount::template sameBaseUnitAs<util::quantities::counts>(),
      "to_mV() must be called with a Counts quantity."
    );
    
    return to_mV(ADC.value());
  }
  else return static_cast<millivolt_t>(ADC * range / ADCrange()); 
} // icarus::ADCsettings::to_mV()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_ADCSETTINGS_H
