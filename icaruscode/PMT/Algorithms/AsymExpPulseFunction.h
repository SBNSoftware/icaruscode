/**
 * @file   icaruscode/PMT/Algorithms/AsymExpPulseFunction.h
 * @brief  Pulse from one photoelectron as sum of two exponential functions.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 20, 2020
 *
 * This library is header only.
 *
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_ASYMEXPPULSEFUNCTION_H
#define ICARUSCODE_PMT_ALGORITHMS_ASYMEXPPULSEFUNCTION_H

// library header
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electronics.h" // counts_f

// C++ standard library
#include <ostream> // std::ostream
#include <string>
#include <utility> // std::forward()
#include <cmath> // std::exp(), std::signbit(), std::isnormal()
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename T> class AsymExpPulseFunction;
}

// -----------------------------------------------------------------------------
/**
 * @brief Describes the waveform from a single photoelectron.
 * @tparam Time type of time unit to be used
 *
 * This functor (class behaving like a function) describes the shape of the
 * response to a single photoelectron as an exponentially raising and then
 * ecponentially falling shape:
 * @f[ A \left[ exp(-(t - t_{0})/\tau_{s}) - exp(-(t - t_{0})/\tau) \right] @f]
 * with @f$ \tau_{s} @f$ the rise time and @f$ \tau @f$ the fall time;
 * @f$ t_{0} @f$ is the start of the signal formation, i.e. the time at which
 * the signal starts raising.
 * 
 * Note that @f$ t_{0} @f$ is not the peak time, and @f$ A @f$ is not the peak
 * amplitude of the signal.
 */
template <typename T>
class icarus::opdet::AsymExpPulseFunction
  : public icarus::opdet::PhotoelectronPulseFunction<T>
{
  using Base_t = icarus::opdet::PhotoelectronPulseFunction<T>;

    public:
  /// Type for ADC counts (floating point).
  using ADCcount = typename Base_t::ADCcount;

  using Time = typename Base_t::Time; ///< Type of time being used.

  /**
   * @brief Constructor: assigns the parameters of the shape.
   * @param amplitude the maximum amplitudes of the shape (at transition)
   * @param peakTime the time of the maximum amplitude of the shape
   * @param raiseTau the raise time (exponential parameter) of the shape
   * @param fallTau the decay time (exponential parameter) of the shape
   *
   * The time parameters (`peakTime`, `raiseTau` and `fallTau`) must be
   * measured in same unit. The `peakTime` defines the position of the shape
   * with respect to time 0.
   * 
   * @note The construction parameters are *not* the ones described by the
   *       functional form documented in the class,
   *       @f$ A \left[ exp(-(t - t_{0})/\tau_{s}) - exp(-(t - t_{0})/\tau) \right] @f$:
   *       the `amplitude` parameter is the actual maximum of the signal
   *       response, and `peakTime` is the actual time at which the signal
   *       peaks.
   *       These parameters are more easily read from a measured waveform.
   *
   */
  AsymExpPulseFunction(
    ADCcount amplitude,
    Time peakTime,
    Time raiseTau,
    Time fallTau
    );

  /// @{
  /// @name Parameter accessors.

//  Time peakTime() const { return fPeakTime; } // inherited
  Time raiseTau() const { return fRaiseTau; }
  Time fallTau() const { return fFallTau; }
  ADCcount amplitude() const { return fAmplitude; }

  /// @}

  /// Returns the value of an exponential: `exp(-t/tau)`.
  static double Exponential(Time t, Time tau) { return std::exp(-t/tau); }

  /// Returns the value of an exponential: `exp(-t/tau)`.
  static double ExponentialDiff(Time t, Time raise, Time fall)
    { return Exponential(t, raise) - Exponential(t, fall); }

    private:
  ADCcount const fAmplitude; ///< Amplitude at peak (transition).
  Time const fPeakTime; ///< Time of maximum signal.
  Time const fRaiseTau; ///< Time constant of signal rise.
  Time const fFallTau; ///< Time constant of signal fall.
  
  Time const fRaiseTime; ///< Time when the signal starts rising above baseline.
  ADCcount const fA; ///< Normalization factor in the functional form.

  /// Returns the time of the peak.
  Time myPeakTime() const { return fPeakTime; }

  /// Returns the amplitude of the pulse from the baseline, including its sign.
  ADCcount myAmplitude() const { return fAmplitude; }


  // --- BEGIN -- Interface implementation -------------------------------------
  /**
   * @brief Evaluates the pulse at the given time.
   * @param time time to evaluate the shape at
   *
   * The scale of the time is defined by the transition time passed
   * at construction.
   */
  virtual ADCcount doEvaluateAt(Time time) const override;

  /// Returns the time at which the first peak is found.
  virtual Time doPeakTime() const override { return myPeakTime(); }

  /// Returns the amplitude of the first peak in ADC counts.
  virtual ADCcount doPeakAmplitude() const override { return myAmplitude(); }

  /**
   * @brief Prints on stream the parameters of this shape.
   * @param out the stream to write into
   * @param indent indentation string, prepended to all lines except first
   * @param indentFirst indentation string prepended to the first line
   */
  virtual void doDump(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const override;

  // --- END -- Interface implementation -------------------------------------

  /// Returns the time value at which `ExponentialDiff()` is maximum.
  static Time expDiffPeak(Time raiseTau, Time fallTau);
  
  
  template <typename V>
  static V round(V value) { return (value.abs() < V{ 1e-6 })? V{ 0.0 }: value; }


}; // class icarus::opdet::AsymExpPulseFunction<>


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename T>
icarus::opdet::AsymExpPulseFunction<T>::AsymExpPulseFunction(
  ADCcount amplitude,
  Time peakTime,
  Time raiseTau,
  Time fallTau
  )
  : fAmplitude(amplitude)
  , fPeakTime(peakTime)
  , fRaiseTau(raiseTau)
  , fFallTau(fallTau)
  , fRaiseTime(fPeakTime - expDiffPeak(fRaiseTau, fFallTau))
  , fA(fAmplitude/ExponentialDiff(fPeakTime - fRaiseTime, fRaiseTau, fFallTau))
{
  // numerical sanity checks (for debugging only)
  assert(std::isnormal(fA.value()));
} // icarus::opdet::AsymExpPulseFunction<>::AsymExpPulseFunction()


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::AsymExpPulseFunction<T>::doEvaluateAt(Time time) const
  -> ADCcount
{
  return (time <= fRaiseTime)
    ? ADCcount{ 0 }
    : round(fA * ExponentialDiff(time - fRaiseTime, fRaiseTau, fFallTau) );
} // icarus::opdet::AsymExpPulseFunction<>::doEvaluateAt()


// -----------------------------------------------------------------------------
template <typename T>
void icarus::opdet::AsymExpPulseFunction<T>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out << firstIndent
    << "Pulse shape: exponential sum with peak at "
      << myPeakTime() << " and amplitude " << amplitude()
      << ", time constants: " << raiseTau() 
      << " (raise) and " << fallTau() << " (fall)"
    << indent << "  start at " << fRaiseTime << ", normalization factor " << fA
    << '\n';
} // icarus::opdet::AsymExpPulseFunction<>::doDump()


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::AsymExpPulseFunction<T>::expDiffPeak
  (Time raiseTau, Time fallTau) -> Time 
{
  // the strange order of the operations makes sure that if Time defines
  // difference (into Time), ratio (into scalar) and multiplication to scalar,
  // the expression is valid.
  return
    raiseTau * (fallTau / (raiseTau - fallTau)) * std::log(raiseTau/fallTau);
} // icarus::opdet::AsymExpPulseFunction<T>::expDiffPeak()


// -----------------------------------------------------------------------------

#endif //  ICARUSCODE_PMT_ALGORITHMS_ASYMEXPPULSEFUNCTION_H
