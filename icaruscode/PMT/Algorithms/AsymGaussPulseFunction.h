/**
 * @file   icaruscode/PMT/Algorithms/AsymGaussPulseFunction.h
 * @brief  Pulse from one photoelectron as two half Gaussian functions.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 4, 2020
 *
 * This library is header only.
 *
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_ASYMGAUSSPULSEFUNCTION_H
#define ICARUSCODE_PMT_ALGORITHMS_ASYMGAUSSPULSEFUNCTION_H

// library header
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electronics.h" // counts_f

// C++ standard library
#include <ostream> // std::ostream
#include <string>
#include <utility> // std::forward()
#include <cmath> // std::exp()


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename T> class AsymGaussPulseFunction;
}

// -----------------------------------------------------------------------------
/**
 * @brief Describes the waveform from a single photoelectron.
 * @tparam Time type of time unit to be used
 *
 * This functor (class behaving like a function) describes the shape of the
 * response to a single photoelectron as an asymmetric Gaussian shape.
 */
template <typename T>
class icarus::opdet::AsymGaussPulseFunction
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
   * @param sigmaLeft the standard deviation of the shape before transition
   * @param sigmaRight the standard deviation of the shape after transition
   *
   * The time parameters (`peakTime`, `sigmaLeft` and `sigmaRight`) must be
   * measured in same unit. The `peakTime` defined the position of the shape
   * with respect to time 0.
   *
   */
  AsymGaussPulseFunction(
    ADCcount amplitude,
    Time peakTime,
    Time sigmaLeft,
    Time sigmaRight
    )
    : fAmplitude(amplitude)
    , fTransitTime(peakTime)
    , fSigmaL(sigmaLeft)
    , fSigmaR(sigmaRight)
    {}

  /// @{
  /// @name Parameter accessors.

//  Time peakTime() const { return fTransitTime; } // inherited
  Time leftSigma() const { return fSigmaL; }
  Time rightSigma() const { return fSigmaR; }
  ADCcount amplitude() const { return fAmplitude; }

  /// @}

  /// Returns the value of normal distribution at specified point.
  static ADCcount Gaussian(Time x, Time mean, Time sigma, ADCcount amplitude)
    { return amplitude * std::exp(-sqr((x - mean)/sigma)/2.0); }

    private:
  ADCcount fAmplitude; ///< Amplitude at peak (transition).
  Time fTransitTime; ///< Time of transition between the two forms of shape.
  Time fSigmaL; ///< RMS parameter of the shape before transition.
  Time fSigmaR; ///< RMS parameter of the shape after transition.

  /// Returns the time of the center of the Gaussian.
  Time myPeakTime() const { return fTransitTime; }

  /// Returns the amplitude of the pulse frome the baseline, including its sign.
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

  template <typename V>
  static V sqr(V value) { return value * value; }

  template <typename V>
  static V round(V value)
    {
      using namespace util::quantities::electronics_literals;
      return (value.abs() < 1e-6_ADCf)? 0.0_ADCf: value;
    }


}; // class icarus::opdet::AsymGaussPulseFunction<>


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::AsymGaussPulseFunction<T>::doEvaluateAt(Time time) const
  -> ADCcount
{
  return round(Gaussian(
    time,
    myPeakTime(),
    ((time < myPeakTime())? leftSigma(): rightSigma()),
    amplitude()
    ));
} // icarus::opdet::AsymGaussPulseFunction<T>::doEvaluateAt()


// -----------------------------------------------------------------------------
template <typename T>
void icarus::opdet::AsymGaussPulseFunction<T>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out
       << firstIndent << "Pulse shape: asymmetric Gaussian with peak at "
          << myPeakTime() << " and amplitude " << amplitude() << ":"
    << '\n' << indent
      << "  (t <  " << myPeakTime() << "): sigma " << leftSigma()
    << '\n' << indent
      << "  (t >= " << myPeakTime() << "): sigma " << rightSigma()
    << '\n';
} // icarus::opdet::AsymGaussPulseFunction<>::doDump()


// -----------------------------------------------------------------------------

#endif //  ICARUSCODE_PMT_ALGORITHMS_ASYMGAUSSPULSEFUNCTION_H
