/**
 * @file   icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h
 * @brief  Abstract interface of shape of a pulse from one photoelectron.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 4, 2020
 *
 * This is a header only library.
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PHOTOELECTRONPULSEFUNCTION_H
#define ICARUSCODE_PMT_ALGORITHMS_PHOTOELECTRONPULSEFUNCTION_H


// LArSoft libraries
#include "lardataalg/Utilities/quantities/electronics.h" // counts_f

// C++ standard library
#include <iosfwd> // std::ostream
#include <string>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::opdet {

  template <typename T> class PhotoelectronPulseFunction;
  template <typename Stream, typename T>
  Stream& operator<< (Stream&& out, PhotoelectronPulseFunction<T> const& pulse);

} // namespace icarus::opdet

// -----------------------------------------------------------------------------
/**
 * @brief Interface for a function describing a pulse from a photoelectron.
 * @tparam T type used to represent time
 *
 * A function implementing this interface expresses a sampled value in ADC
 * counts (`ADCcount` type) as function of time (`Time` type).
 *
 * The time types in LArSoft quantities library can be used as `Time` (`T`):
 * e.g., `PhotoelectronPulseFunction<util::quantities::nanosecond>`.
 *
 */
template <typename T>
class icarus::opdet::PhotoelectronPulseFunction {

    public:
  /// Type for ADC counts (floating point).
  using ADCcount = util::quantities::counts_f;

  using Time = T; ///< Type of time being used.


  /**
    * @brief Evaluates the pulse at the given time.
    * @param time time to evaluate the shape at
    * @return value of the pulse at the specified time
    */
  ADCcount evaluateAt(Time time) const { return doEvaluateAt(time); }

  /// Alias of `evaluateAt()`.
  ADCcount operator() (Time time) const { return evaluateAt(time); }

  /// Returns the time at which the first peak is found.
  Time peakTime() const { return doPeakTime(); }

  /// Returns the amplitude of the first peak in ADC counts.
  ADCcount peakAmplitude() const { return doPeakAmplitude(); }

  /// Returns the baseline of the pulse in ADC counts.
  ADCcount baseline() const { return doBaseline(); }

  /// Returns the polarity of the pulse (`+1`: positive, or `-1`: negative).
  int polarity() const { return doPolarity(); }


  // @{
  /**
    * @brief Prints on stream the parameters of this shape.
    * @param out the stream to write into
    * @param indent indentation string, prepended to all lines except first
    * @param indentFirst indentation string prepended to the first line
    */
  void dump(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const
    { doDump(out, indent, firstIndent); }
  void dump(std::ostream&& out, std::string const& indent = "") const
    { dump(out, indent, indent); }
  // @}


  // @{
  /**
    * @brief Returns the parameters of this shape as a descriptive string.
    * @param indent indentation string, prepended to all lines except first
    * @param indentFirst indentation string prepended to the first line
    * @return a string with the parameters of this shape
    */
  std::string toString
    (std::string const& indent, std::string const& firstIndent) const;
  std::string toString(std::string const& indent = "") const
    { return toString(indent, indent); }
  // @}


    protected:

  // --- BEGIN -- Mandatory customization --------------------------------------
  /// @name Mandatory customization
  /// @{

  /// Implementation of the function evaluation at `time`.
  virtual ADCcount doEvaluateAt(Time time) const = 0;

  /// Returns the time at which the first peak is found.
  virtual Time doPeakTime() const = 0;

  /// @}
  // --- END -- Mandatory customization ----------------------------------------


  // --- BEGIN -- Optional customization ---------------------------------------
  /// @name Optional customization
  /// @{

  /// Returns the amplitude of the first peak in ADC counts.
  virtual ADCcount doPeakAmplitude() const { return evaluateAt(peakTime()); }

  /// Returns the baseline of the pulse.
  virtual ADCcount doBaseline() const { return ADCcount{ 0 }; }

  /// Returns the polarity of the pulse (+1 or -1).
  virtual int doPolarity() const
    { return ((peakAmplitude() - baseline()) >= ADCcount{ 0 })? +1: -1; }


  /**
   * @brief Prints into the stream the parameters of this shape.
   * @param out the C++ output stream to write into
   * @param indent indentation string, prepended to all lines except first
   * @param indentFirst indentation string prepended to the first line
   */
  virtual void doDump(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const {}

  /// @}
  // --- END -- Optional customization -----------------------------------------


}; // class icarus::opdet::PhotoelectronPulseFunction<>


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template <typename T>
std::string icarus::opdet::PhotoelectronPulseFunction<T>::toString
  (std::string const& indent, std::string const& firstIndent) const
{
  std::ostringstream sstr;
  dump(sstr, indent, firstIndent);
  return std::move(sstr).str();
} // icarus::opdet::PhotoelectronPulseFunction<>::toString()


// -----------------------------------------------------------------------------
template <typename Stream, typename T>
Stream& icarus::opdet::operator<<
  (Stream&& out, icarus::opdet::PhotoelectronPulseFunction<T> const& pulse)
  { out << pulse.toString(); return out; }


// -----------------------------------------------------------------------------

#endif //  ICARUSCODE_PMT_ALGORITHMS_PHOTOELECTRONPULSEFUNCTION_H
