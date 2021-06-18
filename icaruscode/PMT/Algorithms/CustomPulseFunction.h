/**
 * @file   icaruscode/PMT/Algorithms/CustomPulseFunction.h
 * @brief  Pulse from one photoelectron fully defined by the configuration.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 18, 2020
 *
 * This library is header only.
 *
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_CUSTOMPULSEFUNCTION_H
#define ICARUSCODE_PMT_ALGORITHMS_CUSTOMPULSEFUNCTION_H

// library header
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electronics.h" // counts_f
#include "larcorealg/CoreUtils/counter.h"

// ROOT
#include "TFormula.h"

// C++ standard library
#include <stdexcept> // std::runtime_error
#include <ostream> // std::ostream
#include <algorithm> // std::set_difference()
#include <vector>
#include <set>
#include <iterator> // std::inserter
#include <string>
#include <utility> // std::forward(), std::pair
#include <cmath> // std::exp()


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename T> class CustomPulseFunction;
}

// -----------------------------------------------------------------------------
/**
 * @brief Describes the waveform from a single photoelectron.
 * @tparam Time type of time unit to be used
 *
 * This functor (class behaving like a function) describes the shape of the
 * response to a single photoelectron as a custom shape.
 * 
 * The shape is a formula interpreted via ROOT's `TFormula`.
 * It is assumed to start at time `0.0` (in `Time` units).
 * 
 */
template <typename T>
class icarus::opdet::CustomPulseFunction
  : public icarus::opdet::PhotoelectronPulseFunction<T>
{
  using Base_t = icarus::opdet::PhotoelectronPulseFunction<T>;

    public:
  /// Type for ADC counts (floating point).
  using ADCcount = typename Base_t::ADCcount;

  using Time = typename Base_t::Time; ///< Type of time being used.

  using ParameterValue_t = double; ///< Type of value for parameters (TFormula).
  
  /// Type of parameter name/value definition.
  using NameAndValue_t = std::pair<std::string, ParameterValue_t>;
  
  /// Type of list of all function parameters.
  using PulseParameters_t = std::vector<NameAndValue_t>;
  
  
  // @{
  /**
   * @brief Constructor: chooses shape and the values of its parameters.
   * @param expression mathematical expression of pulse shape
   * @param parameters the list of named parameters and their values
   * @param peakTime expression of time of the maximum amplitude of the shape
   *
   * The shape mathematical `expression` must be compatible with ROOT 6
   * `TFormula`.
   * The formula must have a single independent variable, time measured in
   * the `Time` scale. The value `0` of that time corresponds to when the
   * photoelectron is emitted from the photocathode.
   * 
   * The peak time argument is used to report when the peak happens.
   * It would be possible to compute it in a general way, but the implementation
   * of that functionality is complicate enough to get it right and fast, that
   * it's better to rely on user's knowledge. The value may be either a number
   * or a string. The number (`Time`) directly represents the time of the
   * peak, while the expression is another `TFormula` expression that can use
   * the same `parameters` as the shape `expression`, so that for example it is
   * possible to indicate the peak time of a shape like `exp(-(x - [mu])**2)`
   * as `[mu]`. In this case, the peak time expression must use no variable.
   *
   */
  CustomPulseFunction(
    std::string const& expression,
    PulseParameters_t const& parameters,
    Time peakTime
    );
  CustomPulseFunction(
    std::string const& expression,
    PulseParameters_t const& parameters,
    std::string const& peakTime
    );
  // @}

  /// @{
  /// @name Parameter accessors.

  /// Returns the value of the parameter specified by `name`.
  ParameterValue_t parameter(std::string const& name) const;

  /// @}


    private:
  
  /// Record of collected information on the pulse shape.
  struct PulseStats_t {
    bool negativePulse;     ///< Whether the pulse is considered negative.
    Time peakTime;          ///< Time of the pulse peak.
    ADCcount peakAmplitude; ///< Pulse amplitude at peak time.
  }; // PulseStats_t


  TFormula fFormula; ///< Formula of the shape.
  
  PulseStats_t fStats; ///< Collected information about the pulse.


  /// Constructor implementation.
  CustomPulseFunction(
    std::string const& expression,
    PulseParameters_t const& parameters
    );
  
  
  // --- BEGIN -- Interface implementation -------------------------------------
  /**
   * @brief Evaluates the pulse at the given time.
   * @param time time to evaluate the shape at
   *
   * The scale of the time is defined by the transition time passed
   * at construction.
   */
  virtual ADCcount doEvaluateAt(Time time) const override
    { return ADCcount::castFrom(fFormula.Eval(static_cast<double>(time))); }

  /// Returns the time at which the first peak is found.
  virtual Time doPeakTime() const override { return fStats.peakTime; }

  /// Returns the amplitude of the first peak in ADC counts.
  virtual ADCcount doPeakAmplitude() const override
    { return fStats.peakAmplitude; }

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
  
  
  /// Extracts statistics from pulse shape.
  PulseStats_t extractStats(Time peakTime) const;
  
  /// Returns a message about parameters in `formula` missing in `parameters`.
  /// @return a message about missing parameters, empty if none
  static std::string checkMissingParameters
    (TFormula const& formula, PulseParameters_t const& parameters);
  
  /// Returns a message about `parameters` which are not in `formula`.
  /// @return a message about excess parameters, empty if none
  static std::string checkExcessParameters
    (TFormula const& formula, PulseParameters_t const& parameters);
  
  /// Sets the value of the `parameters` of `formula`.
  static void setParameters
    (TFormula& formula, PulseParameters_t const& parameters);
  
}; // class icarus::opdet::CustomPulseFunction<>


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename T>
icarus::opdet::CustomPulseFunction<T>::CustomPulseFunction(
  std::string const& expression,
  PulseParameters_t const& parameters,
  Time peakTime
  )
  : CustomPulseFunction(expression, parameters)
{
  fStats = extractStats(peakTime);
} // icarus::opdet::CustomPulseFunction<>::CustomPulseFunction(Time)


// -----------------------------------------------------------------------------
template <typename T>
icarus::opdet::CustomPulseFunction<T>::CustomPulseFunction(
  std::string const& expression,
  PulseParameters_t const& parameters,
  std::string const& peakTime
  )
  : CustomPulseFunction(expression, parameters) 
{
  TFormula peakExpr("CustomPulseFunctionPeak", peakTime.c_str(), false);
  
  // excess parameters are allowed because they might belong to `expression`
  std::string const msg = checkMissingParameters(peakExpr, parameters);
  if (!msg.empty()) throw std::runtime_error("CustomPulseFunction:" + msg);
  
  setParameters(peakExpr, parameters);
  fStats = extractStats(Time{ peakExpr.Eval(0.0) }); // dummy value
  
} // icarus::opdet::CustomPulseFunction<>::CustomPulseFunction(string)


// -----------------------------------------------------------------------------
template <typename T>
icarus::opdet::CustomPulseFunction<T>::CustomPulseFunction(
  std::string const& expression,
  PulseParameters_t const& parameters
  )
  : fFormula{
      "CustomPulseFunction", // fixed name, hope ROOT does not mess it up...
      expression.c_str(),
      false // addToGlobList
    }
{
  std::string const msg = checkMissingParameters(fFormula, parameters)
    + checkExcessParameters(fFormula, parameters);
  if (!msg.empty()) throw std::runtime_error("CustomPulseFunction:" + msg);

  setParameters(fFormula, parameters);
} // icarus::opdet::CustomPulseFunction<>::CustomPulseFunction(impl)


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::CustomPulseFunction<T>::extractStats
  (Time peakTime) const -> PulseStats_t
{
  PulseStats_t stats;
  
  stats.peakTime = peakTime;
  stats.peakAmplitude = (*this)(peakTime);
  stats.negativePulse = stats.peakAmplitude < ADCcount{ 0 };
  
  return stats;
} // icarus::opdet::CustomPulseFunction<>::extractStats()


// -----------------------------------------------------------------------------
template <typename T>
void icarus::opdet::CustomPulseFunction<T>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out
    << firstIndent << "Custom pulse shape: " << fFormula.GetTitle()
    << "\n" << indent << "  peak " << fStats.peakAmplitude
      << " at " << fStats.peakTime
    ;
  if (fFormula.GetNpar() > 0) {
    out << "\n" << indent << "Parameters (" << fFormula.GetNpar() << "):";
    for (auto const iPar: util::counter(fFormula.GetNpar())) {
      out << "\n" << indent << "  [" << fFormula.GetParName(iPar) << "] = "
        << fFormula.GetParameter(iPar);
    } // for
  } // if parameters
  out << '\n';
} // icarus::opdet::CustomPulseFunction<>::doDump()


// -----------------------------------------------------------------------------
template <typename T>
std::string icarus::opdet::CustomPulseFunction<T>::checkMissingParameters
  (TFormula const& formula, PulseParameters_t const& parameters)
{
  std::set<std::string> required;
  for (auto const iPar: util::counter(formula.GetNpar()))
    required.insert(formula.GetParName(iPar));
  std::set<std::string> offered;
  for (auto const& nameAndValue: parameters)
    offered.insert(nameAndValue.first);
  
  std::set<std::string> missing;
  std::set_difference(
    required.cbegin(), required.cend(), offered.cbegin(), offered.cend(),
    std::inserter(missing, missing.begin())
    );
  if (missing.empty()) return {}; // success
  
  std::string msg
    { "\n * " + std::to_string(missing.size()) + " parameters missing:" };
  for (auto const& name: missing) msg += "\n   - '" + name + "'";
  
  return msg;
  
} // icarus::opdet::CustomPulseFunction<>::checkMissingParameters()


// -----------------------------------------------------------------------------
template <typename T>
std::string icarus::opdet::CustomPulseFunction<T>::checkExcessParameters
  (TFormula const& formula, PulseParameters_t const& parameters)
{
  std::set<std::string> offered;
  for (auto const& nameAndValue: parameters)
    offered.insert(nameAndValue.first);
  std::set<std::string> required;
  for (auto const iPar: util::counter(formula.GetNpar()))
    required.insert(formula.GetParName(iPar));
  
  std::set<std::string> excess;
  std::set_difference(
    offered.cbegin(), offered.cend(), required.cbegin(), required.cend(),
    std::inserter(excess, excess.begin())
    );
  if (excess.empty()) return {}; // success
  
  std::string msg
    { "\n * " + std::to_string(excess.size()) + " unused parameters:" };
  for (auto const& name: excess) msg += "\n   - '" + name + "'";
  
  return msg;
  
} // icarus::opdet::CustomPulseFunction<>::checkExcessParameters()


// -----------------------------------------------------------------------------
template <typename T>
void icarus::opdet::CustomPulseFunction<T>::setParameters
  (TFormula& formula, PulseParameters_t const& parameters)
{
  for (auto const& [ name, value ]: parameters) {
    auto const iPar = formula.GetParNumber(name.c_str());
    if (iPar < 0) continue;
    formula.SetParameter(iPar, value);
  } // for
} // icarus::opdet::CustomPulseFunction<>::setParameters()


// -----------------------------------------------------------------------------

#endif //  ICARUSCODE_PMT_ALGORITHMS_CUSTOMPULSEFUNCTION_H
