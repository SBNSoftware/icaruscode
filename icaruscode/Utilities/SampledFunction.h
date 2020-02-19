/**
 * @file   icaruscode/Utilities/SampledFunction.h
 * @brief  Class for a function with precomputed values.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 14, 2020
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_SAMPLEDFUNCTION_H
#define ICARUSCODE_UTILITIES_SAMPLEDFUNCTION_H

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"

// C++ core guideline library
#include "gsl/span"
#include "gsl/gsl_util" // gsl::index

// C++ standard library
#include <vector>
#include <functional> // std::function<>
#include <limits> // std::numeric_limits<>
#include <cmath> // std::isnormal()
#include <cassert>


// ---------------------------------------------------------------------------
namespace util {
  template <typename XType, typename YType> class SampledFunction;
} // namespace util

/**
 * @brief Precomputed discrete sampling of a given function.
 * @tparam XType (default: `double`) type of value accepted by the function
 * @tparam YType (default: as `XType`) type of value returned by the function
 *
 * This object contains the sampling of a specified function at regular values
 * of its variable.
 *
 * If the `size()` of the sampling is requested to be _N_, there will be a
 * sampling of _N_ values covering the specified range in steps of the same
 * length, last value excluded.
 * The sampling happens at the beginning of each step.
 *
 * In addition, subsampling can be requested. If _M_ subsamples are requested,
 * the first step is split in _M_ points and from each one a sampling of _N_
 * steps is started, causing overall _M N_ samples to be computed.
 *
 *
 * Requirements
 * -------------
 *
 * The function must be unary.
 *
 *
 * Technical note
 * ---------------
 *
 * The _M_ subsamples are stored each one contiguously.
 * Therefore a function with _M_ subsamples of size _N_ is different, at least
 * in storage, from a function with a single sampling (no subsamples) of size
 * _M N_.
 *
 */
template <typename XType = double, typename YType = XType>
class util::SampledFunction {
    public:

  using X_t = XType; ///< Type of value accepted by the function.
  using Y_t = YType; ///< Type of value returned by the function.
  using Function_t = std::function<X_t(Y_t)>; ///< Type of sampled function.

  /// Invalid index of sample, returned in case of error.
  static constexpr auto npos = std::numeric_limits<gsl::index>::max();
  
  /// Span of subsample data. Can be forward iterated.
  using SubsampleData_t = gsl::span<Y_t const>;

  // @{
  /**
   * @brief Constructor: samples `function` in the specified range.
   * @param function the function to be sampled
   * @param lower the lower limit of the range to be sampled
   * @param upper the upper limit of the range to be sampled
   * @param nSamples the number (_N_) of samples to be computed
   * @param subsamples (default: `1`) the number (_M_) of subsamples to be
   *                   computed
   *
   * The sampling of `function` is performed on `nSamples` points from `lower`
   * to `upper` (excluded).
   */
  SampledFunction(Function_t const& function,
    X_t lower, X_t upper,
    gsl::index nSamples,
    gsl::index subsamples = 1
    );

  template <typename Func>
  SampledFunction(Func function,
    X_t lower, X_t upper,
    gsl::index nSamples,
    gsl::index subsamples = 1
    )
    : SampledFunction
      (Function_t(function), lower, upper, nSamples, subsamples)
    {}
  // @}



  // --- BEGIN --- Query -------------------------------------------------------
  /// @name Query
  /// @{

  /// Returns the number of samples (in each subsample).
  gsl::index size() const { return fNSamples; }

  /// Returns the number of subsamples.
  gsl::index nSubsamples() const { return fNSubsamples; }

  /// Returns the lower limit of the covered range.
  X_t lower() const { return fLower; }

  /// Returns the upper limit of the covered range.
  X_t upper() const { return fUpper; }

  /// Returns the extension of the covered range.
  X_t rangeSize() const { return upper() - lower(); }

  /// Returns the extension of a step.
  X_t stepSize() const { return fStep; }

  /// Returns the base offset of the subsamples.
  X_t substepSize() const { return stepSize() / nSubsamples(); }

  /// @}
  // --- END --- Query ---------------------------------------------------------


  // --- BEGIN --- Access ------------------------------------------------------
  /// @name Access to the sampled data
  /// @{

  /// Returns the `iSample` value of the subsample with the specified index `n`.
  Y_t value(gsl::index iSample, gsl::index n = 0U) const
    { return subsampleData(n)[iSample]; }
  
  
  /// Returns the data of the subsample with the specified index `n`.
  SubsampleData_t subsample(gsl::index const n) const
    { return { subsampleData(n), static_cast<gsl::index>(fNSamples) }; }

  /**
   * @brief Returns the index of the step including `x`.
   * @param x the argument to the function
   * @param iSubsample the index of the subsample 
   * @return the index of step including `x`, or `npos` if none does
   *
   * This function returns the index of the sample whose step includes `x`.
   * A step includes its lower limit but not its upper limit, which usually
   * belongs to the next step (or it does not belong to any valid step).
   * If there is no step including `x`, the index of the would-be step is 
   * returned (it can be checked e.g. with `isValidStepIndex()`).
   */
  gsl::index stepIndex(X_t const x, gsl::index const iSubsample) const;

  /// Returns whether the specified step index is valid.
  bool isValidStepIndex(gsl::index const index) const
    { return (index >= 0) && (index < size()); }
  
  /**
   * @brief Returns the subsample closest to the value `x`.
   * @param x value to be checked
   * @return the index of the subsample found
   * 
   * The subsample with the bin including `x` whose lower bound is the closest
   * to `x` itself is returned.
   * 
   * For example, assuming bins aligned with 0 and a sampling with steps of
   * size 1 and 5 subsamples, there will be 5 bins contaning the value `x` 3.65:
   * [ 3.0, 4.0 ], [ 3.2, 4.2 ], [ 3.4, 4.4 ], [ 3.6, 4.6 ] and [ 2.8, 3.8 ],
   * one for each subsample: `closestSubsampleIndex(3.65)` will return the
   * sample with the bin [ 3.6, 4.6 ] (that is the fourth one, i.e. subsample
   * number 3), because its lower bound 3.6 is the closest to 3.65.
   * 
   * The value `x` does not need to be in the sampling range. In the example
   * above, the range could have been between 0 and 2, and the result would be
   * the same.
   */
  gsl::index closestSubsampleIndex(X_t x) const;
  
  /// @}
  // --- END --- Access --------------------------------------------------------

  
    private:

  X_t fLower; ///< Lower limit of sampled range.
  X_t fUpper; ///< Upper limit of sampled range.
  gsl::index fNSamples; ///< Number of samples in the range.
  gsl::index fNSubsamples; ///< Number of subsamples.

  X_t fStep; ///< Step size.

  /// All samples, the entire first subsample first.
  std::vector<Y_t> fAllSamples;

  /// Returns the starting point of the subsample `n`.
  X_t subsampleOffset(gsl::index n) const
    { return lower() + substepSize() * n; }

  
  // @{
  /// Start of the block of values for subsample `n` (unchecked).
  Y_t const* subsampleData(gsl::index n) const
    { return fAllSamples.data() + fNSamples * n; }
  Y_t* subsampleData(gsl::index n) { return fAllSamples.data() + fNSamples * n; }
  // @}

  /// Computes the total size of the data.
  std::size_t computeTotalSize() const { return nSubsamples() * size(); }

  /// Samples the `function` and fills the internal caches.
  void fillSamples(Function_t const& function);

  /// Returns `value` made non-negative by adding multiples of `range`.
  template <typename T>
  static T wrapUp(T value, T range);
  
}; // class SampledFunction<>


// =============================================================================
// ===  template implementation
// =============================================================================
template <typename XType, typename YType>
util::SampledFunction<XType, YType>::SampledFunction(
  Function_t const& function,
  X_t lower, X_t upper,
  gsl::index nSamples,
  gsl::index subsamples /* = 1 */
  )
  : fLower(lower)
  , fUpper(upper)
  , fNSamples(nSamples)
  , fNSubsamples(subsamples)
  , fStep(rangeSize() / fNSamples)
  , fAllSamples(computeTotalSize())
{
  assert(nSamples > 0);
  assert(subsamples > 0);
  fillSamples(function);
} // util::SampledFunction<>::SampledFunction()


// -----------------------------------------------------------------------------
template <typename XType, typename YType>
gsl::index util::SampledFunction<XType, YType>::stepIndex
  (X_t const x, gsl::index const iSubsample) const
{
  return static_cast<gsl::index>
    (std::floor((x - subsampleOffset(iSubsample)) / stepSize()));
} // gsl::index util::SampledFunction<XType, YType>::stepIndex()


// -----------------------------------------------------------------------------
template <typename XType, typename YType>
gsl::index util::SampledFunction<XType, YType>::closestSubsampleIndex
  (X_t const x) const
{
  return static_cast<gsl::index>
    (wrapUp(std::fmod(x - lower(), stepSize()), stepSize()) / substepSize());
} // gsl::index util::SampledFunction<XType, YType>::stepIndex()


// -----------------------------------------------------------------------------
template <typename XType, typename YType>
void util::SampledFunction<XType, YType>::fillSamples
  (Function_t const& function)
{

  /*
   * Plan:
   * 0. rely on the currently stored size specifications (range and samples)
   * 1. resize the data structure to the required size
   * 2. fill all the subsamples, in sequence
   *
   */

  //
  // 0. rely on the currently stored size specifications (range and samples)
  //
  std::size_t const dataSize = computeTotalSize();
  assert(dataSize > 0U);
  assert(fLower <= fUpper);
  assert(std::isnormal(fStep));

  //
  // 1. resize the data structure to the required size
  //
  fAllSamples.resize(dataSize);

  //
  // 2. fill all the subsamples, in sequence
  //
  auto iValue = fAllSamples.begin();
  for (gsl::index const iSubsample: util::counter(nSubsamples())) {
    X_t const offset = subsampleOffset(iSubsample);
    for (gsl::index const iStep: util::counter(size())) {
      X_t const x = offset + iStep * stepSize();
      Y_t const y = function(x);
      *iValue++ = y;
    } // for steps
  } // for subsamples

} // util::SampledFunction<>::fillSamples()


// -----------------------------------------------------------------------------
template <typename XType, typename YType>
template <typename T>
T util::SampledFunction<XType, YType>::wrapUp(T value, T range) {
  while (value < T{ 0 }) value += range;
  return value;
} // util::SampledFunction<>::wrapUp()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_SAMPLEDFUNCTION_H
