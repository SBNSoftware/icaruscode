/**
 * @file    icaruscode/Utilities/FastAndPoorGauss.h
 * @brief   Fast approximate Gaussian random translator.
 * @date    February 15, 2020
 */

#ifndef ICARUSCODE_UTILITIES_FASTANDPOORGAUSS_H
#define ICARUSCODE_UTILITIES_FASTANDPOORGAUSS_H

// ROOT
#include "TMath.h" // TMath::ErfInverse()

// C/C++ standard library
#include <array>
#include <type_traits> // std::is_integral_v, std::is_signed_v
#include <cmath> // std::sqrt()
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace util::details {
  
  // ---------------------------------------------------------------------------
  /// Returns whether the integral type `value` is a power of 2.
  template <typename T>
  constexpr bool isPowerOfTwo(T value)
    {
      static_assert(std::is_integral_v<T>);
      if constexpr(std::is_signed_v<T>) if (value < 0) return false;
      do {
        if (value & 1) return (value == 1);
      } while (value >>= 1);
      return false;
    } // isPowerOfTwo()
  
  
  // ---------------------------------------------------------------------------
  
  
} // namespace util::details


// -----------------------------------------------------------------------------
namespace util {
  template <std::size_t N, typename T>
  class FastAndPoorGauss;
  
  template <typename T>
  class GaussianTransformer;
  
  template <typename T>
  class UniformSequence;

} // namespace util

// -----------------------------------------------------------------------------
/**
 * @brief Translates a number _u_ uniformly distributed between 0 and 1
 *        into a Gaussian distributed one _z_.
 * @param N number of possible values returned; it *must* be a power of 2
 * @param T (default: `double`) type of number for _u_ and _z_
 * 
 * The focus of this algorithm is speed, which is paid with some memory and a
 * lot of precision.
 * The algorithm picks one of only `N` precomputed possible values for each
 * input. The mapping of @f$ u \in [ 0, +1 [ @f$ into a real number is a
 * multi-step function with `N` steps. The steps are located at @f$ 1/2N @f$
 * steps through the domain of _u_.
 * 
 * The returned number _z_ is distributed according to a standard normal
 * distribution with mean 0 and standard deviation 1. The number can be turned
 * into one distributed with arbitrary mean @f$ \mu @f$ and standard deviation
 * @f$ \sigma @f$ with the usual transformation @f$ x = \mu + \sigma z @f$
 * (see e.g. `util::GaussianTransformer`).
 * 
 * 
 * Resources
 * ----------
 * 
 * Math is internally performed in `T` precision, except for the initialization
 * that is performed in double precision. The precomputed value table is
 * `sizeof(T) * N` bytes large.
 * 
 * 
 * @note This class is tuned for performance, and it allocates its data on the
 *       stack. Stack overflows have been observed in Linux with a size of
 *       `N` 2^20^. Stack overflows are very puzzling since they do not present
 *       any standard diagnostics and debuggers may not notice them.
 *       Bottom line is: do not overdo with the number of samples, or ask the
 *       author to provide a special version using dynamic memory allocation.
 * 
 */
template <std::size_t N, typename T = double>
class util::FastAndPoorGauss {
  static_assert(
    util::details::isPowerOfTwo(N), "Template parameter N must be a power of 2."
    );
  
    public:
  using Data_t = T; ///< Type of data to deal with.
  
  static constexpr std::size_t NPoints = N; ///< Number of sampled points.
  
  //@{
  /// Returns the Gaussian distributed value corresponding to `u`.
  Data_t transform(Data_t const u) const { return fSamples[indexOf(u)]; }
  Data_t operator() (Data_t const u) const { return transform(u); }
  //@}
  
    private:
  /// Sampled points of inverse Gaussian.
  static std::array<Data_t, N> const fSamples;
  
  /// Returns the index of the precomputed table serving the value `u`.
  std::size_t indexOf(Data_t u) const;
  
  /// Fills the pre-sampling table.
  static std::array<Data_t, NPoints> makeSamples();
  
}; // util::FastAndPoorGauss<>



// -----------------------------------------------------------------------------
/**
 * @brief Transforms a standard normal number into one on a different normal
 *        distribution.
 * @tparam T type of the data
 * 
 * This functor applies the simple mapping of a standard normal variable with
 * mean 0 and standard deviation 1 into one with arbitrary mean and standard
 * deviation: @f$ x = \mu + \sigma z @f$.
 * 
 */
template <typename T>
class util::GaussianTransformer {
  
    public:
  using Data_t = T; ///< Type of data to deal with.
  
  /// Constructor: selects the mean and standard deviation.
  GaussianTransformer(T mean, T stddev): fMean(mean), fStdDev(stddev) {}
  
  //@{
  /// Transforms normal value `z` into the target distribution.
  Data_t transform(Data_t const z) const
    { return transform(z, mean(), stdDev()); }
  Data_t operator() (Data_t const z) const { return transform(z); }
  //@}
  
  /// Returns the mean of the target Gaussian distribution.
  Data_t mean() const { return fMean; } 
  
  /// Returns the standard deviation of the target Gaussian distribution.
  Data_t stdDev() const { return fStdDev; } 
  
  /// Transforms normal value `z` into one distributed with `mean` and `stddev`.
  static Data_t transform
    (Data_t const z, Data_t const mean, Data_t const stddev)
    { return mean + z * stddev; }
  
    private:
  
  Data_t const fMean   = Data_t{ 0.0 };
  Data_t const fStdDev = Data_t{ 1.0 };
  
}; // util::GaussianTransformer<>


//------------------------------------------------------------------------------
/**
 * @brief Samples the interval [ 0, 1 ] in sequence, cyclically.
 * @tparam T (default: `double`) type of number returned
 * 
 * This object `extract()`s in sequence numbers in the interval [ 0 ; 1 ].
 * The interval is split in equally sized bins, and each time `extract()`
 * returns the center of the next bin, starting with the first one.
 * After reaching the last bin, it restarts the cycle.
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::array<double, 100N> seq1, seq2;
 * 
 * util::UniformSequence<> extract { 100 };
 * for (auto& v: seq) v = extract();
 * 
 * std::generate(seq.begin(), seq.end(), util::UniformSequence<>{ 100 });
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * (the two sequences are identical).
 */
template <typename T = double>
class util::UniformSequence {
    public:
  using Data_t = T; ///< Type of data required.
  
  /// Initializes a sequence of N steps.
  UniformSequence(unsigned int N)
    : fN(N), fStep(Data_t{ 1 }/fN), fOffset(fStep/Data_t{ 2 })
    {}
  
  //@{
  /// Returns the next number in the sequence.
  Data_t extract() { return fOffset + next() * fStep; }
  Data_t operator() () { return extract(); }
  //@}
  
  /// Restarts the sequence.
  void reset() { fNext = 0U; }
  
    private:
  unsigned int const fN; ///< The number of points the interval is split into.
  
  Data_t const fStep; ///< Size of each step.
  Data_t const fOffset; ///< Offset of the values (half step).
  
  unsigned int fNext = 0U; ///< The next point to be delivered.
  
  /// Returns the current point and prepares for the next one.
  unsigned int next();
  
}; // class UniformSequence<>


// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
template <std::size_t N, typename T>
std::array<T, N> const util::FastAndPoorGauss<N, T>::fSamples
  = util::FastAndPoorGauss<N, T>::makeSamples();


// -----------------------------------------------------------------------------
template <std::size_t N, typename T>
std::size_t util::FastAndPoorGauss<N, T>::indexOf(Data_t u) const {
  return static_cast<std::size_t>(u * NPoints);
} // util::FastAndPoorGauss<>::indexOf()


// -----------------------------------------------------------------------------
template <std::size_t N, typename T>
auto util::FastAndPoorGauss<N, T>::makeSamples() -> std::array<Data_t, NPoints>
{
  
  std::array<Data_t, NPoints> samples;
  
  double const V2 = std::sqrt(2.0);
  
  util::UniformSequence<Data_t> extract { NPoints };
  
  for(Data_t& value: samples)
    value = static_cast<Data_t>(TMath::ErfInverse(extract() * 2.0 - 1.0) * V2);
  
  return samples;
} // util::FastAndPoorGauss<>::makeSamples()


// -----------------------------------------------------------------------------
// ---  util::UniformSequence
// -----------------------------------------------------------------------------
template <typename T>
unsigned int util::UniformSequence<T>::next() {
  unsigned int i = fNext;
  if (++fNext == fN) reset();
  return i;
} // util::UniformSequence<T>::next()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_FASTANDPOORGAUSS_H
