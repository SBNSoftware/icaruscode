/**
 * @file   icaruscode/PMT/Algorithms/DiscretePhotoelectronPulse.h
 * @brief  Sampling of a photoelectron pulse.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/DiscretePhotoelectronPulse.cxx`
 *
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_DISCRETEPHOTOELECTRONPULSE_H
#define ICARUSCODE_PMT_ALGORITHMS_DISCRETEPHOTOELECTRONPULSE_H


// // ICARUS libraries
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"
#include "icarusalg/Utilities/SampledFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond
#include "lardataalg/Utilities/quantities/frequency.h" // gigahertz
#include "lardataalg/Utilities/quantities/electronics.h" // tick, counts_f

// guidelines library
#include "gsl/gsl_util" // gsl::index

// C++ standard library
#include <string>
#include <utility> // std::forward()
#include <type_traits> // std::is_same_v, std::decay_t
#include <cstdlib> // std::size_t



// -----------------------------------------------------------------------------
namespace icarus::opdet {
  
  using namespace util::quantities::electronics_literals;
  
  class DiscretePhotoelectronPulse;
  
} // namespace icarus::opdet


// -----------------------------------------------------------------------------
/**
 * @brief Precomputed digitized shape of a given function.
 *
 * Multiple samplings ("subsamples") are performed with sub-tick offsets.
 * 
 * The sampled function is of type
 * `PhotoelectronPulseFunction<util::quantities::nanosecond>` (polymorphic
 * implementations are supported).
 *
 * A reference to the sampled function is kept available, so that function needs
 * to be valid for the lifetime of this object.
 * 
 */
class icarus::opdet::DiscretePhotoelectronPulse {
    public:
  using gigahertz = util::quantities::gigahertz;
  using nanoseconds = util::quantities::nanosecond;

  /// Type of shape (times are in nanoseconds).
  using PulseFunction_t = PhotoelectronPulseFunction<nanoseconds>;
  using ADCcount = PulseFunction_t::ADCcount;

    private:
  static_assert(
    std::is_same_v<std::decay_t<PulseFunction_t::Time>, nanoseconds>,
    "The type of single response function does not take nanoseconds."
    );

  /// Internal discretized representation of the sampled pulse shape.
  using SampledFunction_t = util::SampledFunction<nanoseconds, ADCcount>;

    public:
  using Time_t = nanoseconds;
  using Tick_t = util::quantities::tick;

  using SubsampleIndex_t = gsl::index; ///< Type of index of subsample.

  /// Type of subsample data (a sampling of the full range).
  using Subsample_t = SampledFunction_t::SubsampleData_t;

  static_assert(!std::is_same<Time_t, Tick_t>(),
    "Time and tick must be different types!");

  /**
    * @brief Constructor: samples the pulse.
    * @param pulseShape the shape to be pulsed; times in nanoseconds
    * @param samplingFreq frequency of samples [gigahertz]
    * @param nSubsamples (default: `1`) the number of samples within a tick
    * @param samplingThreshold (default: 10^-6^) pulse shape ends when its
    *        value is below this threshold
    *
    * Samples start from time 0, which is the time of the start of the first
    * tick. This time is expected to be the arrival time of the photoelectron.
    *
    * The length of the sampling is determined by the sampling threshold:
    * at the first tick after the peak time where the shape function is below
    * threshold, the sampling ends (that tick under threshold itself is also
    * discarded).
    *
    * The ownership of `pulseShape` is acquired by this object.
    */
  DiscretePhotoelectronPulse(
    PulseFunction_t const& pulseShape,
    gigahertz samplingFreq,
    unsigned int nSubsamples = 1U,
    ADCcount samplingThreshold = 1e-6_ADCf
    );

  /// Returns the length of the sampled pulse in ticks.
  std::size_t pulseLength() const { return fSampledShape.size(); }

  /// Evaluates the shape at the specified time.
  ADCcount operator() (Time_t time) const { return shape()(time); }

  // --- BEGIN -- Access to subsamples ---------------------------------------
  /**
    * @name Access to subsamples.
    *
    * A subsample is represented by a forward-iterable object. For example:
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
    * using namespace util::quantities::time_literals;
    * auto const& subsample = dpp.subsampleFor(5_ns);
    * for (auto sample: subsample) // ...
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *
    * @note Direct iterator access is currently not implemented because the
    *       underlying implementation of the subsample object (`gsl::span`)
    *       does not support "borrowing" the data. See the note in
    *       `util::SampledFunction` for more details. The issue should be
    *       solved with C++20.
    */
  /// @{

  /// Returns the subsample specified by index `i`, undefined if invalid.
  decltype(auto) subsample(SubsampleIndex_t i) const
    { return fSampledShape.subsample(i); }

  /// Returns the index of the subsample whose tick left limit is closest to
  /// `time` (see `util::SampledFunction::closestSubsampleIndex()`).
  decltype(auto) subsampleFor(Time_t time) const
    { return subsample(fSampledShape.closestSubsampleIndex(time)); }

  SubsampleIndex_t nSubsamples() const { return fSampledShape.nSubsamples(); }

  /// @}
  // --- END -- Access to subsamples -----------------------------------------

  // --- BEGIN -- Functional shape -------------------------------------------
  /// @name Functional shape
  /// @{
  /// Returns the function which was sampled.
  PulseFunction_t const& shape() const { return fShape; }

  /// Returns the peak amplitude in ADC counts.
  ADCcount peakAmplitude() const { return shape().peakAmplitude(); }

  /// Returns the time at the peak from the beginning of sampling.
  nanoseconds peakTime() const { return shape().peakTime(); }

  /// Returns the sampling frequency (same units as entered).
  gigahertz samplingFrequency() const { return fSamplingFreq; }

  /// Returns the sampling period (inverse of frequency).
  nanoseconds samplingPeriod() const { return 1.0 / samplingFrequency(); }

  /// Returns the duration of the waveform in time units.
  /// @see `pulseLength()`
  nanoseconds duration() const { return pulseLength() * samplingPeriod(); }

  /// @}
  // --- END -- Functional shape ---------------------------------------------


  // @{
  /**
    * @brief Prints on stream the parameters of this shape.
    * @tparam Stream type of stream to write into
    * @param out the stream to write into
    * @param indent indentation string, prepended to all lines except first
    * @param indentFirst indentation string prepended to the first line
    */
  template <typename Stream>
  void dump(Stream&& out,
    std::string const& indent, std::string const& firstIndent
    ) const;
  template <typename Stream>
  void dump(Stream&& out, std::string const& indent = "") const
    { dump(std::forward<Stream>(out), indent, indent); }
  // @}

  /**
    * @brief Checks that the waveform tails not sampled are negligible.
    * @param limit threshold below which the waveform is considered negligible
    * @param outputCat _(default: empty)_ message facility output category
    *        to use for messages
    * @return whether the two tails are negligible
    *
    * If `outputCat` is empty (default) no message is printed.
    * Otherwise, in case of failure a message is sent to message facility
    * (under category `outputCat`) describing the failure(s).
    */
  bool checkRange(ADCcount limit, std::string const& outputCat) const;

    private:

  /// Analytical shape of the pules.
  PulseFunction_t const& fShape;
  gigahertz fSamplingFreq; ///< Sampling frequency.

  /// Pulse shape, discretized.
  SampledFunction_t fSampledShape;


  /// Builds the sampling cache.
  static SampledFunction_t sampleShape(
    PulseFunction_t const& pulseShape,
    gigahertz samplingFreq, unsigned int nSubsamples,
    ADCcount threshold
    );

}; // class DiscretePhotoelectronPulse<>


// -----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//--- template implementation
//-----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// --- icarus::opdet::DiscretePhotoelectronPulse
// -----------------------------------------------------------------------------
inline icarus::opdet::DiscretePhotoelectronPulse::DiscretePhotoelectronPulse(
  PulseFunction_t const& pulseShape,
  gigahertz samplingFreq, unsigned int nSubsamples, /* = 1U */
  ADCcount samplingThreshold /* = 1e-3_ADCf */
  )
  : fShape(pulseShape)
  , fSamplingFreq(samplingFreq)
  , fSampledShape
    (sampleShape(shape(), fSamplingFreq, nSubsamples, samplingThreshold))
  {}


//-----------------------------------------------------------------------------
template <typename Stream>
void icarus::opdet::DiscretePhotoelectronPulse::dump(Stream&& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out << firstIndent << "Sampled pulse waveform " << pulseLength()
    << " samples long (" << duration()
    << " long, sampled at " << samplingFrequency()
    << ");"
    << "\n" << shape().toString(indent + "  ", indent);
  fSampledShape.dump(out, indent + "  ", indent);
} // icarus::opdet::DiscretePhotoelectronPulse::dump()


//-----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_ALGORITHMS_DISCRETEPHOTOELECTRONPULSE_H

