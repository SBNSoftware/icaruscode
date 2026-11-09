/**
 * @file   icaruscode/PMT/Algorithms/DiscretePhotoelectronPulse.cxx
 * @brief  Sampling of a photoelectron pulse (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/DiscretePhotoelectronPulse.h`
 *
 */

// library header
#include "icaruscode/PMT/Algorithms/DiscretePhotoelectronPulse.h"

// ICARUS libraries
#include "icaruscode/Utilities/quantities_utils.h" // util::...::abs()
#include "icarusalg/Utilities/WaveformOperations.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard library
#include <cmath> // std::abs()
#include <functional> // std::function
#include <utility> // std::cref()


// -----------------------------------------------------------------------------
namespace {
  
  /**
   * Returns `true` after seeing a given number of samples closer to baseline
   * than a specified threshold (either direction).
   */
  template <typename SampleType>
  class CloseToBaselineForClass {
    
    // configuration
    SampleType fMin, fMax;
    unsigned int fCountGoal;
    
    // cache
    unsigned int fCount = 0;
    
      public:
    
    /// Constructor: specifies `threshold` and `count` (baseline is in `ops`).
    template <typename WaveformOperations>
    CloseToBaselineForClass
      (SampleType threshold, unsigned int count, WaveformOperations ops)
      : fMin{ ops.shiftFromBaseline(-threshold) }
      , fMax{ ops.shiftFromBaseline(+threshold) }
      , fCountGoal{ count }
      {
        if (fMin > fMax) std::swap(fMin, fMax);
      }
    
    /// Returns whether there have been `fCountGoal` calls with `sample`
    /// close enough to the baseline (`closeEnough()`).
    template <typename TimeType>
    bool operator() (TimeType, SampleType sample)
      {
        if (closeEnough(sample)) ++fCount;
        else fCount = 0;
        return fCount >= fCountGoal;
      }
    
    /// Returns whether `sample` is close enough to the stored baseline.
    constexpr bool closeEnough(SampleType sample) const noexcept
      { return (sample >= fMin) && (sample <= fMax); }
    
  }; // CloseToBaselineForClass
  
} // local namespace


// -----------------------------------------------------------------------------
// ---  icarus::opdet::DiscretePhotoelectronPulse
// -----------------------------------------------------------------------------
icarus::opdet::DiscretePhotoelectronPulse::DiscretePhotoelectronPulse(
  PulseFunction_t const& pulseShape,
  gigahertz samplingFreq, unsigned int nSubsamples,
  ADCcount samplingThreshold, nanoseconds minTimeBelowThreshold
  )
  : fShape(pulseShape)
  , fSamplingFreq(samplingFreq)
  , fSampledShape(sampleShape(
      shape(), fSamplingFreq, nSubsamples,
      samplingThreshold,
      static_cast<int>(std::ceil(minTimeBelowThreshold * samplingFreq))
      ))
  {}


// -----------------------------------------------------------------------------
auto icarus::opdet::DiscretePhotoelectronPulse::sampleShape(
  PulseFunction_t const& pulseShape,
  gigahertz samplingFreq, unsigned int nSubsamples,
  ADCcount threshold, unsigned int minSamplesBelowThreshold
) -> SampledFunction_t
{
  using namespace util::quantities::time_literals;
  using namespace icarus::waveform_operations;

  // pick the function according to polarity;
  // the pulse polarity is included in the values,
  // the two functions (lambda) are of different type, so they are being wrapped
  // in the common `std::function` type
  
  auto isCloseToBaseline = (pulseShape.polarity() == +1)
    ? CloseToBaselineForClass{
      threshold, minSamplesBelowThreshold,
      PositivePolarityOperations<ADCcount>{ pulseShape.baseline() }
    }
    : CloseToBaselineForClass{
      threshold, minSamplesBelowThreshold,
      NegativePolarityOperations<ADCcount>{ pulseShape.baseline() }
    }
    ;

  return SampledFunction_t{
    std::cref(pulseShape), // function to sample (by reference because abstract)
    0.0_ns,                // sampling start time
    1.0 / samplingFreq,    // tick duration
    isCloseToBaseline,     // when to stop the sampling
    static_cast<gsl::index>(nSubsamples), // how many subsamples per tick
    pulseShape.peakTime() // sample at least until here
    };

} // icarus::opdet::DiscretePhotoelectronPulse::sampleShape()


// -----------------------------------------------------------------------------
bool icarus::opdet::DiscretePhotoelectronPulse::checkRange
  (ADCcount limit, std::string const& outputCat /* = "" */) const
{
  assert(pulseLength() > 0);
  auto const low = *(fSampledShape.subsample(0).begin());
  auto const high
    = *(fSampledShape.subsample(fSampledShape.nSubsamples() - 1).rbegin());

  using std::abs;
  bool const bLowOk = (abs(low) < limit);
  bool const bHighOk = (abs(high) < limit);
  if (bLowOk && bHighOk) return true;
  if (!outputCat.empty()) {
    mf::LogWarning log(outputCat);
    log << "Check on sampled photoelectron waveform template failed!";
    if (!bLowOk) {
      log << "\n => low tail at the starting of sampling is already " << low;
    }
    if (!bHighOk) {
      log
        << "\n => high tail at the end of sampling ("
          << duration() << ") is still at " << high
        ;
    }
    log << "\nShape parameters:" << shape().toString("  ", "");
  } // if writing a message on failure
  return false;
} // icarus::opdet::DiscretePhotoelectronPulse::checkRange()


// -----------------------------------------------------------------------------
