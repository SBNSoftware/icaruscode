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
#include "icaruscode/Utilities/WaveformOperations.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard library
#include <functional> // std::function
#include <utility> // std::cref()


// -----------------------------------------------------------------------------
// ---  icarus::opdet::DiscretePhotoelectronPulse
// -----------------------------------------------------------------------------
auto icarus::opdet::DiscretePhotoelectronPulse::sampleShape(
  PulseFunction_t const& pulseShape,
  gigahertz samplingFreq, unsigned int nSubsamples,
  ADCcount threshold
) -> SampledFunction_t
{
  using namespace util::quantities::time_literals;
  using namespace icarus::waveform_operations;

  // pick the function according to polarity;
  // the pulse polarity is included in the values,
  // the two functions (lambda) are of different type, so they are being wrapped
  // in the common `std::function` type

  auto const isBelowThreshold = (pulseShape.polarity() == +1)
    ? std::function(
      [baseline=pulseShape.baseline(), threshold](nanoseconds, ADCcount s)
        {
          return
            PositivePolarityOperations<ADCcount>::subtractBaseline(s, baseline)
              < threshold
            ;
        }
      )
    : std::function(
      [baseline=pulseShape.baseline(), threshold](nanoseconds, ADCcount s)
        {
          return
            NegativePolarityOperations<ADCcount>::subtractBaseline(s, baseline)
              < threshold
            ;
        }
      )
    ;

  return SampledFunction_t{
    std::cref(pulseShape), // function to sample (by reference because abstract)
    0.0_ns,                // sampling start time
    1.0 / samplingFreq,    // tick duration
    isBelowThreshold,      // when to stop the sampling
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

  bool const bLowOk = (low.abs() < limit);
  bool const bHighOk = (high.abs() < limit);
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
