/**
 * @file   icaruscode/PMT/Algorithms/SampledWaveformFunction.h
 * @brief  Pulse from one photoelectron as a train of samples.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 6, 2022
 *
 * This library is header only.
 *
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_SAMPLEDWAVEFORMFUNCTION_H
#define ICARUSCODE_PMT_ALGORITHMS_SAMPLEDWAVEFORMFUNCTION_H

// library header
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds

// C++ standard library
#include <ostream> // std::ostream
#include <vector>
#include <string>
#include <algorithm> // std::transform()
#include <numeric> // std::reduce()
#include <cmath> // std::round(), std::floor(), ...
#include <cstddef> // std::ptrdiff_t, std::size_t
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  using namespace util::quantities::time_literals; // ""_ns
  template <typename T> class SampledWaveformFunction;
}

// -----------------------------------------------------------------------------
/**
 * @brief Describes the waveform from a single photoelectron.
 * @tparam T type of time unit to be used
 *
 * This functor (class behaving like a function) describes the shape of the
 * response to a single photoelectron in a non-analytical form from a sequence
 * of samples.
 * 
 * The peak time is assigned to the sample with the largest value.
 * 
 * See more details in the constructor.
 * 
 * @note Currently the shape is hard-coded; it is possible to extend this object
 *       to receive the sample data from a text file with proper format
 *       (which needs to include all the information for `WaveformSpecs_t`).
 */
template <typename T>
class icarus::opdet::SampledWaveformFunction
  : public icarus::opdet::PhotoelectronPulseFunction<T>
{
  using Base_t = icarus::opdet::PhotoelectronPulseFunction<T>;

    public:
  /// Type for ADC counts (floating point).
  using ADCcount = typename Base_t::ADCcount;

  using Time = typename Base_t::Time; ///< Type of time being used.

  /// Specifies the waveform shape and features.
  struct WaveformSpecs_t {
    std::string name { "<unknown>" }; ///< Name of this waveform (short).
    std::string description;          ///< Description of this waveform.
    std::string date { "n/a" };  ///< An indication of the date of the waveform.
    unsigned int version { 1 };       ///< Version number.
    std::vector<float> samples;       ///< Samples [mV]
    Time sampleDuration;              ///< Time extension of one sample.
    float gain { 0.0 };               ///< Gain of this waveform.
  };
  
  /**
   * @brief Constructor: initializes from an hard-coded shape.
   * @param waveformSpecs all information on the single photoelectron response
   * @param peakTime time to assign to the peak
   * @param gain the gain of the optical detector
   * 
   * The shape described in waveformSpecs is shifted in time so that evaluation
   * at `peakTime` (`evaluateAt(peakTime)`) returns the peak amplitude; more
   * precisely, `peakTime` is set to match the start of the largest sample.
   * 
   * The `gain` is rescaled starting from the one in the waveform
   * specifications.
   * 
   * The polarity of the waveform is deduced by the value at peak.
   * 
   */
  SampledWaveformFunction(WaveformSpecs_t specs, Time peakTime, float gain);

  /// @{
  /// @name Parameter accessors.

  /// Returns the gain the waveform is representing.
  float gain() const { return fGain; }

  /// @}

    private:
  
  WaveformSpecs_t const fSource; ///< Waveform information.
  
  std::vector<ADCcount> const fSamples; ///< All samples.
  
  float const fGain; ///< The gain this waveform represents.
  
  std::size_t const fPeakSample; ///< The sample with the absolute peak.
  
  Time const fRefTime; ///< The time of the start of sample #0.
  
  /// The duration of each sample.
  Time sampleDuration() const { return fSource.sampleDuration; }
  
  
  // --- BEGIN -- Interface implementation -------------------------------------
  /**
   * @brief Evaluates the pulse at the given time.
   * @param time time to evaluate the shape at
   *
   * The scale of the time is defined by the peak time passed at construction.
   */
  virtual ADCcount doEvaluateAt(Time time) const override;

  /// Returns the time at which the first peak is found.
  virtual Time doPeakTime() const override
    { return fRefTime + fPeakSample * sampleDuration(); }

  /// Returns the amplitude of the first peak in ADC.
  virtual ADCcount doPeakAmplitude() const override
    { return fSamples[fPeakSample]; }

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
  
  
  /// Returns whether a sample with the specified index is within the range.
  bool hasSample(std::ptrdiff_t index) const
    { return (index >= 0) && (std::size_t(index) < fSamples.size()); }
  
  
  /// Returns the integral of the waveform.
  ADCcount integral() const;
  
  /**
   * @brief Transforms the input waveform.
   * @param waveform input waveform (in millivolt and for a known gain)
   * @param targetGain the desired gain
   * @return a sequence of samples in ADC
   * 
   * The returned waveform has the same time domain as the input one, but is
   * expressed in ADC instead of voltage (conversion is perfectly linear,
   * 2 V across 14 bits), and rescaled to meet the target gain.
   */
  std::vector<ADCcount> buildSamples(float targetGain) const;
  
  /// Returns the index of the sample under the peak of the waveform.
  static std::size_t findPeak(std::vector<ADCcount> const& samples);
  
}; // class icarus::opdet::SampledWaveformFunction<>


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename T>
icarus::opdet::SampledWaveformFunction<T>::SampledWaveformFunction
  (WaveformSpecs_t waveformSpecs, Time peakTime, float gain)
  : fSource    { std::move(waveformSpecs) }
  , fSamples   { buildSamples(gain) }
  , fGain      { gain }
  , fPeakSample{ findPeak(fSamples) }
  , fRefTime   { peakTime - fPeakSample * sampleDuration() }
  {}


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::SampledWaveformFunction<T>::doEvaluateAt(Time time) const
  -> ADCcount
{
  std::ptrdiff_t const iSample = static_cast<std::ptrdiff_t>
    (std::floor((time - fRefTime)/sampleDuration()));
  return hasSample(iSample)? fSamples[iSample]: ADCcount{ 0 };
} // icarus::opdet::SampledWaveformFunction<>::doEvaluateAt()


// -----------------------------------------------------------------------------
template <typename T>
void icarus::opdet::SampledWaveformFunction<T>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  
  out << firstIndent
      << "Pulse '" << fSource.name << "' (v. " << fSource.version
      << ", " << fSource.date << "):"
    << "\n" << indent << "  " << fSource.description
    << "\n" << indent
      << "  from " << fSamples.size() << "x " << sampleDuration() << " samples"
      << ", peak at " << Base_t::peakTime()
      << " with amplitude " << Base_t::peakAmplitude()
    << "\n" << indent
      << "  start at " << fRefTime << ", gain " << fGain
      << " (integral: " << integral() << ")"
    << '\n'
    ;
  
  // === BEGIN DEBUG FIXME DELME ===============================================
  out << "\n" << indent << "  samples:";
  for (std::size_t iSample = 0; iSample < fSamples.size(); ++iSample) {
    out << "\n" << indent << "   [" << iSample << "]  " << fSamples[iSample];
  }
  out << "\n";
  // === END DEBUG =============================================================
  
} // icarus::opdet::SampledWaveformFunction<>::doDump()


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::SampledWaveformFunction<T>::integral() const -> ADCcount
  { return std::reduce(fSamples.begin(), fSamples.end()); }


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::opdet::SampledWaveformFunction<T>::buildSamples
  (float targetGain) const -> std::vector<ADCcount>
{
  
  /*
   * Sample conversion
   * ------------------
   * 
   * The waveform is expected in millivolt, and it expresses the response with/
   * the photodetector set at a known gain.
   * Our target is a waveform in ADC and the target gain in argument.
   * The conversion factor is based on the full range of the digitizer,
   * that is 2 V across 14 bit.
   */
  constexpr float VoltageRange = 2'000.0; // millivolt
  constexpr unsigned short int ADCbits = 14;
  
  // 2 V in 14 bits (=> 8.192):
  constexpr float mVtoADC = (1 << ADCbits) / VoltageRange;
  
  // if either the starting gain is unknown or the target gain is not specified,
  // do not scale the gain
  float const gainFactor = (fSource.gain != 0.0 && targetGain != 0.0)
    ? (targetGain / fSource.gain): 1.0;
  float const factor = gainFactor * mVtoADC;
  
  auto voltageToADC = [factor](float mV)
    { return static_cast<ADCcount>(factor * mV); };
  
  std::vector<ADCcount> samples;
  samples.reserve(fSource.samples.size());
  std::transform(fSource.samples.begin(), fSource.samples.end(),
    back_inserter(samples), voltageToADC);
  
  return samples;
  
} // icarus::opdet::SampledWaveformFunction<>::buildSamples()


// -----------------------------------------------------------------------------
template <typename T>
std::size_t icarus::opdet::SampledWaveformFunction<T>::findPeak
  (std::vector<ADCcount> const& samples)
{
  assert(!samples.empty());
  auto const sbegin = samples.begin();
  auto const [ min, max ] = std::minmax_element(sbegin, samples.end());
  // assume baseline 0:
  return ((min->abs() > max->abs())? min: max) - sbegin;
} // icarus::opdet::SampledWaveformFunction<T>::findPeak()


// -----------------------------------------------------------------------------

#endif //  ICARUSCODE_PMT_ALGORITHMS_SAMPLEDWAVEFORMFUNCTION_H
