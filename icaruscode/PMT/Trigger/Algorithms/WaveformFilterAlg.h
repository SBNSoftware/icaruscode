/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WaveformFilterAlg.h
 * @brief  Wrapper class for Fourier transforms and frequency-domain filtering.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 21, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/WaveformFilterAlg.cxx
 *
 * This class was adapted from code generated via Copilot/GPT-5 mini.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WAVEFORMFILTERALG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WAVEFORMFILTERALG_H

// LArSoft libraries
#include "lardataalg/Utilities/quantities/frequency.h" // util::quantities::hertz
#include "lardataalg/Utilities/quantities/spacetime.h" // util::quantities::second
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

// C/C++ standard libraries
#include <complex>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <valarray>


// -----------------------------------------------------------------------------
namespace icarus::trigger { class WaveformFilterAlg; }

/**
 * @brief Handles waveform transformation between time and frequency domain.
 * 
 * `WaveformFilterAlg` provides a wrapper to:
 *  * convert a real-valued waveform to the frequency domain using Fourier
 *    transform (currently, FFTW single-precision API);
 *  * apply an arbitrary complex-valued filter per frequency bin;
 *  * convert the filtered spectrum back to the time domain.
 *
 * This algorithm requires the input time-domain length to be even.
 *
 * Usage example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarus::trigger::WaveformFilterAlg filter{ 500_000_000.0 };
 * 
 * std::valarray<float> waveform = ...; // arbitrary even length > 0
 * 
 * auto Cutoff10MHz = [](float freq_hz) -> std::complex<float>
 *   { return (freq_hz <= 10e6f)? 1.0f: 0.0f; };
 * 
 * std::valarray<float> filtered = filter.applyFilter(waveform, Cutoff10MHz);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Multithreading
 * ---------------
 * 
 * No.
 * 
 * Not yet, at least. This class can make and overwrite buffers at leisure, and
 * plans too. To make it speak multithreading, it should be factored so that
 * 1. FFTW "plans" are cached by size, and the cache object be thread-safe
 * 2. "executor" objects can be spawned, and each should manage their own
 *    buffer; the thread-safety of an executor can be that multiple different
 *    executors can concurrently run, but an executor object may still require
 *    to be executed only on a single thread.
 *
 */
class icarus::trigger::WaveformFilterAlg: private lar::UncopiableClass {
  
    public:
  
  using Sample_t = float; ///< Data type of the samples and harmonic amplitudes.
  
  /**
   * @brief Construct a filter wrapper with the given sample rate.
   * @param sampleRate sample rate
   * @throws std::invalid_argument if `sampleRate` is invalid (non-positive)
   * 
   * A default set of flags is used (i.e. `FFTW_MEASURE`).
   */
  explicit WaveformFilterAlg(util::quantities::hertz sampleRate);

  /**
   * @brief Construct a filter wrapper with the given sample rate.
   * @param sampleRate sample rate
   * @param FFTWflags flags passed to `fftwf_plan_*`
   * @throws std::invalid_argument if `sampleRate` is invalid (non-positive)
   */
  explicit WaveformFilterAlg
    (util::quantities::hertz sampleRate, unsigned int FFTWflags);

  
  ~WaveformFilterAlg(); // need to define after InternalState is defined
  
  
  // Movable
  WaveformFilterAlg(WaveformFilterAlg&&) noexcept = default;
  WaveformFilterAlg& operator=(WaveformFilterAlg&&) noexcept = default;

  /**
   * @brief Apply a frequency-domain filter in-place on `data` (time-domain).
   * @param[in,out] samples time-domain signal; size must be even
   * @param filterFunc function returning the amplification of the frequency in
   *                   argument (in hertz)
   * @throws std::invalid_argument if `data` is empty or of odd size
   * 
   * This method applies a filter function to a discrete (sampled) waveform.
   * 
   * It transforms `samples` (length _N_) into frequency domain, discrete in
   * _N_ / 2 frequency bins, and multiplies the amplitude of each frequency bin
   * by the filter function value for that frequency; and finally transforms
   * the result back into frequency domain.
   * 
   */
  void applyInPlace(
    std::valarray<Sample_t>& samples,
    std::function<std::complex<Sample_t>(Sample_t)> const& filterFunc
    );

  /**
   * @brief Returns the complex spectrum of a discrete waveform.
   * @param samples time-domain discrete waveform (even number _N_ of samples)
   * @return the spectrum as _N_ / 2 + 1 harmonic complex coefficients
   * @throws std::invalid_argument if `samples` length is odd or null
   *
   * This algorithm uses FFTW real-to-complex ("r2c") transform.
   *
   * The spectrum is returned for discrete values of frequency: the frequency
   * for bin _k_ is _k_ / _N_ * `samplingRate()`.
   * The first term ("frequency 0") is a constant term and real. The spectrum
   * includes harmonics from `samplingRate()` / _N_ to `samplingRate()` (for
   * example, a waveform of 10 &micro;s sampled at 500 MHz, with
   * _N_ = 5000 samples, yield a spectrum with harmonics from 100 kHz up to 500
   * MHz).
   */
  std::valarray<std::complex<Sample_t>> waveformToSpectrum
    (std::valarray<Sample_t> const& samples);

  /**
   * @brief Returns the discrete time-domain representation of a spectrum.
   * @param spectrum discrete complex spectrum (length _K_)
   * @return time-domain signal
   * @throws std::invalid_argument if `spectrum` is empty
   * @see `spectrumFrequencies()`
   *
   * A time-domain waveforms with the configured sampling rate is
   * "reconstituted" from its harmonic components specified in `spectrum`.
   * 
   * The spectrum is expected to be represented as complex amplitudes for
   * harmonics of discrete frequency: the amplitude of bin _k_ applies to
   * the frequency _k_ / _N_ * `samplingRate()`.
   * 
   * The length _N_ of the returned waveform is deduced as the canonical even
   * time-domain length from the spectrum length _K_: _N_ = 2*(_K_-1) .
   *
   * This function validates that the spectrum is non-empty and will throw
   * std::invalid_argument if spectrum.empty().
   * 
   * This algorithm uses FFTW complex-to-real ("c2r") transform.
   *
   */
  std::valarray<Sample_t> spectrumToWaveform
    (std::valarray<std::complex<Sample_t>> const& spectrum);
  
  /**
   * @brief Applies a filter to the specified spectrum.
   * @param frequencies the frequencies for the spectrum elements [Hz]
   * @param spectrum amplitudes of the spectrum
   * @param filterFunc a filtering function parametrized on frequency
   * @return the filtered spectrum (same format as `spectrum`)
   * 
   * The filter is applies on `spectrum` and a new, filtered spectrum is
   * returned.
   * 
   * The `spectrum` is expected to have one amplitude corresponding on each
   * frequency in the `frequencies` array.
   * The array of `frequencies` can be obtained from `spectrumFrequencies()`.
   * The filter function object is a functor that takes as argument a frequency
   * (double precision) in hertz, and returns a "dump" factor to be multiplied
   * to the complex amplitude associated to that frequency.
   * 
   * To perform the change "in place", the following approach is supported:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * spectrum = filterAlg.filterSpectrum(frequencies, std::move(spectrum), filter);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which does not involve any copy of `spectrum`.
   * 
   * @note The "same" dump factor is applied to both "positive" and "negative"
   *       frequencies. More precisely, if factor _d_ is applied to frequency
   *       _f_, frequency _-f_ is implicitly applied the complex conjugate of
   *       _d_, so that the resulting spectrum is still complex conjugate and
   *       the transformation back to time domain is still real.
   *       Here in lieu of a frequency we have a phase _i &omega;_ (_i_ being
   *       the imaginary unit) that figures in the exponent of the components
   *       of the waveform, @f$ \sum_{k} a_{k} e^{i\omega_{k}t} @f$, and the
   *       sign of the frequency implied the sign of the exponent.
   * 
   */
  std::valarray<std::complex<Sample_t>> filterSpectrum(
    std::valarray<Sample_t> const& frequencies,
    std::valarray<std::complex<Sample_t>> spectrum,
    std::function<std::complex<Sample_t>(Sample_t)> const& filterFunc
    ) const;
  
  /**
   * @brief Applies a filter to the specified spectrum.
   * @param spectrum amplitudes of the spectrum
   * @param filterFunc a filtering function parametrized on frequency
   * @return the filtered spectrum (same format as `spectrum`)
   * 
   * The filter is applies in place on `spectrum`.
   * 
   * The `spectrum` is expected to have one amplitude corresponding on each
   * frequency in a list of frequencies, which is obtained from
   * `spectrumFrequencies()`.
   * 
   * See the other signature of `filterSpectrum()` for more details.
   */
  std::valarray<std::complex<Sample_t>> filterSpectrum(
    std::valarray<std::complex<Sample_t>> spectrum,
    std::function<std::complex<Sample_t>(Sample_t)> const& filterFunc
    ) const;
  
  /**
   * @brief Returns the frequencies of the spectrum of a `n`-sample waveform.
   * @param nSamples number of samples in the waveform
   * @return the sequence of (`nSamples` / 2 + 1) frequencies, in hertz
   * @see `waveformToSpectrum()`
   * 
   * Given a sequence of `nSamples` samples, the relevant frequencies _f_ in the
   * spectrum analysis are from `0` (continuum) to `samplingRate()`/2 included.
   * in steps of `samplingRate() / nSamples` (inverse of the waveform duration).
   * 
   * This function returns the sequence of the frequency values, in the same
   * order as they apply to the result of `waveformToSpectrum()` function.
   */
  std::valarray<Sample_t> spectrumFrequencies(std::size_t nSamples) const;
  
  
  /// Returns the configured sampling rate.
  util::quantities::hertz samplingRate() const noexcept { return fSampleRate; }

  /// Returns the configured sampling period.
  util::quantities::second samplingPeriod() const noexcept
    { return 1/fSampleRate; }

    private:

  class InternalState; // implementation delegation (pImpl pattern)
  
  // --- BEGIN --- Algorithm configuration -----------------------------------
  
  util::quantities::hertz fSampleRate; ///< Sample rate [Hz]
  
  // ---- END ---- Algorithm configuration -----------------------------------

  std::unique_ptr<InternalState> fState; ///< The internal state.
  
  /// Returns the number of spectrum amplitudes  (including the constant term)
  /// for a waveform with `nSamples` samples.
  static std::size_t spectrumSizeFromSamples(std::size_t nSamples)
    { return nSamples / 2 + 1; }
  
  /// Returns the number of samples of a waveform with a spectrum of
  /// `nAmplitudes` coefficients (including the constant term).
  static std::size_t samplesFromSpectrumSize(std::size_t nAmplitudes)
    { return (nAmplitudes - 1) * 2; }
  
}; // icarus::trigger::WaveformFilterAlg


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WAVEFORMFILTERALG_H
