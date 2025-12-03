/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WaveformFilterAlg.cxx
 * @brief  Wrapper class for Fourier transforms and frequency-domain filtering.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 21, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/WaveformFilterAlg.h
 *
 * Link with -lfftw3f (single-precision FFTW).
 * 
 * This class was adapted from code generated via Copilot/GPT-5 mini.
 */

#undef NDEBUG // FIXME

// library header
#include "icaruscode/PMT/Trigger/Algorithms/WaveformFilterAlg.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"

// numeric libraries
#include <fftw3.h>

// C/C++ standard libraries
#include <cassert>
#include <cstring> // memcpy()
#include <utility> // std::move(), std::exchange(), std::swap()
#include <type_traits> // std::false_type...


// -----------------------------------------------------------------------------
using namespace util::quantities::frequency_literals;


// -----------------------------------------------------------------------------
namespace details {
  
  // some metaprogramming to support most contiguous containers
  // (that is: std::valarray and std::vector)
  
  template <typename T, typename As, typename = void>
  struct const_as { using type = std::remove_const_t<T>; };
  
  template <typename T, typename As>
  struct const_as<T, As, std::enable_if_t<std::is_const_v<As>>>
    { using type = std::add_const_t<T>; };
  
  template <typename Coll, typename = void>
  struct supports_data: std::false_type {};
  
  template <typename Coll>
  struct supports_data<Coll, decltype(std::data(std::declval<Coll>()))>
    : std::true_type {};
  
  template <typename T, typename As>
  using const_as_t = typename details::const_as<T, As>::type;
  
  template <typename Coll>
  constexpr bool supports_data_v = details::supports_data<Coll>::value;
  
} // namespace details

namespace {
  
  template <typename Coll>
  std::enable_if_t
    <details::supports_data_v<Coll>, decltype(std::data(std::declval<Coll>()))>
  dataOf(Coll& coll) { return std::data(coll); }
  
  template <typename Coll>
  std::enable_if_t<
    !details::supports_data_v<Coll>,
    details::const_as_t<typename Coll::value_type, Coll>*
    >
  dataOf(Coll& coll) { return &(coll[0]); }
  
} // local namespace


// -----------------------------------------------------------------------------
// ---  icarus::trigger::WaveformFilterAlg::InternalState
// -----------------------------------------------------------------------------
struct icarus::trigger::WaveformFilterAlg::InternalState {
  
  // --- BEGIN --- FFTW resource wrappers --------------------------------------
  /// @name FFTW resource wrappers
  /// @{

  /**
   * @brief Wrapper for `fftwf_plan`.
   *
   * This wrapper acquires an existing FFTW plan object (which turns out to be
   * an opaque pointer) and manages its lifetime so that its resources will be
   * "automatically" released when the object is destroyed.
   * 
   * ### Multithreading
   * 
   * The functionality of this class is thread-safe only in that it does very
   * little: it does not create a new plan, and it does not execute it either.
   * The destruction of the plan is thread-safe by itself, but it comes with
   * the usual caveat that no other thread should be using the plan after it
   * is destroyed.
   * 
   * However, for real thread-safety, the caller needs to use only thread-safe
   * calls to the plan returned by `get()`, which currently means to call only
   * `execute()`-like functions.
   */
  class Plan: private lar::UncopiableClass {
    
      private:
    
    fftwf_plan fPlan = nullptr;
    
    void destroy();
    
      public:

    Plan() = default;
    
    /// Acquires and manager a plan.
    explicit Plan(fftwf_plan&& p) noexcept;
    
    ~Plan() { destroy(); }
    
    Plan(Plan&& other) noexcept;
    Plan& operator= (Plan&& other) noexcept;
    
    /// Returns whether there is a valid plan.
    explicit operator bool() const noexcept { return fPlan != nullptr; }
    
    /// Returns the actual plan object (pointer).
    fftwf_plan get() const noexcept { return fPlan; }
    
    /// Stops managing the plan and returns it.
    fftwf_plan release() noexcept;

  }; // Plan

  /**
   * @brief FFTW data buffer.
   * @tparam T type of the data in the buffer
   *
   * A buffer has fixed size. Once constructed with a certain size, the size
   * can't be changed any more unless a new buffer is move-assigned into it.
   * 
   * ### Multithreading
   * 
   * This class is thread-safe:
   *  * calls of `const` member functions do not conflict with each other;
   *  * constructor and destructor use global state via `fftwf_malloc()` and
   *    `fftwf_free()`, which are declared to be also thread-safe.
   * 
   * As usual, apart from what stated in the previous list, non-const interface
   * does not offer thread-safety.
   * 
   */
  template <typename T>
  class Buffer: private lar::UncopiableClass {
    
    T* fData = nullptr;     ///< Pointer to the data.
    std::size_t fCount = 0; ///< Number of elements in the buffer.
    
    /// Destroys the resource.
    void destroy();
    
      public:
    
    /// Constructs an empty buffer.
    Buffer() noexcept = default;
    
    /// Allocates a new buffer of the specified size.
    explicit Buffer(std::size_t count);
    
    /// Releases all resources.
    ~Buffer() { destroy(); }
    
    Buffer(Buffer&& other) noexcept;
    Buffer& operator=(Buffer&& other) noexcept;
    
    /// Returns the data buffer.
    T* data() noexcept { return fData; }
    
    /// Returns the data buffer.
    T const* data() const noexcept { return fData; }
    
    /// Returns the size of the data buffer.
    std::size_t size() const noexcept { return fCount; }
    
    /// Returns whether any memory is allocated for the buffer.
    explicit operator bool() const noexcept { return fData != nullptr; }
    
    /// Releases and returns the buffer as is, yielding its lifetime management.
    T* release() noexcept;
    
  }; // Buffer

  /// @}
  // ---- END ---- FFTW resource wrappers --------------------------------------
  
  
  unsigned int FFTWflags; ///< Flags passed to `fftwf_plan_*`.
  
  
  /// Samples (real input buffer, size `planSize()`).
  Buffer<Sample_t> samples;
  
  /// Spectrum (complex output buffer, size `planSize()/2 + 1`).
  Buffer<fftwf_complex> spectrum;
  
  /// FFTW "plan" from real to complex, size `planSize()`.
  Plan samplesToSpectrumPlan;
  
  /// FFTW "plan" from complex to real, size `planSize()/2 + 1`.
  Plan spectrumToSamplesPlan;

  
  /// Constructor: stores the FFTW flags.
  InternalState(unsigned int flags): FFTWflags{ flags } {}
  
  /// Returns the size (in samples) of the currently allocated FFTW "plan".
  std::size_t planSize() const noexcept;
  
  /**
   * @brief Ensure internal buffers and FFTW plans are prepared for length N.
   * @param N number of samples in the waveform (must be > 0 and even).
   * @throws std::invalid_argument if `N` is null or odd.
   * @throws std::bad_alloc on allocation failure.
   * @throws std::runtime_error if plan creation fails.
   *
   * Allocates aligned buffers using `fftwf_malloc` and creates r2c/c2r plans
   * for `N`.
   */
  void ensureBuffersAndPlans(std::size_t N);

  /// Free any allocated FFTW plans and buffers and reset internal state.
  void freePlansAndBuffers();
  
}; // icarus::trigger::WaveformFilterAlg::InternalState


// -----------------------------------------------------------------------------
// ---  icarus::trigger::WaveformFilterAlg::InternalState::Plan
// -----------------------------------------------------------------------------
icarus::trigger::WaveformFilterAlg::InternalState::Plan::Plan
  (fftwf_plan&& p) noexcept
  : fPlan{ std::move(p) }
  {}


// -----------------------------------------------------------------------------
icarus::trigger::WaveformFilterAlg::InternalState::Plan::Plan
  (Plan&& other) noexcept
  : fPlan{ std::exchange(other.fPlan, nullptr) }
  {}


// -----------------------------------------------------------------------------
auto icarus::trigger::WaveformFilterAlg::InternalState::Plan::operator=
  (Plan&& other) noexcept -> Plan&
{
  if (this != &other) {
    destroy();
    std::swap(fPlan, other.fPlan);
  }
  return *this;
}


// -----------------------------------------------------------------------------
fftwf_plan icarus::trigger::WaveformFilterAlg::InternalState::Plan::release()
  noexcept
{
  return std::exchange(fPlan, nullptr);
}


// -----------------------------------------------------------------------------
void icarus::trigger::WaveformFilterAlg::InternalState::Plan::destroy() {
  if (fPlan) fftwf_destroy_plan(fPlan);
  fPlan = nullptr;
}


// -----------------------------------------------------------------------------
// ---  icarus::trigger::WaveformFilterAlg::InternalState::Buffer
// -----------------------------------------------------------------------------
template <typename T>
void icarus::trigger::WaveformFilterAlg::InternalState::Buffer<T>::destroy() {
  if (fData) {
    fftwf_free(fData);
    fData = nullptr;
  }
  fCount = 0;
}


// -----------------------------------------------------------------------------
template <typename T>
icarus::trigger::WaveformFilterAlg::InternalState::Buffer<T>::Buffer
  (std::size_t count)
{
  if (count == 0) return;
  void* p = fftwf_malloc(sizeof(T) * count);
  if (!p) throw std::bad_alloc{};
  fData = static_cast<T*>(p);
  fCount = count;
}


// -----------------------------------------------------------------------------
template <typename T>
icarus::trigger::WaveformFilterAlg::InternalState::Buffer<T>::Buffer
  (Buffer<T>&& other) noexcept
  : fData{ std::exchange(other.fData, nullptr) }
  , fCount{ std::exchange(other.fCount, 0) }
{}


// -----------------------------------------------------------------------------
template <typename T>
auto icarus::trigger::WaveformFilterAlg::InternalState::Buffer<T>::operator=
  (Buffer<T>&& other) noexcept -> Buffer<T>&
{
  if (this != &other) {
    destroy();
    std::swap(fData, other.fData);
    std::swap(fCount, other.fCount);
  }
  return *this;
}


// -----------------------------------------------------------------------------
template <typename T>
T* icarus::trigger::WaveformFilterAlg::InternalState::Buffer<T>::release()
  noexcept
{
  fCount = 0;
  return std::exchange(fData, nullptr);
}


// -----------------------------------------------------------------------------
std::size_t icarus::trigger::WaveformFilterAlg::InternalState::planSize()
  const noexcept
{
  // delegated to the buffer size, under the assumption that when there is a
  // plan there are also samples
  assert(bool(samples) == bool(samplesToSpectrumPlan));
  return samples.size();
}


// -----------------------------------------------------------------------------
void icarus::trigger::WaveformFilterAlg::InternalState::ensureBuffersAndPlans
  (std::size_t N)
{
  
  if (N == 0)
    throw std::invalid_argument{ "Fourier transform length N must be > 0" };
  
  // enforce even length assumption (runtime check)
  if ((N % 2) != 0)
    throw std::invalid_argument("WaveformFilterAlg requires even N");

  if (planSize() == N) return; // already prepared

  // free old resources
  freePlansAndBuffers();

  // allocate aligned buffers with fftwf_malloc via Buffer<T>
  std::size_t const K = spectrumSizeFromSamples(N);
  samples = Buffer<float>(N);
  spectrum = Buffer<fftwf_complex>(K);
  
  // create plans; if plan creation fails, clean up by letting destructors run
  fftwf_plan p_r2c
    = fftwf_plan_dft_r2c_1d(N, samples.data(), spectrum.data(), FFTWflags);
  if (!p_r2c) {
    freePlansAndBuffers();
    throw std::runtime_error{ "WaveformFilterAlg: fftwf_plan_dft_r2c_1d() failed" };
  }
  samplesToSpectrumPlan = Plan(std::move(p_r2c));
  
  fftwf_plan p_c2r
    = fftwf_plan_dft_c2r_1d(N, spectrum.data(), samples.data(), FFTWflags);
  if (!p_c2r) {
    freePlansAndBuffers();
    throw std::runtime_error{ "WaveformFilterAlg: fftwf_plan_dft_c2r_1d() failed" };
  }
  spectrumToSamplesPlan = Plan(std::move(p_c2r));
  
} // icarus::trigger::WaveformFilterAlg::ensureBuffersAndPlans()


// -----------------------------------------------------------------------------
void icarus::trigger::WaveformFilterAlg::InternalState::freePlansAndBuffers() {
  // reset plans and buffers
  samplesToSpectrumPlan = Plan();
  spectrumToSamplesPlan = Plan();
  samples = Buffer<float>();
  spectrum = Buffer<fftwf_complex>();
} // icarus::trigger::WaveformFilterAlg::freePlansAndBuffers()


// -----------------------------------------------------------------------------
// ---  icarus::trigger::WaveformFilterAlg
// -----------------------------------------------------------------------------
icarus::trigger::WaveformFilterAlg::WaveformFilterAlg
  (util::quantities::hertz sampleRate)
  : WaveformFilterAlg{ sampleRate, FFTW_MEASURE }
  {}


// -----------------------------------------------------------------------------
icarus::trigger::WaveformFilterAlg::WaveformFilterAlg
  (util::quantities::hertz sampleRate, unsigned int FFTWflags)
  : fSampleRate{ sampleRate }
  , fState{ std::make_unique<InternalState>(FFTWflags) }
{
  if (fSampleRate <= 0.0_Hz)
    throw std::invalid_argument("Sample rate must be positive.");
}


// -----------------------------------------------------------------------------
icarus::trigger::WaveformFilterAlg::~WaveformFilterAlg() = default;


// -----------------------------------------------------------------------------
void icarus::trigger::WaveformFilterAlg::applyInPlace(
  std::valarray<Sample_t>& samples,
  std::function<std::complex<Sample_t>(Sample_t)> const& filterFunc
) {
  assert(fState);
  
  std::size_t const N = samples.size();
  
  if (N == 0) throw std::invalid_argument{ "Size of data must be non-null." };
  
  if ((N % 2) != 0) {
    throw std::invalid_argument{
      "WaveformFilterAlg::applyInPlace() requires even number of input samples."
      };
  }
  
  fState->ensureBuffersAndPlans(N);
  
  Sample_t* in = fState->samples.data();
  
  // fast copy of contiguous valarray data -> aligned buffer
  std::memcpy(in, dataOf(samples), sizeof(Sample_t) * N);

  // execute forward transform
  // TODO this can be called in parallel, but we need the version with its own I/O buffers for that to be useful
  fftwf_execute(fState->samplesToSpectrumPlan.get());

  fftwf_complex* out = fState->spectrum.data();
  std::size_t const K = spectrumSizeFromSamples(N);

  // apply filter function per bin using in-place arithmetic
  // (avoids constructing std::complex)
  Sample_t const baseFreq
    = fSampleRate.convertInto<util::quantities::hertz>().value() / N;
  for (std::size_t const k: util::counter(K)) {
    
    std::complex<Sample_t> const f = filterFunc(k * baseFreq);
    Sample_t const fr = f.real();
    Sample_t const fi = f.imag();

    Sample_t const ar = out[k][0];
    Sample_t const ai = out[k][1];

    // (x + i y) * (a + i b) = (x*a - y*b) + i(x*b + y*a)
    out[k][0] = ar * fr - ai * fi;
    out[k][1] = ar * fi + ai * fr;
    
  } // for

  // inverse transform (produces unnormalized time-domain signal in `fSamples`)
  fftwf_execute(fState->spectrumToSamplesPlan.get());

  // normalize by dividing by N
  float const n = 1.0f / N;
  for (std::size_t i: util::counter(N)) samples[i] = in[i] * n;
  
  return;
  
  /* FIXME
  std::valarray<std::complex<float>> spectrum = waveformToSpectrum(samples);
  
  Sample_t const baseFreq
    = fSampleRate.convertInto<util::quantities::hertz>().value() / N;
  for (std::size_t const k: util::counter(K))
    spectrum[k] *= filterFunc(k * baseFreq);
  
  samples = spectrumToWaveform(spectrum);
  */
  
} // icarus::trigger::WaveformFilterAlg::applyInPlace()


// -----------------------------------------------------------------------------
auto icarus::trigger::WaveformFilterAlg::waveformToSpectrum
  (std::valarray<Sample_t> const& samples) -> std::valarray<std::complex<Sample_t>>
{
  assert(fState);
  
  std::size_t const N = samples.size();
  
  if (N == 0) throw std::invalid_argument{ "Size of data must be non-null." };
  
  if ((N % 2) != 0) {
    throw std::invalid_argument{
      "WaveformFilterAlg::waveformToSpectrum() requires even number of input samples."
      };
  }
  
  fState->ensureBuffersAndPlans(N);
  
  Sample_t* in = fState->samples.data();

  // fast copy of contiguous valarray data -> aligned buffer
  std::memcpy(in, dataOf(samples), sizeof(Sample_t) * N);

  // execute forward transform
  fftwf_execute(fState->samplesToSpectrumPlan.get());

  fftwf_complex* out = fState->spectrum.data();
  std::size_t const K = spectrumSizeFromSamples(N);

  std::valarray<std::complex<float>> spectrum(K);
  
  // we get our chance to normalize the spectrum amplitudes
  for (std::size_t const k: util::counter(K))
    spectrum[k] = std::complex{ out[k][0], out[k][1] } / Sample_t(N);

  return spectrum;
} // icarus::trigger::WaveformFilterAlg::waveformToSpectrum()


// -----------------------------------------------------------------------------
auto icarus::trigger::WaveformFilterAlg::spectrumToWaveform
  (std::valarray<std::complex<Sample_t>> const& spectrum) -> std::valarray<Sample_t>
{
  assert(fState);
  
  // The inverse transform is unnormalized in FFTW;
  // the implementation divides by N after inverse.
  std::size_t const K = spectrum.size();

  if (K < 2) throw std::invalid_argument{ "Spectrum must not be empty." };

  // Canonical even-length choice: N = 2*(K-1)
  std::size_t const N = samplesFromSpectrumSize(K);

  fState->ensureBuffersAndPlans(N);
  
  fftwf_complex* out = fState->spectrum.data();
  for (std::size_t const k: util::counter(K)) {
    out[k][0] = spectrum[k].real();
    out[k][1] = spectrum[k].imag();
  }

  fftwf_execute(fState->spectrumToSamplesPlan.get());

  std::valarray<Sample_t> samples(N);
  Sample_t* in = fState->samples.data();
  // float const n = 1.0f / N;
  for (std::size_t const i: util::counter(N)) samples[i] = in[i]; // * n;
  
  return samples;
} // icarus::trigger::WaveformFilterAlg::spectrumToWaveform()


// -----------------------------------------------------------------------------
auto icarus::trigger::WaveformFilterAlg::filterSpectrum(
  std::valarray<Sample_t> const& frequencies,
  std::valarray<std::complex<Sample_t>> spectrum,
  std::function<std::complex<Sample_t>(Sample_t)> const& filterFunc
) const -> std::valarray<std::complex<Sample_t>> {
  
  assert(size(frequencies) == size(spectrum));
  
  auto itFreq = begin(frequencies);
  for (std::complex<Sample_t>& amplitude: spectrum) {
    amplitude *= filterFunc(*itFreq);
    ++itFreq;
  }
  assert(itFreq == end(frequencies));
  
  return spectrum;
} // icarus::trigger::WaveformFilterAlg::filterSpectrum(frequencies)


// -----------------------------------------------------------------------------
auto icarus::trigger::WaveformFilterAlg::filterSpectrum(
  std::valarray<std::complex<Sample_t>> spectrum,
  std::function<std::complex<Sample_t>(Sample_t)> const& filterFunc
) const -> std::valarray<std::complex<Sample_t>> {
  
  return filterSpectrum(
    spectrumFrequencies(samplesFromSpectrumSize(spectrum.size())),
    std::move(spectrum), filterFunc
  );
  
} // icarus::trigger::WaveformFilterAlg::filterSpectrum()


// -----------------------------------------------------------------------------
auto icarus::trigger::WaveformFilterAlg::spectrumFrequencies
  (std::size_t nSamples) const -> std::valarray<Sample_t>
{
  
  if (nSamples == 0) throw std::invalid_argument{ "Size of data must be non-null." };
  
  if ((nSamples % 2) != 0) {
    throw std::invalid_argument{
      "WaveformFilterAlg::spectrumFrequencies() requires even number of input samples."
      };
  }
  
  std::size_t const K = spectrumSizeFromSamples(nSamples);
  
  auto const baseFreq = static_cast<Sample_t>
    (samplingRate().convertInto<util::quantities::hertz>().value() / nSamples);
  
  std::valarray<Sample_t> freq(K);
  for (std::size_t k = 0; k < K; ++k) freq[k] = k * baseFreq;
  
  return freq;
  
} // icarus::trigger::WaveformFilterAlg::spectrumFrequencies()


// -----------------------------------------------------------------------------
