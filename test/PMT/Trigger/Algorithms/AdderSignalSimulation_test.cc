/**
 * @file   AdderSignalSimulation_test.cc
 * @brief  Test for `icarus::trigger::AdderChannelSimulator` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 7, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.h
 * 
 */


// library being tested
#include "icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.h"
#include "icaruscode/PMT/Trigger/Algorithms/WaveformFilterAlg.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond...
#include "lardataalg/Utilities/quantities/frequency.h" // hertz
#include "lardataalg/Utilities/StatCollector.h" // lar::util::MinMaxCollector
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// C/C++ standard libraries
#include <cstdint> // std::size_t
#include <array>
#include <fstream>
#include <type_traits> // std::remove_cv_t, ...
#include <utility> // std::move()
#include <valarray>

// Boost libraries
#define BOOST_TEST_MODULE ( TriggerGateData_test )
#include <boost/test/unit_test.hpp>


// -----------------------------------------------------------------------------
// verbatim from https://en.cppreference.com/w/cpp/container/array/to_array
// remove when C++20 is available
namespace detail {
  template<class T, std::size_t N, std::size_t... I>
  constexpr std::array<std::remove_cv_t<T>, N>
  to_array_impl(T (&&a)[N], std::index_sequence<I...>)
    { return {{std::move(a[I])...}}; }
}

template<class T, std::size_t N>
constexpr std::array<std::remove_cv_t<T>, N> to_array(T (&&a)[N])
  { return detail::to_array_impl(std::move(a), std::make_index_sequence<N>{}); }


// -----------------------------------------------------------------------------
// --- global definitions
// ---
using util::quantities::second, util::quantities::microsecond, util::quantities::nanosecond;
using util::quantities::hertz, util::quantities::megahertz;
using namespace util::quantities::time_literals;
using namespace util::quantities::frequency_literals;


// -----------------------------------------------------------------------------
// ---  lower level utilities
// ---

/// Returns the value of `t` in seconds.
template <typename S, typename T>
T toSecond(util::quantities::scaled_second<S, T> t)
  { return t.template convertInto<util::quantities::second_as<T>>().value(); }

/// Use `out << dumpCollection(v)` to print all elements in one line...
template <typename Coll>
struct dumpCollection {
  static constexpr auto OhSoMany = std::numeric_limits<std::size_t>::max();
  Coll const* coll = nullptr;
  std::size_t perLine = OhSoMany;
  std::size_t maxElement = OhSoMany;
  
  dumpCollection(
    Coll const& coll, std::size_t perLine = 0, std::size_t maxElement = OhSoMany
    )
    : coll{ &coll }, perLine{ perLine }, maxElement{ maxElement } {}
  
  template <typename Stream>
  void operator() (Stream& out) const {
    using std::size, std::begin, std::end;
    out << "(" << size(*coll) << " elements) [";
    std::size_t leftInLine = 0;
    for (auto const& [ i, e ]: util::enumerate(*coll)) {
      if ((perLine > 0) && (leftInLine-- == 0)) {
        leftInLine = perLine;
        out << "\n  [" << i << "] ";
      }
      if (i > maxElement) {
        out << " ... and " << (size(*coll) - i) << " more";
        break;
      }
      out << " " << e;
    }
    out << " ]";
  }
}; // struct dumpCollection

template <typename Stream, typename Coll>
Stream& operator<< (Stream&& out, dumpCollection<Coll> const& collDumper) {
  collDumper(out);
  return out;
}


// -----------------------------------------------------------------------------
// ---  waveform managing utilities
// ---

// Description of one component of the waveform.
template <typename T = float>
struct Harmonic_t {
  using value_type = T;    // type of waveform level/voltage/amplitude/samples
  
  value_type amplitude;    // amplitude of this harmonic (must be positive)
  double periods;          // number of periods in the interval
  value_type baseline = 0;
  double phase = 0;        // relative phase in sine wave [rad]
};


/// Returns a sequence of `nSamples` times in seconds, spaced by `samplingPeriod`.
template <typename T>
std::valarray<T> makeTimes(second samplingPeriod, std::size_t nSamples) {
  std::valarray<T> ts(nSamples);
  std::iota(begin(ts), end(ts), 0.0);
  return ts * samplingPeriod.value(); // argument is deliberately in seconds
} // makeTimes()


/**
 * Adjusts duration to be a multiple of `samplingPeriod`.
 * @tparam Time a time duration type (plain number or Quantity)
 * @returns the number of samples
 */
template <typename Time, typename Time2>
std::size_t adjustTimes(Time& duration, Time2 const& samplingPeriod) {
  std::size_t const nSamples = std::round(duration / samplingPeriod);
  duration = nSamples * samplingPeriod; // no time for rounding errors
  return nSamples;
}


/// Returns a collection of the frequencies of `harmonics` (same order).
template <typename Coll, typename Time>
std::vector<megahertz> harmonicFreqs(Coll const& harmonics, Time duration) {
  std::vector<megahertz> freqs(size(harmonics));
  for (std::size_t i = 0; i < freqs.size(); ++i)
    freqs[i] = harmonics[i].periods / duration;
  return freqs;
} // harmonicFreqs()


/// Returns the baseline of the waveform from the specified `harmonics`.
template <typename Coll>
auto harmonicsBaseline(Coll const& harmonics) {
  using value_t = typename Coll::value_type::value_type;
  value_t baseline = 0;
  for (Harmonic_t<value_t> const& harmonic: harmonics)
    baseline += harmonic.baseline;
  return baseline;
} // harmonicsBaseline()


template <typename Components>
auto composeWaveform(Components const& components) {
  std::valarray<typename Components::value_type::value_type> samples
    (0.0, size(components.front()));
  for (auto const& component: components) samples += component;
  return samples;
}


template <typename Coll>
std::valarray<typename Coll::value_type::value_type> generateWaveform(
  std::size_t nSamples, nanosecond samplingPeriod, Coll const& harmonics
) {
  
  using value_t = typename Coll::value_type::value_type;
  using Samples_t = std::valarray<value_t>; // type of waveform samples
  
  microsecond const duration = nSamples * samplingPeriod;
  
  std::vector<megahertz> const freqs = harmonicFreqs(harmonics, duration);
  
  std::cout << "Test filtering a " << duration << " waveform ("
    << nSamples << " samples, " << samplingPeriod << " each) with "
    << harmonics.size() << " harmonics:";
  for (auto const& [ i, harmonic ]: util::enumerate(harmonics)) {
    std::cout << "\n [#" << i << "] A=" << harmonic.amplitude
      << " f=" << freqs[i];
    if (harmonic.phase != 0)
      std::cout << " φ=" << (harmonic.phase / util::pi()) << "π rad";
  }
  std::cout << std::endl;
  
  // time corresponding to each tick; a bit of a waste, but handy.
  std::vector<Samples_t> ComponentSamples;
  std::valarray const times_s = makeTimes<value_t>(samplingPeriod, nSamples);
  for (auto const [ i, harmonic ]: util::enumerate(harmonics)) {
    auto const phases
      = (2 * util::pi() * hertz{ freqs[i] }.value()) * times_s + harmonic.phase;
    ComponentSamples.push_back(harmonic.amplitude * std::cos(phases));
  } // for
  
  value_t const totalBaseline = harmonicsBaseline(harmonics);
  
  return composeWaveform(ComponentSamples) + totalBaseline;

} // generateWaveform()


template <typename T>
std::valarray<T> waveformFromSpectrum(
  std::valarray<T> const& frequencies_hz,
  std::valarray<std::complex<T>> const& spectrum,
  nanosecond samplingPeriod, std::size_t nSamples
) {
  using value_t = T;
  
  // by hand
  std::valarray<T> samples(0.0, nSamples);
  std::valarray const times_s = makeTimes<value_t>(samplingPeriod, nSamples);
  for (auto const [ freq, ampl ]: util::zip(frequencies_hz, spectrum)) {
    
    if (freq == 0) samples += ampl.real(); // constant term comes out already doubled
    else {
      // |ampl| cos (phase + arg(ampl)) = Re(ampl) cos(phase) - Im(ampl) sin(phase)
      // the factor 2 is due to the loop over only half the complex amplitudes
      // (the others, complex conjugate of these ones, are customarily omitted)
      auto const phases = (2 * util::pi() * freq) * times_s + std::arg(ampl);
      samples += (2 * std::abs(ampl)) * std::cos(phases);
    }
  } // for
  
  return samples;
}


// -----------------------------------------------------------------------------
// ---  filter managing utilities
// ---
// -----------------------------------------------------------------------------
template <typename Coll, typename Filter>
Coll filterHarmonics
  (Coll const& harmonics, microsecond duration, Filter&& filter)
{
  using value_t = typename Coll::value_type::value_type;
  Coll filteredHarmonics = harmonics;
  for (Harmonic_t<value_t>& harmonic: filteredHarmonics) {
    hertz const freq = harmonic.periods / duration;
    value_t const dump { filter(freq) };
    harmonic.amplitude *= dump;
  }
  return filteredHarmonics;
} // filterHarmonics()


template <typename T>
struct ThreeLevelHighPassFrequencyFilter {
  
  using value_t = T;
  
  hertz fStartHalve = 0_Hz; ///< Frequencies above this one are halved.
  hertz fStartFull = 0_Hz; ///< Frequencies above this one are untouched.
  
  double factor(hertz freq) const {
    if (freq == 0_Hz) return 1.0; // preserve the baseline
    if (freq >= fStartFull) return 1.0;
    if (freq >= fStartHalve) return 0.5;
    return 0;
    
  }
  
  value_t operator() (hertz freq) const
    { return static_cast<value_t>(factor(freq)); }
  
}; // ThreeLevelHighPassFrequencyFilter


template <typename Coll, typename Time>
ThreeLevelHighPassFrequencyFilter<typename Coll::value_type::value_type>
makeThreeLevelHighPassFrequencyFilter(Coll const& harmonics, Time duration) {
  
  using value_t = typename Coll::value_type::value_type;
  using Filter_t = ThreeLevelHighPassFrequencyFilter<value_t>;
  
  // decide which are the limiting frequencies
  std::size_t const n = size(harmonics);
  
  assert(n > 0);
  
  std::vector<megahertz> const freqs = harmonicFreqs(harmonics, duration);
  
  /*
   * dumped by 50% in a frequency range +/- 50% around the only available frequency
   * 
   *                                     =================>
   * 
   *               ======================
   * 
   *  +============-----------+-----------------------------
   *  0                    freqs[0]
   * 
   */
  if (n == 1) {
    return Filter_t{ freqs[0] * 0.5, freqs[0] * 1.5 };
  }
  
  /*
   * no dump for frequency from the middle between the two frequencies on,
   * and 50% dump from half the first frequency up to the point above
   * 
   *                           ===========================>
   * 
   *          =================
   * 
   *  +=======-------+-----------------+--------------------
   *  0           freqs[0]          freqs[1]
   */
  if (n == 2) {
    return Filter_t{ freqs[0] * 0.5, (freqs[0] + freqs[1]) * 0.5 };
  }
  
  /*
   * split the list of frequencies in three parts with roughly the same number
   * of frequencies in each of them, and use the middle frequency between the
   * borders of those groups as delimiters
   * (in grouping, higher groups are extended first).
   * (illustration with n = 3)
   *                                           ===========>
   * 
   *                           ================
   * 
   *       <=========+=========--------+-------------+-------
   *              freqs[0]          freqs[1]      freqs[2]
   */
  
  return Filter_t{
    (freqs[n/3 - 1] + freqs[n/3]) * 0.5,
    (freqs[2*n/3 - 1] + freqs[2*n/3]) * 0.5,
  };
  
} // makeThreeLevelHighPassFrequencyFilter()


// -----------------------------------------------------------------------------
// ---  testing support
// ---

/// Equality operator: `true` if `a` and `b` are closer than `maxDiff`.
struct EqualWithin {
  double maxDiff;
  template <typename T, typename U>
  bool operator() (T a, U b) const noexcept
    { return std::abs(static_cast<double>(a - b)) <= maxDiff; }
};
template <typename T>
struct EqualTo {
  T ref;
  EqualWithin cmp;
  template <typename U>
  bool operator() (U v) const noexcept { return cmp(static_cast<T>(v), ref); }
};
template<class T> EqualTo(T, EqualWithin)-> EqualTo<T>;


template <typename T, typename U>
T roundAt(T v, U at) {
  return at * static_cast<T>(std::round(v / at));
}

template <typename T, typename U>
std::complex<T> roundAt(std::complex<T> v, U at) {
  return { roundAt(v.real(), at), roundAt(v.imag(), at) };
}

template <typename T>
constexpr T radDiff(T a, T b) {
  T d = a - b;
  while (d >= util::pi<T>()) d -= 2*util::pi<T>();
  while (d < -util::pi<T>()) d += 2*util::pi<T>();
  return d;
}


template <typename Spectrum, typename Coll>
void checkSpectrum(
  Coll const& harmonics,
  microsecond duration,
  std::valarray<typename Coll::value_type::value_type> const& frequencies,
  Spectrum const& spectrum
) {
  /*
   * Compares a spectrum with the harmonics that generated it.
   * 
   * The spectrum is a sequence of complex amplitudes, one per frequency.
   * Checks include the size of the spectrum, the amplitude and phase of each
   * of the generated harmonics, and the absence of all the harmonics
   * that were not generated.
   * Console output details successful and failing checks.
   */
  
  using value_t = typename Coll::value_type::value_type;
  
  BOOST_TEST_REQUIRE(frequencies.size() > 1);
  std::size_t const nSamples = (frequencies.size() - 1) * 2;
  nanosecond const samplingPeriod = duration / nSamples;
  
  std::cout << "Spectrum for " << nSamples << "-sample, " << duration
    << "-long waveform (sampling time: " << samplingPeriod << ")" << std::endl;
  
  BOOST_TEST(spectrum.size() == nSamples / 2 + 1);
  
  // std::cout << "Extracted spectrum with " << spectrum.size() << " frequencies: "
  //   << dumpCollection(spectrum, 16) << std::endl; // FIXME
  
  std::vector<megahertz> const testFreqs = harmonicFreqs(harmonics, duration);
  value_t const testBaseline = harmonicsBaseline(harmonics);
  
  // compute the smallest requested amplitude and base comparison tolerance on it
  value_t minAmplitude = harmonics[0].amplitude;
  for (std::size_t i = 0; i < testFreqs.size(); ++i) {
    if (harmonics[i].amplitude < minAmplitude)
      minAmplitude = harmonics[i].amplitude;
  }
  value_t const amplTol = minAmplitude / 1e3; // 0.1% of the smallest amplitude
  
  megahertz const baseFreq = 1.0 / duration; // e.g. 10 us -> 100 kHz.
  
  unsigned int nFoundFreqs = 0, nUnexpectedFrequencies = 0;
  auto const fbegin = cbegin(testFreqs), fend = cend(testFreqs);
  std::cout << "Detected baseline: " << roundAt(spectrum[0], 1e-4)
    << " (expected: " << testBaseline << ")" << std::endl;
  
  EqualWithin const freqMatch{ baseFreq.convertInto<hertz>().value()/2 };
  auto const rnd = [](auto v){ return roundAt(v, 1e-4); };
      
  for (auto const [ freq, amplitude ]: util::zip(frequencies, spectrum)) {
    if (freq == 0) { // special case: constant term (f=0, baseline)
      BOOST_TEST(amplitude.real() == testBaseline, 0.01% boost::test_tools::tolerance());
      BOOST_TEST(amplitude.imag() == 0, 0.01% boost::test_tools::tolerance());
      continue;
    }
    
    // since the input is N real samples, the N Fourier coefficients are
    // complex conjugate, and the result stores only one of them ("A");
    // so the actual sine/cosine components, cosine's (A + A*) and sine's
    // -i(A - A*), have magnitude twice as the real/imaginary components of A.
    double const modulus = std::abs(amplitude);
    double const phase = std::arg(amplitude);
    // tolerance is relative half the step 
    auto const itFreq
      = std::find_if(fbegin, fend, EqualTo{ hertz{ freq }, freqMatch });
    if (itFreq == fend) { // was not a generated frequency
      if (modulus >= amplTol) {
        ++nUnexpectedFrequencies;
        std::cout << "Frequency: " << (freq/1e6) << " MHz  amplitude: "
          << amplitude << " (modulus: " << modulus
          << ", phase: " << phase/util::pi() << "π rad)"
          << " [SPURIOUS]" << std::endl;
      }
    }
    else { // was a generated frequency
      ++nFoundFreqs;
      
      Harmonic_t<value_t> const& harmonic
        = harmonics.at(std::distance(fbegin, itFreq));
      std::cout << "Frequency: " << (freq/1e6) << " MHz  amplitude: "
        << rnd(amplitude) << " (modulus: " << rnd(modulus)
        << ", phase: " << roundAt(phase/util::pi(), 1e-5) << "π rad → "
        << rnd(2*amplitude.real()) << "·cos ωt "
        << std::showpos << rnd(2*amplitude.imag()) << "·sin ωt)"
        << "; generated as " << rnd(harmonic.amplitude)
        << " (δ=" << rnd(2*modulus - harmonic.amplitude) << ") phase "
        << rnd(harmonic.phase/util::pi()) << "π (δ=" 
        << rnd(radDiff(phase, harmonic.phase)/util::pi())
        << "π rad) = " << rnd(harmonic.amplitude*std::cos(harmonic.phase)) << "·cos ωt "
        << std::showpos << rnd(harmonic.amplitude*std::sin(harmonic.phase)) << "·sin ωt"
        << std::endl;
      BOOST_TEST(2.*modulus == harmonic.amplitude, 0.01% boost::test_tools::tolerance());
      BOOST_TEST(std::abs(radDiff(phase, harmonic.phase)) < 0.001); // rad
    }
  }
  
  // this may be an issue with the testing algorithm, like an improper match
  // tolerance or the choice of frequencies that is too far from the harmonics
  // of the transform.
  BOOST_TEST(nFoundFreqs == testFreqs.size());
  BOOST_TEST(nUnexpectedFrequencies == 0);
  
} // checkSpectrum()


template <typename Waveform>
bool checkWaveform(Waveform const& waveform, Waveform const& generated) {
  
  using value_t = typename Waveform::value_type;
  
  auto findAbsMax = [](Waveform const& samples)
    {
      value_t maxAmpl = 0;
      for (value_t const sample: samples) {
        auto const as = std::abs(sample);
        if (as > maxAmpl) maxAmpl = as;
      }
      return maxAmpl;
    };
  
  value_t const amplTol = findAbsMax(generated) / 1e4;
  lar::util::MinMaxCollector<double> extremes;
  unsigned int nDiff = 0;
  for (auto [ sim, gen ]: util::zip(waveform, generated)) {
    double const d = std::abs(sim - gen);
    extremes.add(d);
    if (d <= amplTol) continue;
    ++nDiff;
  } // for
  
  BOOST_TEST(waveform.size() == generated.size());
  BOOST_TEST(nDiff == 0);
  if (nDiff > 0) {
    std::cout << "Absolute differences ranged between " << extremes.min()
      << " and " << extremes.max() << std::endl;
  }
  
  return nDiff == 0;
} // checkWaveform()


// -----------------------------------------------------------------------------
// --- tests
// -----------------------------------------------------------------------------
template <typename Coll>
void transformTest(
  microsecond duration,
  nanosecond samplingPeriod,
  Coll const& harmonics
) {
  /*
   * This test performs transformation back and forth without any filtering.
   * 
   * Checks include the verification of the spectrum (`checkSpectrum()`)
   * and of the final waveform (`checkWaveform()`).
   */
  
  // type used for data (and sometimes times)
  using value_t = typename Coll::value_type::value_type;
  
  using Samples_t = std::valarray<value_t>; // type of waveform samples
  
  //
  // create the test waveform
  //
  
  /*
   * Reminder: a waveform T with sampling time dt and N samples will host N/2
   *           frequencies from 1/T to N/2T in steps of 1/T, plus a baseline.
   *           For example, a T = 10 us waveform sampled at dt = 2 ns, with
   *           5000 samples will have periods up to 10 us and frequencies from
   *           1/10 us = 100 kHz to 5000 / (2 x 10 us) = 250 MHz (for a total
   *           of 2501: 0 Hz, 100 kHz, 200 kHz, ... , 249.9 MHz, 250 MHz).
   */
  
  std::size_t const nSamples = adjustTimes(duration, samplingPeriod);
  Samples_t const samples = generateWaveform(nSamples, samplingPeriod, harmonics);
  
  //
  // initialize the algorithm
  //
  icarus::trigger::WaveformFilterAlg filterAlg{ 1.0 / samplingPeriod };
  
  // std::cout << "Processing PMT sum with " << samples.size() << " samples: "
  //   << dumpCollection(samples, 16); // FIXME
  auto const frequencies = filterAlg.spectrumFrequencies(nSamples);
  
  //
  // perform spectrum extraction
  //
  
  auto spectrum = filterAlg.waveformToSpectrum(samples);
  BOOST_TEST(spectrum.size() == samples.size() / 2 + 1);
  
  // test of the spectrum components
  checkSpectrum(harmonics, duration, frequencies, spectrum);
  
  // test of the waveform built out of the spectrum
  Samples_t const simSamples
    = waveformFromSpectrum(frequencies, spectrum, samplingPeriod, nSamples);
  checkWaveform(simSamples, samples);
  
  //
  // perform filtering and test
  //
  // double const samplingTime = second{ samplingPeriod }.value(); // -> s
  auto const reshapedSpectrum = spectrum; // TODO
  
  // std::cout << "spectrum reshaped (sampling time: " << samplingTime
  //   << "s, " << reshapedSpectrum.size() << " frequencies)"
  //   << dumpCollection(reshapedSpectrum, 8); // FIXME
  
  //
  // perform reformation of waveform and test
  //
  
  auto reshapedWaveform = filterAlg.spectrumToWaveform(std::move(reshapedSpectrum));
  // std::cout << "Reshaped PMT sum returned to a " << reshapedWaveform.size()
  //   << "-sample waveform: \n" << dumpCollection(reshapedWaveform, 16); // FIXME
  static_assert(std::is_same_v<decltype(reshapedWaveform), Samples_t>);
  
  checkWaveform(reshapedWaveform, samples);

} // transformTest()


template <typename Coll>
void oneStepTransformTest(
  microsecond duration,
  nanosecond samplingPeriod,
  Coll const& harmonics
) {
  /*
   * This test performs transformation back and forth without any filtering.
   * It uses void WaveformFilterAlg::applyInPlace() to perform all the cycle
   * in one call.
   */
  
  // type used for data (and sometimes times)
  using value_t = typename Coll::value_type::value_type;
  
  using Samples_t = std::valarray<value_t>; // type of waveform samples
  
  std::size_t const nSamples = adjustTimes(duration, samplingPeriod);
  Samples_t const samples = generateWaveform(nSamples, samplingPeriod, harmonics);
  
  //
  // initialize the algorithm
  //
  icarus::trigger::WaveformFilterAlg filterAlg{ 1.0 / samplingPeriod };
  
  // std::cout << "Processing PMT sum with " << samples.size() << " samples: "
  //   << dumpCollection(samples, 16); // FIXME
  auto const frequencies = filterAlg.spectrumFrequencies(nSamples);
  
  //
  // perform spectrum extraction, filtering and reconstitution at once:
  //
  auto const identity = [](...) -> value_t { return static_cast<value_t>(1.0); };
  Samples_t reshapedWaveform = samples;
  filterAlg.applyInPlace(reshapedWaveform, identity);
  BOOST_TEST(reshapedWaveform.size() == samples.size());
  
  // std::cout << "Reshaped PMT sum returned to a " << reshapedWaveform.size()
  //   << "-sample waveform: \n" << dumpCollection(reshapedWaveform, 16); // FIXME
  checkWaveform(reshapedWaveform, samples);

} // oneStepTransformTest()


// -----------------------------------------------------------------------------
template <typename Coll>
void threeBandFilterTest(
  microsecond duration,
  nanosecond samplingPeriod,
  Coll const& harmonics
) {
  /*
   * This test performs transformation, filtering and regeneration.
   * 
   * Checks include the verification of the spectrum (`checkSpectrum()`)
   * and of the final waveform (`checkWaveform()`).
   */
  
  // type used for data (and sometimes times)
  using value_t = typename Coll::value_type::value_type;
  
  using Samples_t = std::valarray<value_t>; // type of waveform samples
  
  //
  // create the test waveform
  //
  
  /*
   * Reminder: a waveform T with sampling time dt and N samples will host N/2
   *           frequencies from 1/T to N/2T in steps of 1/T, plus a baseline.
   *           For example, a T = 10 us waveform sampled at dt = 2 ns, with
   *           5000 samples will have periods up to 10 us and frequencies from
   *           1/10 us = 100 kHz to 5000 / (2 x 10 us) = 250 MHz (for a total
   *           of 2501: 0 Hz, 100 kHz, 200 kHz, ... , 249.9 MHz, 250 MHz).
   */
  
  std::size_t const nSamples = adjustTimes(duration, samplingPeriod);
  Samples_t const samples = generateWaveform(nSamples, samplingPeriod, harmonics);
  
  ThreeLevelHighPassFrequencyFilter const filter
    = makeThreeLevelHighPassFrequencyFilter(harmonics, duration);
  auto const filteredHarmonics = filterHarmonics(harmonics, duration, filter);
  Samples_t const expectedReshapedSamples
    = generateWaveform(nSamples, samplingPeriod, filteredHarmonics);
  
  //
  // initialize the algorithm
  //
  icarus::trigger::WaveformFilterAlg filterAlg{ 1.0 / samplingPeriod };
  
  // std::cout << "Processing PMT sum with " << samples.size() << " samples: "
  //   << dumpCollection(samples, 16); // FIXME
  auto const frequencies = filterAlg.spectrumFrequencies(nSamples);
  
  
  //
  // perform spectrum extraction
  //
  
  auto spectrum = filterAlg.waveformToSpectrum(samples);
  BOOST_TEST(spectrum.size() == samples.size() / 2 + 1);
  
  // test of the spectrum components
  checkSpectrum(harmonics, duration, frequencies, spectrum);
  
  // test of the waveform built out of the spectrum
  Samples_t const simSamples
    = waveformFromSpectrum(frequencies, spectrum, samplingPeriod, nSamples);
  checkWaveform(simSamples, samples);
  
  //
  // perform filtering and test
  //
  
  // adapter to the required interface: std::complex<value_t>(*)(value_t)
  auto const filterFunc = [&filter](value_t ampl){ return filter(hertz{ ampl }); };
  
  auto const reshapedSpectrum
    = filterAlg.filterSpectrum(frequencies, spectrum, filterFunc);
  
  // std::cout << "spectrum reshaped (sampling time: " << samplingTime
  //   << "s, " << reshapedSpectrum.size() << " frequencies)"
  //   << dumpCollection(reshapedSpectrum, 8); // FIXME
  
  //
  // perform reformation of waveform and test
  //
  
  auto const reshapedWaveform = filterAlg.spectrumToWaveform(reshapedSpectrum);
  // std::cout << "Reshaped PMT sum returned to a " << reshapedWaveform.size()
  //   << "-sample waveform: \n" << dumpCollection(reshapedWaveform, 16); // FIXME
  static_assert(std::is_same_v<decltype(reshapedWaveform), Samples_t const>);
  
  checkWaveform(reshapedWaveform, expectedReshapedSamples);
  
} // threeBandFilterTest()


// -----------------------------------------------------------------------------
template <typename Coll>
void oneStepThreeBandFilterTest(
  microsecond duration,
  nanosecond samplingPeriod,
  Coll const& harmonics
) {
  /*
   * This test performs transformation back and forth without any filtering.
   * It uses void WaveformFilterAlg::applyInPlace() to perform all the cycle
   * in one call.
   */
  
  // type used for data (and sometimes times)
  using value_t = typename Coll::value_type::value_type;
  
  using Samples_t = std::valarray<value_t>; // type of waveform samples
  
  //
  // create the test waveform
  //
  
  /*
   * Reminder: a waveform T with sampling time dt and N samples will host N/2
   *           frequencies from 1/T to N/2T in steps of 1/T, plus a baseline.
   *           For example, a T = 10 us waveform sampled at dt = 2 ns, with
   *           5000 samples will have periods up to 10 us and frequencies from
   *           1/10 us = 100 kHz to 5000 / (2 x 10 us) = 250 MHz (for a total
   *           of 2501: 0 Hz, 100 kHz, 200 kHz, ... , 249.9 MHz, 250 MHz).
   */
  
  std::size_t const nSamples = adjustTimes(duration, samplingPeriod);
  Samples_t const samples = generateWaveform(nSamples, samplingPeriod, harmonics);
  
  ThreeLevelHighPassFrequencyFilter const filter
    = makeThreeLevelHighPassFrequencyFilter(harmonics, duration);
  auto const filteredHarmonics = filterHarmonics(harmonics, duration, filter);
  Samples_t const expectedReshapedSamples
    = generateWaveform(nSamples, samplingPeriod, filteredHarmonics);
  
  //
  // initialize the algorithm
  //
  icarus::trigger::WaveformFilterAlg filterAlg{ 1.0 / samplingPeriod };
  
  // std::cout << "Processing PMT sum with " << samples.size() << " samples: "
  //   << dumpCollection(samples, 16); // FIXME
  auto const frequencies = filterAlg.spectrumFrequencies(nSamples);
  
  //
  // perform spectrum extraction, filtering and reconstitution at once:
  //
  // adapter to the required interface: std::complex<value_t>(*)(value_t)
  auto const filterFunc
    = [&filter](value_t ampl){ return filter(hertz{ ampl }); };
  
  Samples_t reshapedWaveform = samples;
  filterAlg.applyInPlace(reshapedWaveform, filterFunc);
  BOOST_TEST(reshapedWaveform.size() == samples.size());
  
  // std::cout << "Reshaped PMT sum returned to a " << reshapedWaveform.size()
  //   << "-sample waveform: \n" << dumpCollection(reshapedWaveform, 16); // FIXME
  checkWaveform(reshapedWaveform, expectedReshapedSamples);

} // oneStepThreeBandFilterTest()


// -----------------------------------------------------------------------------
void WaveformFilterAlg_3freq_transform_test() {
  
  /*
   * This test performs filtering on three frequencies, suppressing the lowest
   * one completely, halving the middle one and leaving the third untouched
   * (on good approximation).
   * 
   * The test waveform is 10 us/2 ns = 5000 samples. Its frequencies
   * are chosen so that the lowest has 20 periods in the waveform, the middle
   * one 51 and the highest one 999.
   * 
   */
  std::array const testHarmonics = to_array<Harmonic_t<float>>({
      {  2.0,  20,  3.0 }
    , { 10.0,  51, -2.0, util::pi()/2 }
    , {  0.5, 999,  4.0, util::pi() }
    });
  
  std::cout << "WaveformFilterAlg 3-freq step transform test" << std::endl;
  BOOST_TEST_INFO_SCOPE("WaveformFilterAlg 3-freq step transform test");
  transformTest(10_us, 2_ns, testHarmonics);
  
  std::cout << "WaveformFilterAlg 3-freq monolithic transform test" << std::endl;
  BOOST_TEST_INFO_SCOPE("WaveformFilterAlg 3-freq monolithic transform test");
  oneStepTransformTest(10_us, 2_ns, testHarmonics);
  
} // WaveformFilterAlg_3freq_transform_test()


// -----------------------------------------------------------------------------
void WaveformFilterAlg_3freq_threeBandFilterTest_test() {
  
  /*
   * This test performs filtering on three frequencies, suppressing the lowest
   * one completely, halving the middle one and leaving the third untouched
   * (on good approximation).
   * 
   * The test waveform is 10 us/2 ns = 5000 samples. Its frequencies
   * are chosen so that the lowest has 20 periods in the waveform, the middle
   * one 51 and the highest one 999.
   * 
   */
  std::array const testHarmonics = to_array<Harmonic_t<float>>({
      {  2.0,  20,  3.0 }
    , { 10.0,  51, -2.0, util::pi()/2 }
    , {  0.5, 999,  4.0, util::pi() }
    });
  
  std::cout << "WaveformFilterAlg 3-freq step 3-band filter test" << std::endl;
  BOOST_TEST_INFO_SCOPE("WaveformFilterAlg 3-freq step 3-band filter test");
  threeBandFilterTest(10_us, 2_ns, testHarmonics);
  
  std::cout << "WaveformFilterAlg 3-freq monolithic 3-band filter test" << std::endl;
  BOOST_TEST_INFO_SCOPE("WaveformFilterAlg 3-freq monolithic 3-band filter test");
  oneStepThreeBandFilterTest(10_us, 2_ns, testHarmonics);
  
} // WaveformFilterAlg_3freq_threeBandFilterTest_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(WaveformFilterAlg_testcase) {
  
  WaveformFilterAlg_3freq_transform_test();
  
  WaveformFilterAlg_3freq_threeBandFilterTest_test();
  
} // BOOST_AUTO_TEST_CASE(GateOpeningInfoExtractor_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
