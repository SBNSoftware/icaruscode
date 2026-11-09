/**
 * @file icaruscode/Timing/Tools/PulseStartExtractor.h
 * @brief `icarus::timing::PulseStartExtractor` header.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date  Sun Feb 11 11:37:14 2024
 * Description: Header file for the PulseStartExtractor class implementation. It provides
 *              templated functions to extract the start (rising edge) time of reference
 *              signals (square pulses) using either constant fraction discrimination,
 *              an interpolated constant fraction or a refined logistic function fit.
 **/

#ifndef PULSE_START_EXTRACTOR_H
#define PULSE_START_EXTRACTOR_H

#include <array>
#include <vector>
#include <cstddef>
#include <string>
#include <map>

namespace icarus::timing {

  enum ExtractionMethod
  {
    CONSTANT_FRACTION,
    LOGISTIC_FIT,
    INTERPOLATED_CF
  };

  // map for string to enum conversion
  std::map<std::string, ExtractionMethod> const stringToExtractionMethod = {
    { "LogisticFit", LOGISTIC_FIT },
    { "ConstantFraction", CONSTANT_FRACTION },
    { "InterpolatedCF", INTERPOLATED_CF }
  };

/**
 * @brief Extracts the start time of a pulse from a waveform.
 *
 * This class implements methods for extracting the start time of reference pulses
 * digitized in PMT special channels. These reference waveforms are expected to exhibit
 * a negative-polarity square pulse with a sharp rising edge. Because of the RC coupling,
 * the signal does not hold a plateau: it reaches its minimum shortly after the edge
 * (from a few to several tens of samples, depending on the channel) and then decays back.
 *
 * The extraction can be performed using either a constant-fraction discrimination (CF),
 * an interpolated constant-fraction discrimination or a refined logistic function fit.
 * @see icarus::timing::ExtractionMethod for the available methods.
 *
 * Inputs
 * -------------------------
 *
 * The extraction method is configured via the constructor by specifying one of the available
 * methods (defined by the `icarus::timing::ExtractionMethod` enum) and a threshold value.
 * * `method` (`icarus::timing::ExtractionMethod`): specifies the selected method for the
 *    pulse start extraction. The default is `icarus::timing::CONSTANT_FRACTION`.
 * * `threshold` (`double`): detection threshold in ADC counts to avoid selecting noise if
 *    the reference pulse is missing from its waveform. Default is `500` ADC.
 *
 * The other parameters of the algorithms are defined as constants at the top of the class.
 *
 * Pulse finding
 * ----------------------------------
 *
 * The waveform is expected to contain either one large pulse or only noise.
 * All methods share the same pulse finding:
 * * the absolute minimum of the waveform is found;
 * * the baseline is the median of the waveform from its start up to `kLookBackSamples`
 *   samples before the minimum (at least the first sample is used);
 * * if the baseline-minimum difference (the amplitude) is below threshold, the waveform is
 *   assumed to be noise and the start sample (`0`) is returned;
 * * within the `kLookBackSamples` samples before the minimum, the algorithm searches the
 *   first sample where the signal deviates from the baseline by at least `kCFFraction`
 *   (20%) of the amplitude.
 *
 * ### Constant-fraction discrimination
 *
 * The sample immediately before the 20% crossing is returned as the pulse start time.
 * This algorithm is affected by the 2 ns discretization of the input waveforms as there
 * is no interpolation between the samples.
 *
 * ### Interpolated constant-fraction discrimination
 *
 * The 20% crossing time is linearly interpolated between the two samples across it.
 *
 * ### Logistic Fit Method:
 * The pulse rising edge is fitted using a logistic function: `1/(1+exp(([0]-x)/[1]))*[2]+[3]`.
 * The estimate for the pulse start is then taken to be the inflection point (ie. the steepest point,
 * where the first derivative is maximal) which corresponds to parameter `p0`.
 * It also represents the point at 50% of the amplitude.
 *
 * The procedure is the following:
 * * the initial guess for the pulse start time (`p0`) is the constant-fraction result;
 * * the initial guess for the baseline (`p3`) is the baseline from the pulse finding;
 * * the initial guess for the amplitude (`p2`) is estimated averaging the first
 *   `kFitAmplitudeSamples` samples after the initial guess for the rise time (`p0`);
 * * the logistic function fit is performed over a window of `kFitHalfWindow` samples
 *   around the initial guess for the rise time;
 * * if the fit converges with `p0` within the fit window and a negative amplitude (`p2`),
 *   `p0` is returned as the pulse start time;
 * * otherwise, the fit is retried with the other initial values in `kFitStartOffsets` and
 *   `kFitSlopeSeeds`. If all attempts fail, the interpolated 50% crossing time is returned.
 *
 * This algorithm allows to interpolate between samples, achieving a better timing resolution.
 *
 */
  class PulseStartExtractor
  {
  public:

    /// Samples before the minimum where the start is searched; the baseline region ends there.
    static constexpr std::size_t kLookBackSamples = 100;
    /// Fraction of the amplitude used by the constant-fraction methods.
    static constexpr double kCFFraction = 0.2;
    /// Samples after the initial start guess averaged for the fit amplitude guess.
    static constexpr std::size_t kFitAmplitudeSamples = 50;
    /// Fit window half width around the initial start guess [samples].
    static constexpr std::size_t kFitHalfWindow = 25;
    /// Offsets to the initial start guess tried by the fit [samples].
    static constexpr std::array<double, 3> kFitStartOffsets { 0., -1., 1. };
    /// Initial values of the fit slope parameter (`p1`) tried by the fit [samples].
    static constexpr std::array<double, 4> kFitSlopeSeeds { 0.3, 0.1, 0.5, 1.0 };

    // Constructor: select extraction method and threshold in ADC.
    PulseStartExtractor(ExtractionMethod method = CONSTANT_FRACTION, double threshold= 500);

    // Templated public function to extract the start sample from a waveform.
    // Returns the start sample (can be fractional), or `0` if no pulse is found.
    template <typename T>
    double extractStart(const std::vector<T> &wf) const;

  private:
    ExtractionMethod fMethod;
    double fThreshold;

    /// Pulse found in the waveform.
    struct Pulse
    {
      std::size_t minBin = 0;       ///< sample of the absolute minimum
      std::size_t searchStart = 0;  ///< first sample of the start search window
      double baseline = 0.;         ///< baseline level [ADC]
      double amplitude = 0.;        ///< baseline-minimum difference [ADC]
    };

    // Finds the pulse in the waveform; returns false if only noise is found.
    bool findPulse(std::vector<double> const &wf, Pulse &pulse) const;

    // Returns the crossing time of `fraction` of the pulse amplitude below baseline.
    // If `interpolate` is false, the last sample before the crossing is returned.
    double crossingTime(std::vector<double> const &wf, Pulse const &pulse,
                        double fraction, bool interpolate) const;

    // Logistic-fit extraction seeded with `guess`; returns false if no fit is acceptable.
    bool fitLogistic(std::vector<double> const &wf, Pulse const &pulse, double guess, double &start) const;

  };

} // icarus::timing namespace

#endif // PULSE_START_EXTRACTOR_H
