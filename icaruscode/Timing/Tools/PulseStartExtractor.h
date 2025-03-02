/**
 * @file icaruscode/Timing/Tools/PulseStartExtractor.h
 * @brief `icarus::timing::PulseStartExtractor` header.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date  Sun Feb 11 11:37:14 2024
 * Description: Header file for the PulseStartExtractor class implementation. It provides
 *              templated functions to extract the start (rising edge) time of reference
 *              signals (square pulses) using either constant fraction discrimination or
 *              a refined logistic function fit.
 **/

#ifndef PULSE_START_EXTRACTOR_H
#define PULSE_START_EXTRACTOR_H

#include <vector>
#include <cstddef>
#include <limits>
#include <string>
#include <map>
#include "TF1.h"
#include "TGraph.h"
#include "TFitResultPtr.h"
#include "TRandom.h"

namespace icarus::timing {

  enum ExtractionMethod
  {
    CONSTANT_FRACTION,
    LOGISTIC_FIT
  };

  // map for string to enum conversion 
  std::map<std::string, ExtractionMethod> const stringToExtractionMethod = { 
    { "LogisticFit", LOGISTIC_FIT },
    { "ConstantFraction", CONSTANT_FRACTION }
  };

/**
 * @brief Extracts the start time of a pulse from a waveform.
 * 
 * This class implements methods for extracting the start time of reference pulses
 * digitized in PMT special channels. These reference waveforms are expected to exhibit 
 * a negative-polarity square pulse with a sharp rising edge.
 * 
 * The extraction can be performed using either a legacy constant-fraction discrimination (CF)
 * algorithm or a refined logistic function fit.
 * @see icarus::timing::ExtractionMethod for the available methods.
 * 
 * Inputs
 * -------------------------
 * 
 * The extraction method is configured via the constructor by specifying one of the available 
 * methods (defined by the `icarus::timing::ExtractionMethod` enum) and a threshold value.
 * * `method` (`icarus::timing::ExtractionMethod`): specifies the selected method for the
 *    pulse start extraction. The default is `icarus::timing::CONSTANT_FRACTION`.
 * * `threshold` (`double`): detection threshold in ADC counts to avoid selecting noise of 
 *    the reference pulse is missing from its waveform. Default id `500` ADC.
 * 
 * When using a fit method, the initial parameters guesses are determined using the 
 * constant-fraction algorithm to help ensure fit convergence.
 *
 * Signal timing extraction
 * ----------------------------------
 * 
 * The reference waveform is expected to exhibit a negative-polarity square wave.
 * Its selection is performed simply by looking at the absolute minimum of the waveform.
 * A detection threshold on the amplitude is used to exclude picking up spurios signals
 * in case the large square pulse is missing (which could happen for EW/RWM signals).
 * 
 * ### Constant-fraction discrimination
 * 
 * The procedure is the following:
 * * the absolute minimum of the waveform is found;
 * * an interval starting 20 ticks before that minimum is considered and the baseline level
 *   is defined as the value at the start of that interval;
 * * if the baseline-minimum difference is below a threshold, it is assumed to be noise
 *   and the start sample (`0`) is returned; 
 * * within the defined window, the algorithm searches the first sample where the signal
 *   deviates from the baseline by at least 20% of the full amplitude (baseline-minimum difference);
 *   The sample immediately before is returned as the pulse start time.
 * 
 * This algorithm is affected by the 2ns discretization of the input waveforms as there is no 
 * interpolation between the samples.
 * 
 * ### Logistic Fit Method:
 * The pulse rising edge is fitted using a logistic function: `1/(1+exp(([0]-x)/[1]))*[2]+[3]`.
 * The estimate for the pulse start is then taken to be the inflection point (ie. the steepest point,
 * where the first derivative is maximal) which corresponds to parameter `p0`. 
 * It also represents the point at 50% of the amplitude. 
 * 
 * The procedure is the following:
 * * the constant-fraction method is ran first to obtain an initial guess for the pulse 
 *   start time (`p0`). If no signal is found, the start sample (`0`) is returned;
 * * an initial guess for the baseline (`p3`) is estimated as the median of the full waveform;
 * * an initial guess for the amplitude (`p2`) is estimated averaging the first 50 samples
 *   after the initial guess for the rise time (`p0`).
 * * the logistic function fit is performed over a window of `[-25,25]` samples around the
 *   initial guess for the rise time;
 * * if the fit converges, `p0` is returned as the pulse start time;
 * * if the fit fails to converge, multiple retries with slight random variations around the
 *   initial parameters are attempted. If all retries are exahusted, it falls back to
 *   returning the constant-fraction method’s result.
 *
 * This algorithm allows to interpolate between samples, achieving a better timing resolution.
 * 
 */
  class PulseStartExtractor
  {
  public:


    // Constructor: select extraction method and threshold in ADC (used in the CF method).
    PulseStartExtractor(ExtractionMethod method = CONSTANT_FRACTION, double threshold= 500);

    // Templated public function to extract the start sample from a waveform.
    // Returns the start sample index (std::size_t).
    template <typename T>
    double extractStart(const std::vector<T> &wf) const;

  private:
    ExtractionMethod fMethod;
    double fThreshold;

    // Templated helper: constant-fraction extraction.
    template <typename T>
    std::size_t extractStartSampleCF(const std::vector<T> &wf, double thres) const;

    // Templated helper: finds the index of the minimum element in the range [start, end).
    template <typename T>
    std::size_t findMinBin(const std::vector<T> &wf, std::size_t start, std::size_t end) const;

    // Templated helper: finds the index of the maximum element in the range [start, end).
    template <typename T>
    std::size_t findMaxBin(const std::vector<T> &wf, std::size_t start, std::size_t end) const;

    // Templated helper: computes the median of a vector.
    template <typename T>
    T computeMedian(const std::vector<T> &data) const;

    // Templated helper: logistic-fit extraction; fills par[4] with fit parameters.
    // Returns true if the fit converged.
    template <typename T>
    bool extractStartLogisticFit(const std::vector<T> &wf, std::size_t guess, double par[4]) const;
  };

} // icarus::timing namespace

#endif // PULSE_START_EXTRACTOR_H
