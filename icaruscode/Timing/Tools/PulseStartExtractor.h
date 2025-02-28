/*
 * File: PulseStartExtractor.h
 * Author: M. Vicenzi (mvicenzi@bnl.gov)
 * Description: Header file for the PulseStartExtractor class, which provides methods
 *              templated functions to extract the start (rising edge) time of reference
 *              signals (square pulses) using either constant fraction discrimination or
 *              a refined sigmoid fit.
 */

#ifndef PULSE_START_EXTRACTOR_H
#define PULSE_START_EXTRACTOR_H

#include <vector>
#include <cstddef>
#include <limits>
#include <string>
#include "TF1.h"
#include "TGraph.h"
#include "TFitResultPtr.h"
#include "TRandom.h"

class PulseStartExtractor
{
public:
  enum ExtractionMethod
  {
    CONSTANT_FRACTION,
    SIGMOID_FIT
  };

  // Constructor: select extraction method and threshold (used in the CF method).
  PulseStartExtractor(ExtractionMethod method = CONSTANT_FRACTION, double threshold = 10.0);

  // Templated public function to extract the start sample from a waveform.
  // Returns the start sample index (std::size_t).
  template <typename T>
  double extractStart(const std::vector<T> &wf);

private:
  ExtractionMethod fMethod;
  double fThreshold;

  // Templated helper: constant-fraction extraction.
  template <typename T>
  std::size_t extractStartSampleCF(const std::vector<T> &wf, double thres);

  // Templated helper: finds the index of the minimum element in the range [start, end).
  template <typename T>
  std::size_t findMinBin(const std::vector<T> &wf, std::size_t start, std::size_t end);

  // Templated helper: finds the index of the maximum element in the range [start, end).
  template <typename T>
  std::size_t findMaxBin(const std::vector<T> &wf, std::size_t start, std::size_t end);

  // Templated helper: computes the median of a vector.
  template <typename T>
  T computeMedian(const std::vector<T> &data);

  // Templated helper: sigmoid-fit extraction; fills par[4] with fit parameters.
  // Returns true if the fit converged.
  template <typename T>
  bool extractStartSigmoidFit(const std::vector<T> &wf, std::size_t guess, double par[4]);
};

#endif // PULSE_START_EXTRACTOR_H
