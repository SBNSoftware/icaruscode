/**
 * @file icaruscode/Timing/Tools/PulseStartExtractor.cc
 * @brief `icarus::timing::PulseStartExtractor` source.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date  Sun Feb 11 11:37:14 2024
 * Description: Source file for the PulseStartExtractor class implementation. It provides
 *              templated functions to extract the start (rising edge) time of reference
 *              signals (square pulses) using either constant fraction discrimination or
 *              a refined logistic function fit.
 **/

#include "PulseStartExtractor.h"
#include <algorithm>
#include <iostream>
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
  
  // Constructor
  PulseStartExtractor::PulseStartExtractor(ExtractionMethod method, double threshold)
      : fMethod(method), fThreshold(threshold) {}

  // --------------------------------------------------------------------------------
  template <typename T>
  std::size_t PulseStartExtractor::findMinBin(const std::vector<T> &wf, std::size_t start, std::size_t end) const
  {
    auto minel = std::min_element(wf.begin() + start, wf.begin() + end);
    size_t minsample = std::distance(wf.begin() + start, minel);
    return minsample;
  }

  template <typename T>
  std::size_t PulseStartExtractor::findMaxBin(const std::vector<T> &wf, std::size_t start, std::size_t end) const
  {
    auto maxel = std::max_element(wf.begin() + start, wf.begin() + end);
    size_t maxsample = std::distance(wf.begin() + start, maxel);
    return maxsample;
  }

  // --------------------------------------------------------------------------------
  template <typename T>
  T PulseStartExtractor::computeMedian(const std::vector<T> &data) const
  {
    std::vector<T> copy = data;
    std::nth_element(copy.begin(), copy.begin() + copy.size() / 2, copy.end());
    return copy[copy.size() / 2];
  }

  // --------------------------------------------------------------------------------
  template <typename T>
  double PulseStartExtractor::extractStart(const std::vector<T> &wf) const
  {

    // Try constant-fraction discrimination (default method)
    // if fit method is selected, the result is used as initial guess for the fit
    std::size_t cf_sample = extractStartSampleCF(wf, fThreshold);

    // if no signal is found above threshold, first sample (0) is returned
    // no need to attempt a fit if no signal is there
    if (cf_sample == 0)
      return 0.;

    if (fMethod == LOGISTIC_FIT)
    {
      double par[4] = {0};
      bool ok = extractStartLogisticFit(wf, cf_sample, par);
      if (ok)
        return par[0]; // Rising edge time from the logistic fit (can be non-integer).
      else
        std::cout << "Sigmoid fit did not converge; falling back to constant-fraction discrimation." << std::endl;
    }

// ***** MODOFICATIONS ***** //

if (fMethod == INTERPOLATED_CF) {
      return extractStartInterpolated(wf, 0.2, 2.0); //above 20% fraction (threshold), 2 ns sampling (this is the t0 -t1, t2-t1...fixing it)

	std::cout << "Using interpolated CF method..." << std::endl;

    }

    // Defaulting to constant-fraction discrimation
    return static_cast<double>(cf_sample);
  }

  // --------------------------------------------------------------------------------

  template <typename T>
  std::size_t PulseStartExtractor::extractStartSampleCF(const std::vector<T> &wf, double thres) const
  {
    // Find the position of the minimum along the waveform
    std::size_t minbin = findMinBin(wf, 0, wf.size());

    // Search only a cropped region of the waveform backward from the min (20 samples)
    std::size_t maxbin = (minbin >= 20) ? (minbin - 20) : 0;
    std::size_t startbin = maxbin;

    // check if peak is high enough (remove noise)
    auto delta = wf[maxbin] - wf[minbin];
    if (delta < thres)
      return 0; // return first sample

    // Now we crawl betweem maxbin and minbin and we stop when:
    // maxbin value - bin value > (maxbin value - minbin value )*0.2
    for (std::size_t bin = maxbin; bin < minbin; ++bin)
    {
      auto val = wf[maxbin] - wf[bin];
      if (val >= 0.2 * delta) // 20%
      {
        startbin = bin - 1;
        break;
      }
    }
    if (startbin < maxbin)
      startbin = maxbin;

    return startbin;
  }

  // --------------------------------------------------------------------------------

  template <typename T>
  bool PulseStartExtractor::extractStartLogisticFit(const std::vector<T> &wf, std::size_t guess, double par[4]) const
  {
    std::size_t nsize = wf.size();
    std::vector<double> x(nsize);
    std::vector<double> y(nsize);

    for (std::size_t i = 0; i < nsize; i++)
    {
      x[i] = static_cast<double>(i);
      y[i] = static_cast<double>(wf[i]);
    }

    TGraph *g = new TGraph(nsize, &x[0], &y[0]);

    // prepare initial parameter guesses for the fit
    // baseline: level before signal
    auto baseline = computeMedian(y);

    // amplitude: level after the signal (avg over 50 samples)
    auto amplitude = 0.0;
    std::size_t endIndex = (guess + 50 < nsize) ? guess + 50 : nsize;
    for (std::size_t i = guess; i < endIndex; i++)
      amplitude += y[i];
    amplitude /= (endIndex - guess);
    amplitude -= baseline;

    // crop a window around the initial guess: +/- 25 samples
    std::size_t llim = (guess >= 25) ? guess - 25 : 0;
    std::size_t ulim = (guess + 25 < nsize) ? guess + 25 : nsize - 1;

    std::string logisticModel = "1/(1+exp(([0]-x)/[1]))*[2]+[3]";
    TF1 *fitFunc = new TF1("logistic", logisticModel.c_str(), x[llim], x[ulim]);

    fitFunc->SetParameter(0, x[guess]); // initial guess for rising edge time.
    fitFunc->SetParameter(1, 0.3);      // initial slope.
    fitFunc->SetParameter(2, amplitude);
    fitFunc->SetParameter(3, baseline);

    TFitResultPtr fp = g->Fit(fitFunc, "SRQ", "", x[llim], x[ulim]);
    bool converged = !(bool)(int(fp));

    // if fit did not converge, try again a few times
    // tweak some parameters until it works...
    int maxtries = 100;
    int trial = 0;
    while (!converged && trial < maxtries)
    {
      double newGuess = x[guess] + gRandom->Uniform(-1.5, 1.5);
      double newSlope = 0.3 + gRandom->Uniform(-0.2, 0.2);
      fitFunc->SetParameter(0, newGuess);
      fitFunc->SetParameter(1, newSlope);
      fitFunc->SetParameter(2, amplitude);
      fitFunc->SetParameter(3, baseline);
      fp = g->Fit(fitFunc, "SRQ", "", x[llim], x[ulim]);
      converged = !(bool)(int(fp));
      trial++;
    }

    // return the fit results in the parameters vector.
    // p0 : rising edge as the inflection point (50% of amplitude)
    par[0] = fitFunc->GetParameter(0);
    par[1] = fitFunc->GetParameter(1);
    par[2] = fitFunc->GetParameter(2);
    par[3] = fitFunc->GetParameter(3);

    delete fitFunc;
    delete g;
    return converged;
  }


// ***** MODIFICATIONS ***** //

template <typename T>
double PulseStartExtractor::extractStartInterpolated(const std::vector<T>& wf, double fraction, double dt_ns) const {
  if (wf.size() < 2) return 0.0;

  double baseline = computeMedian(wf);
  auto minIt = std::min_element(wf.begin(), wf.end());
  std::size_t minIndex = std::distance(wf.begin(), minIt);
  double peak = wf[minIndex];

  double threshold = baseline + fraction * (peak - baseline);

  std::size_t maxbin = (minIndex >= 20) ? (minIndex - 20) : 0;

 for (std::size_t i = maxbin; i + 1 < minIndex; ++i) {
    if (wf[i] > threshold && wf[i + 1] <= threshold) {
      // Interpolate between sample i and i+1
      double y1 = wf[i];
      double y2 = wf[i + 1];
      double frac = (y1 - threshold) / (y1 - y2);  // fractional distance
      return static_cast<double>(i) + frac;        // sample index (float)
    }
  } 

  return 0.0;
}
 
  template double PulseStartExtractor::extractStart<short>(const std::vector<short>& wf) const;
  template double PulseStartExtractor::extractStart<int>(const std::vector<int>& wf) const;
  template double PulseStartExtractor::extractStart<float>(const std::vector<float>& wf) const;
  template double PulseStartExtractor::extractStart<double>(const std::vector<double>& wf) const;

  template std::size_t PulseStartExtractor::extractStartSampleCF<short>(const std::vector<short>& wf, double thres) const;
  template std::size_t PulseStartExtractor::extractStartSampleCF<int>(const std::vector<int>& wf, double thres) const;
  template std::size_t PulseStartExtractor::extractStartSampleCF<float>(const std::vector<float>& wf, double thres) const;
  template std::size_t PulseStartExtractor::extractStartSampleCF<double>(const std::vector<double>& wf, double thres) const;

  template std::size_t PulseStartExtractor::findMinBin<short>(const std::vector<short>& wf, std::size_t start, std::size_t end) const;
  template std::size_t PulseStartExtractor::findMinBin<int>(const std::vector<int>& wf, std::size_t start, std::size_t end) const;
  template std::size_t PulseStartExtractor::findMinBin<float>(const std::vector<float>& wf, std::size_t start, std::size_t end) const;
  template std::size_t PulseStartExtractor::findMinBin<double>(const std::vector<double>& wf, std::size_t start, std::size_t end) const;

  template short PulseStartExtractor::computeMedian<short>(const std::vector<short>& data) const;
  template int PulseStartExtractor::computeMedian<int>(const std::vector<int>& data) const;
  template float PulseStartExtractor::computeMedian<float>(const std::vector<float>& data) const;
  template double PulseStartExtractor::computeMedian<double>(const std::vector<double>& data) const;

  template bool PulseStartExtractor::extractStartLogisticFit<short>(const std::vector<short>& wf, std::size_t guess, double par[4]) const;
  template bool PulseStartExtractor::extractStartLogisticFit<int>(const std::vector<int>& wf, std::size_t guess, double par[4]) const;
  template bool PulseStartExtractor::extractStartLogisticFit<float>(const std::vector<float>& wf, std::size_t guess, double par[4]) const;
  template bool PulseStartExtractor::extractStartLogisticFit<double>(const std::vector<double>& wf, std::size_t guess, double par[4]) const;

// ***** MODIFICATIONS ***** //

template double PulseStartExtractor::extractStartInterpolated<short>(const std::vector<short>&, double, double) const;
template double PulseStartExtractor::extractStartInterpolated<int>(const std::vector<int>&, double, double) const;
template double PulseStartExtractor::extractStartInterpolated<float>(const std::vector<float>&, double, double) const;
template double PulseStartExtractor::extractStartInterpolated<double>(const std::vector<double>&, double, double) const;



} // icarus::timing namespace
