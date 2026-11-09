/**
 * @file icaruscode/Timing/Tools/PulseStartExtractor.cc
 * @brief `icarus::timing::PulseStartExtractor` source.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date  Sun Feb 11 11:37:14 2024
 * Description: Source file for the PulseStartExtractor class implementation. It provides
 *              templated functions to extract the start (rising edge) time of reference
 *              signals (square pulses) using either constant fraction discrimination,
 *              an interpolated constant fraction or a refined logistic function fit.
 **/

#include "PulseStartExtractor.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TF1.h"
#include "TGraph.h"
#include "TFitResultPtr.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <vector>

namespace {

  // median of the [begin, end) range (the range is copied)
  double median(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end)
  {
    std::vector<double> copy(begin, end);
    auto const mid = copy.begin() + copy.size() / 2;
    std::nth_element(copy.begin(), mid, copy.end());
    return *mid;
  }

} // local namespace

namespace icarus::timing {

  // Constructor
  PulseStartExtractor::PulseStartExtractor(ExtractionMethod method, double threshold)
      : fMethod(method), fThreshold(threshold) {}

  // --------------------------------------------------------------------------------
  template <typename T>
  double PulseStartExtractor::extractStart(const std::vector<T> &wf) const
  {
    std::vector<double> const y(wf.begin(), wf.end());

    // if no signal is found above threshold, first sample (0) is returned
    Pulse pulse;
    if (!findPulse(y, pulse))
      return 0.;

    switch (fMethod)
    {
      case INTERPOLATED_CF:
        return crossingTime(y, pulse, kCFFraction, true);

      case LOGISTIC_FIT:
      {
        double start = 0.;
        double const guess = crossingTime(y, pulse, kCFFraction, false);
        if (fitLogistic(y, pulse, guess, start))
          return start; // Rising edge time from the logistic fit (can be non-integer).

        mf::LogDebug("PulseStartExtractor") << "Logistic fit did not converge;"
                                            << " falling back to the interpolated 50% crossing.";
        return crossingTime(y, pulse, 0.5, true);
      }

      case CONSTANT_FRACTION:
      default:
        return crossingTime(y, pulse, kCFFraction, false);
    }
  }

  // --------------------------------------------------------------------------------
  bool PulseStartExtractor::findPulse(std::vector<double> const &wf, Pulse &pulse) const
  {
    if (wf.empty())
      return false;

    // Find the position of the minimum along the waveform
    pulse.minBin = std::distance(wf.begin(), std::min_element(wf.begin(), wf.end()));

    // Search only a cropped region of the waveform backward from the min
    pulse.searchStart = (pulse.minBin >= kLookBackSamples) ? (pulse.minBin - kLookBackSamples) : 0;

    // baseline from everything before the search region (at least the first sample)
    std::size_t const baselineEnd = std::max<std::size_t>(pulse.searchStart, 1);
    pulse.baseline = median(wf.begin(), wf.begin() + baselineEnd);
    pulse.amplitude = pulse.baseline - wf[pulse.minBin];

    // check if peak is high enough (remove noise)
    return pulse.amplitude >= fThreshold;
  }

  // --------------------------------------------------------------------------------
  double PulseStartExtractor::crossingTime(std::vector<double> const &wf, Pulse const &pulse,
                                           double fraction, bool interpolate) const
  {
    double const level = pulse.baseline - fraction * pulse.amplitude;

    // first sample at or below the level; the minimum always is
    std::size_t bin = pulse.searchStart;
    while (bin < pulse.minBin && wf[bin] > level)
      ++bin;

    // already below the level at the start of the search window: nothing before to use
    if (bin == pulse.searchStart)
      return static_cast<double>(bin);

    if (!interpolate)
      return static_cast<double>(bin - 1);

    // interpolate between the samples across the level
    double const before = wf[bin - 1];
    double const after = wf[bin];
    return static_cast<double>(bin - 1) + (before - level) / (before - after);
  }

  // --------------------------------------------------------------------------------
  bool PulseStartExtractor::fitLogistic(std::vector<double> const &wf, Pulse const &pulse,
                                        double guess, double &start) const
  {
    std::size_t const nsize = wf.size();
    std::size_t const guessBin = static_cast<std::size_t>(guess);

    // amplitude: level after the signal (avg over kFitAmplitudeSamples samples)
    std::size_t const ampEnd = std::min(guessBin + kFitAmplitudeSamples, nsize);
    double const amplitude = std::accumulate(wf.begin() + guessBin, wf.begin() + ampEnd, 0.)
                             / (ampEnd - guessBin) - pulse.baseline;

    // crop a window around the initial guess
    std::size_t const llim = (guessBin >= kFitHalfWindow) ? guessBin - kFitHalfWindow : 0;
    std::size_t const ulim = std::min(guessBin + kFitHalfWindow, nsize - 1);

    std::vector<double> x(ulim - llim + 1);
    std::iota(x.begin(), x.end(), static_cast<double>(llim));
    TGraph graph(x.size(), x.data(), wf.data() + llim);

    // local function, not registered in the ROOT global list
    auto const logistic = [](double const *xx, double const *p)
      { return p[2] / (1. + std::exp((p[0] - xx[0]) / p[1])) + p[3]; };
    TF1 fitFunc("logistic", logistic, x.front(), x.back(), 4, 1, TF1::EAddToList::kNo);

    // try a fixed set of initial values until the fit is acceptable
    for (double const offset : kFitStartOffsets)
    {
      for (double const slope : kFitSlopeSeeds)
      {
        fitFunc.SetParameters(guess + offset, slope, amplitude, pulse.baseline);
        TFitResultPtr const fp = graph.Fit(&fitFunc, "SRQN");
        if (int(fp) != 0)
          continue;

        // p0 : rising edge as the inflection point (50% of amplitude)
        double const p0 = fitFunc.GetParameter(0);
        bool const inWindow = (p0 >= x.front()) && (p0 <= x.back());
        bool const negativePulse = (fitFunc.GetParameter(2) < 0.);
        if (inWindow && negativePulse)
        {
          start = p0;
          return true;
        }
      }
    }

    return false;
  }

  // --------------------------------------------------------------------------------
  template double PulseStartExtractor::extractStart<short>(const std::vector<short>& wf) const;
  template double PulseStartExtractor::extractStart<int>(const std::vector<int>& wf) const;
  template double PulseStartExtractor::extractStart<float>(const std::vector<float>& wf) const;
  template double PulseStartExtractor::extractStart<double>(const std::vector<double>& wf) const;

} // icarus::timing namespace
