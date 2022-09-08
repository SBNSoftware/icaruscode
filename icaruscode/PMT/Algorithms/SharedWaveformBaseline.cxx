/**
 * @file   icaruscode/PMT/Algorithms/SharedWaveformBaseline.cxx
 * @brief  Extracts and writes PMT waveform baselines (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 5, 2022
 * @see    icaruscode/PMT/Algorithms/SharedWaveformBaseline.h
 */

// library header
#include "icaruscode/PMT/Algorithms/SharedWaveformBaseline.h"

// LArSoft libraries
#include "lardataalg/Utilities/StatCollector.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::nth_element(), std::copy()
#include <iterator> // std::distance(), std::next()
#include <ostream>
#include <cmath> // std::round()
#include <cassert>


// -----------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Extracts the median of the collection between the specified iterators.
  template <typename BIter, typename EIter>
  auto median(BIter begin, EIter end) {
    
    using value_type = typename BIter::value_type;
    
    std::vector<value_type> data{ begin, end };
    assert(!data.empty());
    
    auto const middle = data.begin() + data.size() / 2;
    std::nth_element(data.begin(), middle, data.end());
    
    return *middle;
    
  } // median()
  
  
  /**
   * Returns iterator to the first sample outside `lower`-`upper` range
   * (inclusive) which is the first of at least `maxLower` samples all below
   * lower, or of at least `maxUpper` samples all above `upper`.
   * If all samples are in that range, `end` is returned.
   */
  template <typename Iter>
  Iter findOutOfBoundary(
    Iter begin, Iter end,
    typename Iter::value_type lower, typename Iter::value_type upper,
    unsigned int maxLower, unsigned int maxUpper
  ) {
    assert(lower <= upper);
    assert(maxLower > 0U);
    assert(maxUpper > 0U);

    unsigned int nAbove = 0U;
    unsigned int nBelow = 0U;
    for (auto it = begin; it != end; ++it) {
      
      if (*it > upper) {
        if (++nAbove >= maxUpper) return it - maxUpper + 1;
        nBelow = 0U;
        continue;
      }
      
      if (*it < lower) {
        if (++nBelow >= maxLower) return it - maxLower + 1;
        nAbove = 0U;
        continue;
      }
      
      nAbove = nBelow = 0U;
      
    } // while
    
    return end;
  } // findOutOfBoundary()
  
  
  // ---------------------------------------------------------------------------
  /// Helper for streaming waveforms (`std::cout << waveformIntro(waveform)`).
  struct waveformIntro {
    raw::OpDetWaveform const* waveform;
    waveformIntro(raw::OpDetWaveform const* waveform): waveform(waveform) {}
    waveformIntro(raw::OpDetWaveform const& waveform): waveform(&waveform) {}
  };
  
  
  /// Prints an "introduction" to the specified waveform.
  std::ostream& operator<< (std::ostream& out, waveformIntro wf) {
    if (wf.waveform) {
      raw::OpDetWaveform const& waveform = *(wf.waveform);
      out << "waveform channel="
        << waveform.ChannelNumber() << " timestamp=" << waveform.TimeStamp()
        << " size=" << waveform.size();
    }
    else out << "<no waveform>";
    return out;
  } // operator<< (waveformIntro)
  
  // ---------------------------------------------------------------------------
  
} // local namespace


//------------------------------------------------------------------------------
//---  opdet::SharedWaveformBaseline
//------------------------------------------------------------------------------
auto opdet::SharedWaveformBaseline::operator()
  (std::vector<raw::OpDetWaveform const*> const& waveforms) const
  -> BaselineInfo_t
{
  if (waveforms.empty()) return {};
  
  //
  // first pass: find statistics
  //
  std::vector<raw::ADC_Count_t> samples;
  samples.reserve(fParams.nSample * waveforms.size());
  std::vector<double> RMSs;
  RMSs.reserve(waveforms.size());
  
  for (raw::OpDetWaveform const* waveform: waveforms) {
    
    mf::LogTrace{ fLogCategory } << "Now processing: " << waveformIntro(waveform);
    
    if (waveform->size() < fParams.nSample) {
      mf::LogTrace{ fLogCategory } << waveformIntro(waveform)
        << ": skipped because shorter than " << fParams.nSample
        << " samples";
      continue;
    }
    
    auto const begin = waveform->cbegin();
    auto const end = std::next(begin, fParams.nSample);
    
    lar::util::StatCollector<double> stats;
    for (auto it = begin; it != end; ++it) stats.add(*it);
    RMSs.push_back(stats.RMS());
    
    std::copy(waveform->begin(), waveform->end(), back_inserter(samples));
  } // for
  
  double const medRMS = median(RMSs.cbegin(), RMSs.cend());
  raw::ADC_Count_t const med = median(samples.cbegin(), samples.cend());
  
  mf::LogTrace{ fLogCategory } << "Stats of channel "
    << waveforms.front()->ChannelNumber() << " from "
    << fParams.nSample << " starting samples of " << waveforms.size()
    << " waveforms: median=" << med << " ADC, median RMS of each waveform="
    << medRMS << " ADC";
  
  //
  // collect the samples
  //
  raw::ADC_Count_t const aboveThreshold
    = med + static_cast<raw::ADC_Count_t>(std::round(medRMS * fParams.nRMS));
  raw::ADC_Count_t const belowThreshold
    = med - static_cast<raw::ADC_Count_t>(std::round(medRMS * fParams.nRMS));
  
  auto const sampleOutOfBoundary
    = [this, aboveThreshold, belowThreshold](auto begin, auto end)
    {
      return findOutOfBoundary(begin, end,
                               belowThreshold, aboveThreshold,
                               fParams.nExcessSamples, fParams.nExcessSamples);
    };
  
  lar::util::StatCollector<double> stats;
  unsigned int nUsedWaveforms = 0U;
  for (raw::OpDetWaveform const* waveform: waveforms) {
    
    auto const begin = waveform->cbegin();
    auto const end = std::next(begin, fParams.nSample);
    
    //
    // check whether to use this waveform
    //
    auto const firstExcess = sampleOutOfBoundary(begin, end);
    if (firstExcess != end) {
      
      mf::LogTrace log { fLogCategory };
      log
        << waveformIntro(waveform) << " has " << fParams.nExcessSamples
        << " samples in a row out of [ " << belowThreshold << " ; "
        << aboveThreshold << " ] ADC starting at sample #"
        << (firstExcess - begin) << ":";
      for (
        auto it = firstExcess; it != firstExcess + fParams.nExcessSamples; ++it
      )
        log << " " << *it;
      
      // should we try to recover part of the waveform here? e.g.
      /*
      if (((firstExcess - begin) < fParams.nSample*3/5) || (waveforms.size() >= 5))
        continue;
      end = begin + fParams.nSample / 2;
       */
      
    } // if
    
    //
    // include it
    //
    ++nUsedWaveforms;
    stats.add_unweighted(begin, end);
    
  } // for
  
  return {
      stats.Average()           // baseline
    , medRMS                    // RMS
    , nUsedWaveforms            // nWaveforms
    , (unsigned int) stats.N()  // nSamples
    };
  
} // opdet::SharedWaveformBaseline::operator()


//------------------------------------------------------------------------------
