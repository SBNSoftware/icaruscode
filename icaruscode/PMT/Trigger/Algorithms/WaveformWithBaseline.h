/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WaveformWithBaseline.h
 * @brief  Simple structure to carry around a waveform and its baseline.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WAVEFORMWITHBASELINE_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WAVEFORMWITHBASELINE_H

// SBN libraries
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <tuple>


//------------------------------------------------------------------------------
namespace icarus::trigger { struct WaveformWithBaseline; }
/// Simple object to carry around waveform and its baseline.
struct icarus::trigger::WaveformWithBaseline
  : std::tuple<raw::OpDetWaveform const*, icarus::WaveformBaseline const*>
{
  
  using Waveform_t = raw::OpDetWaveform;
  using Baseline_t = icarus::WaveformBaseline;
  using Base_t = std::tuple<Waveform_t const*, Baseline_t const*>;
  
  using Base_t::Base_t; // inherit constructors
  
  /// Returns a reference to the waveform.
  Waveform_t const& waveform() const { return *waveformPtr(); }
  
  /// Returns a reference to the waveform baseline.
  Baseline_t const& baseline() const { return *baselinePtr(); }
  
  /// Returns a pointer to the waveform.
  Waveform_t const* waveformPtr() const
    { return std::get<Waveform_t const*>(*this); }
  
  /// Returns a pointer to the waveform baseline.
  Baseline_t const* baselinePtr() const
    { return std::get<Baseline_t const*>(*this); }
  
  /// Returns whether the baseline is available for this waveform.
  bool hasBaseline() const { return baselinePtr() != nullptr; }
  
  
  // questionable practices...
  operator Waveform_t const& () const { return waveform(); }
  operator Waveform_t const* () const { return waveformPtr(); }
  operator Baseline_t const& () const { return baseline(); }
  operator Baseline_t const* () const { return baselinePtr(); }
  
}; // icarus::trigger::WaveformWithBaseline


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WAVEFORMWITHBASELINE_H
