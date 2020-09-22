/**
 * @file   icaruscode/PMT/Data/WaveformBaseline.cxx
 * @brief  A baseline for a waveform.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 11, 2020
 * @see    icaruscode/PMT/Data/WaveformBaseline.h
 */

// library header
#include "icaruscode/PMT/Data/WaveformBaseline.h"

// C/C++ standard library
#include <ostream>


// -----------------------------------------------------------------------------
std::ostream& icarus::operator<<
  (std::ostream& out, icarus::WaveformBaseline const& baseline)
  { out << baseline.baseline(); return out; }


// -----------------------------------------------------------------------------
