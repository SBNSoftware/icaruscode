/**
 * @file   icaruscode/PMT/Data/WaveformRMS.cxx
 * @brief  A baseline RMS for a waveform.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2022
 * @see    icaruscode/PMT/Data/WaveformRMS.h
 */

// library header
#include "icaruscode/PMT/Data/WaveformRMS.h"

// C/C++ standard library
#include <ostream>


// -----------------------------------------------------------------------------
std::ostream& icarus::operator<<
  (std::ostream& out, icarus::WaveformRMS const& RMS)
  { out << RMS.RMS(); return out; }


// -----------------------------------------------------------------------------
