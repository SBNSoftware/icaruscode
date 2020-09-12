/**
 * @file   icaruscode/PMT/Data/classes.h
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 11, 2020
 * 
 * Enables dictionary definitions for:
 * 
 * * `icarus::WaveformBaseline`
 *   (and its associations with `raw::OpDetWaveform`)
 * 
 * See also `icaruscode/PMT/Data/classes_def.xml`.
 */

// ICARUS libraries
#include "icaruscode/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

// C++ libraries
#include <ostream>

namespace {
  
  icarus::WaveformBaseline baseline;
  
} // local namespace
