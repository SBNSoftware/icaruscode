/**
 * @file   icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h
 * @brief  Simple type definitions for trigger algorithms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.cxx`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERTYPES_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERTYPES_H


// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // detinfo::timescales
#include "lardataalg/Utilities/quantities/spacetime.h" // util::quantities::microsecond
#include "lardataobj/RawData/OpDetWaveform.h" // raw::ADC_Count_t


namespace icarus::trigger {
  
  //
  // quantities
  //
  using util::quantities::microsecond;
  
  using namespace util::quantities::time_literals;
  
  //
  // times and time scales
  //
  using detinfo::timescales::optical_time;
  using detinfo::timescales::optical_tick;
  using detinfo::timescales::optical_time_ticks;
  
  
  // --- BEGIN Data type definitions -------------------------------------------
  
  /// ADC counts for optical detector readout.
  using ADCCounts_t = util::quantities::counts_as<raw::ADC_Count_t>;
  
  // --- END Data type definitions ---------------------------------------------
  
  
} // namespace icarus::trigger


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_TRIGGERTYPES_H
