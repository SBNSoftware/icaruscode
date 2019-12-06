/**
 * @file   icaruscode/PMT/Trigger/Data/classes.h
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 10, 2019
 * 
 * Enables dictionary definitions for:
 * 
 * * `icarus::trigger::TriggerGateData< TODO >`
 *   (and its associations with `raw::OpDetWaveform`)
 * 
 * See also `icaruscode/PMT/Trigger/Data/classes_def.xml`.
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/TriggerGateData.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

// C++ libraries
#include <ostream>

namespace {
  
  icarus::trigger::OpticalTriggerGate::GateData_t tgd;
  
} // local namespace
