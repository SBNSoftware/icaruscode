/**
 * @file   sbnobj/Common/Trigger/classes.h
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 10, 2019
 * 
 * Enables dictionary definitions for:
 * 
 * * `sbn::ExtraTriggerInfo`
 * 
 * See also `sbnobj/Common/Trigger/classes_def.xml`.
 */

// SBN libraries
// #include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/Decode/DataProducts/ExtraTriggerInfo.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"

// framework libraries
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Wrapper.h"


namespace {
  
  sbn::ExtraTriggerInfo tinfo;
  icarus::TriggerConfiguration tconfig;
  
} // local namespace
