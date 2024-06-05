/**
 * @file   icaruscode/PMT/Trigger/Utilities/TriggerGateDataFormatting.h
 * @brief  Utilities for `TriggerGateData` printout.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 13, 2021
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEDATAFORMATTING_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEDATAFORMATTING_H


// ICARUS libraries
#include "icarusalg/Utilities/IntegerRanges.h" // icarus::makeIntegerRanges()
#include "sbnobj/ICARUS/PMT/Trigger/Data/ReadoutTriggerGate.h"

// C/C++ standard libraries
#include <ostream>


namespace icarus::trigger {
  
  namespace details {
    
    /// Container of a single gate (base class).
    template <typename Gate>
    struct GateWrapper {
      Gate const& gate;
    };
    template <typename Gate>
    GateWrapper(Gate const&) ->GateWrapper<Gate>;
    
    /// Wrapper to format a gate as "compact".
    template <typename Gate>
    struct CompactFormatter: GateWrapper<Gate> {};
    
    /// Implementation of `compactdump()` printout of a gate.
    /// @see `compactdump()`
    template <typename Gate>
    std::ostream& operator<<
      (std::ostream& out, CompactFormatter<Gate> const& wrap);
    
  } // namespace details
  
  
  /**
   * @brief Manipulator-like function for compact format of trigger gates.
   * @tparam Tick template type `Tick` for `ReadoutTriggerGate`
   * @tparam TickInterval template type `TickInterval` for `ReadoutTriggerGate`
   * @tparam ChannelIDType template type `ChannelIDType` for `ReadoutTriggerGate`
   * @param gate the gate to be formatted
   * @return a wrapper with the parameters to format the gate
   * 
   * This helper function is always used in the context of the insertion of a
   * trigger gate data object into an output stream:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * #include "icaruscode/PMT/Trigger/Utilities/TriggerGateDataFormatting.h"
   * #include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
   * #include "larcorealg/CoreUtils/enumerate.h"
   * #include <iostream>
   * 
   * // ...
   * 
   * std::vector<icarus::trigger::OpticalTriggerGate> LVDSgates;
   * 
   * // ...
   * 
   * for (auto const& [ iGate, gate ]: util::enumerate(LVDSgates))
   *   std::cout << "[" << iGate << "] " << compactdump(gate) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The loop will print all the `LVDSgates` in "compact" mode, one per line,
   * prepended by their position in the collection.
   * 
   * Technical note: `compactdump()` in the example is found by the compiler via
   * Koenig lookup.
   * 
   */
  template <typename Tick, typename TickInterval, typename ChannelIDType>
  auto compactdump
    (ReadoutTriggerGate<Tick, TickInterval, ChannelIDType> const& gate)
    -> details::CompactFormatter<ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>>
    { return { details::GateWrapper{ gate } }; }
  
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
// --  template implementation
// -----------------------------------------------------------------------------

template <typename Gate>
std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, CompactFormatter<Gate> const& wrap)
{
  auto const& gate { wrap.gate };
  
  out << "[";
  if (gate.hasChannel()) {
    out << "channel: " << gate.channel();
  }
  else if (gate.hasChannels()) {
    out << gate.nChannels() << " channels: "
      << icarus::makeIntegerRanges(gate.channels());
  }
  else out << "no channel";
  out << "] " << gate.gateLevels();
  
  if (unsigned int const maxTime = gate.findMaxOpen(); maxTime != gate.MaxTick)
  {
    unsigned int const maxOpening = gate.openingCount(maxTime);
    out << " (maximum: " << maxOpening << " at " << maxTime << ")";
  }
  
  return out;
} // icarus::trigger::operator<< (GateWrapper)


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEDATAFORMATTING_H
