/**
 * @file   icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.cxx
 * @brief  Algorithm to produce trigger gates out of optical readout waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.tcc`
 *         `icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h`
 * 
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h"

// framework library
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <algorithm> // std::max()


//------------------------------------------------------------------------------
//--- ManagedTriggerGateBuilder implementation
//------------------------------------------------------------------------------
icarus::trigger::ManagedTriggerGateBuilder::ManagedTriggerGateBuilder
  (Config const& config)
  : Base_t{ config.baseConfig() }
  , fPolarity{ config.Polarity() }
{
  //
  // configuration check
  //
  switch (fPolarity) {
    case util::SignalPolarity::Positive: case util::SignalPolarity::Negative:
      break;
    default:
      throw cet::exception{ "ManagedTriggerGateBuilder" }
        << "Polarity '" << config.Polarity.valueName()
        << "' not supported (only '"
        << config.Polarity.optionName(util::SignalPolarity::Positive)
        << "' and '"
        << config.Polarity.optionName(util::SignalPolarity::Negative)
        << "' are).\n"
        ;
  } // switch
} // icarus::trigger::ManagedTriggerGateBuilder::ManagedTriggerGateBuilder()


//------------------------------------------------------------------------------
