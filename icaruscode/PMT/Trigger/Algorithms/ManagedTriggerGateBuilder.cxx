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
  , fSamplePrescale(std::max(config.SamplePrescale(), std::size_t{ 1 }))
  , fSampleOffset(config.SampleOffset())
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
void icarus::trigger::ManagedTriggerGateBuilder::doDumpConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  Base_t::doDumpConfiguration(out, indent, firstIndent);
  dumpLocalConfiguration(out, indent, indent);
  
} // icarus::trigger::ManagedTriggerGateBuilder::doDumpConfiguration()


//------------------------------------------------------------------------------
void icarus::trigger::ManagedTriggerGateBuilder::dumpLocalConfiguration(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  
  out << firstIndent << " * signal polarity: "
    << util::StandardSelectorFor<util::SignalPolarity>{}.get(fPolarity).name();

  if (fSamplePrescale > 1) {
    out << '\n' << indent
      << " * use only one out of " << fSamplePrescale << " samples";
  }
  if (fSampleOffset > 0) {
    out << '\n' << indent
      << " * skip the first " << fSampleOffset << " samples of each waveform";
  }
  
} // icarus::trigger::ManagedTriggerGateBuilder::dumpLocalConfiguration()


//------------------------------------------------------------------------------
