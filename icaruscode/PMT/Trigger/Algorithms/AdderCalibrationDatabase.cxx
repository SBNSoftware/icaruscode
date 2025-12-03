/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.cxx
 * @brief  Calibration utilities for adder board output and simulation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 5, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.h
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.h"


// -----------------------------------------------------------------------------
icarus::trigger::AdderCalibrationDatabase::RunCalibration::UnknownChannelError::UnknownChannelError
  (icarus::trigger::AdderChannelID channel, std::string const& msg /* = "" */)
  : channel{ channel }
{
  if (msg.empty()) *this << "Unknown adder channel: " << channel;
  else *this << msg;
}


// -----------------------------------------------------------------------------
void icarus::trigger::AdderCalibrationDatabase::doDumpConfig
  (std::ostream& out, details::Indenter& nextLine) const
{
  // no not end the last line.
  out << nextLine << "AdderCalibrationDatabase: no configuration needed.";
}


// -----------------------------------------------------------------------------
