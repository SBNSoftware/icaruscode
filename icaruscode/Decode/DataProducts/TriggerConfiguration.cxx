/**
 * @file   icaruscode/Decode/DataProducts/TriggerConfiguration.cxx
 * @brief  Information from the configuration of a Trigger readout board.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   March 23, 2022
 * @see    icaruscode/Decode/DataProducts/TriggerConfiguration.h
 */

// library header
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"

// C/C++ standard libraries
#include <ostream>
#include <cassert>


//------------------------------------------------------------------------------
static_assert(icarus::TriggerConfiguration::DefaultDumpVerbosity <= icarus::TriggerConfiguration::MaxDumpVerbosity);

//------------------------------------------------------------------------------
void icarus::TriggerConfiguration::dump(std::ostream& out,
  std::string const& indent, std::string const& firstIndent,
  unsigned int verbosity /* = MaxDumpVerbosity */
) const{
  
  do {
    // fake look for easy break: `goto` in disguise;
    // `break` will add a final `'\n'`, `return` will not.
    
    // start a new line with indentation:
    auto const outnl
      = [&out,&indent]() -> std::ostream& { return out << '\n' << indent; };
    
    unsigned int level = 0U;
    
    
    // --- verbosity: 0+ -------------------------------------------------------
    out << firstIndent
      << "\nBasic board configuration:\n "
      << " Use Wr time "    << useWrTime << " \n"
      << " WR time offset " << wrTimeOffset << " ns \n"; 
    
    if (++level > verbosity) break;
    // --- verbosity: 1+ -------------------------------------------------------
    outnl()
      << " FPGA Configuration: \n "
      << " Veto Delay:            "  <<  vetoDelay << " ns,\n"
      << " MajLevelBeamCryoEAST   "  <<  cryoConfig[kEast].majLevelInTime    << " \n"
      << " MajLevelEnableCryoEAST "  <<  cryoConfig[kEast].majLevelDrift  << " \n"
      << " SlidingWindowCryoEAST  "  <<  cryoConfig[kEast].slidingWindow   << " \n"
      << " MajLevelBeamCryoWEST   "  <<  cryoConfig[kWest].majLevelInTime    << " \n"
      << " MajLevelEnableCryoWEST "  <<  cryoConfig[kWest].majLevelDrift  << " \n"
      << " SlidingWindowCryoWEST  "  <<  cryoConfig[kWest].slidingWindow   << " \n"
      << " MajorityTriggerType    "  <<  majorityTriggerType                << " \n"
      << " RunType                "  <<  runType << " ";

    out << "\n";

    if (++level > verbosity) return;
    // --- verbosity: 2+ -------------------------------------------------------
    
    out << indent // absorb the newline from the previous level
      << "\nSPEXI Configuration: \n "
      << "TPCTriggerDelay"         << tpcTriggerDelay << " /400 ns \n "
      << "GateSelection "          << gateSelection << " \n "
      << "BNBBeamWidth "           << gateConfig[kBNB].inTimeWidth << " ns \n "
      << "BNBEnableWidth "         << gateConfig[kBNB].driftWidth  << " ns \n "
      << "NuMIBeamWidth "          << gateConfig[kNuMI].inTimeWidth << " ns \n "
      << "NuMIEnableWidth "        << gateConfig[kNuMI].driftWidth << " ns \n "
      << "PreScaleBNB "            << gateConfig[kBNB].prescaleMinBias << " \n "
      << "PreScaleNuMI"            << gateConfig[kNuMI].prescaleMinBias << "\n"
      << "OffBeamBNBBeamWidth "    << gateConfig[kOffBeamBNB].inTimeWidth<< " ns \n "
      << "OffBeamBNBEnableWidth "  << gateConfig[kOffBeamBNB].driftWidth << " ns \n "
      << "OffBeamNuMIBeamWidth "   << gateConfig[kOffBeamNuMI].inTimeWidth << " ns \n "
      << "OffBeamNuMIEnableWidth " << gateConfig[kOffBeamNuMI].driftWidth << " ns \n "
      << "OffBeamGateRateBNB "     << offBeamGateRate[kBNB] << " \n "
      << "OffBeamGateRateNuMI"     << offBeamGateRate[kNuMI] << " \n "
      << "PreScaleOffBeamBNB "     << gateConfig[kOffBeamBNB].prescaleMinBias << " \n "
      << "PreScaleOffBeamNuMI"     << gateConfig[kOffBeamNuMI].prescaleMinBias << " \n "
      << "ZeroBiasWidth "          << gateConfig[kCalibration].inTimeWidth << " ns \n "
      << "ZeroBiasEnableWidth "    << gateConfig[kCalibration].driftWidth << " ns \n "
      << "ZeroBiasFreq "           << period << " ns \n "
      << "PrescaleZeroBias "       << gateConfig[kCalibration].prescaleMinBias << " \n "
      << "BNBBESOffset "           << earlyWarningOffset[kBNB] << " ns \n "
      << "BNB1DOffset "            << earlyEarlyWarningOffset[kBNB] << " ns \n "
      << "NuMIMIBSOffset "         << earlyWarningOffset[kNuMI] << " ns \n "
      << "NuMIADOffset "           << earlyEarlyWarningOffset[kNuMI] << " ns \n ";

    
    if (++level > verbosity) break;
    // --- verbosity: 3+ -------------------------------------------------------
    
    assert(level == MaxDumpVerbosity + 1U);
    
    // this is more debug information than anything else: verbosity was too high
    outnl() << "No more information available (reached level " << level << ").";
  
  } while(false);
  
  out << "\n";
  
} // icarus::TriggerConfiguration::dump()