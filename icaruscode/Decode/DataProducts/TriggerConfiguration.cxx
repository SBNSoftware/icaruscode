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
      << " Use Wr time "    << UseWrTime << " \n"
      << " WR time offset " << WrTimeOffset << " ns \n"; 
    
    if (++level > verbosity) break;
    // --- verbosity: 1+ -------------------------------------------------------
    outnl()
      << " FPGA Configuration: \n "
      << " Veto Delay:            "  << VetoDelay << " ns,\n"
      << " MajLevelBeamCryoEAST   "  <<  MajLevelBeamCryoEAST   << " \n"
      << " MajLevelEnableCryoEAST "  <<  MajLevelEnableCryoEAST    << " \n"
      << " SlidingWindowCryoEAST  "  <<  SlidingWindowCryoEAST  << " \n"
      << " MajLevelBeamCryoWEST   "  <<  MajLevelBeamCryoWEST   << " \n"
      << " MajLevelEnableCryoWEST "  <<  MajLevelEnableCryoWEST    << " \n"
      << " SlidingWindowCryoWEST  "  <<  SlidingWindowCryoWEST << " \n"
      << " MajorityTriggerType    "  <<  MajorityTriggerType << " \n"
      << " RunType                "  <<  RunType << " ";

    out << "\n";

    if (++level > verbosity) return;
    // --- verbosity: 2+ -------------------------------------------------------
    
    out << indent // absorb the newline from the previous level
      << "\nSPEXI Configuration: \n "
      << "TPCTriggerDelay"         << TPCTriggerDelay << " /400 ns \n"
      << "GateSelection "          << GateSelection << " \n "
      << "BNBBeamWidth "           << BNBBeamWidth << " ns \n "
      << "BNBEnableWidth "         << BNBEnableWidth << " ns \n "
      << "NuMIBeamWidth "          << NuMIBeamWidth << " ns \n "
      << "NuMIEnableWidth "        << NuMIEnableWidth << " ns \n "
      << "PreScaleBNBNuMI "        << PreScaleBNBNuMI << " ns \n "
      << "OffBeamBNBBeamWidth "    << OffBeamBNBBeamWidth << " ns \n "
      << "OffBeamBNBEnableWidth "  << OffBeamBNBEnableWidth << " ns \n "
      << "OffBeamNuMIBeamWidth "   << OffBeamNuMIBeamWidth << " ns \n "
      << "OffBeamNuMIEnableWidth " << OffBeamNuMIEnableWidth << " ns \n "
      << "OffBeamGateRate "        << OffBeamGateRate << " \n "
      << "PreScaleOffBeam "        << PreScaleOffBeam << " \n "
      << "ZeroBiasWidth "          << ZeroBiasWidth << " ns \n "
      << "ZeroBiasEnableWidth "    << ZeroBiasEnableWidth << " ns \n "
      << "ZeroBiasFreq "           << ZeroBiasFreq << " ns \n "
      << "PrescaleZeroBias "       << PrescaleZeroBias << " \n "
      << "BNBBESOffset "           << BNBBESOffset << " ns \n "
      << "BNB1DOffset "            << BNB1DOffset << " ns \n "
      << "NuMIMIBSOffset "         << NuMIMIBSOffset << " ns \n "
      << "NuMIADOffset "           << NuMIADOffset << " ns \n ";

    
    if (++level > verbosity) break;
    // --- verbosity: 3+ -------------------------------------------------------
    
    assert(level == MaxDumpVerbosity + 1U);
    
    // this is more debug information than anything else: verbosity was too high
    outnl() << "No more information available (reached level " << level << ").";
  
  } while(false);
  
  out << "\n";
  
} // icarus::TriggerConfiguration::dump()