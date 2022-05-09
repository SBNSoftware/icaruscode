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
void icarus::TriggerConfiguration::dumpGateConfig
  (std::ostream& out, icarus::TriggerConfiguration::GateConfig const& gateConfig, std::string const& indent) const
{
  out << std::boolalpha
    << indent << "Gate Active:              " << gateConfig.hasGate << "\n"
    << indent << "Drift gate Active         " << gateConfig.hasDriftGate << "\n"
    << indent << "MinBias Gate Active:      " << gateConfig.hasMinBiasGate << "\n"
    << indent << "MinBias Drift Gate Active " << gateConfig.hasMinBiasDriftGate << "\n"
    << indent << "BeamWidth                 " << gateConfig.gateWidth << " ns\n"
    << indent << "EnableWidth               " << gateConfig.driftGateWidth  << " ns\n"
    << indent << "EarlyWarningOffset        " << gateConfig.earlyWarningOffset << " ns\n"
    << indent << "EarlyWarningOffset        " << gateConfig.earlyEarlyWarningOffset << " ns\n"
    << indent << "PreScale                  " << gateConfig.prescaleMinBias << "\n"
    << std::noboolalpha
    ;
}


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
      << "\nBasic trigger configuration:\n ";
      outnl() << " Use Wr time "    << useWrTime << " \n";
      outnl() << " WR time offset " << wrTimeOffset << " ns \n"; 
    
    if (++level > verbosity) break;
    // --- verbosity: 1+ -------------------------------------------------------
    outnl()
      << " FPGA Configuration: \n "
      << " Veto Delay:            "  <<  vetoDelay << " ns,\n"
      << " MajLevelBeamCryoEAST   "  <<  cryoConfig[icarus::trigger::kEast].majLevelInTime    << " \n"
      << " MajLevelEnableCryoEAST "  <<  cryoConfig[icarus::trigger::kEast].majLevelDrift  << " \n"
      << " WindowCryoEAST         "  <<  cryoConfig[icarus::trigger::kEast].slidingWindow   << " \n"
      << " MajLevelBeamCryoWEST   "  <<  cryoConfig[icarus::trigger::kWest].majLevelInTime    << " \n"
      << " MajLevelEnableCryoWEST "  <<  cryoConfig[icarus::trigger::kWest].majLevelDrift  << " \n"
      << " WindowCryoWEST         "  <<  cryoConfig[icarus::trigger::kWest].slidingWindow   << " \n"
      << " MajorityTriggerType    "  <<  majorityTriggerType                << " \n"
      << " RunType                "  <<  runType << " ";

    out << "\n";

    if (++level > verbosity) return;
    // --- verbosity: 2+ -------------------------------------------------------
    
    out << indent // absorb the newline from the previous level
      << "\nSPEXI Configuration: \n ";
      
      out << "TPCTriggerDelay" << tpcTriggerDelay << " /400 ns \n ";
      
      out << indent << " BNB:\n";
      dumpGateConfig( out, gateConfig[icarus::trigger::kBNB], indent+" - " );

      out << indent << " NuMI:\n";
      dumpGateConfig( out, gateConfig[icarus::trigger::kNuMI], indent+" - " );

      out << indent << " Offbeam BNB:\n";
      dumpGateConfig( out, gateConfig[icarus::trigger::kOffBeamBNB], indent + " - " );

      out << indent << " Offbeam NuMI:\n";
      dumpGateConfig(out, gateConfig[icarus::trigger::kOffBeamNuMI], indent + " - " );


      out << indent << " Calibration:\n";
      dumpGateConfig(out, gateConfig[icarus::trigger::kCalibration], indent + " - ");
      out << indent +" - " << "Period                    " << gateConfig[icarus::trigger::kCalibration].period << " ns\n";

    if (++level > verbosity) break;
    // --- verbosity: 3+ -------------------------------------------------------
    
    assert(level == MaxDumpVerbosity + 1U);
    
    // this is more debug information than anything else: verbosity was too high
    outnl() << "No more information available (reached level " << level << ").";
  
  } while(false);
  
  out << "\n";
  
} // icarus::TriggerConfiguration::dump()