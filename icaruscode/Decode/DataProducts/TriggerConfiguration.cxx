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
std::string const icarus::TriggerConfiguration::UnknwonGenerator { "<unknown>" };

//------------------------------------------------------------------------------
void icarus::TriggerConfiguration::dumpGateConfig
  (std::ostream& out, icarus::TriggerConfiguration::GateConfig const& gateConfig, std::string const& indent) const
{
  out << std::boolalpha
    << indent << "Gate Active:               " << gateConfig.hasGate << "\n";
  if (gateConfig.hasGate) {
    
    // gates
    out
      << indent << "Drift gate Active:         " << gateConfig.hasDriftGate << "\n"
      << indent << "MinBias Gate Active:       " << gateConfig.hasMinBiasGate << "\n"
      ;
    if (gateConfig.hasMinBiasGate)
      out << indent << "MinBias Drift Gate Active: " << gateConfig.hasMinBiasDriftGate << "\n";
    else if (gateConfig.hasMinBiasDriftGate)
      out << indent << "MinBias Drift Gate Active: " << gateConfig.hasMinBiasDriftGate << " (ignored)\n";
    
    // gate durations and offsets
    out << indent << "BeamWidth:                 " << (gateConfig.gateWidth - vetoDelay) << " ns";
    if (vetoDelay != 0)
      out << " (plus the veto of " << vetoDelay << " ns)";
    out << "\n";
    if (gateConfig.hasDriftGate)
      out << indent << "EnableWidth:               " << gateConfig.driftGateWidth << " ns\n";
    out << indent << "EarlyWarningOffset:        " << gateConfig.earlyWarningOffset << " ns\n";
    if (gateConfig.hasDriftGate)
      out << indent << "EarlyEarlyWarningOffset:   " << gateConfig.earlyEarlyWarningOffset << " ns\n";
    
    if (gateConfig.hasMinBiasGate)
      out << indent << "MinBias PreScale:          " << gateConfig.prescaleMinBias << "\n";
  }
  else {
    
    if (gateConfig.hasDriftGate)
      out << indent << "Drift gate Active:         " << gateConfig.hasDriftGate << " (ignored)\n";
    if (gateConfig.hasMinBiasGate)
      out << indent << "MinBias Gate Active:       " << gateConfig.hasMinBiasGate << " (ignored)\n";
    if (gateConfig.hasMinBiasDriftGate)
      out << indent << "MinBias Drift Gate Active: " << gateConfig.hasMinBiasDriftGate << " (ignored)\n";
    
  }
  
  out << std::noboolalpha;
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
      << "Basic trigger configuration:";
    outnl() << " Generator name: ";
    if (generator == UnknwonGenerator) out << "unknown";
    else out << "'" << generator << "'";
    outnl() << " Use WR time:    " << std::boolalpha << useWrTime << std::noboolalpha;
    outnl() << " WR time offset: " << wrTimeOffset << " ns";
    
    if (++level > verbosity) break;
    // --- verbosity: 1+ -------------------------------------------------------
    
    auto printWindowMode = [&out](unsigned int slidingWindow){
        switch (slidingWindow) {
          case 0:  out << "fixed";   break;
          case 1:  out << "sliding"; break;
          default: out << "unknown";
        } // switch
        out << " (" << slidingWindow << ")";
      };
    
    outnl() << "FPGA Configuration:";
    outnl() << " Veto Delay:             "  <<  vetoDelay << " ns";
    outnl() << " MajLevelBeamCryoEAST:   "  <<  cryoConfig[icarus::trigger::kEast].majLevelInTime;
    outnl() << " MajLevelEnableCryoEAST: "  <<  cryoConfig[icarus::trigger::kEast].majLevelDrift;
    outnl() << " WindowCryoEAST:         ";
    printWindowMode(cryoConfig[icarus::trigger::kEast].slidingWindow);
    outnl() << " MajLevelBeamCryoWEST:   "  <<  cryoConfig[icarus::trigger::kWest].majLevelInTime;
    outnl() << " MajLevelEnableCryoWEST: "  <<  cryoConfig[icarus::trigger::kWest].majLevelDrift;
    outnl() << " WindowCryoWEST:         ";
    printWindowMode(cryoConfig[icarus::trigger::kWest].slidingWindow);
    outnl() << " MajorityTriggerType:   '"  <<  majorityTriggerType << "'";
    outnl() << " RunType:               '"  <<  runType << "'";

    out << "\n";

    if (++level > verbosity) return;
    // --- verbosity: 2+ -------------------------------------------------------
    
    out << indent // absorb the newline from the previous level
      << "SPEXI Configuration:\n";
      
      out << indent << " TPCTriggerDelay: " << (tpcTriggerDelay * 0.0004) << " ms ("
        << tpcTriggerDelay << " x 400 ns)\n";
      
      out << indent << " BNB:\n";
      dumpGateConfig( out, gateConfig[icarus::trigger::kBNB], indent+"  - " );

      out << indent << " NuMI:\n";
      dumpGateConfig( out, gateConfig[icarus::trigger::kNuMI], indent+"  - " );

      out << indent << " Offbeam BNB:\n";
      dumpGateConfig( out, gateConfig[icarus::trigger::kOffBeamBNB], indent + "  - " );

      out << indent << " Offbeam NuMI:\n";
      dumpGateConfig(out, gateConfig[icarus::trigger::kOffBeamNuMI], indent + "  - " );


      out << indent << " Calibration:\n";
      dumpGateConfig(out, gateConfig[icarus::trigger::kCalibration], indent + "  - ");
      if (gateConfig[icarus::trigger::kCalibration].hasGate) {
        out << indent << "  - " << "Period:                    " << gateConfig[icarus::trigger::kCalibration].period << " ns\n";
      }

    if (++level > verbosity) break;
    // --- verbosity: 3+ -------------------------------------------------------
    
    assert(level == MaxDumpVerbosity + 1U);
    
    // this is more debug information than anything else: verbosity was too high
    outnl() << "No more information available (reached level " << level << ").";
  
  } while(false);
  
  out << "\n";
  
} // icarus::TriggerConfiguration::dump()
