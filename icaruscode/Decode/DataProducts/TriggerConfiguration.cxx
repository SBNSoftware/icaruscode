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
      << " MajLevelBeamCryoEAST   "  <<  cryoConfig[icarus::kEast].majLevelInTime    << " \n"
      << " MajLevelEnableCryoEAST "  <<  cryoConfig[icarus::kEast].majLevelDrift  << " \n"
      << " SlidingWindowCryoEAST  "  <<  cryoConfig[icarus::kEast].slidingWindow   << " \n"
      << " MajLevelBeamCryoWEST   "  <<  cryoConfig[icarus::kWest].majLevelInTime    << " \n"
      << " MajLevelEnableCryoWEST "  <<  cryoConfig[icarus::kWest].majLevelDrift  << " \n"
      << " SlidingWindowCryoWEST  "  <<  cryoConfig[icarus::kWest].slidingWindow   << " \n"
      << " MajorityTriggerType    "  <<  majorityTriggerType                << " \n"
      << " RunType                "  <<  runType << " ";

    out << "\n";

    if (++level > verbosity) return;
    // --- verbosity: 2+ -------------------------------------------------------
    
    out << indent // absorb the newline from the previous level
      << "\nSPEXI Configuration: \n "
      << "TPCTriggerDelay"                        << tpcTriggerDelay << " /400 ns \n "
      << "BNB Gate Active: "                      << gateConfig[icarus::kBNB].hasGate << "\n"
      << "BNB Drift gate Active"                  << gateConfig[icarus::kBNB].hasDriftGate << "\n"
      << "BNB MinBias Gate Active: "              << gateConfig[icarus::kBNB].hasMinBiasGate << "\n"
      << "BNB MinBias Drift Gate Active"          << gateConfig[icarus::kBNB].hasMinBiasDriftGate << "\n"
      << "BNBBeamWidth "                          << gateConfig[icarus::kBNB].gateWidth << " ns \n "
      << "BNBEnableWidth "                        << gateConfig[icarus::kBNB].driftGateWidth  << " ns \n "
      << "BNBBESOffset "                          << gateConfig[icarus::kBNB].earlyWarningOffset << " ns \n "
      << "BNB1DOffset "                           << gateConfig[icarus::kBNB].earlyEarlyWarningOffset << " ns \n "
      << "PreScaleBNB "                           << gateConfig[icarus::kBNB].prescaleMinBias << " \n "
      << "NuMI Gate Active: "                     << gateConfig[icarus::kNuMI].hasGate << "\n"
      << "NuMI Drift gate Active"                 << gateConfig[icarus::kNuMI].hasDriftGate << "\n"
      << "NuMI MinBias Gate Active: "             << gateConfig[icarus::kNuMI].hasMinBiasGate << "\n"
      << "NuMI MinBias Drift Gate Active"         << gateConfig[icarus::kNuMI].hasMinBiasDriftGate << "\n"
      << "NuMIBeamWidth "                         << gateConfig[icarus::kNuMI].gateWidth << " ns \n "
      << "NuMIEnableWidth "                       << gateConfig[icarus::kNuMI].driftGateWidth << " ns \n "
      << "NuMIMIBSOffset "                        << gateConfig[icarus::kNuMI].earlyWarningOffset << " ns \n "
      << "NuMIADOffset "                          << gateConfig[icarus::kNuMI].earlyEarlyWarningOffset << " ns \n "
      << "PreScaleNuMI "                          << gateConfig[icarus::kNuMI].prescaleMinBias << "\n"
      << "OffBeamBNB Gate Active: "               << gateConfig[icarus::kOffBeamBNB].hasGate << "\n"
      << "OffBeamBNB Drift gate Active"           << gateConfig[icarus::kOffBeamBNB].hasDriftGate << "\n"
      << "OffBeamBNB MinBias Gate Active: "       << gateConfig[icarus::kOffBeamBNB].hasMinBiasGate << "\n"
      << "OffBeamBNB MinBias Drift Gate Active"   << gateConfig[icarus::kOffBeamBNB].hasMinBiasDriftGate << "\n"
      << "OffBeamBNBBeamWidth "                   << gateConfig[icarus::kOffBeamBNB].gateWidth<< " ns \n "
      << "OffBeamBNBEnableWidth "                 << gateConfig[icarus::kOffBeamBNB].driftGateWidth << " ns \n "
      << "OffBeamGateRateBNB "                    << gateConfig[icarus::kOffBeamBNB].offBeamGateRate << " \n "
      << "PreScaleOffBeamBNB "                    << gateConfig[icarus::kOffBeamBNB].prescaleMinBias << " \n "
      << "OffBeamNuMI Gate Active: "              << gateConfig[icarus::kOffBeamNuMI].hasGate << "\n"
      << "OffBeamNuMI Drift gate Active"          << gateConfig[icarus::kOffBeamNuMI].hasDriftGate << "\n"
      << "OffBeamNuMI MinBias Gate Active: "      << gateConfig[icarus::kOffBeamNuMI].hasMinBiasGate << "\n"
      << "OffBeamNuMI MinBias Drift Gate Active"  << gateConfig[icarus::kOffBeamNuMI].hasMinBiasDriftGate << "\n"
      << "OffBeamNuMIBeamWidth "                  << gateConfig[icarus::kOffBeamNuMI].gateWidth << " ns \n "
      << "OffBeamNuMIEnableWidth "                << gateConfig[icarus::kOffBeamNuMI].driftGateWidth << " ns \n "
      << "OffBeamGateRateNuMI"                    << gateConfig[icarus::kOffBeamNuMI].offBeamGateRate << " \n "
      << "PreScaleOffBeamNuMI"                    << gateConfig[icarus::kOffBeamNuMI].prescaleMinBias << " \n "
      << "Calibration Gate Active: "              << gateConfig[icarus::kCalibration].hasGate << "\n"
      << "Calibration Drift gate Active"          << gateConfig[icarus::kCalibration].hasDriftGate << "\n"
      << "Calibration MinBias Gate Active: "      << gateConfig[icarus::kCalibration].hasMinBiasGate << "\n"
      << "Calibration MinBias Drift Gate Active"  << gateConfig[icarus::kCalibration].hasMinBiasDriftGate << "\n"
      << "Calibration Width "                     << gateConfig[icarus::kCalibration].gateWidth << " ns \n "
      << "Calibration EnableWidth "               << gateConfig[icarus::kCalibration].driftGateWidth << " ns \n "
      << "Calibration Frequency "                 << gateConfig[icarus::kCalibration].period << " ns \n "
      << "Prescale Calibration "                  << gateConfig[icarus::kCalibration].prescaleMinBias << " \n ";



    
    if (++level > verbosity) break;
    // --- verbosity: 3+ -------------------------------------------------------
    
    assert(level == MaxDumpVerbosity + 1U);
    
    // this is more debug information than anything else: verbosity was too high
    outnl() << "No more information available (reached level " << level << ").";
  
  } while(false);
  
  out << "\n";
  
} // icarus::TriggerConfiguration::dump()