/**
 * @file   icaruscode/Decode/DecoderTools/ICARUSTriggerConfigurationExtractor.cxx
 * @brief  Utility to extract the ICARUSTrigger readout configuration from data. Inspired by the PMTconfigurationExtractor authored by G. Petrillo
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   March 23, 2022
 * @see    icaruscode/Decode/DecoderTools/ICARUSTriggerConfigurationExtractor.h
 */


// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/TriggerConfigurationExtractor.h"

// C/C++ standard libraries
#include <algorithm> // std::find_if()
#include <memory> // std::unique_ptr<>


// -----------------------------------------------------------------------------
// ---  icarus::TriggerConfigurationExtractorBase
// -----------------------------------------------------------------------------
fhicl::ParameterSet
icarus::TriggerConfigurationExtractorBase::convertConfigurationDocuments(
  fhicl::ParameterSet const& container,
  std::string const& configListKey,
  std::initializer_list<std::regex const> components // TODO template with matcher functor
) {
  
  fhicl::ParameterSet const sourceConfig
    = container.get<fhicl::ParameterSet>(configListKey);
  
  fhicl::ParameterSet configDocs;
  for (auto const& key: sourceConfig.get_names()) {
    if (!sourceConfig.is_key_to_atom(key)) continue;
    
    if (!matchKey(key, components.begin(), components.end())) continue;
    
    std::string const psetStr = sourceConfig.get<std::string>(key);
    
    fhicl::ParameterSet pset;
    try {
      pset = fhicl::ParameterSet::make(psetStr);
    }
    catch (cet::exception& e) {
      throw cet::exception{ "convertConfigurationDocuments", "", e }
        << "Error parsing the content of key '" << configListKey << "." << key
        << "'; content was:\n" << psetStr << "\n";
    }
    
    configDocs.put(key, pset);
    
  } // for all main keys
  
  return configDocs;
} // convertConfigurationDocuments()



// -----------------------------------------------------------------------------
// ---  icarus::TriggerConfigurationExtractor
// -----------------------------------------------------------------------------
std::vector<std::regex> const
  icarus::TriggerConfigurationExtractor::ConfigurationNames
  { std::regex{ "icarustrigger.*" } }
  ;

// -----------------------------------------------------------------------------
bool icarus::TriggerConfigurationExtractor::isGoodConfiguration
  (fhicl::ParameterSet const& pset, std::string const& key)
{
  return matchKey(key, ConfigurationNames.begin(), ConfigurationNames.end());
} // icarus::PMTconfigurationExtractor::isGoodConfiguration()


// -----------------------------------------------------------------------------
unsigned long icarus::TriggerConfigurationExtractor::parsePrescaleString
   ( std::string prescaleString, std::size_t source ) const 
{

    /* 
    * The prescale strings represents the hex econding of a 32-bit number
    * The first, least significative 16-bit word is the prescale associated
    * to BNB, while the most significative 16-bits are associated to NuMI.
    * if the word has only 16-bits (4 character string) then it is the 
    * prescale associated to the calibration gate
    */

    using namespace std::string_literals;
    switch( source ){
      case icarus::trigger::kBNB:          return std::stoul( prescaleString.substr(4,8), nullptr, 16);   
      case icarus::trigger::kNuMI:         return std::stoul( prescaleString.substr(0,4), nullptr, 16);
      case icarus::trigger::kOffBeamBNB:   return std::stoul( prescaleString.substr(4,8), nullptr, 16); 
      case icarus::trigger::kOffBeamNuMI:  return std::stoul( prescaleString.substr(0,4), nullptr, 16);
      case icarus::trigger::kCalibration:  return std::stoul( prescaleString, nullptr ,16);
    }
    throw std::runtime_error(" icarus::TriggerConfigurationExtractor::parsePrescaleString triggerSource{"s
    + std::to_string(source) + "}): unknown value"s);

}


// -----------------------------------------------------------------------------
unsigned int icarus::TriggerConfigurationExtractor::parseWindowMode(std::string bitStr) const {
  for (unsigned int bitValue = 0; bitValue < value(sbn::bits::triggerWindowMode::NBits); ++bitValue) {
    auto const bit = static_cast<sbn::bits::triggerWindowMode>(bitValue);
    if (bitStr == bitName(bit)) return bitValue;
  }
  throw cet::exception{ "TriggerConfigurationExtractor" }
    << "Trigger window mode '" << bitStr << "' is not known.\n";
}


// -----------------------------------------------------------------------------
icarus::TriggerConfiguration icarus::TriggerConfigurationExtractor::extract
  (fhicl::ParameterSet const& pset) const
{
  
  std::optional<icarus::TriggerConfiguration> config;
  
  for (std::string const& key: pset.get_names()) {
    
    std::optional<fhicl::ParameterSet> triggerConfigPset
      = readTriggerConfig(pset, key);
    
    if (!triggerConfigPset) continue;
    
    std::optional triggerConfig
      = extractTriggerConfiguration(*triggerConfigPset);
    
    if (config) {
      throw cet::exception{ "TriggerConfigurationExtractor" }
        << "Found multiple configurations for the trigger:"
        << "\n" << std::string(80, '-') << "\n"
        << *config
        << "\n" << std::string(80, '-') << "\n"
        << "and"
        << "\n" << std::string(80, '-') << "\n"
        << *triggerConfig
        << "\n" << std::string(80, '-') << "\n"
        ;
    }
    
    config = std::move(triggerConfig);
    
  } // for
  if (!config) {
    throw cet::exception{ "TriggerConfigurationExtractor" }
      << "No trigger configuration found (fragment type: '"
      << fExpectedFragmentTypeSpec << "').\n";
  }
  
  return *config;
} // icarus::PMTconfigurationExtractor::extract()


// -----------------------------------------------------------------------------
auto icarus::TriggerConfigurationExtractor::extractTriggerConfiguration
  (fhicl::ParameterSet const& pset ) const
  -> icarus::TriggerConfiguration
{
  
  auto const& triggerParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver");

  auto const& fpgaParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver.fpga_init_params");

  auto const& spexiParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver.spexi_init_params");
 
  icarus::TriggerConfiguration rc; // readout config, for friends


  rc.useWrTime                        
    = triggerParams.get<bool>("use_wr_time");
  rc.wrTimeOffset                     
    = triggerParams.get<unsigned int>("wr_time_offset_ns");
  
  rc.vetoDelay                        
    = fpgaParams.get<unsigned int>("Veto.value");

  rc.cryoConfig[icarus::trigger::kEast].majLevelInTime 
    = fpgaParams.get<unsigned int>("MajLevelBeamCryo1.value");   
  rc.cryoConfig[icarus::trigger::kEast].majLevelDrift  
    = fpgaParams.get<unsigned int>("MajLevelEnableCryo1.value");   
  rc.cryoConfig[icarus::trigger::kEast].slidingWindow  
    = parseWindowMode( fpgaParams.get<std::string>("SlidingWindowCryo1.value") );   
  rc.cryoConfig[icarus::trigger::kWest].majLevelInTime 
    = fpgaParams.get<unsigned int>("MajLevelBeamCryo2.value");   
  rc.cryoConfig[icarus::trigger::kWest].majLevelDrift  
    = fpgaParams.get<unsigned int>("MajLevelEnableCryo2.value");   
  rc.cryoConfig[icarus::trigger::kWest].slidingWindow  
    = parseWindowMode( fpgaParams.get<std::string>("SlidingWindowCryo2.value") );   
  rc.majorityTriggerType              
    = fpgaParams.get<std::string>( "MajorityTriggerType.value");
  rc.runType                          
    = fpgaParams.get<std::string>( "RunType.value");

  rc.tpcTriggerDelay        
    = spexiParams.get<unsigned int>("TPCTriggerDelay.value");

  auto gateSelection = sbn::bits::makeMask<sbn::gateSelection>         
    (std::stoul( spexiParams.get<std::string>("GateSelection.value"), nullptr, 16));

  // Read the prescale configuraton as string for now 
  auto prescaleMinBiasBeam = 
    spexiParams.get<std::string>("PreScaleBNBNuMI.value"); 
  auto prescaleMinBiasOffBeam = 
    spexiParams.get<std::string>("PreScaleOffBeam.value"); 
  auto offBeamGateRate = 
    spexiParams.get<std::string>("OffBeamGateRate.value"); 
  auto prescaleMinBiasCalibration =
    spexiParams.get<std::string>("PrescaleZeroBias.value");

  // BNB Full Config
  rc.gateConfig[icarus::trigger::kBNB].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::GateBNB);
  rc.gateConfig[icarus::trigger::kBNB].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::DriftGateBNB);
  rc.gateConfig[icarus::trigger::kBNB].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasGateBNB);
  rc.gateConfig[icarus::trigger::kBNB].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasDriftGateBNB);
  rc.gateConfig[icarus::trigger::kBNB].gateWidth 
    = spexiParams.get<unsigned int>("BNBBeamWidth.value");
  rc.gateConfig[icarus::trigger::kBNB].driftGateWidth          
    = spexiParams.get<unsigned int>("BNBEnableWidth.value");
  rc.gateConfig[icarus::trigger::kBNB].prescaleMinBias 
    = parsePrescaleString(prescaleMinBiasBeam, icarus::trigger::kBNB);
  rc.gateConfig[icarus::trigger::kBNB].earlyWarningOffset
    = spexiParams.get<unsigned int>("BNBBESOffset.value");
  rc.gateConfig[icarus::trigger::kBNB].earlyEarlyWarningOffset
    = spexiParams.get<unsigned int>("BNB1DOffset.value");

  // OffBeamBNB config
  rc.gateConfig[icarus::trigger::kOffBeamBNB].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::GateOffbeamBNB);
  rc.gateConfig[icarus::trigger::kOffBeamBNB].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::DriftGateOffbeamBNB);
  rc.gateConfig[icarus::trigger::kOffBeamBNB].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasGateOffbeamBNB);
  rc.gateConfig[icarus::trigger::kOffBeamBNB].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasDriftGateOffbeamBNB);
  rc.gateConfig[icarus::trigger::kOffBeamBNB].gateWidth            
    = spexiParams.get<unsigned int>("OffBeamBNBBeamWidth.value");
  rc.gateConfig[icarus::trigger::kOffBeamBNB].driftGateWidth             
    = spexiParams.get<unsigned int>("OffBeamBNBEnableWidth.value");
  rc.gateConfig[icarus::trigger::kOffBeamBNB].prescaleMinBias
    = parsePrescaleString(prescaleMinBiasOffBeam, icarus::trigger::kOffBeamBNB);
  rc.gateConfig[icarus::trigger::kOffBeamBNB].offBeamGateRate
    = parsePrescaleString( offBeamGateRate, icarus::trigger::kOffBeamBNB );

  // NuMI Configuration
  rc.gateConfig[icarus::trigger::kNuMI].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::GateNuMI);
  rc.gateConfig[icarus::trigger::kNuMI].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::DriftGateNuMI);
  rc.gateConfig[icarus::trigger::kNuMI].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasGateNuMI);
  rc.gateConfig[icarus::trigger::kNuMI].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasDriftGateNuMI);
  rc.gateConfig[icarus::trigger::kNuMI].gateWidth         
    = spexiParams.get<unsigned int>("NuMIBeamWidth.value");
  rc.gateConfig[icarus::trigger::kNuMI].driftGateWidth          
    = spexiParams.get<unsigned int>("NuMIEnableWidth.value");
  rc.gateConfig[icarus::trigger::kNuMI].prescaleMinBias
    = parsePrescaleString(prescaleMinBiasBeam, icarus::trigger::kNuMI);
  rc.gateConfig[icarus::trigger::kNuMI].earlyWarningOffset
    = spexiParams.get<unsigned int>("NuMIMIBSOffset.value");
  rc.gateConfig[icarus::trigger::kNuMI].earlyEarlyWarningOffset
    = spexiParams.get<unsigned int>("NuMIADOffset.value");

  // OffbeamNuMI config 
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::GateOffbeamNuMI);
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::DriftGateOffbeamNuMI);
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasGateOffbeamNuMI);
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasDriftGateOffbeamNuMI);
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].gateWidth            
    = spexiParams.get<unsigned int>("OffBeamNuMIBeamWidth.value");
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].driftGateWidth             
    = spexiParams.get<unsigned int>("OffBeamNuMIEnableWidth.value");
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].prescaleMinBias 
    = parsePrescaleString(prescaleMinBiasOffBeam, icarus::trigger::kOffBeamNuMI);
  rc.gateConfig[icarus::trigger::kOffBeamNuMI].offBeamGateRate
    = parsePrescaleString( offBeamGateRate, icarus::trigger::kOffBeamNuMI );

  // Calibration configuration
  rc.gateConfig[icarus::trigger::kCalibration].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::GateCalibration);
  rc.gateConfig[icarus::trigger::kCalibration].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::DriftGateCalibration);
  rc.gateConfig[icarus::trigger::kCalibration].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasGateCalibration);
  rc.gateConfig[icarus::trigger::kCalibration].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::MinbiasDriftGateCalibration);
  rc.gateConfig[icarus::trigger::kCalibration].gateWidth 
    = spexiParams.get<unsigned int>("ZeroBiasWidth.value");
  rc.gateConfig[icarus::trigger::kCalibration].driftGateWidth 
    = spexiParams.get<unsigned int>("ZeroBiasEnableWidth.value");
  rc.gateConfig[icarus::trigger::kCalibration].prescaleMinBias
    = parsePrescaleString( prescaleMinBiasCalibration, icarus::trigger::kCalibration );
  rc.gateConfig[icarus::trigger::kCalibration].period 
    = spexiParams.get<unsigned int>("ZeroBiasFreq.value"); //it is actually a period 
  
  rc.generator = triggerParams.get<std::string>("generator");

  return rc;

} // icarus::TriggerConfigurationExtractor::extractTriggerConfiguration()


// -----------------------------------------------------------------------------
std::optional<fhicl::ParameterSet>
icarus::TriggerConfigurationExtractor::readTriggerConfig
  (fhicl::ParameterSet const& pset, std::string const& key) const
{
  std::optional<fhicl::ParameterSet> config;
  
  do { // fake loop for fast exit
    
    // it must be a parameter set
    if (!pset.has_key(key) || !pset.is_key_to_table(key)) break;
    
    auto boardPSet = pset.get<fhicl::ParameterSet>(key);
    
    // its "fragment_type" must be the expected one
    std::string fragmentType;
    if (!boardPSet.get_if_present("daq.fragment_receiver.generator", fragmentType))
      break;
    if (!std::regex_match(fragmentType, fExpectedFragmentType)) break;
    
    config.emplace(std::move(boardPSet)); // success
  } while (false);
  return config;
} // icarus::TriggerConfigurationExtractor::readTriggerConfig()


// -----------------------------------------------------------------------------
// ---  free functions implementation
// -----------------------------------------------------------------------------
icarus::TriggerConfiguration icarus::extractTriggerReadoutConfiguration
  (std::string const& srcFileName, icarus::TriggerConfigurationExtractor extractor)
{
  //
  // TFile::Open() call is needed to support non-local URL
  // (e.g. XRootD URL are not supported by TFile constructor).
  //
  return extractTriggerReadoutConfiguration(*(
    std::unique_ptr<TFile>{ TFile::Open(srcFileName.c_str(), "READ") }
    ),
    std::move(extractor)
    );
} // icarus::extractTriggerReadoutConfiguration(string)


// -----------------------------------------------------------------------------
icarus::TriggerConfiguration icarus::extractTriggerReadoutConfiguration
  (TFile& srcFile, icarus::TriggerConfigurationExtractor extractor)
{
  
  return details::extractTriggerReadoutConfigurationImpl
    (util::readConfigurationFromArtFile(srcFile), std::move(extractor));
  
} // icarus::extractTriggerReadoutConfiguration(TFile)


// -----------------------------------------------------------------------------
