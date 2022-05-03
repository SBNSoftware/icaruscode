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
      case icarus::kBNB:          return std::stoul( prescaleString.substr(4,8), nullptr, 16);   
      case icarus::kNuMI:         return std::stoul( prescaleString.substr(0,4), nullptr, 16);
      case icarus::kOffBeamBNB:   return std::stoul( prescaleString.substr(4,8), nullptr, 16); 
      case icarus::kOffBeamNuMI:  return std::stoul( prescaleString.substr(0,4), nullptr, 16);
      case icarus::kCalibration:  return std::stoul( prescaleString, nullptr ,16);
    }
    throw std::runtime_error(" icarus::TriggerConfigurationExtractor::parsePrescaleString triggerSource{"s
    + std::to_string(source) + "}): unknown value"s);

}



// -----------------------------------------------------------------------------
icarus::TriggerConfiguration icarus::TriggerConfigurationExtractor::extract
  (fhicl::ParameterSet const& pset) const
{
  
  icarus::TriggerConfiguration config;
  
  for (std::string const& key: pset.get_names()) {
    
    std::optional<fhicl::ParameterSet> boardConfig
      = readBoardConfig(pset, key);
    
    if (!boardConfig) continue;
    
    config = extractTriggerConfiguration(*boardConfig);
    
  } // for
  
  return config;
} // icarus::PMTconfigurationExtractor::extract()


// -----------------------------------------------------------------------------
auto icarus::TriggerConfigurationExtractor::extractTriggerConfiguration
  (fhicl::ParameterSet const& pset ) const
  -> icarus::TriggerConfiguration
{
  
  auto const& boardParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver");

  auto const& fpgaParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver.fpga_init_params");

  auto const& spexiParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver.spexi_init_params");
 
  icarus::TriggerConfiguration rc; // readout config, for friends


  rc.useWrTime                        
    = boardParams.get<bool>("use_wr_time");
  rc.wrTimeOffset                     
    = boardParams.get<unsigned int>("wr_time_offset_ns");
  
  rc.vetoDelay                        
    = fpgaParams.get<unsigned int>("Veto.value");

  rc.cryoConfig[icarus::kEast].majLevelInTime 
    = fpgaParams.get<unsigned int>("MajLevelBeamCryo1.value");   
  rc.cryoConfig[icarus::kEast].majLevelDrift  
    = fpgaParams.get<unsigned int>("MajLevelEnableCryo1.value");   
  rc.cryoConfig[icarus::kEast].slidingWindow  
    = fpgaParams.get<std::string>("SlidingWindowCryo1.value");   
  rc.cryoConfig[icarus::kWest].majLevelInTime 
    = fpgaParams.get<unsigned int>("MajLevelBeamCryo2.value");   
  rc.cryoConfig[icarus::kWest].majLevelDrift  
    = fpgaParams.get<unsigned int>("MajLevelEnableCryo2.value");   
  rc.cryoConfig[icarus::kWest].slidingWindow  
    = fpgaParams.get<std::string>("SlidingWindowCryo2.value");   
  rc.majorityTriggerType              
    = fpgaParams.get<std::string>( "MajorityTriggerType.value");
  rc.runType                          
    = fpgaParams.get<std::string>( "RunType.value");

  rc.tpcTriggerDelay        
    = spexiParams.get<unsigned int>("TPCTriggerDelay.value");

  sbn::bits::mask_t<sbn::gateSelection> gateSelection          
    = std::stoul( spexiParams.get<std::string>("GateSelection.value"), nullptr, 16);
  
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
  rc.gateConfig[icarus::kBNB].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::gateBNB);
  rc.gateConfig[icarus::kBNB].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::driftGateBNB);
  rc.gateConfig[icarus::kBNB].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasGateBNB);
  rc.gateConfig[icarus::kBNB].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasDriftGateBNB);
  rc.gateConfig[icarus::kBNB].gateWidth 
    = spexiParams.get<unsigned int>("BNBBeamWidth.value");
  rc.gateConfig[icarus::kBNB].driftGateWidth          
    = spexiParams.get<unsigned int>("BNBEnableWidth.value");
  rc.gateConfig[icarus::kBNB].prescaleMinBias 
    = parsePrescaleString(prescaleMinBiasBeam, icarus::kBNB);
  rc.gateConfig[icarus::kBNB].earlyWarningOffset
    = spexiParams.get<unsigned int>("BNBBESOffset.value");
  rc.gateConfig[icarus::kBNB].earlyEarlyWarningOffset
    = spexiParams.get<unsigned int>("BNB1DOffset.value");

  // OffBeamBNB config
  rc.gateConfig[icarus::kOffBeamBNB].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::gateOffbeamBNB);
  rc.gateConfig[icarus::kOffBeamBNB].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::driftGateOffbeamBNB);
  rc.gateConfig[icarus::kOffBeamBNB].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasGateOffbeamBNB);
  rc.gateConfig[icarus::kOffBeamBNB].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasDriftGateOffbeamBNB);
  rc.gateConfig[icarus::kOffBeamBNB].gateWidth            
    = spexiParams.get<unsigned int>("OffBeamBNBBeamWidth.value");
  rc.gateConfig[icarus::kOffBeamBNB].driftGateWidth             
    = spexiParams.get<unsigned int>("OffBeamBNBEnableWidth.value");
  rc.gateConfig[icarus::kOffBeamBNB].prescaleMinBias
    = parsePrescaleString(prescaleMinBiasOffBeam, icarus::kOffBeamBNB);
  rc.gateConfig[icarus::kOffBeamBNB].offBeamGateRate
    = parsePrescaleString( offBeamGateRate, icarus::kOffBeamBNB );

  // NuMI Configuration
  rc.gateConfig[icarus::kNuMI].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::gateNuMI);
  rc.gateConfig[icarus::kNuMI].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::driftGateNuMI);
  rc.gateConfig[icarus::kNuMI].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasGateBNB);
  rc.gateConfig[icarus::kNuMI].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasDriftGateNuMI);
  rc.gateConfig[icarus::kNuMI].gateWidth         
    = spexiParams.get<unsigned int>("NuMIBeamWidth.value");
  rc.gateConfig[icarus::kNuMI].driftGateWidth          
    = spexiParams.get<unsigned int>("NuMIEnableWidth.value");
  rc.gateConfig[icarus::kNuMI].prescaleMinBias
    = parsePrescaleString(prescaleMinBiasBeam, icarus::kNuMI);
  rc.gateConfig[icarus::kNuMI].earlyWarningOffset
    = spexiParams.get<unsigned int>("NuMIMIBSOffset.value");
  rc.gateConfig[icarus::kNuMI].earlyEarlyWarningOffset
    = spexiParams.get<unsigned int>("NuMIADOffset.value");

  // OffbeamNuMI config 
  rc.gateConfig[icarus::kOffBeamNuMI].hasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::gateOffbeamNuMI);
  rc.gateConfig[icarus::kOffBeamNuMI].hasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::driftGateOffbeamNuMI);
  rc.gateConfig[icarus::kOffBeamNuMI].hasMinBiasGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasGateOffbeamNuMI);
  rc.gateConfig[icarus::kOffBeamNuMI].hasMinBiasDriftGate
    = sbn::bits::hasBitSet(gateSelection, sbn::bits::gateSelection::minbiasDriftGateOffbeamNuMI);
  rc.gateConfig[icarus::kOffBeamNuMI].gateWidth            
    = spexiParams.get<unsigned int>("OffBeamNuMIBeamWidth.value");
  rc.gateConfig[icarus::kOffBeamNuMI].driftGateWidth             
    = spexiParams.get<unsigned int>("OffBeamNuMIEnableWidth.value");
  rc.gateConfig[icarus::kOffBeamNuMI].prescaleMinBias 
    = parsePrescaleString(prescaleMinBiasOffBeam, icarus::kOffBeamNuMI);
  rc.gateConfig[icarus::kOffBeamNuMI].offBeamGateRate
    = parsePrescaleString( offBeamGateRate, icarus::kOffBeamNuMI );

  // Calibration configuration
  rc.gateConfig[icarus::kCalibration].gateWidth 
    = spexiParams.get<unsigned int>("ZeroBiasWidth.value");
  rc.gateConfig[icarus::kCalibration].driftGateWidth 
    = spexiParams.get<unsigned int>("ZeroBiasEnableWidth.value");
  rc.gateConfig[icarus::kCalibration].prescaleMinBias
    = parsePrescaleString( prescaleMinBiasCalibration, icarus::kCalibration );
  rc.gateConfig[icarus::kCalibration].period 
    = spexiParams.get<unsigned int>("ZeroBiasFreq.value"); //it is actually a period 

  return rc;

} // icarus::TriggerConfigurationExtractor::extractTriggerConfiguration()


// -----------------------------------------------------------------------------
std::optional<fhicl::ParameterSet>
icarus::TriggerConfigurationExtractor::readBoardConfig
  (fhicl::ParameterSet const& pset, std::string const& key) const
{
  static std::string const ExpectedFragmentType = "ICARUSTriggerUDP";
  
  std::optional<fhicl::ParameterSet> config;
  
  do { // fake loop for fast exit
    
    // it must be a parameter set
    if (!pset.has_key(key) || !pset.is_key_to_table(key)) break;
    
    auto boardPSet = pset.get<fhicl::ParameterSet>(key);
    
    // its "fragment_type" must be the expected one
    std::string fragmentType;
    if (!boardPSet.get_if_present("daq.fragment_receiver.generator", fragmentType))
      break;
    if (fragmentType != ExpectedFragmentType) break;
    
    config.emplace(std::move(boardPSet)); // success
  } while (false);
  return config;
} // icarus::TriggerConfigurationExtractor::readBoardConfig()


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

