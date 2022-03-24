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


  rc.UseWrTime              = boardParams.get<bool>("use_wr_time");
  rc.WrTimeOffset           = boardParams.get<unsigned int>("wr_time_offset_ns");
  
  rc.VetoDelay              = fpgaParams.get<unsigned int>("Veto.value");
  rc.MajLevelBeamCryoEAST   = fpgaParams.get<unsigned int>("MajLevelBeamCryo1.value");   
  rc.MajLevelEnableCryoEAST = fpgaParams.get<unsigned int>("MajLevelEnableCryo1.value");   
  rc.SlidingWindowCryoEAST  = fpgaParams.get<std::string>("SlidingWindowCryo1.value");   
  rc.MajLevelBeamCryoWEST   = fpgaParams.get<unsigned int>("MajLevelBeamCryo2.value");   
  rc.MajLevelEnableCryoWEST = fpgaParams.get<unsigned int>("MajLevelEnableCryo2.value");   
  rc.SlidingWindowCryoWEST  = fpgaParams.get<std::string>("SlidingWindowCryo2.value");   
  rc.MajorityTriggerType    = fpgaParams.get<std::string>( "MajorityTriggerType.value");
  rc.RunType                = fpgaParams.get<std::string>( "RunType.value");

  rc.TPCTriggerDelay        = spexiParams.get<unsigned int>("TPCTriggerDelay.value");
  rc.GateSelection          = spexiParams.get<std::string>("GateSelection.value");
  rc.BNBBeamWidth           = spexiParams.get<unsigned int>("BNBBeamWidth.value");
  rc.BNBEnableWidth         = spexiParams.get<unsigned int>("BNBEnableWidth.value");
  rc.NuMIBeamWidth          = spexiParams.get<unsigned int>("NuMIBeamWidth.value");
  rc.NuMIEnableWidth        = spexiParams.get<unsigned int>("NuMIEnableWidth.value");
  rc.PreScaleBNBNuMI        = spexiParams.get<std::string>("PreScaleBNBNuMI.value");
  rc.OffBeamBNBBeamWidth    = spexiParams.get<unsigned int>("OffBeamBNBBeamWidth.value");
  rc.OffBeamBNBEnableWidth  = spexiParams.get<unsigned int>("OffBeamBNBEnableWidth.value");
  rc.OffBeamNuMIBeamWidth   = spexiParams.get<unsigned int>("OffBeamNuMIBeamWidth.value");
  rc.OffBeamNuMIEnableWidth = spexiParams.get<unsigned int>("OffBeamNuMIEnableWidth.value");
  rc.OffBeamGateRate        = spexiParams.get<std::string>("OffBeamGateRate.value");
  rc.PreScaleOffBeam        = spexiParams.get<std::string>("PreScaleOffBeam.value");
  rc.ZeroBiasWidth          = spexiParams.get<unsigned int>("ZeroBiasWidth.value");
  rc.ZeroBiasEnableWidth    = spexiParams.get<unsigned int>("ZeroBiasEnableWidth.value");
  rc.ZeroBiasFreq           = spexiParams.get<unsigned int>("ZeroBiasFreq.value");
  rc.PrescaleZeroBias       = spexiParams.get<std::string>("PrescaleZeroBias.value");
  rc.BNBBESOffset           = spexiParams.get<unsigned int>("BNBBESOffset.value");
  rc.BNB1DOffset            = spexiParams.get<unsigned int>("BNB1DOffset.value");
  rc.NuMIMIBSOffset         = spexiParams.get<unsigned int>("NuMIMIBSOffset.value");
  rc.NuMIADOffset           = spexiParams.get<unsigned int>("NuMIADOffset.value") ;    
      
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

