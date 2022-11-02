/**
 * @file   icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.cxx
 * @brief  Utility to extract PMT readout configuration from data.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 18, 2021
 * @see    icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.h
 */


// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.h"

// C/C++ standard libraries
#include <algorithm> // std::find_if()
#include <memory> // std::unique_ptr<>



// -----------------------------------------------------------------------------
// ---  icarus::PMTconfigurationExtractorBase
// -----------------------------------------------------------------------------
fhicl::ParameterSet
icarus::PMTconfigurationExtractorBase::convertConfigurationDocuments(
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
// ---  icarus::PMTconfigurationExtractor
// -----------------------------------------------------------------------------
std::vector<std::regex> const
  icarus::PMTconfigurationExtractor::ConfigurationNames
  { std::regex{ "icaruspmt.*" } }
  ;

// -----------------------------------------------------------------------------
bool icarus::PMTconfigurationExtractor::isGoodConfiguration
  (fhicl::ParameterSet const& pset, std::string const& key)
{
  return matchKey(key, ConfigurationNames.begin(), ConfigurationNames.end());
} // icarus::PMTconfigurationExtractor::isGoodConfiguration()


// -----------------------------------------------------------------------------
sbn::PMTconfiguration icarus::PMTconfigurationExtractor::extract
  (fhicl::ParameterSet const& pset) const
{
  
  sbn::PMTconfiguration config;
  
  for (std::string const& key: pset.get_names()) {
    
    std::optional<fhicl::ParameterSet> boardConfig
      = readBoardConfig(pset, key);
    
    if (!boardConfig) continue;
    
    config.boards.push_back(extractV1730configuration(*boardConfig, key));
    
  } // for
  
  return config;
} // icarus::PMTconfigurationExtractor::extract()


// -----------------------------------------------------------------------------
auto icarus::PMTconfigurationExtractor::extractV1730configuration
  (fhicl::ParameterSet const& pset, std::string const& boardName) const
  -> sbn::V1730Configuration
{
  
  auto const& boardParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver");
  
  sbn::V1730Configuration rc; // readout config, for friends
  rc.boardName = boardName;
  
  rc.boardID = boardParams.get<unsigned int>("board_id");
  rc.fragmentID = boardParams.get<unsigned int>("fragment_id");
  
  rc.bufferLength = boardParams.get<int>("recordLength");
  rc.postTriggerFrac = boardParams.get<float>("postPercent") / 100.0f;
  
  rc.nChannels = boardParams.get<unsigned int>("nChannels");
  
  rc.useTimeTagForTimeStamp = boardParams.get("UseTimeTagForTimeStamp", false);
  
  for (unsigned short int channelNo = 0; channelNo < 16; ++channelNo)
    rc.channels.push_back(extractChannelConfiguration(boardParams, channelNo));
  
  return rc;
} // icarus::PMTconfigurationExtractor::extractV1730configuration()


// -----------------------------------------------------------------------------
sbn::V1730channelConfiguration
icarus::PMTconfigurationExtractor::extractChannelConfiguration
  (fhicl::ParameterSet const& boardPSet, unsigned short int channelNo) const
{
  std::string const ChannelStr = std::to_string(channelNo);
  
  sbn::V1730channelConfiguration channel;
  
  channel.channelNo = channelNo;
  
  channel.baseline = boardPSet.get<short signed int>
    ("BaselineCh" + std::to_string(channelNo + 1));
  channel.threshold = boardPSet.get<short signed int>
    ("triggerThreshold" + ChannelStr);
  channel.enabled = boardPSet.get<bool>("channelEnable" + ChannelStr);
  
  return channel;
} // icarus::PMTconfigurationExtractor::extractChannelConfiguration()


// ---------------------------------------------------------------------------
sbn::PMTconfiguration icarus::PMTconfigurationExtractor::finalize
  (sbn::PMTconfiguration config) const
{

  for (sbn::V1730Configuration& readoutBoardConfig: config.boards) {
    auto const fragmentID = readoutBoardDBfragmentID(readoutBoardConfig);
    if (!fChannelMap->hasPMTDigitizerID(fragmentID)) {
      mf::LogWarning("PMTconfigurationExtractor")
        << "No entry found in PMT channel mapping database for board '"
        << readoutBoardConfig.boardName << "' (fragment ID: "
        << readoutBoardConfig.fragmentID << " => " << std::hex << fragmentID
        << ")\n";
      continue;
    }
  
    icarusDB::DigitizerChannelChannelIDPairVec const& digitizerChannelVec
      = fChannelMap->getChannelIDPairVec(fragmentID);
    
    // finds the channel ID matching the specified channel number of this board
    auto const toChannelID = [&channelIDs=digitizerChannelVec]
      (short unsigned int channelNo)
      {
        auto const it = std::find_if(channelIDs.begin(), channelIDs.end(),
          [channelNo](auto const& p){ return std::get<0U>(p) == channelNo; });
        return (it != channelIDs.end())
          ? std::get<1U>(*it)
          : sbn::V1730channelConfiguration::NoChannelID
          ;
      };
    
    for (auto& channelInfo: readoutBoardConfig.channels)
      channelInfo.channelID = toChannelID(channelInfo.channelNo);
    
  } // for boards
  
  return config;
} // icarus::PMTconfigurationExtractor::finalize()


// -----------------------------------------------------------------------------
std::optional<fhicl::ParameterSet>
icarus::PMTconfigurationExtractor::readBoardConfig
  (fhicl::ParameterSet const& pset, std::string const& key) const
{
  static std::string const ExpectedFragmentType = "CAENV1730";
  
  std::optional<fhicl::ParameterSet> config;
  
  do { // fake loop for fast exit
    
    // it must be a parameter set
    if (!pset.has_key(key) || !pset.is_key_to_table(key)) break;
    
    auto boardPSet = pset.get<fhicl::ParameterSet>(key);
    
    // its "fragment_type" must be the expected one
    std::string fragmentType;
    if (!boardPSet.get_if_present("daq.fragment_receiver.fragment_type", fragmentType))
      break;
    if (fragmentType != ExpectedFragmentType) break;
    
    config.emplace(std::move(boardPSet)); // success
  } while (false);
  return config;
} // icarus::PMTconfigurationExtractor::readBoardConfig()


// -----------------------------------------------------------------------------
// ---  free functions implementation
// -----------------------------------------------------------------------------
sbn::PMTconfiguration icarus::extractPMTreadoutConfiguration
  (std::string const& srcFileName, icarus::PMTconfigurationExtractor extractor)
{
  //
  // TFile::Open() call is needed to support non-local URL
  // (e.g. XRootD URL are not supported by TFile constructor).
  //
  return extractPMTreadoutConfiguration(*(
    std::unique_ptr<TFile>{ TFile::Open(srcFileName.c_str(), "READ") }
    ),
    std::move(extractor)
    );
} // icarus::extractPMTreadoutConfiguration(string)


// -----------------------------------------------------------------------------
sbn::PMTconfiguration icarus::extractPMTreadoutConfiguration
  (TFile& srcFile, icarus::PMTconfigurationExtractor extractor)
{
  
  return details::extractPMTreadoutConfigurationImpl
    (util::readConfigurationFromArtFile(srcFile), std::move(extractor));
  
} // icarus::extractPMTreadoutConfiguration(TFile)


// ---------------------------------------------------------------------------
unsigned int icarus::PMTconfigurationExtractor::readoutBoardDBfragmentID
  (sbn::V1730Configuration const& boardConfig)
{
  return boardConfig.fragmentID & 0xFF; // secret recipe
} // icarus::PMTconfigurationExtractor::readoutBoardDBfragmentID()


// -----------------------------------------------------------------------------

