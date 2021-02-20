/**
 * @file   icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.cxx
 * @brief  Utility to extract PMT readout configuration from data.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 18, 2021
 * @see    icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.h
 */


// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.h"


// -----------------------------------------------------------------------------
icarus::PMTconfiguration icarus::PMTconfigurationExtractor::extract
  (fhicl::ParameterSet const& pset) const
{
  
  icarus::PMTconfiguration config;
  
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
  -> icarus::V1730Configuration
{
  
  auto const& boardParams
    = pset.get<fhicl::ParameterSet>("daq.fragment_receiver");
  
  icarus::V1730Configuration rc; // readout config, for friends
  rc.boardName = boardName;
  
  rc.boardID = boardParams.get<unsigned int>("board_id");
  rc.fragmentID = boardParams.get<unsigned int>("fragment_id");
  
  rc.bufferLength = boardParams.get<int>("recordLength");
  rc.postTriggerFrac = boardParams.get<float>("postPercent") / 100.0f;
  
  rc.nChannels = boardParams.get<unsigned int>("nChannels");
  
  for (unsigned short int channelNo = 0; channelNo < 16; ++channelNo)
    rc.channels.push_back(extractChannelConfiguration(boardParams, channelNo));
  
  return rc;
} // icarus::PMTconfigurationExtractor::extractV1730configuration()


// -----------------------------------------------------------------------------
icarus::V1730channelConfiguration
icarus::PMTconfigurationExtractor::extractChannelConfiguration
  (fhicl::ParameterSet const& boardPSet, unsigned short int channelNo) const
{
  std::string const ChannelStr = std::to_string(channelNo);
  
  icarus::V1730channelConfiguration channel;
  
  channel.channelNo = channelNo;
  
  channel.baseline = boardPSet.get<short signed int>
    ("BaselineCh" + std::to_string(channelNo + 1));
  channel.threshold = boardPSet.get<short signed int>
    ("triggerThreshold" + ChannelStr);
  channel.enabled = boardPSet.get<bool>("channelEnable" + ChannelStr);
  
  return channel;
} // icarus::PMTconfigurationExtractor::extractChannelConfiguration()


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

