/**
 *  @file   PMTDecoder_tool.cc
 *
 *  @brief  This tool provides "standard" 3D hits built (by this tool) from 2D hits
 * 
 *  This code provided by Andrea Scarpelli
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

#include "sbnobj/ICARUS/PMT/Data/PMTconfiguration.h"
#include "sbnobj/ICARUS/PMT/Data/V1730Configuration.h"
#include "sbnobj/ICARUS/PMT/Data/V1730channelConfiguration.h"
#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/DecoderTools/PMTconfigurationExtractor.h"
#include "icaruscode/Utilities/ReadArtConfiguration.h" // util::readConfigurationFromArtFile()

// ROOT
#include "TFile.h"

// std includes
#include <regex>
#include <initializer_list>
#include <string>
#include <optional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  PMTDecoder class definiton
 */
class PMTDecoder : virtual public IDecoder
{
public:
    
    struct Config {
      
      fhicl::Atom<bool> DiagnosticOutput {
        fhicl::Name("DiagnosticOutput"),
        fhicl::Comment("enable additional console output"),
        false // default
        };
      
      fhicl::Atom<bool> ExtractConfiguration {
        fhicl::Name("ExtractConfiguration"),
        fhicl::Comment("creates a data product with the PMT configuration"),
        false // default
        };
      
    }; // Config

  
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit PMTDecoder(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~PMTDecoder();

    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::ProducesCollector&) override;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) override;

    /**
     *  @brief Initialize any data products the tool will output
     *
     */
    virtual void initializeDataProducts() override;

    
    /// Extracts the PMT readout board configuration from the specified `file`
    /// and keeps it in the object.
    virtual void decodeConfigurationFromFile(TFile& file) override;
    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param fragment            The artdaq fragment to process
     *  @param rawDigitColllection The output RawDigits
     */
    virtual void process_fragment(const artdaq::Fragment &fragment) override;

    /**
     *  @brief Output the data products to the event store
     * 
     *  @param event The event store objects
     */
    virtual void outputDataProducts(art::Event& event) override;

private:

    bool                               fDiagnosticOutput;       ///< If true will spew endless messages to output
    bool                               fExtractConfig; ///< Parse PMT readout configuration.

    using OpDetWaveformCollection    = std::vector<raw::OpDetWaveform>;
    using OpDetWaveformCollectionPtr = std::unique_ptr<OpDetWaveformCollection>;

    OpDetWaveformCollectionPtr         fOpDetWaveformCollection;  ///< The output data collection pointer
    const geo::Geometry*               fGeometry = nullptr;       ///< pointer to the Geometry service
    const icarusDB::IICARUSChannelMap* fChannelMap = nullptr;
    
    
    /**
     * @brief Returns a parameter set with the content of
     *        `configuration_documents` key from `container`.
     * @param container parameter set including a `configuration_documents`
     * @param components keys to be converted (as regular expressions)
     * @return a parameter set
     * 
     * The `configuration_documents` element of `container` is processed: for
     * each of its keys which match at least one of the `components` regular
     * expression patterns (`std::regex_match()`), the associated string value
     * is parsed with FHiCL parser, and the result is set as a FHiCL table in
     * the output parameter set.
     * For example, if the `components` are
     * `{ std::regex{".*pmt.*"}, std::regex{".*trigger.*"} }`, the returned
     * value is a parameter set that may have keys like `icaruspmtee01`,
     * `icaruspmtew02`, `icarustrigger` etc., each one with a FHiCL table as
     * `value.
     */
    static fhicl::ParameterSet convertConfigurationDocuments(
      fhicl::ParameterSet const& container,
      std::initializer_list<std::regex const> components
      );
    
};

PMTDecoder::PMTDecoder(fhicl::ParameterSet const &pset)
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PMTDecoder::~PMTDecoder()
{
}

void PMTDecoder::produces(art::ProducesCollector& collector)
{
    collector.produces<OpDetWaveformCollection>();
}

//------------------------------------------------------------------------------------------------------------------------------------------
void PMTDecoder::configure(fhicl::ParameterSet const &pset)
{
    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);
    fExtractConfig = pset.get<bool>("ExtractReadoutConfig", false);

    fGeometry   = art::ServiceHandle<geo::Geometry const>{}.get();
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    return;
}

void PMTDecoder::initializeDataProducts()
{
    fOpDetWaveformCollection = OpDetWaveformCollectionPtr(new OpDetWaveformCollection);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------
fhicl::ParameterSet PMTDecoder::convertConfigurationDocuments(
  fhicl::ParameterSet const& container,
  std::initializer_list<std::regex const> components
) {
  static std::string const ConfigListKey { "configuration_documents" };
  
  fhicl::ParameterSet const sourceConfig
    = container.get<fhicl::ParameterSet>(ConfigListKey);
  
  fhicl::ParameterSet configDocs;
  for (auto const& key: sourceConfig.get_names()) {
    if (!sourceConfig.is_key_to_atom(key)) continue;
    
    // filter by key
    if (components.size() > 0U) {
      bool keep = false;
      for (auto const& pattern: components) {
        if (!std::regex_match(key, pattern)) continue;
        keep = true;
        break;
      }
      if (!keep) continue;
    } // if filtering
    
    std::string const psetStr = sourceConfig.get<std::string>(key);
    
    fhicl::ParameterSet pset;
    try {
      fhicl::make_ParameterSet(psetStr, pset);
    }
    catch (cet::exception& e) {
      throw cet::exception{ "PMTDecoder", "", e }
        << "Error parsing the content of key '" << ConfigListKey << "." << key
        << "'; content was:\n" << psetStr << "\n";
    }
    
    configDocs.put(key, pset);
    
  } // for all main keys
  
  return configDocs;
} // PMTDecoder::convertConfigurationDocuments()


void PMTDecoder::decodeConfigurationFromFile(TFile& file) {
  
  if (!fExtractConfig) return;
  
  /*
   * The plan is to look in all the FHiCL configuration fragments we can find
   * in the input file, and find all the useful configuration therein.
   * Given that there may be multiple input files, there may also be multiple
   * configurations for the same detector components.
   * In that case, we will extract parameters from each and every one of the
   * configurations, and throw an exception if they are not all consistent.
   * 
   * Consistency is tested only for the extracted parameters, not for the whole
   * FHiCL configuration fragment.
   */
  icarus::PMTconfigurationExtractor extractor;
  
  auto const& globalConfigColl = util::readConfigurationFromArtFile(file);
  
  std::optional<icarus::PMTconfiguration> config;
  
  // look in the global configuration for all parameter sets which contain
  // `configuration_documents` as a (direct) name;
  for (auto const& [ id, pset ]: globalConfigColl) {
    if (!pset.has_key("configuration_documents")) continue;
    
    fhicl::ParameterSet const configDocs
      = convertConfigurationDocuments(pset, { std::regex{ "icaruspmt.*" } });
    
    icarus::PMTconfiguration candidateConfig = extractor.extract(configDocs);
    if (config) {
      if (config.value() == candidateConfig) continue;
      mf::LogError log("PMTDecoder");
      log << "Found two candidate configurations differring:"
        "\nFirst:\n" << config.value()
        << "\nSecond:\n" << candidateConfig
        ;
      throw cet::exception("PMTDecoder")
        << "PMTDecoder::decodeConfigurationFromFile() found inconsistent configurations.\n";
    } // if incompatible configurations
    
    config.emplace(std::move(candidateConfig));
  } // for all configuration documents
  
  if (!config) {
    throw cet::exception("PMTDecoder")
      << "PMTDecoder::decodeConfigurationFromFile() could not find a suitable configuration.\n";
  }
  
  for (icarus::V1730Configuration& readoutBoardConfig: config->boards) {
    if (!fChannelMap->hasPMTDigitizerID(readoutBoardConfig.fragmentID))
      continue;
    icarusDB::DigitizerChannelChannelIDPairVec const& digitizerChannelVec
      = fChannelMap->getChannelIDPairVec(readoutBoardConfig.fragmentID);
    
    // finds the channel ID matching the specified channel number of this board
    auto const toChannelID = [&channelIDs=digitizerChannelVec]
      (short unsigned int channelNo)
      {
        auto const it = std::find_if(channelIDs.begin(), channelIDs.end(),
          [channelNo](auto const& p){ return p.first == channelNo; });
        return (it != channelIDs.end())
          ? it->second
          : icarus::V1730channelConfiguration::NoChannelID
          ;
      };
    
    for (auto& channelInfo: readoutBoardConfig.channels)
      channelInfo.channelID = toChannelID(channelInfo.channelNo);
    
  } // for boards
  
  mf::LogInfo("PMTDecoder")
    << "PMT readout: " << config.value();
  
} // PMTDecoder::decodeConfigurationFromFile()

//------------------------------------------------------------------------------------------------------------------------------------------
void PMTDecoder::process_fragment(const artdaq::Fragment &artdaqFragment)
{
    size_t fragment_id = artdaqFragment.fragmentID();

    // convert fragment to Nevis fragment
    sbndaq::CAENV1730Fragment         fragment(artdaqFragment);
    sbndaq::CAENV1730FragmentMetadata metafrag = *fragment.Metadata();
    sbndaq::CAENV1730Event            evt      = *fragment.Event();
    sbndaq::CAENV1730EventHeader      header   = evt.Header;

    size_t nChannelsPerBoard  = metafrag.nChannels; //fragment.nChannelsPerBoard();

    //std::cout << "Decoder:channels_per_board: " << nChannelsPerBoard << std::endl;

    uint32_t ev_size_quad_bytes         = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(sbndaq::CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes     = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t nSamplesPerChannel         = data_size_double_bytes/nChannelsPerBoard;

    raw::TimeStamp_t time_tag =  header.triggerTimeTag;

    size_t boardId = nChannelsPerBoard * fragment_id;

    if (fDiagnosticOutput)
    {
        std::cout << "----> PMT Fragment ID: " << fragment_id << ", boardID: " << boardId << ", nChannelsPerBoard: " << nChannelsPerBoard << ", nSamplesPerChannel: " << nSamplesPerChannel << std::endl;
        std::cout << "      size: " << ev_size_quad_bytes << ", data size: " << data_size_double_bytes << ", samples/channel: " << nSamplesPerChannel << ", time: " << time_tag << std::endl;
    }

    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
    const uint16_t* value_ptr  =  data_begin;
    uint16_t        value      = 0;
    size_t          ch_offset  = 0;

    // Temporary? 
    time_tag = 0;

    // Recover the information for this fragment
    if (fChannelMap->hasPMTDigitizerID(fragment_id))
    {
        const icarusDB::DigitizerChannelChannelIDPairVec& digitizerChannelVec = fChannelMap->getChannelIDPairVec(fragment_id);

        // Allocate the vector outside the loop just since we'll resuse it over and over... 
        std::vector<uint16_t> wvfm(nSamplesPerChannel);

        for(const auto& digitizerChannelPair : digitizerChannelVec)
        {
            size_t         digitizerChannel = digitizerChannelPair.first;
            raw::Channel_t channelID        = digitizerChannelPair.second;

            ch_offset = digitizerChannel * nSamplesPerChannel;

            for(size_t i_t=0; i_t < nSamplesPerChannel; ++i_t)
            {
                value_ptr = data_begin + ch_offset + i_t; /*pointer arithmetic*/
                value     = *(value_ptr);
  	            wvfm[i_t] = value;
            }

            fOpDetWaveformCollection->emplace_back(time_tag, channelID, wvfm);
        } 
    }
    else std::cout << "*** PMT could not find channel information for fragment: " << fragment_id << std::endl;

    if (fDiagnosticOutput) std::cout << "      - size of output collection: " << fOpDetWaveformCollection->size() << std::endl;

    return;
}

void PMTDecoder::outputDataProducts(art::Event& event)
{
    // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
    std::sort(fOpDetWaveformCollection->begin(),fOpDetWaveformCollection->end(),[](const auto& left,const auto&right){return left.ChannelNumber() < right.ChannelNumber();});

    // Now transfer ownership to the event store
    event.put(std::move(fOpDetWaveformCollection));

    return;
}

DEFINE_ART_CLASS_TOOL(PMTDecoder)
} // namespace lar_cluster3d
