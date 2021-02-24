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
    
    /// Configuration of PMT readout as extracted from the FHiCL configuration.
    std::optional<icarus::PMTconfiguration> fPMTconfig;
    
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

void PMTDecoder::decodeConfigurationFromFile(TFile& file) {
  
  // there should be a way to read this information from a data product instead;
  // the module for the data product is `PMTconfigurationExtraction`.
  
  if (!fExtractConfig) return;
  
  icarus::PMTconfigurationExtractor extractor { *fChannelMap };
  
  fPMTconfig.emplace(extractPMTreadoutConfiguration(file, extractor));
  
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
