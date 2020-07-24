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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "icaruscode/Decode/TPCChannelmapping.h"

// std includes
#include <string>
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

    bool                                     fDiagnosticOutput;       //< If true will spew endless messages to output

    using OpDetWaveformCollection    = std::vector<raw::OpDetWaveform>;
    using OpDetWaveformCollectionPtr = std::unique_ptr<OpDetWaveformCollection>;

    OpDetWaveformCollectionPtr              fOpDetWaveformCollection;  //< The output data collection pointer
    database::FragmentToDigitizerChannelMap fFragmentToDigitizerMap;   //< Channel mapping from daq to LArSoft
    const geo::Geometry*                    fGeometry;                 //< pointer to the Geometry service
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

    fGeometry = art::ServiceHandle<geo::Geometry const>{}.get();

    // Do the channel mapping initialization
    if (database::BuildFragmentToDigitizerChannelMap(fFragmentToDigitizerMap))
    {
        throw cet::exception("PMTDecoder") << "Cannot recover the Fragment ID channel map from the database \n";
    }
    else if (fDiagnosticOutput)
    {
        std::cout << "FragmentID to Readout ID map has " << fFragmentToDigitizerMap.size() << " Fragment IDs";
        for(const auto& pair : fFragmentToDigitizerMap) std::cout << "   Frag: " << std::hex << pair.first << ", # pairs: " << std::dec << pair.second.size() << std::endl;
    }

    return;
}

void PMTDecoder::initializeDataProducts()
{
    fOpDetWaveformCollection = OpDetWaveformCollectionPtr(new OpDetWaveformCollection);

    return;
}

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
        std::cout << "----> Fragment ID: " << fragment_id << ", boardID: " << boardId << ", nChannelsPerBoard: " << nChannelsPerBoard << ", nSamplesPerChannel: " << nSamplesPerChannel << std::endl;
        std::cout << "      size: " << ev_size_quad_bytes << ", data size: " << data_size_double_bytes << ", samples/channel: " << nSamplesPerChannel << ", time: " << time_tag << std::endl;
    }

    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
    const uint16_t* value_ptr  =  data_begin;
    uint16_t        value      = 0;
    size_t          ch_offset  = 0;

    // Temporary? 
    time_tag = 0;

    // Recover the information for this fragment
    if (fFragmentToDigitizerMap.find(fragment_id) != fFragmentToDigitizerMap.end())
    {
        const database::DigitizerChannelChannelIDPairVec& digitizerChannelVec = fFragmentToDigitizerMap[fragment_id];

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
    else std::cout << "*** could not find channel information for fragment: " << fragment_id << std::endl;

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
