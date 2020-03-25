/**
 *  @file   TPCDecoder_tool.cc
 *
 *  @brief  This tool provides "standard" 3D hits built (by this tool) from 2D hits
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/TPC/Decode/DecoderTools/IDecoder.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  TPCDecoder class definiton
 */
class TPCDecoder : virtual public IDecoder
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit TPCDecoder(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TPCDecoder();

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

    using RawDigitCollectionPtr = std::unique_ptr<std::vector<raw::RawDigit>>;

    uint32_t                                    fFragment_id_offset;    //< Allow offset for id

    RawDigitCollectionPtr                       fRawDigitCollection;    //< The output data collection pointer

    const geo::Geometry*                        fGeometry;              //< pointer to the Geometry service
    const detinfo::DetectorProperties*          fDetector;              //< Pointer to the detector properties
};

TPCDecoder::TPCDecoder(fhicl::ParameterSet const &pset)
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TPCDecoder::~TPCDecoder()
{
}

void TPCDecoder::produces(art::ProducesCollector& collector)
{
    collector.produces< std::vector<raw::RawDigit>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------
void TPCDecoder::configure(fhicl::ParameterSet const &pset)
{
    fFragment_id_offset = pset.get<uint32_t>("fragment_id_offset");

    fGeometry = art::ServiceHandle<geo::Geometry const>{}.get();
    fDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();

    return;
}

void TPCDecoder::initializeDataProducts()
{
    fRawDigitCollection = RawDigitCollectionPtr(new std::vector<raw::RawDigit>);

    return;
}

void TPCDecoder::process_fragment(const artdaq::Fragment &fragment)
{
    size_t fragment_id = fragment.fragmentID() - fFragment_id_offset;
  
    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment physCrateFragment(fragment);
    
    size_t nBoardsPerFragment = physCrateFragment.nBoards();
    size_t nChannelsPerBoard  = physCrateFragment.nChannelsPerBoard();

    //int channel_count=0;
    for(size_t board = 0; board < physCrateFragment.nBoards(); board++)
    {
//        size_t event_number = physCrateFragment.BoardEventNumber(i_b);
//        size_t timestamp    = physCrateFragment.BoardTimeStamp(board);

        size_t boardId = nChannelsPerBoard * (nBoardsPerFragment * fragment_id + board);

        // Get the pointer to the start of this board's block of data
        const icarus::A2795DataBlock::data_t* dataBlock = physCrateFragment.BoardData(board);

        //A2795DataBlock const& block_data = *(crate_data.BoardDataBlock(i_b));
        for(size_t channel = 0; channel < physCrateFragment.nChannelsPerBoard(); channel++)
        {
            //raw::ChannelID_t channel_num = (i_ch & 0xff ) + (i_b << 8);
            raw::ChannelID_t           channel_num = boardId + channel;
            raw::RawDigit::ADCvector_t wvfm(physCrateFragment.nSamplesPerChannel());

            // It seems that the data is read from each channel for each tick so the 
            // loop indices below are chosen to pick out the "right" ticks for a given channel
            for(size_t tick = 0; tick < physCrateFragment.nSamplesPerChannel(); tick++)
                wvfm[tick] = dataBlock[channel + tick * physCrateFragment.nChannelsPerBoard()];
        
            fRawDigitCollection->emplace_back(channel_num,physCrateFragment.nSamplesPerChannel(),wvfm);
        }//loop over channels
    }//loop over boards

    return;
}

void TPCDecoder::outputDataProducts(art::Event& event)
{
    // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
    std::sort(fRawDigitCollection->begin(),fRawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});

    // Now transfer ownership to the event store
    event.put(std::move(fRawDigitCollection));

    return;
}

DEFINE_ART_CLASS_TOOL(TPCDecoder)
} // namespace lar_cluster3d
