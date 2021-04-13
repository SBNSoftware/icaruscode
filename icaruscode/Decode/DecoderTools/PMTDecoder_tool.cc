/**
 *  @file   icaruscode/Decode/DecoderTools/PMTDecoder_tool.cc
 *
 *  @brief  This tool provides "standard" 3D hits built (by this tool) from 2D hits
 * 
 *  @author Andrea Scarpelli (ascarpell@bnl.gov)
 *
 */

// Framework Includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
// #include "cetlib/cpu_timer.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

// std includes
#include <string>
#include <memory>


//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 * @brief Turns PMT readout fragments from DAQ into LArSoft data products.
 * 
 * The tool can read fragments from CAEN V1730 readout boards delivered by
 * artDAQ.
 * 
 * This decoder must support both a off-line mode (for storage and downstream
 * processing) and a on-line mode (for monitoring).
 * In particular, the on-line workflow is such that it may not be possible to
 * access the FHiCL configuration of the job and therefore the PMT configuration
 * data (see `icarus::PMTconfigurationExtraction` module).
 * 
 * Configuration
 * --------------
 * 
 * The set of supported parameters can be seen on command line by running
 * `lar --print-description PMTDecoder`.
 * 
 * Description of the configuration parameters:
 * * `DiagnosticOutput` (flag, default: `false`): enables additional console
 *     output, including dumping of the fragments (that is huge output).
 * * `LogCategory` (string, default: `PMTDecoder`): name of the message facility
 *     category where the output is sent.
 * 
 */
class PMTDecoder: public IDecoder
{
public:
    
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<bool> DiagnosticOutput {
        Name("DiagnosticOutput"),
        Comment("enable additional console output"),
        false // default
        };
      
      fhicl::Atom<std::string> LogCategory {
        Name("LogCategory"),
        Comment("name of the category for message stream"),
        "PMTDecoder" // default
        };
      
    }; // Config
    
    using Parameters = art::ToolConfigTable<Config>;
  
    /**
     *  @brief  Constructor
     *
     *  @param  params configuration parameter set
     */
    explicit PMTDecoder(Parameters const& params);


    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::ProducesCollector&) override;

    /// Reconfiguration is not supported: all configuration at construction time.
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

    // --- BEGIN -- Configuration parameters -----------------------------------
    bool const        fDiagnosticOutput; ///< If true will spew endless messages to output.
    
    std::string const fLogCategory; ///< Message facility category.
    
    // --- END ---- Configuration parameters -----------------------------------

    
    // --- BEGIN -- Services ---------------------------------------------------
    
    geo::GeometryCore const&           fGeometry; ///< Geometry service provider.
    icarusDB::IICARUSChannelMap const& fChannelMap; ///< Fragment/channel mapping database.

    // --- END ---- Services ---------------------------------------------------

    using OpDetWaveformCollection    = std::vector<raw::OpDetWaveform>;
    using OpDetWaveformCollectionPtr = std::unique_ptr<OpDetWaveformCollection>;

    OpDetWaveformCollectionPtr         fOpDetWaveformCollection;  ///< The output data collection pointer
    
}; // class PMTDecoder


//------------------------------------------------------------------------------------------------------------------------------------------
PMTDecoder::PMTDecoder(Parameters const& params)
  : fDiagnosticOutput{ params().DiagnosticOutput() }
  , fLogCategory{ params().LogCategory() }
  , fGeometry{ *(lar::providerFrom<geo::Geometry const>()) }
  , fChannelMap{ *(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}) }
{
}

//------------------------------------------------------------------------------------------------------------------------------------------


void PMTDecoder::produces(art::ProducesCollector& collector)
{
    collector.produces<OpDetWaveformCollection>();
}

//------------------------------------------------------------------------------------------------------------------------------------------
void PMTDecoder::configure(fhicl::ParameterSet const&) {
  // Configuration all happens during construction
  throw cet::exception("PMTDecoder")
    << "This tool does not support reconfiguration.\n"; 
} // PMTDecoder::configure()


//------------------------------------------------------------------------------------------------------------------------------------------
void PMTDecoder::initializeDataProducts()
{
    fOpDetWaveformCollection = OpDetWaveformCollectionPtr(new OpDetWaveformCollection);
}

void PMTDecoder::process_fragment(const artdaq::Fragment &artdaqFragment)
{
    size_t fragment_id = artdaqFragment.fragmentID() & 0x0fff;

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
        mf::LogVerbatim(fLogCategory)
          << "----> PMT Fragment ID: " << fragment_id << ", boardID: " << boardId
            << ", nChannelsPerBoard: " << nChannelsPerBoard
            << ", nSamplesPerChannel: " << nSamplesPerChannel
          << "\n      size: " << ev_size_quad_bytes
            << ", data size: " << data_size_double_bytes
            << ", samples/channel: " << nSamplesPerChannel
            << ", time: " << time_tag
          ;
    }

    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));
    const uint16_t* value_ptr  =  data_begin;
    uint16_t        value      = 0;
    size_t          ch_offset  = 0;

    // Temporary? 
    time_tag = 0;

    // Recover the information for this fragment
    if (fChannelMap.hasPMTDigitizerID(fragment_id))
    {
        const icarusDB::DigitizerChannelChannelIDPairVec& digitizerChannelVec
          = fChannelMap.getChannelIDPairVec(fragment_id);

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
    else {
      mf::LogError(fLogCategory)
        << "*** PMT could not find channel information for fragment: "
          << fragment_id;
    }

    if (fDiagnosticOutput) {
      mf::LogVerbatim(fLogCategory)
        << "      - size of output collection: " << fOpDetWaveformCollection->size();
    }
    
} // PMTDecoder::process_fragment()

void PMTDecoder::outputDataProducts(art::Event& event)
{
    // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
    std::sort(fOpDetWaveformCollection->begin(),fOpDetWaveformCollection->end(),[](const auto& left,const auto&right){return left.ChannelNumber() < right.ChannelNumber();});

    // Now transfer ownership to the event store
    event.put(std::move(fOpDetWaveformCollection));

}

DEFINE_ART_CLASS_TOOL(PMTDecoder)
} // namespace daq
