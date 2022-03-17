/**
 *  @file   icaruscode/Decode/DecoderTools/PMTLaserDecoder_tool.cc
 *
 *  @brief  Simple Decoder for data collected using the Laser system
 * 
 *  @author Andrea Scarpelli (ascarpell@bnl.gov)
 *
 */

// Framework Includes
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
// #include "cetlib/cpu_timer.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/intervals_fhicl.h" // for nanoseconds in FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/ExternalTrigger.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

#include "icaruscode/Decode/DecoderTools/details/PMTDecoderUtils.h"
#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSLaserMap.h"
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()
#include "icarusalg/Utilities/BinaryDumpUtils.h" // icarus::ns::util::bin()
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration

// ROOT libraries
#include "TTree.h"

// std includes
#include <ostream>
#include <algorithm> // std::lower_bound()
#include <string>
#include <vector>
#include <tuple>
#include <optional>
#include <memory>
#include <cmath> // std::floor()
#include <cassert>

using namespace util::quantities::time_literals;

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq { class PMTSpareChannelDecoder; }

// Doxygen bookeeping 

class daq::PMTSpareChannelDecoder: public IDecoder
{

    public:

        struct Config {
      
            using Name = fhicl::Name;
            using Comment = fhicl::Comment;

            fhicl::Atom<bool> GetTiming {
                Name("GetTiming"),
                Comment("Decode PMT channel 16 signal and extract timing information")
            };


        };

        using Parameters = art::ToolConfigTable<Config>;

        // --- END ---- FHiCL configuration ----------------------------------------
    
        explicit PMTSpareChannelDecoder(Parameters const& params);
    
        virtual void produces(art::ProducesCollector&) override;

        virtual void configure(const fhicl::ParameterSet&) override;
    
        virtual void initializeDataProducts() override;

        virtual void process_fragment(const artdaq::Fragment &fragment) override;

        virtual void outputDataProducts(art::Event& event) override;


    private: 

        bool fGetTiming;

        //nanoseconds const fOpticalTick;

        using OpDetWaveformCollection    = std::vector<raw::OpDetWaveform>;
        using OpDetWaveformCollectionPtr = std::unique_ptr<OpDetWaveformCollection>;

        OpDetWaveformCollectionPtr         fOpDetWaveformCollection;  ///< The output data collection pointer

        template <std::size_t NBits, typename T>
            constexpr std::pair<std::array<std::size_t, NBits>, std::size_t>
                setBitIndices(T value) noexcept;

};


daq::PMTSpareChannelDecoder::PMTSpareChannelDecoder(Parameters const& params)
    : fGetTiming{ params().GetTiming() }
    //, fChannelMap{ *(art::ServiceHandle<icarusDB::IICARUSLaserMap const>{}) }
    //, fOpticalTick{ fDetTimings.OpticalClockPeriod() }
{}


template <std::size_t NBits, typename T>
constexpr std::pair<std::array<std::size_t, NBits>, std::size_t>
daq::PMTSpareChannelDecoder::setBitIndices(T value) noexcept {
    
    std::pair<std::array<std::size_t, NBits>, std::size_t> res;
    auto& [ indices, nSetBits ] = res;
    for (std::size_t& index: indices) {
        index = (value & 1)? nSetBits++: NBits;
        value >>= 1;
    } // for
    return res;
    
} 


void daq::PMTSpareChannelDecoder::produces(art::ProducesCollector& collector){
    collector.produces<OpDetWaveformCollection>();
}


void daq::PMTSpareChannelDecoder::configure(fhicl::ParameterSet const&) {
  // Configuration all happens during construction
  throw cet::exception("PMTDecoder")
    << "This tool does not support reconfiguration.\n"; 
} 



void daq::PMTSpareChannelDecoder::initializeDataProducts() {
    fOpDetWaveformCollection = OpDetWaveformCollectionPtr(new OpDetWaveformCollection);
}


void daq::PMTSpareChannelDecoder::process_fragment( const artdaq::Fragment &artdaqFragment ) {

    size_t const fragment_id = artdaqFragment.fragmentID();
    size_t const eff_fragment_id = fragment_id & 0x0fff;

    sbndaq::CAENV1730Fragment fragment(artdaqFragment);
    sbndaq::CAENV1730FragmentMetadata metafrag = *fragment.Metadata();
    sbndaq::CAENV1730Event evt = *fragment.Event();
    sbndaq::CAENV1730EventHeader header = evt.Header;
    size_t nChannelsPerBoard = metafrag.nChannels;

  
    uint32_t ev_size_quad_bytes         = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(sbndaq::CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes     = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t nSamplesPerChannel         = data_size_double_bytes/nChannelsPerBoard;
    uint16_t const enabledChannels      = header.ChannelMask();

    //artdaq::Fragment::timestamp_t const fragmentTimestamp = artdaqFragment.timestamp(); 
    unsigned int const time_tag =  header.triggerTimeTag*8;

    auto const [ chDataMap, nEnabledChannels ] = setBitIndices<16U>(enabledChannels);
    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));

    std::vector<uint16_t> wvfm(nSamplesPerChannel);

    std::size_t const channelPosInData = chDataMap[15];

    if (channelPosInData >= nEnabledChannels) return; // not enabled
    
    std::size_t const ch_offset = channelPosInData * nSamplesPerChannel;

    std::copy_n(data_begin + ch_offset, nSamplesPerChannel, wvfm.begin());

    //uint32_t spareChWaveformStartTime = fragmentTimestamp;

    fOpDetWaveformCollection->emplace_back(time_tag, eff_fragment_id, wvfm);

}


void daq::PMTSpareChannelDecoder::outputDataProducts(art::Event& event) {
    
    // Now transfer ownership to the event store
    event.put(std::move(fOpDetWaveformCollection));

}

DEFINE_ART_CLASS_TOOL(daq::PMTSpareChannelDecoder)




