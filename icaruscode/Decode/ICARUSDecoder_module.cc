///////////////////////////////////////////////////////////////////////
// Class:       ICARUSTPCDecoder
// Plugin Type: producer (art v2_09_06)
// File:        ICARUSTPCDecoder.cxx
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Principal/Run.h"
//#include "art/Framework/Principal/SubRun.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h" 

#include "artdaq-core/Data/Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"

#include "lardataobj/RawData/RawDigit.h"

#include "TTree.h"

//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

namespace daq
{

/*
  * The Decoder module takes as input "NevisTPCFragments" and
  * outputs raw::RawDigits. It also handles in and all issues
  * with the passed in header and fragments (or at least it will).
*/
class ICARUSDecoder : public art::EDProducer 
{
public:
    explicit ICARUSDecoder(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    // Plugins should not be copied or assigned.
    ICARUSDecoder(ICARUSDecoder const &) = delete;
    ICARUSDecoder(ICARUSDecoder &&) = delete;

    ICARUSDecoder & operator = (ICARUSDecoder const &) = delete;

    ICARUSDecoder & operator = (ICARUSDecoder &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // get checksum from a Icarus fragment
    //static uint32_t compute_checksum(sbnddaq::NevisTPCFragment &fragment);
    static uint32_t compute_checksum(icarus::PhysCrateFragment &fragment);
private:
    class Config {
    public:
        int      wait_sec;
        int      wait_usec;
        bool     produce_header;
        bool     produce_metadata;
        bool     baseline_calc;
        unsigned n_mode_skip;
        bool     subtract_pedestal;
        unsigned channel_per_slot;
        unsigned min_slot_no;
        // for converting frame time into timestamp
        unsigned timesize;
        double   frame_to_dt;

        Config(fhicl::ParameterSet const & p);
    };

    // process an individual fragment inside an art event
    void process_fragment(art::Event &event, const artdaq::Fragment &frag, std::unique_ptr<std::vector<raw::RawDigit>> &product_collectionf);
    // Gets the WIRE ID of the channel. This wire id can be then passed
    // to the Lariat geometry.
    //raw::ChannelID_t get_wire_id(const sbnddaq::NevisTPCHeader *header, uint16_t nevis_channel_id);
    // whether the given nevis readout channel is mapped to a wire
    //bool is_mapped_channel(const sbnddaq::NevisTPCHeader *header, uint16_t nevis_channel_id);


    // Tools for decoding fragments depending on type
    std::vector<std::unique_ptr<IDecoder>>  fDecoderToolVec;  ///< Decoder tools

    art::InputTag _tag;
    Config        _config;
    // keeping track of incrementing numbers
    uint32_t _last_event_number;
    uint32_t _last_trig_frame_number;

    uint32_t _fragment_id_offset;

    TTree* _header_ana_tree; 

    unsigned int _event_number;
    unsigned int _fragment_id;
    unsigned int _board_id;
    unsigned int _timestamp;

};

DEFINE_ART_MODULE(ICARUSDecoder)

ICARUSDecoder::ICARUSDecoder(fhicl::ParameterSet const & params): art::EDProducer{params},
    _tag(params.get<std::string>("raw_data_label", "daq"),params.get<std::string>("fragment_type_label", "PHYSCRATEDATA")),
    _config(params),
    _last_event_number(0),
    _last_trig_frame_number(0),
    _fragment_id_offset(params.get<uint32_t>("fragment_id_offset"))
{
  
    // produce stuff
//    produces<std::vector<raw::RawDigit>>();
    // if (_config.produce_header) {
    //   produces<std::vector<tpcAnalysis::HeaderData>>();
    // }

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& decoderTools = params.get<fhicl::ParameterSet>("DecoderToolVec");
    
    fDecoderToolVec.resize(decoderTools.get_pset_names().size());
    
    for(const std::string& decoderTool : decoderTools.get_pset_names())
    {
        const fhicl::ParameterSet& decoderToolParamSet = decoderTools.get<fhicl::ParameterSet>(decoderTool);
        
        // Get instance of tool
        fDecoderToolVec.push_back(art::make_tool<IDecoder>(decoderToolParamSet));

        // Announce any producers
        fDecoderToolVec.back()->produces(producesCollector());
    }

    art::ServiceHandle<art::TFileService> tfs;
    
    _header_ana_tree = tfs->make<TTree>("_header_ana_tree","Header Ana Tree");

    _header_ana_tree->Branch("fragment_id",&_fragment_id,"fragment_id/i");
    _header_ana_tree->Branch("board_id",&_board_id,"board_id/i");
    _header_ana_tree->Branch("event",&_event_number,"event/i");
    _header_ana_tree->Branch("timestamp",&_timestamp,"timestamp/i");

    return;
}

ICARUSDecoder::Config::Config(fhicl::ParameterSet const & param) 
{
    // amount of time to wait in between processing events
    // useful for debugging redis
    double wait_time = param.get<double>("wait_time", -1 /* units of seconds */);
    wait_sec = (int) wait_time;
    wait_usec = (int) (wait_time / 1000000.);

    // whether to calcualte the pedestal (and set it in SetPedestal())
    baseline_calc = param.get<bool>("baseline_calc", false);

    // whether to put headerinfo in the art root file
    produce_header = param.get<bool>("produce_header", false);

    // how many adc values to skip in mode/pedestal finding
    n_mode_skip = param.get<unsigned>("n_mode_skip", 1);

    // whether to subtract pedestal
    subtract_pedestal = param.get<bool>("subtract_pedestal", false);

    // icarus readout window length
    timesize = param.get<unsigned>("timesize", 1);

    // icarus tick length (for timestamp)
    // should be 1/(2.5MHz) = 0.4mus
    frame_to_dt = param.get<double>("frame_to_dt", 1);

    // number of channels in each slot
    channel_per_slot = param.get<unsigned>("channel_per_slot", 0);

    // index of 0th slot
    min_slot_no = param.get<unsigned>("min_slot_no", 0);

    // Initialize the tool.

    return;
}

void ICARUSDecoder::produce(art::Event & event)
{
    if (_config.wait_sec > 0)
      std::this_thread::sleep_for(std::chrono::seconds(_config.wait_sec) + std::chrono::microseconds(_config.wait_usec));

    auto const& daq_handle = event.getValidHandle<artdaq::Fragments>(_tag);
    
    // storage for waveform
//    std::unique_ptr<std::vector<raw::RawDigit>> product_collection(new std::vector<raw::RawDigit>);
    fDecoderToolVec.back()->initializeDataProducts();

    // storage for header info
   // std::unique_ptr<std::vector<tpcAnalysis::HeaderData>> header_collection(new std::vector<tpcAnalysis::HeaderData>);
    for (auto const &rawFrag: *daq_handle) 
    {
      //process_fragment(event, rawfrag, product_collection, header_collection);
//      process_fragment(event, rawfrag, product_collection);
        fDecoderToolVec.back()->process_fragment(rawFrag);
    }

//    event.put(std::move(product_collection));
    fDecoderToolVec.back()->outputDataProducts(event);

    return;
}


void ICARUSDecoder::process_fragment(art::Event                                     &event, 
                                        const artdaq::Fragment                      &frag, 
                                        std::unique_ptr<std::vector<raw::RawDigit>> &product_collection) 
{
    _fragment_id = frag.fragmentID() - _fragment_id_offset;
  
    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment fragment(frag);
    
    size_t nBoardsPerFragment = fragment.nBoards();
    size_t nChannelsPerBoard  = fragment.nChannelsPerBoard();

    //int channel_count=0;
    for(size_t i_b=0; i_b < fragment.nBoards(); i_b++)
    {
        _board_id     = i_b;
        _event_number = fragment.BoardEventNumber(i_b);
        _timestamp    = fragment.BoardTimeStamp(i_b);
  
        _header_ana_tree->Fill();

        std::cout << "-- Fragment ID: " << _fragment_id << ", board ID: " << _board_id << ", time: " << _timestamp << std::endl;

        size_t boardId = nChannelsPerBoard * (nBoardsPerFragment * _fragment_id + i_b);

        //A2795DataBlock const& block_data = *(crate_data.BoardDataBlock(i_b));
        for(size_t i_ch = 0; i_ch < fragment.nChannelsPerBoard(); ++i_ch)
        {
            //raw::ChannelID_t channel_num = (i_ch & 0xff ) + (i_b << 8);
            raw::ChannelID_t           channel_num = boardId + i_ch;
            raw::RawDigit::ADCvector_t wvfm(fragment.nSamplesPerChannel());

            for(size_t i_t=0; i_t < fragment.nSamplesPerChannel(); ++i_t) {
  	            wvfm[i_t] = fragment.adc_val(i_b,i_ch,i_t);
            }
        
            //   product_collection->emplace_back(channel_count++,fragment.nSamplesPerChannel(),wvfm);
            product_collection->emplace_back(channel_num,fragment.nSamplesPerChannel(),wvfm);
        }//loop over channels
    }//loop over boards

    return;
}

} // end namespace

