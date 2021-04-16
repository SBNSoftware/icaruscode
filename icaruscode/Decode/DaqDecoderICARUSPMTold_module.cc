////////////////////////////////////////////////////////////////////////
// Class:       DaqDecoderIcarus
// Plugin Type: producer (art v2_09_06)
// File:        DaqDecoderIcarus.cxx
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "artdaq-core/Data/Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

//#include "sbnddaq-datatypes/Overlays/NevisTPCFragment.hh"

//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

namespace daq 
{

class DaqDecoderIcarusPMTold : public art::EDProducer 
{
public:
    explicit DaqDecoderIcarusPMTold(fhicl::ParameterSet const & p);

    // Required functions.
    void produce(art::Event & e) override;

    // get checksum from a Icarus fragment
    //static uint32_t compute_checksum(sbnddaq::NevisTPCFragment &fragment);
//    static uint32_t compute_checksum(icarus::PhysCrateFragment &fragment);

private:
    class Config 
    {
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

    art::InputTag _tag;
    Config        _config;
    // keeping track of incrementing numbers
    // commented out while unused
    //uint32_t      _last_event_number;
    //uint32_t      _last_trig_frame_number;
};


DEFINE_ART_MODULE(daq::DaqDecoderIcarusPMTold)

DaqDecoderIcarusPMTold::DaqDecoderIcarusPMTold(fhicl::ParameterSet const & param): 
    art::EDProducer{param},
    _tag(param.get<art::InputTag>("FragmentsLabel", "daq:CAENV1730")),
    _config(param)//,
    //_last_event_number(0),
    //_last_trig_frame_number(0)
{
  
    // produce stuff
    produces<std::vector<raw::OpDetWaveform>>();
    // if (_config.produce_header) {
    //   produces<std::vector<tpcAnalysis::HeaderData>>();
    // }
}

DaqDecoderIcarusPMTold::Config::Config(fhicl::ParameterSet const & param) 
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
}

void DaqDecoderIcarusPMTold::produce(art::Event & event)
{
    if (_config.wait_sec >= 0) 
    {
        std::this_thread::sleep_for(std::chrono::seconds(_config.wait_sec) + std::chrono::microseconds(_config.wait_usec));
    }
    //auto const& daq_handle = event.getValidHandle<artdaq::Fragments>(_tag);
  
    // storage for waveform
    std::unique_ptr<std::vector<raw::OpDetWaveform>> product_collection(new std::vector<raw::OpDetWaveform>());
    // storage for header info
    // std::unique_ptr<std::vector<tpcAnalysis::HeaderData>> header_collection(new std::vector<tpcAnalysis::HeaderData>);

    /************************************************************************************************/
    art::Handle< std::vector<artdaq::Fragment> > rawFragHandle;
    event.getByLabel(_tag, rawFragHandle);
  
    if (rawFragHandle.isValid()) {

    std::cout << "######################################################################\n";
    std::cout << "Run " << event.run() << ", subrun " << event.subRun() << std::endl;

    for (size_t idx = 0; idx < rawFragHandle->size(); ++idx) /*loop over the fragments*/
    {
        //--use this fragment as a reference to the same data
        const auto& frag((*rawFragHandle)[idx]); 
        sbndaq::CAENV1730Fragment bb(frag);
        auto const* md = bb.Metadata();
        sbndaq::CAENV1730Event const* event_ptr = bb.Event();
        sbndaq::CAENV1730EventHeader header = event_ptr->Header;

        std::cout << "\tFragment ID: " << frag.fragmentID() << ", type: " << frag.typeString() << ", type: " << frag.type() << ", boardID: " << header.boardID << ", msk lo: " << header.channelMask_lo << ", hi: " << header.channelMask_hi << std::endl;
        std::cout << "\tFrom header, event counter is "  << header.eventCounter   << "\n";
        std::cout << "\tFrom header, triggerTimeTag is " << header.triggerTimeTag << "\n";
        std::vector< std::vector<uint16_t> >  fWvfmsVec;
        size_t nChannels = md->nChannels;
        std::cout <<"\tFrom header , no of channel is " << nChannels << "\n";
      
        uint32_t ev_size_quad_bytes = header.eventSize;
        //   std::cout << "Event size in quad bytes is: " << ev_size_quad_bytes << "\n";
        uint32_t evt_header_size_quad_bytes = sizeof(sbndaq::CAENV1730EventHeader)/sizeof(uint32_t);
        uint32_t data_size_double_bytes = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
        uint32_t wfm_length = data_size_double_bytes/nChannels;

        //--note, needs to take into account channel mask
        //   std::cout << "Channel waveform length = " << wfm_length << "\n";

        const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(frag.dataBeginBytes() 
                                   + sizeof(sbndaq::CAENV1730EventHeader));

        const uint16_t* value_ptr =  data_begin;
        size_t ch_offset = 0;
        uint16_t value;

        // loop over channels
        for (size_t i_ch=0; i_ch<nChannels; ++i_ch)
        {
            //fWvfmsVec[i_ch].resize(wfm_length);
            ch_offset = i_ch * wfm_length;
      
            raw::OpDetWaveform my_wf(0.00, i_ch, wfm_length);
            my_wf.resize(wfm_length);
            // Loop over waveform
      
            for(size_t i_t=0; i_t<wfm_length; ++i_t)
            { 
                //fTicksVec[i_t] = t0*Ttt_DownSamp + i_t;   /*timestamps, event level*/
                value_ptr = data_begin + ch_offset + i_t; /*pointer arithmetic*/
	              //value_ptr = (i_t%2 == 0) ? (index+1) : (index-1); 
                value = *(value_ptr);
                my_wf[i_t] = value;
	              //std::cout<<"Value is" << value << "and" << my_wf[i_t] << "\n";
	              //if (i_ch == 0 && fEvent == 0) h_wvfm_ev0_ch0->SetBinContent(i_t,value);
                //fWvfmsVec[i_ch][i_t] = value;
            }  //--end loop samples
 
            product_collection->push_back(my_wf);
          }// end loop over channels
        // std::cout<<"product collection"<<product_collection->back().size()<<"\n";
      } // end loop over fragments
    }
    event.put(std::move(product_collection));
}


/*
bool daq::DaqDecoderIcarusPMTold::is_mapped_channel(const sbnddaq::NevisTPCHeader *header, uint16_t nevis_channel_id) {
  // TODO: make better
  return true;
}

raw::ChannelID_t daq::DaqDecoderIcarusPMTold::get_wire_id(const sbnddaq::NevisTPCHeader *header, uint16_t nevis_channel_id) {
  // TODO: make better
  return (header->getSlot() - _config.min_slot_no) * _config.channel_per_slot + nevis_channel_id;
}
*/
void daq::DaqDecoderIcarusPMTold::process_fragment(art::Event &event, const artdaq::Fragment &frag, std::unique_ptr<std::vector<raw::RawDigit>> &product_collection) 
{
    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment fragment(frag);
	  std::cout << " n boards " << fragment.nBoards() << std::endl;
    //int channel_count=0;

    for(size_t i_b=0; i_b < fragment.nBoards(); i_b++)
    {
	      //A2795DataBlock const& block_data = *(crate_data.BoardDataBlock(i_b));

	      for(size_t i_ch=0; i_ch < fragment.nChannelsPerBoard(); ++i_ch)
        {
	          //raw::ChannelID_t channel_num = (i_ch & 0xff ) + (i_b << 8);
            raw::ChannelID_t channel_num = i_ch+i_b*64;
	          raw::RawDigit::ADCvector_t wvfm(fragment.nSamplesPerChannel());
	  
            for(size_t i_t=0; i_t < fragment.nSamplesPerChannel(); ++i_t) 
            {
	              wvfm[i_t] = fragment.adc_val(i_b,i_ch,i_t);
                // if(channel_num==1855) std::cout << " sample " << i_t << " wave " << wvfm[i_t] << std::endl;
            }
            //   product_collection->emplace_back(channel_count++,fragment.nSamplesPerChannel(),wvfm);
            product_collection->emplace_back(channel_num,fragment.nSamplesPerChannel(),wvfm);
            //std::cout << " channel " << channel_num << " waveform size " << fragment.nSamplesPerChannel() << std::endl;

	      }//loop over channels

    }//loop over boards
    std::cout << "Total number of channels added is " << product_collection->size() << std::endl;
}


// Computes the checksum, given a nevis tpc header
// Ideally this would be in sbnddaq-datatypes, but it's not and I can't
// make changes to it, so put it here for now
//
// Also note that this only works for uncompressed data
/*
uint32_t daq::DaqDecoderIcarusPMTold::compute_checksum(icarus::PhysCrateFragment &fragment) {
  uint32_t checksum = 0;

  const sbnddaq::NevisTPC_ADC_t* data_ptr = fragment.data();
  // RETURN VALUE OF getADCWordCount IS OFF BY 1
  size_t n_words = fragment.header()->getADCWordCount() + 1;

  for (size_t word_ind = 0; word_ind < n_words; word_ind++) {
    const sbnddaq::NevisTPC_ADC_t* word_ptr = data_ptr + word_ind;
    checksum += *word_ptr;
  }
  // only first 6 bytes of checksum are used
  return checksum & 0xFFFFFF;

}
*/

} // end of namespace


 

