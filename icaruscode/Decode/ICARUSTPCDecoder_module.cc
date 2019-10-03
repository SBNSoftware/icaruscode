///////////////////////////////////////////////////////////////////////
// Class:       ICARUSTPCDecoder
// Plugin Type: producer (art v2_09_06)
// File:        ICARUSTPCDecoder.cxx
//
////////////////////////////////////////////////////////////////////////
#include "ICARUSTPCDecoder.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <memory>
#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <thread>
#include "art/Framework/Core/ModuleMacros.h"
#include "artdaq-core/Data/Fragment.hh"
#include "lardataobj/RawData/RawDigit.h"
//#include "sbnddaq-datatypes/Overlays/NevisTPCFragment.hh"
//#include "sbnddaq-datatypes/NevisTPC/NevisTPCTypes.hh"
//#include "sbnddaq-datatypes/NevisTPC/NevisTPCUtilities.hh"
//#include "../HeaderData.hh"
//#include "../../Mode/Mode.hh"
DEFINE_ART_MODULE(daq::ICARUSTPCDecoder)
// constructs a header data object from a icarus header
// construct from a icarus header
/*tpcAnalysis::HeaderData daq::ICARUSTPCDecoder::Fragment2HeaderData(art::Event &event, const artdaq::Fragment &frag) {

      icarus::PhysCrateFragment fragment(frag);

  const sbnddaq::NevisTPCHeader *raw_header = fragment.header();
  tpcAnalysis::HeaderData ret;

  ret.crate = raw_header->getFEMID();
  ret.slot = raw_header->getSlot();
  ret.event_number = raw_header->getEventNum();
  // ret.frame_number = raw_header->getFrameNum();
  ret.checksum = raw_header->getChecksum();
  
  ret.adc_word_count = raw_header->getADCWordCount();
  // ret.trig_frame_number = raw_header->getTrigFrame();
  
  // formula for getting unix timestamp from nevis frame number:
  // timestamp = frame_number * (timesize + 1) + trigger_sample
  ret.timestamp = (raw_header->getFrameNum() * (_config.timesize + 1) + raw_header->get2mhzSample()) * _config.frame_to_dt;

  ret.index = raw_header->getSlot() - _config.min_slot_no;

  return ret;

}
*/
daq::ICARUSTPCDecoder::ICARUSTPCDecoder(fhicl::ParameterSet const & param): 
  art::EDProducer{param},
  _tag(param.get<std::string>("raw_data_label", "daq"),param.get<std::string>("fragment_type_label", "PHYSCRATEDATA")),
  _config(param),
  _last_event_number(0),
  _last_trig_frame_number(0)
 {
  
  // produce stuff
  produces<std::vector<raw::RawDigit>>();
 // if (_config.produce_header) {
 //   produces<std::vector<tpcAnalysis::HeaderData>>();
 // }
}
daq::ICARUSTPCDecoder::Config::Config(fhicl::ParameterSet const & param) {
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
void daq::ICARUSTPCDecoder::produce(art::Event & event)
{
  if (_config.wait_sec >= 0) {
    std::this_thread::sleep_for(std::chrono::seconds(_config.wait_sec) + std::chrono::microseconds(_config.wait_usec));
  }
  auto const& daq_handle = event.getValidHandle<artdaq::Fragments>(_tag);
  
  // storage for waveform
  std::unique_ptr<std::vector<raw::RawDigit>> product_collection(new std::vector<raw::RawDigit>);
  // storage for header info
 // std::unique_ptr<std::vector<tpcAnalysis::HeaderData>> header_collection(new std::vector<tpcAnalysis::HeaderData>);
  for (auto const &rawfrag: *daq_handle) {
    //process_fragment(event, rawfrag, product_collection, header_collection);
process_fragment(event, rawfrag, product_collection);
  }
  event.put(std::move(product_collection));
  //if (_config.produce_header) {
 //   event.put(std::move(header_collection));
 // }
}
/*
bool daq::ICARUSTPCDecoder::is_mapped_channel(const sbnddaq::NevisTPCHeader *header, uint16_t nevis_channel_id) {
  // TODO: make better
  return true;
}

raw::ChannelID_t daq::ICARUSTPCDecoder::get_wire_id(const sbnddaq::NevisTPCHeader *header, uint16_t nevis_channel_id) {
  // TODO: make better
  return (header->getSlot() - _config.min_slot_no) * _config.channel_per_slot + nevis_channel_id;
}
*/
void daq::ICARUSTPCDecoder::process_fragment(art::Event &event, const artdaq::Fragment &frag, 
  std::unique_ptr<std::vector<raw::RawDigit>> &product_collection) {
  // convert fragment to Nevis fragment
  icarus::PhysCrateFragment fragment(frag);
        std::cout << " n boards " << fragment.nBoards() << std::endl;
//int channel_count=0;
for(size_t i_b=0; i_b < fragment.nBoards(); i_b++){
        //A2795DataBlock const& block_data = *(crate_data.BoardDataBlock(i_b));
        for(size_t i_ch=0; i_ch < fragment.nChannelsPerBoard(); ++i_ch){
          //raw::ChannelID_t channel_num = (i_ch & 0xff ) + (i_b << 8);
           raw::ChannelID_t channel_num = i_ch+i_b*64;
          raw::RawDigit::ADCvector_t wvfm(fragment.nSamplesPerChannel());
          for(size_t i_t=0; i_t < fragment.nSamplesPerChannel(); ++i_t) {
            wvfm[i_t] = fragment.adc_val(i_b,i_ch,i_t);
           // if(channel_num==1855) std::cout << " sample " << i_t << " wave " << wvfm[i_t] << std::endl;
}
     //   product_collection->emplace_back(channel_count++,fragment.nSamplesPerChannel(),wvfm);
      product_collection->emplace_back(channel_num,fragment.nSamplesPerChannel(),wvfm);
 //std::cout << " channel " << channel_num << " waveform size " << fragment.nSamplesPerChannel() << std::endl;
        }//loop over channels
      }//loop over boards
      std::cout << "Total number of channels added is " << product_collection->size() << std::endl;
  /*std::unordered_map<uint16_t,sbnddaq::NevisTPC_Data_t> waveform_map;
  size_t n_waveforms = fragment.decode_data(waveform_map);
  (void)n_waveforms;

  if (_config.produce_header) {
    auto header_data = Fragment2HeaderData(event, frag);
    if (_config.produce_header || _config.produce_metadata) {
      // Construct HeaderData from the Nevis Header and throw it in the collection
      header_collection->push_back(header_data);
    }
  }

  for (auto waveform: waveform_map) {
    // ignore channels that aren't mapped to a wire
    if (!is_mapped_channel(fragment.header(), waveform.first)) continue;

    std::vector<int16_t> raw_digits_waveform;
    raw::ChannelID_t wire_id = get_wire_id(fragment.header(), waveform.first);
    // TODO: is this too expensive an operation?
    for (auto digit: waveform.second) {
      raw_digits_waveform.push_back( (int16_t) digit);
    }  
    // calculate the mode and set it as the pedestal
    if (_config.baseline_calc || _config.subtract_pedestal) {
      int16_t mode = Mode(raw_digits_waveform, _config.n_mode_skip);
      if (_config.subtract_pedestal) {
        for (unsigned i = 0; i < raw_digits_waveform.size(); i++) {
          raw_digits_waveform[i] -= mode;
        }
      }

      // construct the next RawDigit object
      product_collection->emplace_back(wire_id, raw_digits_waveform.size(), raw_digits_waveform);

      if (_config.baseline_calc) {
        (*product_collection)[product_collection->size() - 1].SetPedestal( mode);
      }
    }
    // just push back
    else {
      product_collection->emplace_back(wire_id, raw_digits_waveform.size(), raw_digits_waveform);

    }
  }
*/
}
// Computes the checksum, given a nevis tpc header
// Ideally this would be in sbnddaq-datatypes, but it's not and I can't
// make changes to it, so put it here for now
//
// Also note that this only works for uncompressed data
/*
uint32_t daq::ICARUSTPCDecoder::compute_checksum(icarus::PhysCrateFragment &fragment) {
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

