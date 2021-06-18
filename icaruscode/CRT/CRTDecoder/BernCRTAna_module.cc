#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"

//#include "art/Framework/Services/Optional/TFileService.h"
#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TNtuple.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>

namespace sbndaq {
  class BernCRTAna;
}

class sbndaq::BernCRTAna : public art::EDAnalyzer {

public:
  explicit BernCRTAna(fhicl::ParameterSet const & pset);
  virtual ~BernCRTAna();
  
  virtual void analyze(art::Event const & evt);
  

private:
  void analyze_fragment(artdaq::Fragment & frag);
  
  TTree * events;

//data payload
  uint8_t mac5; //last 8 bits of FEB mac5 address
  uint16_t flags;
  uint16_t lostcpu;
  uint16_t lostfpga;
  uint32_t ts0;
  uint32_t ts1;
  uint16_t adc[32];
  uint32_t coinc;

//metadata
  uint64_t  run_start_time;
  uint64_t  this_poll_start;
  uint64_t  this_poll_end;
  uint64_t  last_poll_start;
  uint64_t  last_poll_end;
  int32_t   system_clock_deviation;
  uint32_t  feb_events_per_poll;
  uint32_t  feb_event_number;

//information from fragment header
  uint32_t  sequence_id;
  uint64_t  fragment_timestamp;
};

//Define the constructor
sbndaq::BernCRTAna::BernCRTAna(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{

  art::ServiceHandle<art::TFileService> tfs; //pointer to a file named tfs

  events = tfs->make<TTree>("events", "FEB events");

//event data
  events->Branch("mac5",        &mac5,          "mac5/b");
  events->Branch("flags",       &flags,         "flags/s");
  events->Branch("lostcpu",     &lostcpu,       "lostcpu/s");
  events->Branch("lostfpga",    &lostfpga,      "lostfpga/s");
  events->Branch("ts0",         &ts0,           "ts0/i");
  events->Branch("ts1",         &ts1,           "ts1/i");
  events->Branch("adc",         &adc,           "adc[32]/s");
  events->Branch("coinc",       &coinc,         "coinc/i");

//metadata
  events->Branch("run_start_time",            &run_start_time,              "run_start_time/l");
  events->Branch("this_poll_start",           &this_poll_start,             "this_poll_start/l");
  events->Branch("this_poll_end",             &this_poll_end,               "this_poll_end/l");
  events->Branch("last_poll_start",           &last_poll_start,             "last_poll_start/l");
  events->Branch("last_poll_end",             &last_poll_end,               "last_poll_end/l");
  events->Branch("system_clock_deviation",    &system_clock_deviation,      "system_clock_deviation/I");
  events->Branch("feb_events_per_poll",       &feb_events_per_poll,         "feb_events_per_poll/i");
  events->Branch("feb_event_number",          &feb_event_number,            "feb_event_number/i");

  events->Branch("sequence_id",               &sequence_id,                 "sequence_id/i");
  events->Branch("fragment_timestamp",        &fragment_timestamp,          "fragment_timestamp/l");
}

sbndaq::BernCRTAna::~BernCRTAna()
{
}

void sbndaq::BernCRTAna::analyze_fragment(artdaq::Fragment & frag) {

  BernCRTFragment bern_fragment(frag);

  fragment_timestamp        = frag.timestamp();
  sequence_id               = frag.sequenceID();

  //event data
  BernCRTEvent const* bevt = bern_fragment.eventdata();

//  TLOG(TLVL_INFO)<<*bevt;

  mac5     = bevt->MAC5();
  flags    = bevt->flags;
  lostcpu  = bevt->lostcpu;
  lostfpga = bevt->lostfpga;
  ts0      = bevt->Time_TS0();
  ts1      = bevt->Time_TS1();
  coinc    = bevt->coinc;

  for(int ch=0; ch<32; ch++) adc[ch] = bevt->ADC(ch);

  //metadata
  const BernCRTFragmentMetadata* md = bern_fragment.metadata();
//  TLOG(TLVL_INFO)<<*md;

  run_start_time            = md->run_start_time();
  this_poll_start           = md->this_poll_start();
  this_poll_end             = md->this_poll_end();
  last_poll_start           = md->last_poll_start();
  last_poll_end             = md->last_poll_end();
  system_clock_deviation    = md->system_clock_deviation();
  feb_events_per_poll       = md->feb_events_per_poll();
  feb_event_number          = md->feb_event_number();

  events->Fill();
}


void sbndaq::BernCRTAna::analyze(art::Event const & evt)
{
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  evt.getManyByType(fragmentHandles);
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;
    
    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      //Container fragment
      for (auto cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() != sbndaq::detail::FragmentType::BERNCRT)
          continue;
        for (size_t ii = 0; ii < contf.block_count(); ++ii)
          analyze_fragment(*contf[ii].get());
      }
    }
    else {
      //normal fragment
      if (handle->front().type() != sbndaq::detail::FragmentType::BERNCRT) continue;
      for (auto frag : *handle) 
        analyze_fragment(frag);
    }
  }
} //analyze

DEFINE_ART_MODULE(sbndaq::BernCRTAna)

