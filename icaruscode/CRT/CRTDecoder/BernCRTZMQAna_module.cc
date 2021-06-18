////////////////////////////////////////////////////////////////////////
// Class:       BernCRTZMQAna
// Module Type: analyzer
// File:        BernCRTZMQAna_module.cc
// Description: Makes a tree with waveform information.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTZMQFragment.hh"
#include "artdaq-core/Data/Fragment.hh"

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
  class BernCRTZMQAna;
}

class sbndaq::BernCRTZMQAna : public art::EDAnalyzer {

public:
  explicit BernCRTZMQAna(fhicl::ParameterSet const & pset); // explicit doesn't allow for copy initialization
  virtual ~BernCRTZMQAna();
  
  virtual void analyze(art::Event const & evt);
  

private:
  
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
sbndaq::BernCRTZMQAna::BernCRTZMQAna(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{

  //this is how you setup/use the TFileService
  //I do this here in the constructor to setup the histogram just once
  //but you can call the TFileService and make new objects from anywhere
  art::ServiceHandle<art::TFileService> tfs; //pointer to a file named tfs

  events = tfs->make<TTree>("events", "FEB events");

//event data
  events->Branch("mac5",        &mac5,          "mac5/b");
  events->Branch("flags",       &flags,         "flags/s");
  events->Branch("lostcpu",     &lostcpu,       "lostcpu/s");
  events->Branch("lostfpga",    &lostcpu,       "lostfpga/s");
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

sbndaq::BernCRTZMQAna::~BernCRTZMQAna()
{
}

void sbndaq::BernCRTZMQAna::analyze(art::Event const & evt)
{


  //can get the art event number
  //art::EventNumber_t eventNumber = evt.event();
  
  //we get a 'handle' to the fragments in the event
  //this will act like a pointer to a vector of artdaq fragments
  art::Handle< std::vector<artdaq::Fragment> > rawFragHandle;
  
  //we fill the handle by getting the right data according to the label
  //the module label will be 'daq', while the instance label (second argument) matches the type of fragment
  evt.getByLabel("daq","BERNCRTZMQ", rawFragHandle);

  //this checks to make sure it's ok
  if (!rawFragHandle.isValid()) return;


  for (size_t idx = 0; idx < rawFragHandle->size(); ++idx) { // loop over the fragments of an event

    const auto& frag((*rawFragHandle)[idx]); // use this fragment as a reference to the same data

    //this applies the 'overlay' to the fragment, to tell it to treat it like a BernCRTZMQ fragment
    BernCRTZMQFragment bern_fragment(frag);
    
    fragment_timestamp        = frag.timestamp();
    sequence_id               = frag.sequenceID();

    //event data
    BernCRTZMQEvent const* bevt = bern_fragment.eventdata();

    mac5     = bevt->MAC5();
    flags    = bevt->flags;
    lostcpu  = bevt->lostcpu;
    lostfpga = bevt->lostfpga;
    ts0      = bevt->Time_TS0();
    ts1      = bevt->Time_TS1();
    coinc    = bevt->coinc;

    for(int ch=0; ch<32; ch++) adc[ch] = bevt->ADC(ch);

    //metadata
    const BernCRTZMQFragmentMetadata* md = bern_fragment.metadata();

    run_start_time            = md->run_start_time();
    this_poll_start           = md->this_poll_start();
    this_poll_end             = md->this_poll_end();
    last_poll_start           = md->last_poll_start();
    last_poll_end             = md->last_poll_end();
    system_clock_deviation    = md->system_clock_deviation();
    feb_events_per_poll       = md->feb_events_per_poll();
    feb_event_number          = md->feb_event_number();

    events->Fill();

  }//end loop over fragments
}

DEFINE_ART_MODULE(sbndaq::BernCRTZMQAna)
//this is where the name is specified
