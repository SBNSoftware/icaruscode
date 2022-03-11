#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTTranslator.hh"

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
  TTree * hits;
  bool     IsSideCRT(icarus::crt::BernCRTTranslator & hit);
  bool     IsTopCRT(icarus::crt::BernCRTTranslator & hit);
//data payload
  uint8_t mac5; //last 8 bits of FEB mac5 address
  uint16_t flags;
  uint16_t lostcpu;
  uint16_t lostfpga;
  uint32_t ts0;
  uint32_t ts1;
  uint16_t adc[32];
  uint32_t coinc;
  int subSys;  // Top CRT 0, Side CRT 1, Bottom 2


//metadata
  uint64_t  run_start_time;
  uint64_t  this_poll_start;
  uint64_t  this_poll_end;
  uint64_t  last_poll_start;
  uint64_t  last_poll_end;
  int32_t   system_clock_deviation;
  uint32_t  feb_hits_in_poll;
  uint32_t  feb_hit_number;

//information from fragment header
  uint32_t  sequence_id;
  uint64_t  fragment_timestamp;
  uint64_t  last_accepted_timestamp;
  uint16_t  lost_hits;
  uint16_t  hits_in_fragment;
};

//Define the constructor
sbndaq::BernCRTAna::BernCRTAna(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{

  art::ServiceHandle<art::TFileService> tfs; //pointer to a file named tfs

  hits = tfs->make<TTree>("hits", "FEB hits");

  hits->Branch("mac5",        &mac5,          "mac5/b");
  hits->Branch("flags",       &flags,         "flags/s");
  hits->Branch("lostcpu",     &lostcpu,       "lostcpu/s");
  hits->Branch("lostfpga",    &lostfpga,      "lostfpga/s");
  hits->Branch("ts0",         &ts0,           "ts0/i");
  hits->Branch("ts1",         &ts1,           "ts1/i");
  hits->Branch("adc",         &adc,           "adc[32]/s");
  hits->Branch("coinc",       &coinc,         "coinc/i");
  hits->Branch("subSys",         &subSys,           "subSys/i");

  hits->Branch("run_start_time",            &run_start_time,              "run_start_time/l");
  hits->Branch("this_poll_start",           &this_poll_start,             "this_poll_start/l");
  hits->Branch("this_poll_end",             &this_poll_end,               "this_poll_end/l");
  hits->Branch("last_poll_start",           &last_poll_start,             "last_poll_start/l");
  hits->Branch("last_poll_end",             &last_poll_end,               "last_poll_end/l");
  hits->Branch("system_clock_deviation",    &system_clock_deviation,      "system_clock_deviation/I");
  hits->Branch("feb_hits_in_poll",          &feb_hits_in_poll,            "feb_hits_in_poll/i");
  hits->Branch("feb_hit_number",            &feb_hit_number,              "feb_hit_number/i");

  hits->Branch("sequence_id",               &sequence_id,                 "sequence_id/i");
  hits->Branch("fragment_timestamp",        &fragment_timestamp,          "fragment_timestamp/l");

  hits->Branch("last_accepted_timestamp",   &last_accepted_timestamp,      "last_accepted_timestamp/l");
  hits->Branch("lost_hits",                 &lost_hits,                    "lost_hits/s");
  hits->Branch("hits_in_fragment",          &hits_in_fragment,             "hits_in_fragment/s");

}

bool sbndaq::BernCRTAna::IsSideCRT(icarus::crt::BernCRTTranslator & hit) {
  /**
 *    * Fragment ID described in SBN doc 16111
 *       */
 return (hit.fragment_ID & 0x3100) == 0x3100;
}
bool sbndaq::BernCRTAna::IsTopCRT(icarus::crt::BernCRTTranslator & hit) {
  /**
 *  *    * Fragment ID described in SBN doc 16111
 *   *       */
  return (hit.fragment_ID & 0x3200) == 0x3200;
}

sbndaq::BernCRTAna::~BernCRTAna()
{
}

void sbndaq::BernCRTAna::analyze(art::Event const & evt) {

  //WK 09/02/21. Update to BernCRTTranslator in sbndaq_artdaq_core
  std::vector<icarus::crt::BernCRTTranslator> hit_vector;

  auto fragmentHandles = evt.getMany<artdaq::Fragments>();
  for (auto  handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;

    auto this_hit_vector = icarus::crt::BernCRTTranslator::getCRTData(*handle);

    hit_vector.insert(hit_vector.end(),this_hit_vector.begin(),this_hit_vector.end());

  }

  for(auto & hit : hit_vector) {
   // TLOG(TLVL_INFO)<<hit;
    
    if(IsSideCRT(hit)){
      subSys    = 1;
    }
    else if(IsTopCRT(hit)){
      subSys    = 0;
    }
    else{
      subSys    = 2;
    }
    fragment_timestamp        = hit.timestamp;
    sequence_id               = hit.sequence_id;

    mac5     = hit.mac5;
    flags    = hit.flags;
    lostcpu  = hit.lostcpu;
    lostfpga = hit.lostfpga;
    ts0      = hit.ts0;
    ts1      = hit.ts1;
    coinc    = hit.coinc;

    for(int ch=0; ch<32; ch++) adc[ch] = hit.adc[ch];

    run_start_time            = hit.run_start_time;
    this_poll_start           = hit.this_poll_start;
    this_poll_end             = hit.this_poll_end;
    last_poll_start           = hit.last_poll_start;
    last_poll_end             = hit.last_poll_end;
    system_clock_deviation    = hit.system_clock_deviation;
    feb_hits_in_poll          = hit.hits_in_poll;
    feb_hit_number            = hit.feb_hit_number; 

    last_accepted_timestamp   = hit.last_accepted_timestamp;
    lost_hits                 = hit.lost_hits;
    hits_in_fragment          = hit.hits_in_fragment;

    hits->Fill();
  }
} //analyze


DEFINE_ART_MODULE(sbndaq::BernCRTAna)

