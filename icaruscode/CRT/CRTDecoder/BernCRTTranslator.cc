//#include "art/Framework/Core/EDAnalyzer.h"
//#include "art/Framework/Core/ModuleMacros.h"
//#include "art/Framework/Principal/Event.h"
//#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTZMQFragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"

#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"

//#include <algorithm>
//#include <cassert>
//#include <cmath>
//#include <fstream>
//#include <iomanip>
#include <vector>
#include <iostream>

#include "BernCRTTranslator.hh"

/* void sbndaq::BernCRTAna::analyze_fragment(artdaq::Fragment & frag) {

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

} */

std::vector<icarus::crt::BernCRTTranslator> icarus::crt::BernCRTTranslator::getCRTData(art::Event const & evt) {
  std::vector<icarus::crt::BernCRTTranslator> out;

  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  evt.getManyByType(fragmentHandles);
/*  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;
    
    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      //Container fragment
      for (auto cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() != sbndaq::detail::FragmentType::BERNCRT)
          continue;
        for (size_t ii = 0; ii < contf.block_count(); ++ii) {
        //  analyze_fragment(*contf[ii].get());
        }
      }
    }
    else {
      //normal fragment
      if (handle->front().type() != sbndaq::detail::FragmentType::BERNCRT) continue;
      for (auto frag : *handle) {
      //  analyze_fragment(frag);
      }
    }
  } */

  return out;
}


