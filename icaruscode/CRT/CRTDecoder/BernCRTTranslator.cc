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
#include <iostream>

#include "BernCRTTranslator.hh"

icarus::crt::BernCRTTranslator icarus::crt::BernCRTTranslator::analyze_BernCRTZMQFragment(artdaq::Fragment & frag) {
  icarus::crt::BernCRTTranslator out;

  const sbndaq::BernCRTZMQFragment bern_fragment(frag);

  out.timestamp        = frag.timestamp();
  out.sequence_id      = frag.sequenceID();

  //event data
  const sbndaq::BernCRTZMQEvent* bevt = bern_fragment.eventdata();

//  TLOG(TLVL_INFO)<<*bevt;

  out.mac5     = bevt->MAC5();
  out.flags    = bevt->flags;
  out.lostcpu  = bevt->lostcpu;
  out.lostfpga = bevt->lostfpga;
  out.ts0      = bevt->Time_TS0();
  out.ts1      = bevt->Time_TS1();
  out.coinc    = bevt->coinc;

  for(int ch=0; ch<32; ch++) out.adc[ch] = bevt->ADC(ch);

  //metadata
  const sbndaq::BernCRTZMQFragmentMetadata* md = bern_fragment.metadata();
//  TLOG(TLVL_INFO)<<*md;

  out.run_start_time            = md->run_start_time();
  out.this_poll_start           = md->this_poll_start();
  out.this_poll_end             = md->this_poll_end();
  out.last_poll_start           = md->last_poll_start();
  out.last_poll_end             = md->last_poll_end();
  out.system_clock_deviation    = md->system_clock_deviation();
  out.hits_in_poll              = md->feb_events_per_poll();
  out.feb_hit_number            = md->feb_event_number();
  out.lost_hits                 = md->omitted_events();
  out.last_accepted_timestamp   = md->last_accepted_timestamp();
  return out;
}

icarus::crt::BernCRTTranslator icarus::crt::BernCRTTranslator::analyze_BernCRTFragment(artdaq::Fragment & frag) {
  icarus::crt::BernCRTTranslator out;

  const sbndaq::BernCRTFragment bern_fragment(frag);

  out.timestamp        = frag.timestamp();
  out.sequence_id      = frag.sequenceID();

  //event data
  const sbndaq::BernCRTEvent* bevt = bern_fragment.eventdata();

//  TLOG(TLVL_INFO)<<*bevt;

  out.mac5     = bevt->MAC5();
  out.flags    = bevt->flags;
  out.lostcpu  = bevt->lostcpu;
  out.lostfpga = bevt->lostfpga;
  out.ts0      = bevt->Time_TS0();
  out.ts1      = bevt->Time_TS1();
  out.coinc    = bevt->coinc;

  for(int ch=0; ch<32; ch++) out.adc[ch] = bevt->ADC(ch);

  //metadata
  const sbndaq::BernCRTFragmentMetadata* md = bern_fragment.metadata();
//  TLOG(TLVL_INFO)<<*md;

  out.run_start_time            = md->run_start_time();
  out.this_poll_start           = md->this_poll_start();
  out.this_poll_end             = md->this_poll_end();
  out.last_poll_start           = md->last_poll_start();
  out.last_poll_end             = md->last_poll_end();
  out.system_clock_deviation    = md->system_clock_deviation();
  out.hits_in_poll              = md->feb_events_per_poll();
  out.feb_hit_number            = md->feb_event_number();
  out.lost_hits                 = md->omitted_events();
  out.last_accepted_timestamp   = md->last_accepted_timestamp();
  return out;
}

std::vector<icarus::crt::BernCRTTranslator> icarus::crt::BernCRTTranslator::analyze_BernCRTFragmentV2(artdaq::Fragment & frag) {
  
  const sbndaq::BernCRTFragmentV2 bern_fragment(frag);
  const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
//  TLOG(TLVL_INFO)<<*md;

  unsigned int feb_hits_in_fragment      = md->hits_in_fragment();
  std::vector<icarus::crt::BernCRTTranslator> OUT(feb_hits_in_fragment);

  for(unsigned int iHit = 0; iHit < feb_hits_in_fragment; iHit++) {
    icarus::crt::BernCRTTranslator out;

    const sbndaq::BernCRTHitV2* bevt = bern_fragment.eventdata(iHit);
    //  TLOG(TLVL_INFO)<<*bevt;

    out.flags                     = bevt->flags;
    out.lostcpu                   = bevt->lostcpu;
    out.lostfpga                  = bevt->lostfpga;
    out.ts0                       = bevt->ts0;
    out.ts1                       = bevt->ts1;
    out.coinc                     = bevt->coinc;
    out.feb_hit_number            = bevt->feb_hit_number;
    out.timestamp                 = bevt->timestamp;
    out.last_accepted_timestamp   = bevt->last_accepted_timestamp;
    out.lost_hits                 = bevt->lost_hits;

    for(int ch=0; ch<32; ch++) out.adc[ch] = bevt->adc[ch];

    out.sequence_id               = frag.sequenceID();

    //metadata
    out.mac5                      = md->MAC5();
    out.run_start_time            = md->run_start_time();
    out.this_poll_start           = md->this_poll_start();
    out.this_poll_end             = md->this_poll_end();
    out.last_poll_start           = md->last_poll_start();
    out.last_poll_end             = md->last_poll_end();
    out.system_clock_deviation    = md->system_clock_deviation();
    out.hits_in_poll              = md->hits_in_poll();

    OUT.push_back(out);
  }
  return OUT;
} 

std::vector<icarus::crt::BernCRTTranslator> icarus::crt::BernCRTTranslator::getCRTData(art::Event const & evt) {
  std::vector<icarus::crt::BernCRTTranslator> out;

  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  evt.getManyByType(fragmentHandles);
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;
    
    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      //Container fragment
      for (auto cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        switch(contf.fragment_type()) {
          case sbndaq::detail::FragmentType::BERNCRTZMQ:
            out.reserve(out.size() + contf.block_count());
            for (size_t ii = 0; ii < contf.block_count(); ++ii) {
              out.push_back(analyze_BernCRTZMQFragment(*contf[ii].get()));
            }
            break;
          case sbndaq::detail::FragmentType::BERNCRT:
            out.reserve(out.size() + contf.block_count());
            for (size_t ii = 0; ii < contf.block_count(); ++ii) {
              out.push_back(analyze_BernCRTFragment(*contf[ii].get()));
            }
            break;
          case sbndaq::detail::FragmentType::BERNCRTV2:
            for (size_t ii = 0; ii < contf.block_count(); ++ii) {
              std::vector<icarus::crt::BernCRTTranslator> v = analyze_BernCRTFragmentV2(*contf[ii].get());
              out.insert( out.end(), v.begin(), v.end() );
            }
            break; 
        }
      }
    }
    else {
      //normal fragment
      switch(handle->front().type()) {
        case sbndaq::detail::FragmentType::BERNCRTZMQ:
          out.reserve(out.size() + handle->size());
          for (auto frag : *handle) {
            out.push_back(analyze_BernCRTZMQFragment(frag));
          }
          break;
        case sbndaq::detail::FragmentType::BERNCRT:
          out.reserve(out.size() + handle->size());
          for (auto frag : *handle) {
            out.push_back(analyze_BernCRTFragment(frag));
          }
          break;
        case sbndaq::detail::FragmentType::BERNCRTV2:
          out.reserve(out.size() + handle->size());
          for (auto frag : *handle) {
            std::vector<icarus::crt::BernCRTTranslator> v = analyze_BernCRTFragmentV2(frag);
            out.insert( out.end(), v.begin(), v.end() );
          }
          break; 
      }
    }
  } 

  return out;
}


