#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTZMQFragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragment.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"

#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"

#include <iostream>
#include <bitset>

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
  
  out.hits_in_fragment          = 1;

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

  out.hits_in_fragment          = 1;

  return out;
}

std::vector<icarus::crt::BernCRTTranslator> icarus::crt::BernCRTTranslator::analyze_BernCRTFragmentV2(artdaq::Fragment & frag) {
  
  const sbndaq::BernCRTFragmentV2 bern_fragment(frag);
  const sbndaq::BernCRTFragmentMetadataV2* md = bern_fragment.metadata();
//  TLOG(TLVL_INFO)<<*bern_fragment;

  const unsigned int nhits      = md->hits_in_fragment();
  std::vector<icarus::crt::BernCRTTranslator> OUT;
  OUT.reserve(nhits);

  for(unsigned int iHit = 0; iHit < nhits; iHit++) {
    icarus::crt::BernCRTTranslator out;

    const sbndaq::BernCRTHitV2* bevt = bern_fragment.eventdata(iHit);

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
 
    out.hits_in_fragment          = nhits;

    OUT.push_back(out);
  }
  return OUT;
} 

std::vector<icarus::crt::BernCRTTranslator> icarus::crt::BernCRTTranslator::getCRTData(art::Event const & evt) {
  std::vector<icarus::crt::BernCRTTranslator> out;

  //std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  //evt.getManyByType(fragmentHandles);
  auto fragmentHandles = evt.getMany<artdaq::Fragments>();
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

std::ostream & icarus::crt::operator << (std::ostream & os, icarus::crt::BernCRTTranslator const& t){
  os <<"\n\tMAC5:              0x" << std::hex << (int)t.mac5 << std::dec << " (" << (int)t.mac5 << ")"
     << "\n\tRun start time:    " << sbndaq::BernCRTFragment::print_timestamp(t.run_start_time)
     << "\n\tThis poll start:   " << sbndaq::BernCRTFragment::print_timestamp(t.this_poll_start)
     << "\n\tThis poll finish:  " << sbndaq::BernCRTFragment::print_timestamp(t.this_poll_end)
     << "\n\tLast poll start:   " << sbndaq::BernCRTFragment::print_timestamp(t.last_poll_start)
     << "\n\tLast poll finish:  " << sbndaq::BernCRTFragment::print_timestamp(t.last_poll_end)
     << "\n\tClock deviation:   " << t.system_clock_deviation<<" ns"
     << "\n\tFEB hits/poll:     " << t.hits_in_poll
     << "\n\tFEB hits/fragment: " << t.hits_in_fragment
     << "\n\tFlags:             " << std::bitset<8>(t.flags)
     <<(t.IsOverflow_TS0() ?" [T0 overflow]" :"")
     <<(t.IsOverflow_TS1() ?" [T1 overflow]" :"")
     <<(t.IsReference_TS0()?" [T0 reference]":"")
     <<(t.IsReference_TS1()?" [T1 reference]":"")
     << "\n\tLostCPU:           " << t.lostcpu
     << "\n\tLostFPGA:          " << t.lostfpga
     << "\n\tTime1 (TS0):       " << sbndaq::BernCRTFragment::print_timestamp(t.ts0)
     << "\n\tTime2 (TS1):       " << sbndaq::BernCRTFragment::print_timestamp(t.ts1);

  os << "\n\t[#ch]: ADC  ";
  for(size_t i_c=0; i_c<32; ++i_c) {
    if(!(i_c % 8)) os<<"\n\t";
    os << " ["<<std::setw(2)<<i_c<<"]: " <<std::setw(4)<< t.adc[i_c];
  }

  os << "\n\tCoincidence:       " << std::bitset<32>(t.coinc)
     << "\n\tTimestamp:         " << sbndaq::BernCRTFragment::print_timestamp(t.timestamp)
     << "\n\tFEB hit number:    " << t.feb_hit_number
     << "\n\tLost hits:         " << t.lost_hits
     << "\n\tLast timestamp:    " << sbndaq::BernCRTFragment::print_timestamp(t.last_accepted_timestamp)
     << "\n";

  return os;
}

