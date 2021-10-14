////////////////////////////////////////////////////////////////////////
// Class:       DecoderICARUSCRT
// Plugin Type: producer (art v3_06_03)
// File:        DecoderICARUSCRT_module.cc
//
// Generated at Sat May  1 20:19:33 2021 by Biswaranjan Behera using cetskelgen
// from cetlib version v3_11_01.
//
// Thanks to Gianluca Petrillo for helping me on improving the decoder 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTTranslator.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

#include "sbnobj/ICARUS/CRT/CRTData.hh"

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "TH1F.h"
#include "TNtuple.h"

#include <memory>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include<stdlib.h>
#include <map>

namespace crt {
  class DecoderICARUSCRT;
}


class crt::DecoderICARUSCRT : public art::EDProducer {
public:
  explicit DecoderICARUSCRT(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DecoderICARUSCRT(DecoderICARUSCRT const&) = delete;
  DecoderICARUSCRT(DecoderICARUSCRT&&) = delete;
  DecoderICARUSCRT& operator=(DecoderICARUSCRT const&) = delete;
  DecoderICARUSCRT& operator=(DecoderICARUSCRT&&) = delete;

  // Required functions.
  void produce(art::Event& evt) override;

private:
  uint64_t CalculateTimestamp(icarus::crt::BernCRTTranslator& hit);

  // Declare member data here.
  const icarusDB::IICARUSChannelMap* fChannelMap = nullptr;

  std::map<uint8_t, int32_t> FEB_delay; //<mac5, delay in ns>
};


crt::DecoderICARUSCRT::DecoderICARUSCRT(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
  produces< std::vector<icarus::crt::CRTData> >();

  std::vector<std::vector<int32_t> > delays =  p.get<std::vector<std::vector<int32_t> > >("FEB_delay");
  for(auto & feb : delays) {
    int32_t & mac = feb[0];
    int32_t & d   = feb[1];
    FEB_delay[mac] = d;
    std::cout<<"Read delay for mac5 "<<std::setw(3)<<(int)mac<<": "<<std::setw(4)<<d<<" ns\n";
  }
}

uint64_t crt::DecoderICARUSCRT::CalculateTimestamp(icarus::crt::BernCRTTranslator& hit) {
  /**
   * Calculate timestamp based on nanosecond from FEB and poll times measured by server
   * see: https://sbn-docdb.fnal.gov/cgi-bin/private/DisplayMeeting?sessionid=7783
   */
  int32_t ts0  = hit.ts0; //must be signed int

  //add PPS cable length offset modulo 1s
  ts0 = (ts0 + FEB_delay.at(hit.mac5)) % (1'000'000'000);
  if(ts0 < 0) ts0 += 1000'000'000; //just in case the cable offset is negative (should be positive normally)

  uint64_t mean_poll_time = hit.last_poll_start/2 + hit.this_poll_end/2;
  int mean_poll_time_ns = mean_poll_time % (1000'000'000); 
  
  return mean_poll_time - mean_poll_time_ns + ts0
    + (ts0 - mean_poll_time_ns < -500'000'000) * 1000'000'000
    - (ts0 - mean_poll_time_ns >  500'000'000) * 1000'000'000;
}

void crt::DecoderICARUSCRT::produce(art::Event& evt)
{

  // Implementation of required member function here.
  //  std::unique_ptr< std::vector<icarus::crt::CRTData> > crtdata( new std::vector<icarus::crt::CRTData>);
  auto crtdata = std::make_unique<std::vector<icarus::crt::CRTData>>();
  
  //WK 09/02/21. Update to BernCRTTranslator in sbndaq_artdaq_core
  std::vector<icarus::crt::BernCRTTranslator> hit_vector;

  auto fragmentHandles = evt.getMany<artdaq::Fragments>();
  for (auto  handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;

    auto this_hit_vector = icarus::crt::BernCRTTranslator::getCRTData(*handle);

    hit_vector.insert(hit_vector.end(),this_hit_vector.begin(),this_hit_vector.end());

  }

  for (auto & hit : hit_vector){ 

    icarus::crt::CRTData data;
    data.fMac5  = fChannelMap->getSimMacAddress(hit.mac5);
    data.fTs0   = CalculateTimestamp(hit);
    data.fTs1   = hit.ts1;
    // data.fEntry = hit.entry;
    //data.coinc    = hit.coinc;

    for(int ch=0; ch<32; ch++) {

      // East-center wall
      if (hit.mac5 == 27 || hit.mac5 == 28
          || hit.mac5 == 33 || hit.mac5 == 34){
        if (ch > 19) {
          data.fAdc[ch] = 0.;
        }else {data.fAdc[ch] = hit.adc[ch+2];}

	// West-center wall
      }else if (hit.mac5 == 21 || hit.mac5 == 22
                || hit.mac5 == 13 || hit.mac5 == 14){
        if (ch > 19) {
          data.fAdc[ch] = 0.;
        }else {data.fAdc[ch] = hit.adc[ch+2];}
	// All side-crt walls
      } else  if (ch > 29){
	data.fAdc[ch] = 0.;
      }else {
	data.fAdc[ch] = hit.adc[ch+2];
      }
    }
    crtdata->push_back(std::move(data));
  }
  
  evt.put(std::move(crtdata));
 
}

DEFINE_ART_MODULE(crt::DecoderICARUSCRT)
