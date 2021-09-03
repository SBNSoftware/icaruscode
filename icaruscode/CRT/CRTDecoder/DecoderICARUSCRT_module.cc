////////////////////////////////////////////////////////////////////////
// Class:       DecoderICARUSCRT
// Plugin Type: producer (art v3_06_03)
// File:        DecoderICARUSCRT_module.cc
//
// Generated at Sat May  1 20:19:33 2021 by Biswaranjan Behera using cetskelgen
// from cetlib version v3_11_01.
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

  // Declare member data here.
  const icarusDB::IICARUSChannelMap* fChannelMap = nullptr;
};


crt::DecoderICARUSCRT::DecoderICARUSCRT(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
// More initializers here.

{
  fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
  // Call appropriate produces<>() functions here.
  produces< std::vector<icarus::crt::CRTData> >();
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void crt::DecoderICARUSCRT::produce(art::Event& evt)
{
  // Implementation of required member function here.
  //  std::unique_ptr< std::vector<icarus::crt::CRTData> > crtdata( new std::vector<icarus::crt::CRTData>);
  auto crtdata = std::make_unique<std::vector<icarus::crt::CRTData>>();
  
  //WK 09/02/21. Update to BernCRTTranslator in sbndaq_artdaq_core
  std::vector<icarus::crt::BernCRTTranslator> hit_vector;

  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  evt.getManyByType(fragmentHandles);
  for (auto  handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;

    auto this_hit_vector = icarus::crt::BernCRTTranslator::getCRTData(*handle);

    hit_vector.insert(hit_vector.end(),this_hit_vector.begin(),this_hit_vector.end());

  }

  for (auto & hit : hit_vector){ 

    icarus::crt::CRTData data;
    data.fMac5  = fChannelMap->getSimMacAddress(hit.mac5);
    data.fTs0   = hit.ts0;
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
