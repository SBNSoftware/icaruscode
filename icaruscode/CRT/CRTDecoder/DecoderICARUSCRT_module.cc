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
  : EDProducer{p}
{
  fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
  produces< std::vector<icarus::crt::CRTData> >();
  
}

void crt::DecoderICARUSCRT::produce(art::Event& evt)
{

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
  
  struct Recipe_t {
    
    unsigned int destMac5;
    unsigned int firstSourceChannel;
    unsigned int lastSourceChannel;
    unsigned int firstDestChannel;
    unsigned int lastDestChannel;
    
    int direction; // +1 or -1
    
  };

  // vector: Mac5 -> its CRT data
  std::vector<icarus::crt::CRTData> allCRTdata ( 305 + 1); // TODO size this correctly!
  
  for (auto & hit : hit_vector){ 

    std::array<Recipe_t, 3U> allRecipes;
    //
    // fill the recipe
    //
    if (!((hit.mac5 >= 88 && hit.mac5 <= 91)
	  || hit.mac5 == 96 || hit.mac5 == 97
	  || hit.mac5 ==  1 || hit.mac5 ==  3
	  || hit.mac5 ==  6 || hit.mac5 ==  7)) { // look for FEB those are not between 88 to 91 
      
      int const destMac5 = fChannelMap->getSimMacAddress(hit.mac5);
      
      Recipe_t recipe;
      
      // 
      // first block of 10 channels from source
      // 
      recipe.destMac5 = destMac5;
      
      recipe.firstSourceChannel =  2;
      recipe.lastSourceChannel  = 11;
      
      recipe.firstDestChannel   =  0;
      recipe.lastDestChannel    =  9;
      recipe.direction          = +1;
      allRecipes[0] = recipe;
      
      // 
      // second block of 10 channels from source
      // 
      recipe.destMac5 = destMac5;
      recipe.firstSourceChannel = 12;
      recipe.lastSourceChannel  = 21;
      
      recipe.firstDestChannel   = 10;
      recipe.lastDestChannel    = 19;
      recipe.direction          = +1;
      allRecipes[1] = recipe;
      
      // 
      // third block of 10 channels from source
      // 
      recipe.destMac5 = destMac5;
      recipe.firstSourceChannel  = 22;
      recipe.lastSourceChannel   = 31;
      
      recipe.firstDestChannel    = 20;
      recipe.lastDestChannel     = 29;
      recipe.direction           = +1;
      allRecipes[2] = recipe;
      
      
    } // "normal assignment"
    else if (hit.mac5 ==  97) { // south wall - east side top horizontal module channels are reversed
      
      int const destMac5 = fChannelMap->getSimMacAddress(hit.mac5);
      
      Recipe_t recipe;
      
      // 
      // first block of 10 channels from source
      // 
      recipe.destMac5 = destMac5;
      
      recipe.firstSourceChannel =  2;
      recipe.lastSourceChannel  = 11;
      
      recipe.firstDestChannel   =  0;
      recipe.lastDestChannel    =  9;
      recipe.direction          = +1;
      allRecipes[0] = recipe;
      
      // 
      // second block of 10 channels from source
      // 
      recipe.destMac5 = destMac5;
      recipe.firstSourceChannel = 12;
      recipe.lastSourceChannel  = 21;
      
      recipe.firstDestChannel   = 10;
      recipe.lastDestChannel    = 19;
      recipe.direction          = +1;
      allRecipes[1] = recipe;
      
      // 
      // third block of 10 channels from source
      // 
      recipe.destMac5 = destMac5;
      recipe.firstSourceChannel  = 22;
      recipe.lastSourceChannel   = 31;
      
      recipe.firstDestChannel    = 29;
      recipe.lastDestChannel     = 20;
      recipe.direction           = -1;
      allRecipes[2] = recipe;
      
      
    } 
    else if (hit.mac5 == 1 || hit.mac5 == 3 ||
	     hit.mac5 == 6 || hit.mac5 == 7 ||
	     hit.mac5 == 96) { // north wall inner layer and south wall west side top three horizontal layer orientation is reversed 

      int const destMac5 = fChannelMap->getSimMacAddress(hit.mac5);

      Recipe_t recipe;

      //
      // first block of 10 channels from source
      //
      recipe.destMac5 = destMac5;

      recipe.firstSourceChannel  =  2;
      recipe.lastSourceChannel   = 11;

      recipe.firstDestChannel    =  9;
      recipe.lastDestChannel     =  0;
      recipe.direction           = -1;
      allRecipes[0] = recipe;

      //
      // second block of 10 channels from source
      //
      recipe.destMac5 = destMac5;
      recipe.firstSourceChannel = 12;
      recipe.lastSourceChannel  = 21;

      recipe.firstDestChannel   = 19;
      recipe.lastDestChannel    = 10;
      recipe.direction          = -1;
      allRecipes[1] = recipe;

      //
      // third block of 10 channels from source: special mapping
      //
      recipe.destMac5 = destMac5;
      recipe.firstSourceChannel = 22;
      recipe.lastSourceChannel  = 31;

      recipe.firstDestChannel   = 29;
      recipe.lastDestChannel    = 20;
      recipe.direction          = -1;
      allRecipes[2] = recipe;

    }
    else if (hit.mac5 == 88) {
      
      int const destMac5 = fChannelMap->getSimMacAddress(hit.mac5);
      
      Recipe_t recipe;
      
      // 
      // first block of 10 channels from source
      // 
      recipe.destMac5 = 79; 
      
      recipe.firstSourceChannel  =  2;
      recipe.lastSourceChannel   = 11;
      
      recipe.firstDestChannel    =  29;
      recipe.lastDestChannel     =  20;
      recipe.direction           =  -1;
      allRecipes[0] = recipe;
      
      // 
      // second block of 10 channels from source
      // 
      recipe.destMac5 = destMac5; 
      recipe.firstSourceChannel = 12;
      recipe.lastSourceChannel  = 21;
      
      recipe.firstDestChannel   = 10;
      recipe.lastDestChannel    = 19;
      recipe.direction          = +1;
      allRecipes[1] = recipe;
      
      // 
      // third block of 10 channels from source: special mapping
      // 
      recipe.destMac5 = destMac5; 
      recipe.firstSourceChannel = 22;
      recipe.lastSourceChannel  = 31;
      
      recipe.firstDestChannel   =  0;
      recipe.lastDestChannel    =  9;
      recipe.direction          = +1;
      allRecipes[2] = recipe;
      
    } 
    else if ((hit.mac5 >= 89) && (hit.mac5 <= 91)) {
      
      int const destMac5 = fChannelMap->getSimMacAddress(hit.mac5); 
      
      Recipe_t recipe;
      
      // 
      // first block of 10 channels from source
      // 
      recipe.destMac5 = destMac5; 
      
      recipe.firstSourceChannel  =  2;
      recipe.lastSourceChannel   = 11;
      
      recipe.firstDestChannel    = 19;
      recipe.lastDestChannel     = 10;
      recipe.direction           = -1;
      allRecipes[0] = recipe;
      
      // 
      // second block of 10 channels from source
      // 
      recipe.destMac5 = destMac5; 
      recipe.firstSourceChannel = 12;
      recipe.lastSourceChannel  = 21;
      
      recipe.firstDestChannel   =  9;
      recipe.lastDestChannel    =  0;
      recipe.direction          = -1;
      allRecipes[1] = recipe;
      
      // 
      // third block of 10 channels from source: special mapping
      // 
      recipe.destMac5 = destMac5 - 1; 
      recipe.firstSourceChannel  = 22;
      recipe.lastSourceChannel   = 31;
      
      recipe.firstDestChannel    = 29;
      recipe.lastDestChannel     = 20;
      recipe.direction           = -1;
      
      allRecipes[2] = recipe;
      
    } // if not 88
    
      
    //
    // cook the crtdata
    //
    for (Recipe_t const& recipe: allRecipes) {
      if (recipe.firstSourceChannel == recipe.lastSourceChannel) continue;
      
      icarus::crt::CRTData& data = allCRTdata.at(recipe.destMac5);
      data.fMac5  = recipe.destMac5;
      //data.fTs0   = hit.ts0;
      data.fTs0   = hit.timestamp;
      data.fTs1   = hit.ts1;
      

      unsigned destCh = recipe.firstDestChannel;
      for (unsigned srcCh = recipe.firstSourceChannel; srcCh <= recipe.lastSourceChannel; ++srcCh) {

        data.fAdc[destCh] = hit.adc[srcCh];
        destCh += recipe.direction; // increase or decrease the source

      }
      
    } // for all recipes
    
  } // for all input data
  
  // move the data which is actually present in the final data product
  auto crtdata = std::make_unique<std::vector<icarus::crt::CRTData>>();
  for (icarus::crt::CRTData& crtDataElem: allCRTdata) {
    if (crtDataElem.fMac5 == 0) continue; // not a valid Mac5, data is not present
    crtdata->push_back(std::move(crtDataElem));
  }
    
  evt.put(std::move(crtdata));
 
}

DEFINE_ART_MODULE(crt::DecoderICARUSCRT)
