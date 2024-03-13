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

#include "icaruscode/Utilities/ArtDataProductSelectors.h"
#include "icaruscode/Utilities/ArtHandleTrackerManager.h"
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
#include <stdlib.h>
#include <map>

//Bottom specific includes
#include "TGraph.h"
#include "Bottom/Fragments/CRTFragment.hh"
#include "Bottom/Data/CRTTrigger.h"

using namespace CRT;
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
  void     CorrectForCableDelay(icarus::crt::BernCRTTranslator & hit);
  bool     IsSideCRT(icarus::crt::BernCRTTranslator & hit);
  bool     IsTopCRT(icarus::crt::BernCRTTranslator & hit);
  

  // Declare member data here.
  const icarusDB::IICARUSChannelMap* fChannelMap = nullptr;
  
  /// Selector object detecting all the suitable data products.
  util::RegexDataProductSelector const fInputTagPatterns;
  bool fDropRawDataAfterUse; ///< Clear fragment data product cache after use.
  
  std::map<uint8_t, int32_t> FEB_delay_side; //<mac5, delay in ns>
  std::map<uint8_t, int32_t> FEB_delay_top;  //<mac5, delay in ns>
  //Bottom variables and functions
  const art::InputTag fFragTag;
  const bool fLookForContainer;
  uint64_t fEarliestTime;
  std::vector<size_t> fChannelMapB; //Simple map from raw data module number to offline module number.
  // Compartmentalize internal functionality so that I can reuse it with both regular Fragments and "container" Fragments
  void FragmentToTriggers(const artdaq::Fragment& artFrag, std::unique_ptr<std::vector<CRT::Trigger>>& triggers);
  // For the first Event of every job, I want to set fEarliestTime to the earliest time in that Event
  void SetEarliestTime(const artdaq::Fragment& frag);
  //Counting variables for sub module event. for debugging
  uint64_t      topCount = 0; 
  uint64_t      sideCount = 0;
  uint64_t      bottomCount = 0;
  uint64_t      totalCount = 0;
};


crt::DecoderICARUSCRT::DecoderICARUSCRT(fhicl::ParameterSet const& p): EDProducer{p},
    fInputTagPatterns{
    util::RegexDataProductSelector::makePatterns(p.get(
      "FragmentTagPatterns",
      std::vector<std::string>{ "daq:(Container)?BERNCRT.*" }
      ))
    }, 
    fDropRawDataAfterUse( p.get<bool>("DropRawDataAfterUse", true) ),
    fFragTag(p.get<std::string>("RawDataTag","hello")), 
    fLookForContainer(p.get<bool>("LookForContainer", false)),    
    fEarliestTime(std::numeric_limits<decltype(fEarliestTime)>::max())

{
std::cout<<"RawDataTag: "<<fFragTag<<'\n';  
fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
  produces< std::vector<icarus::crt::CRTData> >();
  //mayConsumeMany<artdaq::Fragments>();
  consumes<std::vector<artdaq::Fragment>>(fFragTag);
  {
    std::vector<std::vector<int32_t> > delays =  p.get<std::vector<std::vector<int32_t> > >("FEB_delay_side");
    for(auto & feb : delays) {
      int32_t & mac = feb[0];
      int32_t & d   = feb[1];
      FEB_delay_side[mac] = d;
    }
  }

  {
    std::vector<std::vector<int32_t> > delays =  p.get<std::vector<std::vector<int32_t> > >("FEB_delay_top");
    for(auto & feb : delays) {
      int32_t & mac = feb[0];
      int32_t & d   = feb[1];
      FEB_delay_top[mac] = d;
    }
  }
}

bool crt::DecoderICARUSCRT::IsSideCRT(icarus::crt::BernCRTTranslator & hit) {
  /**
   * Fragment ID described in SBN doc 16111
   */
  return (hit.fragment_ID & 0x3100) == 0x3100;
}
bool crt::DecoderICARUSCRT::IsTopCRT(icarus::crt::BernCRTTranslator & hit) {
  /**
   * Fragment ID described in SBN doc 16111
   */
  return (hit.fragment_ID & 0x3200) == 0x3200;
}
void crt::DecoderICARUSCRT::CorrectForCableDelay(icarus::crt::BernCRTTranslator & hit) {
  if(!hit.IsReference_TS0() && !hit.IsReference_TS1()) { //don't correct reference T0 and T1 hits for cable length
    try {
      int32_t delay;
      if(IsSideCRT(hit)) {
        delay = FEB_delay_side.at(hit.mac5);
      }
      else {
        delay = FEB_delay_top.at(hit.mac5);
      }
      hit.ts0 += delay;
      hit.ts0 %= 1'000'000'000;
      if(hit.ts0 < 0) hit.ts0 += 1000'000'000; //just in case the cable offset is negative (should be positive normally)
      hit.ts1 += delay;
    } catch(const std::out_of_range & e) {
      TLOG(TLVL_ERROR)<<"CRT MAC "<<(int)(hit.mac5)<<" not found in the FEB_delay array!!! Please update FEB_delay FHiCL file";
      throw cet::exception("DecoderICARUSCRT")
        << "CRT MAC "<<(int)(hit.mac5)<<" not found in the FEB_delay array!!! Please update FEB_delay FHiCL file\n";
    }
  }
}

uint64_t crt::DecoderICARUSCRT::CalculateTimestamp(icarus::crt::BernCRTTranslator& hit) {
  /**
   * Calculate timestamp based on nanosecond from FEB and poll times measured by server
   * see: https://sbn-docdb.fnal.gov/cgi-bin/private/DisplayMeeting?sessionid=7783
   */
  int32_t ts0  = hit.ts0; //must be signed int

  uint64_t mean_poll_time = hit.last_poll_start/2 + hit.this_poll_end/2;
  int mean_poll_time_ns = mean_poll_time % (1000'000'000); 
  
  return mean_poll_time - mean_poll_time_ns + ts0
    + (ts0 - mean_poll_time_ns < -500'000'000) * 1000'000'000
    - (ts0 - mean_poll_time_ns >  500'000'000) * 1000'000'000;
}

//Bottom functions

void crt::DecoderICARUSCRT::FragmentToTriggers(const artdaq::Fragment& artFrag, std::unique_ptr<std::vector<CRT::Trigger>>& triggers)
  {
    CRT::Fragment frag(artFrag);                                                                                                                                                   
    TLOG(TLVL_DEBUG,"CRT") << "Is this Fragment good?  " << ((frag.good_event())?"true":"false") << "\n";                                                                                                                                                   
    std::vector<CRT::Hit> hits;
    //Make Channel Mapping array
    size_t map_array[64] = {};
    int count = 1;
    bool alternativeMapping = false;
    if (alternativeMapping){
     //Basic alternative map filling.
     for (int i = 64; i>= 1; i--){ 
 	     map_array[i] = count;
	     count++; 
     }
    }
    else{
    //Apropriate map filling.
     int even = 64;
     int odd = 63;
     for (int i = 1; i<= 64; i++){
       if (count <= 4) {
       map_array [i] = even;
       even-=2;
       count++;
       }
       else if (count <= 8){
        map_array[i] = odd;
        odd-=2;
        if (count != 8){
        count++;
        }
        else{
        count = 1;
        }
       }
     }
    }
  
    //Make a CRT::Hit from each non-zero ADC value in this Fragment
    for(size_t hitNum = 0; hitNum < frag.num_hits(); ++hitNum)
    {
      const auto hit = *(frag.hit(hitNum)); 
      size_t offline_channel = map_array[hit.channel]; 
      hits.emplace_back(offline_channel, hit.adc);      
    }
                                                                                                                                                   
    try
    {  
      triggers->emplace_back(fChannelMapB.at(frag.module_num()), frag.fifty_mhz_time(), std::move(hits)); 
    }
    catch(const std::out_of_range& e)
    {
       TLOG(TLVL_WARNING,"CRT")<< "Got CRT channel number " << frag.module_num() << " that is greater than the number of boards"
                                        << " in the channel map: " << fChannelMapB.size() << ".  Throwing out this Trigger.\n";
    }    

  }

void crt::DecoderICARUSCRT::SetEarliestTime(const artdaq::Fragment& frag)
  {
    CRT::Fragment crt(frag);
    if(crt.fifty_mhz_time() < fEarliestTime) fEarliestTime = crt.fifty_mhz_time();
  }
  

void crt::DecoderICARUSCRT::produce(art::Event& evt)
{
  
  util::LocalArtHandleTrackerManager dataCacheRemover
    (evt, fDropRawDataAfterUse);
  //Side and Top Fragment Handles:
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles
     = evt.getMany<artdaq::Fragments>(fInputTagPatterns);
  
  {
    mf::LogDebug log { "DecoderICARUSCRT" };
    log << "Found " << fragmentHandles.size()
      << " CRT fragment data product candidates:";
    for (auto const& handle: fragmentHandles)
      log << "\n - '" << handle.provenance()->inputTag().encode() << '"';
  }
  
  std::vector<icarus::crt::BernCRTTranslator> hit_vector;

  for (auto const& handle : fragmentHandles) {
    if (!handle.isValid()) continue;
    
    dataCacheRemover.registerHandle(handle);
    
    if (handle->empty()) continue;

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
  std::vector<icarus::crt::CRTData> allCRTdata;
  
  for (auto & hit : hit_vector){
    CorrectForCableDelay(hit);  //add PPS cable length

    if(IsSideCRT(hit)) {
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

        icarus::crt::CRTData data;
        data.fMac5  = recipe.destMac5;
        data.fTs0   = CalculateTimestamp(hit);
        data.fTs1   = hit.ts1;
        data.fFlags                   = hit.flags;
        data.fThisPollStart           = hit.this_poll_start;
        data.fLastPollStart           = hit.last_poll_start;
        data.fHitsInPoll              = hit.hits_in_poll;
        data.fCoinc                   = hit.coinc;
        data.fLastAcceptedTimestamp   = hit.last_accepted_timestamp;
        data.fLostHits                = hit.lost_hits;

        unsigned destCh = recipe.firstDestChannel;
        for (unsigned srcCh = recipe.firstSourceChannel; srcCh <= recipe.lastSourceChannel; ++srcCh) {

          data.fAdc[destCh] = hit.adc[srcCh];
          destCh += recipe.direction; // increase or decrease the source

        }
        allCRTdata.push_back(data);
      } // for all recipes
      sideCount++;
    }
    else if (IsTopCRT(hit)) { //not side CRT, therefore check if top CRT
      //this code needs review by the TOP CRT group!!!
        icarus::crt::CRTData data;
        data.fMac5  = fChannelMap->gettopSimMacAddress(hit.mac5); 
        data.fTs0   = CalculateTimestamp(hit);
        data.fTs1   = hit.ts1;
        data.fFlags                   = hit.flags;
        data.fThisPollStart           = hit.this_poll_start;
        data.fLastPollStart           = hit.last_poll_start;
        data.fHitsInPoll              = hit.hits_in_poll;
        data.fCoinc                   = hit.coinc;
        data.fLastAcceptedTimestamp   = hit.last_accepted_timestamp;
        data.fLostHits                = hit.lost_hits;

        memcpy(data.fAdc, hit.adc, 32*sizeof(hit.adc[0]));
        
        allCRTdata.push_back(data);
        topCount++;
    }
    totalCount++;
  } // loop over all hits in an event


//Bottom pieces:
auto triggers = std::make_unique<std::vector<CRT::Trigger>>();
//Next two lines map each index to itself, taken for begin_job from bottom decoder
const size_t nModules = 56; //
for(size_t module = 0; module < nModules; ++module) fChannelMapB.push_back(module);
try
    {      
      std::cout<<"In the bottom try"<<'\n';  
      //-----
      //Try to get artdaq::Fragments produced from CRT data.  The following line is the reason for 
      //this try-catch block.  I don't expect anything else to throw a cet::Exception.
      const auto& fragHandle = evt.getValidHandle<std::vector<artdaq::Fragment>>(fFragTag);
      //std::cout<<"In the try for produce function"<<'\n';//AC:Added to see where just running the decoder goes. 
      if(fLookForContainer)
      {
        //If this is the first event, set fEarliestTime
        if(fEarliestTime == std::numeric_limits<decltype(fEarliestTime)>::max())
        {
	  //Loop every fragment in the Event for each fragment create a containerFragment and loop over each block. Calling SetEarliestTime for each one. 
          for(const auto& frag: *fragHandle) 
          {
            artdaq::ContainerFragment container(frag);
            for(size_t pos = 0; pos < container.block_count(); ++pos) SetEarliestTime(*container[pos]);
          }
        }
        
        for(const auto& artFrag: *fragHandle)
        {
          artdaq::ContainerFragment container(artFrag);
          for(size_t pos = 0; pos < container.block_count(); ++pos) FragmentToTriggers(*container[pos], triggers); 
        }
      }
      else
      {
        //If this is the first event, set fEarliestTime
        if(fEarliestTime == std::numeric_limits<decltype(fEarliestTime)>::max())
        {
          for(const auto& frag: *fragHandle) SetEarliestTime(frag);
        }

        //Convert each fragment into a CRT::Trigger.
        for(const auto& artFrag: *fragHandle) FragmentToTriggers(artFrag, triggers);
      }
    }//End of bottom try to make Triggers vector equivalent to the hit_vector

catch(const cet::exception& exc) //If there are no artdaq::Fragments in this Event, just add an empty container of CRT::Triggers.
    {
       //TLOG(TLVL_WARNING,"CRT") << "No artdaq::Fragments produced by " << fFragTag << " in this event, so "
       //                             << "not doing anything.\n";
        std::cout << "No artdaq::Fragments produced by " << fFragTag << " in this event, so "
                              << "not doing anything.\n";
    }

for(const auto& trigger: *triggers)
       {
        std::cout<<"In the bottom triggers loop"<<'\n';
        //template<typename T, int N>
        //using raw_array = T[N];
        //auto adc_array =  raw_array<decltype(auto),64> {};
        icarus::crt::CRTData data;
          //trigger.dump(outputStream);
          data.fMac5 = trigger.Channel();
          data.fTs0 = trigger.Timestamp();	          
          for(const auto& hit: trigger.Hits()){
            //adc_array[hit.Channel()] = hit.ADC();    
            //data.fAdc = adc_array;
            std::cout<<"ADC: "<<hit.ADC()<<'\n';
            bottomCount++;          
            totalCount++;   
            allCRTdata.push_back(data);
          }           
       }


// move the data which is actually present in the final data product
auto crtdata = std::make_unique<std::vector<icarus::crt::CRTData>>();
  for (icarus::crt::CRTData& crtDataElem: allCRTdata) {
    if (crtDataElem.fMac5 == 0) continue; // not a valid Mac5, data is not present
    crtdata->push_back(std::move(crtDataElem));
  }
evt.put(std::move(crtdata));


//Testing counts
  std::cout<<"Hit counts for each subdetector. Total count: "<<totalCount<<'\n';
  std::cout<<"Top count: "<<topCount<<'\n';
  std::cout<<"Side count: "<<sideCount<<'\n';
  std::cout<<"Bottom count: "<<bottomCount<<'\n';
}//End of produce

DEFINE_ART_MODULE(crt::DecoderICARUSCRT)

