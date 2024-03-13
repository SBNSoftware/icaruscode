////////////////////////////////////////////////////////////////////////
// Class:       CRTRawDecoder
// File:        CRTRawDecoder_module.cc
//
// Edited from Andrew Olivier original by Axel Campos
// Stripped down version of decoder for ICARUS Bottom CRT
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
//dune-artdaq includes
#include "artdaq-core/Data/ContainerFragment.hh"
//Log system
//#include "TRACE/tracemf.h"
#include "Bottom/Fragments/CRTFragment.hh"
#include "Bottom/Data/CRTTrigger.h"

//#include "/exp/icarus/app/users/bcaxel/testdir/srcs/icaruscode/icaruscode/CRT/CRTDecoder/Bottom/Fragments/CRTFragment.hh"

//ROOT includes
#include "TGraph.h"

//c++ includes
#include <memory>

using namespace CRT;

//A CRTRawDecoder takes artdaq::Fragments made from Cosmic Ray Tagger input 
//and produces a CRT::Trigger for each time a CRT module triggered in this Event.  
namespace CRT
{
  class CRTRawDecoder : public art::EDProducer 
  {
    public:
      explicit CRTRawDecoder(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.
    
      // Plugins should not be copied or assigned.
      CRTRawDecoder(CRTRawDecoder const &) = delete;
      CRTRawDecoder(CRTRawDecoder &&) = delete;
      CRTRawDecoder & operator = (CRTRawDecoder const &) = delete;
      CRTRawDecoder & operator = (CRTRawDecoder &&) = delete;
    
      // Required functions.
      void produce(art::Event & e) override;
    
      // Selected optional functions.  I might get services that I expect to change here.  So far, seems like only service I might need is timing.
      void beginJob() override;
      //void beginRun(art::Run & r) override;
      //void beginSubRun(art::SubRun & sr) override;
    
    private:
    
      // Configuration data
      const art::InputTag fFragTag; //Label and product instance name (try eventdump.fcl on your input file) of the module that produced 
                                      //artdaq::Fragments from CRT data that I will turn into CRT::Triggers. Usually seems to be "daq" from artdaq.  

      const bool fLookForContainer; //Should this module look for "container" artdaq::Fragments instead of "single" artdaq::Fragments?  
                                    //If you ran the CRT board reader with the "request_mode" FHICL parameter set to "window", you need 
                                    //to look for "container" Fragments and probably won't get any useful Fragments otherwise.  If you 
                                    //ran the CRT board reader with "request_mode" set to "ignore", then set fLookForContainer to "false" 
                                    //which is the default.  
                                    //You probably want this set to "true" because the CRT should be run in "window" request_mode for 
                                    //production.  You might set it to "false" if you took debugging data in which you wanted to ignore 
                                    //timestamp matching with other detectors.  
      const bool fMatchOfflineMapping; //Should the hardware channel mapping match the order in which AuxDetGeos are sorted in the offline 
                                       //framework?  Default is true.  Please set to false for online monitoring so that CRT experts can 
                                       //understand problems more quickly.  
      std::vector<size_t> fChannelMap; //Simple map from raw data module number to offline module number.  Initialization depends on 
                                       //fMatchOfflineMapping above.

      // Compartmentalize internal functionality so that I can reuse it with both regular Fragments and "container" Fragments
      void FragmentToTriggers(const artdaq::Fragment& artFrag, std::unique_ptr<std::vector<CRT::Trigger>>& triggers);

      // For the first Event of every job, I want to set fEarliestTime to the earliest time in that Event
      void SetEarliestTime(const artdaq::Fragment& frag);

      //Sync diagnostic plots
      
      

      uint64_t fEarliestTime; //Earliest time in clock ticks      

  };
  
  
  CRTRawDecoder::CRTRawDecoder(fhicl::ParameterSet const & p): EDProducer{p}, fFragTag(p.get<std::string>("RawDataTag")), 
                                                               fLookForContainer(p.get<bool>("LookForContainer", false)),
                                                               fMatchOfflineMapping(p.get<bool>("MatchOfflineMapping", true)),
                                                               fEarliestTime(std::numeric_limits<decltype(fEarliestTime)>::max())
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<CRT::Trigger>>();
    consumes<std::vector<artdaq::Fragment>>(fFragTag);
    //Initialize Metric Manager
    //if (p.has_key("metrics")) {
    //sbndaq::InitializeMetricManager(p.get<fhicl::ParameterSet>("metrics"));
    // }
    //if (p.has_key("metric_config")) {
    //  sbndaq::GenerateMetricConfig(p.get<fhicl::ParameterSet>("metric_config"));
    //}
    //sbndaq::InitializeMetricManager(pset.get<fhicl::ParameterSet>("metrics")); //This causes the error for no "metrics" at the beginning or the end
    //sbndaq::GenerateMetricConfig(p.get<fhicl::ParameterSet>("metric_channel_config"));
    //sbndaq::GenerateMetricConfig(p.get<fhicl::ParameterSet>("metric_board_config"));  //This line cauess the code to not be able to compile
    //art::ServiceHandle<art::TFileService> tfs;
    //tfs->registerFileSwitchCallback(this, &CRTRawDecoder::createSyncPlots);
  }

  void CRTRawDecoder::FragmentToTriggers(const artdaq::Fragment& artFrag, std::unique_ptr<std::vector<CRT::Trigger>>& triggers)
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
      //Metrics per channel            
      //std::cout<<"Metrics Sent: " << "Channel: "<<std::to_string(hit.channel)<<"ADC: "<< std::to_string(hit.adc) <<'\n'; 
      //sbndaq::sendMetric("CRT_channel_bottom", std::to_string(hit.channel), "ADC", std::to_string(hit.adc), 0, artdaq::MetricMode::Average);
      //MF_LOG_DEBUG("CRTRaw") << "Channel: " << (int)(hit.channel) << "\n"
      //                    << "ADC: " << hit.adc << "\n";
      //Determine the offline channel number for each strip      
      size_t offline_channel = map_array[hit.channel]; 
      hits.emplace_back(offline_channel, hit.adc);      
    }
                                                                                                                                                   
    try
    {  
      triggers->emplace_back(fChannelMap.at(frag.module_num()), frag.fifty_mhz_time(), std::move(hits)); 
    }
    catch(const std::out_of_range& e)
    {
       TLOG(TLVL_WARNING,"CRT")<< "Got CRT channel number " << frag.module_num() << " that is greater than the number of boards"
                                        << " in the channel map: " << fChannelMap.size() << ".  Throwing out this Trigger.\n";
    }    

  }

  void CRTRawDecoder::SetEarliestTime(const artdaq::Fragment& frag)
  {
    CRT::Fragment crt(frag);
    if(crt.fifty_mhz_time() < fEarliestTime) fEarliestTime = crt.fifty_mhz_time();
  }
  
  //Read artdaq::Fragments produced by fFragTag, and use CRT::Fragment to convert them to CRT::Triggers.  
  void CRTRawDecoder::produce(art::Event & e)
  {
    //Create an empty container of CRT::Triggers.  Any Triggers in this container will be put into the 
    //event at the end of produce.  I will try to fill this container, but just not produce any CRT::Triggers 
    //if there are no input artdaq::Fragments.  
    auto triggers = std::make_unique<std::vector<CRT::Trigger>>();

    try
    {
      //Try to get artdaq::Fragments produced from CRT data.  The following line is the reason for 
      //this try-catch block.  I don't expect anything else to throw a cet::Exception.
      const auto& fragHandle = e.getValidHandle<std::vector<artdaq::Fragment>>(fFragTag);
      
      if(fLookForContainer)
      {
        //If this is the first event, set fEarliestTime
        if(fEarliestTime == std::numeric_limits<decltype(fEarliestTime)>::max())
        {
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
    }
    catch(const cet::exception& exc) //If there are no artdaq::Fragments in this Event, just add an empty container of CRT::Triggers.
    {
       TLOG(TLVL_WARNING,"CRT") << "No artdaq::Fragments produced by " << fFragTag << " in this event, so "
                                    << "not doing anything.\n";
    }

    //Put a vector of CRT::Triggers into this Event for other modules to read.
    e.put(std::move(triggers));
  }
  
  void CRT::CRTRawDecoder::beginJob()
  {
    const size_t nModules = 56; //14; 
    for(size_t module = 0; module < nModules; ++module) fChannelMap.push_back(module); //Map each index to itself    
  }

}

DEFINE_ART_MODULE(CRT::CRTRawDecoder)

