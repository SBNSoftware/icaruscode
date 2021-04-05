////////////////////////////////////////////////
//   
//    File: TriggerDecoder_tool.cc
//       
//    Description: Starting point for extracting ICARUS trigger fragment information into LArSoft object TBD 
//
//    Author: Jacob Zettlemoyer, FNAL
//
///////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/ExternalTrigger.h" //JCZ: TBD, placeholder for now to represent the idea

#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerUDPFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"

#include <cstdlib>
#include <iostream>
#include <memory>

namespace daq 
{
  
  class TriggerDecoder : virtual public IDecoder
  {
  public:
    explicit TriggerDecoder(fhicl::ParameterSet const &pset);
    ~TriggerDecoder();
    
    virtual void produces(art::ProducesCollector&) override;
    virtual void configure(const fhicl::ParameterSet&) override;
    virtual void initializeDataProducts() override;
    virtual void process_fragment(const artdaq::Fragment &fragment) override;
    virtual void outputDataProducts(art::Event &event) override;
   
  private: 
    using TriggerCollection = std::vector<raw::ExternalTrigger>;
    using TriggerPtr = std::unique_ptr<TriggerCollection>;
    TriggerPtr fTrigger;
    bool fDiagnosticOutput; //< Produces large number of diagnostic messages, use with caution!
    bool fDebug; //< Use this for debugging this tool
    //Add in trigger data member information once it is selected, current LArSoft object likely not enough as is
  };

  TriggerDecoder::TriggerDecoder(fhicl::ParameterSet const &pset)
  {
    this->configure(pset);
  }

  TriggerDecoder::~TriggerDecoder() {}
  
  void TriggerDecoder::produces(art::ProducesCollector& collector) 
  {
    collector.produces<TriggerCollection>();
  }

  void TriggerDecoder::configure(fhicl::ParameterSet const &pset) 
  {
    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);
    fDebug = pset.get<bool>("Debug", false);
    return;
  }
  
  void TriggerDecoder::initializeDataProducts()
  {
    //use until different object chosen 
    //fTrigger = new raw::Trigger();
    fTrigger = TriggerPtr(new TriggerCollection);
    return;
  }

  void TriggerDecoder::process_fragment(const artdaq::Fragment &fragment)
  {
    size_t fragmentID = fragment.fragmentID();
    icarus::ICARUSTriggerUDPFragment frag(fragment);
    //artdaq::Fragment::timestamp_t ts = frag.timestamp();
    icarus::ICARUSTriggerUDPFragmentMetadata meta = *frag.Metadata();
    if(fDiagnosticOutput)
      std::cout << meta << std::endl;
    int local_trig = frag.getName();
    int local_event_no = frag.getEventNo();
    int local_seconds = frag.getSeconds();
    long local_nsec = frag.getNanoSeconds();
    int wr_trig = frag.getWRName();
    int wr_event_no = frag.getWREventNo();
    long wr_seconds = frag.getWRSeconds();
    long wr_nsec = frag.getWRNanoSeconds();
    uint64_t trigger_ts = wr_seconds*1e9 + wr_nsec;
    fTrigger->emplace_back(wr_event_no, trigger_ts);
    if(fDiagnosticOutput || fDebug)
    {
      std::cout << "Fragment ID: " << fragmentID <<  " Local Trigger ID: " << local_trig << " Event Number: " << local_event_no << " Local Seconds: " << local_seconds << " Local Nanoseconds: " << local_nsec << std::endl;
      std::cout << "White Rabbit Trigger ID: " << wr_trig << " Event Number: " << wr_event_no << " White Rabbit POSIX Seconds: " << std::fixed << wr_seconds << " White Rabbit Timestamp Nanoseconds: " << wr_nsec << std::endl;
      std::cout << "Full Timestamp = " << std::fixed << trigger_ts << std::endl;
    }
    
    //Once we have full trigger data object, set up and place information into there
    return;
  }

  void TriggerDecoder::outputDataProducts(art::Event &event)
  {
    //Place trigger data object into raw data store 
    event.put(std::move(fTrigger));
    return;
  }

  DEFINE_ART_CLASS_TOOL(TriggerDecoder)

}




  
