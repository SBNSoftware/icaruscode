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

#include "lardataalg/Utilities/quantities/spacetime.h" // util::quantities::nanosecond
#include "lardataobj/RawData/ExternalTrigger.h" //JCZ: TBD, placeholder for now to represent the idea
#include "lardataobj/Simulation/BeamGateInfo.h" //JCZ:, another placeholder I am unsure if this and above will be combined at some point into a dedicated object 

#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerUDPFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"

#include <cstdlib>
#include <iostream>
#include <memory>

namespace daq 
{
  
  class TriggerDecoder : public IDecoder
  {
    using nanoseconds = util::quantities::nanosecond;
  public:
    explicit TriggerDecoder(fhicl::ParameterSet const &pset);
    
    virtual void produces(art::ProducesCollector&) override;
    virtual void configure(const fhicl::ParameterSet&) override;
    virtual void initializeDataProducts() override;
    virtual void process_fragment(const artdaq::Fragment &fragment) override;
    virtual void outputDataProducts(art::Event &event) override;
   
  private: 
    using TriggerCollection = std::vector<raw::ExternalTrigger>;
    using TriggerPtr = std::unique_ptr<TriggerCollection>;
    using BeamGateInfoCollection = std::vector<sim::BeamGateInfo>;
    using BeamGateInfoPtr = std::unique_ptr<BeamGateInfoCollection>;
    TriggerPtr fTrigger;
    TriggerPtr fPrevTrigger;
    BeamGateInfoPtr fBeamGateInfo; 
    bool fDiagnosticOutput; //< Produces large number of diagnostic messages, use with caution!
    bool fDebug; //< Use this for debugging this tool
    int fOffset; //< Use this to determine additional correction needed for TAI->UTC conversion from White Rabbit timestamps. Needs to be changed if White Rabbit firmware is changed and the number of leap seconds changes! 
    //Add in trigger data member information once it is selected, current LArSoft object likely not enough as is
    uint64_t fLastTimeStamp = 0;
    long fLastEvent = 0;
    
    /// Name of the data product instance for the current trigger.
    static std::string const CurrentTriggerInstanceName;
    
    /// Name of the data product instance for the previous trigger.
    static std::string const PreviousTriggerInstanceName;
    
    static constexpr nanoseconds BNBgateDuration { 1600. };
    static constexpr nanoseconds NuMIgateDuration { 9500. };
  };


  std::string const TriggerDecoder::CurrentTriggerInstanceName {};
  std::string const TriggerDecoder::PreviousTriggerInstanceName { "previous" };
  

  TriggerDecoder::TriggerDecoder(fhicl::ParameterSet const &pset)
  {
    this->configure(pset);
  }

  
  void TriggerDecoder::produces(art::ProducesCollector& collector) 
  {
    collector.produces<TriggerCollection>(CurrentTriggerInstanceName);
    collector.produces<TriggerCollection>(PreviousTriggerInstanceName);
    collector.produces<BeamGateInfoCollection>();
  }
    

  void TriggerDecoder::configure(fhicl::ParameterSet const &pset) 
  {
    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);
    fDebug = pset.get<bool>("Debug", false);
    fOffset = pset.get<int>("TAIOffset", 2);
    return;
  }
  
  void TriggerDecoder::initializeDataProducts()
  {
    //use until different object chosen 
    //fTrigger = new raw::Trigger();
    fTrigger = std::make_unique<TriggerCollection>();
    fPrevTrigger = std::make_unique<TriggerCollection>();
    fBeamGateInfo = BeamGateInfoPtr(new BeamGateInfoCollection);
    return;
  }

  void TriggerDecoder::process_fragment(const artdaq::Fragment &fragment)
  {
    uint64_t artdaq_ts = fragment.timestamp();
    icarus::ICARUSTriggerUDPFragment frag(fragment);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerString(buffer);
    uint64_t wr_ts = datastream_info.getNanoseconds_since_UTC_epoch() + fOffset;
    int gate_type = datastream_info.gate_type;
    long delta_gates_bnb = frag.getDeltaGatesBNB();
    long delta_gates_numi = frag.getDeltaGatesOther(); //this is actually NuMI due to abrupt changes in trigger board logic
    long delta_gates_other = frag.getDeltaGatesNuMI();
    uint64_t fLastTrigger = 0;
    if(fDiagnosticOutput || fDebug)
    {
      std::cout << "Full Timestamp = " << wr_ts << std::endl;
      double cross_check = wr_ts - artdaq_ts;
      if(abs(cross_check) > 1000)
        std::cout << "Loaded artdaq TS and fragment data TS are > 1 ms different! They are " << cross_check << " nanoseconds different!" << std::endl;
      std::cout << delta_gates_bnb << " " << delta_gates_numi << " " << delta_gates_other << std::endl; // nonsensical print statement to avoid error that I don't use these...until we have an object to store them in...    
    }
    if(gate_type == 1) //BNB
    {
      
      fTrigger->emplace_back(datastream_info.wr_event_no, wr_ts);
      if(datastream_info.wr_event_no == 1)
      {
        fLastEvent = 0;
      }
      else 
      {
        fLastEvent = datastream_info.wr_event_no - 1;
      }
      fLastTrigger = frag.getLastTimestampBNB();
      fPrevTrigger->emplace_back(fLastEvent, fLastTrigger);
      fBeamGateInfo->emplace_back
        (wr_ts, BNBgateDuration.convertInto<nanoseconds>().value(), sim::kBNB);
    }
    if(gate_type == 2) //NuMI
    {
      fTrigger->emplace_back(datastream_info.wr_event_no, wr_ts);
      if(datastream_info.wr_event_no == 1)
        {
          fLastEvent = 0;
        }
      else
        fLastEvent = datastream_info.wr_event_no - 1;
      fLastTrigger = frag.getLastTimestampOther(); //actually NuMI for now
      fPrevTrigger->emplace_back(fLastEvent, fLastTrigger);
      fBeamGateInfo->emplace_back
        (wr_ts, NuMIgateDuration.convertInto<nanoseconds>().value(), sim::kNuMI);
    }
    //Once we have full trigger data object, set up and place information into there
    return;
  }

  void TriggerDecoder::outputDataProducts(art::Event &event)
  {
    //Place trigger data object into raw data store 
    event.put(std::move(fTrigger), CurrentTriggerInstanceName);
    event.put(std::move(fPrevTrigger), PreviousTriggerInstanceName);
    event.put(std::move(fBeamGateInfo));
    return;
  }

  DEFINE_ART_CLASS_TOOL(TriggerDecoder)

}





  
