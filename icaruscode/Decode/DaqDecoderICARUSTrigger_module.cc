///////////////////////////////////////////////////////
///// Class: DaqDecoderICARUSTrigger
///// Plugin Type: producer
///// File: DaqDecoderICARUSTrigger.cc
//////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "artdaq-core/Data/Fragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"

#include <iostream>
#include <cstdlib>
#include <vector>

namespace daq 
{
  class DaqDecoderICARUSTrigger: public art::EDProducer
  {
  public:
    
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> FragmentsLabel {
        Name("FragmentsLabel"),
        Comment("input tag for the DAQ fragment of the trigger"),
        "daq:ICARUSTriggerUDP"
        };
      
      fhicl::DelegatedParameter DecoderTool {
        Name("DecoderTool"),
        Comment("configuration for the trigger decoding tool")
        };
      
      
    }; // Config
    
    using Parameters = art::EDProducer::Table<Config>;
    
    explicit DaqDecoderICARUSTrigger(Parameters const & params);

    void produce(art::Event & e) override;
    
  private:
    std::unique_ptr<IDecoder> fDecoderTool;
    art::InputTag const fInputTag;

  };

  DEFINE_ART_MODULE(DaqDecoderICARUSTrigger)
  
  DaqDecoderICARUSTrigger::DaqDecoderICARUSTrigger(Parameters const & params)
    : art::EDProducer{params}
    , fDecoderTool{
      art::make_tool<IDecoder>(params().DecoderTool.get<fhicl::ParameterSet>())
      }
    , fInputTag{ params().FragmentsLabel() }
  {
    fDecoderTool->produces(producesCollector());
  }
  
  void DaqDecoderICARUSTrigger::produce(art::Event & event)
  {
    fDecoderTool->initializeDataProducts();
    auto const & daq_handle = event.getValidHandle<artdaq::Fragments>(fInputTag);
    if(daq_handle.isValid() && daq_handle->size() > 0)
    {
      for (auto const & rawFrag: *daq_handle) fDecoderTool->process_fragment(rawFrag);
    }
    else
      std::cout << "No Trigger Fragment Information Found!" << std::endl;

    fDecoderTool->outputDataProducts(event);
    return;

  }
      
} //end namespace
