///////////////////////////////////////////////////////////////////////
// Class:       DaqDecoderIcarusPMT
// Plugin Type: producer (art v2_09_06)
// File:        DaqDecoderIcarusPMT.cxx
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Principal/Event.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "artdaq-core/Data/Fragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"

//some standard C++ includes
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

namespace daq
{

/*
  * The Decoder module takes as input "NevisTPCFragments" and
  * outputs raw::RawDigits. It also handles in and all issues
  * with the passed in header and fragments (or at least it will).
*/
class DaqDecoderICARUSPMT: public art::EDProducer 
{
public:
    
    /// FHiCL configuration of the module.
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> FragmentsLabel {
        Name("FragmentsLabel"),
        Comment("data product with the PMT fragments from DAQ"),
        "daq:CAEN1730" // default
        };
      
      fhicl::DelegatedParameter DecoderTool {
        Name("DecoderTool"),
        Comment
          ("configuration of decoding tool (daq::IDecoder) to be loaded and used")
        };
      
    }; // struct Config
    
    using Parameters = art::EDProducer::Table<Config>;
  
  
    explicit DaqDecoderICARUSPMT(Parameters const & params);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    // Plugins should not be copied or assigned.
    DaqDecoderICARUSPMT(DaqDecoderICARUSPMT const &) = delete;
    DaqDecoderICARUSPMT(DaqDecoderICARUSPMT &&) = delete;

    DaqDecoderICARUSPMT & operator = (DaqDecoderICARUSPMT const &) = delete;

    DaqDecoderICARUSPMT & operator = (DaqDecoderICARUSPMT &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

private:

    // Tools for decoding fragments depending on type
    std::unique_ptr<IDecoder>  fDecoderTool;  ///< Decoder tools
    
    
    art::InputTag const fInputTag;

};

DEFINE_ART_MODULE(DaqDecoderICARUSPMT)

DaqDecoderICARUSPMT::DaqDecoderICARUSPMT(Parameters const & params)
  : art::EDProducer{params}
  , fInputTag{ params().FragmentsLabel() }
{
    fDecoderTool = art::make_tool<IDecoder>
      (params().DecoderTool.get<fhicl::ParameterSet>());
    // Recover the PMT decoder tool

    consumes<artdaq::Fragments>(fInputTag);
    
    // Announce what the tool will do
    fDecoderTool->produces(producesCollector());

    return;
}

void DaqDecoderICARUSPMT::produce(art::Event & event)
{
    // storage for waveform
    fDecoderTool->initializeDataProducts();

    // Protect for runs with no PMT info
    try
    {
        // Recover the data fragments for the PMT 
        auto const& fragments = event.getByLabel<artdaq::Fragments>(fInputTag);
    
        // Make sure data available
        if (!fragments.empty())
        {
            for (auto const &rawFrag: fragments)  fDecoderTool->process_fragment(rawFrag);
        }
    }
    catch(art::Exception const& e) {
        std::cout << "DaqDecoderICARUSPMT: Did not find daq data products to decode:"
          << "\n" << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout << "DaqDecoderICARUSPMT: Did not find daq data products to decode" << std::endl;
    }

    fDecoderTool->outputDataProducts(event);

    return;
}


} // end namespace

