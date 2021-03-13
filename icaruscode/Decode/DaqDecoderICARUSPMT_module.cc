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
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
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
    explicit DaqDecoderICARUSPMT(fhicl::ParameterSet const & p);
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

    art::InputTag fInputTag;
};

DEFINE_ART_MODULE(DaqDecoderICARUSPMT)

DaqDecoderICARUSPMT::DaqDecoderICARUSPMT(fhicl::ParameterSet const & params): art::EDProducer{params},
    fInputTag(params.get<std::string>("FragmentsLabel", "daq:CAEN1730"))
{
    // Recover the PMT decoder tool
    fDecoderTool = art::make_tool<IDecoder>(params.get<fhicl::ParameterSet>("DecoderTool"));

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
        auto const& daq_handle = event.getValidHandle<artdaq::Fragments>(fInputTag);
    
        // Make sure data available
        if (daq_handle.isValid() && daq_handle->size() > 0)
        {
            for (auto const &rawFrag: *daq_handle)  fDecoderTool->process_fragment(rawFrag);
        }
    }
    catch(...)
    {
        std::cout << "DaqDecoderICARUSPMT: Did not find daq data products to decode" << std::endl;
    }

    fDecoderTool->outputDataProducts(event);

    return;
}


} // end namespace

