// icaruscode includes
#include "icaruscode/TPC/Compression/PhysCrateCompressedFragment.cc"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// C++ includes
#include <string>
#include <vector>

//-------------------------------------------------------------------------------------------------------------------------

namespace icarus
{
  namespace tpc
  {
    class ProduceCompressedFragment : public art::EDProducer
    {
      public:
        explicit ProduceCompressedFragment(fhicl::ParameterSet const& pset);
        ~ProduceCompressedFragment() override = default;
        void reconfigure(fhicl::ParameterSet const & p);
        void produce(art::Event & e) override;
      private:
        art::InputTag fFragmentsLabel;
    };

    //.............................................................................
    ProduceCompressedFragment::ProduceCompressedFragment(fhicl::ParameterSet const& pset)
    : art::EDProducer::EDProducer(pset)
    {
      this->reconfigure(pset);
      produces<std::vector<artdaq::Fragment>>();
    }

    //.............................................................................
    void ProduceCompressedFragment::reconfigure(fhicl::ParameterSet const& pset)
    {
      fFragmentsLabel = pset.get<art::InputTag>("FragmentsLabel", "daq:PHYSCRATEDATA");
    }

    //.............................................................................
    void ProduceCompressedFragment::produce(art::Event & e)
    {
      // make a unique_ptr to store what we are producing
      std::unique_ptr<std::vector<artdaq::Fragment>> producedFragments(new std::vector<artdaq::Fragment>);

      // get the uncompressed fragments from the event
      art::Handle<std::vector<artdaq::Fragment>> uncompressedFragmentsHandle;
      std::vector<art::Ptr<artdaq::Fragment>>    uncompressedFragments;

      if (e.getByLabel(fFragmentsLabel, uncompressedFragmentsHandle))
        art::fill_ptr_vector(uncompressedFragments, uncompressedFragmentsHandle);

      // loop over the uncompressed fragments and make the compressed fragments
      for (auto const uncompFrag : uncompressedFragments)
      {
        // art::Ptrs are kind of annoying
        // so to get the Fragment we need to `.get()` the Ptr and dereference
        artdaq::Fragment newFrag = icarus::PhysCrateCompressedFragment::compressArtdaqFragment(*(uncompFrag.get()));
        producedFragments->emplace_back(newFrag);
      }

      e.put(std::move(producedFragments));
    }
    
    DEFINE_ART_MODULE(ProduceCompressedFragment)
  }
}

