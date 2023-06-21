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
#include <chrono>
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
      //art::Handle<std::vector<artdaq::Fragment>> uncompressedFragmentsHandle;
      //std::vector<art::Ptr<artdaq::Fragment>>    uncompressedFragments;
      //if (e.getByLabel(fFragmentsLabel, uncompressedFragmentsHandle))
      //  art::fill_ptr_vector(uncompressedFragments, uncompressedFragmentsHandle);
      std::vector<artdaq::Fragment> uncompressedFragments = e.getProduct<std::vector<artdaq::Fragment>>(fFragmentsLabel);

      // set up averages to keep track of how long the (de)compression takes and how effective it is
      double avgCompTime   = 0;
      double avgDecompTime = 0;
      double avgCompSize   = 0;
      double avgDecompSize = 0;

      // loop over the uncompressed fragments and make the compressed fragments
      for (auto const uncompFrag : uncompressedFragments)
      {
        // time the creation of the compressed fragment
        std::chrono::steady_clock::time_point compT1 = std::chrono::steady_clock::now();
        // art::Ptrs are kind of annoying
        // so to get the Fragment we need to `.get()` the Ptr and dereference
        //artdaq::Fragment newFrag = icarus::PhysCrateCompressedFragment::compressArtdaqFragment(*(uncompFrag.get()));
        artdaq::Fragment newFrag = icarus::PhysCrateCompressedFragment::compressArtdaqFragment(uncompFrag);
        std::chrono::steady_clock::time_point compT2 = std::chrono::steady_clock::now();
        producedFragments->emplace_back(newFrag);
        double compTime = std::chrono::duration_cast<std::chrono::nanoseconds>(compT2 - compT1).count();
        double compSize = newFrag.dataSizeBytes();
        MF_LOG_VERBATIM("ProduceCompressedFragment")
          << "It took "
          << compTime
          << " ns to compress the artdq::Fragment" << '\n'
          << "Fragment is " << compSize << " bytes.";
        avgCompTime += compTime;
        avgCompSize += compSize;

        // for the sake of completeness time the decompression
        std::chrono::steady_clock::time_point decompT1 = std::chrono::steady_clock::now();
        artdaq::Fragment decompFrag = icarus::PhysCrateCompressedFragment::decompressArtdaqFragment(newFrag);
        std::chrono::steady_clock::time_point decompT2 = std::chrono::steady_clock::now();
        double decompTime = std::chrono::duration_cast<std::chrono::nanoseconds>(decompT2 - decompT1).count();
        double decompSize = decompFrag.dataSizeBytes();
        MF_LOG_VERBATIM("ProduceCompressedFragment")
          << "It took "
          << decompTime
          << " ns to decompress the artdq::Fragment" << '\n'
          << "Fragment is " << decompSize << " bytes.";
        avgDecompTime += decompTime;
        avgDecompSize += decompSize;
      }
      
      avgCompTime   /= uncompressedFragments.size();
      avgCompSize   /= uncompressedFragments.size();
      avgDecompTime /= uncompressedFragments.size();
      avgDecompSize /= uncompressedFragments.size();
      MF_LOG_VERBATIM("ProduceCompressedFragment")
        << "Average compression time: "   << avgCompTime << " ns"     << '\n'
        << "Average decompression time: " << avgDecompTime << " ns"   << '\n'
        << "Average compressed size: "    << avgCompSize << " bytes"  << '\n'
        << "Average decompressed size: "  << avgDecompSize << " bytes";

      e.put(std::move(producedFragments));
    }
    
    DEFINE_ART_MODULE(ProduceCompressedFragment)
  }
}

