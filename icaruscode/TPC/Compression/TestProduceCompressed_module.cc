// std inlcudes
#include <string>
#include <vector>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// wiremod testing
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

//namespace
namespace reprocessRaw
{

  class TestProduceCompressed : public art::EDProducer
  {
    public:
      explicit TestProduceCompressed(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

    private:
      art::InputTag fFragmentLabel; // which Fragments are we pulling in?

  }; // end TestProduceCompressed class

  //------------------------------------------------
  TestProduceCompressed::TestProduceCompressed(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void TestProduceCompressed::reconfigure(fhicl::ParameterSet const& pset)
  {
    // What are we procesing
    fFragmentLabel  = pset.get<art::InputTag>("FragmentLabel", "");

    MF_LOG_VERBATIM("TestProduceCompressed")
      << "Getting fragments " << fFragmentLabel.label() << ":" << fFragmentLabel.instance();

    // we make these things
    produces<std::vector<artdaq::Fragment>>();
  }

  //------------------------------------------------
  void TestProduceCompressed::produce(art::Event& evt)
  {
    // get the old fragments
    art::Handle<std::vector<artdaq::Fragment>> fragHandle;
    evt.getByLabel(fFragmentLabel, fragHandle);
    auto const& old_fragments(*fragHandle);

    // make a vector to put the new fragments in
    std::unique_ptr<std::vector<artdaq::Fragment>> new_fragments(new std::vector<artdaq::Fragment>());

    // fill the new fragments vector
    for (auto const& fragment : old_fragments)
      new_fragments->emplace_back(icarus::PhysCrateFragment::compressArtdaqFragment(fragment));

    // put the new fragments into the event
    evt.put(std::move(new_fragments));
  }
  DEFINE_ART_MODULE(TestProduceCompressed)
} // end namespace
