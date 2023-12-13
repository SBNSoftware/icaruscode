// std inlcudes
#include <string>
#include <vector>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// overlay testing
#include "icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h"
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
    fFragmentLabel = pset.get<art::InputTag>("FragmentLabel", "");

    MF_LOG_VERBATIM("TestProduceCompressed")
      << "Getting fragments " << fFragmentLabel.label() << ":" << fFragmentLabel.instance();

    // we make these things
    // try passing it the input label
    produces<std::vector<artdaq::Fragment>>(fFragmentLabel.instance());
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
    for (auto const& old_fragment : old_fragments)
    {
      artdaq::Fragment new_fragment = icarus::PhysCrateFragment::fragmentSwitch(old_fragment, true);
      //artdaq::Fragment new_fragment = icarus::PhysCrateFragment::compressArtdaqFragment(old_fragment);
      MF_LOG_DEBUG("TestProduceCompressed")
        << "-------------------------" << '\n'
        << "  size of old_fragment (in bytes) = " << old_fragment.sizeBytes() << '\n'
        << "    sizd of payload: " << old_fragment.dataSizeBytes() << '\n'
        << "  size of new_fragment (in byets) = " << new_fragment.sizeBytes() << '\n'
        << "    sizd of payload: " << new_fragment.dataSizeBytes();
      new_fragments->emplace_back(new_fragment);

      // check adc vals
      icarus::PhysCrateFragment old_overlay(old_fragment);
      icarus::PhysCrateFragment new_overlay(new_fragment);
      for (size_t b = 0; b < old_overlay.nBoards(); ++b)
      {
        for (size_t s = 0; s < old_overlay.nSamplesPerChannel(); ++s)
        {
          bool sampleWrong = false;
          std::array<bool, 16> sampleKey = new_overlay.CompressionKey(b, s);
          size_t nCompressed = 0;
          for (auto const& bit : sampleKey)
            nCompressed += bit;
          //MF_LOG_VERBATIM("TestProduceCompressed")
          //  << "Sample key for board " << b << ", sample " << s << ": " << sampleKey[0]
          //                                                              << sampleKey[1]
          //                                                              << sampleKey[2]
          //                                                              << sampleKey[3]
          //                                                              << sampleKey[4]
          //                                                              << sampleKey[5]
          //                                                              << sampleKey[6]
          //                                                              << sampleKey[7]
          //                                                              << sampleKey[8]
          //                                                              << sampleKey[9]
          //                                                              << sampleKey[10]
          //                                                              << sampleKey[11]
          //                                                              << sampleKey[12]
          //                                                              << sampleKey[13]
          //                                                              << sampleKey[14]
          //                                                              << sampleKey[15] << '\n'
          //  << " there are " << nCompressed << " compressions in the sample.";
          for (size_t c = 0; c < old_overlay.nChannelsPerBoard(); ++c)
          {
            if (old_overlay.adc_val(b, c, s) != new_overlay.adc_val(b, c, s) && not sampleWrong)
            {
              MF_LOG_VERBATIM("TestProduceCompressed")
                << "  MISMATCH ON BOAD " << b << ", CHANNEL " << c << ", SAMPLE " << s << '\n'
                << "    OLD ADC: " << old_overlay.adc_val(b, c, s) << '\n'
                << "    NEW ADC: " << new_overlay.adc_val(b, c, s);
              sampleWrong = true;
            }
          }
          if (sampleWrong)
            break;
        }
      }

      // dump the Fragments
      MF_LOG_VERBATIM("TestProduceCompressed")
        << "---------------------------------" << '\n'
        << "DUMP FRAGMENT — OLD"               << '\n'
        << sbndaq::dumpFragment(old_fragment)  << '\n'
        << "#################################" << '\n'
        << "DUMP FRAGMENT — NEW"               << '\n'
        << sbndaq::dumpFragment(new_fragment)  << '\n'
        << "---------------------------------";
    }


    // put the new fragments into the event
    evt.put(std::move(new_fragments), fFragmentLabel.instance());
  }
  DEFINE_ART_MODULE(TestProduceCompressed)
} // end namespace
