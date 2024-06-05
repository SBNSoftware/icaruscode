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
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

//namespace
namespace reprocessRaw
{

  class ProduceCompressed : public art::EDProducer
  {
    public:
      explicit ProduceCompressed(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produceForLabel(art::Event& evt, art::InputTag fFragmentLabel);
      void produce(art::Event& evt) override;

    private:
      std::vector<art::InputTag> fFragmentLabelVec; // which Fragments are we pulling in?
      bool fDebug = false;
      const icarusDB::IICARUSChannelMap* fChannelMap;

  }; // end ProduceCompressed class

  //------------------------------------------------
  ProduceCompressed::ProduceCompressed(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void ProduceCompressed::reconfigure(fhicl::ParameterSet const& pset)
  {
    // What are we procesing
    fFragmentLabelVec = pset.get<std::vector<art::InputTag>>("FragmentLabelVec");
    fDebug = pset.get<bool>("RunDebugMode", false);
    if (fDebug)
      fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();
    
    for (auto const& fragLabel : fFragmentLabelVec)
    {
      MF_LOG_VERBATIM("ProduceCompressed")
        << "Getting fragments " << fragLabel.label() << ":" << fragLabel.instance();

      // we make these things
      produces<std::vector<artdaq::Fragment>>(fragLabel.instance());
    }
  }

  //------------------------------------------------
  void ProduceCompressed::produceForLabel(art::Event& evt, art::InputTag fFragmentLabel)
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
      if (fDebug)
        MF_LOG_VERBATIM("ProduceCompressed")
          << "-------------------------" << '\n'
          << "  size of old_fragment (in bytes) = " << old_fragment.sizeBytes() << '\n'
          << "    size of payload: " << old_fragment.dataSizeBytes() << '\n'
          << "    64-bit words in payload: " << old_fragment.dataSize() << '\n'
          << "    compression scheme: " << old_fragment.metadata<icarus::PhysCrateFragmentMetadata>()->compression_scheme() << '\n'
          << "  size of new_fragment (in byets) = " << new_fragment.sizeBytes() << '\n'
          << "    size of payload: " << new_fragment.dataSizeBytes() << '\n'
          << "    64-bit words in payload: " << new_fragment.dataSize() << '\n'
          << "    compression scheme: " << new_fragment.metadata<icarus::PhysCrateFragmentMetadata>()->compression_scheme() << '\n'
          << "-------------------------";
      new_fragments->emplace_back(new_fragment);

      if (fDebug)
      {
        // check adc vals
        icarus::PhysCrateFragment old_overlay(old_fragment);
        icarus::PhysCrateFragment new_overlay(new_fragment);
        artdaq::detail::RawFragmentHeader::fragment_id_t fragmentID = new_fragment.fragmentID();
        const std::string& crateName = fChannelMap->getCrateName(fragmentID);
        MF_LOG_VERBATIM("ProduceCompressed")
          << "########################################################" << '\n'
          << "Checking compression for FragmentID " << std::hex << fragmentID << std::dec << ", Crate " << crateName;
        bool adcsGood = true;
        for (size_t b = 0; b < old_overlay.nBoards(); ++b)
        {
          for (size_t s = 0; s < old_overlay.nSamplesPerChannel(); ++s)
          {
            for (size_t c = 0; c < old_overlay.nChannelsPerBoard(); ++c)
            {
              uint16_t old_adc_val = old_overlay.adc_val(b, c, s);
              uint16_t new_adc_val = new_overlay.adc_val(b, c, s);
              int16_t adc_difference = (old_adc_val - new_adc_val);
              if (old_adc_val != new_adc_val)
              {
                adcsGood = false;
                MF_LOG_VERBATIM("ProduceCompressed")
                  << "ADC mismatch in board " << b << ", channel "<< c << ", sample " << s << '\n'
                  << "  Old value: " << old_adc_val << '\n'
                  << "  New value: " << new_adc_val << '\n'
                  << "    Difference: " << adc_difference;
              }
            }
          }
        }
        if (adcsGood)
        {
          MF_LOG_VERBATIM("ProduceCompressed")
            << "Fragments for FragmentID " << std::hex << fragmentID << std::dec << ", Crate " << crateName << " agree!" << '\n'
            << "########################################################";
        } else {
          MF_LOG_VERBATIM("ProduceCompressed")
            << "Fragments for FragmentID " << std::hex << fragmentID << std::dec << ", Crate " << crateName << " DO NOT MATCH" << '\n'
            << "  Checking fragment word-by-word:";
          uint64_t dataOffset = 0;
          for (size_t b = 0; b < old_overlay.nBoards(); ++b)
          {
            uint32_t boardDataTileSize = sizeof(icarus::PhysCrateDataTileHeader) + sizeof(icarus::A2795DataBlock::Header);
            for (size_t boardHeaderWord = 0; boardHeaderWord < (boardDataTileSize / sizeof(icarus::A2795DataBlock::data_t)); ++boardHeaderWord)
            {
              MF_LOG_VERBATIM("ProduceCompressed")
                << "    " << std::bitset<16>(icarus::PhysCrateFragment::getA2795Word(new_fragment, dataOffset)) << " | board " << b << " header word (" << dataOffset << ")";
              ++dataOffset;
            }
            std::vector<uint16_t> prevADCVals(64, 0);
            for (size_t s = 0; s < old_overlay.nSamplesPerChannel(); ++s)
            {
              bool oddCompressions = false;
              size_t c = 0;
              while (c < 64)
              {
                uint16_t word = icarus::PhysCrateFragment::getA2795Word(new_fragment, dataOffset);
                if ((word & 0xF000) == 0x8000)
                {
                  int16_t delta = (word & 0x0FFF) + (0xF000 * ((word & 0x0800) == 0x0800));
                  MF_LOG_VERBATIM("ProduceCompressed")
                    << "    " << std::bitset<16>(word) << " | uncompressed word (" << dataOffset << ")" << '\n'
                    << "      board " << b << ", channel " << c << ", sample " << s << " has ADC value " << prevADCVals[c] << " + " << delta << " = " << (prevADCVals[c] + delta) << '\n'
                    << "        ADC value from the old fragment is " << old_overlay.adc_val(b, c, s);
                  prevADCVals[c] += delta;
                  ++c;
                } else {
                  oddCompressions = (not oddCompressions);
                  int16_t delta0 = ((word >> (4*0)) & 0x000F) + (0xFFF0 * ((word & 0x0008) == 0x0008));
                  int16_t delta1 = ((word >> (4*1)) & 0x000F) + (0xFFF0 * ((word & 0x0080) == 0x0080));
                  int16_t delta2 = ((word >> (4*2)) & 0x000F) + (0xFFF0 * ((word & 0x0800) == 0x0800));
                  int16_t delta3 = ((word >> (4*3)) & 0x000F) + (0xFFF0 * ((word & 0x8000) == 0x8000));
                  MF_LOG_VERBATIM("ProduceCompressed")
                    << "    " << std::bitset<16>(word) << " | compressed word (" << dataOffset << ")" << '\n'
                    << "      board " << b << ", channel " << c + 0 << ", sample " << s << " has ADC value " << prevADCVals[c + 0] << " + " << delta0 << " = " << (prevADCVals[c + 0] + delta0) << '\n'
                    << "        ADC value from the old fragment is " << old_overlay.adc_val(b, c + 0, s) << '\n'
                    << "      board " << b << ", channel " << c + 1 << ", sample " << s << " has ADC value " << prevADCVals[c + 1] << " + " << delta1 << " = " << (prevADCVals[c + 1] + delta1) << '\n'
                    << "        ADC value from the old fragment is " << old_overlay.adc_val(b, c + 1, s) << '\n'
                    << "      board " << b << ", channel " << c + 2 << ", sample " << s << " has ADC value " << prevADCVals[c + 2] << " + " << delta2 << " = " << (prevADCVals[c + 2] + delta2) << '\n'
                    << "        ADC value from the old fragment is " << old_overlay.adc_val(b, c + 2, s) << '\n'
                    << "      board " << b << ", channel " << c + 3 << ", sample " << s << " has ADC value " << prevADCVals[c + 3] << " + " << delta3 << " = " << (prevADCVals[c + 3] + delta3) << '\n'
                    << "        ADC value from the old fragment is " << old_overlay.adc_val(b, c + 3, s);
                  prevADCVals[c + 0] += delta0;
                  prevADCVals[c + 1] += delta1;
                  prevADCVals[c + 2] += delta2;
                  prevADCVals[c + 3] += delta3;
                  c += 4;
                }
                ++dataOffset;
              }
              if (oddCompressions)
              {
                MF_LOG_VERBATIM("ProduceCompressed")
                  << "    " << std::bitset<16>(icarus::PhysCrateFragment::getA2795Word(new_fragment, dataOffset)) << " | board " << b << " offset word for odd compressions in sample " << s << "(" << dataOffset << ")";
                ++dataOffset;
              }
            }
            for (size_t trailerWord = 0; trailerWord < 4; ++trailerWord)
            {
              MF_LOG_VERBATIM("ProduceCompressed")
                << "    " << std::bitset<16>(icarus::PhysCrateFragment::getA2795Word(new_fragment, dataOffset)) << " | board " << b << " trailer word (" << dataOffset << ")";
              ++dataOffset;
            }
          }
          MF_LOG_VERBATIM("ProduceCompressed")
            << "########################################################";
        }
      }
    }

    // put the new fragments into the event
    evt.put(std::move(new_fragments), fFragmentLabel.instance());
  }

  //------------------------------------------------
  void ProduceCompressed::produce(art::Event& evt)
  {
    for (auto const& fragLabel : fFragmentLabelVec)
      this->produceForLabel(evt, fragLabel);
  }

  DEFINE_ART_MODULE(ProduceCompressed)
} // end namespace
