//std includes
#include <vector>

//ROOT includes
#include "TFile.h"
#include "TH1.h"

//Framework includes
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

//icarus includes
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

//sbndaq includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.cc"

namespace tcpCompression {
  class ValidateCompression : public art::ResultsProducer {
  public:
    explicit ValidateCompression(fhicl::ParameterSet const& pset);
    ~ValidateCompression() override = default;
    void event(art::Event const& evt) override; 
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    art::InputTag fFragmentsLabel;
    bool          fCheckOldFragments;
    bool          fDumpADCs;
    bool          fFoundFirstUncompressed = false;
    bool          fFoundFirstCompressed   = false;
    TH1F*         fCompHist;
    TH1F*         fDcmpHist;
    const icarusDB::IICARUSChannelMap* fChannelMap;
  };// end ValidateCompression class

  //------------------------------------------------------------------
  ValidateCompression::ValidateCompression(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------
  void ValidateCompression::reconfigure(fhicl::ParameterSet const& pset)
  {
    fFragmentsLabel    = pset.get<art::InputTag>("FragmentsLabel"   , "daq:PHYSCRATEDATA");
    fCheckOldFragments = pset.get<bool>         ("CheckOldFragments", false);
    fDumpADCs          = pset.get<bool>         ("DumpADCs"         , false);

    art::ServiceHandle<art::TFileService> tfs;
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    fCompHist = tfs->make<TH1F>(("CompHist_"+fFragmentsLabel.instance()).c_str(), ";Size (bytes);Number of Fragments", 100, 0, 10000000);
    fDcmpHist = tfs->make<TH1F>(("DcmpHist_"+fFragmentsLabel.instance()).c_str(), ";Size (bytes);Number of Fragments", 100, 0, 10000000);
  }

  //------------------------------------------------------------------
  void ValidateCompression::event(art::Event const& evt)
  {
    std::vector<artdaq::Fragment> originalFragments = evt.getProduct<std::vector<artdaq::Fragment>>(fFragmentsLabel);

    // if we want to make a histogram for the waveform
    art::ServiceHandle<art::TFileService> tfs;

    // this is useful
    size_t fragNumber = 0;
    for (auto const& frag : originalFragments)
    {

      // put the fragment in a compression compliant overlay
      icarus::PhysCrateFragment fragOverlay(frag);
      std::string fragCrateName = fChannelMap->getCrateName(frag.fragmentID());

      if (true)
      {
        // iterate through the fragment which wasn't compressed and check which boards are or aren't working
        MF_LOG_VERBATIM("ValidateCompression")
          << "****************************" << '\n'
          << "Crate " << fragCrateName;
        size_t cumulativePrevBlockSize = 0;
        icarus::A2795DataBlock::data_t const* boardData = new icarus::A2795DataBlock::data_t[frag.dataSize()];
        for (size_t board = 0; board < fragOverlay.nBoards(); ++board)
        {
          bool boardCompressed = true;
          size_t word = 0;
          boardData = reinterpret_cast<icarus::A2795DataBlock::data_t const*>
                      ( frag.dataBeginBytes()
                      + (1+board)*sizeof(icarus::PhysCrateDataTileHeader)
                      + (1+board)*sizeof(icarus::A2795DataBlock::Header)
                      + 4*board*sizeof(uint16_t)
                      + cumulativePrevBlockSize                     );

          for (size_t sample = 0; sample < fragOverlay.nSamplesPerChannel(); ++sample)
          {
            bool oddCompressions = false;
            std::vector<icarus::A2795DataBlock::data_t> prevSample(fragOverlay.nChannelsPerBoard(), 0);
            for (size_t cSet = 0; cSet < fragOverlay.nChannelsPerBoard() / 4; ++cSet)
            {
              // get word
              icarus::A2795DataBlock::data_t dataWord = boardData[word];
              
              // check if the word is interpretable as a compressed set of 4 ADC differences
              bool wordNotTagged = ((dataWord & 0xF000) != 0x8000);
              if (not wordNotTagged)
              {
                // if the word is tagged, then it's not compressed, but should exist in a compressed fragment
                // check the next three words are also tagged, if they are then we are good (might need to check the delta sizes)
                // otherwise something has gone wrong
                bool wordTagged_plus1 = (boardData[word + 1] & 0xF000) == 0x8000;
                bool wordTagged_plus2 = (boardData[word + 2] & 0xF000) == 0x8000;
                bool wordTagged_plus3 = (boardData[word + 3] & 0xF000) == 0x8000;
                if (wordTagged_plus1 && wordTagged_plus2 && wordTagged_plus3)
                {
                  if (not boardCompressed)
                  {
                    MF_LOG_VERBATIM("ValidateCompression")
                      << "SOMEHOW FOUND COMPRESSED SEQUENCE IN UNCOMPRESSED BOARD";
                  }
                  // it's good
                  word += 4;
                  cumulativePrevBlockSize += 8;
                  for (size_t cInSet = 0; cInSet < 4; ++cInSet)
                  {
                    icarus::A2795DataBlock::data_t twelveBitDiff = (boardData[word + 1] & 0x0FFF);
                    bool isNeg = (sample != 0) && (twelveBitDiff >> 11);
                    prevSample[4*cSet + cInSet] += isNeg*0xF000 + twelveBitDiff;
                  }
                } else {
                  // something has gone wrong
                  MF_LOG_VERBATIM("ValidateCompression")
                    << "TAGGED BUT NOT FROM COMPRESSED FRAGMENT";
                  boardCompressed = false;
                  word += 4;
                  cumulativePrevBlockSize += 8;
                  for (size_t cInSet = 0; cInSet < 4; ++cInSet)
                  {
                    prevSample[4*cSet + cInSet] = (boardData[word + cInSet] & (~(1<<(fragOverlay.metadata()->num_adc_bits()+1))));
                  }
                }
              } else {
                if (sample == 0)
                {
                  boardCompressed = false;
                }
                if (boardCompressed)
                {
                  oddCompressions = (not oddCompressions);
                  // if we don't have the tag check if it makes more sense as one value or 4 differences
                  // first try adding to previous samples and see if the results are reasonable
                  for (size_t cInSet = 0; cInSet < 4; ++cInSet)
                  {
                    icarus::A2795DataBlock::data_t fourBitDiff = (dataWord >> (4*cInSet)) & 0x000F;
                    bool isNeg = (fourBitDiff >> 3);
                    prevSample[4*cSet + cInSet] += (isNeg*0xFFF0 + fourBitDiff);
                  }
                  ++word;
                  cumulativePrevBlockSize += 2;
                } else {
                  word += 4;
                  cumulativePrevBlockSize += 8;
                  for (size_t cInSet = 0; cInSet < 4; ++cInSet)
                  {
                    prevSample[4*cSet + cInSet] = (boardData[word + cInSet] & (~(1<<(fragOverlay.metadata()->num_adc_bits()+1))));
                  }
                }
              }
            }
            if (oddCompressions)
            {
              ++word;
              cumulativePrevBlockSize += 2;
            }
          }
          std::string boardCompStr = (boardCompressed) ? " is compressed" : " is UNCOMPRESSED";
          MF_LOG_VERBATIM("ValidateCompression")
            << "  Board " << board << boardCompStr;
        }
        uint32_t totalBytes = fragOverlay.nBoards()*sizeof(icarus::PhysCrateDataTileHeader)
                            + fragOverlay.nBoards()*sizeof(icarus::A2795DataBlock::Header)
                            + 4*fragOverlay.nBoards()*sizeof(uint16_t)
                            + cumulativePrevBlockSize;
        MF_LOG_DEBUG("ValidateCompression")
          << "****************************" << '\n'
          << "Crate " << fragCrateName << '\n'
          << "  Estimated data bytes " << totalBytes << '\n'
          << "  should be            " << frag.dataSizeBytes() << '\n'
          << "****************************";
      }

      ++fragNumber;
    }
  }

  //------------------------------------------------------------------
  void ValidateCompression::writeResults(art::Results& r)
  {
    MF_LOG_VERBATIM("ValidateCompression")
      << "Not writing results (yet)";
  }

  //------------------------------------------------------------------
  void ValidateCompression::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(ValidateCompression)
}// end tcpCompression namespace
