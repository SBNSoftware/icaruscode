//std includes
#include <vector>
#include <chrono>
#include <limits>
#include <tuple>

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

//icaruscode includes
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

//sbndaq-artdaq-core includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.cc"

namespace tcpCompression
{
  class FragmentPlots : public art::EDAnalyzer
  {
  public:
    explicit FragmentPlots(fhicl::ParameterSet const& pset);
    ~FragmentPlots() override = default;
    void analyze(art::Event const& evt) override; 
    void reconfigure(fhicl::ParameterSet const& p);
    std::tuple<uint16_t, double> timeADC(const icarus::PhysCrateFragment& frag,
                                         const size_t& board, const size_t& channel, const size_t& sample,
                                         const bool& key = false);

  private:
    art::InputTag fFragmentsLabel;
    bool fInitiallyCompressed;
    std::unique_ptr<icarusDB::IICARUSChannelMap> fChannelMap;
    std::unique_ptr<TH1D> adcDiffs;
    std::unique_ptr<TH1D> fragSizeRatio;
    std::unique_ptr<TH1D> fragSize_compressed;
    std::unique_ptr<TH1D> fragSize_uncompressed;
    std::array<std::unique_ptr<TH3D>, 9> randReads_keyed_arr;
    std::array<std::unique_ptr<TH3D>, 9> randReads_compressed_arr;
    std::array<std::unique_ptr<TH3D>, 9> randReads_uncompressed_arr;
  };// end FragmentPlots class

  //------------------------------------------------------------------
  std::tuple<uint16_t, double> FragmentPlots::timeADC(const icarus::PhysCrateFragment& frag,
                                                      const size_t& board, const size_t& channel, const size_t& sample,
                                                      const bool& key = false)
  {
    uint16_t adc = 0;
    double msTiming = std::numeric_limits<double>::max();
    if ((not key) && (frag.isCompressed()))
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      size_t s = 0;
      const uint16_t* boardZeroData = frag.BoardData(0);
      size_t nWordsIntoFrag = 0;
      for (size_t b = 0; b < board; ++b)
      {
        // header
        nWordsIntoFrag += (sizeof(PhysCrateDataTileHeader) + sizeof(A2795DataBlock::Header)) / sizeof(uint16_t);
        // channel data
        for (size_t s = 0; s <= frag.nSamplesPerChannel(); ++s)
        {
          bool oddCompressions = true;
          for (size_t c = 0; c < frag.nChannels(); c += 4)
          {
            uint16_t nextWord = *(boardZeroData + nWordsIntoFrag);
            bool tagged = ((nextWord & 0xF000) == 0x8000);
            nWordsIntoFrag += (tagged) ? 4 : 1;
            oddCompressions = (tagged) ? (not oddCompressions) : oddCompressions;
          }
          // compression offset
          if (oddCompressions)
            ++nWordsIntoFrag;
        }
        // trailer
        nWordsIntoFrag += 4;
      }
      // header for our board
      nWordsIntoFrag += (sizeof(PhysCrateDataTileHeader) + sizeof(A2795DataBlock::Header)) / sizeof(uint16_t);
      for (size_t s = 0; s <= sample; ++s)
      {
        bool oddCompressions = true;
        for (size_t c = 0; c < frag.nChannels(); c += 4)
        {
          uint16_t nextWord = *(boardZeroData + nWordsIntoFrag);
          bool tagged = ((nextWord & 0xF000) == 0x8000);
          nWordsIntoFrag += (tagged) ? 4 : 1;
          oddCompressions = (tagged) ? (not oddCompressions) : oddCompressions;
          if ((c / 4) == (channel / 4))
          {
            if (tagged)
            {
              uint16_t twelveBitDiff = nextWord & 0x0FFF;
              bool isNeg = (fourBitDiff >> 11);
              adc += isNeg*0xF000 + twelveBitDiff;
            }
            else
            {
               uint16_t fourBitDiff = (nextWord >> 4*(channel % 4)) & 0x000F;
               bool isNeg = (fourBitDiff >> 3);
               adc += isNeg*0xFFF0 + fourBitDiff;
            }
          }
          if (s == sample)
            break;
        }
        if (s == sample)
          break;
        // compression offset
        if (oddCompressions)
          ++nWordsIntoFrag;
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::micro> msDur(t1 - t0);
      msTiming = msDur.count();
    }
    else if (frag.isCompressed())
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      const uint16_t* boardData = frag.BoardData(board);
      size_t prevSampleWords = 0;
      adc = ( *(boardData + channel) & (~(1<<(frag.metadata()->num_adc_bits()+1))) );
      for (size_t s = 1; s <= sample; ++s)
      {
        prevSampleWords += (frag.SampleBytesFromKey(frag.CompressionKey(board, s - 1))) / sizeof(uint16_t);
        size_t curSampleWords = 0;
        for (size_t cSet = 0; cSet < (channel / 4) - 1; ++cSet)
        {
          curSampleWords += (frag.CompressionKey(board, s)[cSet]) ? 1 : 4;
        }
        if (frag.CompressionKey(board, s)[channel / 4])
        {
          uint16_t word = *(boardData + prevSampleWords + curSampleWords);
          uint16_t fourBitDiff = (word >> 4*(channel % 4)) & 0x000F;
          bool isNeg = (fourBitDiff >> 3);
          adc += isNeg*0xFFF0 + fourBitDiff;
        }
        else
        {
          uint16_t twelveBitDiff = *(boardData + prevSampleWords + curSampleWords + (channel % 4)) & 0x0FFF;
          bool isNeg = (fourBitDiff >> 11);
          adc += isNeg*0xF000 + twelveBitDiff;
        }
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::micro> msDur(t1 - t0);
      msTiming = msDur.count();
    }
    else
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      adc = frag.adc_val(board, channel, sample);
      auto t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::micro> msDur(t1 - t0);
      msTiming = msDur.count();
    }
    return std::tie(adc, msTiming);
  }

  //------------------------------------------------------------------
  FragmentPlots::FragmentPlots(fhicl::ParameterSet const& pset) : EDAnalyzer(config)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------
  void FragmentPlots::reconfigure(fhicl::ParameterSet const& pset)
  {
    fFragmentsLabel = pset.get<art::InputTag>("FragmentsLabel", "daq:PHYSCRATEDATATPC??");

    adcDiffs.reset(
      tfs->make<TH1D>("adcDiffs", ";Difference in Subsequent ADC Samples;Number of Differences;", 21, -10.5, 10.5));
   
    fragSizeRatio.reset(
      tfs->make<TH1D>("fragSizeRatio", ";Compression Factor;Number of Fragments", 100, 0, 1));

    std::string fragSize_name = "fragSize";
    std::string fragSize_title = ";Fragment Size (GB);Fragments";
    int nBins = 100;
    double binMin = 0;
    double binMax = 10;
    fragSize_compressed  .reset(tfs->make<TH1D>((fragSize_name + "_compressed").c_str(), fragSize_title, nBins, binMin, binMax));
    fragSize_uncompressed.reset(tfs->make<TH1D>((fragSize_name+"_uncompressed").c_str(), fragSize_title, nBins, binMin, binMax));
    
    std::string randReads_name = "randReads";
    std::string randReads_title = ";Channel;Sample;Read Time (#mus)";
    int nBins_randReads = 10000;
    double binMin_randReads = 0;
    double binMax_randReads = 10000;
    for (size_t b = 0; b < 9; ++b)
    {
      randReads_keyed_arr[b]       .reset(tfs->make<TH3D>((randReads_name +      "_keyed_"+std::to_string(b)).c_str(),
                                                          randReads_title,
                                                          64, 0, 63,
                                                          4096, 0, 4095,
                                                          nBins, binMin, binMax));
      randReads_compressed_arr[b]  .reset(tfs->make<TH3D>((randReads_name + "_compressed_"+std::to_string(b)).c_str(),
                                                          randReads_title,
                                                          64, 0, 63,
                                                          4096, 0, 4095,
                                                          nBins, binMin, binMax));
      randReads_uncompressed_arr[b].reset(tfs->make<TH3D>((randReads_name+"_uncompressed_"+std::to_string(b)).c_str(),
                                                          randReads_title,
                                                          64, 0, 63,
                                                          4096, 0, 4095,
                                                          nBins, binMin, binMax));
    }
  }

  //------------------------------------------------------------------
  void FragmentPlots::analyze(art::Event const& evt)
  {
    // get the fragments
    art::Handle<std::vector<artdaq::Fragment>> fragmentsHandle;
    evt.getByLabel(fFragmentLabel, fragmentsHandle);

    // loop over the fragments
    for (size_t fragIdx = 0; fragIdx < fragmentsHandle->size(); ++fragIdx)
    {
      // instantiate the overlays for compress/uncompressed fragments
      icarus::PhysCrateFragment   compressedFragment(fragmentsHandle->at(fragIdx));
      icarus::PhysCrateFragment uncompressedFragment = compressedFragment.makeUncompressedFragment();
     
      // get size (in GB) of the payloads
      // DataPayloadSize returns the size in bytes as a size_t, so do some converting
      double   compressedGB =  0.000000001 * std::static_cast<double>(  compressedFragment.DataPayloadSize());
      double uncompressedGB =  0.000000001 * std::static_cast<double>(uncompressedFragment.DataPayloadSize());
      fragSize_compressed  ->Fill(  compressedGB);
      fragSize_uncompressed->Fill(uncompressedGB);

      // get adcs
      // store the times for each board/sample/channel
      // keep track of the differences as well
      for (size_t board = 0; board < compressedFragment.nBoards(); ++boards)
      {
        for (size_t channel = 0; channel < compressedFragment.nChannelsPerBoard(); ++channel)
        {
          uint16_t prevSampleADC = 0;
          for (size_t sample = 0; sample < compressedFragment.nSamplesPerChannel(); ++sample)
          {
            auto [sampleADC_keyed,        timing_keyed       ] = this->timeADC(  compressedFragment,
                                                                               board, channel, sample, true);
            auto [sampleADC_compressed,   timing_compressed  ] = this->timeADC(  compressedFragment,
                                                                               board, channel, sample);
            auto [sampleADC_uncompressed, timing_uncompressed] = this->timeADC(uncompressedFragment,
                                                                               board, channel, sample);
            // assert that the ADCs are all identical
            if ((sampleADC_keyed != sampleADC_compressed) || (sampleADC_keyed != sampleADC_uncompressed))
            {
              std::string whatStr = fFragmentsLabel.encode()
                                  + " event " + std::to_string(evt.id().event())
                                  + ", Board "   + std::to_string(board)
                                  + ", channel " + std::to_string(channel)
                                  + ", sample "  + std::to_string(sample)
                                  + " has inconsistent ADC counts:\n"
                                  + "         Keyed: " + std::to_string(sampleADC_keyed) + '\n'
                                  + "    Compressed: " + std::to_string(sampleADC_compressed) + '\n'
                                  + "  Uncompressed: " + std::to_string(sampleADC_uncompressed);
              throw std::runtime_error(whatStr);
            }
            randReads_keyed_arr       [board]->Fill(channel, sample, timing_keyed       );
            randReads_compressed_arr  [board]->Fill(channel, sample, timing_compressed  );
            randReads_uncompressed_arr[board]->Fill(channel, sample, timing_uncompressed);
            if (sample != 0)
            {
              adcDiff->Fill(sampleADC_keyed - prevSampleADC);
            }
            prevSampleADC = sampleADC_keyed;
          }
        }
      }
    }
  }

  DEFINE_ART_MODULE(FragmentPlots)
}// end tcpCompression namespace
