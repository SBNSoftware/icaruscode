//std includes
#include <vector>
#include <chrono>
#include <limits>
#include <tuple>

//ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TTree.h"

//Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
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

//sbndaq-artdaq-core includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

namespace tcpCompression
{
  class FragmentPlots : public art::EDAnalyzer
  {
  public:
    explicit FragmentPlots(fhicl::ParameterSet const& pset);
    ~FragmentPlots() override = default;
    void analyze(art::Event const& evt) override; 
    void reconfigure(fhicl::ParameterSet const& p);
    void endJob() override;
    void fillADCCalls(const icarus::PhysCrateFragment& frag,
                      const size_t& board, const size_t& channel,
                      const bool& key = false);
    std::tuple<uint16_t, double> timeADC(const icarus::PhysCrateFragment& frag,
                                         const size_t& board, const size_t& channel, const size_t& sample,
                                         const bool& key = false);

  private:
    bool fDebug;
    art::InputTag fFragmentLabel;
    std::unique_ptr<TH1D> adcDiffs;
    std::unique_ptr<TH1D> fragSizeRatio;
    std::unique_ptr<TH1D> fragSize_compressed;
    std::unique_ptr<TH1D> fragSize_uncompressed;
    size_t fBoard;
    size_t fChannel;
    size_t fSample;
    std::unique_ptr<TTree> randReads_keyed;
    double readTime_keyed;
    std::unique_ptr<TTree> randReads_compressed;
    double readTime_compressed;
    std::unique_ptr<TTree> randReads_uncompressed;
    double readTime_uncompressed;
  };// end FragmentPlots class

  //------------------------------------------------------------------
  void FragmentPlots::fillADCCalls(const icarus::PhysCrateFragment& frag,
                                   const size_t& board, const size_t& channel,
                                   const bool& key)
  {
    // only fill ADC diffs on keyed calls
    // initialize some value stores
    uint16_t prevADC = 0;
    uint16_t currADC = 0;
    double msTiming = 0;
    fBoard = board;
    fChannel = channel;

    // read differently base on the fragment/key use
    //
    if ((not key) && (frag.isCompressed()))
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      const uint16_t* boardZeroData = frag.BoardData(0); // starts after the header for board 0
      size_t nWordsUntilBoardStart = 0;
      for (size_t b = 0; b < board; ++b)
      {
        // header
        if (b != 0)
        {
          nWordsUntilBoardStart += (sizeof(icarus::PhysCrateDataTileHeader) + sizeof(icarus::A2795DataBlock::Header))
                                   / sizeof(uint16_t);
        }
        // channel data
        for (size_t sample = 0; sample < frag.nSamplesPerChannel(); ++sample)
        {
          bool oddCompressions = false;
          for (size_t cSet = 0; cSet < (frag.nChannelsPerBoard() / 4); ++cSet)
          {
            uint16_t nextWord = *(boardZeroData + nWordsUntilBoardStart);
            bool tagged = ((nextWord & 0xF000) == 0x8000);
            nWordsUntilBoardStart += (tagged) ? 4 : 1;
            oddCompressions = (tagged) ? (not oddCompressions) : oddCompressions;
          }
          // compression offset
          if (oddCompressions)
            ++nWordsUntilBoardStart;
        }
        // trailer
        nWordsUntilBoardStart += 4;
      }
      // header for our board
      if (board != 0)
      {
        nWordsUntilBoardStart += (sizeof(icarus::PhysCrateDataTileHeader) + sizeof(icarus::A2795DataBlock::Header))
                                 / sizeof(uint16_t);
      }
      size_t nWordsIntoBoard = 0;
      for (size_t sample = 0; sample < frag.nSamplesPerChannel(); ++sample)
      {
        bool oddCompressions = false;
        for (size_t cSet = 0; cSet < (frag.nChannelsPerBoard() / 4); ++cSet)
        {
          uint16_t nextWord = *(boardZeroData + nWordsUntilBoardStart + nWordsIntoBoard);
          bool tagged = ((nextWord & 0xF000) == 0x8000);
          if ((channel / 4) == cSet)
          {
            if (tagged)
            {
              uint16_t twelveBitDiff = *(boardZeroData + nWordsUntilBoardStart + nWordsIntoBoard + (channel % 4)) & 0x0FFF;
              bool isNeg = (twelveBitDiff >> 11) && (sample != 0);
              currADC += isNeg*0xF000 + twelveBitDiff;
            }
            else
            {
              uint16_t fourBitDiff = (nextWord >> 4*(channel % 4)) & 0x000F;
              bool isNeg = (fourBitDiff >> 3);
              currADC += isNeg*0xFFF0 + fourBitDiff;
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::micro> msDur(t1 - t0);
            msTiming += msDur.count();
            fSample = sample;
            readTime_compressed = msTiming;
            randReads_compressed->Fill();
            prevADC = currADC;
            t0 = std::chrono::high_resolution_clock::now();
          }
          nWordsIntoBoard += (tagged) ? 4 : 1;
          oddCompressions = (tagged) ? (not oddCompressions) : oddCompressions;
        }
        // compression offset
        if (oddCompressions)
          ++nWordsIntoBoard;
      }
    }
    else if (frag.isCompressed())
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      const uint16_t* boardData = frag.BoardData(board);
      size_t prevSampleWords = 0;
      currADC = ( *(boardData + channel) & 0x0FFF );
      auto t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::micro> msDur(t1 - t0);
      msTiming += msDur.count();
      fSample = 0;
      readTime_keyed = msTiming;
      randReads_keyed->Fill();
      t0 = std::chrono::high_resolution_clock::now();
      for (size_t sample = 1; sample < frag.nSamplesPerChannel(); ++sample)
      {
        prevSampleWords += (frag.SampleBytesFromKey(frag.CompressionKey(board, sample - 1))) / sizeof(uint16_t);
        size_t curSampleWords = 0;
        for (size_t cSet = 0; cSet < (channel / 4); ++cSet)
        {
          curSampleWords += (frag.CompressionKey(board, sample).at(cSet)) ? 1 : 4;
        }
        if (frag.CompressionKey(board, sample).at(channel / 4))
        {
          uint16_t word = *(boardData + prevSampleWords + curSampleWords);
          uint16_t fourBitDiff = (word >> 4*(channel % 4)) & 0x000F;
          bool isNeg = (fourBitDiff >> 3);
          currADC += isNeg*0xFFF0 + fourBitDiff;
        }
        else
        {
          uint16_t twelveBitDiff = *(boardData + prevSampleWords + curSampleWords + (channel % 4)) & 0x0FFF;
          bool isNeg = (twelveBitDiff >> 11);
          currADC += isNeg*0xF000 + twelveBitDiff;
        }
        t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> msDur(t1 - t0);
        msTiming += msDur.count();
        fSample = sample;
        readTime_keyed = msTiming;
        randReads_keyed->Fill();
        adcDiffs->Fill(currADC - prevADC);
        prevADC = currADC;
        t0 = std::chrono::high_resolution_clock::now();
      }
    }
    else
    {
      for (size_t sample = 0; sample < frag.nSamplesPerChannel(); ++sample)
      {
        fSample = sample;
        prevADC = currADC;
        auto t0 = std::chrono::high_resolution_clock::now();
        currADC = frag.adc_val(board, channel, sample);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> msDur(t1 - t0);
        readTime_uncompressed = msDur.count();
        randReads_uncompressed->Fill();
      }
    }
  }

  //------------------------------------------------------------------
  std::tuple<uint16_t, double> FragmentPlots::timeADC(const icarus::PhysCrateFragment& frag,
                                                      const size_t& board, const size_t& channel, const size_t& sample,
                                                      const bool& key)
  {
    uint16_t adc = 0;
    double msTiming = std::numeric_limits<double>::max();
    if ((not key) && (frag.isCompressed()))
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      const uint16_t* boardZeroData = frag.BoardData(0); // starts after the header for board 0
      size_t nWordsUntilBoardStart = 0;
      for (size_t b = 0; b < board; ++b)
      {
        // header
        if (b != 0)
        {
          nWordsUntilBoardStart += (sizeof(icarus::PhysCrateDataTileHeader) + sizeof(icarus::A2795DataBlock::Header))
                                   / sizeof(uint16_t);
        }
        // channel data
        for (size_t s = 0; s < frag.nSamplesPerChannel(); ++s)
        {
          bool oddCompressions = false;
          for (size_t cSet = 0; cSet < (frag.nChannelsPerBoard() / 4); ++cSet)
          {
            uint16_t nextWord = *(boardZeroData + nWordsUntilBoardStart);
            bool tagged = ((nextWord & 0xF000) == 0x8000);
            nWordsUntilBoardStart += (tagged) ? 4 : 1;
            oddCompressions = (tagged) ? (not oddCompressions) : oddCompressions;
          }
          // compression offset
          if (oddCompressions)
            ++nWordsUntilBoardStart;
        }
        // trailer
        nWordsUntilBoardStart += 4;
      }
      // header for our board
      if (board != 0)
      {
        nWordsUntilBoardStart += (sizeof(icarus::PhysCrateDataTileHeader) + sizeof(icarus::A2795DataBlock::Header))
                                 / sizeof(uint16_t);
      }
      size_t nWordsIntoBoard = 0;
      for (size_t s = 0; s <= sample; ++s)
      {
        bool oddCompressions = false;
        for (size_t cSet = 0; cSet < (frag.nChannelsPerBoard() / 4); ++cSet)
        {
          uint16_t nextWord = *(boardZeroData + nWordsUntilBoardStart + nWordsIntoBoard);
          bool tagged = ((nextWord & 0xF000) == 0x8000);
          if ((channel / 4) == cSet)
          {
            if (tagged)
            {
              uint16_t twelveBitDiff = *(boardZeroData + nWordsUntilBoardStart + nWordsIntoBoard + (channel % 4)) & 0x0FFF;
              bool isNeg = (twelveBitDiff >> 11) && (s != 0);
              adc += isNeg*0xF000 + twelveBitDiff;
            }
            else
            {
              uint16_t fourBitDiff = (nextWord >> 4*(channel % 4)) & 0x000F;
              bool isNeg = (fourBitDiff >> 3);
              adc += isNeg*0xFFF0 + fourBitDiff;
            }
            if (s == sample)
              break;
          }
          nWordsIntoBoard += (tagged) ? 4 : 1;
          oddCompressions = (tagged) ? (not oddCompressions) : oddCompressions;
        }
        if (s == sample)
          break;
        // compression offset
        if (oddCompressions)
          ++nWordsIntoBoard;
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
      adc = ( *(boardData + channel) & 0x0FFF );
      for (size_t s = 1; s <= sample; ++s)
      {
        prevSampleWords += (frag.SampleBytesFromKey(frag.CompressionKey(board, s - 1))) / sizeof(uint16_t);
        size_t curSampleWords = 0;
        for (size_t cSet = 0; cSet < (channel / 4); ++cSet)
        {
          curSampleWords += (frag.CompressionKey(board, s).at(cSet)) ? 1 : 4;
        }
        if (frag.CompressionKey(board, s).at(channel / 4))
        {
          uint16_t word = *(boardData + prevSampleWords + curSampleWords);
          uint16_t fourBitDiff = (word >> 4*(channel % 4)) & 0x000F;
          bool isNeg = (fourBitDiff >> 3);
          adc += isNeg*0xFFF0 + fourBitDiff;
        }
        else
        {
          uint16_t twelveBitDiff = *(boardData + prevSampleWords + curSampleWords + (channel % 4)) & 0x0FFF;
          bool isNeg = (twelveBitDiff >> 11);
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
  FragmentPlots::FragmentPlots(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------------
  void FragmentPlots::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDebug         = pset.get<bool>         ("Debug",          false);
    fFragmentLabel = pset.get<art::InputTag>("FragmentLabel", "daq:PHYSCRATEDATATPC??");
    mf::LogVerbatim("FragmentPlots")
      << "Running with" << '\n'
      << "  fDebug: " << std::boolalpha << fDebug << '\n'
      << "  fFragmentLabel: " << fFragmentLabel;

    art::ServiceHandle<art::TFileService> tfs;

    adcDiffs.reset(
      tfs->make<TH1D>("adcDiffs", ";Difference in Subsequent ADC Samples;Number of Differences;", 129, -64.5, 64.5));
   
    fragSizeRatio.reset(
      tfs->make<TH1D>("fragSizeRatio", ";Compression Factor;Number of Fragments", 100, 0, 1));

    std::string fragSize_name = "fragSize";
    std::string fragSize_title = ";Fragment Size (GB);Fragments";
    int nBins = 100;
    double binMin = 0;
    double binMax = 10;
    fragSize_compressed  .reset(tfs->make<TH1D>((fragSize_name + "_compressed").c_str(), fragSize_title.c_str(),
                                                nBins, binMin, binMax));
    fragSize_uncompressed.reset(tfs->make<TH1D>((fragSize_name+"_uncompressed").c_str(), fragSize_title.c_str(),
                                                nBins, binMin, binMax));
    
    randReads_keyed.reset(tfs->make<TTree>("randReads_keyed", "Random Reads - Keyed"));
    randReads_keyed->Branch("board",    &fBoard);
    randReads_keyed->Branch("channel",  &fChannel);
    randReads_keyed->Branch("sample",   &fSample);
    randReads_keyed->Branch("readTime", &readTime_keyed);

    randReads_compressed.reset(tfs->make<TTree>("randReads_compressed", "Random Reads - Compressed"));
    randReads_compressed->Branch("board",    &fBoard);
    randReads_compressed->Branch("channel",  &fChannel);
    randReads_compressed->Branch("sample",   &fSample);
    randReads_compressed->Branch("readTime", &readTime_compressed);

    randReads_uncompressed.reset(tfs->make<TTree>("randReads_uncompressed", "Random Reads - Uncompressed"));
    randReads_uncompressed->Branch("board",    &fBoard);
    randReads_uncompressed->Branch("channel",  &fChannel);
    randReads_uncompressed->Branch("sample",   &fSample);
    randReads_uncompressed->Branch("readTime", &readTime_uncompressed);

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
      if (fDebug)
        mf::LogVerbatim("FragmentPlots")
          << "Processing " <<fFragmentLabel << " fragment " << (fragIdx + 1) << " of " << fragmentsHandle->size(); 

      // instantiate the overlays for compress/uncompressed fragments
      icarus::PhysCrateFragment   compressedFragment(fragmentsHandle->at(fragIdx));
      artdaq::Fragment bareUncompFrag = icarus::PhysCrateFragment::decompressArtdaqFragment(fragmentsHandle->at(fragIdx));
      icarus::PhysCrateFragment uncompressedFragment(bareUncompFrag);
     
      // get size (in GB) of the payloads
      // DataPayloadSize returns the size in bytes as a size_t, so do some converting
      double   compressedMB =  0.000001 * static_cast<double>(  compressedFragment.DataPayloadSize());
      double uncompressedMB =  0.000001 * static_cast<double>(uncompressedFragment.DataPayloadSize());
      fragSize_compressed  ->Fill(  compressedMB);
      fragSize_uncompressed->Fill(uncompressedMB);
      fragSizeRatio->Fill(compressedMB / uncompressedMB);

      // get adcs
      // store the times for each board/sample/channel
      // keep track of the differences as well
      for (size_t board = 0; board < compressedFragment.nBoards(); ++board)
      {
        fBoard = board;
        for (size_t channel = 0; channel < compressedFragment.nChannelsPerBoard(); ++channel)
        {
          fChannel = channel;
          this->fillADCCalls(  compressedFragment, board, channel, true);
          this->fillADCCalls(  compressedFragment, board, channel);
          this->fillADCCalls(uncompressedFragment, board, channel);
        }
      }
    }
  }

  //------------------------------------------------------------------
  void FragmentPlots::endJob()
  {
    // release unique_ptrs so the file services can make the hists
    adcDiffs              .release();
    fragSizeRatio         .release();
    fragSize_compressed   .release();
    fragSize_uncompressed .release();
    randReads_keyed       .release();
    randReads_compressed  .release();
    randReads_uncompressed.release();
  }

  DEFINE_ART_MODULE(FragmentPlots)
}// end tcpCompression namespace
