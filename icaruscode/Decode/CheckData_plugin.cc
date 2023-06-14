// Framework includes
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"
#include "icarusalg/Utilities/BinaryDumpUtils.h"
#include "icaruscode/Decode/DecoderTools/IDecoderFilter.h"
#include "icaruscode/Decode/DataProducts/CompressedFragment.h"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_arena.h"
#include "tbb/spin_mutex.h"
#include "tbb/concurrent_vector.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/RawDigit.h"

#include "icarus_signal_processing/ICARUSSigProcDefs.h"
#include "icarus_signal_processing/WaveformTools.h"

// ROOT
#include "TH1.h"
#include "TH2.h"

// std includes
#include <iostream>
#include <memory>
#include <string>
#include <bitset>
#include <chrono>

//-------------------------------------------------------------------------------------------------------

namespace daq
{
  class CheckData : public art::ResultsProducer
  //class CheckData : public art::EDProducer
  {
    public:

      // Copnstructor & destructor
      explicit CheckData(fhicl::ParameterSet const & pset);
      ~CheckData() override = default;

      // overrides
      void beginJob() override;
      void event(art::Event const& e) override;
      //void produce(art::Event& e) override;
      void writeResults(art::Results& r) override;
      void clear() override;
      void endJob() override;

      // the interesting functions
      void reconfigure(fhicl::ParameterSet const& pset);
      void checkLabel(art::Handle<std::vector<artdaq::Fragment>>& label);

      // methods for reading compressed fragment
      size_t  countADCDiffs(artdaq::RawDataType const& word); // count how many ADC Diffs are in the word
      int16_t getADCDiffFromWord(size_t const& n, artdaq::RawDataType const& word); // get nth diff in word
      std::pair<size_t, artdaq::RawDataType*>
              advanceNDiffs(size_t const& n, artdaq::RawDataType* loc); // move forward n diffs
                                                                       // return word and loc of ACD diff
                                                                       // in word
      int16_t getADCVal(size_t const& board, size_t const& channel, size_t const& sample); // get ADC value

    private:
  };

  // inputs
  std::vector<art::InputTag> fFragmentsLabelVec; // vector for multiple fragments

  // outputs
  std::map<art::InputTag,
           std::unique_ptr<artdaq::Fragments>> fFragMap; // the DAQ fragments for each tag
  std::map<art::InputTag, TH1D*>               fHistMap; // A histogram of ADC diffs for each tag
  std::map<art::InputTag, TH1D*>               fAbsMap;  // A histogram of ADC abs diffs for each tag
  std::map<art::InputTag, TH1D*>               fCompMap; // Histograms of fragment compressions
  TH1D*                                        fOldTime; // Old Fragment time vs sample
  TH1D*                                        fNewTime; // New Fragment time vs sample
  TH1D*                                        fAltTime; // New (without keys) Fragment time vs sample
  TH1D*                                        fRecTime; // New (with recursion) Fragment time vs sample
  TH1D*                                        fNewDTime;// New - Old time
  TH1D*                                        fAltDTime;// Alt - Old time
  TH1D*                                        fRecDTime;// Rec - Old time

  //......................................................................
  size_t CheckData::countADCDiffs(artdaq::RawDataType const& word)
  {
    // size of word in bits
    size_t wordBits = 8*sizeof(word);
    
    // count the ADC diffs stored in the word
    // check bits until all are used
    size_t nDiffs = 0;
    size_t bitsChecked = 0;

    while (bitsChecked < wordBits)
    {
      if (((word >> (12 + bitsChecked)) & 0xF) == 0x8)
      {
        nDiffs += 1;
      } else {
        nDiffs += 4;
      }
      bitsChecked += 16;
    }
    
    return nDiffs;
  }

  //......................................................................
  CheckData::CheckData(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
    //produces<artdaq::Fragments>();
  }

  //......................................................................
  void CheckData::reconfigure(fhicl::ParameterSet const & pset)
  {
    fFragmentsLabelVec = pset.get<std::vector<art::InputTag>>("FragmentsLabelVec", std::vector<art::InputTag>()={"daq:PHYSCRATEDATA"});

    art::ServiceHandle<art::TFileService>  tfs;
    for (const auto& fragmentLabel : fFragmentsLabelVec)
    {
      std::string histName = fragmentLabel.instance();
      fHistMap[fragmentLabel] = tfs->make<TH1D>((histName+"_adc").c_str(), 
                                                (histName+";ADC Differences;Counts").c_str(),
                                                1001, -500, 500);
      fAbsMap [fragmentLabel] = tfs->make<TH1D>((histName+"_absDiffs").c_str(),
                                                (histName+";|#DeltaADC|;Counts").c_str(),
                                                10001, 0, 10000);
      fCompMap[fragmentLabel] = tfs->make<TH1D>((histName+"_compression").c_str(), 
                                                (histName+";Compression Ratio;Fragments").c_str(),
                                                10000, 0, 1);
      fFragMap[fragmentLabel] = std::make_unique<artdaq::Fragments>();
    }

    // come back to these once you know the actual times
    fOldTime = tfs->make<TH1D>("TimeVsSample_old",
                               "Old Fragment;Time (ns);",
                               1000, 0, 100000);
    fNewTime = tfs->make<TH1D>("TimeVsSample_new",
                               "New Fragment;Time (ns);",
                               1000, 0, 100000);
    fAltTime = tfs->make<TH1D>("TimeVsSample_alt",
                               "New Fragment without keys;Time (ns);",
                               1000, 0, 100000);
    fRecTime = tfs->make<TH1D>("TimeVsSample_rec",
                               "New Fragment with recursion;Time (ns);",
                               1000, 0, 100000);

    fNewDTime = tfs->make<TH1D>("DTimeVsSample_new",
                                "New Fragment;#DeltaTime (ns);",
                                1000, 0, 100000);
    fAltDTime = tfs->make<TH1D>("DTimeVsSample_alt",
                                "New Fragment without keys;#DeltaTime (ns);",
                                1000, 0, 100000);
    fRecDTime = tfs->make<TH1D>("DTimeVsSample_rec",
                                "New Fragment with recursion;#DeltaTime (ns);",
                                1000, 0, 100000);
  }

  //......................................................................
  // Method to clear out the collections of data products after the writeResults method is called.
  void CheckData::clear()
  {
  }

  //......................................................................
  void CheckData::checkLabel(art::Handle<std::vector<artdaq::Fragment>>& handle)
  {
    if (handle.isValid())
    {
      MF_LOG_VERBATIM("CheckData")
        << "Label Is Valid.";
    } else {
      MF_LOG_VERBATIM("CheckData")
        << "Label Is Not Valid.";
    }
  }

  //......................................................................
  void CheckData::beginJob()
  {
  }

  //......................................................................
  void CheckData::event(art::Event const& e)
  //void CheckData::produce(art::Event& e)
  {
    for(const auto& fragmentLabel : fFragmentsLabelVec)
    {
      art::Handle<std::vector<artdaq::Fragment>> daq_handle;
      e.getByLabel(fragmentLabel, daq_handle);

      for (auto const& fragment : *daq_handle)
      {
        icarus::PhysCrateFragment physCrateFragment(fragment); // convert to a physCrateFragment

        size_t nBoardsPerFragment = physCrateFragment.nBoards();
        size_t nChannelsPerBoard  = physCrateFragment.nChannelsPerBoard();
        size_t nSamplesPerChannel = physCrateFragment.nSamplesPerChannel();

        size_t headerSize              = fragment.headerSizeBytes();
        size_t uncompressedPayloadSize = fragment.dataSizeBytes();
        size_t metadataSize            = fragment.sizeBytes() - headerSize - uncompressedPayloadSize;
        size_t compressedPayloadSize   = 0;

        // setup waveform vector which is nBoards x nChannels x nSamples
        // note we'll index by b*nSamples*nChannels + s*nChannels + c
        std::vector<uint16_t> waveform(nBoardsPerFragment * nChannelsPerBoard * nSamplesPerChannel);

        // setup the ADC differences vector
        std::vector<int16_t> adcDiffs(nBoardsPerFragment * nChannelsPerBoard * nSamplesPerChannel);

        // a vector of bools storing if the adc diff is to be compressed
        // we compress in blocks of 4 channels, so only 1/4th the channel size needed
        // note we'll index by (b*nSamples*nChannels + s*nChannels + c)/4
        std::vector<bool> isCompressed(nBoardsPerFragment * (nChannelsPerBoard / 4) * nSamplesPerChannel);

        for (size_t board = 0; board < nBoardsPerFragment; board++)
        {
          //size_t boardIdx = nChannelsPerBoard * (nBoardsPerFragment * fragment_id + board);

          const icarus::A2795DataBlock::data_t* dataBlock = physCrateFragment.BoardData(board);
            
          // a counter for how many blocks of 4 channels are compressed
          size_t nBlocks = nChannelsPerBoard * nSamplesPerChannel / 4;
          size_t nCompressedBlocks = 0;
            
          for (size_t sample = 0; sample < nSamplesPerChannel; sample++)
          {

            for(size_t channel = 0; channel < nChannelsPerBoard; channel++)
            {
              //raw::ChannelID_t           channel_num = boardIdx + channel;

              size_t index = board * nSamplesPerChannel * nChannelsPerBoard
                           + sample * nChannelsPerBoard
                           + channel;
              size_t blockIndex = index / 4;

              waveform[index] = dataBlock[sample * nChannelsPerBoard + channel];

              if (sample == 0)
              {
                // compression takes sample 0 as the referent
                adcDiffs[index] = waveform[index];

                // we don't compress sample 0
                // only have to set this once every 4 channels
                if ((channel % 4) == 3)
                  isCompressed[blockIndex] = false;
              } else {
                size_t prevSampleIndex = board * nSamplesPerChannel * nChannelsPerBoard
                                       + (sample - 1) * nChannelsPerBoard
                                       + channel;
                adcDiffs[index] = waveform[index] - waveform[prevSampleIndex];
                fHistMap[fragmentLabel]->Fill(adcDiffs[index]);
                fAbsMap [fragmentLabel]->Fill(std::abs(adcDiffs[index]));

                // compression is done on blocks of 4 channels
                // if all differences are less than 7 in magnitude each adc difference is compressed
                if ((channel % 4) == 3)
                {
                  bool isCmp =   (std::abs(adcDiffs[index - 3]) < 8)
                              && (std::abs(adcDiffs[index - 2]) < 8)
                              && (std::abs(adcDiffs[index - 1]) < 8)
                              && (std::abs(adcDiffs[index    ]) < 8);
                  
                  isCompressed[blockIndex] = isCmp;

                  nCompressedBlocks += isCmp;
                }
              }
            }
          }
          // each board has a header
          compressedPayloadSize += sizeof(icarus::PhysCrateDataTileHeader)
                                 + sizeof(icarus::A2795DataBlock::Header);

          // calculate size, in bytes, needed for compressed data
          // each uncompressed adcDiff takes 16 bits = 2 bytes
          // each compressed adcDiff is 4 bits = 0.5 bytes
          // each block is 4 adcDiffs, so compression takes us from 8 to 2 bytes per block
          // if there are an odd number of uncompressed blocks there is a 2 byte overflow
          compressedPayloadSize += 8 *  (nBlocks - nCompressedBlocks) 
                                 + 2 *   nCompressedBlocks
                                 + 2 * ((nBlocks - nCompressedBlocks) % 2);

          // each board has a trailer(?)
          compressedPayloadSize += 4*sizeof(uint16_t);
        }

        size_t totalCompressedSize = compressedPayloadSize + headerSize + metadataSize;
        MF_LOG_DEBUG("CheckData")
          << "Compressed Fragment Size: " << totalCompressedSize;

        artdaq::Fragment newFrag(fragment);
        newFrag.resize(compressedPayloadSize);
        auto oldBegin = reinterpret_cast<const uint16_t* const>(fragment.dataBegin()  );
        auto newBegin = reinterpret_cast<      uint16_t* const>(newFrag .dataAddress());

        // loop add data
        size_t nWord_ofOld = 0;
        size_t nWord_ofNew = 0;
        for (size_t board = 0; board < nBoardsPerFragment; board++)
        {
          size_t nHeaderWords = (sizeof(icarus::PhysCrateDataTileHeader) + sizeof(icarus::A2795DataBlock::Header)) / sizeof(uint16_t);
          for (size_t headerWord = 0; headerWord < nHeaderWords; headerWord++)
          {
            *(newBegin + nWord_ofNew) = *(oldBegin + nWord_ofOld);
            nWord_ofOld++;
            nWord_ofNew++;
          }
          for (size_t sample = 0; sample < nSamplesPerChannel; sample++)
          {
            size_t channel = 0;
            size_t nUncompressed = 16;
            while (channel < nChannelsPerBoard)
            {
              size_t index = board  * nSamplesPerChannel * nChannelsPerBoard
                           + sample * nChannelsPerBoard
                           + channel;
              size_t blockIndex = index/4;

              if (not isCompressed[blockIndex])
              {
                *(newBegin + nWord_ofNew) = (adcDiffs[index] & 0x0FFF) + 0x8000;
                nWord_ofOld++;
                nWord_ofNew++;
                channel++;
              } else {
                *(newBegin + nWord_ofNew) =  (adcDiffs[index    ] & 0x000F)
                                          + ((adcDiffs[index + 1] & 0x000F) <<  4)
                                          + ((adcDiffs[index + 2] & 0x000F) <<  8)
                                          + ((adcDiffs[index + 3] & 0x000F) << 12);
                nWord_ofOld += 4;
                nWord_ofNew++;
                channel += 4;
                nUncompressed--;
              }
            }
            // add offset when odd number of uncompressed blocks
            if ((nUncompressed % 2) == 1)
            {
              *(newBegin + nWord_ofNew) = 0;
              nWord_ofNew++;
            }
          }
        }
        
        daq::compressedFragment compPhysFrag(newFrag);
        double fracCompress = (double) newFrag.sizeBytes() / fragment.sizeBytes();
        fCompMap[fragmentLabel]->Fill(fracCompress);
        
        // check compressed fragment
        // we want to time the average adc pull and plot things
        //for (size_t board = 0; board < nBoardsPerFragment; board++)
        //for (size_t board = 0; board < 1; board++)
        //{
        //  for (size_t sample = 0; sample < nSamplesPerChannel; sample++)
        //  {
        //    //for (size_t channel = 0; channel < nChannelsPerBoard; channel++)
        //    for (size_t channel = 32; channel < 33; channel++)
        //    {
        //      
        //      // time how long each call takes
        //      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        //      uint16_t oldVal = physCrateFragment.adc_val          (board, channel, sample);
        //      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

        //      std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
        //      uint16_t newVal = compPhysFrag     .adc_val          (board, channel, sample);
        //      std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();

        //      std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
        //      uint16_t altVal = compPhysFrag     .adc_val_unkeyed  (board, channel, sample);
        //      std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();

        //      std::chrono::steady_clock::time_point t7 = std::chrono::steady_clock::now();
        //      uint16_t recVal = compPhysFrag     .adc_val_recursive(board, channel, sample);
        //      std::chrono::steady_clock::time_point t8 = std::chrono::steady_clock::now();

        //      if ((newVal != oldVal) || (altVal != oldVal) || (recVal != oldVal))
        //        MF_LOG_VERBATIM("CheckData")
        //          << "*******************" << '\n'
        //          << "ADC value mismatch!" << '\n'
        //          << "Old: " << oldVal     << '\n'
        //          << "New: " << newVal     << '\n'
        //          << "Alt: " << altVal     << '\n'
        //          << "Rec: " << recVal     << '\n'
        //          << "*******************";

        //      // time for old in ns
        //      double oldTime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        //      fOldTime->Fill(oldTime);
        //      // time for new in ns
        //      double newTime = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();
        //      fNewTime->Fill(newTime);
        //      // time for alt in ns
        //      double altTime = std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5).count();
        //      fAltTime->Fill(altTime);
        //      // time for rec in ns
        //      double recTime = std::chrono::duration_cast<std::chrono::nanoseconds>(t8 - t7).count();
        //      fRecTime->Fill(recTime);

        //      fNewDTime->Fill(newTime - oldTime);
        //      fAltDTime->Fill(altTime - oldTime);
        //      fRecDTime->Fill(recTime - oldTime);
        //    }
        //  }
        //}

        //fFragMap[fragmentLabel]->push_back(newFrag);
      }
      //e.put(std::move(fFragMap[fragmentLabel]));
    }
  }
  //......................................................................
  //void CheckData::produce(art::Event& e)
  //{
  //  for (const auto& fragmentLabel : fFragmentsLabelVec)
  //  {
  //    e.put(std::move(fFragMap[fragmentLabel]));
  //  } 
  //}

  //......................................................................
  void CheckData::endJob()
  {
  }

  //----------------------------------------------------------------------------
  void CheckData::writeResults(art::Results& r)
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(CheckData)

}
