// std inlcudes
#include <arpa/inet.h>
#include <array>
#include <string>
#include <vector>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "artdaq-core/Data/Fragment.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

//namespace
namespace reprocessRaw
{
  struct TileHeader
  {
    uint32_t token;
    uint32_t info1;
    uint32_t info2;
    uint32_t info3;
    uint32_t timeinfo;
    uint32_t pkt_fmt_ver  : 8;
    uint32_t crate_id     : 8;
    uint32_t board_id     : 8;
    uint32_t board_status : 8;
    uint32_t packSize;

    TileHeader(){};
  }; // end of TileHeader struct

  class MetaData
  {
    public:
      MetaData(){}
      MetaData(uint32_t run_number,
               uint32_t n_boards,
               uint32_t channels_per_board,
               uint32_t samples_per_channel,
               uint32_t adcs_per_sample,
               uint32_t compression,
               std::vector<uint32_t> const& idvec)
      {
        _run_number = run_number;
        _num_boards = n_boards;
        _channels_per_board = channels_per_board;
        _samples_per_channel = samples_per_channel;
        _num_adc_bits = adcs_per_sample;
        _compression_scheme = compression;
        SetBoardIDs(idvec);
      }

      uint32_t const& run_number() const { return _run_number; }
      uint32_t const& samples_per_channel() const { return _samples_per_channel; }
      uint32_t const& num_adc_bits() const { return _num_adc_bits; }
      uint32_t const& channels_per_board() const { return _channels_per_board; }
      uint32_t const& num_boards() const { return _num_boards; }
      uint32_t const& compression_scheme() const { return _compression_scheme; }
      uint32_t const& board_id(size_t i) const { return _board_ids[i]; }

      void SetBoardID(size_t i,id_t id) { _board_ids[i] = id; }
      void SetBoardIDs(std::vector<uint32_t> const& idvec) {_board_ids = idvec; }
      void SetCompressionScheme(uint32_t scheme) { _compression_scheme = scheme; }
    
    private:
      uint32_t _run_number;
      uint32_t _samples_per_channel;
      uint32_t _num_adc_bits;
      uint32_t _channels_per_board;
      uint32_t _num_boards;
      uint32_t _compression_scheme;
      std::vector<uint32_t> _board_ids;
  }; // end MetaData class

  class ICARUSProduceCompressed : public art::EDProducer
  {
    public:
      explicit ICARUSProduceCompressed(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produceForLabel(art::Event& evt, art::InputTag fFragmentLabel);
      void produce(art::Event& evt) override;

      uint64_t getFragmentWord(artdaq::Fragment const& f, size_t nWord)
      {
        if (nWord > f.dataSize())
          throw cet::exception("PhysCrateFragment:getFragmentWord")
           << "Asked for 64-bit word " << nWord << " of the fragment, but there are only " << f.dataSize() << " in the payload" << '\n';

        return *(f.dataBegin() + nWord);
      }
      uint16_t getA2795Word(artdaq::Fragment const& f, size_t nWord)
      {
        size_t nFragmentWord = nWord / 4;
        size_t shift = 16*(nWord % 4);
        uint64_t longWord = getFragmentWord(f, nFragmentWord);
        uint16_t shortWord = (longWord >> shift) & 0xFFFF;
        return shortWord;
      }
      void setFragmentWord(artdaq::Fragment& f, size_t nWord, uint64_t value)
      {
        *(f.dataBegin() + nWord) = value;
      }
      void setA2795Word(artdaq::Fragment& f, size_t nWord, uint16_t value)
      {
        uint64_t value64 = value;
        size_t nFragmentWord = nWord / 4;
        size_t shift = 16*(nWord % 4);
        uint64_t oldWord = getFragmentWord(f, nFragmentWord);
        uint64_t filterBlock = 0xFFFF;
        uint64_t filter = ~(filterBlock << shift);
        uint64_t newWord = (oldWord & filter) + ((value64 << shift) & ~filter);
        setFragmentWord(f, nWord / 4, newWord);
      }
      size_t SampleBytesFromKey(std::array<bool, 16> const& key)
      {
        size_t nCompressed = 0;
        for (auto const& bit : key)
          nCompressed += bit;

        return 128 - 6*nCompressed + 2*(nCompressed % 2);
      }
      uint16_t adc_val(artdaq::Fragment const& f, size_t b, size_t c, size_t s)
      {
        size_t nChannels = f.metadata<MetaData>()->channels_per_board();
        size_t nSamples  = f.metadata<MetaData>()->samples_per_channel();
        uint16_t const* boardData = reinterpret_cast<uint16_t const*>(f.dataBeginBytes() + (b+1)*28 + (b+1)*8 + b*(nChannels*nSamples + 4)*sizeof(uint16_t));
        return *(boardData + s*nChannels + c) & (~(1<<(f.metadata<MetaData>()->num_adc_bits()+1)));
      }
      artdaq::Fragment compressArtdaqFragment(artdaq::Fragment const & f);

    private:
      std::vector<art::InputTag> fFragmentLabelVec; // which Fragments are we pulling in?
      bool fDebug; // run with debugging?

  }; // end ICARUSProduceCompressed class

  //------------------------------------------------
  ICARUSProduceCompressed::ICARUSProduceCompressed(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void ICARUSProduceCompressed::reconfigure(fhicl::ParameterSet const& pset)
  {
    // What are we procesing
    fFragmentLabelVec = pset.get<std::vector<art::InputTag>>("FragmentLabelVec");
    fDebug = pset.get<bool>("RunDebugger", false);
    
    for (auto const& fragLabel : fFragmentLabelVec)
    {
      MF_LOG_VERBATIM("ICARUSProduceCompressed")
        << "Getting fragments " << fragLabel.label() << ":" << fragLabel.instance();

      // we make these things
      produces<std::vector<artdaq::Fragment>>(fragLabel.instance());
    }
  }

  //------------------------------------------------
  artdaq::Fragment ICARUSProduceCompressed::compressArtdaqFragment(artdaq::Fragment const & f)
  {
    // get boards, channels, samples in Fragment
    size_t nBoards   = f.metadata<MetaData>()->num_boards();
    size_t nChannels = f.metadata<MetaData>()->channels_per_board();
    size_t nSamples  = f.metadata<MetaData>()->samples_per_channel();

    // initialize the output fragment
    artdaq::Fragment compressed_fragment(f);

    // go through board/channel/sample and load up the data
    // data is stored in channel->sample->board
    // keep track of where we are reading from/writing to
    size_t   compressedDataOffset = 0;
    size_t uncompressedDataOffset = 0;
    std::vector<std::array<bool, 16>> compressionKeys;
    compressionKeys.resize(nBoards * nSamples);
    size_t totalDataTileSize = 0;
    for (size_t board = 0; board < nBoards; ++board)
    {
      // each board has a header...
      uint32_t boardDataTileSize = 36;
      for (size_t boardHeaderWord = 0; boardHeaderWord < 18; ++boardHeaderWord)
      {
        setA2795Word(compressed_fragment, compressedDataOffset, getA2795Word(f, uncompressedDataOffset));
        ++compressedDataOffset;
        ++uncompressedDataOffset;
      } // done with board header

      // since the 0th sample is used as a reference, store it as is
      compressionKeys[board*nSamples] = {};
      for (size_t channel = 0; channel < nChannels; ++channel)
      {
        uint16_t sampleADC  = (adc_val(f, board, channel, 0) & 0x0FFF);
        sampleADC += 0x8000;
        setA2795Word(compressed_fragment, compressedDataOffset, sampleADC);
        ++compressedDataOffset;
        ++uncompressedDataOffset;
      } // done with channel in sample 0
      boardDataTileSize += SampleBytesFromKey({});

      // from here on we store the difference between samples for each channel
      for (size_t sample = 1; sample < nSamples; ++sample)
      {
        bool oddCompressions = false;
        // work in blocks of 4 channels
        for (size_t channelBlock = 0; channelBlock < nChannels / 4; ++channelBlock)
        {
          int16_t adcDiff_0 = adc_val(f, board, 4*channelBlock + 0, sample) - adc_val(f, board, 4*channelBlock + 0, sample - 1);
          int16_t adcDiff_1 = adc_val(f, board, 4*channelBlock + 1, sample) - adc_val(f, board, 4*channelBlock + 1, sample - 1);
          int16_t adcDiff_2 = adc_val(f, board, 4*channelBlock + 2, sample) - adc_val(f, board, 4*channelBlock + 2, sample - 1);
          int16_t adcDiff_3 = adc_val(f, board, 4*channelBlock + 3, sample) - adc_val(f, board, 4*channelBlock + 3, sample - 1);

          bool isComp = (std::abs(adcDiff_0) < 8) &&
                        (std::abs(adcDiff_1) < 8) &&
                        (std::abs(adcDiff_2) < 8) &&
                        (std::abs(adcDiff_3) < 8);

          uncompressedDataOffset += 4;

          if (isComp)
          {
            compressionKeys[board*nSamples + sample][channelBlock] = true;
            oddCompressions = (not oddCompressions);
            uint16_t new_word  = (adcDiff_3 & 0x000F);
            new_word <<= 4;
                     new_word += (adcDiff_2 & 0x000F);
            new_word <<= 4;
                     new_word += (adcDiff_1 & 0x000F);
            new_word <<= 4;
                     new_word += (adcDiff_0 & 0x000F);
            setA2795Word(compressed_fragment, compressedDataOffset, new_word);
            ++compressedDataOffset;
          } else {
            compressionKeys[board*nSamples + sample][channelBlock] = false;
            setA2795Word(compressed_fragment, compressedDataOffset + 0, (adcDiff_0 & 0x0FFF) + 0x8000);
            setA2795Word(compressed_fragment, compressedDataOffset + 1, (adcDiff_1 & 0x0FFF) + 0x8000);
            setA2795Word(compressed_fragment, compressedDataOffset + 2, (adcDiff_2 & 0x0FFF) + 0x8000);
            setA2795Word(compressed_fragment, compressedDataOffset + 3, (adcDiff_3 & 0x0FFF) + 0x8000);
            compressedDataOffset += 4;
          }
        } // done with channels

        // if there are an odd number of words in the difference then there needs to be a spacer added
        if (oddCompressions)
        {
          setA2795Word(compressed_fragment, compressedDataOffset, 0);
          ++compressedDataOffset;
        }

        // add up the size of the differences
        boardDataTileSize += SampleBytesFromKey(compressionKeys[board*nSamples + sample]);
      } // done with samples

      // ...and each board has a trailer
      boardDataTileSize += 4*sizeof(uint16_t);
      for (size_t boardTrailerWord = 0; boardTrailerWord < 4; ++boardTrailerWord)
      {
        setA2795Word(compressed_fragment, compressedDataOffset, getA2795Word(f, uncompressedDataOffset));
        ++compressedDataOffset;
        ++uncompressedDataOffset;
      }

      // now that we know the size for the data tile, store that in the tile header
      // endianness is weird for this, hence the htonl...
      TileHeader* boardHeader = reinterpret_cast<TileHeader*>(compressed_fragment.dataBeginBytes() + totalDataTileSize);
      boardHeader->packSize = htonl(boardDataTileSize);
      totalDataTileSize += boardDataTileSize;
    } // done with boards

    // resize the fragment down to the compressed size
    compressed_fragment.resizeBytes(totalDataTileSize);

    // updated the metadata to reflect the compression
    compressed_fragment.metadata<MetaData>()->SetCompressionScheme(1);

    return compressed_fragment;
  }
  //------------------------------------------------
  void ICARUSProduceCompressed::produceForLabel(art::Event& evt, art::InputTag fFragmentLabel)
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
      size_t nBoards   = old_fragment.metadata<MetaData>()->num_boards();
      size_t nChannels = old_fragment.metadata<MetaData>()->channels_per_board();
      size_t nSamples  = old_fragment.metadata<MetaData>()->samples_per_channel();

      artdaq::Fragment new_fragment = compressArtdaqFragment(old_fragment);
      new_fragments->emplace_back(new_fragment);

      if (fDebug)
      {
        // store ADC vals
        std::vector<uint16_t> adcValues(nBoards*nChannels*nSamples, std::numeric_limits<uint16_t>::min());

        // loop to fill
        uint64_t fragWord = 0;
        for (size_t b = 0; b < nBoards; b++)
        {
          // skip board headers...
          fragWord += ((28 + 8) / sizeof(uint16_t));
          for (size_t s = 0; s < nSamples; s++)
          {
            size_t keyCount = 0;
            for (size_t bit = 0; bit < nChannels/4; bit++)
            {
              uint16_t word = getA2795Word(new_fragment, fragWord);
              bool isCompressed = (((word & 0xF000) != 0x8000));
              for (size_t cInSet = 0; cInSet < 4; ++ cInSet)
              {
                size_t c = 4*bit + cInSet;
                if (not isCompressed)
                {
                  uint16_t tempWord = getA2795Word(new_fragment, fragWord + cInSet);
                  int16_t twelveBitDiff = (tempWord & 0x0FFF);
                  bool isNeg = (twelveBitDiff >> 11) && (s != 0);
                  uint16_t prevSample = (s != 0) ? adcValues[b*nSamples*nChannels + c*nSamples + s - 1] : 0;
                  uint16_t currSample = (isNeg*0xF000 + twelveBitDiff + prevSample);
                  adcValues.at(b*nSamples*nChannels + c*nSamples + s) = currSample;
                } else {
                  int16_t fourBitDiff = (word >> (4*cInSet)) & 0x000F;
                  bool isNeg = (fourBitDiff >> 3);
                  uint16_t prevSample = (s != 0) ? adcValues[b*nSamples*nChannels + c*nSamples + s - 1] : 0;
                  uint16_t currSample = (isNeg*0xFFF0 + fourBitDiff + prevSample);
                  adcValues.at(b*nSamples*nChannels + c*nSamples + s) = currSample;
                }
              }
              fragWord += (isCompressed) ? 1 : 4;
              keyCount += isCompressed;
            }
            if ((keyCount % 2) == 1)
              fragWord += 1;
          }
          // ...and skip board trailers
          fragWord += 4;
        }
        for (size_t board = 0; board < nBoards; ++board)
          for (size_t sample = 0; sample < nSamples; ++sample)
            for (size_t channel = 0; channel < nChannels; ++channel)
            {
              uint16_t oldADC = adc_val(old_fragment, board, channel, sample);
              uint16_t newADC = adcValues.at(board*nSamples*nChannels + channel*nSamples + sample);
              if (oldADC != newADC)
                MF_LOG_VERBATIM("ICARUSProduceCompressed")
                  << "ERROR - ADC Mismatch in board " << board << ", sample " << sample << ", channel " << channel << '\n'
                  << "  old: " << oldADC << '\n'
                  << "  new: " << newADC;
            }
      }
    }

    // put the new fragments into the event
    evt.put(std::move(new_fragments), fFragmentLabel.instance());
  }

  //------------------------------------------------
  void ICARUSProduceCompressed::produce(art::Event& evt)
  {
    for (auto const& fragLabel : fFragmentLabelVec)
      this->produceForLabel(evt, fragLabel);
  }

  DEFINE_ART_MODULE(ICARUSProduceCompressed)
} // end namespace
