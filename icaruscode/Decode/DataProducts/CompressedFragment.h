// LArSoft includes
#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

// std includes
#include <bitset>
#include <arpa/inet.h>

//-------------------------------------------------------------------------------------------------------

namespace daq
{
  class compressedFragment : public icarus::PhysCrateFragment
  {
  public:
    // we need to generate the keys in the constructor
    compressedFragment(artdaq::Fragment const& f) : icarus::PhysCrateFragment(f),
                                                    comp_frag_(f),
                                                    //compressionKeys(GenerateKeys(f))
                                                    compressionKeys(f)
                                                    {
                                                    }
   
    size_t BoardBlockSize(size_t board);
    icarus::A2795DataBlock::data_t const* BoardData(uint16_t b);
    icarus::A2795DataBlock::data_t adc_val(size_t b, size_t c, size_t s);
    icarus::A2795DataBlock::data_t adc_val_unkeyed(size_t b, size_t c, size_t s);
    icarus::A2795DataBlock::data_t adc_val_recursive(size_t b, size_t c, size_t s);
    icarus::A2795DataBlock::data_t adc_val_debugger(size_t b, size_t c, size_t s);
    uint16_t const& CompressionKey(size_t b, size_t s) const 
                                    {
                                      size_t index = b*this->nSamplesPerChannel() + s;
                                      return compressionKeys._keys[index];
                                    };
  
  private:
    // because this is inheriting from icarus::PhysCrateFragment we don't have the raw fragment
    // make a copy because that's nonsense
    artdaq::Fragment const&      comp_frag_;
    //std::vector<uint16_t> const& compressionKeys;
    struct Keys
    {
      Keys(artdaq::Fragment const& f) : _keys(compressedFragment::GenerateKeys(f)) {};
      std::vector<uint16_t> const _keys;
    } compressionKeys;
    static size_t                SampleBytesFromKey(uint16_t const& key);
    size_t cumulativeSampleSize(size_t b, size_t s, size_t runningTotal = 0)
      {
        uint16_t const& key = this->CompressionKey(b, s);

        if (s == 0)
          return runningTotal + this->SampleBytesFromKey(key);

        return this->cumulativeSampleSize(b, s - 1, runningTotal + this->SampleBytesFromKey(key));
      };
    size_t cumulativeBoardSize(size_t b, size_t runningTotal = 0)
      {
        size_t nSamples = this->nSamplesPerChannel();

        if (b == 0)
          return this->cumulativeSampleSize(0, nSamples - 1, runningTotal);

        return this->cumulativeBoardSize(b - 1,
                                         this->cumulativeSampleSize(b, nSamples - 1, runningTotal));
      };
    static std::vector<uint16_t> GenerateKeys(artdaq::Fragment const& f);
    std::pair<icarus::A2795DataBlock::data_t, const icarus::A2795DataBlock::data_t*>
      adc_val_recursive_helper(size_t b, size_t c, size_t s, size_t sTarget,
        std::pair<icarus::A2795DataBlock::data_t, const icarus::A2795DataBlock::data_t*> pair);
  };

  //......................................................................
  size_t compressedFragment::SampleBytesFromKey(uint16_t const& key)
  {
    // each uncompressed set of channels uses 8 bytes
    // compressed sets use 2 bytes
    size_t nCompressed = std::bitset<16>(key).count();
    size_t blockSizeBytes = 128 - 6*nCompressed;
    
    // due to offsets, if there are an odd number of uncompressed sets
    // we need to have 2 bytes of filler
    if ((nCompressed % 2) == 1)
      blockSizeBytes += 2;

    return blockSizeBytes;
  }

  //......................................................................
  std::vector<uint16_t> compressedFragment::GenerateKeys(artdaq::Fragment const& f)
  {
    size_t nBoards   = f.metadata<icarus::PhysCrateFragmentMetadata>()->num_boards();
    size_t nSamples  = f.metadata<icarus::PhysCrateFragmentMetadata>()->samples_per_channel();

    std::vector<uint16_t> keys(nBoards*nSamples, 0);

    size_t cumulativePrevBlockSize = 0;
    for (size_t b = 0; b < nBoards; b++)
    {
      icarus::A2795DataBlock::data_t const* boardData
                                              = reinterpret_cast<icarus::A2795DataBlock::data_t const*>
                                                         ( f.dataBeginBytes()
                                                         + sizeof(icarus::PhysCrateDataTileHeader)
                                                         + sizeof(icarus::A2795DataBlock::Header)
                                                         + b*sizeof(icarus::PhysCrateDataTileHeader)
                                                         + b*4*sizeof(uint16_t)
                                                         + cumulativePrevBlockSize                   );

      size_t nWord = 0;
      for (size_t s = 0; s < nSamples; s++)
      {
        uint16_t key = 0;
        for (size_t bit = 0; bit < 16; bit++)
        {
          icarus::A2795DataBlock::data_t word = boardData[nWord];
          bool isCompressed = ((word & 0xF000) != 0x8000);
          key += (isCompressed << bit);
          nWord += (isCompressed) ? 1 : 4;
        }
        keys[b*nSamples + s] = key;
        nWord += (std::bitset<16>(key).count() % 2);
        cumulativePrevBlockSize += compressedFragment::SampleBytesFromKey(key);
      }
    }

    return keys;
  }

  //......................................................................
  size_t compressedFragment::BoardBlockSize(size_t board)
  {
    metadata()->BoardExists(board);

    size_t boardDataSize = 0;
    for (size_t sample = 0; sample < this->nSamplesPerChannel(); sample++)
    {
      boardDataSize += SampleBytesFromKey(this->CompressionKey(board, sample));
    }

    return boardDataSize;
  }

  //......................................................................
  icarus::A2795DataBlock::data_t const* compressedFragment::BoardData(uint16_t b)
  {
    metadata()->BoardExists(b);

    size_t cumulativePrevBlockSize = (b == 0) ? 0 : cumulativeBoardSize(b - 1);

    return ( reinterpret_cast< icarus::A2795DataBlock::data_t const*>
              ( this->comp_frag_.dataBeginBytes()
              + sizeof(icarus::PhysCrateDataTileHeader)
              + sizeof(icarus::A2795DataBlock::Header)
              + b*sizeof(icarus::PhysCrateDataTileHeader)
              + b*4*sizeof(uint16_t)
              + cumulativePrevBlockSize                  )           );
  }

  //......................................................................
  icarus::A2795DataBlock::data_t compressedFragment::adc_val(size_t b, size_t c, size_t s)
  {
    // initialie out output
    icarus::A2795DataBlock::data_t adcVal = 0;

    const size_t cSet = c / 4;
    const size_t channelInSet = c % 4;

    size_t sIndex = 0;
    for (size_t sPrime = 0; sPrime < s + 1; sPrime++)
    {
      if (sPrime != 0)
      {
        uint16_t prevSampleKey = this->CompressionKey(b, sPrime - 1);
        sIndex += (SampleBytesFromKey(prevSampleKey) / sizeof(icarus::A2795DataBlock::data_t));
      }

      uint16_t sampleKey = this->CompressionKey(b, sPrime);
      //size_t cIndex = 0;
      //for (size_t cSetPrime = 0; cSetPrime < cSet; cSetPrime++)
      //  cIndex += 4 - 3*((sampleKey >> cSetPrime) & 0x0001);

      size_t cIndex = 4*cSet - 3*std::bitset<16>(sampleKey % (1 << cSet)).count();

      // dealing with 4-bit and 12 bit numbers
      // we need to convert these to 16 bits
      // check the most significant bit to get sign
      // if it's negative, fill the empty bits with 1s
      bool isChannelSetCompressed = ((sampleKey >> cSet) & 0x0001) == 0x0001;
      if (isChannelSetCompressed)
      {
        int16_t dataSlice = *(BoardData(b)+sIndex+cIndex);
        int16_t fourBitDiff = (dataSlice >> (4*channelInSet)) & 0x000F;
        bool isNeg = (fourBitDiff >> 3);
        adcVal += isNeg*0xFFF0 + fourBitDiff; 
      } else {
        int16_t dataSlice = *(BoardData(b)+sIndex+cIndex+channelInSet);
        int16_t twelveBitDiff = dataSlice & 0x0FFF;
        bool isNeg = (sPrime != 0) && (twelveBitDiff >> 11);
        adcVal += isNeg*0xF000 + twelveBitDiff;
      }
    }
 
    return adcVal;
  }

  //......................................................................
  icarus::A2795DataBlock::data_t compressedFragment::adc_val_recursive(size_t b, size_t c, size_t s)
  {
    // call the helper function to do this recursively
    // but we only need the first value in the return
    return adc_val_recursive_helper(b, c, 0, s, std::make_pair(static_cast<icarus::A2795DataBlock::data_t>(0),BoardData(b))).first;
  }

  //......................................................................
  std::pair<icarus::A2795DataBlock::data_t, const icarus::A2795DataBlock::data_t*>
    compressedFragment::adc_val_recursive_helper(size_t b, size_t c, size_t s, size_t sTarget,
     std::pair<icarus::A2795DataBlock::data_t, const icarus::A2795DataBlock::data_t*> pair)
  {
    icarus::A2795DataBlock::data_t runningVal = pair.first;
    const icarus::A2795DataBlock::data_t* loc = pair.second;

    const size_t cSet = c / 4;
    const size_t channelInSet = c % 4;
    uint16_t sampleKey = this->CompressionKey(b, s);

    size_t cIndex = 4*cSet - 3*std::bitset<16>(sampleKey % (1 << cSet)).count();
    bool isChannelSetCompressed = ((sampleKey >> cSet) & 0x0001) == 0x0001;

    if (isChannelSetCompressed)
    {
      int16_t dataSlice = *(loc+cIndex);
      int16_t fourBitDiff = (dataSlice >> (4*channelInSet)) & 0x000F;
      bool isNeg = (fourBitDiff >> 3);
      runningVal += isNeg*0xFFF0 + fourBitDiff; 
    } else {
      int16_t dataSlice = *(loc+cIndex+channelInSet);
      int16_t twelveBitDiff = dataSlice & 0x0FFF;
      bool isNeg = (s != 0) && (twelveBitDiff >> 11);
      runningVal += isNeg*0xF000 + twelveBitDiff;
    }

    size_t keyCount = std::bitset<16>(sampleKey).count();
    size_t increment = 64 - 3*keyCount + (keyCount % 2);

    if (s == sTarget)
      return std::make_pair(runningVal, loc + increment);

    return this->adc_val_recursive_helper(b, c, s + 1, sTarget, std::make_pair(runningVal, loc + increment));
  }

  //......................................................................
  icarus::A2795DataBlock::data_t compressedFragment::adc_val_unkeyed(size_t b, size_t c, size_t s)
  {
    // this method is to gauge how long it takes to get the ADC values without keys
    // so loop over the whole block and count until you're there

    // initialie out output
    icarus::A2795DataBlock::data_t adcVal = 0;

    // what word are we on?
    size_t nWord = 0;

    for (size_t sPrime = 0; sPrime < s + 1; sPrime++)
    {
      size_t nUncompressed = 0;
      size_t cPrime = c % 4;
      while (cPrime < 64)
      {
        int16_t dataSlice = *(BoardData(b) + nWord);
        bool isCompressed = ((dataSlice & 0xF000) != 0x8000);
        bool isNeg;
        int16_t difference;
        if (isCompressed)
        {
          size_t bitsToShift = 4*(c % 4);
          int16_t fourBitDiff = (dataSlice >> bitsToShift) & 0x000F;
          isNeg = (fourBitDiff >> 3);
          difference = isNeg*0xFFF0 + fourBitDiff;
          nWord += 1;
        } else {
          nUncompressed++;
          dataSlice = *(BoardData(b) + nWord + (c % 4));
          int16_t twelveBitDiff = dataSlice & 0x0FFF;
          isNeg = (sPrime != 0) && (twelveBitDiff >> 11);
          difference = isNeg*0xF000 + twelveBitDiff;
          nWord += 4;
        }
        if (cPrime == c)
          adcVal += difference;
        cPrime += 4;
      }
      if ((nUncompressed % 2) == 1)
        nWord++;
    }
 
    return adcVal;
  }

  //......................................................................
  icarus::A2795DataBlock::data_t compressedFragment::adc_val_debugger(size_t b, size_t c, size_t s)
  {
    MF_LOG_VERBATIM("compressedFragment") << "----------------------------------------------------------------------------------------------------" << '\n'
      << "Getting ADC value for board " << b << ", channel " << c << ", sample " << s << '\n';

    // initialie out output
    icarus::A2795DataBlock::data_t adcVal = 0;

    size_t cSet = c / 4;
    size_t channelInSet = c % 4;

    MF_LOG_VERBATIM("compressedFragment") << "Channel set is " << cSet << ", with this channel being number " << channelInSet << " of the set" << '\n';

    size_t sIndex = 0;
    for (size_t sPrime = 0; sPrime < s + 1; sPrime++)
    {
      if (sPrime > 0)
      {
        uint16_t prevSampleKey = this->CompressionKey(b, sPrime - 1);
        sIndex += (SampleBytesFromKey(prevSampleKey) / sizeof(icarus::A2795DataBlock::data_t));
        if (sPrime == s) MF_LOG_VERBATIM("compressedFragment") << "Key for sample " << sPrime - 1 << " is " << std::bitset<16>(prevSampleKey) << '\n'
          << "This gives us " << SampleBytesFromKey(prevSampleKey) << " / "
                              << sizeof(icarus::A2795DataBlock::data_t) << " = "
                              << (SampleBytesFromKey(prevSampleKey) / sizeof(icarus::A2795DataBlock::data_t)) << " words for sample" << '\n';
      }
      if (sPrime == s)MF_LOG_VERBATIM("compressedFragment") << "For sPrime " << sPrime << ", sIndex = " << sIndex << '\n';

      uint16_t sampleKey = this->CompressionKey(b, sPrime);
      if (sPrime == s)MF_LOG_VERBATIM("compressedFragment") << "Key for sample " << sPrime << " is " << std::bitset<16>(sampleKey) << '\n';
      size_t cIndex = 0;
      for (size_t cSetPrime = 0; cSetPrime < cSet; cSetPrime++)
      { 
        cIndex += 4 - 3*((sampleKey >> cSetPrime) & 0x0001);
        if (sPrime == s)MF_LOG_VERBATIM("compressedFragment") << "Bit " << cSetPrime << " of key is " << ((sampleKey >> cSetPrime) & 0x0001)
          << " which means we increment cIndex by " << 4 - 3*((sampleKey >> cSetPrime) & 0x0001) << '\n';
      }
      if (sPrime == s)MF_LOG_VERBATIM("compressedFragment") << "cIndex is " << cIndex << '\n';

      // dealing with 4-bit and 12 bit numbers
      // we need to convert these to 16 bits
      // check the most significant bit to get sign
      // if it's negative, fill the empty bits with 1s
      int16_t difference;
      bool isChannelSetCompressed = ((sampleKey >> cSet) & 0x0001) == 0x0001;
      if (sPrime == s)MF_LOG_VERBATIM("compressedFragment")
        << "Looking at bit " << cSet << " of " << std::bitset<16>(sampleKey) << " are we compressed? "
        << (((sampleKey >> cSet) & 0x0001) == 0x0001) << '\n';
      if (isChannelSetCompressed)
      {
        int16_t dataSlice = *(BoardData(b)+sIndex+cIndex);
        int16_t  fourBitDiff = (dataSlice >> (4*channelInSet)) & 0x000F;
        bool isNeg = (fourBitDiff >> 3);
        difference = isNeg*0xFFF0 + fourBitDiff; 
        if (sPrime == s)MF_LOG_VERBATIM("compressedFragment")
          << "We are compressed." << '\n'
          << "The data slice is " << std::bitset<16>(dataSlice) << '\n'
          << "of which the difference is from section " << cSet << ": " << std::bitset<4>(fourBitDiff) << '\n'
          << "The difference is thus " << std::bitset<16>(isNeg*0xFFF0 + fourBitDiff) << " = " << difference << '\n';
      } else {
        int16_t dataSlice = *(BoardData(b)+sIndex+cIndex+channelInSet);
        int16_t twelveBitDiff = dataSlice & 0x0FFF;
        bool isNeg = (twelveBitDiff >> 11);
        difference = isNeg*0xF000 + twelveBitDiff;
        if (sPrime == s)MF_LOG_VERBATIM("compressedFragment")
          << "We are not compressed." << '\n'
          << "The data slice is " << std::bitset<16>(dataSlice) << '\n'
          << "The twelve bit differnce is then " << std::bitset<12>(twelveBitDiff)
          << "The difference is thus " << std::bitset<16>(isNeg*0xF000 + twelveBitDiff) << " = " << difference << '\n';
      }
      
      MF_LOG_VERBATIM("compressedFragment") << "Our running ADC value is " << adcVal << " + " << difference << " = " << adcVal + difference << '\n';
      adcVal += difference;
    }
 
    MF_LOG_VERBATIM("compressedFragment")
      << "The totak ADC value is thus " << adcVal << '\n'
      << "----------------------------------------------------------------------------------------------------";
    return adcVal;
  }
}
