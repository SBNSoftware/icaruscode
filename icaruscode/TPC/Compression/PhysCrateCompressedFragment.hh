#ifndef sbndaq_artdaq_core_Overlays_ICARUS_PhysCrateCompressedFragment_hh
#define sbndaq_artdaq_core_Overlays_ICARUS_PhysCrateCompressedFragment_hh

#include "artdaq-core/Data/Fragment.hh"
#include "cetlib_except/exception.h"

#include <bitset>
#include <iostream>
#include <vector>

#include "sbndaq-artdaq-core/Overlays/ICARUS/structures.h"

// Implementation of "PhysCrateCompressedFragment", an artdaq::Fragment overlay class
// used for ICARUS DAQ

namespace icarus {
  class PhysCrateCompressedFragment;
  std::ostream & operator << (std::ostream &, PhysCrateCompressedFragment const &);

  class PhysCrateCompressedFragmentMetadata;
  std::ostream & operator << (std::ostream &, PhysCrateCompressedFragmentMetadata const&);

  struct PhysCrateCompressedDataTileHeader;
  std::ostream & operator << (std::ostream &, struct PhysCrateCompressedDataTileHeader const&);

  struct A2795CompressedDataBlock;
  std::ostream & operator << (std::ostream &, struct A2795CompressedDataBlock const&);
}

//same structure that is found in structures.h
struct icarus::PhysCrateCompressedDataTileHeader{
  
  //'DATA'
  uint32_t token;
  
  //info1
  //this is an update based on looking at info from May2019
  //uint32_t stop_addr : 16;
  //uint32_t mode      : 16;
  //this is an update based on looking at info from May2019

  uint32_t info1;
  /*
  uint32_t ctrl_bufforg  : 4;
  uint32_t ctrl_status   : 4; //(CTRL_INHIBIT_TRG=0x8 | CTRL_ACQRUN=0x1)
  uint32_t ctrl_datamode : 4; //(CTRL_TTL_PROPAG=0x4 | CTRL_SRAM_TEST=0x2 | CTRL_TEST_PATTERN=0x1)
  uint32_t ctrl_unused1  : 4;
  uint32_t ctrl_unused2  : 4;
  uint32_t ctrl_unused3  : 4;
  uint32_t ctrl_unused4  : 4;
  uint32_t ctrl_testmode : 4; //( CTRL_TEST_EVENT=0x8 | CTRL_TEST_MODE=0x1)
  */

  //info2
  //uint32_t dead_peak_mask;
  //this is an update based on looking at info from May2019

  uint32_t info2;
  /*
  uint32_t slot_id : 4;
  uint32_t status  : 4; // (GBUSY=0x8 | BUSY=0x4 | DRDY=0x2 | RUNNING=0x1)
  uint32_t unused1 : 24;
  */

  //info3
  //uint32_t tv_trcr  : 16;
  //uint32_t abs_time : 16;
  //this is an update based on looking at info from May2019
  
  uint32_t info3;
  //uint32_t nev_stored;

  //timeinfo
  uint32_t timeinfo;
  
  //chID
  uint32_t pkt_fmt_ver  : 8;
  uint32_t crate_id     : 8;
  uint32_t board_id     : 8;
  uint32_t board_status : 8;
  
  //packSize
  uint32_t packSize;

  PhysCrateCompressedDataTileHeader(){};
  PhysCrateCompressedDataTileHeader( struct DataTile::Header const& );


  uint32_t CtrlStatus_BuffOrg()     const { return (info1 & 0x0000000F); }
  uint32_t CtrlStatus_AcqRun()      const { return (info1 & 0x00000010) >> 4; }
  uint32_t CtrlStatus_InhibitTrg()  const { return (info1 & 0x00000080) >> 4; }
  uint32_t CtrlStatus_TestPattern() const { return (info1 & 0x00000100) >> 8; }
  uint32_t CtrlStatus_SRAMTest()    const { return (info1 & 0x00000200) >> 9; }
  uint32_t CtrlStatus_TTLPropag()   const { return (info1 & 0x00000400) >> 10; }
  uint32_t CtrlStatus_TestMode()    const { return (info1 & 0x10000000) >> 28; }
  uint32_t CtrlStatus_TestEvent()   const { return (info1 & 0x80000000) >> 31; }  

  uint32_t StatusReg_SlotID()       const { return (info2 & 0x0000000F); }
  uint32_t StatusReg_Running()      const { return (info2 & 0x00000010) >> 4; }
  uint32_t StatusReg_DataReady()    const { return (info2 & 0x00000020) >> 5; }
  uint32_t StatusReg_Busy()         const { return (info2 & 0x00000040) >> 6; }
  uint32_t StatusReg_GBusy()        const { return (info2 & 0x00000080) >> 7; }


};
struct icarus::A2795CompressedDataBlock{

  typedef uint32_t header_t;
  typedef uint16_t data_t;

  typedef struct {
    header_t event_number : 24;
    header_t unused1      :  8;
    header_t time_stamp;
  } Header;

  Header  header;
  data_t* data;
};

class icarus::PhysCrateCompressedFragmentMetadata {
  
public:
  
  typedef uint32_t data_t; //fundamental unit of metadata
  typedef data_t   id_t; //type for board IDs
  
  PhysCrateCompressedFragmentMetadata(){}
  PhysCrateCompressedFragmentMetadata(data_t run_number,
			    data_t n_boards,
			    data_t channels_per_board,
			    data_t samples_per_channel,
			    data_t adcs_per_sample,
			    data_t compression,
			   std::vector<id_t> const& idvec)
  { 
    _run_number = run_number;
    _num_boards = n_boards;
    _channels_per_board = channels_per_board;
    _samples_per_channel = samples_per_channel;
    _num_adc_bits = adcs_per_sample;
    _compression_scheme = compression;
    SetBoardIDs(idvec);
  }

  data_t const& run_number() const { return _run_number; }
  data_t const& samples_per_channel() const { return _samples_per_channel; }
  data_t const& num_adc_bits() const { return _num_adc_bits; }
  data_t const& channels_per_board() const { return _channels_per_board; }
  data_t const& num_boards() const { return _num_boards; }
  data_t const& compression_scheme() const { return _compression_scheme; }

  id_t   const& board_id(size_t i) const
  { BoardExists(i); return _board_ids[i]; }

  void  SetBoardID(size_t i,id_t id)
  { BoardExists(i); _board_ids[i] = id; }
  void  SetBoardIDs(std::vector<id_t> const& idvec)
  { CheckNBoards(idvec.size()); _board_ids = idvec; } 

  void BoardExists(size_t i) const;
  void CheckNBoards(size_t i) const;

  size_t ExpectedDataSize() const 
  { return _num_boards*(sizeof(PhysCrateCompressedDataTileHeader) +
			sizeof(A2795CompressedDataBlock::Header) +
			_channels_per_board*_samples_per_channel*sizeof(A2795CompressedDataBlock::data_t) +
			2*sizeof(uint32_t)); }


private:

  data_t _run_number;
  data_t _samples_per_channel;  
  data_t _num_adc_bits;
  data_t _channels_per_board;
  data_t _num_boards;
  data_t _compression_scheme;
  std::vector<id_t> _board_ids;
  
  void UpdateBoardIDSize(){ _board_ids.resize(_num_boards); }

};

class icarus::PhysCrateCompressedFragment {

  public:

  PhysCrateCompressedFragment(artdaq::Fragment const & f) : artdaq_Fragment_(f), 
                                                  compressionKeys_(f) {}

  PhysCrateCompressedFragment(artdaq::Fragment const & f, bool const & compressionSwitch)
    : PhysCrateCompressedFragment(this->fragmentSwitch(f, compressionSwitch)) {}

  PhysCrateCompressedFragmentMetadata const * metadata() const { return artdaq_Fragment_.metadata<PhysCrateCompressedFragmentMetadata>(); }

  size_t RunNumber() const { return metadata()->run_number(); }
  size_t nBoards() const { return metadata()->num_boards(); }
  size_t nChannels() const { return metadata()->num_boards()*metadata()->channels_per_board(); }
  size_t nSamplesPerChannel() const { return metadata()->samples_per_channel(); }
  size_t nChannelsPerBoard() const { return metadata()->channels_per_board(); }
  size_t CompressionScheme() const { return metadata()->compression_scheme(); }

  bool   isCompressed() const { return (CompressionScheme()!=0); } // check that 0 in the md is compressed

  size_t DataPayloadSize() const { return artdaq_Fragment_.dataSizeBytes(); }

  PhysCrateCompressedDataTileHeader const * DataTileHeader(uint16_t b=0) const;
  size_t DataTileHeaderLocation(uint16_t b=0) const;

  size_t BoardBlockSize() const
  { return sizeof(A2795CompressedDataBlock::Header)+nChannelsPerBoard()*nSamplesPerChannel()*sizeof(A2795CompressedDataBlock::data_t); }

  A2795CompressedDataBlock           const* BoardDataBlock(uint16_t b=0) const;
  A2795CompressedDataBlock::Header   const& BoardHeader(uint16_t b=0) const;
  A2795CompressedDataBlock::header_t        BoardEventNumber(uint16_t b=0) const;
  A2795CompressedDataBlock::header_t        BoardTimeStamp(uint16_t b=0) const;
  A2795CompressedDataBlock::data_t   const* BoardData(uint16_t b=0) const;

  A2795CompressedDataBlock::data_t adc_val(size_t b,size_t c,size_t s) const;

  bool Verify() const;

  uint16_t const& CompressionKey(size_t b, size_t s) const
                                  {
                                    size_t index = b*this->nSamplesPerChannel() + s;
                                    return compressionKeys_.keys_[index];
                                  }

  typedef std::pair<A2795CompressedDataBlock::data_t, const A2795CompressedDataBlock::data_t*> recursionPair;

  PhysCrateCompressedFragment makeCompressedFragment()   const { return PhysCrateCompressedFragment(  compressArtdaqFragment(artdaq_Fragment_)); }
  PhysCrateCompressedFragment makeUncompressedFragment() const { return PhysCrateCompressedFragment(decompressArtdaqFragment(artdaq_Fragment_)); }

  static artdaq::Fragment   compressArtdaqFragment(artdaq::Fragment const & f);
  static artdaq::Fragment decompressArtdaqFragment(artdaq::Fragment const & f);
  static artdaq::Fragment fragmentSwitch(artdaq::Fragment const & f, bool const & compressionSwitch)
  {
    return (compressionSwitch) ? compressArtdaqFragment(f) : decompressArtdaqFragment(f);
  }

private:

  artdaq::Fragment const & artdaq_Fragment_;

  void   throwIfCompressed() const;

  // here are things helpful for the comrpessed fragments
  struct Keys
  {
     Keys(artdaq::Fragment const& f) : keys_(GenerateKeys(f)) {};
     std::vector<uint16_t> const keys_;
  } compressionKeys_;

  static size_t SampleBytesFromKey(uint16_t const& key)
  {
    size_t nCompressed = std::bitset<16>(key).count();
    return 128 - 6*nCompressed + ((nCompressed % 2) == 1)*2;
  }

  size_t cumulativeSampleSize(size_t b, size_t s, size_t runningTotal = 0) const
  {
    // tail recursive function to total the bytes used for each sample
    // requires the compression keys to function
    uint16_t const& key = this->CompressionKey(b, s);

    if (s == 0)
      return runningTotal + this->SampleBytesFromKey(key);

    return this->cumulativeSampleSize(b, s - 1, runningTotal + this->SampleBytesFromKey(key));
  }

  size_t cumulativeBoardSize(size_t b, size_t runningTotal = 0) const
  {
    // tail recursive function to total the bytes used for each board
    // each board is made up of nSamples samples
    // requires the compression keys to function
    size_t nSamples = this->nSamplesPerChannel();
    
    if (b == 0)
      return this->cumulativeSampleSize(0, nSamples - 1, runningTotal);

    return this->cumulativeBoardSize(b - 1, this->cumulativeSampleSize(b, nSamples - 1, runningTotal));
  }

  static std::vector<uint16_t> GenerateKeys(artdaq::Fragment const& f);

  recursionPair adc_val_recursive_helper(size_t b, size_t c, size_t s, size_t sTarget, recursionPair pair) const;
};


#endif /* sbndaq_artdaq_core_Overlays_ToyFragment_hh */
