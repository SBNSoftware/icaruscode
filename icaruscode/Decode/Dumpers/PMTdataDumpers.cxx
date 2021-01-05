/**
 * @file   icaruscode/Decode/Dumpers/PMTdataDumpers.cxx
 * @brief  Functions to dump the content of PMT (V1730) data fragments to console
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 22, 2020
 * @see    icaruscode/Decode/Dumpers/PMTdataDumpers.h
 */

#ifndef ICARUSCODE_DECODE_DUMPERS_PMTdataDumpers_H
#define ICARUSCODE_DECODE_DUMPERS_PMTDUMPINGUTILS_H

// library header
#include "icaruscode/Decode/Dumpers/PMTdataDumpers.h"

// ICARUS libraries
#include "icaruscode/Decode/Dumpers/BinaryDumpUtils.h"

// framework libraries
#include <ostream>
#include <iomanip>
// #include <type_traits> // std::is_integral_v


// -----------------------------------------------------------------------------
std::ostream& sbndaq::operator<<
  (std::ostream& out, sbndaq::CAENV1730FragmentMetadata const& meta)
{
  using namespace sbn::data::dump;
  out << "Metadata: @" << ((void const*) &meta)
    << "\nChannels:  " << meta.nChannels
    << "\nSamples:   " << meta.nSamples
    << "\nTimeStamp: "
      << meta.timeStampSec << "." << zeropad(meta.timeStampNSec, 9U)
    << "\nTemps:    "
    ;
  for (std::size_t i = 0U; i < CAEN_V1730_MAX_CHANNELS; ++i)
    out << " [" << i << "] " << meta.chTemps[i];
  
  return out;
} // sbndaq::operator<< (sbndaq::CAENV1730FragmentMetadata)


// -----------------------------------------------------------------------------
std::ostream& sbndaq::operator<<
  (std::ostream& out, sbndaq::CAENV1730EventHeader const& header)
{
  
  using namespace sbn::data;
  out
    <<   "event header:    @" << ((void const*) &header)
    << "\nmarker:          " << std::hex << header.marker << std::dec;
  if (header.marker != 0x0A)
    out << " [UNEXPECTED (should be 0xA), DATA MAY BE CORRUPTED]";
  out
    << "\nsize:            " << header.eventSize
    << "\nVME64X board ID: " << std::hex << header.boardID << std::dec
    << "\nevent counter:   " << header.eventCounter
    << "\nchannel mask:    " << std::hex << header.ChannelMask() << std::dec
      << " " << dump::bin<16U>(header.ChannelMask())
    ;
  if constexpr (true) {
    // this branch of code is for with extended trigger time tag (ETTT) disabled
    out << "\ntime tag:        " << header.triggerTime() << " x 8ns";
    if (header.triggerTimeRollOver()) out << " (rollover!)";
    out << "\nLVDS pattern:    " << std::hex << header.pattern << std::dec
      << " " << dump::bin<16U>(header.pattern);
  }
  else {
    // this branch of code is for with extended trigger time tag (ETTT) enabled
    out << "\ntime tag (long): " << header.extendedTriggerTime() << " x 8ns";
  }
  out
    << "\nformat:          " << header.eventFormat
    << "\nfail:            " << std::boolalpha << static_cast<bool>(header.boardFail)
    << "\n(reserved):      " << header.reserved
    ;
  auto const excessData = header.eventSize - 4 * sizeof(uint32_t);
  if (excessData != 0)
    out << "\n  (more unknown data: " << excessData << " bytes!)";
  out << "\nheader dump:"
    << dump::hexdump(reinterpret_cast<uint16_t const*>(&header), sizeof(header) / sizeof(uint16_t));
  return out;
} // sbndaq::operator<<(sbndaq::CAENV1730EventHeader)


// -----------------------------------------------------------------------------
std::ostream& sbndaq::operator<<
  (std::ostream& out, sbndaq::CAENV1730Event const& event)
{
  
  out << "Event: @" << ((void const*) &event)
    << ", data block " << event.DataBlock << ", header:\n" << event.Header;
  return out;
  
} // sbndaq::operator<<(sbndaq::CAENV1730Event)


// -----------------------------------------------------------------------------
std::ostream& sbndaq::operator<<
  (std::ostream& out, sbndaq::CAENV1730Fragment const& frag)
{
  
  out << "Fragment @" << ((void const*) &frag)
    << ", data size: " << frag.DataPayloadSize() << " bytes ("
    << frag.ExpectedDataSize() << " expected)";
  
  if (frag.Verify()) out << " [verified]";
  else out << " [VERIFICATION FAILED]";
  
  sbndaq::CAENV1730FragmentMetadata const* meta = frag.Metadata();
  out << "\nMetadata: ";
  if (meta) out << *meta;
  else out << "n/a";

  sbndaq::CAENV1730Event const* event = frag.Event();
  out << "\n";
  if (event) {
    out << *event;
    
    std::size_t const dataSize
      = (frag.DataPayloadSize() - sizeof(event->Header)) / sizeof(uint16_t);
    
    // interface does not help much here, we need to rely on insider knowledge:
    // start from just after the event header
    auto const data = reinterpret_cast<uint16_t const*>(&(event->Header) + 1U);
    out << "\nData:" << sbn::data::dump::hexdump(data, dataSize);
    
  }
  else out << "Event: n/a";
  
  return out << std::endl;
} // sbndaq::operator<< (CAENV1730Fragment)


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_DUMPERS_PMTDUMPINGUTILS_H
