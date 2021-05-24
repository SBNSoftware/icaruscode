/**
 * @file   icaruscode/Decode/DecoderTools/details/PMTDecoderUtils.h
 * @brief  Some helpers for PMT decoder tool.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 13, 2021
 */


// library header
#include "icaruscode/Decode/DecoderTools/details/PMTDecoderUtils.h"

// C/C++ standard libraries
#include <ostream>
#include <cassert>


// -----------------------------------------------------------------------------
std::ostream& daq::details::operator<<
  (std::ostream& out, BoardInfoLookup const& db)
{
  
  out << "Information on " << db.nBoardInfo() << " boards recorded:";
  for (BoardInfoLookup::BoardInfo_t const& boardInfo: db.allBoardInfo()) {
    assert(boardInfo.setup);
    out << "\n  board \"" << boardInfo.setup->name
      << "\" (fragment ID " << std::hex << boardInfo.fragmentID << std::dec
      << "): trigger delay " << boardInfo.setup->triggerDelay
      << ", TTT reset delay " << boardInfo.setup->TTTresetDelay
      << ", pre-trigger buffer length " << boardInfo.facts.preTriggerTime;
    if (boardInfo.config) {
      out << ", buffer " << boardInfo.config->bufferLength
        << " tick long, board ID " << boardInfo.config->boardID
        << " with " << boardInfo.config->nChannels << " channels";
    }
    else {
      out << " (no PMT configuration)";
    }
  } // for board info
  out << '\n';
  return out;
} // daq::details::operator<< (daq::details::BoardInfoLookup)


// -----------------------------------------------------------------------------
