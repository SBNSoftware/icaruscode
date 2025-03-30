/**
 * @file   icaruscode/TPC/Utilities/TPCchannelStatusDBdumper.cxx
 * @brief  Provides the `TPCchannelStatusDBdumper` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanfird.edu)
 * @date   March 11, 2025
 * @see    icaruscode/TPC/Utilities/TPCchannelStatusDBdumper.h
 * 
 */

// library header
#include "icaruscode/TPC/Utilities/TPCchannelStatusDBdumper.h"

// LArSoft libraries
#include "larevt/CalibrationDBI/IOVData/IOVDataError.h"

#include "larcorealg/CoreUtils/counter.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// framework libraries
#include "canvas/Utilities/Exception.h"

// C/C++ standard libraries
#include <string>


// -----------------------------------------------------------------------------
// ---  lariov::TPCchannelStatusDBdumper
// -----------------------------------------------------------------------------
lariov::TPCchannelStatusDBdumper::TPCchannelStatusDBdumper(
  Config const& config,
  lariov::ChannelStatusProvider& channelStatus, unsigned int nChannels
)
  : fPrintChannels{ config.PrintIndividualChannels() }
  , fPrintSummary { config.PrintSummary() }
  , fNChannels    { nChannels }
  , fChannelStatus{ channelStatus }
{
  if (!fPrintSummary && !fPrintChannels) {
    throw art::Exception{ art::errors::Configuration }
      << "All printing options have been disabled.";
  }
}


// -----------------------------------------------------------------------------
void lariov::TPCchannelStatusDBdumper::dump(std::ostream& out) const {
  unsigned int nErrors = 0;
  if (fPrintChannels) {
    for (raw::ChannelID_t const channel: util::counter(fNChannels)) {
      
      out << "\nCH=" << channel << ": ";
      try {
        if     (!fChannelStatus.IsPresent(channel)) out << "NOT PRESENT";
        else if (fChannelStatus.IsBad    (channel)) out << "BAD";
        else if (fChannelStatus.IsNoisy  (channel)) out << "NOISY";
        else if (fChannelStatus.IsGood   (channel)) out << "good";
        else
          out << "UNKNOWN STATUS (code=" << fChannelStatus.Status(channel) << ")";
      }
      catch(lariov::IOVDataError const& e) {
        out << "ERROR: " << e.what();
        ++nErrors;
      }
    }
  } // if print channels
  
  if (fPrintSummary) {
    if (fPrintChannels) out << "\n";
    
    try {
      lariov::ChannelStatusProvider::ChannelSet_t const& badChannels
        = fChannelStatus.BadChannels();
      out << "\nCounting " << badChannels.size() << " BAD channels";
      if (badChannels.empty()) out << ".";
      else {
        out << ":";
        for (raw::ChannelID_t const channel: badChannels) out << " " << channel;
      }
    }
    catch(lariov::IOVDataError const&) {
      out << "\nCount of NOISY channels not available.";
    }
    
    try {
      lariov::ChannelStatusProvider::ChannelSet_t const& noisyChannels
        = fChannelStatus.NoisyChannels();
      out << "\nCounting " << noisyChannels.size() << " NOISY channels";
      if (noisyChannels.empty()) out << ".";
      else {
        out << ":";
        for (raw::ChannelID_t const channel: noisyChannels) out << " " << channel;
      }
    }
    catch(lariov::IOVDataError const&) {
      out << "\nCount of BAD channels not available.";
    }
    
    try {
      lariov::ChannelStatusProvider::ChannelSet_t const& goodChannels
        = fChannelStatus.GoodChannels();
      out << "\nCounting " << goodChannels.size() << " good channels.";
    }
    catch(lariov::IOVDataError const&) {
      out << "\nCount of good channels not available.";
    }
    
    if (nErrors > 0) {
      out << "\nEncountered errors querying " << nErrors << "/" << fNChannels
        << " channels.";
    }
  } // if print summary
  
} // TPCchannelStatusDBdumper::dumpCurrentTimestamp()


// -----------------------------------------------------------------------------
std::ostream& lariov::operator<<
  (std::ostream& out, TPCchannelStatusDBdumper::dumpStruct dumpInfo)
{
  dumpInfo.dump(out);
  return out;
}

// -----------------------------------------------------------------------------
