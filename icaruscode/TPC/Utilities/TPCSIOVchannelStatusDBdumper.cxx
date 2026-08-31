/**
 * @file   icaruscode/TPC/Utilities/TPCSIOVchannelStatusDBdumper.cxx
 * @brief  Provides the `TPCSIOVchannelStatusDBdumper` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanfird.edu)
 * @date   March 11, 2025
 * @see    icaruscode/TPC/Utilities/TPCSIOVchannelStatusDBdumper.h
 * 
 */

// library header
#include "icaruscode/TPC/Utilities/TPCSIOVchannelStatusDBdumper.h"

// LArSoft libraries
#include "larevt/CalibrationDBI/IOVData/IOVDataError.h"

// C/C++ standard libraries
#include <iomanip>
#include <string>


// -----------------------------------------------------------------------------
// ---  lariov::TPCSIOVchannelStatusDBdumper
// -----------------------------------------------------------------------------
lariov::TPCSIOVchannelStatusDBdumper::TPCSIOVchannelStatusDBdumper(
  Config const& config,
  lariov::SIOVChannelStatusProvider& channelStatus, unsigned int nChannels
)
  : fChannelStatus{ channelStatus }
  , fDumper{ config.TPCchannelStatusDBdumperConfig(), channelStatus, nChannels }
{
}


// -----------------------------------------------------------------------------
void lariov::TPCSIOVchannelStatusDBdumper::dumpTimestamp
  (std::ostream& out, DBTimeStamp_t timestamp) const
{
  fChannelStatus.UpdateTimeStamp(timestamp);
  
  out <<   "=== BEGIN TIMESTAMP: "
    << std::setw(19) << timestamp << " " << std::string(39, '=');
  
  dumpCurrentTimestamp(out);
  
  out << "\n=== END   TIMESTAMP: "
    << std::setw(19) << timestamp << " " << std::string(39, '=');
  
} // TPCSIOVchannelStatusDBdumper::dumpTimestamp()


// -----------------------------------------------------------------------------
void lariov::TPCSIOVchannelStatusDBdumper::dumpCurrentTimestamp
  (std::ostream& out) const
{
  
  fDumper.dump(out);
  
} // TPCSIOVchannelStatusDBdumper::dumpCurrentTimestamp()


// -----------------------------------------------------------------------------
std::ostream& lariov::operator<<
  (std::ostream& out, TPCSIOVchannelStatusDBdumper::timestampDump dumpInfo)
{
  dumpInfo.dump(out);
  return out;
}


// -----------------------------------------------------------------------------
