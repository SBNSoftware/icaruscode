/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderChannelID.cxx
 * @brief  Simple object representing an ID for an adder channel.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 6, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderChannelID.h
 * 
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.h"

// C/C++ standard libraries
#include <ios> // std::dec, std::hex
#include <sstream>
#include <ostream>

// -----------------------------------------------------------------------------
namespace {
  
  template <typename Stream>
  decltype(auto) formatAdderChannel
    (Stream&& out, icarus::trigger::AdderChannelID const& channel)
  {
    std::ostringstream buffer;
    if (channel >= 0x1000)
      out << "0x" << std::hex << std::uppercase << raw::Channel_t{ channel } << std::dec;
    else
      out << raw::Channel_t{ channel };
    return std::forward<Stream>(out);
  } // formatAdderChannel()
  
} // local namespace


// -----------------------------------------------------------------------------
icarus::trigger::AdderChannelID::operator std::string() const {
  
  return formatAdderChannel(std::ostringstream{}, *this).str();
  
}


// -----------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, AdderChannelID channel)
{
  return formatAdderChannel(out, channel);
}


// -----------------------------------------------------------------------------
