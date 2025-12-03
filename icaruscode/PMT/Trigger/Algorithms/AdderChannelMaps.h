/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.h
 * @brief  Utilities for assigning IDs to adder channels.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 6, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELMAPS_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELMAPS_H

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelID.h"
#include "icarusalg/Utilities/mfLoggingClass.h"

// LArSoft and framework libraries
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t
#include "cetlib_except/exception.h" // cet::exception

// C/C++ standard libraries
#include <iosfwd>
#include <limits>
#include <map>

// -----------------------------------------------------------------------------
// forward declarations
namespace icarus::trigger {
  
  struct AdderChannelInfo_t;
  
  class AdderChannelMapBuilder;
  
} // namespace icarus::trigger

namespace icarusDB { class IICARUSChannelMapProvider; }


// -----------------------------------------------------------------------------
struct icarus::trigger::AdderChannelInfo_t {
  using PMTchannelList_t = std::vector<raw::Channel_t>;
  
  /// Special value to indicate an invalid channel.
  static constexpr raw::Channel_t InvalidID
    = std::numeric_limits<raw::Channel_t>::max();
  
  AdderChannelID channel = InvalidID; ///< Adder channel ID.
  PMTchannelList_t PMTs; ///< List of PMT channels covered by this adder.
  
  // TODO: add board information as needed
  
}; // icarus::trigger::AdderChannelInfo_t

namespace icarus::trigger {
  
  /// Map of adder channel ID to adder channel information.
  using AdderChannelMap = std::map<AdderChannelID, AdderChannelInfo_t>;
  

  // ---------------------------------------------------------------------------
  /**
   * @brief Algorithm to build a adder channel map.
   * 
   * The current implementation requires that there is one adder channel per
   * PMT digitizer board, and assigns a channel ID according to which that board
   * is.
   * 
   * 
   * Channel numbering scheme
   * -------------------------
   * 
   * The convention is that the channel ID is the sequence number of the adder
   * board, defined based on the PMT the adder covers and following the LArSoft
   * convention where the PMT with lower _x_, then _z_, then _y_ get lower ID
   * and based on the lowest PMT covered
   * 
   * Mapping conventions
   * --------------------
   * 
   * As written above, the current implementation assumes that each adder
   * channel covers all the PMT channels on the same digitizer board.
   * If this is not the case, an exception will be thrown.
   * 
   */
  class AdderChannelMapBuilder: protected icarus::ns::util::mfLoggingClass {
    
      public:
    
    AdderChannelMapBuilder(std::string logCategory = "AdderChannelMapBuilder");
    
    /**
     * @brief Builds and returns the map.
     * @param PMTchannelMap PMT channel mapping provider
     * @param run run number or period describing the map; for diagnostics only
     * @return a map from adder channel to adder channel information
     * @throw cet::exception if inconsistency is detected
     */
    AdderChannelMap build
      (icarusDB::IICARUSChannelMapProvider const& PMTchannelMap, int run = 0)
      const;
    
  }; // AdderChannelMapBuilder

} // namespace icarus::trigger


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELMAPS_H
