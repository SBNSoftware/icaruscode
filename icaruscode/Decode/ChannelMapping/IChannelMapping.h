/**
 * @file icaruscode/Decode/ChannelMapping/IChannelMapping.h
 * @brief Interface definition for helpers handling the channel mapping.
 * @author T. Usher (usher@slac.stanford.edu)
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_ICHANNELMAPPING_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_ICHANNELMAPPING_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h" // icarusDB::RunPeriod
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapDataTypes.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// C/C++ standard libraries
#include <algorithm> // std::transform()
#include <array>
#include <string>
#include <cctype> // std::tolower()


// -----------------------------------------------------------------------------
namespace icarusDB { class IChannelMapping; }
/**
 * @brief Interface for helper classes extracting channel mapping from database.
 * 
 * This is currently used by both SQLite and PostgreSQL database backends,
 * and the service provider implementations take advantage of the similar
 * interface.
 */
class icarusDB::IChannelMapping {
    public:
  /**
   *  @brief  Virtual Destructor
   */
  virtual ~IChannelMapping() noexcept = default;
  
  /**
   *  @brief Prepares the tool to serve information about the specified run period.
   *  @param period the run period to be served
   *  @return whether the served values are going to be different than before
   * 
   * If this method returns `false`, the values from the previous call to it
   * are still current.
   * 
   * Run periods are defined in `icarusDB::RunPeriods`.
   */
  virtual bool SelectPeriod(icarusDB::RunPeriod period) = 0;
  
  
  // --- BEGIN --- TPC mapping -------------------------------------------------
  /// @name TPC mapping.
  /// @{
  
  /// Fill mapping between TPC Fragment IDs and the related crate and readout
  /// information.

  virtual int BuildTPCFragmentIDToReadoutIDMap
    (TPCFragmentIDToReadoutIDMap&) const = 0;
  
  /// Fill mapping between TPC readout boards and the channel information.
  virtual int BuildTPCReadoutBoardToChannelMap
    (TPCReadoutBoardToChannelMap&) const = 0;

  /// @}
  // --- END ----- TPC mapping -------------------------------------------------
  
  
  // --- BEGIN --- PMT mapping -------------------------------------------------
  /// @name PMT mapping.
  /// @{
  
  /// Fill the mapping between PMT Fragment IDs and the related crate and
  /// readout information.
  virtual int BuildPMTFragmentToDigitizerChannelMap
    (FragmentToDigitizerChannelMap&) const = 0;

  /// @}
  // --- END ----- PMT mapping -------------------------------------------------
  
  
  // --- BEGIN --- CRT mapping -------------------------------------------------
  /// @name CRT mapping.
  /// @{
  
  /// Fill mapping between side CRT hardware mac_address and the simulated
  /// mac_address.
  virtual int BuildCRTChannelIDToHWtoSimMacAddressPairMap
    (CRTChannelIDToHWtoSimMacAddressPairMap&) const = 0;

  /// Fill mapping between top CRT channel ID and the simulated mac_address.
  virtual int BuildTopCRTHWtoSimMacAddressPairMap
    (TopCRTHWtoSimMacAddressPairMap&) const = 0;
  
  /// Fill CRT calibration information.
  virtual int BuildSideCRTCalibrationMap
    (SideCRTChannelToCalibrationMap&) const = 0;
  
  /// @}
  // --- END ----- CRT mapping -------------------------------------------------
  
  
  // --- BEGIN --- Utility functions for the implementation --------------------
  
  static constexpr unsigned int ChannelsPerTPCreadoutBoard = 64;
  
  /// Turns and returns a string into lowercase.
  template <typename String>
  static void makeLowerCase(String& s) {
    using std::cbegin, std::cend, std::begin, std::end;
    std::transform(cbegin(s), cend(s), begin(s),
      [](unsigned char c){ return std::tolower(c); });
  }
  
  /// Returns the plane identified by `planeStr`: `0` to `2`; `3` if invalid.
  static unsigned int TPCplaneIdentifierToPlane(std::string planeStr)
    {
      using namespace std::string_literals;
      makeLowerCase(planeStr); // Make sure lower case... (sigh...)
      static std::array const PlaneNames
        = { "induction 1"s, "induction 2"s, "collection"s };
      for (auto const& [ plane, planeName ]: util::enumerate(PlaneNames))
        if (planeStr.find(planeName) != std::string::npos) return plane;
      return PlaneNames.size();
    }
  
  // --- END ----- Utility functions for the implementation --------------------
  
}; // icarusDB::IChannelMapping


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_ICHANNELMAPPING_H
