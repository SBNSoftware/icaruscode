/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapDataTypes.h
 * @brief  Data types used for caching channel maps.
 * @author T. Usher (refactoring by G. Petrillo, petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPDATATYPES_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPDATATYPES_H

// C/C++ standard libraries
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include <limits>
#include <cstddef> // std::size_t


namespace icarusDB {
  
  // --- BEGIN --- ICARUS channel mapping data definitions ---------------------
  /// @name ICARUS channel mapping data definitions.
  /// @{

  
  /// @name Data structures for mapping between TPC fragment IDs and the related
  ///       crate and readout information.
  /// @{
  
  /// Sequence of TPC readout board ID.
  using ReadoutIDVec = std::vector<unsigned int>;
  
  /// Name of the TPC readout crate, and list of IDs of hosted readout boards.
  using CrateNameReadoutIDPair = std::pair<std::string, ReadoutIDVec>;
  
  /// Map from fragment ID to crate information (name and readout board IDs).
  using TPCFragmentIDToReadoutIDMap
    = std::map<unsigned int, CrateNameReadoutIDPair>;
  
  /// @}
  
  
  /// @name Data structures for mapping between TPC readout board IDs and the
  ///       related crate and readout information.
  /// @{
  
  /// A TPC channel ID/TPC plane number pair.
  using ChannelPlanePair = std::pair<unsigned int, unsigned int>;
  
  /// A sequence of TPC channel IDs paired with their TPC plane number.
  using ChannelPlanePairVec = std::vector<ChannelPlanePair>;

  /// TPC readout board slot number and all channels in the board with plane.
  using SlotChannelVecPair = std::pair<unsigned int, ChannelPlanePairVec>;
  
  /// Map from TPC readout board ID to board information (slot and channels).
  using TPCReadoutBoardToChannelMap
    = std::map<unsigned int, SlotChannelVecPair>;
  
  /// @}
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /// @name Data structures for PMT channel mapping.
  /// @{
  
  /// Collection of information pertaining a single PMT.
  struct PMTChannelInfo_t {
    static constexpr int NoDigitizerChannel
      = std::numeric_limits<unsigned int>::max();
    static constexpr unsigned int NoChannelID
      = std::numeric_limits<unsigned int>::max();
    static constexpr unsigned int NoLaserChannel
      = std::numeric_limits<unsigned int>::max();
    
    /// Number of the channel within its digitizer.
    unsigned int digitizerChannelNo = NoDigitizerChannel;
    
    /// LArSoft channel ID.
    unsigned int channelID = NoChannelID;
    
    /// Number of laser channel shining into this PMT.
    unsigned int laserChannelNo = NoLaserChannel;
    
    /// Sorting by channel ID.
    constexpr bool operator< (PMTChannelInfo_t const& other) const noexcept
      { return channelID < other.channelID; }
    
  }; // PMTChannelInfo_t
  
  /// A sequence of PMT channel information for a single digitizer.
  using PMTdigitizerInfoVec = std::vector<PMTChannelInfo_t>;

  /// Map from PMT fragment ID to information on all its digitizer channels.
  using PMTFragmentToDigitizerChannelMap
    = std::map<unsigned int, PMTdigitizerInfoVec>;
  
  /// @}
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /// @name Data structures for CRT channel mapping.
  /// @{
  
  /// Hardware and Sim Mac addresses.
  using CRTHWtoSimMacAddressPair = std::pair<unsigned int, unsigned int>;
  
  /// Map of channel ID to Mac address pair for side CRT.
  using CRTChannelIDToHWtoSimMacAddressPairMap
    = std::map<unsigned int, CRTHWtoSimMacAddressPair>;

  /// Map of hardware Mac address to Sim Mac address for top CRT.
  using TopCRTHWtoSimMacAddressPairMap = std::map<unsigned int, unsigned int>;
  
  // @}
  
  
  /// @name Data structures for CRT calibration information.
  /// @{
  
  /// Pair of channel ID and hardware Mac address for CRT.
  using SideCRTMac5ToChannelPair = std::pair<unsigned int, unsigned int>;
  
  /// CRT channel calibration information: gain and pedestal.
  using SideCRTGainToPedPair = std::pair<double, double>;
  
  /// Map of channel/hardware address to its calibration information.
  using SideCRTChannelToCalibrationMap
    = std::map<SideCRTMac5ToChannelPair, SideCRTGainToPedPair>;
  
  /// @}
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  /// @}
  // --- END ----- ICARUS channel mapping data definitions ---------------------
  
} // namespace icarusDB


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPDATATYPES_H
