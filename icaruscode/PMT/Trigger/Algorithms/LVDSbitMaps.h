/**
 * @file   icaruscode/PMT/Trigger/Algorithms/LVDSbitMaps.h
 * @brief  Utilities for mapping CAEN V1730B LVDS channels to PMT channels.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 10, 2024
 * @see    icaruscode/PMT/Trigger/Algorithms/LVDSbitMaps.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_LVDSBITMAPS_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_LVDSBITMAPS_H

// LArSoft and framework libraries
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t
#include "cetlib_except/coded_exception.h"

// C/C++ standard libraries
#include <ostream>
#include <iomanip> // std::setw(), std::setfill()
#include <initializer_list>
#include <vector>
#include <string>
#include <array>
#include <limits>

// -----------------------------------------------------------------------------
// forward declarations
namespace icarus::trigger {
  
  struct PMTpairBitID;
  struct AdderBitID;
  struct LVDSHWbitID;
  struct LVDSinfo_t;
  struct AdderInfo_t;
  
  class LVDSbitMaps;
  
  std::ostream& operator<< (std::ostream& out, PMTpairBitID const& bit);
  std::ostream& operator<< (std::ostream& out, LVDSHWbitID const& bit);
  
} // namespace icarus::trigger

namespace icarusDB { class IICARUSChannelMapProvider; }


// -----------------------------------------------------------------------------
/// Identifier of a "logical" LVDS bit for PMT pairs
/// following `sbn::ExtraTriggerInfo` convention.
namespace icarus::trigger::details {
  template <std::size_t MaxBits = 64> struct PMTwallBitID;
  template <std::size_t MaxBits>
  std::ostream& operator<<
    (std::ostream& out, PMTwallBitID<MaxBits> const& bit);
}
template <std::size_t MaxBits /* = 64 */>
struct icarus::trigger::details::PMTwallBitID {
  
  using BitID_t = PMTwallBitID<MaxBits>; ///< This ID type.
  using StatusBit_t = unsigned short int; ///< Type for bit in the status word.
  using Index_t = unsigned short int; ///< Type for bit ID index.
  
  /// How many cryostats in the detector.
  static constexpr std::size_t NCryostats = 2;
  
  /// How many PMT walls per cryostat. Each has one status word.
  static constexpr std::size_t NPMTwallsPerCryostat = 2;
  
  /// How many bits in a status word.
  static constexpr unsigned short int NStatusWordBits = MaxBits;
  
  
  /// Special value to identify the absence of cryostat information.
  static constexpr std::size_t NoCryostat
    = std::numeric_limits<std::size_t>::max();
  
  /// Special value to identify the absence of PMT wall information.
  static constexpr std::size_t NoPMTwall
    = std::numeric_limits<std::size_t>::max();
  
  /// Special value to identify the absence of a status bit.
  static constexpr StatusBit_t NoStatusBit
    = std::numeric_limits<StatusBit_t>::max();
  
  /// Special value for an invalid bit index.
  static constexpr Index_t InvalidIndex = std::numeric_limits<Index_t>::max();
  
  
  // --- BEGIN ---  Data members  --------------------------------------------
  
  /// The cryostat the bit comes from (`sbn::ExtraTriggerInfo` convention).
  std::size_t cryostat = NoCryostat;
  
  /// The PMT wall the bit comes from (`sbn::ExtraTriggerInfo` convention).
  std::size_t PMTwall = NoPMTwall;
  
  /// The number of bit in the status word (`0` is the least significant one).
  StatusBit_t statusBit = NoStatusBit;
  
  // --- END -----  Data members  --------------------------------------------
  
  /// Returns if the cryostat is set; ID is valid even if this is `false`.
  constexpr bool hasCryostat() const { return cryostat != NoCryostat; }
  
  /// Returns if the PMT wall is set; ID is valid even if this is `false`.
  constexpr bool hasPMTwall() const { return PMTwall != NoPMTwall; }
  
  /// Returns an "unique" index for this information.
  /// Supports the absence of PMT wall or cryostat (not of bit).
  constexpr Index_t index() const
    {
      if (!isValid()) return InvalidIndex;
      std::size_t const nWalls = hasPMTwall()
        ? (PMTwall + (hasCryostat()? cryostat: 0) * NPMTwallsPerCryostat): 0;
      return statusBit + nWalls * NStatusWordBits;
    }
  
  
  /// Returns whether this ID is valid
  /// (does not require `hasCryostat()` nor `hasPMTwall()`).
  constexpr bool isValid() const { return statusBit != NoStatusBit; }
  
  /// Returns whether this ID is valid
  /// (does not require `hasCryostat()` nor `hasPMTwall()`).
  constexpr operator bool() const { return isValid(); }
  
  
  /// Rebuilds an ID object from its index.
  static constexpr BitID_t fromIndex(Index_t index);
  
  /// Returns the number of possible indices (`0` to `NIndices() - 1`).
  static constexpr Index_t NIndices()
    { return NCryostats * NPMTwallsPerCryostat * NStatusWordBits; }
  
}; // icarus::trigger::details::PMTwallBitID


// -----------------------------------------------------------------------------
/// Identifier of a "logical" LVDS bit for PMT pairs
/// following `sbn::ExtraTriggerInfo` convention.
struct icarus::trigger::PMTpairBitID: public details::PMTwallBitID<64> {
  
  using BaseID_t = details::PMTwallBitID<64>;
  
  /// Rebuilds an ID object from its index.
  static PMTpairBitID fromIndex(Index_t index)
    { return { BaseID_t::fromIndex(index) }; }
  
}; // icarus::trigger::PMTpairBitID


// -----------------------------------------------------------------------------
/// Identifier of a "logical" LVDS bit for adders
/// following `sbn::ExtraTriggerInfo` convention.
struct icarus::trigger::AdderBitID: public details::PMTwallBitID<16> {
  
  using BaseID_t = details::PMTwallBitID<16>;
  
  /// Rebuilds an ID object from its index.
  static AdderBitID fromIndex(Index_t index)
    { return { BaseID_t::fromIndex(index) }; }
  
}; // icarus::trigger::AdderBitID


// -----------------------------------------------------------------------------
/// Identifier of a "hardware" LVDS bit featured in the PMT channel mapping
/// databases.
struct icarus::trigger::LVDSHWbitID {
  
  using Connector_t = unsigned short int; ///< Type for connector number.
  using ConnectorBit_t = unsigned short int; ///< Type for connector bit.
  using Index_t = unsigned short int; ///< Type for LVDS bit ID index.
  
  /// How many cryostats in the detector.
  static constexpr unsigned short int NCryostats = PMTpairBitID::NCryostats;
  
  /// How many 32-bit connectors come out of each cryostat.
  static constexpr unsigned short int NConnectorsPerCryostat = 4;
  
  /// Special value to identify the absence of cryostat information.
  static constexpr unsigned short int NoCryostat
    = std::numeric_limits<unsigned short int>::max();
  
  /// Special value to identify the absence of a connector.
  static constexpr Connector_t NoConnector
    = std::numeric_limits<Connector_t>::max();
  
  /// Special value to identify the absence of a connector bit.
  static constexpr ConnectorBit_t NoConnectorBit
    = std::numeric_limits<ConnectorBit_t>::max();
  
  /// Special value for an invalid bit index.
  static constexpr Index_t InvalidIndex = std::numeric_limits<Index_t>::max();
  
  
  // --- BEGIN ---  Data members  --------------------------------------------
  
  /// The cryostat the bit comes from (`sbn::ExtraTriggerInfo` convention).
  unsigned short int cryostat = NoCryostat;
  
  /// The number of connector within the appropriate wall.
  Connector_t connector = NoConnector;
  
  /// The bit index within the 32-bit connector word (0 = least significant).
  ConnectorBit_t bit = NoConnectorBit;
  
  // --- END -----  Data members  --------------------------------------------
  
  /// Returns if the cryostat is set; ID is valid even if this is `false`.
  constexpr bool hasCryostat() const { return cryostat != NoCryostat; }
  
  /// Returns an "unique" index for this information.
  constexpr Index_t index() const
    {
      return isValid()
        ? ((hasCryostat()? (cryostat * NConnectorsPerCryostat * 32): 0)
          + connector * 32 + bit)
        : InvalidIndex;
    }
  
  
  /// Returns whether this ID is valid (does not require `hasCryostat()`).
  constexpr bool isValid() const
    { return (connector != NoConnector) && (bit != NoConnectorBit); }
  
  /// Returns whether this ID is valid (does not require `hasCryostat()`).
  constexpr operator bool() const { return isValid(); }
  
  
  /// Rebuilds an ID object from its index.
  static LVDSHWbitID fromIndex(Index_t index);
  
  /// Returns the number of possible indices (`0` to `NIndices() - 1`).
  static constexpr Index_t NIndices()
    { return NCryostats * NConnectorsPerCryostat * 32; }
  
}; // icarus::trigger::LVDSHWbitID


// -----------------------------------------------------------------------------
/// Information about a PMT pair: source bit and connected channels.
struct icarus::trigger::LVDSinfo_t {
  
  PMTpairBitID ID; ///< Identifier of this bit.
  
  LVDSHWbitID source; ///< Identifier of the hardware LVDS bit feeding this one.
  
  std::vector<raw::Channel_t> channels; ///< PMT channels covered by this bit.
  
  /// Returns whether this information is fully valid (both channel and ID).
  operator bool() const { return source; }
  
  /// Order operator: sorted after the ID.
  bool operator< (LVDSinfo_t const& other) const { return ID < other.ID; }
  
}; // icarus::trigger::LVDSinfo_t


// -----------------------------------------------------------------------------
/// Description of a (logical) adder bit and additional information.
struct icarus::trigger::AdderInfo_t {
  
  AdderBitID ID; ///< Identifier of this bit.
  
  LVDSHWbitID source; ///< Identifier of the hardware LVDS bit feeding this one.
  
  std::vector<raw::Channel_t> channels; ///< PMT channels covered by this bit.
  
  /// Returns whether this information is fully valid (both channel and ID).
  operator bool() const { return source; }
  
  /// Order operator: sorted after the source index.
  bool operator< (AdderInfo_t const& other) const { return ID < other.ID; }
  
}; // icarus::trigger::AdderInfo_t


// -----------------------------------------------------------------------------
/**
 * @brief Utility to build and maintain ICARUS trigger LVDS bit maps.
 * 
 * This object caches the relations between the LVDS bits. The maps include:
 *  * from logic LVDS bit ID to hardware LVDS bit ID (see below).
 * 
 * Upon request, the following mappings may be added:
 *  * from PMT channel to LVDS bit ID (logic), and vice versa.
 * 
 * 
 * Maps are build based on the information from the PMT hardware channel
 * mapping.
 * 
 * 
 * 
 * Mapping conventions
 * --------------------
 * 
 * The logic mapping is described in `sbn::ExtraTriggerInfo::LVDSstatus`.
 * Compared with that mapping, the PMT pair bit ID are assigned so that
 *  * the index of the ID, modulo 64 bits, matches the `LVDSstatus` words;
 *  * cryostats and PMT walls of the IDs are honoured, so that they can be used
 *    to directly address `LVDSstatus`.
 * 
 * Likewise, the adder bit ID are assigned so that
 *  * the index of the ID, modulo 6 bits, matches the `sectorStatus` words;
 *  * cryostats and PMT walls of the IDs are honoured, so that they can be used
 *    to directly address `sectorStatus`.
 * 
 * This results in the following mapping, "first bit" meaning the most
 * significant.
 * 
 * PMT pairs:
 * * east cryostat, east wall: `0x00aAbBcC'00dDeEfF`, channel `0` the most
 *   significant bit in `a` (bit #55) and channel `45` the one in `d` (bit 23);
 * * east cryostat, west wall: `0x00aAbBcC'00dDeEfF`, channel `90` the most
 *   significant bit in `a` (bit #55) and channel `135` the one in `d` (bit 23);
 * * west cryostat, east wall: `0x00aAbBcC'00dDeEfF`, channel `180` the most
 *   significant bit in `a` (bit #55) and channel `225` the one in `d` (bit 23);
 * * west cryostat, west wall: `0x00aAbBcC'00dDeEfF`, channel `270` the most
 *   significant bit in `a` (bit #55) and channel `315` the one in `d` (bit 23).
 * 
 * Adders:
 * * east cryostat, east wall: `0b00000000'00abcdef`, with `a` to `f` covering
 *   adders starting with channel `0`, `15`, 30`, `45`, 60` and `75`;
 * * east cryostat, west wall: `0b00000000'00abcdef`, with `a` to `f` covering
 *   adders starting with channel `90`, `105`, 120`, `135`, 150` and `165`;
 * * west cryostat, east wall: `0b00000000'00abcdef`, with `a` to `f` covering
 *   adders starting with channel `180`, `195`, 210`, `225`, 240` and `255`;
 * * west cryostat, west wall: `0b00000000'00abcdef`, with `a` to `f` covering
 *   adders starting with channel `270`, `285`, 300`, `315`, 330` and `345`.
 * 
 * 
 * Naming conventions
 * -------------------
 * 
 * Each PMT has a channel ID associated with a location; the relation between
 * the channel number value and the position of the PMT is fixed by a
 * convention (from LArSoft). This is a "logical" channel number.
 * On the other side, each PMT has a unique identifier that was assigned by
 * the ICARUS PMT Working Group and is used as key in an hardware database.
 * This is the PMT hardware channel number. These two channel numbers may be
 * (and in fact, are) different.
 * A similar distinction is made of the LVDS signals.
 * 
 * * Eight LVDS signals with discriminated PMT pair information are emitted by
 *   each CAEN V1730B PMT readout board. These bits are collected into four
 *   64-bit words, one for each of the four PMT walls. The word values are
 *   stored in `sbn::ExtraTriggerInfo`, which also defines the convention
 *   relating each of the bits in the 64-bit words with the physical location of 
 *   the PMT feeding the bit. This is equivalent to the logical PMT channel
 *   number described above, and here we call it a "logical LVDS bit ID".
 *   Analyzers will typically deal with this logical ID only. It is described by
 *   `icarus::trigger::PMTpairBitID` data type.
 * * The assembly of the bits in hardware includes the combination of 8-bit
 *   information from groups of three PMT readout boards, physically using LVDS
 *   format and for that often just dubbed "the LVDS signals", via custom LVDS
 *   boxes into 32-bit cables and connectors, with the excess 8 bits in each of
 *   them partly used for other information.
 *   Like in the case of single PMT, we define a "hardware LVDS bit ID"
 *   representing the position of a bit in the 32-bit word in the connector,
 *   and identifying the connector. The relation of these bits with the actual
 *   PMT channels is encoded in the hardware database; analyzers typically don't
 *   have to deal with them, because the trigger decoding software will
 *   translate them into the `sbn::ExtraTriggerInfo` convention (using, in fact,
 *   this `LVDSbitMaps` object for the mapping). The data type devoted to this
 *   concept is `LVDSHWbitID`. Note that the NIM FPGA boards combine two 32-bit
 *   connector words into 64-bit words for transportation. These 64-bit words
 *   are stored in the trigger DAQ fragment. At first look, these 64-bit words
 *   are similar  to the ones described in the previous bullet, but their bits
 *   are in a different relation with the actual readout boards and PMT.
 *   Also note that the decoding of these 64-bit words from the FPGA require,
 *   on top oft he information from this mapping, the information of where each
 *   connector is located within the word (hint: connectors 0 and 2 are stored
 *   as the most significant 32-bit word, 1 and 3 as the least significant one).
 *
 */
class icarus::trigger::LVDSbitMaps {
  
    public:
  
  // --- BEGIN --- Hard-coded geometry -----------------------------------------
  /// @name Good old hard-coded ICARUS geometry.
  /// @{
  
  /// Total number of PMT channels in the detector.
  static constexpr unsigned int NChannels = 360;
  
  /// Number of PMT channels in each cryostat.
  static constexpr unsigned int NChannelsPerCryostat = 180;
  
  /// Number of cryostats in the detector.
  static constexpr std::size_t NCryostats = NChannels / NChannelsPerCryostat;
  
  ///@}
  // --- END ----- Hard-coded geometry -----------------------------------------
  
  /// Error codes from this object. See `errorMessage()` code for descriptions.
  enum class Errors {
    NoError,               ///< No error condition.
    Error,                 ///< An unspecified error condition.
    NoDatabase,            ///< No channel mapping database access.
    NoPMTfragment,         ///< Missing a required PMT fragment information.
    NoPairBitInformation,  ///< No PMT pair bit information in the database.
    NoAdderBitInformation, ///< No adder bit information in the database.
    Unknown                ///< Unknown error.
  };
  
  /// Supported maps.
  enum class Map {
    PMTpairs,  ///< Map of PMT pairs ("LVDS").
    Adders,    ///< Map of adder signals.
  };
  
  /// Returns a predefined message describing the specified error `code`.
  static std::string errorMessage(Errors code);
  
  /// Type of exception thrown by this object.
  using Exception = cet::coded_exception<Errors, &errorMessage>;
  
  /// Default constructor: maps empty until `rebuild()` is called.
  LVDSbitMaps() = default;
  
  /// Constructor: builds all the maps based on the specified PMT channel map.
  LVDSbitMaps(icarusDB::IICARUSChannelMapProvider const& channelMap)
    { buildAllMaps(channelMap); }
  
  /// Constructor: builds the specified maps based on the specified PMT channel
  /// map.
  LVDSbitMaps(
    icarusDB::IICARUSChannelMapProvider const& channelMap,
    std::initializer_list<Map> maps
  )
    { buildMaps(channelMap, std::move(maps)); }
  
  
  // --- BEGIN -----  Query interface  -----------------------------------------
  /// @name Query interface.
  /// @{
  
  /// Returns whether the specified map is present.
  bool hasMap(Map map) const;
  
  /**
   * @brief Returns the source bit for the specified encoded PMT pair bit.
   * @param bit ID of the logical PMT pair bit to get the source information of
   * @return the number of hardware connector and bit where to find the value
   * @throws std::out_of_range if the specified bit is not in the map
   * 
   * The output `connector`/`bit` bits are sorted according to the prescription
   * of `sbn::ExtraTriggerInfo`, which requires the least significant bit to
   * have the most upstream PMT pair (i.e. the lowest channels).
   * 
   * The returned information includes:
   *  * the cryostat the bit comes from (the same as the one in `bit`);
   *  * the hardware connector number the bit comes from;
   *  * the hardware bit (wire) in the connector which carries the `bit` value;
   *  * the PMT channels covered by the bit.
   * 
   * See `buildPMTpairSourceMap()` for details of the mapping definitions.
   */
  LVDSinfo_t bitSource(PMTpairBitID const& bit) const;
  
  
  /**
   * @brief Returns the source bit for the specified encoded adder bit.
   * @param bit ID of the logical adder bit to get the source information of
   * @return the number of hardware connector and bit where to find the value
   * @throws std::out_of_range if the specified bit is not in the map
   * 
   * The output `connector`/`bit` bits are sorted according to the prescription
   * of `sbn::ExtraTriggerInfo`, which requires the least significant bit to
   * have the most upstream adder signal (i.e. the lowest channels).
   * 
   * The returned information includes:
   *  * the cryostat the bit comes from (the same as the one in `bit`);
   *  * the hardware connector number the bit comes from;
   *  * the hardware bit (wire) in the connector which carries the `bit` value;
   *  * the PMT channels covered by the bit.
   * 
   * See `buildAdderSourceMap()` for details of the mapping definitions.
   */
  AdderInfo_t bitSource(AdderBitID const& bit) const;
  
  
  /**
   * @brief Returns a list of channels paired in a single LVDS.
   * @return a list of vectors of one or two PMT channel ID
   * 
   * This method recreates the pairing of PMT channels into LVDS bits.
   * The order of the "pairs" follows the order of the logic bits.
   * These "pairs" can actually contain either two or one ID; empty pairs
   * (from invalid bits) are not included.
   */
  std::vector<std::vector<raw::Channel_t>> channelPairs() const;
  
  /// @}
  // --- END -------  Query interface  -----------------------------------------
  
  
  // --- BEGIN -----  Map updates  ---------------------------------------------
  
  /// Rebuilds all the LVDS maps based on the specified channel map status.
  void rebuild(icarusDB::IICARUSChannelMapProvider const& channelMap)
    { buildAllMaps(channelMap); }
  
  // --- END -------  Map updates  ---------------------------------------------
  
  
    private:
  
  /// Description of a LVDS bit and association to a PMT channel.
  /// @note This record does not include the logical LVDS bit ID it provides
  ///       information of.
  struct LVDSchannel_t {
    
    /// Special value identifying no channel set.
    static constexpr raw::Channel_t NoChannel
      = std::numeric_limits<raw::Channel_t>::max();
    
    raw::Channel_t channel = NoChannel; ///< ID of an associated PMT channel.
    
    LVDSHWbitID PMTpairSource; ///< Identifier of the PMT pair bit source.
    LVDSHWbitID adderSource; ///< Identifier of the adder bit source.
    
    /// Returns an "unique" index for the LVDS source.
    constexpr unsigned short int index() const { return PMTpairSource.index(); }
    
    /// Returns whether this information is fully valid (both channel and ID).
    operator bool() const
      { return (channel != NoChannel) && PMTpairSource && adderSource; }
    
    /// Order operator: sorted after the channel ID only.
    bool operator< (LVDSchannel_t const& other) const
      { return channel < other.channel; }
    
  }; // LVDSchannel_t
  
  
  /// Cached information of each PMT (sorted by PMT channel).
  std::vector<LVDSchannel_t> fPMTchannelInfo;
  
  /// For each PMT pair bit index, the HW LVDS ID and information are stored.
  std::vector<LVDSinfo_t> fLVDSinfoMap;
  
  /// For each adder bit ID, adder information is stored.
  std::vector<AdderInfo_t> fAdderInfoMap;
  
  
  /// Builds all the maps
  void buildAllMaps(icarusDB::IICARUSChannelMapProvider const& channelMap);
  
  /// Builds the specified maps. If any error occurs, an `Exception` is thrown.
  void buildMaps(
    icarusDB::IICARUSChannelMapProvider const& channelMap,
    std::initializer_list<Map> maps
    );
  
  /**
   * @brief Fills a list of PMT channels and their relevant information
   * @param channelMap service provider to associate each PMT channel to LVDS
   * @return the list, sorted by PMT channel ID
   */
  std::vector<LVDSchannel_t> buildChannelInfo
    (icarusDB::IICARUSChannelMapProvider const& channelMap) const;
  
  /**
   * @brief Builds the LVDS map connecting each bit to its hardware source.
   * @param channelMap service provider to associate each PMT channel to LVDS
   * @return a map
   * 
   * The map is associating the unique PMT pair bit ID (`PMTpairBitID::index()`)
   * and the information of the LVDS ID matching that ID.
   * The information includes the source in the hardware connector of the bit
   * and the involved channels (see `LVDSinfo_t` for more).
   */
  std::vector<LVDSinfo_t> buildPMTpairSourceMap
    (icarusDB::IICARUSChannelMapProvider const& channelMap) const;
  
  /**
   * @brief Builds the adder map connecting each bit to its hardware source.
   * @param channelMap service provider to associate each PMT channel to LVDS
   * @return a map
   * 
   * The map is associating the unique adder bit ID (`AdderBitID::index()`)
   * and the information of the hardware LVDS ID matching that ID.
   * The information includes the source in the hardware connector of the bit
   * and the involved channels (see `AdderInfo_t` for more).
   * 
   */
  std::vector<AdderInfo_t> buildAdderSourceMap
    (icarusDB::IICARUSChannelMapProvider const& channelMap) const;

  
}; // icarus::trigger::LVDSbitMaps


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
// ---  icarus::trigger::PMTwallBitID
// -----------------------------------------------------------------------------
template <std::size_t MaxBits>
constexpr auto icarus::trigger::details::PMTwallBitID<MaxBits>::fromIndex
  (Index_t index) -> BitID_t
{
  auto const [ cryostat, wb ]
    = std::div(index, NPMTwallsPerCryostat * NStatusWordBits);
  auto const [ wall, bit ] = std::div(wb, NStatusWordBits);
  return { (std::size_t) cryostat, (std::size_t) wall, (StatusBit_t) bit };
} // icarus::trigger::PMTwallBitID<>::fromIndex()


// -----------------------------------------------------------------------------
template <std::size_t MaxBits>
std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, PMTwallBitID<MaxBits> const& bit)
{
  static constexpr char Side[] = "EW";
  auto const fillch = out.fill(); // width is reset to 0, fill character is not
  out << (bit.hasCryostat()? Side[bit.cryostat]: '?')
    << '/' << (bit.hasPMTwall()? Side[bit.PMTwall]: '?')
    << '/' << std::setw(2) << std::setfill('0') << bit.statusBit
    << std::setfill(fillch);
  return out;
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_LVDSBITMAPS_H
