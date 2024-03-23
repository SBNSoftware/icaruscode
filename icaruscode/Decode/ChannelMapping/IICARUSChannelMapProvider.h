/**
 * @file   icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h
 * @brief  Interface class for hardware/software channel mapping for ICARUS.
 * @author T. Usher
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAPPROVIDER_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAPPROVIDER_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapDataTypes.h"
#include "icaruscode/Utilities/CacheCounter.h"

// C/C++ standard libraries
#include <vector>
#include <map>
#include <string>
#include <tuple>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------

namespace icarusDB { class IICARUSChannelMapProvider; }
/**
 * @brief Interface of ICARUS channel mapping service provider.
 * 
 * The service providers implementing this interface provide access to ICARUS
 * channel mapping.
 * 
 * The mapping organization reflects what was needed at the time of inception.
 * 
 * The service provider is expected to have status and to be compatible with
 * multithreading only within the same run period. In fact, a current run
 * (`forRun()`) or run period (`forPeriod()`) must be selected and all the
 * provided information will pertain the selected run period.
 * 
 * Also, given that the interface yields complex structures by reference, it is
 * expected that implementations cache results in their final form and provide
 * direct access to that cache.
 * 
 * The interface includes cache tracking support via `util::CacheCounter`.
 * Implementations are expected to maintain the tracking of the cache with the
 * default (empty) name; by default the cache will be always considered
 * outdated.
 */
class icarusDB::IICARUSChannelMapProvider: public util::CacheCounter {
    
    public:
  
  virtual ~IICARUSChannelMapProvider() noexcept = default;
  
  
  /// --- BEGIN --- Data period selection --------------------------------------
  /// @name Data period selection
  /// @{
  
  /// Loads the mapping for `run`, returns whether a new mapping was loaded.
  virtual bool forRun(int run) = 0;
  
  /// Loads the mapping for `period`, returns whether a new mapping was loaded.
  virtual bool forPeriod(icarusDB::RunPeriod period) = 0;
  
  /// @}
  /// --- END ----- Data period selection --------------------------------------
  
  
  /// --- BEGIN --- TPC information --------------------------------------------
  /// @name TPC information
  /// @{
  
  /// --- BEGIN - - TPC fragment information - - - - - - - - - - - - - - - - - -
  /// @name TPC fragment information
  /// @{
  
  /// Returns whether the specified `ID` is a known TPC fragment ID.
  virtual bool hasFragmentID(unsigned int ID) const = 0;

  /// Returns the number of TPC fragment IDs known to the service.
  virtual unsigned int nTPCfragmentIDs() const = 0;
  
  /// Returns the name of the crate served by the specified `fragmentID`.
  virtual std::string const& getCrateName(unsigned int fragmentID) const = 0;
  
  /// Returns the list of board IDs included in the specified `fragmentID`.
  virtual ReadoutIDVec const& getReadoutBoardVec(unsigned int fragmentID) const
    = 0;
  
  /// @}
  /// --- END - - - TPC fragment information - - - - - - - - - - - - - - - - - -
  
  /**
   * @brief Returns the full TPC channel mapping.
   * 
   * The returned mapping is an associative container, associating each board
   * ID (`boardID`) to the following information (in a `std::tuple`):
   *  * `[0]` slot number the board is on (like `getBoardSlot(boardID)`)
   *  * `[1]` list of channels served by the board, and their plane
   *          (like `getChannelPlanePair(boardID)`)
   * 
   */
  virtual TPCReadoutBoardToChannelMap const& getReadoutBoardToChannelMap() const
    = 0;
  
  /// --- BEGIN - - TPC board information  - - - - - - - - - - - - - - - - - - -
  /// @name TPC board information
  /// @{
  
  /// Returns whether there is a board with the specified `ID`.
  virtual bool hasBoardID(unsigned int ID) const = 0;
  
  /// Returns the number of TPC board IDs known to the service.
  virtual unsigned int nTPCboardIDs() const = 0;
  
  /// Returns the number of slot the `boardID` is on.
  virtual unsigned int getBoardSlot(unsigned int boardID) const = 0;
  
  /// Returns a list of channels served by the `boardID` and for each the plane
  /// it is on (`0` to `2`).
  virtual ChannelPlanePairVec const& getChannelPlanePair
    (unsigned int boardID) const = 0;

  /// @}
  /// --- END - - - TPC board information  - - - - - - - - - - - - - - - - - - -
  
  
  /// --- BEGIN --- PMT information --------------------------------------------
  /// @name PMT information
  /// @{
  
  /// Returns whether the specified fragment `ID` is known to the mapping.
  virtual bool hasPMTDigitizerID(unsigned int ID) const = 0;
  
  /// Returns the number of PMT fragment IDs known to the mapping.
  virtual unsigned int nPMTfragmentIDs() const = 0;
  
  /// Returns records on all the PMT channels covered by the fragment `ID`.
  virtual PMTdigitizerInfoVec const& getPMTchannelInfo(unsigned int ID) const
    = 0;
  
  /// @}
  /// --- END ----- PMT information --------------------------------------------


  /// --- BEGIN --- CRT information --------------------------------------------
  /// @name CRT information
  /// @{
  
  /// Returns the Sim Mac address corresponding to the specified side CRT
  /// hardware address.
  virtual unsigned int getSimMacAddress(unsigned int hwmacaddress) const
    = 0;
  
  /// Returns the sim Mac address corresponding to the specified top CRT
  /// hardware address.
  virtual unsigned int gettopSimMacAddress(unsigned int) const = 0;

  /// Returns the gain and pedestal for Side CRT.
  virtual std::pair<double, double> getSideCRTCalibrationMap
    (int mac5, int chan) const = 0;
  
  /// @}
  /// --- END ----- CRT information --------------------------------------------

}; // icarusDB::IICARUSChannelMapProvider


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAPPROVIDER_H
