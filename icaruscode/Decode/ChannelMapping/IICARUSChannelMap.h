/**
 * @file   icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h
 * @brief  Interface class for hardware/software channel mapping for ICARUS.
 * @author T. Usher
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAP_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAP_H

#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include <vector>
#include <map>
#include <string>
#include <tuple>

namespace icarusDB {

using ReadoutIDVec                     = std::vector<unsigned int>;
using ChannelPlanePair                 = std::pair<unsigned int, unsigned int>;
using ChannelPlanePairVec              = std::vector<ChannelPlanePair>;

using DigitizerChannelChannelIDPair    = std::tuple<size_t,size_t, size_t>;
using DigitizerChannelChannelIDPairVec = std::vector<DigitizerChannelChannelIDPair>;

using ChannelPlanePair                 = std::pair<unsigned int, unsigned int>;
using ChannelPlanePairVec              = std::vector<ChannelPlanePair>;
using SlotChannelVecPair               = std::pair<unsigned int, ChannelPlanePairVec>;
using TPCReadoutBoardToChannelMap      = std::map<unsigned int, SlotChannelVecPair>;

class IICARUSChannelMap //: private lar::EnsureOnlyOneSchedule
{
public:
    virtual ~IICARUSChannelMap() noexcept = default;
    
    
    /// --- BEGIN --- Data period selection ------------------------------------
    /// @name Data period selection
    /// @{
    
    /// Loads the mapping for `run`, returns whether a new mapping was loaded.
    virtual bool                                    forRun(int run)                               = 0;
    
    /// Loads the mapping for `period`, returns whether a new mapping was loaded.
    virtual bool                                    forPeriod(icarusDB::RunPeriod period)         = 0;
    
    /// @}
    /// --- END ----- Data period selection ------------------------------------
    
    
    /// --- BEGIN --- TPC information ------------------------------------------
    /// @name TPC information
    /// @{
    
    /// --- BEGIN - - TPC fragment information - - - - - - - - - - - - - - - - -
    /// @name TPC fragment information
    /// @{
    
    /// Returns whether the specified `ID` is a known TPC fragment ID.
    virtual bool                                    hasFragmentID(const unsigned int ID)       const = 0;
    /// Returns the number of TPC fragment IDs known to the service.

    /// Returns the number of known TPC fragments.
    virtual unsigned int                            nTPCfragmentIDs() const = 0;
    
    /// Returns the name of the crate served by the specified `fragmentID`.
    virtual const std::string&                      getCrateName(const unsigned int fragmentID) const = 0;
    
    /// Returns the list of board IDs included in the specified `fragmentID`.
    virtual const ReadoutIDVec&                     getReadoutBoardVec(const unsigned int fragmentID) const = 0;
    
    /// @}
    /// --- END - - - TPC fragment information - - - - - - - - - - - - - - - - -
    
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
    virtual const TPCReadoutBoardToChannelMap&      getReadoutBoardToChannelMap()           const = 0;
    
    /// --- BEGIN - - TPC board information  - - - - - - - - - - - - - - - - - -
    /// @name TPC board information
    /// @{
    
    /// Returns whether there is a board with the specified `ID`.
    virtual bool                                    hasBoardID(const unsigned int ID)          const = 0;
    /// Returns the number of TPC board IDs known to the service.
    virtual unsigned int                            nTPCboardIDs() const = 0;
    /// Returns the number of slot the `boardID` is on.
    virtual unsigned int                            getBoardSlot(const unsigned int boardID) const = 0;
    /// Returns a list of channels served by the `boardID` and for each the plane it is on (`0` to `2`).
    virtual const ChannelPlanePairVec&              getChannelPlanePair(const unsigned int boardID) const = 0;

    /// @}
    /// --- END - - - TPC board information  - - - - - - - - - - - - - - - - - -
    
    
    /// --- BEGIN --- PMT information ------------------------------------------
    /// @name PMT information
    /// @{
    
    /// Returns whether the specified fragment `ID` is known to the mapping.
    virtual bool                                    hasPMTDigitizerID(const unsigned int ID) const = 0;
    
    /// Returns the number of PMT fragment IDs known to the mapping.
    virtual unsigned int                            nPMTfragmentIDs() const = 0;
    
    /// Returns a list of triplets: 
    virtual const DigitizerChannelChannelIDPairVec& getChannelIDPairVec(const unsigned int) const = 0;
    
    /// @}
    /// --- END ----- PMT information ------------------------------------------


    /// --- BEGIN --- CRT information ------------------------------------------
    /// @name CRT information
    /// @{
    
    /// Returns the sim Mac address corresponding to the specified side CRT hardware address.
    virtual unsigned int                            getSimMacAddress(const unsigned int hwmacaddress)    const = 0;
    /// Returns the sim Mac address corresponding to the specified top CRT hardware address.
    virtual unsigned int                         gettopSimMacAddress(const unsigned int)    const = 0;

    /// Returns the Gain and Pedestal for Side CRT.
    virtual std::pair<double, double>          getSideCRTCalibrationMap(int mac5, int chan) const = 0;
    
    /// @}
    /// --- END ----- CRT information ------------------------------------------

}; // IICARUSChannelMap

} // namespace icarusDB

DECLARE_ART_SERVICE_INTERFACE(icarusDB::IICARUSChannelMap, SHARED)


#endif // ICARUSCODE_DECODE_CHANNELMAPPING_IICARUSCHANNELMAP_H
