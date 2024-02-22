/**
 *  @file   IChannelMapping.h
 *
 *  @brief  This provides an art tool interface definition for tools handle the channel mapping
 *          The idea is to be able to switch between postgres and sqlite implementations
 *
 *  @author usher@slac.stanford.edu
 *
 */
#ifndef IChannelMapping_h
#define IChannelMapping_h

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// C/C++ standard libraries
#include <algorithm> // std::transform()
#include <array>
#include <vector>
#include <map>
#include <string>
#include <tuple>
#include <cctype> // std::tolower()


//------------------------------------------------------------------------------------------------------------------------------------------

namespace icarusDB
{
  /**
   *  @brief  IChannelMapping interface class definiton
   */
  class IChannelMapping
  {
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
    
    /**
     *  @brief Define the returned data structures for a mapping between TPC Fragment IDs
     *         and the related crate and readout information. 
     *         Then define the function interface to fill these data structures 
     */
    using ReadoutIDVec                = std::vector<unsigned int>;
    using CrateNameReadoutIDPair      = std::pair<std::string,ReadoutIDVec>;
    using TPCFragmentIDToReadoutIDMap = std::map<unsigned int, CrateNameReadoutIDPair>;
    
    virtual int BuildTPCFragmentIDToReadoutIDMap(TPCFragmentIDToReadoutIDMap&) const = 0;
    
    /**
     *  @brief Define the returned data structures for a mapping between TPC readout boards
     *         and the channel information 
     *         Then define the function interface to fill these data structures 
     */
    using ChannelPlanePair            = std::pair<unsigned int, unsigned int>;
    using ChannelPlanePairVec         = std::vector<ChannelPlanePair>;
    using SlotChannelVecPair          = std::pair<unsigned int, ChannelPlanePairVec>;
    using TPCReadoutBoardToChannelMap = std::map<unsigned int, SlotChannelVecPair>;

    virtual int BuildTPCReadoutBoardToChannelMap(TPCReadoutBoardToChannelMap&) const = 0;

    /**
     *  @brief Define the returned data structures for a mapping between PMT Fragment IDs
     *         and the related crate and readout information. 
     *         Then define the function interface to fill these data structures 
     */
    using DigitizerChannelChannelIDPair    = std::tuple<size_t,size_t,size_t>; // std::tuple<DigitizerChannel, ChannelID, LaserChannel>
    using DigitizerChannelChannelIDPairVec = std::vector<DigitizerChannelChannelIDPair>;
    using FragmentToDigitizerChannelMap    = std::map<size_t, DigitizerChannelChannelIDPairVec>;

    virtual int BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap&) const = 0;

    /**
     *  @brief Define the returned data structures for a mapping between CRT hardware mac_address
     *         to the simulated mac_address. 
     *         Then define the function interface to fill these data structures
     */

    using CRTHWtoSimMacAddressPair                  = std::pair<unsigned int, unsigned int>;
    using CRTChannelIDToHWtoSimMacAddressPairMap    = std::map<unsigned int, CRTHWtoSimMacAddressPair>;

    virtual int BuildCRTChannelIDToHWtoSimMacAddressPairMap(CRTChannelIDToHWtoSimMacAddressPairMap&) const = 0;


    using TopCRTHWtoSimMacAddressPairMap    = std::map<unsigned int, unsigned int>;

    virtual int BuildTopCRTHWtoSimMacAddressPairMap(TopCRTHWtoSimMacAddressPairMap&) const = 0;

    
    using SideCRTMac5ToChannelPair = std::pair<unsigned int, unsigned int>;
    using SideCRTGainToPedPair = std::pair<double, double>;
    using SideCRTChannelToCalibrationMap = std::map< SideCRTMac5ToChannelPair, SideCRTGainToPedPair >;

    virtual int BuildSideCRTCalibrationMap(SideCRTChannelToCalibrationMap&) const = 0;

    // --- BEGIN --- Utility functions for the implementation ------------------
    
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
    
    // --- BEGIN --- Utility functions for the implementation ------------------
    
};

} // namespace lar_cluster3d
#endif
