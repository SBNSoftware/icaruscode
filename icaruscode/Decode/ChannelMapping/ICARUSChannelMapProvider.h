////////////////////////////////////////////////////////////////////////
/// \file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h
/// \author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
/// \see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.cxx
////////////////////////////////////////////////////////////////////////

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPROVIDER_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPROVIDER_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <string>
#include <memory> // std::unique_ptr<>


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMapProvider; }
/**
 * @brief Interface to the channel mapping databases of TPC, CRT and PMT subdetectors.
 * 
 * Before retrieving the information, a run period (`forPeriod()`) or a run
 * (`forRun()`) must be selected.
 */
class icarusDB::ICARUSChannelMapProvider: public IICARUSChannelMap
{
public:
    
    // Constructor, destructor.
    ICARUSChannelMapProvider(const fhicl::ParameterSet& pset);
    
    // Section to access fragment to board mapping
    bool                                    hasFragmentID(const unsigned int)       const override;

    /// Returns the number of TPC fragment IDs known to the service.
    unsigned int                            nTPCfragmentIDs() const override;
    const std::string&                      getCrateName(const unsigned int)        const override;
    const ReadoutIDVec&                     getReadoutBoardVec(const unsigned int)  const override;
    const TPCReadoutBoardToChannelMap&      getReadoutBoardToChannelMap()           const override;

    // Section to access channel information for a given board
    bool                                    hasBoardID(const unsigned int)          const override;

    /// Returns the number of TPC board IDs known to the service.
    unsigned int                            nTPCboardIDs() const override;
    unsigned int                            getBoardSlot(const unsigned int)        const override;
    const ChannelPlanePairVec&              getChannelPlanePair(const unsigned int) const override;

    // Section for PMT channel mapping
    bool                                    hasPMTDigitizerID(const unsigned int)   const override;

    /// Returns the number of PMT fragment IDs known to the service.
    unsigned int                            nPMTfragmentIDs() const override;
    const DigitizerChannelChannelIDPairVec& getChannelIDPairVec(const unsigned int) const override;

    // Section for CRT channel mapping    
    unsigned int                            getSimMacAddress   (const unsigned int)    const override;
    unsigned int                            gettopSimMacAddress(const unsigned int)    const override;

    /// Returns the Gain and Pedestal for Side CRT 
    std::pair<double, double>               getSideCRTCalibrationMap(int mac5, int chan) const override;    

    /// Returns the channel mapping database key for the specified PMT fragment ID.
    static constexpr unsigned int PMTfragmentIDtoDBkey(unsigned int fragmentID);
    
    /// Returns the PMT fragment ID for the specified channel mapping database key.
    static constexpr unsigned int DBkeyToPMTfragmentID(unsigned int DBkey);

  
    /// Loads the mapping for `run`, returns whether a new mapping was loaded.
    bool                                    forRun(int run)                         override;
    
    /// Loads the mapping for `period`, returns whether a new mapping was loaded.
    bool                                    forPeriod(icarusDB::RunPeriod period)   override;

private:
    
    bool fDiagnosticOutput;
      
    IChannelMapping::TPCFragmentIDToReadoutIDMap   fFragmentToReadoutMap;
      
    IChannelMapping::TPCReadoutBoardToChannelMap   fReadoutBoardToChannelMap;

    IChannelMapping::FragmentToDigitizerChannelMap fFragmentToDigitizerMap; 

    IChannelMapping::CRTChannelIDToHWtoSimMacAddressPairMap fCRTChannelIDToHWtoSimMacAddressPairMap;

    IChannelMapping::TopCRTHWtoSimMacAddressPairMap fTopCRTHWtoSimMacAddressPairMap;

    IChannelMapping::SideCRTChannelToCalibrationMap fSideCRTChannelToCalibrationMap;

    std::unique_ptr<IChannelMapping>               fChannelMappingTool;


    /// Has the channel mapping tool fill the mapping caches.
    void readFromDatabase();

    
    /// Returns the list of board channel-to-PMT channel ID mapping within the specified fragment.
    /// @returns a pointer to the mapping list, or `nullptr` if invalid fragment
    DigitizerChannelChannelIDPairVec const* findPMTfragmentEntry
      (unsigned int fragmentID) const;
    
    

}; // icarusDB::ICARUSChannelMapProvider


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPROVIDER_H
