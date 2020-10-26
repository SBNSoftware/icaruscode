////////////////////////////////////////////////////////////////////////
/// \file   ICARUSChannelMap_service.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()

#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/TPCChannelmapping.h"

#include <fstream>

namespace icarusDB
{

class ICARUSChannelMap : virtual public IICARUSChannelMap
{
public:
    
    // Constructor, destructor.
    ICARUSChannelMap(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);
    
    // Section to access fragment to board mapping
    bool                                    hasFragmentID(const unsigned int)       const override;
    const std::string&                      getCrateName(const unsigned int)        const override;
    const ReadoutIDVec&                     getReadoutBoardVec(const unsigned int)  const override;

    // Section to access channel information for a given board
    bool                                    hasBoardID(const unsigned int)          const override;
    unsigned int                            getBoardSlot(const unsigned int)        const override;
    const ChannelPlanePairVec&              getChannelPlanePair(const unsigned int) const override;

    // Section for PMT channel mapping
    bool                                    hasPMTDigitizerID(const unsigned int)   const override;
    const DigitizerChannelChannelIDPairVec& getChannelIDPairVec(const unsigned int) const override;
    
private:
    
    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    bool fDiagnosticOutput;
      
    database::TPCFragmentIDToReadoutIDMap   fFragmentToReadoutMap;
      
    database::TPCReadoutBoardToChannelMap   fReadoutBoardToChannelMap;

    database::FragmentToDigitizerChannelMap fFragmentToDigitizerMap; 

};

//----------------------------------------------------------------------
// Constructor.
ICARUSChannelMap::ICARUSChannelMap(const fhicl::ParameterSet& pset, art::ActivityRegistry& /* reg */)
{
    reconfigure(pset);
}

//----------------------------------------------------------------------
// Reconfigure method.
void ICARUSChannelMap::reconfigure(const fhicl::ParameterSet& pset)
{
    mf::LogInfo("ICARUSChannelMap") << "Building the channel mapping" ;

    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);
    
    cet::cpu_timer theClockFragmentIDs;

    theClockFragmentIDs.start();

    if (database::BuildTPCFragmentIDToReadoutIDMap(fFragmentToReadoutMap))
    {
        throw cet::exception("ICARUSChannelMap") << "Cannot recover the Fragment ID channel map from the database \n";
    }
    else if (fDiagnosticOutput)
    {
        std::cout << "FragmentID to Readout ID map has " << fFragmentToReadoutMap.size() << " elements";
        for(const auto& pair : fFragmentToReadoutMap) std::cout << "   Frag: " << std::hex << pair.first << ", Crate: " << pair.second.first << ", # boards: " << std::dec << pair.second.second.size() << std::endl;
    }

    theClockFragmentIDs.stop();

    double fragmentIDsTime = theClockFragmentIDs.accumulated_real_time();

    cet::cpu_timer theClockReadoutIDs;

    theClockReadoutIDs.start();

    if (database::BuildTPCReadoutBoardToChannelMap(fReadoutBoardToChannelMap))
    {
        std::cout << "******* FAILED TO CONFIGURE CHANNEL MAP ********" << std::endl;
        throw cet::exception("ICARUSChannelMap") << "POS didn't read the F'ing database again \n";
    }

    // Do the channel mapping initialization
    if (database::BuildFragmentToDigitizerChannelMap(fFragmentToDigitizerMap))
    {
        throw cet::exception("PMTDecoder") << "Cannot recover the Fragment ID channel map from the database \n";
    }
    else if (fDiagnosticOutput)
    {
        std::cout << "FragmentID to Readout ID map has " << fFragmentToDigitizerMap.size() << " Fragment IDs";
        for(const auto& pair : fFragmentToDigitizerMap) std::cout << "   Frag: " << std::hex << pair.first << ", # pairs: " << std::dec << pair.second.size() << std::endl;
    }

    theClockReadoutIDs.stop();

    double readoutIDsTime = theClockReadoutIDs.accumulated_real_time();


    mf::LogInfo("ICARUSChannelMap") << "==> FragmentID map time: " << fragmentIDsTime << ", Readout IDs time: " << readoutIDsTime << std::endl;
    
    return;
}

bool ICARUSChannelMap::hasFragmentID(const unsigned int fragmentID) const 
{
    return fFragmentToReadoutMap.find(fragmentID) != fFragmentToReadoutMap.end();
}

const std::string&  ICARUSChannelMap::getCrateName(const unsigned int fragmentID) const
{
    database::TPCFragmentIDToReadoutIDMap::const_iterator fragToReadoutItr = fFragmentToReadoutMap.find(fragmentID);

    if (fragToReadoutItr == fFragmentToReadoutMap.end())
        throw cet::exception("ICARUSChannelMap") << "Fragment ID " << fragmentID << " not found in lookup map when looking up crate name \n";

    return fragToReadoutItr->second.first;
}

const ReadoutIDVec& ICARUSChannelMap::getReadoutBoardVec(const unsigned int fragmentID) const
{
    database::TPCFragmentIDToReadoutIDMap::const_iterator fragToReadoutItr = fFragmentToReadoutMap.find(fragmentID);

    if (fragToReadoutItr == fFragmentToReadoutMap.end())
        throw cet::exception("ICARUSChannelMap") << "Fragment ID " << fragmentID << " not found in lookup map when looking up board vector \n";

    return fragToReadoutItr->second.second;

}

bool ICARUSChannelMap::hasBoardID(const unsigned int boardID)  const
{
    return fReadoutBoardToChannelMap.find(boardID) != fReadoutBoardToChannelMap.end();
}

unsigned int ICARUSChannelMap::getBoardSlot(const unsigned int boardID)  const
{
    database::TPCReadoutBoardToChannelMap::const_iterator readoutBoardItr = fReadoutBoardToChannelMap.find(boardID);

    if (readoutBoardItr == fReadoutBoardToChannelMap.end())
        throw cet::exception("ICARUSChannelMap") << "Board ID " << boardID << " not found in lookup map when looking up board slot \n";

    return readoutBoardItr->second.first;
}

 const ChannelPlanePairVec& ICARUSChannelMap::getChannelPlanePair(const unsigned int boardID) const
{
    database::TPCReadoutBoardToChannelMap::const_iterator readoutBoardItr = fReadoutBoardToChannelMap.find(boardID);

    if (readoutBoardItr == fReadoutBoardToChannelMap.end())
        throw cet::exception("ICARUSChannelMap") << "Board ID " << boardID << " not found in lookup map when looking up channel/plane pair \n";

    return readoutBoardItr->second.second;

}

bool ICARUSChannelMap::hasPMTDigitizerID(const unsigned int fragmentID)   const
{
    return fFragmentToDigitizerMap.find(fragmentID) != fFragmentToDigitizerMap.end();
}

const DigitizerChannelChannelIDPairVec& ICARUSChannelMap::getChannelIDPairVec(const unsigned int fragmentID) const
{
    database::FragmentToDigitizerChannelMap::const_iterator digitizerItr = fFragmentToDigitizerMap.find(fragmentID);

    if (digitizerItr == fFragmentToDigitizerMap.end())
        throw cet::exception("ICARUSChannelMap") << "Fragment ID " << fragmentID << " not found in lookup map when looking for PMT channel info \n";

    return digitizerItr->second;
      
}


} // end namespace

DECLARE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap, SHARED)
DEFINE_ART_SERVICE_INTERFACE_IMPL(icarusDB::ICARUSChannelMap, icarusDB::IICARUSChannelMap)

