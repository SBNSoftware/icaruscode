////////////////////////////////////////////////////////////////////////
/// \file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.cxx
/// \author T. Usher (factorized by Gianluca Petrillo, petrillo@slac.stanford.edu)
/// \see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h
////////////////////////////////////////////////////////////////////////

// library header
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h"

#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"

#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
//#include "icaruscode/Decode/ChannelMapping/TPCChannelmapping.h"
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"

#include <string>
#include <iostream>

namespace icarusDB
{


//----------------------------------------------------------------------
// Constructor.
ICARUSChannelMapProvider::ICARUSChannelMapProvider(const fhicl::ParameterSet& pset) {

    mf::LogInfo("ICARUSChannelMapProvider") << "Building the channel mapping" ;

    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet&channelMappingParams = pset.get<fhicl::ParameterSet>("ChannelMappingTool");

    // Get instance of the mapping tool (allowing switch between database instances)
    fChannelMappingTool = art::make_tool<IChannelMapping>(channelMappingParams);

    cet::cpu_timer theClockFragmentIDs;

    theClockFragmentIDs.start();

    if (fChannelMappingTool->BuildTPCFragmentIDToReadoutIDMap(fFragmentToReadoutMap))
    {
        throw cet::exception("ICARUSChannelMapProvider") << "Cannot recover the Fragment ID channel map from the database \n";
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

    if (fChannelMappingTool->BuildTPCReadoutBoardToChannelMap(fReadoutBoardToChannelMap))
    {
        std::cout << "******* FAILED TO CONFIGURE CHANNEL MAP ********" << std::endl;
        throw cet::exception("ICARUSChannelMapProvider") << "POS didn't read the F'ing database again \n";
    }

    // Do the channel mapping initialization
    if (fChannelMappingTool->BuildFragmentToDigitizerChannelMap(fFragmentToDigitizerMap))
      {
	throw cet::exception("ICARUSChannelMapProvider") << "Cannot recover the Fragment ID channel map from the database \n";
      }
    else if (fDiagnosticOutput)
      {
	std::cout << "FragmentID to Readout ID map has " << fFragmentToDigitizerMap.size() << " Fragment IDs";
        for(const auto& pair : fFragmentToDigitizerMap) std::cout << "   Frag: " << std::hex << pair.first << ", # pairs: " << std::dec << pair.second.size() << std::endl;
      }
    
    // Do the channel mapping initialization for CRT
    if (fChannelMappingTool->BuildCRTChannelIDToHWtoSimMacAddressPairMap(fCRTChannelIDToHWtoSimMacAddressPairMap))
      {
        throw cet::exception("CRTDecoder") << "Cannot recover the HW MAC Address  from the database \n";
      }
    else if (fDiagnosticOutput)
      {
	std::cout << "ChannelID to MacAddress map has " << fCRTChannelIDToHWtoSimMacAddressPairMap.size() << " Channel IDs";
        for(const auto& pair : fCRTChannelIDToHWtoSimMacAddressPairMap) std::cout <<" ChannelID: "<< pair.first
                                                                                  << ", hw mac address: " << pair.second.first
                                                                                  <<", sim mac address: " << pair.second.second << std::endl;
      }
    
    
    theClockReadoutIDs.stop();

    double readoutIDsTime = theClockReadoutIDs.accumulated_real_time();


    mf::LogInfo("ICARUSChannelMapProvider") << "==> FragmentID map time: " << fragmentIDsTime << ", Readout IDs time: " << readoutIDsTime << std::endl;
    
    return;
}

bool ICARUSChannelMapProvider::hasFragmentID(const unsigned int fragmentID) const 
{
    return fFragmentToReadoutMap.find(fragmentID) != fFragmentToReadoutMap.end();
}


unsigned int ICARUSChannelMapProvider::nTPCfragmentIDs() const {
  return fFragmentToReadoutMap.size();
}


const std::string&  ICARUSChannelMapProvider::getCrateName(const unsigned int fragmentID) const
{
    IChannelMapping::TPCFragmentIDToReadoutIDMap::const_iterator fragToReadoutItr = fFragmentToReadoutMap.find(fragmentID);

    if (fragToReadoutItr == fFragmentToReadoutMap.end())
        throw cet::exception("ICARUSChannelMapProvider") << "Fragment ID " << fragmentID << " not found in lookup map when looking up crate name \n";

    return fragToReadoutItr->second.first;
}

const ReadoutIDVec& ICARUSChannelMapProvider::getReadoutBoardVec(const unsigned int fragmentID) const
{
    IChannelMapping::TPCFragmentIDToReadoutIDMap::const_iterator fragToReadoutItr = fFragmentToReadoutMap.find(fragmentID);

    if (fragToReadoutItr == fFragmentToReadoutMap.end())
        throw cet::exception("ICARUSChannelMapProvider") << "Fragment ID " << fragmentID << " not found in lookup map when looking up board vector \n";

    return fragToReadoutItr->second.second;

}

bool ICARUSChannelMapProvider::hasBoardID(const unsigned int boardID)  const
{
    return fReadoutBoardToChannelMap.find(boardID) != fReadoutBoardToChannelMap.end();
}


unsigned int ICARUSChannelMapProvider::nTPCboardIDs() const {
  return fReadoutBoardToChannelMap.size();
}


unsigned int ICARUSChannelMapProvider::getBoardSlot(const unsigned int boardID)  const
{
    IChannelMapping::TPCReadoutBoardToChannelMap::const_iterator readoutBoardItr = fReadoutBoardToChannelMap.find(boardID);

    if (readoutBoardItr == fReadoutBoardToChannelMap.end())
        throw cet::exception("ICARUSChannelMapProvider") << "Board ID " << boardID << " not found in lookup map when looking up board slot \n";

    return readoutBoardItr->second.first;
}

 const ChannelPlanePairVec& ICARUSChannelMapProvider::getChannelPlanePair(const unsigned int boardID) const
{
    IChannelMapping::TPCReadoutBoardToChannelMap::const_iterator readoutBoardItr = fReadoutBoardToChannelMap.find(boardID);

    if (readoutBoardItr == fReadoutBoardToChannelMap.end())
        throw cet::exception("ICARUSChannelMapProvider") << "Board ID " << boardID << " not found in lookup map when looking up channel/plane pair \n";

    return readoutBoardItr->second.second;

}

bool ICARUSChannelMapProvider::hasPMTDigitizerID(const unsigned int fragmentID)   const
{
    return fFragmentToDigitizerMap.find(fragmentID) != fFragmentToDigitizerMap.end();
}


unsigned int ICARUSChannelMapProvider::nPMTfragmentIDs() const {
  return fFragmentToDigitizerMap.size();
}


const DigitizerChannelChannelIDPairVec& ICARUSChannelMapProvider::getChannelIDPairVec(const unsigned int fragmentID) const
{
    IChannelMapping::FragmentToDigitizerChannelMap::const_iterator digitizerItr = fFragmentToDigitizerMap.find(fragmentID);

    if (digitizerItr == fFragmentToDigitizerMap.end())
        throw cet::exception("ICARUSChannelMapProvider") << "Fragment ID " << fragmentID << " not found in lookup map when looking for PMT channel info \n";

    return digitizerItr->second;
      
}

  unsigned int ICARUSChannelMapProvider::getSimMacAddress(const unsigned int hwmacaddress)  const
  {
    unsigned int   simmacaddress = -99999;
    for(const auto& pair : fCRTChannelIDToHWtoSimMacAddressPairMap){
      if (pair.second.first == hwmacaddress)
        simmacaddress = pair.second.second;
    }
    return simmacaddress;
  }
  
} // end namespace

