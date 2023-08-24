////////////////////////////////////////////////////////////////////////
/// \file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.cxx
/// \author T. Usher (factorized by Gianluca Petrillo, petrillo@slac.stanford.edu)
/// \see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h
////////////////////////////////////////////////////////////////////////

// library header
#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProvider.h"

#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"

#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"

#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"

#include <string>
#include <iostream>
#include <cassert>

namespace icarusDB
{


//----------------------------------------------------------------------
// Constructor.
ICARUSChannelMapProvider::ICARUSChannelMapProvider(const fhicl::ParameterSet& pset) {

    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet&channelMappingParams = pset.get<fhicl::ParameterSet>("ChannelMappingTool");

    // Get instance of the mapping tool (allowing switch between database instances)
    fChannelMappingTool = art::make_tool<IChannelMapping>(channelMappingParams);

}

bool ICARUSChannelMapProvider::forRun(int run) {
  
  RunPeriod period = RunPeriod::NPeriods;
  try {
    period = RunPeriods::withRun(run);
  } catch (...) {
    mf::LogError("ICARUSChannelMapProvider")
      << "Failed to associate a run period number to run " << run;
    throw;
  }
  
  if (period == RunPeriod::NPeriods) {
    throw std::runtime_error(
      "ICARUSChannelMapProvider::forRun(): can't determine the period of run "
      + std::to_string(run) + "!\n"
      );
  }
  
  return forPeriod(period);
  
}


bool ICARUSChannelMapProvider::forPeriod(icarusDB::RunPeriod period) {
  
  // if cache is not invalidated, we don't refresh it
  if (!fChannelMappingTool->SelectPeriod(period)) return false;
  
  readFromDatabase();
  return true;
}


void ICARUSChannelMapProvider::readFromDatabase() {

    mf::LogInfo("ICARUSChannelMapProvider") << "Building the channel mapping" ;

    cet::cpu_timer theClockFragmentIDs;

    theClockFragmentIDs.start();

    fFragmentToReadoutMap.clear();
    if (fChannelMappingTool->BuildTPCFragmentIDToReadoutIDMap(fFragmentToReadoutMap))
    {
        throw cet::exception("ICARUSChannelMapProvider") << "Cannot recover the Fragment ID channel map from the database \n";
    }
    else if (fDiagnosticOutput)
    {
        std::cout << "FragmentID to Readout ID map has " << fFragmentToReadoutMap.size() << " elements";
	for(const auto& pair : fFragmentToReadoutMap) std::cout << "   Frag: " << std::hex << pair.first << ", Crate: " 
								<< pair.second.first << ", # boards: " << std::dec << pair.second.second.size() << std::endl;
	
    }

    theClockFragmentIDs.stop();

    double fragmentIDsTime = theClockFragmentIDs.accumulated_real_time();

    cet::cpu_timer theClockReadoutIDs;

    theClockReadoutIDs.start();

    fReadoutBoardToChannelMap.clear();
    if (fChannelMappingTool->BuildTPCReadoutBoardToChannelMap(fReadoutBoardToChannelMap))
    {
        std::cout << "******* FAILED TO CONFIGURE CHANNEL MAP ********" << std::endl;
        throw cet::exception("ICARUSChannelMapProvider") << "POS didn't read the F'ing database again \n";
    }

    // Do the channel mapping initialization
    fFragmentToDigitizerMap.clear();
    if (fChannelMappingTool->BuildFragmentToDigitizerChannelMap(fFragmentToDigitizerMap))
      {
	throw cet::exception("ICARUSChannelMapProvider") << "Cannot recover the Fragment ID channel map from the database \n";
      }
    else if (fDiagnosticOutput)
      {
	std::cout << "FragmentID to Readout ID map has " << fFragmentToDigitizerMap.size() << " Fragment IDs";
	 for(const auto& pair : fFragmentToDigitizerMap) std::cout << "   Frag: " << std::hex << pair.first << ", # pairs: " 
								   << std::dec << pair.second.size() << std::endl;
      }
    
    // Do the channel mapping initialization for CRT
    fCRTChannelIDToHWtoSimMacAddressPairMap.clear();
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
    
    
    // Do the channel mapping initialization for top CRT
    fTopCRTHWtoSimMacAddressPairMap.clear();
    if (fChannelMappingTool->BuildTopCRTHWtoSimMacAddressPairMap(fTopCRTHWtoSimMacAddressPairMap))
      {
        throw cet::exception("CRTDecoder") << "Cannot recover the Top CRT HW MAC Address  from the database \n";
      }
    else if (fDiagnosticOutput)
      {
	std::cout << "Top CRT MacAddress map has " << fTopCRTHWtoSimMacAddressPairMap.size() << " rows";
        for(const auto& pair : fTopCRTHWtoSimMacAddressPairMap) std::cout << ", hw mac address: " << pair.first
									  <<", sim mac address: " << pair.second << std::endl;
      }
    

    // Do the CRT Charge Calibration initialization
    fSideCRTChannelToCalibrationMap.clear();
    if (fChannelMappingTool->BuildSideCRTCalibrationMap(fSideCRTChannelToCalibrationMap))
      {
	std::cout << "******* FAILED TO CONFIGURE CRT Calibration  ********" << std::endl;
        throw cet::exception("ICARUSChannelMapProvider") << "Cannot recover the charge calibration information from the database \n";
      }
    else if (fDiagnosticOutput)
      {
	std::cout << "side crt calibration map has " << fSideCRTChannelToCalibrationMap.size() << " list of rows \n";

	for(const auto& pair : fSideCRTChannelToCalibrationMap) std::cout <<" mac5: "<< pair.first.first
									  << ", chan: " << pair.first.second
									  << ", Gain: " << pair.second.first
									  << ", Pedestal: " << pair.second.second << std::endl;

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

const TPCReadoutBoardToChannelMap& ICARUSChannelMapProvider::getReadoutBoardToChannelMap() const
{
    return fReadoutBoardToChannelMap;
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
    return findPMTfragmentEntry(fragmentID) != nullptr;
}


unsigned int ICARUSChannelMapProvider::nPMTfragmentIDs() const {
  return fFragmentToDigitizerMap.size();
}


const DigitizerChannelChannelIDPairVec& ICARUSChannelMapProvider::getChannelIDPairVec(const unsigned int fragmentID) const
{
    mf::LogTrace{ "ICARUSChannelMapProvider" }
      << "Call to: ICARUSChannelMapProvider::getChannelIDPairVec(" << fragmentID << ")";
    DigitizerChannelChannelIDPairVec const* digitizerPair = findPMTfragmentEntry(fragmentID);
    
    if (digitizerPair) return *digitizerPair;
    throw cet::exception("ICARUSChannelMapProvider") << "Fragment ID " << fragmentID << " not found in lookup map when looking for PMT channel info \n";
    
}

  unsigned int ICARUSChannelMapProvider::getSimMacAddress(const unsigned int hwmacaddress)  const
  {
    unsigned int   simmacaddress = 0;

    for(const auto& pair : fCRTChannelIDToHWtoSimMacAddressPairMap){
      if (pair.second.first == hwmacaddress)
	simmacaddress = pair.second.second;
    }

    return simmacaddress;
  }
  
  unsigned int ICARUSChannelMapProvider::gettopSimMacAddress(const unsigned int hwmacaddress)  const
  {
    unsigned int   simmacaddress = 0;

    for(const auto& pair : fTopCRTHWtoSimMacAddressPairMap){
      if (pair.first == hwmacaddress)
	simmacaddress = pair.second;
    }

    return simmacaddress;
    
    /* untested:
    auto const it = fTopCRTHWtoSimMacAddressPairMap.find(hwmacaddress);
    return (it == fTopCRTHWtoSimMacAddressPairMap.end())? 0: it->second;
    */
  }
   
  std::pair<double, double> ICARUSChannelMapProvider::getSideCRTCalibrationMap(int mac5, int chan) const
  {
    auto const itGainAndPedestal = fSideCRTChannelToCalibrationMap.find({ mac5, chan });
    return (itGainAndPedestal == fSideCRTChannelToCalibrationMap.cend())
      ? std::pair{ -99., -99. }: itGainAndPedestal->second;
  }

auto ICARUSChannelMapProvider::findPMTfragmentEntry(unsigned int fragmentID) const
  -> DigitizerChannelChannelIDPairVec const*
{
  auto it = fFragmentToDigitizerMap.find(PMTfragmentIDtoDBkey(fragmentID));
  return (it == fFragmentToDigitizerMap.end())? nullptr: &(it->second);
}


constexpr unsigned int ICARUSChannelMapProvider::PMTfragmentIDtoDBkey
  (unsigned int fragmentID)
{
  /*
   * PMT channel mapping database stores the board number (0-23) as key.
   * Fragment ID are currently in the pattern 0x20xx, with xx the board number.
   */
  
  // protest if this is a fragment not from the PMT;
  // but make an exception for old PMT fragment IDs (legacy)
  assert(((fragmentID & ~0xFF) == 0x00) || ((fragmentID & ~0xFF) == 0x20));
  
  return fragmentID & 0xFF;
  
} // ICARUSChannelMapProvider::PMTfragmentIDtoDBkey()


constexpr unsigned int ICARUSChannelMapProvider::DBkeyToPMTfragmentID
  (unsigned int DBkey)
{
  /*
   * PMT channel mapping database stores the board number (0-23) as key.
   * Fragment ID are currently in the pattern 0x20xx, with xx the board number.
   */
  
  // protest if this is a fragment not from the PMT;
  // but make an exception for old PMT fragment IDs (legacy)
  assert((DBkey & 0xFF) < 24);
  
  return (DBkey & 0xFF) | 0x2000;
  
} // ICARUSChannelMapProvider::PMTfragmentIDtoDBkey()


} // end namespace

