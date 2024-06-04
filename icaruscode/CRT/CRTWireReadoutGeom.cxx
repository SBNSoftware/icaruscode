///////////////////////////////////////////////////////////////////////////////
/// \file CRTWireReadoutGeom.cxx
/// \brief Algorithm class for ICARUS auxiliary detector channel mapping
///
/// Originally ported from AuxDetChannelMapLArIATAlg.cxx (Author: brebel@fnal.gov)
/// and modified for SBND (Author: mastbaum@uchicago.edu) then ICARUS
/// \version $Id: 0.0 $
/// \author chilge@rams.colostate.edu
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "icaruscode/CRT/CRTWireReadoutGeom.h"
#include "TVector3.h"
#include <iostream>
#include <ostream>
#include <vector>

//namespace icarus{
//namespace crt {
namespace geo {


  //---------------------------------------------------------------------------
  CRTWireReadoutGeom::CRTWireReadoutGeom(
      fhicl::ParameterSet const& p)
    : fSorter(geo::CRTGeoObjectSorter(p)) {}

  //---------------------------------------------------------------------------
  void CRTWireReadoutGeom::Initialize(AuxDetGeometryData_t& geodata) {
    Uninitialize();

    std::vector<geo::AuxDetGeo>& adgeo = geodata.auxDets;
    std::vector<geo::AuxDetGeo*> adgeo_copy;
    for (size_t i=0; i<adgeo.size(); i++) adgeo_copy.push_back(&adgeo[i]);

    // Sort the AuxDetGeo objects and map them to names of the detectors
    // (SBN workshop 19 March) seemed to cause some problems for SBND folks --comment out?
    //    fSorter.SortAuxDets(adgeo);

    //    for (size_t a=0; a<adgeo.size(); a++){
    // std::cout << "module number:  " << a << "   volname: " << adgeo[a].TotalVolume()->GetName() << std::endl;     
    // }
    
    fSorter.SortCRTs(adgeo);

    /*
    for (size_t i=0; i<adgeo.size(); i++) {
      std::cout << "i: " << i
      << " | name1: " << adgeo[i].TotalVolume()->GetName()
      << " | name2: " << adgeo_copy[i]->TotalVolume()->GetName() << std::endl;
    }
    */


    // Map the AuxDetGeo names to their position in the sorted vector
    //
    // The ICARUS CRT is composed of three subsystems
    //   (1) The top and 'eaves' portion is new construction from CERN/Bologna.
    //       Each module consists of two layers of 8 strips in an X-Y configuration.
    //       Each strip contains 2 fibers read out seperately for a total of 32 channels/module.
    //       Each module is read out by a single CAEN front-end board (same as SBND, uBooNE).
    //   (2) The sides are made up of MINOS modules. Each module consists of 1 layer of 20 strips.
    //       Each strip contains a single fiber read out at both ends. 2 strips are optically
    //       ganged on a single SiPM on each end for a total of 20 channels/module. A single end
    //       of 3 modules are read out by a single CAEN front-end board.
    //   (3) The bottom is made up of Double Chooz veto modules. Each module consists of two layers of
    //       32 strips in an X-X configuration with the layers laterally offset by 0.5*(strip width).
    //       Each strip contains a single fiber read out at a single end by an M64 PMT for a total of
    //       64 channels / module. The front-end electronics are different from the subsystems 1 & 2.
    //
    // In the geometry, CRT modules are AuxDets and CRT strips are the
    // AuxDetSensitives. Each strip has one-two SiPM channels depending on module
    // type, one per optical fiber (not in the geometry).
    //
    // Top:                84 modules, 1344 strips, 2688 channels
    // Front, back eaves:  6 modules, 96 strips, 192 channels
    // Left, right eaves:  13 modules, 208 strips, 416 channels
    // Front, back sides:  20 modules, 400 strips, 400 channels
    // Left, right sides:  54 modules, 1080 strips, 1080 channels
    // Bottom:             14 modules, 896 strips, 896 channels
    //
    // Total: 284 modules, 5808 strips, 7760 channels
    //    CERN type:  122 modules, 122 FEBs, 3904 channels: 4052-7956
    //    MINOS type: 148 modules, 100 FEBs, 2960 channels: 0-3155 (2 channels/FEB not used)
    //    DblCh type:  14 modules,  14 FEBs,  896 channels: 3156-4051

    fADGeoToName.clear();
    fADGeoToChannelAndSV.clear();


        
    int chID = 0;

    //loop over modules
    for (size_t a=0; a<adgeo.size(); a++){

      std::string volName(adgeo[a].TotalVolume()->GetName());

      // define a list of pair of sensitive volume ID and sensitive volume
      // we are making a pair because the sorted strip number is necessary
      // while assigning each strip number to a channelID.
      std::vector<std::pair<int, geo::AuxDetSensitiveGeo> > adslist;

      //get number of strips in this module
      size_t nsv = adgeo[a].NSensitiveVolume();

      // make a pair of sensitive volume ID and sensitive volume
      for(int isv = 0; isv < (int)nsv; isv++ ){
	adslist.push_back(std::make_pair(isv,adgeo[a].SensitiveVolume(isv)));
	//std::cout << isv << "\t"  << std::endl;
      }
      //      std::cout << " before sort................. \t"  << std::endl;
      //for (size_t a=0; a<adslist.size(); a++){
      //std::cout << adslist[a].first << "\t"  << std::endl;
      //}
      // sort the sensitive volume inside a module
      fSorter.SortCRTSensitive(adslist);
      //std::cout << " after sort..................... \t"  << std::endl;
      //for (size_t a=0; a<adslist.size(); a++){
      //std::cout << adslist[a].first << "\t"  << std::endl;
      //}

      if (nsv != 20 && nsv !=16 && nsv != 64) {
        throw cet::exception("CRTChannelMap")
        << "Wrong number of sensitive volumes for CRT volume "
        << volName << " (got " << nsv << ", expected 20, 16 or 64)" << std::endl;
      }//if strip number correct

      fADGeoToName[a] = volName;
      fNameToADGeo[volName] = a;


      //-------------------------------------------------------
      
      // Example string
      std::string str =  volName;
      
      // For atoi, the input string has to start with a digit, so lets search for the first digit
      size_t i = 0;
      for ( ; i < str.length(); i++ ){ if ( std::isdigit(str[i]) ) break; }
      
      // remove the first chars, which aren't digits
      str = str.substr(i, str.length() - i );
      
      // convert the remaining text to an integer
      int id = std::atoi(str.c_str());
      
      //-------------------------------------------------------


      //if volume is a module, not a strip (search doesn't hit end)
      if (volName.find("volAuxDet") != std::string::npos) {
        //loop over strips
        //for (size_t svID=0; svID<nsv; svID++) {
      //std::cout<< "hello << -------------------------- >> " << std::endl;
          //CERN modules
          if (nsv==16){
            for (size_t svID=0; svID<16; svID++) {
              //each strip has 2 fibers, 1 SiPM per fiber at the same strip end
              //ch id in range (0,31)
              for (size_t ich=0; ich<2; ich++) {
                fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, adslist[svID].first));
		chID++;
		//std::cout << " " << a << " \t" << volName << "\t "<< id << " \t" << adslist[svID].first << " \t" << chID << std::endl;
	      }
            }
          }
          //DC modules
          if (nsv==64){
            //1 fiber per strip read out at one end
            for (size_t svID=0; svID<64; svID++) {
	      fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, adslist[svID].first));
	      chID++;
	      //    std::cout << "module index | " << a << " volName: | " << volName << "channelno: | " << chID << std::endl;      
	      //      std::cout << " " << a << " \t" << volName << "\t "<< id << " \t" << adslist[svID].first << " \t" << chID << std::endl;
            }
          }
          //MINOS modules
          if (nsv==20){

	      for (size_t svID=0; svID<20; svID+=2) {
		//  size_t chID = svID/2 + ich+31; //2 fibers on one SiPM, at both ends, w/different FEBs
		fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, adslist[svID].first));
		fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, adslist[svID+1].first));
		chID++;
		//std::cout << "module index | " << a << " volName: | " << volName <<  "\t " << adslist[svID].first <<  "\t " << adslist[svID+1].first<<" channelno: | " << chID << std::endl;      
		//	std::cout << " " << a << " \t" << volName << "\t "<< id << " \t" << adslist[svID].first << " \t" << chID << std::endl;

		if (adslist[svID].second.Length() == 800){
		  fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, adslist[svID].first));
		  fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, adslist[svID+1].first));
		  chID++;

		  //std::cout << " " << a << " \t" << volName << "\t "<< id << " \t" << adslist[svID].first << " \t" << chID << std::endl;

		  // std::cout << "module index :  800 cm | " << a << " volName: | " << volName << "\t " << adslist[svID].first <<  "\t " << adslist[svID+1].first << " channelno: | " << chID << std::endl;      
		}
	      }


	  } // if side crt strips
	  //extractIntegerWords(volName);



	  // print the result
	  // std::cout << id << std::endl;

	  //	  std::cout << "module index  " << a << " volName  " << volName << " module no "<< id << " channelID:  " << chID << std::endl;      
	  std::cout << " " << a << " \t" << volName << "\t "<< id << " \t" << chID << std::endl;      
      } //if correct name
      else {
	throw cet::exception("CRTChannelMap")
	  << "Unrecognized volume found: " << volName << std::endl;
      }
    }//for modules
    
    /*
    for (auto& csvItr : fADGeoToChannelAndSV){

      for (auto& pairItr : csvItr.second)
	std::cout <<  "module index | "<< csvItr.first << " channel: | "<< pairItr.first << " sv: | "<< pairItr.second <<std::endl;
    }
    */

  }//Initialize
  
  //----------------------------------------------------------------------------
  void CRTWireReadoutGeom::Uninitialize() {}

  //----------------------------------------------------------------------------
  uint32_t CRTWireReadoutGeom::PositionToAuxDetChannel(
      geo::Point_t const& worldLoc,
      std::vector<geo::AuxDetGeo> const& auxDets,
      size_t& ad,
      size_t& sv) const {

    // Set the default to be that we don't find the position in any AuxDet
     uint32_t channel = UINT_MAX;
     //uint32_t channel;

    // Figure out which detector we are in
    ad = 0;
    sv = this->NearestSensitiveAuxDet(worldLoc, auxDets, ad);


    //    std::cout << ad << "\t"<< auxDets[ad].TotalVolume()->GetName() << std::endl;
    // for (auto& adtoname : fADGeoToName) std::cout << "ad: "<<adtoname.first << " | volname: " << adtoname.second<<std::endl;
    // Check to see which AuxDet this position corresponds to
    auto gnItr = fADGeoToName.find(ad);
    if (gnItr != fADGeoToName.end()){
    // Get the vector of channel and sensitive volume pairs
    
      auto csvItr = fADGeoToChannelAndSV.find(ad);
      
      if (csvItr == fADGeoToChannelAndSV.end()) {
	throw cet::exception("CRTWireReadoutGeom")
	  << "No entry in channel and sensitive volume map for AuxDet index "
	  << ad;
      }
      //*/
      /*
      // N.B. This is the ID on the nth channel, and the strip has n and n+1
    auto pairItr = std::find_if((csvItr->second).begin(), (csvItr->second).end(),
				[sv](const std::pair<uint32_t,size_t>& element){ return element.second == sv;} );
    channel = pairItr->first;
    */
    // }

    

      for (auto& pairItr : csvItr->second)
	if (pairItr.second == sv) channel = pairItr.first;
    }
      
    /*if (channel == UINT_MAX) {
      throw cet::exception("CRTWireReadoutGeom")
      << "position ("
      << worldLoc[0] << "," << worldLoc[1] << "," << worldLoc[2]
      << ") does not correspond to any AuxDet";
    }*/

    return channel;
  }

  //----------------------------------------------------------------------------
  geo::Point_t CRTWireReadoutGeom::AuxDetChannelToPosition(
      uint32_t const channel,
      std::string const& auxDetName,
      std::vector<geo::AuxDetGeo> const& auxDets) const {

    // Figure out which detector we are in
    size_t ad = UINT_MAX;
    if (fNameToADGeo.count(auxDetName) > 0) {
      ad = fNameToADGeo.find(auxDetName)->second;
    }
    else {
      throw cet::exception("CRTWireReadoutGeom")
      << "No AuxDetGeo with name " << auxDetName;
    }

    // Get the vector of channel and sensitive volume pairs
    auto csvItr = fADGeoToChannelAndSV.find(ad);

    if (csvItr == fADGeoToChannelAndSV.end()) {
      throw cet::exception("CRTWireReadoutGeom")
      << "No entry in channel and sensitive volume"
      << " map for AuxDet index " << ad << " bail";
    }

    // Loop over the vector of channel and sensitive volumes to determine the
    // sensitive volume for this channel. Then get the origin of the sensitive
    // volume in the world coordinate system.
    for (auto csv : csvItr->second) {
      if (csv.first == channel) {
        // Get the center of the sensitive volume for this channel
        return auxDets[ad].SensitiveVolume(csv.second).GetCenter();
      }
    }

    return {};
  }

}  // namespace crt
//} //namespace icarus
