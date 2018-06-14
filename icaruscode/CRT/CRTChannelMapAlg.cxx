///////////////////////////////////////////////////////////////////////////////
/// \file CRTChannelMapAlg.cxx
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
#include "icaruscode/CRT/CRTChannelMapAlg.h"
#include "TVector3.h"
#include <iostream>
#include <ostream>

namespace geo {

  //---------------------------------------------------------------------------
  CRTChannelMapAlg::CRTChannelMapAlg(
      fhicl::ParameterSet const& p)
    : fSorter(geo::CRTGeoObjectSorter(p)) {}

  //---------------------------------------------------------------------------
  void CRTChannelMapAlg::Initialize(AuxDetGeometryData_t& geodata) {
    Uninitialize();

    std::vector<geo::AuxDetGeo*>& adgeo = geodata.auxDets;

    // Sort the AuxDetGeo objects and map them to names of the detectors
    // (SBN workshop 19 March) seemed to cause some problems for SBND folks --comment out?
    fSorter.SortAuxDets(adgeo);

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

    //loop over modules
    for (size_t a=0; a<adgeo.size(); a++){
      std::string volName(adgeo[a]->TotalVolume()->GetName());

      //get number of strips in this module
      size_t nsv = adgeo[a]->NSensitiveVolume();
      if (nsv != 20 && nsv !=16 && nsv != 64) {
        throw cet::exception("CRTChannelMap")
        << "Wrong number of sensitive volumes for CRT volume "
        << volName << " (got " << nsv << ", expected 20, 16 or 64)" << std::endl;
      }//if strip number correct

      fADGeoToName[a] = volName;
      fNameToADGeo[volName] = a;

      //if volName contains "Module" (search doesn't hit end)
      if (volName.find("_Module_") != std::string::npos) {
	//loop over strips
        for (size_t svID=0; svID<nsv; svID++) {
          size_t chID=UINT_MAX;
          if (nsv==16){
            chID = 32 * a + 2 * svID + 0 - 1132;
            fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, svID));
            chID = 32 * a + 2 * svID + 1 - 1132;
            fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, svID));
          }
	  if (nsv==64){
	    chID = 64 * a + svID - 6316;
	    fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, svID));
	  }
          if (nsv==20){
            chID = 32 * (a/3) + svID/2 + 10*(a % 3);
            fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, svID));
            chID = 32 * (a/3) + svID/2 + 10*(a % 3) + 1578;
            fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, svID));
          }
        }//for strips
      } //if correct name

      else {
        throw cet::exception("CRTChannelMap")
        << "Unrecognized volume found: " << volName << std::endl;
      }

    }//for modules
  }//Initialize

  //----------------------------------------------------------------------------
  void CRTChannelMapAlg::Uninitialize() {}

  //----------------------------------------------------------------------------
  uint32_t CRTChannelMapAlg::PositionToAuxDetChannel(
      double const worldLoc[3],
      std::vector<geo::AuxDetGeo*> const& auxDets,
      size_t& ad,
      size_t& sv) const {

    // Set the default to be that we don't find the position in any AuxDet
    uint32_t channel = UINT_MAX;

    // Figure out which detector we are in
    ad = 0;
    sv = this->NearestSensitiveAuxDet(worldLoc, auxDets, ad);

    // Get the origin of the sensitive volume in the world coordinate system
    double svOrigin[3] = {0, 0, 0};
    double localOrigin[3] = {0, 0, 0};

    auxDets[ad]->SensitiveVolume(sv).LocalToWorld(localOrigin, svOrigin);

    // Check to see which AuxDet this position corresponds to
    auto gnItr = fADGeoToName.find(ad);
    if (gnItr != fADGeoToName.end()){
      // Get the vector of channel and sensitive volume pairs
      auto csvItr = fADGeoToChannelAndSV.find(ad);

      if (csvItr == fADGeoToChannelAndSV.end()) {
        throw cet::exception("CRTChannelMapAlg")
        << "No entry in channel and sensitive volume map for AuxDet index "
        << ad;
      }

      // N.B. This is the ID on the nth channel, and the strip has n and n+1
      auto pairItr = std::find_if((csvItr->second).begin(), (csvItr->second).end(),
        [sv](const std::pair<uint32_t,size_t>& element){ return element.second == sv;} );
      channel = pairItr->first;
    }

    /*if (channel == UINT_MAX) {
      throw cet::exception("CRTChannelMapAlg")
      << "position ("
      << worldLoc[0] << "," << worldLoc[1] << "," << worldLoc[2]
      << ") does not correspond to any AuxDet";
    }*/

    return channel;
  }

  //----------------------------------------------------------------------------
  const TVector3 CRTChannelMapAlg::AuxDetChannelToPosition(
      uint32_t const& channel,
      std::string const& auxDetName,
      std::vector<geo::AuxDetGeo*> const& auxDets) const {
    double x = 0;
    double y = 0;
    double z = 0;

    // Figure out which detector we are in
    size_t ad = UINT_MAX;
    if (fNameToADGeo.count(auxDetName) > 0) {
      ad = fNameToADGeo.find(auxDetName)->second;
    }
    else {
      throw cet::exception("CRTChannelMapAlg")
      << "No AuxDetGeo with name " << auxDetName;
    }

    // Get the vector of channel and sensitive volume pairs
    auto csvItr = fADGeoToChannelAndSV.find(ad);

    if (csvItr == fADGeoToChannelAndSV.end()) {
      throw cet::exception("CRTChannelMapAlg")
      << "No entry in channel and sensitive volume"
      << " map for AuxDet index " << ad << " bail";
    }

    // Loop over the vector of channel and sensitive volumes to determine the
    // sensitive volume for this channel. Then get the origin of the sensitive
    // volume in the world coordinate system.
    double svOrigin[3] = {0, 0, 0};
    double localOrigin[3] = {0, 0, 0};
    for (auto csv : csvItr->second) {
      if (csv.first == channel) {
        // Get the center of the sensitive volume for this channel
        auxDets[ad]->SensitiveVolume(csv.second).LocalToWorld(localOrigin,
                                                              svOrigin);

        x = svOrigin[0];
        y = svOrigin[1];
        z = svOrigin[2];

        break;
      }
    }

    return TVector3(x, y, z);
  }

}  // namespace geo

