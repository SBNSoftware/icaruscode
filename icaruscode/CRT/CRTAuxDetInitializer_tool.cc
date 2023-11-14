#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetReadoutGeom.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "icaruscode/CRT/CRTGeoObjectSorter.h"
#include "icaruscode/CRT/compareCRTs.h"

#include "art/Utilities/ToolMacros.h"

#include "TVector3.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

  class CRTAuxDetInitializer : public geo::AuxDetInitializer {
  public:
    explicit CRTAuxDetInitializer(fhicl::ParameterSet const&) {}

  private:
    geo::AuxDetReadoutInitializers
    initialize(std::vector<geo::AuxDetGeo> const& adgeo) const override
    {
      geo::AuxDetReadoutInitializers result;
      auto& [ADGeoToName, NameToADGeo, ADGeoToChannelAndSV] = result;

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
          adslist.emplace_back(isv,adgeo[a].SensitiveVolume(isv));
        }

        // sort the sensitive volume inside a module
        std::sort(adslist.begin(), adslist.end(),
                  [](auto const& ads1, auto const& ads2) {
                    return icarus::compareCRTs(ads1.second, ads2.second);
                  });

        if (nsv != 20 && nsv !=16 && nsv != 64) {
          throw cet::exception("CRTChannelMap")
            << "Wrong number of sensitive volumes for CRT volume "
            << volName << " (got " << nsv << ", expected 20, 16 or 64)" << std::endl;
        }//if strip number correct

        ADGeoToName[a] = volName;
        NameToADGeo[volName] = a;

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
        if (volName.find("volAuxDet_") != std::string::npos) {
          //loop over strips
          //CERN modules
          if (nsv==16){
            for (size_t svID=0; svID<16; svID++) {
              //each strip has 2 fibers, 1 SiPM per fiber at the same strip end
              //ch id in range (0,31)
              for (size_t ich=0; ich<2; ich++) {
                ADGeoToChannelAndSV[a].emplace_back(chID, adslist[svID].first);
                chID++;
              }
            }
          }
          //DC modules
          if (nsv==64){
            //1 fiber per strip read out at one end
            for (size_t svID=0; svID<64; svID++) {
              ADGeoToChannelAndSV[a].emplace_back(chID, adslist[svID].first);
              chID++;
            }
          }
          //MINOS modules
          if (nsv==20){

            for (size_t svID=0; svID<20; svID+=2) {
              ADGeoToChannelAndSV[a].emplace_back(chID, adslist[svID].first);
              ADGeoToChannelAndSV[a].emplace_back(chID, adslist[svID+1].first);
              chID++;

              if (adslist[svID].second.Length() == 800){
                ADGeoToChannelAndSV[a].emplace_back(chID, adslist[svID].first);
                ADGeoToChannelAndSV[a].emplace_back(chID, adslist[svID+1].first);
                chID++;
              }
            }
          } // if side crt strips

          // print the result
          std::cout << " " << a << " \t" << volName << "\t "<< id << " \t" << chID << std::endl;
        } //if correct name
        else {
          throw cet::exception("CRTChannelMap")
            << "Unrecognized volume found: " << volName << std::endl;
        }
      } //for modules

      return result;
    }
  };

}

DEFINE_ART_CLASS_TOOL(CRTAuxDetInitializer)
