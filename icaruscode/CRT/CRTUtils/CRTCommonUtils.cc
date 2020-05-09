#ifndef IC_CRTCOMMONUTILS_CC
#define IC_CRTCOMMONUTILS_CC

#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

namespace icarus {
 namespace crt {

    //given an AuxDetGeo object, returns name of the CRT subsystem to which it belongs
    char CRTCommonUtils::GetAuxDetType(geo::AuxDetGeo const& adgeo) {
        std::string volName(adgeo.TotalVolume()->GetName());
        if (volName.find("MINOS") != std::string::npos) return 'm';
        if (volName.find("CERN")  != std::string::npos) return 'c';
        if (volName.find("DC")    != std::string::npos) return 'd';

        mf::LogError("CRT") << "AuxDetType not found!" << '\n';
        return 'e';

    }

    //------------------------------------------------------------------------------------
    //given an AuxDetGeo object, returns name of the CRT region to which it belongs
    std::string CRTCommonUtils::GetAuxDetRegion(geo::AuxDetGeo const& adgeo) {
        char type = CRTCommonUtils::GetAuxDetType(adgeo);
        std::string base = "volAuxDet_", region="";
        switch ( type ) {
          case 'c' : base+= "CERN"; break;
          case 'd' : base+= "DC"; break;
          case 'm' : base+= "MINOS"; break;
        }
        base+="_module_###_";
        std::string volName(adgeo.TotalVolume()->GetName());
      
        //module name has 2 possible formats
        //  volAuxDet_<subsystem>_module_###_<region>
        //  volAuxDet_<subsystem>_module_###_cut###_<region>
      
        region = volName.substr(base.length(),volName.length());
        if( region.find("_")==std::string::npos)
          return region;
      
        else
              return region.substr(region.find("_")+1,region.length());
      
    }

    //------------------------------------------------------------------------------
    int CRTCommonUtils::GetAuxDetRegionNum(std::string reg)
    {
        if(reg == "Top")        return 30;
        if(reg == "RimWest")    return 31;
        if(reg == "RimEast")    return 32;
        if(reg == "RimSouth")   return 33;
        if(reg == "RimNorth")   return 34;
        if(reg == "WestSouth")  return 40;
        if(reg == "WestCenter") return 41;
        if(reg == "WestNorth")  return 42;
        if(reg == "EastSouth")  return 43;
        if(reg == "EastCenter") return 44;
        if(reg == "EastNorth")  return 45;
        if(reg == "South")      return 46;
        if(reg == "North")      return 47;
        if(reg == "Bottom")     return 50;
        mf::LogError("CRT") << "region not found!" << '\n';
        return INT_MAX;
    }//GetAuxDetRegionNum()

    //---------------------------------------------------------------------------------
    std::string CRTCommonUtils::GetRegionNameFromNum(int num) {
        switch(num) {
            case 30 :
                return "Top";
            case 31 :
                return "RimWest";
            case 32 :
                return "RimEast";
            case 33 :
                return "RimSouth";
            case 34 :
                return "RimNorth";
            case 40 :
                return "WestSouth";
            case 41 :
                return "WestCenter";
            case 42 :
                return "WestNorth";
            case 43 :
                return "EastSouth";
            case 44 :
                return "EastCenter";
            case 45 :
                return "EastNorth";
            case 46 :
                return "South";
            case 47 :
                return "North";
            case 50 :
                return "Bottom";
        }
    
        return "Region not found!";
    }

    //--------------------------------------------------------------------------------------
    int CRTCommonUtils::MacToRegion(int mac){

        if(mac>=107 && mac<=190) return 30; //top
        if(mac>=191 && mac<=204) return 31; //rim west
        if(mac>=205 && mac<=218) return 32; //rim east
        if(mac>=219 && mac<=224) return 33; //rim south
        if(mac>=225 && mac<=230) return 34; //rim north
        if(            mac<=12 ) return 40; //west side, south stack
        if(mac>=13  && mac<=24 ) return 41; //west side, center stack
        if(mac>=25  && mac<=36 ) return 42; //west side, north stack
        if(mac>=37  && mac<=48 ) return 43; //east side, south stack
        if(mac>=49  && mac<=60 ) return 44; //east side, center stack
        if(mac>=61  && mac<=72 ) return 45; //east side, north stack
        if(mac>=73  && mac<=84 ) return 46; //south
        if(mac>=85  && mac<=92 ) return 47; //north
        if(mac>=93 && mac<=106) return 50; //bottom

        std::cout << "ERROR in CRTHitRecoAlg::MacToRegion: unknown mac address " << mac << std::endl;
        return 0;
    }

    //--------------------------------------------------------------------------------------
    char CRTCommonUtils::MacToType(int mac)
    {     
       int reg = MacToRegion(mac);
       if(reg>=30&&reg<40) return 'c';
       if(reg>=40&&reg<50) return 'm';
       if(reg==50) return 'd';
       std::cout << "ERROR in CRTHitRecoAlg::MacToType: type not set!" << std::endl;
       return 'e';
    }
    //--------------------------------------------------------------------------------------

    std::string CRTCommonUtils::MacToRegionName(int mac)
    {
        int reg = MacToRegion(mac);
        switch(reg) {
            case 30 : return "top";
            case 31 : return "rimWest";
            case 32 : return "rimEast";
            case 33 : return "rimSouth";
            case 34 : return "rimNorth";
            case 40 : return "westSouth";
            case 41 : return "westCenter";
            case 42 : return "westNorth";
            case 43 : return "eastSouth";
            case 44 : return "eastCenter";
            case 45 : return "eastNorth";
            case 46 : return "south";
            case 47 : return "north";
            case 50 : return "bottom";
        }
        return "";
    }

    //--------------------------------------------------------------------------------------
    
    int CRTCommonUtils::ChannelToAuxDetSensitiveID(int mac, int chan) {
      char type = MacToType(mac);
      if (type=='d') return chan;
      if (type=='c') return chan/2;
      if (type=='m') return (chan % 10)*2;
  
      return INT_MAX;
    }

    //--------------------------------------------------------------------------------------
    
    int CRTCommonUtils::MacToAuxDetID(int mac, int chan)
    {
	auto febmap = GetFebMap();

        char type = MacToType(mac);
        if (type == 'e') return INT_MAX;
  
        int pos=1;
        if(type=='m')
            pos = chan/10 + 1;
  
        for(const auto&  p : febmap) {
            if(p.second[0].first == mac && p.second[0].second==pos)
                return (uint32_t)p.first;
            if(p.second.size()==2)
                if(p.second[1].first==mac && p.second[1].second==pos)
                    return (uint32_t)p.first;
        }
  
  
      std::cout << "ERROR in CRTHitRecoAlg::MacToAuxDetID: auxDetID not set!" << std::endl;
      return INT_MAX;
    }


    //--------------------------------------------------------------------------------------
    //reads a file generated by CRT geometry generation script and
    // returns a map modID->FEB(s), FEB channel subset(s). can be 2 if MINOS module (not cut)
    std::map<int,std::vector<std::pair<int,int>>> CRTCommonUtils::GetFebMap() {

	std::map<int,std::vector<std::pair<int,int>>> febMap;

        std::string fullFileName;
        cet::search_path searchPath("FW_SEARCH_PATH");
        searchPath.find_file("feb_map.txt",fullFileName);
        std::ifstream fin;
        fin.open(fullFileName,std::ios::in);

        if(fin.good()) 
            std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
        else
            throw cet::exception("CRTDetSim::FillFebMap") 
              << "Unable to find/open file 'feb_map.txt'" << std::endl;

        std::vector<std::string> row;
        std::string line, word;
        while(getline(fin,line)) {
            row.clear();
            std::stringstream s(line);
            int mod;
            while (std::getline(s, word, ',')) {
                row.push_back(word);
            }
            mod = std::stoi(row[0]);
            febMap[mod].push_back(std::make_pair(std::stoi(row[1]),std::stoi(row[2])));
            if(row.size()>3)
                febMap[mod].push_back(std::make_pair(std::stoi(row[3]),std::stoi(row[4])));
        }
        std::cout << "filled febMap with " << febMap.size() << " entries" << std::endl;
        fin.close();

	return febMap;

    }

 }//namespace crt
}//namespace icarus


#endif
