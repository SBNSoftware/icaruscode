#ifndef IC_CRTCOMMONUTILS_CC
#define IC_CRTCOMMONUTILS_CC

#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

using namespace icarus::crt;

//given an AuxDetGeo object, returns name of the CRT subsystem to which it belongs
char CRTCommonUtils::GetAuxDetType(geo::AuxDetGeo const& adgeo) {
    string volName(adgeo.TotalVolume()->GetName());
    if (volName.find("MINOS") != string::npos) return 'm';
    if (volName.find("CERN")  != string::npos) return 'c';
    if (volName.find("DC")    != string::npos) return 'd';

    mf::LogError("CRT") << "AuxDetType not found!" << '\n';
    return 'e';

}

//------------------------------------------------------------------------------------
//given an AuxDetGeo object, returns name of the CRT region to which it belongs
string CRTCommonUtils::GetAuxDetRegion(geo::AuxDetGeo const& adgeo) {
    char type = GetAuxDetType(adgeo);
    string base = "volAuxDet_", region="";
    switch ( type ) {
        case 'c' : base+= "CERN"; break;
        case 'd' : base+= "DC"; break;
        case 'm' : base+= "MINOS"; break;
    }
    base+="_module_###_";
    string volName(adgeo.TotalVolume()->GetName());
  
    //module name has 2 possible formats
    //  volAuxDet_<subsystem>_module_###_<region>
    //  volAuxDet_<subsystem>_module_###_cut###_<region>
  
    region = volName.substr(base.length(),volName.length());
    if( region.find("_")==string::npos)
        return region;
  
    else
        return region.substr(region.find("_")+1,region.length());
  
}

//------------------------------------------------------------------------------
int CRTCommonUtils::GetAuxDetRegionNum(string reg)
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
string CRTCommonUtils::GetRegionNameFromNum(int num) {
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

//-------------------------------------------------------------------------------------------
char CRTCommonUtils::GetRegTypeFromRegName(string name) {

    //CERN modules
    if(name=="Top"||name=="RimWest"||name=="RimEast"||name=="RimSouth"||name=="RimNorth")
        return 'c';
    //MINOS modules
    if(name=="WestSouth"||name=="WestCenter"||name=="WestNorth"||
       name=="EastSouth"||name=="EastCenter"||name=="EastNorth"||
       name=="North"||name=="South")
        return 'm';
    //DC modules
    if(name=="Bottom")
        return 'd';

    //error
    return 'e';
}


//---------------------------------------------------------------------------------------------
//for C- and D-modules, 1 mac5 address
//three M-modules / FEB, full-length modules read out at both ends (2 FEBs)
//  cut modules read out at one end (1 FEB)
//  numbering convention is module from FEB i 
//  return pair<FEB i,FEB i> (C-, D-, Cut M- module)
//  return pair<FEB i,FEB j> (Full-length M-modules)
pair<uint8_t,uint8_t> CRTCommonUtils::ADToMac(const map<int,vector<pair<uint8_t,int>>>& febMap, uint32_t adid) {

    for(auto const& p : febMap) {
        if((uint32_t)p.first!=adid)
            continue;
        if(p.second.size()==2)
            return std::make_pair(p.second[0].first,p.second[1].first);
        else
            return std::make_pair(p.second[0].first,p.second[0].first);
    }
    return std::make_pair(UINT8_MAX,UINT8_MAX);
}


//--------------------------------------------------------------------------------------
int CRTCommonUtils::MacToRegion(uint8_t mac){

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
char CRTCommonUtils::MacToType(uint8_t mac)
{     
   int reg = MacToRegion(mac);
   if(reg>=30&&reg<40) return 'c';
   if(reg>=40&&reg<50) return 'm';
   if(reg==50) return 'd';
   std::cout << "ERROR in CRTHitRecoAlg::MacToType: type not set!" << std::endl;
   return 'e';
}
//--------------------------------------------------------------------------------------

string CRTCommonUtils::MacToRegionName(uint8_t mac)
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

int CRTCommonUtils::ChannelToAuxDetSensitiveID(uint8_t mac, int chan) {
  char type = MacToType(mac);
  if (type=='d') return chan;
  if (type=='c') return chan/2;
  if (type=='m') return (chan % 10)*2;

  return INT_MAX;
}

//--------------------------------------------------------------------------------------

int CRTCommonUtils::MacToAuxDetID(uint8_t mac, int chan)
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

//-------------------------------------------------------------------------------
int CRTCommonUtils::ModToTypeCode(geo::AuxDetGeo const& adgeo) {
   size_t nstrips = adgeo.NSensitiveVolume();
   if (nstrips==16) return 0; //'c'
   if (nstrips==20) return 1; //'m'
   if (nstrips==64) return 2; //'d'
   return INT_MAX;
}

//--------------------------------------------------------------------------------------
//reads a file generated by CRT geometry generation script and
// returns a map modID->vector<pair<FEB, FEB channel subset>,+1 dual readoutfor MINOS module>
// channel subset = 1,2,or 3 (always =1 for c or d modules)
map<int,vector<pair<uint8_t,int>>> CRTCommonUtils::GetFebMap() {

    map<int,vector<pair<uint8_t,int>>> febMap;

    string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file("feb_map.txt",fullFileName);
    std::ifstream fin;
    fin.open(fullFileName,std::ios::in);

    if(fin.good()) 
        std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else
        throw cet::exception("CRTDetSim::FillFebMap") 
          << "Unable to find/open file 'feb_map.txt'" << std::endl;

    vector<string> row;
    string line, word;
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

//-----------------------------------------------------------------------
//returns average 4-position in the scintillator strip
//ROOT::Math::XYZTVector
TLorentzVector CRTCommonUtils::AvgIDEPoint(sim::AuxDetIDE ide){
    double x = 0.5*(ide.entryX + ide.exitX);
    double y = 0.5*(ide.entryY + ide.exitY);
    double z = 0.5*(ide.entryZ + ide.exitZ);
    double t = 0.5*(ide.entryT + ide.exitT);
    //ROOT::Math::XYZTVector
    TLorentzVector v(x,y,z,t);

    return v;
}

//-----------------------------------------------------------------------
//returns the path length in the scintillator strip
double CRTCommonUtils::LengthIDE(sim::AuxDetIDE ide) {
    double dx=0., dy=0., dz=0.;
    dx = ide.entryX - ide.exitX;
    dy = ide.entryY - ide.exitY;
    dz = ide.entryZ - ide.exitZ;
    return sqrt(dx*dx+dy*dy+dz*dz);
}

//----------------------------------------------------------------------
int CRTCommonUtils::GetLayerID(geo::GeometryCore const* geoService, sim::AuxDetSimChannel const& adsc){
    int layer = -1;

    auto const& adGeo = geoService->AuxDet(adsc.AuxDetID());
    auto const& adsGeo = adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());
    int region = GetAuxDetRegionNum(GetAuxDetRegion(adGeo));
    int type = ModToTypeCode(adGeo);

    std::set<string> volNames = { adsGeo.TotalVolume()->GetName() };
    vector<vector<TGeoNode const*> > paths = geoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
            path += "/";
        }
    }
    TGeoManager* manager = geoService->ROOTGeoManager();
    manager->cd(path.c_str());
    TGeoNode* nodeStrip = manager->GetCurrentNode();
    TGeoNode* nodeInner = manager->GetMother(1);
    TGeoNode* nodeModule = manager->GetMother(2);
    double origin[3] = {0, 0, 0};
    double modulePosMother[3]; //position in CRT region volume
    double stripPosMother[3]; // strip position in module frame
    double stripPosModule[3];

    nodeModule->LocalToMaster(origin, modulePosMother);
    nodeStrip->LocalToMaster(origin, stripPosMother);
    nodeInner->LocalToMaster(stripPosMother,stripPosModule);

    //if 'c' or 'd' type
    if ( type == 0 || type == 2 )
        layer = (stripPosModule[1] > 0);

    // if 'm' type
    if ( type == 1 ) {
        // if east or west stacks (6 in total)
        if ( region >=40 && region <=45 ) {
            layer = ( modulePosMother[0]>0 );
        }
        // if north stack
        if ( region == 47) {
            layer = ( modulePosMother[2]> 0 );
        }
        // if south stack
        if( region == 46) {
            auto const& adsGeo = adGeo.SensitiveVolume(0);
            if(adsGeo.Length() < 500) //is cut module?
                layer = 1;
            else 
                layer = 0;
        }
    }

    return layer;

}

//----------------------------------------------------------------------
int CRTCommonUtils::GetLayerID(geo::GeometryCore const* geoService, const art::Ptr<sim::AuxDetSimChannel> adsc){
    int layer = -1;

    auto const& adGeo = geoService->AuxDet(adsc->AuxDetID());
    auto const& adsGeo = adGeo.SensitiveVolume(adsc->AuxDetSensitiveID());
    int region = GetAuxDetRegionNum(GetAuxDetRegion(adGeo));
    int type = ModToTypeCode(adGeo);

    std::set<string> volNames = { adsGeo.TotalVolume()->GetName() };
    vector<vector<TGeoNode const*> > paths = geoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
            path += "/";
        }
    }
    TGeoManager* manager = geoService->ROOTGeoManager();
    manager->cd(path.c_str());
    TGeoNode* nodeStrip = manager->GetCurrentNode();
    TGeoNode* nodeInner = manager->GetMother(1);
    TGeoNode* nodeModule = manager->GetMother(2);
    double origin[3] = {0, 0, 0};
    double modulePosMother[3]; //position in CRT region volume
    double stripPosMother[3]; // strip position in module frame
    double stripPosModule[3];

    nodeModule->LocalToMaster(origin, modulePosMother);
    nodeStrip->LocalToMaster(origin, stripPosMother);
    nodeInner->LocalToMaster(stripPosMother,stripPosModule);

    //if 'c' or 'd' type
    if ( type == 0 || type == 2 )
        layer = (stripPosModule[1] > 0);

    // if 'm' type
    if ( type == 1 ) {
        // if east or west stacks (6 in total)
        if ( region >=40 && region <=45 ) {
            layer = ( modulePosMother[0]>0 );
        }
        // if north stack
        if ( region == 47) {
            layer = ( modulePosMother[2]> 0 );
        }
        // if south stack
        if( region == 46) {
            if(adsGeo.Length() < 500) //is cut module?
                layer = 1;
            else 
                layer = 0;
        }
    }

    return layer;

}

//----------------------------------------------------------------------------------------------------------
int CRTCommonUtils::GetMINOSLayerID(geo::GeometryCore const* geoService, geo::AuxDetGeo const& adGeo) {
    int layer = -1;

    int region = GetAuxDetRegionNum(GetAuxDetRegion(adGeo));
    int type = ModToTypeCode(adGeo);
    if(type!=1) {
        mf::LogError("CRTCommonUtils") << "non-MINOS module provided to GetMINOSLayerID";
        return layer;
    }

    std::set<string> volNames = { adGeo.TotalVolume()->GetName() };
    vector<vector<TGeoNode const*> > paths = geoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
            path += "/";
        }
    }
    TGeoManager* manager = geoService->ROOTGeoManager();
    manager->cd(path.c_str());
    TGeoNode* nodeModule = manager->GetCurrentNode();
    double origin[3] = {0, 0, 0};
    double modulePosMother[3]; //position in CRT region volume

    nodeModule->LocalToMaster(origin, modulePosMother);

    // if east or west stacks (6 in total)
    if ( region >=40 && region <=45 ) {
        layer = ( modulePosMother[0]>0 );
    }
    // if north stack
    if ( region == 47) {
        layer = ( modulePosMother[2]> 0 );
    }
    // if south stack
    if( region == 46) {
        auto const& adsGeo = adGeo.SensitiveVolume(0);
        if(adsGeo.Length() < 500) //is cut module?
            layer = 1;
        else
            layer = 0;
    }

    return layer;


}

//-----------------------------------------------------------------------------------------------------------
// given mac address and mac channel, return CRT strip center in module coordinates (w.r.t. module center)
TVector3 CRTCommonUtils::ChanToLocalCoords(geo::GeometryCore const* geoService, uint8_t mac, int chan) {

    TVector3 coords(0.,0.,0.);
    int adid  = MacToAuxDetID(mac,chan); //CRT module ID
    auto const& adGeo = geoService->AuxDet(adid); //CRT module
    int adsid = ChannelToAuxDetSensitiveID(mac,chan); //CRT strip ID
    auto const& adsGeo = adGeo.SensitiveVolume(adsid); //CRT strip

    double origin[3] = {0,0,0};
    double stripPosWorld[3], modPos[3];

    adsGeo.LocalToWorld(origin,stripPosWorld);
    adGeo.WorldToLocal(stripPosWorld,modPos);

    coords.SetXYZ(modPos[0],modPos[1],modPos[2]);
    return coords;
}

//-----------------------------------------------------------------------------------------------------------
// given mac address and mac channel, return CRT strip center in World coordinates (w.r.t. LAr active volume center)
TVector3 CRTCommonUtils::ChanToWorldCoords(geo::GeometryCore const* geoService, uint8_t mac, int chan) {

    TVector3 coords(0.,0.,0.);
    int adid  = MacToAuxDetID(mac,chan); //CRT module ID
    auto const& adGeo = geoService->AuxDet(adid); //CRT module
    int adsid = ChannelToAuxDetSensitiveID(mac,chan); //CRT strip ID
    auto const& adsGeo = adGeo.SensitiveVolume(adsid); //CRT strip

    double origin[3] = {0,0,0};
    double stripPosWorld[3];

    adsGeo.LocalToWorld(origin,stripPosWorld);

    coords.SetXYZ(stripPosWorld[0],stripPosWorld[1],stripPosWorld[2]);
    return coords;
}

#endif
