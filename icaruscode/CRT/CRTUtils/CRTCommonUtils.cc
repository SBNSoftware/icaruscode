#ifndef IC_CRTCOMMONUTILS_CC
#define IC_CRTCOMMONUTILS_CC

#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

using namespace icarus::crt;

CRTCommonUtils::CRTCommonUtils() {
  // fGeoService = art::ServiceHandle<geo::AuxDetGeometry const>()->GetProviderPtr();
  fGeoService  = lar::providerFrom<geo::Geometry>();
  FillFebMap();
  FillAuxDetMaps();
}

//given an AuxDetGeo object, returns name of the CRT subsystem to which it belongs
char CRTCommonUtils::GetAuxDetType(size_t adid) {
    if(fAuxDetIdToType.find(adid)==fAuxDetIdToType.end()) {
        throw cet::exception("CRTCommonUtils::GetAuxDetType")
          << "unknown AuxDetID passed to function";
    }
    return fAuxDetIdToType[adid];
}

//------------------------------------------------------------------------------------
int CRTCommonUtils::GetAuxDetTypeCode(size_t adid) {
    char type = GetAuxDetType(adid);
    switch(type){
        case 'c': return 0;
        case 'm': return 1;
        case 'd': return 2;
        default:  return -1;
    }
}

//------------------------------------------------------------------------------------
//given an AuxDetGeo object, returns name of the CRT region to which it belongs
string CRTCommonUtils::GetAuxDetRegion(size_t adid) {
    if(fAuxDetIdToRegion.find(adid)==fAuxDetIdToRegion.end()) {
        throw cet::exception("CRTCommonUtils::GetAuxDetRegion")
          << "unknown AuxDetID passed to function";
    }
    return fAuxDetIdToRegion[adid];  
}

//------------------------------------------------------------------------------
int CRTCommonUtils::AuxDetRegionNameToNum(string reg)
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
    throw cet::exception("CRTCommonUtils::GetRegTypeFromRegName")
          << "passed region name not recognized!" << '\n';
}

//-------------------------------------------------------------------------------------------
int CRTCommonUtils::GetTypeCodeFromRegion(string name) {
    char type = GetRegTypeFromRegName(name);
    switch(type){
        case 'c': return 0;
        case 'm': return 1;
        case 'd': return 2;
        default:  return -1;
    }
}

//---------------------------------------------------------------------------------------------
//for C- and D-modules, 1 mac5 address
//three M-modules / FEB, full-length modules read out at both ends (2 FEBs)
//  cut modules read out at one end (1 FEB)
//  numbering convention is module from FEB i 
//  return pair<FEB i,FEB i> (C-, D-, Cut M- module)
//  return pair<FEB i,FEB j> (Full-length M-modules)
pair<uint8_t,uint8_t> CRTCommonUtils::ADToMac(size_t adid) {

    if(fAuxDetIdToFeb.find(adid)==fAuxDetIdToFeb.end()) {
        throw cet::exception("CRTCommonUtils::ADToMac")
          << "unknown AuxDetID passed to function";
    }
    vector<pair<uint8_t,int>> febs = fAuxDetIdToFeb[adid];
    if(febs.size()==2)
        return std::make_pair(febs[0].first,febs[1].first);
    else
        return std::make_pair(febs[0].first,febs[0].first);
}

//---------------------------------------------------------------------------------
int CRTCommonUtils::ADToChanGroup(size_t adid){
    if(fAuxDetIdToChanGroup.find(adid)==fAuxDetIdToChanGroup.end()){
        throw cet::exception("CRTCommonUtils::ADMacToChanGroup")
          << "unknown AuxDetID passed to function";
    }
    return fAuxDetIdToChanGroup[adid];
}

//-------------------------------------------------------------------------------------
int CRTCommonUtils::NFeb(size_t adid){
    if(fAuxDetIdToFeb.find(adid)==fAuxDetIdToFeb.end()) {
        throw cet::exception("CRTCommonUtils::NFeb")
          << "unknown AuxDetID passed to function";
    }
    return(fAuxDetIdToFeb[adid]).size();
}

//--------------------------------------------------------------------------------------
string CRTCommonUtils::MacToRegion(uint8_t mac){
    if(fFebToAuxDetId.find(mac)==fFebToAuxDetId.end()) {
        throw cet::exception("CRTCommonUtils::MacToRegion")
          << "unknown mac passed to function";
    }
    return GetAuxDetRegion(fFebToAuxDetId[mac][0]);
}

//--------------------------------------------------------------------------------------
char CRTCommonUtils::MacToType(uint8_t mac)
{ 
    if(fFebToAuxDetId.find(mac)==fFebToAuxDetId.end()) {
        throw cet::exception("CRTCommonUtils::MacToType")
          << "unknown mac passed to function";
    }
    return GetAuxDetType(fFebToAuxDetId[mac][0]);
}

//--------------------------------------------------------------------------------------
int CRTCommonUtils::MacToTypeCode(uint8_t mac)
{
    if(fFebToAuxDetId.find(mac)==fFebToAuxDetId.end()) {
        throw cet::exception("CRTCommonUtils::MacToType")
          << "unknown mac passed to function";
    }
    char type = GetAuxDetType(fFebToAuxDetId[mac][0]);
    switch(type){
        case 'c': return 0;
        case 'm': return 1;
        case 'd': return 2;
        default:  return -1;
    }
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

size_t CRTCommonUtils::MacToAuxDetID(uint8_t mac, int chan)
{
    char type = MacToType(mac);
    int pos=1;
    

    /*
    if(type=='m'){
      if (chan >=0 && chan < 12) pos = 1;
      else if (chan >=12 && chan < 22) pos = 2;
      else if (chan >=22 && chan < 32) pos = 3;
      //else pos = 1; 
    }
    */
    if(type=='m') pos = chan/10 + 1;


    
    if(fFebToAuxDetId.find(mac)==fFebToAuxDetId.end()) {
      throw cet::exception("CRTCommonUtils::MacToAuxDetID")
	<< "unknown mac passed to function";
    }
    
    
    for(auto const& adid : fFebToAuxDetId[mac]) {
      
      if(fAuxDetIdToChanGroup.find(adid)==fAuxDetIdToChanGroup.end()){
	throw cet::exception("CRTCommonUtils::ADMacToChanGroup")
	  << "unknown AuxDetID passed to function";
      }
      

      //std::cout << "mac: "  << (int)mac << " ,module number " << adid <<  ", pos " << fAuxDetIdToChanGroup[adid] 
      //	<<  ", chan " << chan <<", pos to match " << pos << std::endl;
      //std::cout << " good upto Step 1" << std::endl;
      /*
      if (mac == 49 && pos == 3){ 
	//	std::cout << " good upto Step 2" << std::endl;
	continue;
      } else if (fAuxDetIdToChanGroup[adid]==pos){
	//std::cout << " good upto Step 3" << std::endl;
	return adid;
      } else {
	//std::cout << " good upto Step 4" << std::endl;
	return 0;
      }
      
      */

      if (fAuxDetIdToChanGroup[adid]==pos){
	return adid;
      }
      //	 std::cout << "mac: "  << (int)mac << ", chan " << chan << ", pos " << pos << ", adid " << adid<< std::endl;
    }
    
    //    std::cout << "------------------------- mac: "  << (int)mac << ", chan " << chan << ", pos " << pos << std::endl;
    
    throw cet::exception("CRTCommonUtils::MacToAuxDetID")
      << "AuxDetID not found!";
    
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
int CRTCommonUtils::GetLayerID(sim::AuxDetSimChannel const& adsc){
    int layer = -1;

    auto const& adGeo = fGeoService->AuxDet(adsc.AuxDetID());
    auto const& adsGeo = adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());
    int region = AuxDetRegionNameToNum(GetAuxDetRegion(adsc.AuxDetID()));
    char type = GetAuxDetType(adsc.AuxDetID());

    std::set<string> volNames = { adsGeo.TotalVolume()->GetName() };
    vector<vector<TGeoNode const*> > paths = fGeoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
            path += "/";
        }
    }
    TGeoManager* manager = fGeoService->ROOTGeoManager();
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
    if ( type == 'c' || type == 'd' )
        layer = (stripPosModule[1] > 0);

    // if 'm' type
    if ( type == 'm' ) {
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
	  //            if(adsGeo.Length() < 500) //is cut module?
	  if(adsGeo.Length() == 400 || adsGeo.Length() == 485.15) //is cut module?
	    layer = 1;
	  else 
	    layer = 0;
        }
    }

    return layer;

}

//----------------------------------------------------------------------
int CRTCommonUtils::GetLayerID(const art::Ptr<sim::AuxDetSimChannel> adsc){
    int layer = -1;

    auto const& adGeo = fGeoService->AuxDet(adsc->AuxDetID());
    auto const& adsGeo = adGeo.SensitiveVolume(adsc->AuxDetSensitiveID());
    int region = AuxDetRegionNameToNum(GetAuxDetRegion(adsc->AuxDetID()));
    char type = GetAuxDetType(adsc->AuxDetID());

    std::set<string> volNames = { adsGeo.TotalVolume()->GetName() };
    vector<vector<TGeoNode const*> > paths = fGeoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
            path += "/";
        }
    }
    TGeoManager* manager = fGeoService->ROOTGeoManager();
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
    if ( type == 'c' || type == 'd' )
        layer = (stripPosModule[1] > 0);

    // if 'm' type
    if ( type == 'm' ) {
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
	  //if(adsGeo.Length() < 500) //is cut module?
	  if(adsGeo.Length() == 400 || adsGeo.Length() == 485.15) //is cut module?
                layer = 1;
            else 
                layer = 0;
        }
    }

    return layer;

}

//--------------------------------------------------------------------------------------------------
int CRTCommonUtils::GetMINOSLayerID(size_t adid) {
    int layer = -1;

    int region = AuxDetRegionNameToNum(GetAuxDetRegion(adid));
    char type = GetAuxDetType(adid);
    auto const& adGeo = fGeoService->AuxDet(adid);
    if(type!='m') {
        mf::LogError("CRTCommonUtils") << "non-MINOS module provided to GetMINOSLayerID";
        return layer;
    }

    std::set<string> volNames = { adGeo.TotalVolume()->GetName() };
    vector<vector<TGeoNode const*> > paths = fGeoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
            path += "/";
        }
    }
    TGeoManager* manager = fGeoService->ROOTGeoManager();
    manager->cd(path.c_str());
    TGeoNode* nodeModule = manager->GetCurrentNode(); //Mother(2);
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
	if(adsGeo.Length() == 400 || adsGeo.Length() == 485.15) //is cut module?
	//        if(adsGeo.Length() < 500) //is cut module?
            layer = 1;
        else
            layer = 0;
    }

    if(layer==-1)
        mf::LogError("CRTCommonUtils::GetMINOSLayerID")
           << "layer ID not set!";

    return layer;


}

//--------------------------------------------------------------------------------------------------
// given mac address and mac channel, return CRT strip center in module coordinates (w.r.t. module center)
TVector3 CRTCommonUtils::ChanToLocalCoords(uint8_t mac, int chan) {

    TVector3 coords(0.,0.,0.);
    size_t adid  = MacToAuxDetID(mac,chan); //CRT module ID
    auto const& adGeo = fGeoService->AuxDet(adid); //CRT module
    int adsid = ChannelToAuxDetSensitiveID(mac,chan); //CRT strip ID
    auto const& adsGeo = adGeo.SensitiveVolume(adsid); //CRT strip

    double origin[3] = {0,0,0};
    double stripPosWorld[3], modPos[3];

    adsGeo.LocalToWorld(origin,stripPosWorld);
    adGeo.WorldToLocal(stripPosWorld,modPos);

    coords.SetXYZ(modPos[0],modPos[1],modPos[2]);
    return coords;
}

//--------------------------------------------------------------------------------------------------
// given mac address and mac channel, return CRT strip center in World coordinates (w.r.t. LAr active volume center)
TVector3 CRTCommonUtils::ChanToWorldCoords(uint8_t mac, int chan) {

    TVector3 coords(0.,0.,0.);
    int adid  = MacToAuxDetID(mac,chan); //CRT module ID
    auto const& adGeo = fGeoService->AuxDet(adid); //CRT module
    int adsid = ChannelToAuxDetSensitiveID(mac,chan); //CRT strip ID
    auto const& adsGeo = adGeo.SensitiveVolume(adsid); //CRT strip

    double origin[3] = {0,0,0};
    double stripPosWorld[3];

    adsGeo.LocalToWorld(origin,stripPosWorld);

    coords.SetXYZ(stripPosWorld[0],stripPosWorld[1],stripPosWorld[2]);
    return coords;
}

//--------------------------------------------------------------------------------------
//reads a file generated by CRT geometry generation script and
// fills a map modID->vector<pair<FEB, FEB channel subset>,+1 dual readoutfor MINOS module>
// channel subset = 1,2,or 3 (always =1 for c or d modules)
void CRTCommonUtils::FillFebMap() {

    string fullFileName;// = "/icarus/app/users/chilgenb/ana_icaruscode_v08_52_00/feb_map.txt";
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file("feb_map.txt",fullFileName);
    std::ifstream fin;
    fin.open(fullFileName,std::ios::in);

    if(!fin.good())
//        std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
//    else
        throw cet::exception("CRTDetSim::FillFebMap")
          << "Unable to find/open file 'feb_map.txt'" << std::endl;

    vector<string> row;
    string line, word;
    //each line has pattern 
    //  auxDetID, mac5, chan pos, '\n' (CERN, DC, cut MINOS) OR
    //  auxDetID, mac5, chan pos, mac5', chan pos', '\n' (full length MINOS)
    while(getline(fin,line)) {
        row.clear();
        std::stringstream s(line);
        while (std::getline(s, word, ',')) { //parse line
            row.push_back(word);
        }
        int mod = (size_t)std::stoi(row[0]);       //auxDetID
        uint8_t mac5 = (uint8_t)std::stoi(row[1]); //febID
        int pos = std::stoi(row[2]);               //feb channel block
        fAuxDetIdToFeb[mod].push_back(std::make_pair(mac5,pos));
        fAuxDetIdToChanGroup[mod]=pos;
        fFebToAuxDetId[mac5].push_back(mod);
//        std::cout << "mod: " << mod << ", mac: " << (int)mac5 << ", pos: " << pos;
        if(row.size()>3) { //if dual ended readout MINOS module
            mac5 = (uint8_t)std::stoi(row[3]);
            fAuxDetIdToFeb[mod].push_back(std::make_pair(mac5,pos));
            fFebToAuxDetId[mac5].push_back(mod);
//            if(pos!=std::stoi(row[4])) //feb channel block same on both febs
//              std::cout << "WARNING in CRTComUtil: 2 unique chan groups for ADId!" << std::endl;
//            std::cout << ", mac: " << (int)mac5;
        }
//        std::cout << std::endl;
    }
//    std::cout << "filled febMap with " << fAuxDetIdToFeb.size() << " entries" << std::endl;
    fin.close();
}

//------------------------------------------------------------------------
void CRTCommonUtils::FillAuxDetMaps() {

    for(auto const& ad : fAuxDetIdToFeb){
        auto const& adGeo = fGeoService->AuxDet(ad.first);

        //AuxDetType
        switch(adGeo.NSensitiveVolume()) {
            case 16:
                fAuxDetIdToType[ad.first] = 'c';
                break;
            case 20:
                fAuxDetIdToType[ad.first] = 'm';
                break;
            case 64:
                fAuxDetIdToType[ad.first] = 'd';
                break;
        }

        //AuxDet region
        string name = adGeo.TotalVolume()->GetName();
        fNameToAuxDetId[name] = ad.first;
        fAuxDetIdToRegion[ad.first] = AuxDetNameToRegion(name);
    }

}

//--------------------------------------------------------------------
string CRTCommonUtils::AuxDetNameToRegion(string name) {

    string base("volAuxDet_");
    const char type = fAuxDetIdToType[fNameToAuxDetId[name]];

    switch( type ){
      case 'c' : base+= "CERN"; break;
      case 'd' : base+= "DC";   break;
      case 'm' : base+= "MINOS"; break;
      default  : 
          throw cet::exception("CRTCommonUtils::Constructor::AuxDetNameToRegion")
                << "AuxDet type not set!";
    }
    base+="_module_###_";

    //module name has 2 possible formats
    //  volAuxDet_<subsystem>_module_###_<region>
    //  volAuxDet_<subsystem>_module_###_cut###_<region>
    string region(name.substr(base.length(),name.length()));
    if( region.find("_")==string::npos) 
        return region;
    
    else 
        return region.substr(region.find("_")+1,region.length());
}

//--------------------------------------------------------------------------
TVector3 CRTCommonUtils::WorldToModuleCoords(TVector3 point, size_t adid) {

    char type = GetAuxDetType(adid);
    auto const& adGeo = fGeoService->AuxDet(adid);
    double world[3], local[3];
    world[0] = point.X();
    world[1] = point.Y();
    world[2] = point.Z();

    adGeo.WorldToLocal(world,local);
    TVector3 localpoint(local[0],local[1],local[2]);

    if(type=='c') {
        localpoint.SetX( 8.*(localpoint.X()/adGeo.HalfWidth1()+1.));
        localpoint.SetY(0);
        localpoint.SetZ( 8.*(localpoint.Z()/adGeo.HalfWidth1()+1.));
        return localpoint;
    }
    if(type=='m') { //needs fixing
        localpoint.SetX(0);
        localpoint.SetY(ADToChanGroup(adid)*10.*(localpoint.Y()/adGeo.HalfHeight()+1.));
        localpoint.SetZ(localpoint.Z());
        return localpoint;
    }
    if(type=='d') {
        localpoint.SetX(32.5*(localpoint.X()/adGeo.HalfWidth1()+1.));
        localpoint.SetY(0);
        localpoint.SetZ(localpoint.Z());
        return localpoint;
    }

    return localpoint;
}

// Simple distance of closest approach between infinite track and centre of hit
double CRTCommonUtils::SimpleDCA(sbn::crt::CRTHit hit, TVector3 start, TVector3 direction){
  
  TVector3 pos (hit.x_pos, hit.y_pos, hit.z_pos);
  TVector3 end = start + direction;
  double denominator = direction.Mag();
  double numerator = (pos - start).Cross(pos - end).Mag();
  return numerator/denominator;

}

// Minimum distance from infinite track to CRT hit assuming that hit is a 2D square
double CRTCommonUtils::DistToCrtHit(sbn::crt::CRTHit hit, TVector3 start, TVector3 end){

  // Check if track goes inside hit
  TVector3 min (hit.x_pos - hit.x_err, hit.y_pos - hit.y_err, hit.z_pos - hit.z_err);
  TVector3 max (hit.x_pos + hit.x_err, hit.y_pos + hit.y_err, hit.z_pos + hit.z_err);
  if(CubeIntersection(min, max, start, end).first.X() != -99999) return 0;

  // Calculate the closest distance to each edge of the CRT hit
  // Assume min error is the fixed position of tagger
  TVector3 vertex1 (hit.x_pos, hit.y_pos - hit.y_err, hit.z_pos - hit.z_err);
  TVector3 vertex2 (hit.x_pos, hit.y_pos + hit.y_err, hit.z_pos - hit.z_err);
  TVector3 vertex3 (hit.x_pos, hit.y_pos - hit.y_err, hit.z_pos + hit.z_err);
  TVector3 vertex4 (hit.x_pos, hit.y_pos + hit.y_err, hit.z_pos + hit.z_err);
  if(hit.y_err < hit.x_err && hit.y_err < hit.z_err){
    vertex1.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos, hit.z_pos - hit.z_err);
    vertex2.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos, hit.z_pos - hit.z_err);
    vertex3.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos, hit.z_pos + hit.z_err);
    vertex4.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos, hit.z_pos + hit.z_err);
  }
  if(hit.z_err < hit.x_err && hit.z_err < hit.y_err){
    vertex1.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos - hit.y_err, hit.z_pos);
    vertex2.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos - hit.y_err, hit.z_pos);
    vertex3.SetXYZ(hit.x_pos - hit.x_err, hit.y_pos + hit.y_err, hit.z_pos);
    vertex4.SetXYZ(hit.x_pos + hit.x_err, hit.y_pos + hit.y_err, hit.z_pos);
  }

  double dist1 = LineSegmentDistance(vertex1, vertex2, start, end);
  double dist2 = LineSegmentDistance(vertex1, vertex3, start, end);
  double dist3 = LineSegmentDistance(vertex4, vertex2, start, end);
  double dist4 = LineSegmentDistance(vertex4, vertex3, start, end);

  return std::min(std::min(dist1, dist2), std::min(dist3, dist4));

}


// Distance between infinite line (2) and segment (1)
// http://geomalgorithms.com/a07-_distance.html
double CRTCommonUtils::LineSegmentDistance(TVector3 start1, TVector3 end1, TVector3 start2, TVector3 end2){

  double smallNum = 0.00001;

  // 1 is segment
  TVector3 direction1 = end1 - start1;
  // 2 is infinite line
  TVector3 direction2 = end2 - start2;

  TVector3 u = direction1;
  TVector3 v = direction2;
  TVector3 w = start1 - start2;

  double a = u.Dot(u);
  double b = u.Dot(v);
  double c = v.Dot(v);
  double d = u.Dot(w);
  double e = v.Dot(w);
  double D = a * c - b * b;
  double sc, sN, sD = D; // sc = sN/sD
  double tc, tN, tD = D; // sc = sN/sD

  // Compute the line parameters of the two closest points
  if(D < smallNum){ // Lines are almost parallel
    sN = 0.0;
    sD = 1.0;
    tN = e;
    tD = c;
  }
  else{
    sN = (b * e - c * d)/D;
    tN = (a * e - b * d)/D;
    if(sN < 0.){ // sc < 0, the s = 0 edge is visible
      sN = 0.;
      tN = e;
      tD = c;
    }
    else if(sN > sD){ // sc > 1, the s = 1 edge is visible
      sN = sD;
      tN = e + b;
      tD = c;
    } 
  }

  sc = (std::abs(sN) < smallNum ? 0.0 : sN / sD);
  tc = (std::abs(tN) < smallNum ? 0.0 : tN / tD);
  // Get the difference of the two closest points
  TVector3 dP = w + (sc * u) - (tc * v);

  return dP.Mag();

}

// Intersection between axis-aligned cube and infinite line
// (https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection)
std::pair<TVector3, TVector3> CRTCommonUtils::CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end){

  TVector3 dir = (end - start);
  TVector3 invDir (1./dir.X(), 1./dir.Y(), 1/dir.Z());

  double tmin, tmax, tymin, tymax, tzmin, tzmax;

  TVector3 enter (-99999, -99999, -99999);
  TVector3 exit (-99999, -99999, -99999);

  // Find the intersections with the X plane
  if(invDir.X() >= 0){
    tmin = (min.X() - start.X()) * invDir.X();
    tmax = (max.X() - start.X()) * invDir.X();
  }
  else{
    tmin = (max.X() - start.X()) * invDir.X();
    tmax = (min.X() - start.X()) * invDir.X();
  }

  // Find the intersections with the Y plane
  if(invDir.Y() >= 0){
    tymin = (min.Y() - start.Y()) * invDir.Y();
    tymax = (max.Y() - start.Y()) * invDir.Y();
  }
  else{
    tymin = (max.Y() - start.Y()) * invDir.Y();
    tymax = (min.Y() - start.Y()) * invDir.Y();
  }

  // Check that it actually intersects
  if((tmin > tymax) || (tymin > tmax)) return std::make_pair(enter, exit);

  // Max of the min points is the actual intersection
  if(tymin > tmin) tmin = tymin;

  // Min of the max points is the actual intersection
  if(tymax < tmax) tmax = tymax;

  // Find the intersection with the Z plane
  if(invDir.Z() >= 0){
    tzmin = (min.Z() - start.Z()) * invDir.Z();
    tzmax = (max.Z() - start.Z()) * invDir.Z();
  }
  else{
    tzmin = (max.Z() - start.Z()) * invDir.Z();
    tzmax = (min.Z() - start.Z()) * invDir.Z();
  }

  // Check for intersection
  if((tmin > tzmax) || (tzmin > tmax)) return std::make_pair(enter, exit);

  // Find final intersection points
  if(tzmin > tmin) tmin = tzmin;

  // Find final intersection points
  if(tzmax < tmax) tmax = tzmax;

  // Calculate the actual crossing points
  double xmin = start.X() + tmin * dir.X();
  double xmax = start.X() + tmax * dir.X();
  double ymin = start.Y() + tmin * dir.Y();
  double ymax = start.Y() + tmax * dir.Y();
  double zmin = start.Z() + tmin * dir.Z();
  double zmax = start.Z() + tmax * dir.Z();

  // Return pair of entry and exit points
  enter.SetXYZ(xmin, ymin, zmin);
  exit.SetXYZ(xmax, ymax, zmax);
  return std::make_pair(enter, exit);

}


#endif
