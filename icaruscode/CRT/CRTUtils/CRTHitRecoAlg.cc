#include "CRTHitRecoAlg.h"

namespace icarus {

  std::map<int,std::vector<std::pair<int,int>>> CRTHitRecoAlg::fFebMap;

  CRTHitRecoAlg::CRTHitRecoAlg(const Config& config){
    this->reconfigure(config);
  
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTrigClock = fDetectorClocks->TriggerClock();
    //fAuxDetGeo = &(*fAuxDetGeoService);
    //fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();
  }

  CRTHitRecoAlg::CRTHitRecoAlg(){
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTrigClock = fDetectorClocks->TriggerClock();
    //fAuxDetGeo = &(*fAuxDetGeoService);
    //fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();
  }


  CRTHitRecoAlg::~CRTHitRecoAlg(){

  }

  void CRTHitRecoAlg::reconfigure(const Config& config){
    fUseReadoutWindow = config.UseReadoutWindow(); 
    fQPed = config.QPed();
    fQSlope = config.QSlope();
    fPropDelay = config.PropDelay();
    return;
  }

  std::vector<std::pair<crt::CRTHit, std::vector<int>>> CRTHitRecoAlg::CreateCRTHits(std::vector<art::Ptr<crt::CRTData>> crtList) {

    //std::cout << "inside CreateCRTHits" << std::endl;
    std::vector<std::pair<crt::CRTHit, std::vector<int>>> returnHits;
    std::vector<int> dataIds;
    float pes =0.0;
    std::vector<uint8_t> tfeb_id;// = {0};
    std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap; //maps febID->vector<pair<channel,pes>>
    //tpesmap[0] = {std::make_pair(0,0)};

    uint16_t nHitMiss = 0, nHitC = 0, nHitD = 0, nHitM = 0;
    if (fVerbose) mf::LogInfo("CRT") << "Found " << crtList.size() << " FEB events" << '\n';

    double origin[3] = {0,0,0};
    std::map<int,int> regCounts;
    std::set<int> regs;

    //for (auto const &febdat : (*crtList)) {
    for (size_t febdat_i=0; febdat_i<crtList.size(); febdat_i++) {

      // all FEB entries have at least one strip hit
      int mac = crtList[febdat_i]->Mac5(); //febdat.Mac5();
      char type = MacToType(mac);
      int region = MacToRegion(mac);
      std::string regionName = MacToRegionName(mac);
      std::pair<int,int> trigpair = crtList[febdat_i]->TrigPair(); //pair of strips that provided the trigger (for c and d mods)
      std::pair<int,int> macPair = crtList[febdat_i]->MacPair(); // pair of FEBs that provided validation (for m mods)

      dataIds.clear();
      tfeb_id.clear();
      tpesmap.clear();

      dataIds.push_back(febdat_i);
      tfeb_id.push_back(mac);

      if ((regs.insert(region)).second) regCounts[region] = 1;
      else regCounts[region]++;

      int adid  = MacToAuxDetID(mac,trigpair.first); //module ID
      //std::cout << "auxDetID: " << adid << std::endl;
      auto const& adGeo = fGeometryService->AuxDet(adid); //module
      int adsid1 = ChannelToAuxDetSensitiveID(mac,trigpair.first); //trigger strip ID
      //std::cout << "auxDetSensitiveID: " << adsid1 << std::endl;
      auto const& adsGeo1 = adGeo.SensitiveVolume(adsid1); //trigger strip

      double stripPosWorld1[3], stripPosWorld2[3], stripPosWorld3[3], stripPosWorld4[3];
      double modPos1[3], modPos2[3], hitLocal[3];

      adsGeo1.LocalToWorld(origin,stripPosWorld1);
      adGeo.WorldToLocal(stripPosWorld1,modPos1);

      double length = adsGeo1.Length(); //strip length
      double width = adsGeo1.HalfWidth1()*2; //strip width
      double height = adsGeo1.HalfHeight()*2; //strip height

      double tCor1=0, tCor2=0;
      double t01=0, t02=0, t01corr=0, t02corr=0;
      double t11=0, t12=0, t11corr=0, t12corr=0;
      double t0=0, t0corr=0, t1=0, t1corr=0;

      double hitPoint[3];
      double hitPointErr[3];
      std::vector<crt::CRTChannelData> chandat = crtList[febdat_i]->ChanData();
      //std::set<int> hitTrack;// = chandat[0].TrackID()[0];
      //int hitMod = adid;
      //int hitStrip = adsid1;
      double distToReadoutX=0;
      double distToReadoutY=0;

      //CERN and DC modules have intermodule coincidence
      if ( type == 'c' || type == 'd' ) {
          //time and charge information for each channel above threshold
          //std::vector<icarus::crt::CRTChannelData> chandat = febdat.ChanData();
          hitLocal[1] = 0.0;

          //2nd strip within the module that provided the coincidence
          int  adsid2 = ChannelToAuxDetSensitiveID(mac,trigpair.second);
          auto const& adsGeo2 = adGeo.SensitiveVolume(adsid2);

          adsGeo2.LocalToWorld(origin,stripPosWorld2); //fetch origin of strip in global coords
          adGeo.WorldToLocal(stripPosWorld2,modPos2); //coincidence strip pos in module coords

          //CERN modules have X-Y configuration
          if ( type == 'c' ) {

              if ((adsid1<8&&adsid2<8)||(adsid1>7&&adsid2>7))
                  mf::LogInfo("CRT") << "incorrect CERN trig pair!!" << '\n';

              //if first strip along z-direction, second strip along x-dir
              if (adsid1 < 8 ) {
                  hitLocal[0] = modPos1[0];
                  hitLocal[2] = modPos2[2];
                  distToReadoutX = abs( adsGeo1.HalfLength() - modPos1[2]);
                  distToReadoutY = abs( adsGeo1.HalfLength() - modPos2[0]);
                  tCor1 = distToReadoutX*fPropDelay;
                  tCor2 = distToReadoutY*fPropDelay;
              }
              else {
                  hitLocal[0] = modPos2[0];
                  hitLocal[2] = modPos1[2];
                  distToReadoutX = abs( adsGeo1.HalfLength() - modPos2[2]);
                  distToReadoutY = abs( adsGeo1.HalfLength() - modPos1[0]);
                  tCor1 = distToReadoutY*fPropDelay;
                  tCor2 = distToReadoutX*fPropDelay;

              }

              adGeo.LocalToWorld(hitLocal,hitPoint); //tranform from module to world coords

              if(regionName=="top") {
                  hitPointErr[0] = width/sqrt(12);
                  hitPointErr[1] = height/sqrt(12);
                  hitPointErr[2] = width/sqrt(12);
              }
              else if(regionName=="rimWest" || regionName=="rimEast" {
                  hitPointErr[0] = height/sqrt(12);
                  hitPointErr[1] = width/sqrt(12);
                  hitPointErr[2] = width/sqrt(12);
              }
              else if(regionName=="rimSouth" || regionName=="rimNorth"){
                  hitPointErr[0] = width/sqrt(12);
                  hitPointErr[1] = width/sqrt(12);;
                  hitPointErr[2] = height/sqrt(12); 
              }
              else {
                  hitPointErr[0] = 0.0;
                  hitPointErr[1] = 0.0;
                  hitPointErr[2] = 0.0;
		  std::cout << "ERROR: no valid region name found for c-module. Hit errors set to 0!" << std::endl;
              }
                  
              nHitC++;

          }//if c type


          if ( type == 'd' ) {

             if ((adsid1<32&&adsid2<32)||(adsid1>31&&adsid2>31))
                  mf::LogInfo("CRT") << "incorrect DC trig pair!!" << '\n';

             hitPoint[0] = 0.5 * ( stripPosWorld1[0] + stripPosWorld2[0] );
             hitPoint[1] = 0.5 * ( stripPosWorld1[1] + stripPosWorld2[1] );
             hitPoint[2] = 0.5 * ( stripPosWorld1[2] + stripPosWorld2[2] );
             
             hitPointErr[1] = height/sqrt(12);
             //if strips are parallel to z-dir
             if(abs(stripPosWorld1[2])>300) {
                 hitPointErr[0] = 1.5*width/sqrt(12);
                 hitPointErr[2] = length/sqrt(12);
             }
             //else they are parallel to x-dir
             else {
                 hitPointErr[0] = length/sqrt(12);
                 hitPointErr[2] = 1.5*width/sqrt(12);
             }
             nHitD++;

             distToReadoutX = abs( adsGeo1.HalfLength() - modPos1[2]);
             distToReadoutY = abs( adsGeo1.HalfLength() - modPos2[2]);
             tCor1 = length*0.5*fPropDelay;
             tCor2 = tCor1;
          } // if d type

          bool cdp1=false, cdp2=false;
          //std::cout << "looping over chandat" << std::endl;
          for(auto chan : chandat) {

              //for ( auto const& trk : chan.TrackID() )
              //    hitTrack.insert(trk);

              if (chan.Channel() == trigpair.first) {
                  t01 = chan.T0();
                  t11 = chan.T1();
                  t01corr = t01-tCor1;
                  t11corr = t11-tCor1;
                  pes+=0.5*chan.ADC();
                  tpesmap[mac].push_back(std::make_pair(chan.Channel(),chan.ADC()));
                  cdp1=true;
              }
              if (chan.Channel() == trigpair.second) {
                  t02 = chan.T0();
                  t12 = chan.T1();
                  t02corr = t02-tCor2;
                  t12corr = t12-tCor2;
                  pes+=0.5*chan.ADC();
                  tpesmap[mac].push_back(std::make_pair(chan.Channel(),chan.ADC()));
                  cdp2=true;
              }
              if (t01!=0 && t02!=0) {
                  t0 = 0.5 * ( t01 + t02 );
                  t1 = 0.5 * ( t11 + t12 );
                  t0corr = 0.5 * ( t01corr + t02corr );
                  t1corr = 0.5 * ( t11corr + t12corr );
                  break;
              }
          }//for ChannelData
          if(!cdp1||!cdp2)
              std::cout << "c or d trig pair data not found!" << '\n'
                        << "  trigPair: (" << trigpair.first << ", " << trigpair.second << ")" << std::endl;


          if(fVerbose) mf::LogInfo("CRT") << " CRT HIT PRODUCED!" << '\n'
                  << "   AuxDetID: " << adid << '\n'
                  << "   AuxDetSID1 / AuxDetSID2 : " << adsid1 << " / " << adsid2 << '\n'
                  << "   x: " << hitPoint[0] << " , y: " << hitPoint[1] << " , z: " << hitPoint[2] << '\n'
                  << "   xerr: " << hitPointErr[0] << " , yerr: " << hitPointErr[1] << " , zerr: " << hitPointErr[2] << '\n'
                  << "   t0: " << t0 << " , t0corr: " << t0corr << " , t1: " << t1 << " , t1corr: " << t1corr << '\n';


          // Create a CRT hit
          crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pes, t0corr, t1corr, 0, hitPoint[0], hitPointErr[0], 
                                          hitPoint[1], hitPointErr[1], hitPoint[2], hitPointErr[2], regionName);

          returnHits.push_back(std::make_pair(crtHit, dataIds));


      }//if c or d type

      if ( type == 'm' ) {

          double ttrig1 = crtList[febdat_i]->TTrig();
          bool coinModFound = false;

          for (size_t febdat_j=0; febdat_j<crtList.size(); febdat_j++) {
              if(febdat_i==febdat_j) continue;

              int mac2 = crtList[febdat_j]->Mac5();
              if (macPair.second != mac2) continue;
              //if (MacToRegion(mac2) != region) continue;

              //std::cout << " found M mac pair: " << mac << ", " << mac2 << '\n'
              //          << "   in region " << region << std::endl;
              auto const& adsGeo3 = adGeo.SensitiveVolume(adsid1+1);
              adsGeo3.LocalToWorld(origin,stripPosWorld3);

              auto const& trigpair2 = crtList[febdat_j]->TrigPair();
              int adid2  = MacToAuxDetID(mac2,trigpair2.first);
              auto const& adGeo2 = fGeometryService->AuxDet(adid2);
              int  adsid2 = ChannelToAuxDetSensitiveID(mac2,trigpair2.first);
              auto const& adsGeo2 = adGeo2.SensitiveVolume(adsid2);
              adsGeo2.LocalToWorld(origin,stripPosWorld2);
              auto const& adsGeo4 = adGeo2.SensitiveVolume(adsid2+1);
              adsGeo4.LocalToWorld(origin,stripPosWorld4);

              //std::cout << "   giving ad1/ads1: " << adid  << " / " << adsid1 << '\n'
              //          << "          ad2/ads2: " << adid2 << " / " << adsid2 << std::endl;

              double ttrig2 = crtList[febdat_j]->TTrig();
              std::vector<crt::CRTChannelData> chandat2 = crtList[febdat_j]->ChanData();
              double pes1=0.0, pes2=0.0;
              int chan1=-1, chan2=-1;
              if (util::absDiff(ttrig2,ttrig1)<50) {

                  for(auto chan : chandat) {
                      //for ( auto const& trk : chan.TrackID() )
                          //hitTrack.insert(trk);
                      if (chan.Channel()== trigpair.first) {
                          pes1=chan.ADC();
                          chan1=chan.Channel();
                      }
                  }
                  for(auto chan : chandat2) {
                      //for ( auto const& trk : chan.TrackID() )
                          //hitTrack.insert(trk);
                      if (chan.Channel()== trigpair.second) {
                          pes2=chan.ADC();
                          chan2=chan.Channel();
                      }
                  }

                    //if east, west, or north regions (X-X)
                    if ((region >=40&&region<=45) || region==47) {
                        //average over 2 channels (4strips)
                        hitPoint[0] = 0.25*(stripPosWorld1[0]+stripPosWorld2[0]+stripPosWorld3[0]+stripPosWorld4[0]);
                        hitPoint[1] = 0.25*(stripPosWorld1[1]+stripPosWorld2[1]+stripPosWorld3[1]+stripPosWorld4[1]);
                        hitPoint[2] = 0.25*(stripPosWorld1[2]+stripPosWorld2[2]+stripPosWorld3[2]+stripPosWorld4[2]);

                        hitPointErr[0] = abs(stripPosWorld1[0] - stripPosWorld2[0])/sqrt(12);
                        hitPointErr[1] = abs(stripPosWorld1[1] - stripPosWorld2[1])/sqrt(12);
                        hitPointErr[2] = length/sqrt(12);

                        t0 = 0.5*(ttrig1+ttrig2);
                        t1 = t0;
                        t0corr = t0-length*0.5*fPropDelay;
                        t1corr = t0corr;

                        coinModFound = true;
                    }// if left or right

                    //if south region (X-Y)
                    if (region == 46 ) {
                        if  (adid<104 || adid>130)
                            mf::LogInfo("CRT") << "WARNING: unexpected MINOS mod ID for south wall" << '\n';

                        if (adid>=122&&adid<=130) { //if X module
                            hitPoint[0] = stripPosWorld1[1];
                            hitPoint[1] = stripPosWorld2[0];
                            hitPoint[2] = 0.5 * ( stripPosWorld1[2] + stripPosWorld2[2] );
                        }
                        else {
                            hitPoint[1] = stripPosWorld1[1];
                            hitPoint[0] = stripPosWorld2[0];
                            hitPoint[2] = 0.5 * ( stripPosWorld1[2] + stripPosWorld2[2] );
                        }

                        hitPointErr[0] = width/sqrt(12);
                        hitPointErr[1] = width/sqrt(12);
                        hitPointErr[2] = 0.5*util::absDiff( stripPosWorld1[2],stripPosWorld2[2]);

                        t0 = 0.5*(ttrig1+ttrig2);
                        t1 = t0;
                        t0corr = t0-length*0.5*fPropDelay;
                        t1corr = t0corr;

                        coinModFound = true;
                    }//if front or back

              }//if trigger pair module found

              if (coinModFound) {
                  nHitM++;

                  tfeb_id.push_back(mac2);
                  tpesmap[mac].push_back(std::make_pair(chan1,pes1));
                  tpesmap[mac2].push_back(std::make_pair(chan2,pes2));
                  pes = 0.5*(pes1+pes2);

                  // Create a CRT hit
                  crt::CRTHit crtHit = FillCrtHit(tfeb_id, tpesmap, pes, t0corr, t1corr, 0, hitPoint[0], hitPointErr[0],
                                          hitPoint[1], hitPointErr[1], hitPoint[2], hitPointErr[2], regionName);
                  dataIds.push_back(febdat_j);

                  returnHits.push_back(std::make_pair(crtHit, dataIds));
                  break;
              }

          }// for febdat         

          if (!coinModFound) nHitMiss++;

      } //if m type

    }//loop over CRTData products

    if(fVerbose) {
          mf::LogInfo("CRT") << returnHits.size() << " CRT hits produced!" << '\n'
              << "  nHitC: " << nHitC << " , nHitD: " << nHitD << " , nHitM: " << nHitM << '\n'
              << "    " << nHitMiss << " CRT hits missed!" << '\n';
          std::map<int,int>::iterator cts = regCounts.begin();
          mf::LogInfo("CRT") << " CRT Hits by region" << '\n';
          while (cts != regCounts.end()) {
              mf::LogInfo("CRT") << "reg: " << (*cts).first << " , hits: " << (*cts).second << '\n';
              cts++;
          }
    }//if Verbose


    return returnHits;

  }

  // Function to make filling a CRTHit a bit faster
  icarus::crt::CRTHit CRTHitRecoAlg::FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                              std::vector<std::pair<int,float>>> tpesmap, float peshit, double time0, double time1, int plane, 
                              double x, double ex, double y, double ey, double z, double ez, std::string tagger){
      icarus::crt::CRTHit crtHit;
      crtHit.feb_id      = tfeb_id;
      crtHit.pesmap      = tpesmap;
      crtHit.peshit      = peshit;
      crtHit.ts0_s_corr  = time0*1e-9; 
      crtHit.ts0_ns      = time0;
      crtHit.ts0_ns_corr = time0; 
      crtHit.ts1_ns      = time1;
      crtHit.ts0_s       = time0 * 1e-9;
      crtHit.plane       = plane;
      crtHit.x_pos       = x;
      crtHit.x_err       = ex;
      crtHit.y_pos       = y; 
      crtHit.y_err       = ey;
      crtHit.z_pos       = z;
      crtHit.z_err       = ez;
      crtHit.tagger      = tagger;

      //std::cout << "produced hit with t1 = " << crtHit.ts1_ns << std::endl;
 
      return crtHit;

  } // CRTHitRecoAlg::FillCrtHit()

  void CRTHitRecoAlg::FillFEBMap() 
  {
    if(!(this->fFebMap).empty())
        return;

    std::string dir = "/icarus/app/users/chilgenb/dev_areas/v08_22_00_prof/srcs/icaruscode/icaruscode/Geometry/gdml/";
    std::ifstream fin;
    fin.open(dir+"feb_map.txt",std::ios::in);

    if(fin.good()) std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else std::cout << "could not open file 'feb_map.txt' for reading!" << std::endl;

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
        (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[1]),std::stoi(row[2])));
        if(row.size()>3){
            (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[3]),std::stoi(row[4])));
        }             
    }
    std::cout << "filled febMap with " << (this->fFebMap).size() << " entries" << std::endl;
    fin.close();
  }

  int CRTHitRecoAlg::MacToRegion(int mac){

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

  char CRTHitRecoAlg::MacToType(int mac) 
  {

      int reg = MacToRegion(mac);
      if(reg>=30&&reg<40) return 'c';
      if(reg>=40&&reg<50) return 'm';
      if(reg==50) return 'd';
      std::cout << "ERROR in CRTHitRecoAlg::MacToType: type not set!" << std::endl;
      return 'e';
  }

  std::string CRTHitRecoAlg::MacToRegionName(int mac)
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

  int CRTHitRecoAlg::ChannelToAuxDetSensitiveID(int mac, int chan) {
    char type = MacToType(mac);
    if (type=='d') return chan;
    if (type=='c') return chan/2;
    if (type=='m') return (chan % 10)*2;

    return INT_MAX;
  }

  int CRTHitRecoAlg::MacToAuxDetID(int mac, int chan)
  {
      FillFEBMap();
      char type = MacToType(mac);
      if (type == 'e') return INT_MAX;

      int pos=1;
      if(type=='m')
          pos = chan/10 + 1;

      if((this->fFebMap).empty()) 
          std::cout << "ERROR in CRTHitRecoAlg::MacToAuxDetID: FEBMap is empty!" << std::endl;

      for(const auto&  p : this->fFebMap) {
          if(p.second[0].first == mac && p.second[0].second==pos)
              return (uint32_t)p.first;
          if(p.second.size()==2)
              if(p.second[1].first==mac && p.second[1].second==pos)
                  return (uint32_t)p.first;
      }


    std::cout << "ERROR in CRTHitRecoAlg::MacToAuxDetID: auxDetID not set!" << std::endl;
    return INT_MAX;
  }

} //end namespace icarus
