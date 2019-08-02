#include "CRTHitRecoAlg.h"

namespace icarus {

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
      auto const& adGeo = fGeometryService->AuxDet(adid); //module
      int adsid1 = ChannelToAuxDetSensitiveID(mac,trigpair.first); //trigger strip ID
      auto const& adsGeo1 = adGeo.SensitiveVolume(adsid1); //trigger strip

      double stripPosWorld1[3], stripPosWorld2[3], stripPosWorld3[3], stripPosWorld4[3];
      double modPos1[3], modPos2[3], hitLocal[3];

      adsGeo1.LocalToWorld(origin,stripPosWorld1);
      adGeo.WorldToLocal(stripPosWorld1,modPos1);

      double length = adsGeo1.Length(); //strip length
      double width = adsGeo1.HalfWidth1()*2; //strip width

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
                  hitPointErr[1] = 2*adGeo.HalfHeight()/sqrt(12);
                  hitPointErr[2] = width/sqrt(12);
              }
              else if(regionName=="rimWest" || regionName=="rimEast" {
                  hitPointErr[0] = 2*adGeo.HalfHeight()/sqrt(12);
                  hitPointErr[1] = width/sqrt(12);
                  hitPointErr[2] = width/sqrt(12);
              }
              else if(regionName=="rimSouth" || regionName=="rimNorth"){
                  hitPointErr[0] = width/sqrt(12);
                  hitPointErr[1] = adGeo.HalfHeight();
                  hitPointErr[2] = 2*adGeo.HalfHeight()/sqrt(12); 
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
             hitPointErr[1] = 2*adGeo.HalfHeight()/sqrt(12);
             if(abs(stripPosWorld1[2])>300) {
                 hitPointErr[0] = abs(stripPosWorld1[0] - stripPosWorld2[0])/sqrt(12);
                 hitPointErr[2] = length/sqrt(12);
             }
             nHitD++;

             distToReadoutX = abs( adsGeo1.HalfLength() - modPos1[2]);
             distToReadoutY = abs( adsGeo1.HalfLength() - modPos2[2]);
             tCor1 = length*0.5*fPropDelay;
             tCor2 = tCor1;
          } // if d type


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
              }
              if (chan.Channel() == trigpair.second) {
                  t02 = chan.T0();
                  t12 = chan.T1();
                  t02corr = t02-tCor2;
                  t12corr = t12-tCor2;
                  pes+=0.5*chan.ADC();
                  tpesmap[mac].push_back(std::make_pair(chan.Channel(),chan.ADC()));
              }
              if (t01!=0 && t02!=0) {
                  t0 = 0.5 * ( t01 + t02 );
                  t1 = 0.5 * ( t11 + t12 );
                  t0corr = 0.5 * ( t01corr + t02corr );
                  t1corr = 0.5 * ( t11corr + t12corr );
                  break;
              }
          }//for ChannelData


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

                    //if left or right regions (X-X)
                    if (region == 50 || region == 54) {
                        //hitPoint[0] = 0.5*(stripPosWorld1[0]+stripPosWorld2[0]);
                        //hitPoint[1] = 0.5*(stripPosWorld1[1]+stripPosWorld2[1]);
                        //hitPoint[2] = 0.5*(stripPosWorld1[2]+stripPosWorld2[2]);

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

                    //if front or back regions (X-Y)
                    if (region == 44 || region == 42) {
                        if  ((adid>107&&adid<119&&adid2>107&&adid2<119)||
                             (adid>118&&adid<128&&adid2>118&&adid2<128)||
                             (adid>127&&adid<139&&adid2>127&&adid2<139)||
                             (adid>138&&adid<148&&adid2>138&&adid2<148) )
                            mf::LogInfo("CRT") << "incorrect MINOS mac pair!!!" << '\n'
                                               << " mac 1/2: " << mac << ", " << mac2 << "\n"
                                               << " auxDetID 1/2: " << adid << ", " << adid2 << "\n";

                        if (adid>107&&adid<119) { //if X module
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

  char CRTHitRecoAlg::MacToType(int mac) {

      if(mac<=99) return 'm';
      if(mac>=148 && mac<=161) return 'd';
      if(mac>=162&&mac<=283) return 'c';

      return 'e';
  }

  int CRTHitRecoAlg::MacToRegion(int mac){

      if(mac>=162 && mac<=245) return 38; //top
      if(mac>=246 && mac<=258) return 52; //slope left
      if(mac>=259 && mac<=271) return 56; //slope right
      if(mac>=272 && mac<=277) return 48; //slope front
      if(mac>=278 && mac<=283) return 46; //slope back
      if(            mac<=17 ) return 50; //left
      if(mac>=50  && mac<=67 ) return 50; //left
      if(mac>=18  && mac<=35 ) return 54; //right
      if(mac>=68  && mac<=85 ) return 54; //right
      if(mac>=36  && mac<=42 ) return 44; //front
      if(mac>=86  && mac<=92 ) return 44; //front
      if(mac>=43  && mac<=49 ) return 42; //back
      if(mac>=93  && mac<=99 ) return 42; //back
      if(mac>=148 && mac<=161) return 58; //bottom

      return 0;
  }


  std::string CRTHitRecoAlg::MacToRegionName(int mac){

      if(mac>=162 && mac<=245) return "top";
      if(mac>=246 && mac<=258) return "slope left";
      if(mac>=259 && mac<=271) return "slope right";
      if(mac>=272 && mac<=277) return "slope front";
      if(mac>=278 && mac<=283) return "slope back";
      if(            mac<=17 ) return "left";
      if(mac>=50  && mac<=67 ) return "left";
      if(mac>=18  && mac<=35 ) return "right";
      if(mac>=68  && mac<=85 ) return "right";
      if(mac>=36  && mac<=42 ) return "front";
      if(mac>=86  && mac<=92 ) return "front";
      if(mac>=43  && mac<=49 ) return "back";
      if(mac>=93  && mac<=99 ) return "back";
      if(mac>=148 && mac<=161) return "bottom";

      return "";
  }

  int CRTHitRecoAlg::ChannelToAuxDetSensitiveID(int mac, int chan) {
    char type = MacToType(mac);
    if (type=='d') return chan;
    if (type=='c') return chan/2;
    if (type=='m') return (chan % 10)*2;

    return INT_MAX;
  }

  int CRTHitRecoAlg::MacToAuxDetID(int mac, int chan){
    char type = MacToType(mac);
    if (type == 'e') return INT_MAX;
    if (type == 'c' || type == 'd') return mac;
    if (type == 'm') {
      if (mac>49) mac+=-50;
      if (mac<40) {
          if (            chan<=9 ) return mac*3;
          if (chan>=10 && chan<=19) return mac*3 + 1;
          if (chan>=20 && chan<=29) return mac*3 + 2;
      }
      else if (mac<47) {
          if (            chan<=9 ) return mac*3 - 1;
          if (chan>=10 && chan<=19) return mac*3;
          if (chan>=20 && chan<=29) return mac*3 + 1;
      }
      else {
          if (            chan<=9 ) return mac*3 - 2;
          if (chan>=10 && chan<=19) return mac*3 - 1;
          if (chan>=20 && chan<=29) return mac*3;
      }

    }

    return INT_MAX;
  }

} //end namespace icarus
