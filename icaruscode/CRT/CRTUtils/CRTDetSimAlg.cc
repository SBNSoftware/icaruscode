#ifndef IC_CRTDETSIMALG_CC
#define IC_CRTDETSIMALG_CC

#include "icaruscode/CRT/CRTUtils/CRTDetSimAlg.h"

namespace icarus{
 namespace crt {

    /*struct Tagger {
      char type;
      int modid;
      std::string reg; //crt region where FEB is located
      std::set<int> layerid; //keep track of layers hit accross whole event window
      std::map<int,int> chanlayer; //map chan # to layer
      //std::pair<int,int> macPair; //which two FEBs provided coincidence (applies to m mods only)
      std::vector<icarus::crt::CRTChannelData> data; //time and charge info for each channel > thresh
      std::vector<int> ichan;
    };//Tagger*/

    bool TimeOrderCRTData(icarus::crt::CRTChannelData crtdat1, icarus::crt::CRTChannelData crtdat2) {
        return ( crtdat1.T0() < crtdat2.T0() );
    }//TimeOrderCRTData()

    //constructor 
    CRTDetSimAlg::CRTDetSimAlg(fhicl::ParameterSet const & p, CLHEP::HepRandomEngine& engine) :
        fNsim_m(0), fNsim_d(0), fNsim_c(0), fNchandat_m(0), fNchandat_d(0), fNchandat_c(0),
        fNmissthr_c(0), fNmissthr_d(0), fNmissthr_m(0), fNmiss_strcoin_c(0), fNdual_m(0),
        fRandEngine(engine), fFebMap(CRTCommonUtils::GetFebMap())
    {

        this->reconfigure(p);
        fRegCounts.clear();
	fRegions.clear();
        fTaggers.clear();
    }

    //getting parameter values from FHiCL
    void CRTDetSimAlg::reconfigure(fhicl::ParameterSet const & p) {
       fVerbose = p.get<bool>("Verbose");
       fUltraVerbose = p.get<bool>("UltraVerbose");
       fGlobalT0Offset = p.get<double>("GlobalT0Offset");
       fTDelayNorm = p.get<double>("TDelayNorm");
       fTDelayShift = p.get<double>("TDelayShift");
       fTDelaySigma = p.get<double>("TDelaySigma");
       fTDelayOffset = p.get<double>("TDelayOffset");
       fTDelayRMSGausNorm = p.get<double>("TDelayRMSGausNorm");
       fTDelayRMSGausShift = p.get<double>("TDelayRMSGausShift");
       fTDelayRMSGausSigma = p.get<double>("TDelayRMSGausSigma");
       fTDelayRMSExpNorm = p.get<double>("TDelayRMSExpNorm");
       fTDelayRMSExpShift = p.get<double>("TDelayRMSExpShift");
       fTDelayRMSExpScale = p.get<double>("TDelayRMSExpScale");
       fPropDelay = p.get<double>("PropDelay");
       fPropDelayError = p.get<double>("PropDelayError");
       fTResInterpolator = p.get<double>("TResInterpolator");
       fUseEdep = p.get<bool>("UseEdep");
       fQ0 = p.get<double>("Q0");
       fQPed = p.get<double>("QPed");
       fQSlope = p.get<double>("QSlope");
       fQRMS = p.get<double>("QRMS");
       fQThresholdC = p.get<double>("QThresholdC");
       fQThresholdM = p.get<double>("QThresholdM");
       fQThresholdD = p.get<double>("QThresholdD");
       fStripCoincidenceWindow = p.get<double>("StripCoincidenceWindow");
       fApplyStripCoinC = p.get<bool>("ApplyStripCoincidenceC");
       fApplyCoincidenceC = p.get<bool>("ApplyCoincidenceC");
       fApplyCoincidenceM = p.get<bool>("ApplyCoincidenceM");
       fApplyCoincidenceD = p.get<bool>("ApplyCoincidenceD");
       fLayerCoincidenceWindowC = p.get<double>("LayerCoincidenceWindowC");
       fLayerCoincidenceWindowM = p.get<double>("LayerCoincidenceWindowM");
       fLayerCoincidenceWindowD = p.get<double>("LayerCoincidenceWindowD");
       //fAbsLenEffC = p.get<double>("AbsLenEffC");
       //fAbsLenEffM = p.get<double>("AbsLenEffM");
       //fAbsLenEffD = p.get<double>("AbsLenEffD");
       fDeadTime = p.get<double>("DeadTime");
       fBiasTime = p.get<double>("BiasTime");
    }//CRTDetSim::reconfigure()

    //-----------------------------------------------------------------------------------------------------------
    vector<pair<CRTData,vector<int>>> CRTDetSimAlg::CreateData() //vector<art::Ptr<sim::AuxDetSimChannel>> channels)
    {
        if(fTaggers.size()==0)
            throw cet::exception("CRTDetSimAlg") << "CreateData() called with empty taggers map!";

        vector<pair<CRTData, vector<int>>> dataCol;
	int eve=1;

        // Services: Geometry, DetectorClocks, RandomNumberGenerator
        //art::ServiceHandle<geo::Geometry> geoService;
        //art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
        //detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

        // A list of hit taggers, before any coincidence requirement
        //std::map<int, Tagger> taggers;

        /*int nsim_m=0, nsim_d=0, nsim_c=0; //number of strips in each subsystem with deposited energy
        int nchandat_m=0, nchandat_d=0, nchandat_c=0; //number of SiPM channel signals above threshold
        int nmissthr_c = 0, nmissthr_d = 0, nmissthr_m = 0; //number of channel signals below threshold
        int nmiss_strcoin_c = 0; //number of channel signals missed due to no fiber-fiber coincidence in a cern strip
        int ndual_m = 0; //number of energy deposits producing signals above threshold at both ends of a minos strip

        std::map<int,int> regCounts;
        std::set<int> regions;

        // Loop through truth AuxDetChannel objects
        // Each channel contains the module and strip IDs and
        //  collection of position/time/deposited E for all MCParticles
        //  (AuxDetIDE) for all tracks in the event
        for (int ichan=0; ichan<(int)channels.size(); ichan++) {

            const int adid = channels[ichan]->AuxDetID(); //CRT module ID number (from gdml)
            const int adsid = channels[ichan]->AuxDetSensitiveID(); //CRT strip ID number (from gdml)
            const geo::AuxDetGeo& adGeo = geoService->AuxDet(adid); //pointer to module object

            //check stripID is consistent with number of sensitive volumes
            if( (int)adGeo.NSensitiveVolume() < adsid){
                mf::LogError("CRT") << "adsID out of bounds! Skipping..." << "\n"
                          << "   " << adGeo.Name()  << " / modID "   << adid
                          << " / stripID " << adsid << '\n';
                continue;
            }

            const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsid); //pointer to strip object
            const char auxDetType = CRTCommonUtils::GetAuxDetType(adGeo); //CRT module type (c, d, or m)
            if (auxDetType=='e') mf::LogError("CRT") << "COULD NOT GET AD TYPE!" << '\n';
            const std::string region = CRTCommonUtils::GetAuxDetRegion(adGeo); //CRT region

            int layid = INT_MAX; //set to 0 or 1 if layerid determined
            int mac5=INT_MAX, mac5dual=INT_MAX; //front-end board ID, dual for MINOS modules (not cut)

            // Find the path to the strip geo node, to locate it in the hierarchy
            std::set<std::string> volNames = { adsGeo.TotalVolume()->GetName() };
            std::vector<std::vector<TGeoNode const*> > paths = geoService->FindAllVolumePaths(volNames);

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

            // Module position in parent (tagger) frame
            double modulePosMother[3]; //position in CRT region volume
            nodeModule->LocalToMaster(origin, modulePosMother);

            // strip position in module frame
            double stripPosMother[3];
            double stripPosModule[3];
            nodeStrip->LocalToMaster(origin, stripPosMother);
            nodeInner->LocalToMaster(stripPosMother,stripPosModule);

            // Determine layid
            //  for C modules, two diff. lay thick
            //    1cm for top (y>0) and 1.5cm for bottom (y<0)
            if (auxDetType == 'c' || auxDetType == 'd')
                layid = (stripPosModule[1] > 0);

            if (auxDetType == 'm') {
              //lateral stacks (6 total, 3 per side)
              if ( region.find("West")!=std::string::npos || region.find("East")!=std::string::npos ) {
                  layid = ( modulePosMother[0]>0 );
              }
              //longitudinal walls
              if ( region=="South" || region=="North" ) {
                  layid = ( modulePosMother[2]> 0 );
              }
            }

            if(layid==INT_MAX) mf::LogError("CRT") << "layid NOT SET!!!" << '\n'
                                       << "   ADType: " << auxDetType << '\n'
                                       << "   ADRegion: " << region << '\n';

            // Simulate the CRT response for each hit in this strip
            for (auto ide : channels[ichan]->AuxDetIDEs()) {
              // count true number of energy deposits for each strip
              if (auxDetType=='c') nsim_c++;
              if (auxDetType=='d') nsim_d++;
              if (auxDetType=='m') nsim_m++;

              //track ID of MC paritcle depositing energy for truth matching
              std::vector<int> trkid;
              trkid.push_back(ide.trackID);

              // What is the distance from the hit (centroid of the entry
              // and exit points) to the readout end?
              double x = (ide.entryX + ide.exitX) / 2;
              double y = (ide.entryY + ide.exitY) / 2;
              double z = (ide.entryZ + ide.exitZ) / 2;
              double world[3] = {x, y, z};
              double svHitPosLocal[3];
              double modHitPosLocal[3];
              adsGeo.WorldToLocal(world, svHitPosLocal); //position in strip frame  (origin at center)
              adGeo.WorldToLocal(world, modHitPosLocal); //position in module frame (origin at center)

              //check hit point is contained within the strip according to geometry
              if ( abs(svHitPosLocal[0])>adsGeo.HalfWidth1()+0.001 ||
                   abs(svHitPosLocal[1])>adsGeo.HalfHeight()+0.001 ||
                   abs(svHitPosLocal[2])>adsGeo.HalfLength()+0.001)
                 mf::LogWarning("CRT") << "HIT POINT OUTSIDE OF SENSITIVE VOLUME!" << '\n'
                                    << "  AD: " << adid << " , ADS: " << adsid << '\n'
                                    << "  Local position (x,y,z): ( " << svHitPosLocal[0]
                                    << " , " << svHitPosLocal[1] << " , " << svHitPosLocal[2] << " )" << '\n';

              // The expected number of PE, using a quadratic model for the distance
              // dependence, and scaling linearly with deposited energy.
              double qr = fUseEdep ? ide.energyDeposited / fQ0 : 1.0;
              if (auxDetType == 'c'&& layid==0) qr *= 1.5; //c bottom layer strips 50% thicker

              //longitudinal distance (m) along the strip for fiber atten. calculation
              //assuming SiPM is on +z end (also -z for m modules)
              double distToReadout = abs( adsGeo.HalfLength() - svHitPosLocal[2])*0.01;
              double distToReadout2 = abs(-adsGeo.HalfLength() - svHitPosLocal[2])*0.01;

              double npeExpected = GetLongAtten(distToReadout) * qr;
              double npeExpected2 = GetLongAtten(distToReadout2) * qr;

              //Attenuation factor for transverse propegation in the bulk (c modules only)
              double abs0=1.0, abs1=1.0;
              if(auxDetType=='c'){
                  std::pair<double,double> tmp = GetTransAtten(svHitPosLocal[0]);
                  abs0 = tmp.first;
                  abs1 = tmp.second;
              }

              //most probable # photons arriving at SiPM
              double npeExp0 = npeExpected * abs0;;
              double npeExp1 = npeExpected * abs1;
              double npeExp0Dual = npeExpected2 * abs0;

              //sanity check on simulated light output
              if (npeExp0<0||npeExp1<0||npeExp0Dual<0) mf::LogError("CRT") << "NEGATIVE PE!!!!!" << '\n';

              // Observed PE (Poisson-fluctuated)
              int npe0 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0);
              int npe1 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp1);
              int npe0Dual = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0Dual);

              // Time relative to trigger [ns], accounting for propagation delay and 'walk'
              // for the fixed-threshold discriminator
              double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
              double t0 = \
                GetChannelTriggerTicks(trigClock, tTrue, npe0, distToReadout*100);
              double t1 = \
                GetChannelTriggerTicks(trigClock, tTrue, npe1, distToReadout*100);
              double t0Dual = \
                GetChannelTriggerTicks(trigClock, tTrue, npe0Dual, distToReadout2*100);

              // Time relative to PPS: Random for now! (FIXME)
              int ppsTicks = \
                CLHEP::RandFlat::shootInt(&fRandEngine, trigClock.Frequency() * 1e6);

              // SiPM and ADC response: Npe to ADC counts
              int q0 = \
                CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
              int q1 = \
                CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));
              int q0Dual = \
                CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0Dual, fQRMS * sqrt(npe0Dual));

              if (q0<0||q1<0||q0Dual<0) mf::LogError("CRT") << "NEGATIVE ADC!!!!!" << '\n';

              // Adjacent channels on a strip are numbered sequentially.
              //
              // In the AuxDetChannelMapAlg methods, channels are identified by an
              // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
              // module, and a channel number from 0 to 32.

              int channel0ID=0, channel1ID=0;

              switch (auxDetType){
                  case 'c' :
                      mac5 =fFebMap[adid][0].first;
                      channel0ID = 2 * adsid + 0;
                      channel1ID = 2 * adsid + 1;
                      if(mac5<107||mac5>230)
                          std::cout << "WARNING: mac5 out of bounds for c-type!" << std::endl;
                      if(channel0ID<0 || channel0ID > 31 || channel1ID<0 || channel1ID>31)
                          std::cout << "WARNING: channel out of bounds for c-type!" << std::endl;
                      break;
                  case 'd' :
                      mac5 = fFebMap[adid][0].first;
                      channel0ID = adsid;
                      if(mac5<93||mac5>106)
                          std::cout << "WARNING: mac5 out of bounds for d-type!" << std::endl;
                      if(channel0ID<0 || channel0ID > 63)
                          std::cout << "WARNING: channel out of bounds for d-type!" << std::endl;
                      break;
                  case 'm' :
                      mac5 = fFebMap[adid][0].first;
                      channel0ID = adsid/2 + 10*(fFebMap[adid][0].second-1);
                      if(mac5<1||mac5>92)
                          std::cout << "WARNING: mac5 out of bounds for m-type!" << std::endl;
                      if(channel0ID<0 || channel0ID > 31)
                          std::cout << "WARNING: channel out of bounds for m-type!" << std::endl;
                      if (fFebMap[adid].size()==2)  {
                          mac5dual = fFebMap[adid][1].first;
                          if(mac5dual<1||mac5dual>92)
                              std::cout << "WARNING: mac5dual out of bounds for m-type!" << std::endl;
                      }
                      break;

              }

              if (mac5==INT_MAX) std::cout << "mac addrs not set!" << std::endl; //'\n';

              // Apply ADC threshold and strip-level coincidence (both fibers fire)
              //if (auxDetType=='c' && q0 > fQThresholdC && q1 > fQThresholdC && util::absDiff(t0, t1) < fStripCoincidenceWindow) {
              if (auxDetType=='c') {
                  if ((fApplyStripCoinC && q0>fQThresholdC && q1>fQThresholdC
                      && util::absDiff(t0,t1)<fStripCoincidenceWindow)||
                      (!fApplyStripCoinC && (q0>fQThresholdC || q1>fQThresholdC)) )
                  {
                      Tagger& tagger = taggers[mac5];
                      tagger.layerid.insert(layid);
                      tagger.chanlayer[channel0ID] = layid;
                      tagger.chanlayer[channel1ID] = layid;
                      tagger.reg = region;
                      tagger.type = 'c';
                      tagger.modid = adid;
                      if (q0>fQThresholdC) {
                          tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
                          tagger.ichan.push_back(ichan);
                          nchandat_c++;
                      }
                      else nmissthr_c++;
                      if (q1>fQThresholdC) {
                          tagger.data.push_back(icarus::crt::CRTChannelData(channel1ID,t1,ppsTicks,q1,trkid));
                          tagger.ichan.push_back(ichan);
                          nchandat_c++;
                      }
                      else nmissthr_c++;
                      //nchandat_c+=2;
                  }
              }//if fiber-fiber coincidence

              if (auxDetType=='d' && q0 > fQThresholdD) {
                      Tagger& tagger = taggers[mac5];
                      tagger.layerid.insert(layid);
                      tagger.chanlayer[channel0ID] = layid;
                      tagger.reg = region;
                      tagger.type = 'd';
                      tagger.modid = adid;
                      tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
                      tagger.ichan.push_back(ichan);
                      nchandat_d++;
              }//if one strip above threshold

              if (auxDetType=='m') {
                      if(q0 > fQThresholdM) {
                        if(mac5>300) std::cout << "WARNING: filling tagger with bad mac5!" << std::endl;
                        Tagger& tagger = taggers[mac5];
                        tagger.layerid.insert(layid);
                        tagger.chanlayer[channel0ID] = layid;
                        tagger.reg = region;
                        tagger.type = 'm';
                        tagger.modid = adid;
                        tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
                        tagger.ichan.push_back(ichan);
                        nchandat_m++;
                      }
                      if(q0Dual > fQThresholdM && fFebMap[adid].size()==2) {
                        if(mac5dual>300) std::cout << "WARNING: filling tagger with bad mac5dual!" << std::endl;
                        Tagger& tagger = taggers[mac5dual];
                        tagger.layerid.insert(layid);
                        tagger.chanlayer[channel0ID] = layid;
                        tagger.reg = region;
                        tagger.type = 'm';
                        tagger.modid = adid;
                        tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0Dual,ppsTicks,q0Dual,trkid));
                        tagger.ichan.push_back(ichan);
                        nchandat_m++;
                      }
                      if(q0 > fQThresholdM && q0Dual > fQThresholdM) ndual_m++;
              }//if one strip above threshold at either end

              //counting losses
              if (auxDetType == 'c') {
                  //if (q0 < fQThresholdC || q1 < fQThresholdC) nmissthr_c++;
                  //if ( util::absDiff(t0,t1) >= fStripCoincidenceWindow ) nmiss_strcoin_c++;
                  if ( fApplyStripCoinC && util::absDiff(t0,t1) >= fStripCoincidenceWindow ) nmiss_strcoin_c++;
              }
              if (auxDetType == 'd' && q0 < fQThresholdD) nmissthr_d++;
              if (auxDetType == 'm') {
                  if( q0 < fQThresholdM) nmissthr_m++;
                  if( q0Dual < fQThresholdM && fFebMap[adid].size()==2) nmissthr_m++;
              }

              //print detsim info (if enabled)
              if (fUltraVerbose&&
                 ( //(auxDetType=='c' && q0>fQThresholdC && q1>fQThresholdC) ||
                   (auxDetType=='c' && q0>fQThresholdC && q1>fQThresholdC && util::absDiff(t0, t1) < fStripCoincidenceWindow) ||
                   (auxDetType=='d' && q0>fQThresholdD ) ||
                   (auxDetType=='m' && (q0>fQThresholdM || q0Dual>fQThresholdM)) ))
                std::cout << '\n'
                << "CRT HIT VOL " << (adGeo.TotalVolume())->GetName() << "\n"
                << "CRT HIT SENSITIVE VOL " << (adsGeo.TotalVolume())->GetName() << "\n"
                << "CRT HIT AuxDetID " <<  channels[ichan]->AuxDetID() << " / AuxDetSensitiveID " << channels[ichan]->AuxDetSensitiveID() << "\n"
                << "CRT module type: " << auxDetType << " , CRT region: " << region << '\n'
                << "CRT channel: " << channel0ID << " , mac5: " << mac5 << '\n'
                << "CRT HIT POS (world coords) " << x << " " << y << " " << z << "\n"
                << "CRT STRIP POS (module coords) " << svHitPosLocal[0] << " " << svHitPosLocal[1] << " " << svHitPosLocal[2] << "\n"
                << "CRT MODULE POS (region coords) " << modHitPosLocal[0] << " " << modHitPosLocal[1] << " "<< modHitPosLocal[2] << " " << "\n"
                << "CRT layer ID: " << layid << "\n"
                << "CRT distToReadout: " << distToReadout << ", distToReadout2: " << distToReadout2 << ", qr = " << qr << '\n'
                << "CRT abs0: " << abs0 << " , abs1: " << abs1 << '\n'
                << "CRT npeExpected: " << npeExpected << " , npeExpected2: " << npeExpected2 << '\n'
                //<< "CRT npeExp0: " << npeExp0 << " , npeExp1: " << npeExp1 << " , npeExp0Dual: " << npeExp0Dual << '\n'
                << "CRT npeSiPM0: " << npe0 << " , npeSiPM1: " << npe1 << " , npeSiPM0Dual: " << npe0Dual << '\n'
                << "CRT charge q0: " << q0 << ", q1: " << q1 << '\n'
                //<< "CRT timing: tTrue: " << tTrue << ", t0: " << t0 << ", t1: " << t1 << ", dt: " << util::absDiff(t0,t1) << '\n'
                << "CRT timing: tTrue: " << tTrue << ", t0: " << t0 << ", t1: " << t1 << ", |t0-t1|: " << util::absDiff(t0,t1) << '\n'
                //<< " recoT-trueT = " << t0-tTrue << std::endl; 
                << " recoT0-trueT = " << t0-tTrue << ", recoT1-trueT = " << t1-tTrue << ", recoT0Dual-trueT = " << t0Dual-tTrue << ", recoT0-recoT0Dual = " << t0-t0Dual << std::endl;

            }//for AuxDetIDEs 
          }//for AuxDetChannels */

          //----------------------------------------------------------------------------------
          //if(fUltraVerbose) std::cout << "outside of AD loop" << std::endl;

          // Apply coincidence trigger requirement

          int ncombined_c=0, ncombined_m=0, ncombined_d=0; //channel signals close in time, biasing first signal entering into track and hold circuit
          int nmiss_lock_c=0, nmiss_lock_d=0, nmiss_lock_m=0; //channel signals missed due to channel already locked, but before readout (dead time)
          int nmiss_dead_c=0, nmiss_dead_d=0, nmiss_dead_m=0; //track losses from deadtime
          int nmiss_opencoin_c = 0, nmiss_opencoin_d = 0; //track losses where no coincidence possible
          int nmiss_coin_c = 0;
          int nmiss_coin_d = 0;
          int nmiss_coin_m = 0;
          int event = 0;      //FEB entry number (can have multiple per art event)
          int nhit_m=0, nhit_c=0, nhit_d=0; //channels with signals contained in readout of FEB
          int neve_m=0, neve_c=0, neve_d=0; //no. readouts of FEBs for each subsystem

          // Loop over all FEBs (key for taggers) with a hit and check coincidence requirement.
          // For each FEB, find channel providing trigger and determine if
          //  other hits are in concidence with the trigger (keep) 
          //  or if hits occur during R/O (dead time) (lost)
          //  or if hits are part of a different event (keep for now)
          // First apply dead time correction, biasing effect if configured to do so.
          // Front-end logic: For CERN or DC modules require at least one hit in each X-X layer.
          if (fUltraVerbose) std::cout << '\n' << "about to loop over taggers (size " << fTaggers.size() << " )" << std::endl;

          if(fTaggers.find(2147483647)!=fTaggers.end()) std::cout << "WARNING FOUND BAD MAC5 BEFORE ENTERING TAGGER LOOP!" << std::endl;

          for (auto trg : fTaggers) {

              event = 0;
              icarus::crt::CRTChannelData *chanTrigData(nullptr);
              icarus::crt::CRTChannelData *chanTmpData(nullptr);
              std::set<int> trackNHold = {}; //set of channels close in time to triggered readout above threshold
              std::set<int> layerNHold = {}; //layers with channels above threshold (used for checking layer-layer coincidence)
              if(trg.first==INT_MAX) std::cout << "WARNING: bad mac5 found!" << std::endl;
              std::pair<int,int> macPair = std::make_pair(trg.first,trg.first); //FEB-FEB validation pair
              std::pair<int,int> tpair;  //pair of channels making the coincidence
              bool minosPairFound = false, istrig=false;
              std::vector<icarus::crt::CRTChannelData> passingData; //data to be included in "readout" of FEB
              std::vector<int> passingIchan;
              double ttrig=0.0, ttmp=0.0;  //time stamps on trigger channel, channel considered as part of readout
              int ichantmp =0;
              size_t trigIndex = 0;

              //check "open" coincidence (just check if coincdence possible w/hit in both layers) 
              if (trg.second.type=='c' && fApplyCoincidenceC && trg.second.layerid.size()<2) {
                  nmiss_opencoin_c++;
                  continue;
              }
              if (trg.second.type=='d' && fApplyCoincidenceD && trg.second.layerid.size()<2) {
                  nmiss_opencoin_d++;
                  continue;
              }

              //time order ChannalData objects in this FEB by T0
              std::sort((trg.second.data).begin(),(trg.second.data).end(),TimeOrderCRTData);

              if (fUltraVerbose) std::cout << "processing data for FEB " << trg.first << " with "
                                      << trg.second.data.size() << " entries..." << '\n'
                                      << "    type: " <<  trg.second.type << '\n'
                                      << "    region: " <<  trg.second.reg << '\n'
                                      << "    layerID: " << *trg.second.layerid.begin() << '\n' << std::endl;

              //outer (primary) loop over all data products for this FEB
              for ( size_t i=0; i< trg.second.data.size(); i++ ) {

                //get data for earliest entry
                if(i==0) {
                  chanTrigData = &(trg.second.data[0]);
                  ttrig = chanTrigData->T0(); //time stamp [ns]
                  trackNHold.insert(chanTrigData->Channel());
                  layerNHold.insert(trg.second.chanlayer[chanTrigData->Channel()]);
                  //tpair = std::make_pair(chanTrigData->Channel(),chanTrigData->Channel());
                  passingData.push_back(*chanTrigData);
                  passingIchan.push_back(trg.second.ichan[0]);
                  ttmp = ttrig; 
                  ichantmp = passingIchan.back();
                  if(trg.second.type!='m')
                      continue;
                }

                else {
                  chanTmpData = &(trg.second.data[i]);
                  ttmp = chanTmpData->T0();
                  ichantmp = trg.second.ichan[i];
                }

                if(i!=0 && i==trigIndex) std::cout << "WARNING: tmp index same as trig index!" << std::endl;
                //check that time sorting works
                if ( ttmp < ttrig ) mf::LogError("CRT") << "SORTING OF DATA PRODUCTS FAILED!!!"<< "\n";

                //for C and D modules only and coin. enabled, if assumed trigger channel has no coincidence
                // set trigger channel to tmp channel and try again
                if ( layerNHold.size()==1 &&
                     ( (trg.second.type=='c' && fApplyCoincidenceC && ttmp-ttrig>fLayerCoincidenceWindowC) ||
                       (trg.second.type=='d' && fApplyCoincidenceD && ttmp-ttrig>fLayerCoincidenceWindowD )))
                {
                     trigIndex++;
                     chanTrigData = &(trg.second.data[trigIndex]);
                     i = trigIndex; //+1
                     ttrig = chanTrigData->T0();
                     trackNHold.clear();
                     layerNHold.clear();
                     passingData.clear();
                     passingIchan.clear();
                     trackNHold.insert(chanTrigData->Channel());
                     layerNHold.insert(trg.second.chanlayer[chanTrigData->Channel()]);
                     passingData.push_back(*chanTrigData);
                     passingIchan.push_back(trg.second.ichan[trigIndex]);
                     if(trg.second.type=='c') nmiss_coin_c++;
                     if(trg.second.type=='d') nmiss_coin_d++;
                     continue;
                }

                if(fUltraVerbose) std::cout << "checking for coincidence..." << std::endl;
                //check if coincidence condtion met
                //for c and d modules, just need time stamps within tagger obj
                //for m modules, need to check coincidence with other tagger objs
                if (trg.second.type=='m' && !minosPairFound && fApplyCoincidenceM) {
                    for (auto trg2 :fTaggers) {

                        if(trg2.first==INT_MAX) std::cout << "WARNING: bad mac5 found in inner tagger loop!" << std::endl;

                        if( trg2.second.type!='m' || //is other mod 'm' type
                          trg.second.modid == trg2.second.modid || //other mod not same as this one
                          trg.second.reg != trg2.second.reg || //other mod is in same region
                          *trg2.second.layerid.begin() == *trg.second.layerid.begin()) //modules are in adjacent layers
                            continue;

                        //time sort all data for this FEB (all times for this event)
                        std::sort((trg2.second.data).begin(),(trg2.second.data).end(),TimeOrderCRTData);

                        //find entry within coincidence window starting with this FEB's
                        //triggering channel in coincidence candidate's FEB
                        for ( size_t j=0; j< trg2.second.data.size(); j++ ) {
                            double t2tmp = trg2.second.data[j].T0(); //in us
                            if ( util::absDiff(t2tmp,ttrig) < fLayerCoincidenceWindowM) {
                                minosPairFound = true;
                                macPair = std::make_pair(trg.first,trg2.first);
                                tpair = std::make_pair(chanTrigData->Channel(),trg2.second.data[j].Channel());
                                break;
                            }
                        }
                        //we found a valid pair so move on to next step
                        if (minosPairFound) {
                            if(fUltraVerbose) std::cout << "MINOS pair found! (" << macPair.first
                                                   << "," << macPair.second << ")" << std::endl;
                            break;
                        }
                    }//inner loop over febs (taggers)

                    //if no coincidence pairs found, reinitialize and move to next FEB
                    if(!minosPairFound) {
                        if(fVerbose) std::cout << "MINOS pair NOT found! Skipping to next FEB..." << std::endl;
                        if(trg.second.data.size()==1) continue;
                        trigIndex++;
                        chanTrigData = &(trg.second.data[trigIndex]);
                        i = trigIndex;//+1;
                        ttrig = chanTrigData->T0();
                        trackNHold.clear();
                        layerNHold.clear();
                        passingData.clear();
                        tpair = {};
                        trackNHold.insert(chanTrigData->Channel());
                        layerNHold.insert(trg.second.chanlayer[chanTrigData->Channel()]);
                        passingData.push_back(*chanTrigData);
			passingIchan.push_back(trg.second.ichan[trigIndex]);
                        nmiss_coin_m++;
                        continue;
                    }
                    else
                        istrig = true;
                }//if minos module and no pair yet found

                if(fUltraVerbose) std::cout << "done checking coinceidence...moving on to latency effects..." << std::endl;

                if(!minosPairFound && trg.second.type=='m') std::cout << "WARNING: !minosPairFound...should not see this message!" << std::endl;

                int adctmp = 0;
                std::vector<int> combined_trackids;
                //currently assuming layer coincidence window is same as track and hold window (FIX ME!)
                if (i>0 && ((trg.second.type=='c' && ttmp < ttrig + fLayerCoincidenceWindowC) || 
                    (trg.second.type=='d' && ttmp < ttrig + fLayerCoincidenceWindowD) ||
                    (trg.second.type=='m' && ttmp < ttrig + fLayerCoincidenceWindowM)) )
                {

                    //if channel not locked
                    if ((trackNHold.insert(chanTmpData->Channel())).second) {

                        //channel added to vector of ChannelData for current FEB readout
                        passingData.push_back(*chanTmpData);
			passingIchan.push_back(ichantmp);

                        //if not m module, check to see if strip is first in time in adjacent layer w.r.t. trigger strip
                        if (layerNHold.insert(trg.second.chanlayer[chanTmpData->Channel()]).second
                           && trg.second.type != 'm')
                        {
                            //flagging strips which produce triggering condition (layer-layer coincidence)
                            tpair =std::make_pair(chanTrigData->Channel(),chanTmpData->Channel());
                            istrig = true;

                            //make sure trigger pair strips come from adjacent layers, not same layer
                            if (trg.second.type=='c'&&
                              ((tpair.first<16&&tpair.second<16)||(tpair.first>15&&tpair.second>15)) )
                                mf::LogError("CRT")<< "incorrect CERN trigger pair!!!" << '\n'
                                                   << "  " << tpair.first << ", " << tpair.second << "\n";

                            if (trg.second.type=='d'&&
                              ((tpair.first<32&&tpair.second<32)||(tpair.first>31&&tpair.second>31)) )
                                mf::LogError("CRT")<< "incorrect DC trigger pair!!!" << '\n'
                                                   << "  " << tpair.first << ", " << tpair.second << "\n";
                        }
                    } //end if channel not locked

                    //check for signal biasing
                    else if (ttmp < ttrig + fBiasTime) {
                        combined_trackids.clear();
                        for(size_t dat=0; dat<passingData.size(); dat++) {
                            if(passingData[dat].Channel()!=chanTmpData->Channel())
                                continue;
                            adctmp = passingData[dat].ADC();
                            adctmp += chanTmpData->ADC();
                            passingData[dat].SetADC(adctmp);
		            passingIchan.push_back(ichantmp);

                            for ( auto const& ids : passingData[dat].TrackID())
                                combined_trackids.push_back(ids);

                            combined_trackids.push_back(chanTmpData->TrackID()[0]);
                            passingData[dat].SetTrackID(combined_trackids);
                            break;

                        }

                        switch (trg.second.type) {
                            case 'c' : ncombined_c++; break;
                            case 'd' : ncombined_d++; break;
                            case 'm' : ncombined_m++; break;
                        }
                    } //channel is locked but hit close in time to bias pulse height

                    else switch (trg.second.type) {
                        case 'c' : nmiss_lock_c++; break;
                        case 'd' : nmiss_lock_d++; break;
                        case 'm' : nmiss_lock_m++; break;
                    } //data is discarded (not writeable)

                }//if hits inside track and hold window

                else if ( i>0 && ttmp <= ttrig + fDeadTime ) {
                    switch (trg.second.type) {
                        case 'c' : nmiss_dead_c++; break;
                        case 'd' : nmiss_dead_d++; break;
                        case 'm' : nmiss_dead_m++; break;
                    }
                } // hits occuring during digitization lost (dead time)

                //"read out" data for this event, first hit after dead time as next trigger channel
                else if ( ttmp > ttrig + fDeadTime) {

                  if(istrig)
                    {int regnum = CRTCommonUtils::GetAuxDetRegionNum(trg.second.reg);
                    if( (fRegions.insert(regnum)).second) fRegCounts[regnum] = 1;
                    else fRegCounts[regnum]++;

                    if (fUltraVerbose) {
                        std::cout << "creating CRTData product just after deadtime" << '\n'
                                  << "  event:          " << eve << '\n'
                                  << "  mac5:           " << trg.first << '\n'
                                  << "  FEB entry:      " << event << '\n'
                                  << "  trig time:      " << ttrig << '\n'
                                  << "  trig channel:   " << chanTrigData->Channel() << '\n'
                                  << "  trigger pair:    (" << tpair.first << "," << tpair.second << ")" << '\n'
                                  << "  mac pair:        (" << macPair.first << "," << macPair.second << ")" << '\n'
                                  << "  passing data:   " << std::endl;
                        for(size_t i=0; i<passingData.size(); i++) {
                            std::cout
                                  << "     index:          " << i << '\n'
                                  << "     channel:        " << passingData[i].Channel() << '\n'
                                  << "     t0:             " << passingData[i].T0() << '\n'
                                  << "     t1:             " << passingData[i].T1() << '\n'
                                  << "     adc:            " << passingData[i].ADC() << '\n'
                                  << "     trackIDs:       " << std::endl;
                                  for(auto const& id:passingData[i].TrackID())
                                  std::cout
                                  << "         trackID:         " << id << std::endl;
                        }
                    }

                    if(tpair.first!=chanTrigData->Channel())
                        std::cout << "trigChan mismatch with trigPair!" << '\n'
                                  << "  trigChan: " << chanTrigData->Channel() << ", trigPair[0]: " << tpair.first << std::endl;

                    dataCol.push_back(std::make_pair(
                      CRTData(eve, trg.first,event,ttrig,chanTrigData->Channel(),tpair,macPair,passingData),
                      passingIchan) );

                    if (fUltraVerbose) std::cout << " ...success!" << std::endl;
                    event++;
                    if (trg.second.type=='c') {neve_c++; nhit_c+=passingData.size(); }
                    if (trg.second.type=='d') {neve_d++; nhit_d+=passingData.size(); }
                    if (trg.second.type=='m') {neve_m++; nhit_m+=passingData.size(); } }
                    trigIndex = i;
                    ttrig = ttmp;
                    chanTrigData = chanTmpData;
                    passingData.clear();
                    trackNHold.clear();
                    layerNHold.clear();
                    passingData.push_back(*chanTrigData);
                    trackNHold.insert(chanTrigData->Channel());
                    layerNHold.insert(trg.second.chanlayer[chanTmpData->Channel()]);
                    tpair = {};
                    minosPairFound = false;
                    istrig = false;
                }

                if (!(ttmp > ttrig + fDeadTime) && i==trg.second.data.size()-1 && istrig) {

                    if (fUltraVerbose) {
                        std::cout << "creating CRTData product at end of FEB events..." << '\n'
                                  << "  event:          " << eve << '\n'
                                  << "  mac5:           " << trg.first << '\n'
                                  << "  FEB entry:      " << event << '\n'
                                  << "  trig time:      " << ttrig << '\n'
                                  << "  trig channel:   " << chanTrigData->Channel() << '\n'
                                  << "  trigger pair:    (" << tpair.first << "," << tpair.second << ")" << '\n'
                                  << "  mac pair:        (" << macPair.first << "," << macPair.second << ")" << '\n'
                                  << "  passing data:   " << std::endl;
                        for(size_t i=0; i<passingData.size(); i++) {
                            std::cout
                                  << "     index:          " << i << '\n'
                                  << "     channel:        " << passingData[i].Channel() << '\n'
                                  << "     t0:             " << passingData[i].T0() << '\n'
                                  << "     t1:             " << passingData[i].T1() << '\n'
                                  << "     adc:            " << passingData[i].ADC() << '\n'
                                  << "     trackIDs:       " << std::endl;
                                  for(auto const& id:passingData[i].TrackID())
                                  std::cout
                                  << "         trackID:         " << id << std::endl;
                        }
                    }

                    if(tpair.first!=chanTrigData->Channel())
                        std::cout << "trigChan mismatch with trigPair!" << '\n'
                                  << "  trigChan: " << chanTrigData->Channel() << ", trigPair[0]: " << tpair.first << std::endl;

                    int regnum = CRTCommonUtils::GetAuxDetRegionNum(trg.second.reg);
                    if( (fRegions.insert(regnum)).second) fRegCounts[regnum] = 1;
                    else fRegCounts[regnum]++;

                    dataCol.push_back( std::make_pair(
                      CRTData(eve,trg.first,event,ttrig,chanTrigData->Channel(),tpair,macPair,passingData),
                      passingIchan) );
                    if (fUltraVerbose) std::cout << " ...success!" << std::endl;
                    event++;
                    if (trg.second.type=='c') {neve_c++; nhit_c+=passingData.size(); }
                    if (trg.second.type=='d') {neve_d++; nhit_d+=passingData.size(); }
                    if (trg.second.type=='m') {neve_m++; nhit_m+=passingData.size(); }

                }//if last event and not already written 
              }//for data entries (hits)

              if(fUltraVerbose) std::cout << " outside loop over FEB data entries...moving on to next FEB..." << std::endl;

          } // for taggers

          if (fVerbose) {
            std::cout << "|---------------------- CRTDetSim Summary ----------------------|" << '\n'
             << " - - - EDeps from AuxDetIDE's - - -" << '\n'
             << "  CERN  sim strip hits: " << fNsim_c << '\n'
             << "  MINOS sim strip hits: " << fNsim_m << '\n'
             << "  DC    sim strip hits: " << fNsim_d << '\n' << '\n'
             << " - - - Single Channel Threshold - - -" << '\n'
             << " Pass:" << '\n'
             << "   CERN channel signals  > thresh: " << fNchandat_c << '\n'
             << "   MINOS channel signals > thresh: " << fNchandat_m << '\n'
             << "               both ends > thresh: " << fNdual_m << '\n'
             << "   DC channel signals    > thresh: " << fNchandat_d << '\n'
             << " Lost:" << '\n'
             << "   CERN channel signals  < thresh: " << fNmissthr_c << '\n'
             << "   MINOS channel signals < thresh: " << fNmissthr_m << '\n'
             << "   DC channel signals    < thresh: " << fNmissthr_d << '\n' << '\n'
             << " - - - System Specifc Coincidence Loss- - -" << '\n'
             << "  CERN      fiber-fiber: " << fNmiss_strcoin_c << " (" << fNmiss_strcoin_c*2 << " channels signals lost)" << '\n'
             << "  CERN open layer-layer: " << nmiss_opencoin_c << '\n'
             << "  DC   open layer-layer: " << nmiss_opencoin_d << '\n'
             << "  CERN      layer-layer: " << nmiss_coin_c << " (" << 100.0*nmiss_coin_c/fNchandat_c << "%)" << '\n'
             << "  MINOS     layer-layer: " << nmiss_coin_m << " (" << 100.0*nmiss_coin_m/fNchandat_m << "%)" << '\n'
             << "  DC        layer-layer: " << nmiss_coin_d << " (" << 100.0*nmiss_coin_d/fNchandat_d << "%)" << '\n' << '\n'
             << " - - - Front-End Electronics Effects - - -" << '\n'
             << "  CERN  trackNHold loss: " << nmiss_lock_c << " (" << 100.0*nmiss_lock_c/fNchandat_c << "%)" << '\n'
             << "  MINOS trackNHold loss: " << nmiss_lock_m << " (" << 100.0*nmiss_lock_m/fNchandat_m << "%)" << '\n'
             << "  DC    trackNHold loss: " << nmiss_lock_d << " (" << 100.0*nmiss_lock_d/fNchandat_d << "%)" << '\n'
             << "  CERN   combined:  " << ncombined_c  << " (" << 100.0*ncombined_c/fNchandat_c  << "%)" << '\n'
             << "  MINOS  combined:  " << ncombined_m  << " (" << 100.0*ncombined_m/fNchandat_m  << "%)" << '\n'
             << "  DC     combined:  " << ncombined_d  << " (" << 100.0*ncombined_d/fNchandat_d  << "%)" << '\n'
             << "  CERN    deadTime loss: " << nmiss_dead_c << " (" << 100.0*nmiss_dead_c/fNchandat_c << "%)" << '\n'
             << "  MINOS   deadTime loss: " << nmiss_dead_m << " (" << 100.0*nmiss_dead_m/fNchandat_m << "%)" << '\n'
             << "  DC      deadTime loss: " << nmiss_dead_d << " (" << 100.0*nmiss_dead_d/fNchandat_d << "%)" << '\n' << '\n'
             << " - - - Passing Hits Pushed To Event - - -" << '\n'
             << "  CERN   strip hits remaining: " << nhit_c+ncombined_c << " (" << 100.0*(nhit_c+ncombined_c)/fNchandat_c << "%)" << '\n'
             << "  MINOS  strip hits remaining: " << nhit_m+ncombined_m << " (" << 100.0*(nhit_m+ncombined_m)/fNchandat_m << "%)" << '\n'
             << "  DC     strip hits remaining: " << nhit_d+ncombined_d << " (" << 100.0*(nhit_d+ncombined_d)/fNchandat_d << "%)" << '\n'
             << "  CERN   channel signals generated: " << nhit_c << '\n'
             << "  MINOS  channel signals generated: " << nhit_m << '\n'
             << "  DC     channel signals generated: " << nhit_d << '\n'
             << "  CERN  events: " << neve_c << '\n'
             << "  MINOS events: " << neve_m << '\n'
             << "  DC    events: " << neve_d << '\n'
             << "  Total pushes to event: " << dataCol.size() << std::endl;

            std::map<int,int>::iterator it = fRegCounts.begin();
            std::cout << '\n' << "FEB events per CRT region: " << '\n' << std::endl;

            do {
                std::cout << CRTCommonUtils::GetRegionNameFromNum((*it).first) << ": , events: " << (*it).second << '\n' << std::endl;
                it++;
            }
            while ( it != fRegCounts.end() );

          } //if verbose

	return dataCol;
    }

    //-----------------------------------------------------------------------------
    void CRTDetSimAlg::FillTaggers(const uint32_t adid, const uint32_t adsid, const std::unique_ptr<vector<sim::AuxDetIDE>>& ides, const int nide) {

        art::ServiceHandle<geo::Geometry> geoService;
        art::ServiceHandle<detinfo::DetectorClocksService> detClocks;
        detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

        const geo::AuxDetGeo& adGeo = geoService->AuxDet(adid); //pointer to module object

        //check stripID is consistent with number of sensitive volumes
        if( adGeo.NSensitiveVolume() < adsid){
            mf::LogError("CRT") << "adsID out of bounds! Skipping..." << "\n"
                      << "   " << adGeo.Name()  << " / modID "   << adid
                      << " / stripID " << adsid << '\n';
            //continue;
        }

        const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsid); //pointer to strip object
        const char auxDetType = CRTCommonUtils::GetAuxDetType(adGeo); //CRT module type (c, d, or m)
        if (auxDetType=='e') mf::LogError("CRT") << "COULD NOT GET AD TYPE!" << '\n';
        const std::string region = CRTCommonUtils::GetAuxDetRegion(adGeo); //CRT region

        int layid = INT_MAX; //set to 0 or 1 if layerid determined
        int mac5=INT_MAX, mac5dual=INT_MAX; //front-end board ID, dual for MINOS modules (not cut)

        // Find the path to the strip geo node, to locate it in the hierarchy
        std::set<std::string> volNames = { adsGeo.TotalVolume()->GetName() };
        std::vector<std::vector<TGeoNode const*> > paths = geoService->FindAllVolumePaths(volNames);

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
        // Module position in parent (tagger) frame
        double modulePosMother[3]; //position in CRT region volume
        nodeModule->LocalToMaster(origin, modulePosMother);

        // strip position in module frame
        double stripPosMother[3];
        double stripPosModule[3];
        nodeStrip->LocalToMaster(origin, stripPosMother);
        nodeInner->LocalToMaster(stripPosMother,stripPosModule);

        // Determine layid
        //  for C modules, two diff. lay thick
        //    1cm for top (y>0) and 1.5cm for bottom (y<0)
        if (auxDetType == 'c' || auxDetType == 'd')
            layid = (stripPosModule[1] > 0);

        if (auxDetType == 'm') {
          //lateral stacks (6 total, 3 per side)
          if ( region.find("West")!=std::string::npos || region.find("East")!=std::string::npos ) {
              layid = ( modulePosMother[0]>0 );
          }
          //longitudinal walls
          if ( region=="South" || region=="North" ) {
              layid = ( modulePosMother[2]> 0 );
          }
        }

        if(layid==INT_MAX) mf::LogError("CRT") << "layid NOT SET!!!" << '\n'
                                   << "   ADType: " << auxDetType << '\n'
                                   << "   ADRegion: " << region << '\n';

        // Simulate the CRT response for each hit in this strip
        for (int ide=ides->size()-nide; ide<(int)ides->size(); ide++) {
          // count true number of energy deposits for each strip
          if (auxDetType=='c') fNsim_c++;
          if (auxDetType=='d') fNsim_d++;
          if (auxDetType=='m') fNsim_m++;

          //track ID of MC paritcle depositing energy for truth matching
          std::vector<int> trkid;
          trkid.push_back(ides->at(ide).trackID);

          // What is the distance from the hit (centroid of the entry
          // and exit points) to the readout end?
          double x = (ides->at(ide).entryX + ides->at(ide).exitX) / 2;
          double y = (ides->at(ide).entryY + ides->at(ide).exitY) / 2;
          double z = (ides->at(ide).entryZ + ides->at(ide).exitZ) / 2;
          double world[3] = {x, y, z};
          double svHitPosLocal[3];
          double modHitPosLocal[3];
          adsGeo.WorldToLocal(world, svHitPosLocal); //position in strip frame  (origin at center)
          adGeo.WorldToLocal(world, modHitPosLocal); //position in module frame (origin at center)

          //check hit point is contained within the strip according to geometry
          if ( abs(svHitPosLocal[0])>adsGeo.HalfWidth1()+0.001 ||
               abs(svHitPosLocal[1])>adsGeo.HalfHeight()+0.001 ||
               abs(svHitPosLocal[2])>adsGeo.HalfLength()+0.001)
             mf::LogWarning("CRT") << "HIT POINT OUTSIDE OF SENSITIVE VOLUME!" << '\n'
                                << "  AD: " << adid << " , ADS: " << adsid << '\n'
                                << "  Local position (x,y,z): ( " << svHitPosLocal[0]
                                << " , " << svHitPosLocal[1] << " , " << svHitPosLocal[2] << " )" << '\n';

          // The expected number of PE, using a quadratic model for the distance
          // dependence, and scaling linearly with deposited energy.
          double qr = fUseEdep ? ides->at(ide).energyDeposited / fQ0 : 1.0;
          if (auxDetType == 'c'&& layid==0) qr *= 1.5; //c bottom layer strips 50% thicker

          //longitudinal distance (m) along the strip for fiber atten. calculation
          //assuming SiPM is on +z end (also -z for m modules)
          double distToReadout = abs( adsGeo.HalfLength() - svHitPosLocal[2])*0.01;
          double distToReadout2 = abs(-adsGeo.HalfLength() - svHitPosLocal[2])*0.01;

          double npeExpected = GetLongAtten(distToReadout) * qr;
          double npeExpected2 = GetLongAtten(distToReadout2) * qr;

          //Attenuation factor for transverse propegation in the bulk (c modules only)
          double abs0=1.0, abs1=1.0;
          if(auxDetType=='c'){
              std::pair<double,double> tmp = GetTransAtten(svHitPosLocal[0]);
              abs0 = tmp.first;
              abs1 = tmp.second;
          }

          //most probable # photons arriving at SiPM
          double npeExp0 = npeExpected * abs0;;
          double npeExp1 = npeExpected * abs1;
          double npeExp0Dual = npeExpected2 * abs0;

          //sanity check on simulated light output
          if (npeExp0<0||npeExp1<0||npeExp0Dual<0) mf::LogError("CRT") << "NEGATIVE PE!!!!!" << '\n';

          // Observed PE (Poisson-fluctuated)
          int npe0 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0);
          int npe1 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp1);
          int npe0Dual = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0Dual);

          // Time relative to trigger [ns], accounting for propagation delay and 'walk'
          // for the fixed-threshold discriminator
          double tTrue = (ides->at(ide).entryT + ides->at(ide).exitT) / 2 + fGlobalT0Offset;
          double t0 = \
            GetChannelTriggerTicks(trigClock, tTrue, npe0, distToReadout*100);
          double t1 = \
            GetChannelTriggerTicks(trigClock, tTrue, npe1, distToReadout*100);
          double t0Dual = \
            GetChannelTriggerTicks(trigClock, tTrue, npe0Dual, distToReadout2*100);

          // Time relative to PPS: Random for now! (FIXME)
          int ppsTicks = \
            CLHEP::RandFlat::shootInt(&fRandEngine, trigClock.Frequency() * 1e6);

          // SiPM and ADC response: Npe to ADC counts
          int q0 = \
            CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
          int q1 = \
            CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));
          int q0Dual = \
            CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0Dual, fQRMS * sqrt(npe0Dual));

          if (q0<0||q1<0||q0Dual<0) mf::LogError("CRT") << "NEGATIVE ADC!!!!!" << '\n';

         // Adjacent channels on a strip are numbered sequentially.
          //
          // In the AuxDetChannelMapAlg methods, channels are identified by an
          // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
          // module, and a channel number from 0 to 32.

          int channel0ID=0, channel1ID=0;

          switch (auxDetType){
              case 'c' :
                  mac5 =fFebMap[adid][0].first;
                  channel0ID = 2 * adsid + 0;
                  channel1ID = 2 * adsid + 1;
                  if(mac5<107||mac5>230)
                      std::cout << "WARNING: mac5 out of bounds for c-type!" << std::endl;
                  if(channel0ID<0 || channel0ID > 31 || channel1ID<0 || channel1ID>31)
                      std::cout << "WARNING: channel out of bounds for c-type!" << std::endl;
                  break;
              case 'd' :
                  mac5 = fFebMap[adid][0].first;
                  channel0ID = adsid;
                  if(mac5<93||mac5>106)
                      std::cout << "WARNING: mac5 out of bounds for d-type!" << std::endl;
                  if(channel0ID<0 || channel0ID > 63)
                      std::cout << "WARNING: channel out of bounds for d-type!" << std::endl;
                  break;
              case 'm' :
                  mac5 = fFebMap[adid][0].first;
                  channel0ID = adsid/2 + 10*(fFebMap[adid][0].second-1);
                  if(mac5<1||mac5>92)
                      std::cout << "WARNING: mac5 out of bounds for m-type!" << std::endl;
                  if(channel0ID<0 || channel0ID > 31)
                      std::cout << "WARNING: channel out of bounds for m-type!" << std::endl;
                  if (fFebMap[adid].size()==2)  {
                      mac5dual = fFebMap[adid][1].first;
                      if(mac5dual<1||mac5dual>92)
                          std::cout << "WARNING: mac5dual out of bounds for m-type!" << std::endl;
                  }
                  break;

          }

         if (mac5==INT_MAX) /*mf::LogError("CRT")*/ std::cout << "mac addrs not set!" << std::endl; //'\n';

          // Apply ADC threshold and strip-level coincidence (both fibers fire)
          //if (auxDetType=='c' && q0 > fQThresholdC && q1 > fQThresholdC && util::absDiff(t0, t1) < fStripCoincidenceWindow) {
          if (auxDetType=='c') {
              if ((fApplyStripCoinC && q0>fQThresholdC && q1>fQThresholdC
                  && util::absDiff(t0,t1)<fStripCoincidenceWindow)||
                  (!fApplyStripCoinC && (q0>fQThresholdC || q1>fQThresholdC)) )
              {
                  Tagger& tagger = fTaggers[mac5];
                  tagger.layerid.insert(layid);
                  tagger.chanlayer[channel0ID] = layid;
                  tagger.chanlayer[channel1ID] = layid;
                  tagger.reg = region;
                  tagger.type = 'c';
                  tagger.modid = adid;
                  if (q0>fQThresholdC) {
                      tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
                      tagger.ichan.push_back(ide);
                      fNchandat_c++;
                  }
                  else fNmissthr_c++;
                  if (q1>fQThresholdC) {
                      tagger.data.push_back(icarus::crt::CRTChannelData(channel1ID,t1,ppsTicks,q1,trkid));
                      tagger.ichan.push_back(ide);
                      fNchandat_c++;
                  }
                  else fNmissthr_c++;
                  //nchandat_c+=2;
              }
          }//if fiber-fiber coincidence

          if (auxDetType=='d' && q0 > fQThresholdD) {
                  Tagger& tagger = fTaggers[mac5];
                  tagger.layerid.insert(layid);
                  tagger.chanlayer[channel0ID] = layid;
                  tagger.reg = region;
                  tagger.type = 'd';
                  tagger.modid = adid;
                  tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
                  tagger.ichan.push_back(ide);
                  fNchandat_d++;
          }//if one strip above threshold

          if (auxDetType=='m') {
                  if(q0 > fQThresholdM) {
                    if(mac5>300) std::cout << "WARNING: filling tagger with bad mac5!" << std::endl;
                    Tagger& tagger = fTaggers[mac5];
                    tagger.layerid.insert(layid);
                    tagger.chanlayer[channel0ID] = layid;
                    tagger.reg = region;
                    tagger.type = 'm';
                    tagger.modid = adid;
                    tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0,ppsTicks,q0,trkid));
                    tagger.ichan.push_back(ide);
                    fNchandat_m++;
                  }
                  if(q0Dual > fQThresholdM && fFebMap[adid].size()==2) {
                    if(mac5dual>300) std::cout << "WARNING: filling tagger with bad mac5dual!" << std::endl;
                    Tagger& tagger = fTaggers[mac5dual];
                    tagger.layerid.insert(layid);
                    tagger.chanlayer[channel0ID] = layid;
                    tagger.reg = region;
                    tagger.type = 'm';
                    tagger.modid = adid;
                    tagger.data.push_back(icarus::crt::CRTChannelData(channel0ID,t0Dual,ppsTicks,q0Dual,trkid));
                    tagger.ichan.push_back(ide);
                    fNchandat_m++;
                  }
                  if(q0 > fQThresholdM && q0Dual > fQThresholdM) fNdual_m++;
          }//if one strip above threshold at either end

          //counting losses
          if (auxDetType == 'c') {
              //if (q0 < fQThresholdC || q1 < fQThresholdC) nmissthr_c++;
              //if ( util::absDiff(t0,t1) >= fStripCoincidenceWindow ) nmiss_strcoin_c++;
              if ( fApplyStripCoinC && util::absDiff(t0,t1) >= fStripCoincidenceWindow ) fNmiss_strcoin_c++;
          }
          if (auxDetType == 'd' && q0 < fQThresholdD) fNmissthr_d++;
          if (auxDetType == 'm') {
              if( q0 < fQThresholdM) fNmissthr_m++;
              if( q0Dual < fQThresholdM && fFebMap[adid].size()==2) fNmissthr_m++;
          }

          //print detsim info (if enabled)
          if (fUltraVerbose&&
             ( //(auxDetType=='c' && q0>fQThresholdC && q1>fQThresholdC) ||
               (auxDetType=='c' && q0>fQThresholdC && q1>fQThresholdC && util::absDiff(t0, t1) < fStripCoincidenceWindow) ||
               (auxDetType=='d' && q0>fQThresholdD ) ||
               (auxDetType=='m' && (q0>fQThresholdM || q0Dual>fQThresholdM)) ))
            std::cout << '\n'
            << "CRT HIT VOL " << (adGeo.TotalVolume())->GetName() << "\n"
            << "CRT HIT SENSITIVE VOL " << (adsGeo.TotalVolume())->GetName() << "\n"
            << "CRT HIT AuxDetID " <<  adid << " / AuxDetSensitiveID " << adsid << "\n"
            << "CRT module type: " << auxDetType << " , CRT region: " << region << '\n'
            << "CRT channel: " << channel0ID << " , mac5: " << mac5 << '\n'
            << "CRT HIT POS (world coords) " << x << " " << y << " " << z << "\n"
            << "CRT STRIP POS (module coords) " << svHitPosLocal[0] << " " << svHitPosLocal[1] << " " << svHitPosLocal[2] << "\n"
            << "CRT MODULE POS (region coords) " << modHitPosLocal[0] << " " << modHitPosLocal[1] << " "<< modHitPosLocal[2] << " " << "\n"
            << "CRT layer ID: " << layid << "\n"
            << "CRT distToReadout: " << distToReadout << ", distToReadout2: " << distToReadout2 << ", qr = " << qr << '\n'
            << "CRT abs0: " << abs0 << " , abs1: " << abs1 << '\n'
            << "CRT npeExpected: " << npeExpected << " , npeExpected2: " << npeExpected2 << '\n'
            //<< "CRT npeExp0: " << npeExp0 << " , npeExp1: " << npeExp1 << " , npeExp0Dual: " << npeExp0Dual << '\n'
            << "CRT npeSiPM0: " << npe0 << " , npeSiPM1: " << npe1 << " , npeSiPM0Dual: " << npe0Dual << '\n'
            << "CRT charge q0: " << q0 << ", q1: " << q1 << '\n'
            //<< "CRT timing: tTrue: " << tTrue << ", t0: " << t0 << ", t1: " << t1 << ", dt: " << util::absDiff(t0,t1) << '\n'
            << "CRT timing: tTrue: " << tTrue << ", t0: " << t0 << ", t1: " << t1 << ", |t0-t1|: " << util::absDiff(t0,t1) << '\n'
            //<< " recoT-trueT = " << t0-tTrue << std::endl; 
            << " recoT0-trueT = " << t0-tTrue << ", recoT1-trueT = " << t1-tTrue << ", recoT0Dual-trueT = " << t0Dual-tTrue << ", recoT0-recoT0Dual = " << t0-t0Dual << std::endl;

        }//for AuxDetIDEs 

    }

    //---------------------------------------------------------------
    std::pair<double,double> CRTDetSimAlg::GetTransAtten(const double pos) {
        //attenuation coefficients for each fiber
        double abs0=0.0, abs1=0.0;

        //coefficiencts for transverse attenuation in CERN modules 
        //from early beta-source tranverse scan data
        double at0_c = 0.682976;
        double at1_c = -0.0204477;
        double at2_c = -0.000707564;
        double at3_c = 0.000636617;
        double at4_c = 0.000147957;
        double at5_c = -3.89078e-05;

        double at0_r = 0.139941;
        double at1_r = 0.168238;
        double at2_r = -0.0198199;
        double at3_r = 0.000781752;

        double at0_l = 8.78875;
        double at1_l = 3.54602;
        double at2_l = 0.595592;
        double at3_l = 0.0449169;
        double at4_l = 0.00127892;

        //hit between both fibers 
        if ( abs(pos) <= 5.5 ) {
            abs0 = at5_c*pow(pos,5) + at4_c*pow(pos,4) + at3_c*pow(pos,3)
                    + at2_c*pow(pos,2) + at1_c*pos + at0_c; 
            abs1 = -1*at5_c*pow(pos,5) + at4_c*pow(pos,4) - at3_c*pow(pos,3)
                    + at2_c*pow(pos,2) - at1_c*pos + at0_c;
        }

        //hit to right of both fibers
        if ( pos > 5.5 ) {
            abs0 = at3_r*pow(pos,3) + at2_r*pow(pos,2) + at1_r*pos + at0_r;
            abs1 = at4_l*pow(pos,4) - at3_l*pow(pos,3)
              + at2_l*pow(pos,2) - at1_l*pos + at0_l;
        }

        //hit to left of both fibers
        if ( pos < -5.5 ) {
             abs0 = at4_l*pow(pos,4) + at3_l*pow(pos,3) \
                    + at2_l*pow(pos,2) + at1_l*pos + at0_l;
             abs1 = -1*at3_r*pow(pos,3) + at2_r*pow(pos,2) - at1_r*pos + at0_r;
        }

        return std::make_pair(abs0,abs1);

    }

    //---------------------------------------------------------------
    double CRTDetSimAlg::GetLongAtten(const double dist) {

        double pes=0.0;

        //coefficients for quadratic fit to MINOS test data w/S14
        //obtained for normally incident cosmic muons
        //p2_m * x^2 + p1_m * x + p0_m, x is distance from readout [m]
        double p0_m = 36.5425; //initial light yield (pe) before any attenuation in m scintillator
        double p1_m = -6.3895;
        double p2_m =  0.3742;

        pes = p2_m * pow(dist,2) + p1_m * dist + p0_m;

        return pes;

    }

    //function for simulating time response
    //  takes true hit time, LY(PE) observed, and intitudinal distance from readout
    //  uses 12 FHiCL configurable parameters
    //  returns simulated time in units of clock ticks
    double CRTDetSimAlg::GetChannelTriggerTicks(detinfo::ElecClock& clock,
                                             float t0, float npeMean, float r)
    {
        // Hit timing, with smearing and NPE dependence
        double tDelayMean =
          fTDelayNorm *
            exp(-0.5 * pow((npeMean - fTDelayShift) / fTDelaySigma, 2)) +
          fTDelayOffset;
    
        double tDelayRMS =
          fTDelayRMSGausNorm *
            exp(-pow(npeMean - fTDelayRMSGausShift, 2) / fTDelayRMSGausSigma) +
          fTDelayRMSExpNorm *
            exp(-(npeMean - fTDelayRMSExpShift) / fTDelayRMSExpScale);
    
        double tDelay = CLHEP::RandGauss::shoot(&fRandEngine, tDelayMean, tDelayRMS);
    
        // Time resolution of the interpolator
        tDelay += CLHEP::RandGauss::shoot(&fRandEngine, 0, fTResInterpolator);
    
        // Propagation time
        double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;
    
        double t = t0 + tProp + tDelay;
    
        // Get clock ticks
        clock.SetTime(t / 1e3);  // SetTime takes microseconds
    
        /*if (fUltraVerbose) mf::LogInfo("CRT")
          << "CRT TIMING: t0=" << t0
          << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
          << ", tDelay=" << tDelay << ", tDelay(interp)="
          << tDelay << ", tProp=" << tProp << ", t=" << t << "\n";//, ticks=" << clock.Ticks() << "\n"; 
        */
        return t;//clock.Ticks();
    }//CRTDetSim::GetChannelTriggerTicks()

    void CRTDetSimAlg::ClearTaggers() {
        fTaggers.clear();
    }
    
 }//namespace crt
}//namespace icarus

#endif
