#ifndef IC_CRTDETSIMALG_CC
#define IC_CRTDETSIMALG_CC

#include "icaruscode/CRT/CRTUtils/CRTDetSimAlg.h"

namespace icarus{
 namespace crt {

    bool TimeOrderCRTData(std::pair<ChanData, AuxDetIDE> crtdat1, 
                          std::pair<ChanData, AuxDetIDE> crtdat2) {
        return ( crtdat1.first.ts < crtdat2.first.ts );
    }//TimeOrderCRTData()

    //-------------------------------------------------------------------------------------------
    //constructor 
    CRTDetSimAlg::CRTDetSimAlg(fhicl::ParameterSet const & p, CLHEP::HepRandomEngine& engine) :
        fNsim_m(0), fNsim_d(0), fNsim_c(0), fNchandat_m(0), fNchandat_d(0), fNchandat_c(0),
        fNmissthr_c(0), fNmissthr_d(0), fNmissthr_m(0), fNmiss_strcoin_c(0), fNdual_m(0),
        fHasFilledTaggers(false), fRandEngine(engine)
    {

        this->reconfigure(p);
        fRegCounts.clear();
	fRegions.clear();
        fTaggers.clear();
        fCrtutils = new CRTCommonUtils();
    }
    //----------------------------------------------------------------------
    //getting parameter values from FHiCL (called by constructor)
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
       fQMax = p.get<uint16_t>("QMax");
       fQThresholdC = p.get<uint16_t>("QThresholdC");
       fQThresholdM = p.get<uint16_t>("QThresholdM");
       fQThresholdD = p.get<uint16_t>("QThresholdD");
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
       fUseBirks = p.get<bool>("UseBirks");
       fKbirks   = p.get<double>("Kbirks");
    }//CRTDetSim::reconfigure()

    //-------------------------------------------------------------------------------------------------------------
    // function called after loop over AuxDetSimChannels where FillTaggers was used to perform first detsim step.
    // this function applies trigger logic, deadtime, and close-in-time signal biasing effects. it produces the
    // "triggered data" products which make it into the event
    vector<pair<CRTData,vector<sim::AuxDetIDE>>> CRTDetSimAlg::CreateData()
    {        
        vector<pair<CRTData, vector<AuxDetIDE>>> dataCol;

        if(!fHasFilledTaggers)
            throw cet::exception("CRTDetSimAlg") << "CreateData() called with empty taggers map!";
        if(fTaggers.size()==0)
            return dataCol;

	int eve=1;

        int ncombined_c=0, ncombined_m=0, ncombined_d=0; //channel signals close in time, biasing first signal entering into track and hold circuit
        int nmiss_lock_c=0, nmiss_lock_d=0, nmiss_lock_m=0; //channel signals missed due to channel already locked, but before readout (dead time)
        int nmiss_dead_c=0, nmiss_dead_d=0, nmiss_dead_m=0; //track losses from deadtime
        int nmiss_opencoin_c = 0, nmiss_opencoin_d = 0; //track losses where no coincidence possible
        int nmiss_coin_c = 0;
        int nmiss_coin_d = 0;
        int nmiss_coin_m = 0;
        int event = 0;                    //FEB entry number (can have multiple per art event)
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

        for (auto trg : fTaggers) {
            //if(trg.second.data.size()!=trg.second.ide.size())
            //    std::cout << "WARNING DATA AND INDEX VECTOR SIZE MISMATCH!" << std::endl;

            event = 0;
            ChanData *chanTrigData(nullptr);
            ChanData *chanTmpData(nullptr);
            set<int> trackNHold = {}; //set of channels close in time to triggered readout above threshold
            set<int> layerNHold = {}; //layers with channels above threshold (used for checking layer-layer coincidence)
            if(trg.first==INT_MAX) std::cout << "WARNING: bad mac5 found!" << std::endl;
            bool minosPairFound = false, istrig=false;
            vector<ChanData> passingData; //data to be included in "readout" of FEB
            vector<AuxDetIDE> passingIDE;
            uint64_t ttrig=0.0, ttmp=0.0;  //time stamps on trigger channel, channel considered as part of readout
            AuxDetIDE idetmp;
            size_t trigIndex = 0;
            uint16_t adc[64];

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

            if (fUltraVerbose) std::cout << "processing data for FEB " << (int)trg.first << " with "
                                    << trg.second.data.size() << " entries..." << '\n'
                                    << "    type: " <<  trg.second.type << '\n'
                                    << "    region: " <<  trg.second.reg << '\n'
                                    << "    layerID: " << *trg.second.layerid.begin() << '\n' << std::endl;

            //outer (primary) loop over all data products for this FEB
            for ( size_t i=0; i< trg.second.data.size(); i++ ) {

              //get data for earliest entry
              if(i==0) {
                chanTrigData = &(trg.second.data[0].first);
                ttrig = chanTrigData->ts; //time stamp [ns]
                trackNHold.insert(chanTrigData->channel);
                layerNHold.insert(trg.second.chanlayer[chanTrigData->channel]);
                passingData.push_back(*chanTrigData);
                passingIDE.push_back(trg.second.data[0].second);
                ttmp = ttrig; 
                idetmp = passingIDE.back();
                if(trg.second.type!='m')
                    continue;
              }

              else {
                chanTmpData = &(trg.second.data[i].first);
                ttmp = chanTmpData->ts;
                idetmp = trg.second.data[i].second;
              }

              //for C and D modules only and coin. enabled, if assumed trigger channel has no coincidence
              // set trigger channel to tmp channel and try again
              if ( layerNHold.size()==1 &&
                   ( (trg.second.type=='c' && fApplyCoincidenceC && ttmp-ttrig>fLayerCoincidenceWindowC) ||
                     (trg.second.type=='d' && fApplyCoincidenceD && ttmp-ttrig>fLayerCoincidenceWindowD )))
              {
                   trigIndex++;
                   chanTrigData = &(trg.second.data[trigIndex].first);
                   i = trigIndex; //+1
                   ttrig = chanTrigData->ts;
                   trackNHold.clear();
                   layerNHold.clear();
                   passingData.clear();
                   passingIDE.clear();
                   trackNHold.insert(chanTrigData->channel);
                   layerNHold.insert(trg.second.chanlayer[chanTrigData->channel]);
                   passingData.push_back(*chanTrigData);
                   passingIDE.push_back(trg.second.data[trigIndex].second);
                   if(trg.second.type=='c') nmiss_coin_c++;
                   if(trg.second.type=='d') nmiss_coin_d++;
                   continue;
              }

              if(fUltraVerbose)
                  std::cout << "checking for coincidence..." << std::endl;
              //check if coincidence condtion met
              //for c and d modules, just need time stamps within tagger obj
              //for m modules, need to check coincidence with other tagger objs
              if (trg.second.type=='m' && !minosPairFound && fApplyCoincidenceM) {
                  for (auto trg2 :fTaggers) {

                      if( trg2.second.type!='m' || //is other mod 'm' type
                        trg.second.modid == trg2.second.modid || //other mod not same as this one
                        trg.second.reg != trg2.second.reg || //other mod is in same region
                        *trg2.second.layerid.begin() == *trg.second.layerid.begin()) //modules are in adjacent layers
                          continue;

                      //find entry within coincidence window starting with this FEB's
                      //triggering channel in coincidence candidate's FEB
                      for ( size_t j=0; j< trg2.second.data.size(); j++ ) {
                          uint64_t t2tmp = trg2.second.data[j].first.ts; //in us
                          if ( util::absDiff(t2tmp,ttrig) < fLayerCoincidenceWindowM) {
                              minosPairFound = true;
                              break;
                          }
                      }
                      //we found a valid pair so move on to next step
                      if (minosPairFound)
                          break;
                  }//inner loop over febs (taggers)

                  //if no coincidence pairs found, reinitialize and move to next FEB
                  if(!minosPairFound) {
                      if(fUltraVerbose) 
                          std::cout << "MINOS pair NOT found! Skipping to next FEB..." << std::endl;
                      if(trg.second.data.size()==1) continue;
                      trigIndex++;
                      chanTrigData = &(trg.second.data[trigIndex].first);
                      i = trigIndex;//+1;
                      ttrig = chanTrigData->ts;
                      trackNHold.clear();
                      layerNHold.clear();
                      passingData.clear();
                      passingIDE.clear();
                      trackNHold.insert(chanTrigData->channel);
                      layerNHold.insert(trg.second.chanlayer[chanTrigData->channel]);
                      passingData.push_back(*chanTrigData);
	      	      passingIDE.push_back(trg.second.data[trigIndex].second);
                      nmiss_coin_m++;
                      continue;
                  }
                  else
                      istrig = true;
              }//if minos module and no pair yet found

              if(fUltraVerbose) 
                  std::cout << "done checking coinceidence...moving on to latency effects..." << std::endl;
              if(!minosPairFound && trg.second.type=='m') 
                  std::cout << "WARNING: !minosPairFound...should not see this message!" << std::endl;

              int adctmp = 0;
              //currently assuming layer coincidence window is same as track and hold window (FIX ME!)
              if (i>0 && ((trg.second.type=='c' && ttmp < ttrig + fLayerCoincidenceWindowC) || 
                  (trg.second.type=='d' && ttmp < ttrig + fLayerCoincidenceWindowD) ||
                  (trg.second.type=='m' && ttmp < ttrig + fLayerCoincidenceWindowM)) )
              {

                  //if channel not locked
                  if ((trackNHold.insert(chanTmpData->channel)).second) {

                      //channel added to vector of ChannelData for current FEB readout
                      passingData.push_back(*chanTmpData);
	              passingIDE.push_back(idetmp);

                      //if not m module, check to see if strip is first in time in adjacent layer w.r.t. trigger strip
                      if (layerNHold.insert(trg.second.chanlayer[chanTmpData->channel]).second
                         && trg.second.type != 'm')
                      {
                          //flagging strips which produce triggering condition (layer-layer coincidence)
                          istrig = true;
                      }
                  } //end if channel not locked

                  //check for signal biasing
                  else if (ttmp < ttrig + fBiasTime) {
                      for(size_t dat=0; dat<passingData.size(); dat++) {
                          if(passingData[dat].channel!=chanTmpData->channel)
                              continue;
                          adctmp = passingData[dat].adc;
                          adctmp += chanTmpData->adc;
                          if(adctmp>fQMax) adctmp = fQMax;
                          passingData[dat].adc = adctmp;
                          passingIDE.push_back(idetmp);
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
                  {int regnum = fCrtutils->AuxDetRegionNameToNum(trg.second.reg);
                  if( (fRegions.insert(regnum)).second) fRegCounts[regnum] = 1;
                  else fRegCounts[regnum]++;

                  if (fUltraVerbose) {
                      std::cout << "creating CRTData product just after deadtime" << '\n'
                                << "  event:          " << eve << '\n'
                                << "  mac5:           " << trg.first << '\n'
                                << "  FEB entry:      " << event << '\n'
                                << "  trig time:      " << ttrig << '\n'
                                << "  trig channel:   " << chanTrigData->channel << '\n'
                                << "  passing data:   " << std::endl;
                      for(size_t i=0; i<passingData.size(); i++) {
                          std::cout
                                << "     index:          " << i << '\n'
                                << "     channel:        " << passingData[i].channel << '\n'
                                << "     t0:             " << passingData[i].ts << '\n'
                                << "     adc:            " << passingData[i].adc << std::endl;
                      }
                  }

                  if(passingData.size()>passingIDE.size()) 
                      std::cout << "data/IDE size mismatch!" <<  passingData.size()-passingIDE.size() << std::endl;
                  FillAdcArr(passingData,adc);
                  dataCol.push_back(std::make_pair(
                    FillCRTData(trg.first,event,ttrig,ttrig,adc),
                    passingIDE) );

                  if (fUltraVerbose) std::cout << " ...success!" << std::endl;
                  event++;
                  if (trg.second.type=='c') {neve_c++; nhit_c+=passingData.size(); }
                  if (trg.second.type=='d') {neve_d++; nhit_d+=passingData.size(); }
                  if (trg.second.type=='m') {neve_m++; nhit_m+=passingData.size(); } }
                  trigIndex = i;
                  ttrig = ttmp;
                  chanTrigData = chanTmpData;
                  passingData.clear();
                  passingIDE.clear();
                  trackNHold.clear();
                  layerNHold.clear();
                  passingData.push_back(*chanTrigData);
                  passingIDE.push_back(idetmp);
                  trackNHold.insert(chanTrigData->channel);
                  layerNHold.insert(trg.second.chanlayer[chanTmpData->channel]);
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
                                << "  trig channel:   " << chanTrigData->channel << '\n'
                                << "  passing data:   " << std::endl;
                      for(size_t i=0; i<passingData.size(); i++) {
                          std::cout
                                << "     index:          " << i << '\n'
                                << "     channel:        " << passingData[i].channel << '\n'
                                << "     t0:             " << passingData[i].ts << '\n'
                                << "     adc:            " << passingData[i].adc << std::endl;
                      }
                  }

                  int regnum = fCrtutils->AuxDetRegionNameToNum(trg.second.reg);
                  if( (fRegions.insert(regnum)).second) fRegCounts[regnum] = 1;
                  else fRegCounts[regnum]++;

                  if(passingData.size()>passingIDE.size()) 
                      std::cout << "data/IDE size mismatch! " << passingData.size()-passingIDE.size() << std::endl;
                  FillAdcArr(passingData,adc);
                  dataCol.push_back( std::make_pair(
                    FillCRTData(trg.first,event,ttrig,ttrig,adc),
                    passingIDE) );
                  if (fUltraVerbose) 
                      std::cout << " ...success!" << std::endl;
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

          if(dataCol.size()>0) {
              std::map<int,int>::iterator it = fRegCounts.begin();
              std::cout << '\n' << "FEB events per CRT region: " << '\n' << std::endl;

              do {
                  std::cout << fCrtutils->GetRegionNameFromNum((*it).first) << ": , events: " << (*it).second << '\n' << std::endl;
                  it++;
              }
              while ( it != fRegCounts.end() );
          }//if any CRTData
        } //if verbose

	return dataCol;

    }//end CreateData()

    //-----------------------------------------------------------------------------
    // intented to be called within loop over AuxDetChannels and provided the 
    // AuxDetChannelID, AuxDetSensitiveChannelID, vector of AuxDetIDEs and 
    // the number of ides from the end of the vector to include in the detector sim.
    // handles deposited energy -> light output at SiPM (including attenuation) 
    // -> PEs from the SiPM (or PMT in case of bottom CRT) with associated time stamps
    void CRTDetSimAlg::FillTaggers(detinfo::DetectorClocksData const& clockData,
                                   const uint32_t adid, const uint32_t adsid, const vector<sim::AuxDetIDE>& ides) {

        art::ServiceHandle<geo::Geometry> geoService;
        detinfo::ElecClock trigClock = clockData.TriggerClock();
        fHasFilledTaggers = true;

        const geo::AuxDetGeo& adGeo = geoService->AuxDet(adid); //pointer to module object
        const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsid); //pointer to strip object
        const char auxDetType = fCrtutils->GetAuxDetType(adid); //CRT module type (c, d, or m)
        const string region = fCrtutils->GetAuxDetRegion(adid); //CRT region

        int layid = INT_MAX; //set to 0 or 1 if layerid determined
        uint8_t mac5=UINT8_MAX, mac5dual=UINT8_MAX; //front-end board ID, dual for MINOS modules (not cut)

        // Find the path to the strip geo node, to locate it in the hierarchy
        set<string> volNames = { adsGeo.TotalVolume()->GetName() };
        vector<vector<TGeoNode const*> > paths = geoService->FindAllVolumePaths(volNames);

        string path = "";
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

        if(layid==INT_MAX) 
            mf::LogError("CRT") 
                << "layid NOT SET!!!" << '\n'
                << "   ADType: " << auxDetType << '\n'
                << "   ADRegion: " << region;

        // ------------- loop over AuxDetIDEs ---------
        // Simulate the CRT response for each hit in this strip
        for (auto const& ide : ides ) {
            // count true number of energy deposits for each strip
            if (auxDetType=='c') fNsim_c++;
            if (auxDetType=='d') fNsim_d++;
            if (auxDetType=='m') fNsim_m++;

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
               mf::LogError("CRT") << "HIT POINT OUTSIDE OF SENSITIVE VOLUME!" << '\n'
                                  << "  AD: " << adid << " , ADS: " << adsid << '\n'
                                  << "  Local position (x,y,z): ( " << svHitPosLocal[0]
                                  << " , " << svHitPosLocal[1] << " , " << svHitPosLocal[2] << " )";

            //calculate Birks' correction factor
            double dl = sqrt(pow(ide.entryX-ide.exitX,2)+pow(ide.entryY-ide.exitY,2)+pow(ide.entryZ-ide.exitZ,2));
            double birksCorr = 1./(1+fKbirks*ide.energyDeposited/dl);

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
                pair<double,double> tmp = GetTransAtten(svHitPosLocal[0]);
                abs0 = tmp.first;
                abs1 = tmp.second;
            }

            //# photons arriving at SiPM
            double npeExp0 = npeExpected * abs0;;
            double npeExp1 = npeExpected * abs1;
            double npeExp0Dual = npeExpected2 * abs0;

            //apply Birks' quenching if requested
            if(fUseBirks){
                npeExp0 *= birksCorr;
                npeExp1 *= birksCorr;
                npeExp0Dual *= birksCorr;
            }

            //sanity check on simulated light output
            if(npeExp0<0||npeExp1<0||npeExp0Dual<0) 
                mf::LogError("CRT") << "NEGATIVE PE!!!!!";

            // Observed PE (Poisson-fluctuated)
            int npe0 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0);
            int npe1 = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp1);
            int npe0Dual = CLHEP::RandPoisson::shoot(&fRandEngine, npeExp0Dual);

            // Time relative to trigger [ns], accounting for propagation delay and 'walk'
            // for the fixed-threshold discriminator
            double tTrue = (ide.entryT + ide.exitT) / 2 + fGlobalT0Offset;
            if(tTrue<0) 
                mf::LogError("CRTDetSim")
                  << "negative true time passed to time stamp generator!";
            uint64_t t0 = 
              GetChannelTriggerTicks(trigClock, tTrue, npe0, distToReadout*100);
            uint64_t t1 = 
              GetChannelTriggerTicks(trigClock, tTrue, npe1, distToReadout*100);
            uint64_t t0Dual = 
              GetChannelTriggerTicks(trigClock, tTrue, npe0Dual, distToReadout2*100);

            // Time relative to PPS: Random for now! (FIXME)
            //int ppsTicks = 
            //  CLHEP::RandFlat::shootInt(&fRandEngine, trigClock.Frequency() * 1e6);

            // SiPM and ADC response: Npe to ADC counts
            uint16_t q0 = 
              CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0, fQRMS * sqrt(npe0));
            uint16_t q1 = 
              CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe1, fQRMS * sqrt(npe1));
            uint16_t q0Dual = 
              CLHEP::RandGauss::shoot(&fRandEngine, fQPed + fQSlope * npe0Dual, fQRMS * sqrt(npe0Dual));

            if(q0>fQMax) q0 = fQMax;
            if(q1>fQMax) q1 = fQMax;
            if(q0Dual>fQMax) q0Dual = fQMax;

           // Adjacent channels on a strip are numbered sequentially.
            //
            // In the AuxDetChannelMapAlg methods, channels are identified by an
            // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
            // module, and a channel number from 0 to 32.

            int channel0ID=0, channel1ID=0;
            pair<uint8_t,uint8_t> macs = fCrtutils->ADToMac(adid);
            mac5 = macs.first;
            int changroup = fCrtutils->ADToChanGroup(adid);

            switch (auxDetType){
                case 'c' :
                    channel0ID = 2 * adsid + 0;
                    channel1ID = 2 * adsid + 1;
                    break;
                case 'd' :
                    channel0ID = adsid;
                    break;
                case 'm' :
                    channel0ID = adsid/2 + 10*(changroup-1);
                    if (fCrtutils->NFeb(adid)==2)  {
                        mac5dual = macs.second;
                    }
                    break;
            }

           if (mac5==UINT8_MAX) mf::LogError("CRT") << "mac addrs not set!";

            // Apply ADC threshold and strip-level coincidence (both fibers fire)
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
                        tagger.data.push_back(
                          std::make_pair(FillChanData(channel0ID,q0,t0),ide));
                        fNchandat_c++;
                    }
                    else fNmissthr_c++;
                    if (q1>fQThresholdC) {
                        tagger.data.push_back(
                          std::make_pair(FillChanData(channel1ID,q1,t1),ide));
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
                    tagger.data.push_back(
                      std::make_pair(FillChanData(channel0ID,q0,t0),ide));
                    fNchandat_d++;
            }//if one strip above threshold

            if (auxDetType=='m') {
                    if(q0 > fQThresholdM) {
                      Tagger& tagger = fTaggers[mac5];
                      tagger.layerid.insert(layid);
                      tagger.chanlayer[channel0ID] = layid;
                      tagger.reg = region;
                      tagger.type = 'm';
                      tagger.modid = adid;
                      tagger.data.push_back(
                        std::make_pair(FillChanData(channel0ID,q0,t0),ide));
                      fNchandat_m++;
                    }
                    if(q0Dual > fQThresholdM && fCrtutils->NFeb(adid)==2) {
                      Tagger& tagger = fTaggers[mac5dual];
                      tagger.layerid.insert(layid);
                      tagger.chanlayer[channel0ID] = layid;
                      tagger.reg = region;
                      tagger.type = 'm';
                      tagger.modid = adid;
                      tagger.data.push_back(
                        std::make_pair(FillChanData(channel0ID,q0Dual,t0Dual),ide));
                      fNchandat_m++;
                    }
                    if(q0 > fQThresholdM && q0Dual > fQThresholdM) fNdual_m++;
            }//if one strip above threshold at either end

            //counting losses
            if (auxDetType == 'c') {
                if ( fApplyStripCoinC && util::absDiff(t0,t1) >= fStripCoincidenceWindow ) fNmiss_strcoin_c++;
            }
            if (auxDetType == 'd' && q0 < fQThresholdD) fNmissthr_d++;
            if (auxDetType == 'm') {
                if( q0 < fQThresholdM) fNmissthr_m++;
                if( q0Dual < fQThresholdM && fCrtutils->NFeb(adid)==2) fNmissthr_m++;
            }

            //print detsim info (if enabled)
            if (fUltraVerbose&&
               (
                 (auxDetType=='c' && q0>fQThresholdC && q1>fQThresholdC && util::absDiff(t0, t1) < fStripCoincidenceWindow) ||
                 (auxDetType=='d' && q0>fQThresholdD ) ||
                 (auxDetType=='m' && (q0>fQThresholdM || q0Dual>fQThresholdM)) ))
              std::cout << '\n'
              << "CRT HIT VOL " << (adGeo.TotalVolume())->GetName() << "\n"
              << "CRT HIT SENSITIVE VOL " << (adsGeo.TotalVolume())->GetName() << "\n"
              << "CRT HIT AuxDetID " <<  adid << " / AuxDetSensitiveID " << adsid << "\n"
              << "CRT module type: " << auxDetType << " , CRT region: " << region << '\n'
              << "CRT channel: " << channel0ID << " , mac5: " << (int)mac5 << '\n'
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

    } //end FillTaggers

    //---------------------------------------------------------------
    // for CERN modules only. Given position in local strip coordinates, 
    // calculate the transverse attentuation factor
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
    uint64_t CRTDetSimAlg::GetChannelTriggerTicks(detinfo::ElecClock& clock,
                                             double t0, float npeMean, float r)
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
        clock = clock.WithTime(t / 1e3);  // WithTime takes microseconds

        if (fUltraVerbose) mf::LogInfo("CRT")
          << "CRT TIMING: t0=" << t0
          << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
          << ", tDelay=" << tDelay << ", tDelay(interp)="
          << tDelay << ", tProp=" << tProp << ", t=" << t 
          << ", ticks=" << clock.Ticks() << "\n"; 
        
        return (uint64_t)t;//clock.Ticks();
    }//CRTDetSim::GetChannelTriggerTicks()

    //---------------------------------------------------------------
    // function to clear member data at beginning of each art::event
    void CRTDetSimAlg::ClearTaggers() {

        fTaggers.clear();
        fHasFilledTaggers = false;

        fNsim_m = 0;
        fNsim_d = 0; 
        fNsim_c = 0;
        fNchandat_m = 0;
        fNchandat_d = 0;
        fNchandat_c = 0;
        fNmissthr_c = 0;
        fNmissthr_d = 0;
        fNmissthr_m = 0;
        fNmiss_strcoin_c = 0;
        fNdual_m = 0;

        fRegions.clear();
        fRegCounts.clear();
    }

    //----------------------------------------------------------------
    // function to make fill CRTData products a bit easer
    CRTData CRTDetSimAlg::FillCRTData(uint8_t mac, uint32_t entry, uint64_t t0, uint64_t t1, uint16_t adc[64]){
        CRTData dat;
        dat.fMac5 = mac;
        dat.fEntry = entry;
        dat.fTs0 = t0;
        dat.fTs1 = t1;
        for(size_t chan=0; chan<64; chan++) {
            dat.fAdc[chan] = adc[chan];
        }

        return dat;
    }
 
    //------------------------------------------------------------------
    ChanData CRTDetSimAlg::FillChanData(int channel, uint16_t adc, uint64_t ts) {
        ChanData dat;
        dat.channel = channel;
        dat.adc = adc;
        dat.ts = ts;

        return dat;
    }

    //---------------------------------------------------------------------------
    void CRTDetSimAlg::FillAdcArr(const vector<ChanData>& data, uint16_t arr[64]) {
       for(int chan=0; chan<64; chan++)
           arr[chan] = 0;
       for(auto const& dat : data) {
           if(dat.channel<0 || dat.channel>63){
               std::cout << "ChanData channel out of bounds!" << std::endl;
               return;
           }
           arr[dat.channel] = dat.adc;
       }
       return;
    }
    
 }//namespace crt
}//namespace icarus

#endif
