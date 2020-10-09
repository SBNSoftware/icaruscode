////////////////////////////////////////////////////////////////////////
// Class:       AnaProducer
// Module Type: analyzer
// File:        AnaProducer_module.cc
// Description: Makes a tree with waveform information.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"

//#include "art/Framework/Services/Optional/TFileService.h"
#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TNtuple.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <iostream>

#include "icaruscode/CRT/CRTDecoder/CrtCal.h"
#include "icaruscode/CRT/CRTDecoder/CrtCalTree.h"
#include "icaruscode/CRT/CRTDecoder/CRTTiming.h"
#include "icaruscode/CRT/CRTDecoder/CRTRawTree.h"

using std::vector;
using std::map;

namespace icarus {
 namespace crt {
  class AnaProducer;
 }
}

class icarus::crt::AnaProducer : public art::EDAnalyzer {

public:
   struct Config {

       // Save some typing:
       using Name = fhicl::Name;
       using Comment = fhicl::Comment;

       fhicl::Atom<float> PeThresh {
         Name("peThreshold"),
         Comment("used for determining if a channel is above the discriminator threshold")
         };

       // One Atom for each parameter
       fhicl::Atom<std::string> CalFile {
         Name("calFile"),
         Comment("Full path to calibration file")
         };
       // One Atom for each parameter
       fhicl::Atom<bool> Calibrate {
         Name("calibrate"),
         Comment("calibrate the data (true) or read from cal file (false)")
         };

       fhicl::Atom<uint8_t> Region {
         Name("region"),
         Comment("CRT region e.g. north (n) or west (w)")
         };
       fhicl::Atom<uint8_t> Layer {
         Name("layer"),
         Comment("CRT layer inner (i) or outer (o)")
         };

       fhicl::Sequence<uint8_t> MacsNI {
         Name("macsNI"),
         Comment("mac addresses for each FEB in the north, inner daisychain")
         };
       fhicl::Sequence<uint8_t> MacsNO {
         Name("macsNO"),
         Comment("mac addresses for each FEB in the north, outer daisychain")
         };
       fhicl::Sequence<uint8_t> MacsWI {
         Name("macsWI"),
         Comment("mac addresses for each FEB in the west, inner daisychain")
         };
       fhicl::Sequence<uint8_t> MacsWO {
         Name("macsWO"),
         Comment("mac addresses for each FEB in the west, outer daisychain")
         };

   };

   using Parameters = art::EDAnalyzer::Table<Config>;

   explicit AnaProducer(Parameters const& config); // explicit doesn't allow for copy initialization
   virtual ~AnaProducer();
 
   virtual void beginJob();
   virtual void analyze(art::Event const & evt);
   virtual void endJob();
  

private:
   void analyze_fragment(artdaq::Fragment &);

   string  pFile;
   bool    pCalibrate;
   vector<uint8_t> pMacs;
   float   pPeThresh;

   map<uint8_t,vector<TH1F*>*> fMacToHistos;
   map<uint8_t,vector<TH1F*>*> fMacToPEHistos;

   TTree *fRawTree, *fCalTree, *fAnaTree;

   //rawTree
     //data payload
   uint8_t   fMac5; //last 8 bits of FEB mac5 address
   uint16_t  fFlags;
   uint16_t  fLostcpu;
   uint16_t  fLostfpga;
   uint32_t  fTs0;
   uint32_t  fTs1;
   uint16_t  fAdc[32];
   uint32_t  fCoinc;

     //metadata
   uint64_t  fRun_start_time;
   uint64_t  fThis_poll_start;
   uint64_t  fThis_poll_end;
   uint64_t  fLast_poll_start;
   uint64_t  fLast_poll_end;
   int32_t   fSystem_clock_deviation;
   uint32_t  fFeb_per_poll;
   uint32_t  fFeb_event_number;

      //information from fragment header
   uint32_t  fSequence_id;
   uint64_t  fFragment_timestamp;

   //anaTree
   bool    fIsNoise;
   uint8_t fMaxChan;
   float   fMaxPE;
   float   fTotPE;
   float   fPE[32];
   uint8_t fNChanAbove;
   bool    fAbove[32];
   uint64_t fT0;
   int     fRegion;
   int     fLayer;
   float   fPollRate;
   float   fInstRate;

   //calTree
   bool    fActive[32];
   float   fGain[32];
   float   fGainErr[32];
   float   fGainXsqr[32];
   short   fGainNdf[32];
   float   fGainPed[32];
   float   fGainPedErr[32];
   short   fNpeak[32];
   float   fPed[32];
   float   fPedErr[32];
   float   fPedXsqr[32];
   short   fPedNdf[32];
   float   fPedSigma[32];
   float   fPedSigmaErr[32];
   float   fPedNorm[32];
   float   fPedNormErr[32];
   int     fThreshADC[32];
   float   fThreshPE[32];
   int     fNabove[32];
   float   fPeakNorm[32][5];
   float   fPeakNormErr[32][5];
   float   fPeakSigma[32][5];
   float   fPeakSigmaErr[32][5];
   float   fPeakMean[32][5];
   float   fPeakMeanErr[32][5];
   float   fPeakXsqr[32][5];
   short   fPeakNdf[32][5];
  
};

//Define the constructor
icarus::crt::AnaProducer::AnaProducer(Parameters const& config)
  : EDAnalyzer(config), 
    pFile(config().CalFile()), 
    pCalibrate(config().Calibrate()), 
    //pMacs(config().Macs()), 
    pPeThresh(config().PeThresh()),
    fRegion(config().Region()),
    fLayer(config().Layer())
{

  if     (fRegion==0 && fLayer==1) pMacs = config().MacsNI();
  else if(fRegion==0 && fLayer==0) pMacs = config().MacsNO();
  else if(fRegion==1 && fLayer==1) pMacs = config().MacsWI();
  else if(fRegion==1 && fLayer==0) pMacs = config().MacsWO();
  else std::cout << "ERROR in AnaProducer::AnaProducer: bad region and/or layer codes!" << std::endl;

  //this is how you setup/use the TFileService
  //I do this here in the constructor to setup the histogram just once
  //but you can call the TFileService and make new objects from anywhere
  art::ServiceHandle<art::TFileService> tfs; //pointer to a file named tfs

  fCalTree = 0;
  if(pCalibrate)
	fCalTree = tfs->make<TTree>("calTree","gain, pedestal, and fit statistics");
  fRawTree = tfs->make<TTree>("rawTree", "fields from DAQ framqments with metadata");
  fAnaTree = tfs->make<TTree>("anaTree","calibrated charge, low-level ana");

  //rawTree
  fRawTree->Branch("mac5",                   &fMac5,                   "mac5/b");
  fRawTree->Branch("flags",                  &fFlags,                  "flags/s");
  fRawTree->Branch("lostcpu",                &fLostcpu,                "lostcpu/s");
  fRawTree->Branch("lostfpga",               &fLostfpga,               "lostfpga/s");
  fRawTree->Branch("ts0",                    &fTs0,                    "ts0/i");
  fRawTree->Branch("ts1",                    &fTs1,                    "ts1/i");
  fRawTree->Branch("adc",                    &fAdc,                    "adc[32]/s");
  fRawTree->Branch("coinc",                  &fCoinc,                  "coinc/i");
  fRawTree->Branch("run_start_time",         &fRun_start_time,         "run_start_time/l");
  fRawTree->Branch("this_poll_start",        &fThis_poll_start,        "this_poll_start/l");
  fRawTree->Branch("this_poll_end",          &fThis_poll_end,          "this_poll_end/l");
  fRawTree->Branch("last_poll_start",        &fLast_poll_start,        "last_poll_start/l");
  fRawTree->Branch("last_poll_end",          &fLast_poll_end,          "last_poll_end/l");
  fRawTree->Branch("system_clock_deviation", &fSystem_clock_deviation, "system_clock_deviation/I");
  fRawTree->Branch("feb_per_poll",           &fFeb_per_poll,           "feb_per_poll/i");
  fRawTree->Branch("feb_event_number",       &fFeb_event_number,       "feb_event_number/i");
  fRawTree->Branch("sequence_id",            &fSequence_id,            "sequence_id/i");
  fRawTree->Branch("fragment_timestamp",     &fFragment_timestamp,     "fragment_timestamp/l");

  //anaTree
  fAnaTree->Branch("mac5",           &fMac5,                "mac5/b");
  fAnaTree->Branch("pe",             fPE,                   "pe[32]/F");
  fAnaTree->Branch("active",         fActive,               "active[32]/O");
  fAnaTree->Branch("maxChan",        &fMaxChan,             "maxChan/b");
  fAnaTree->Branch("maxPE",          &fMaxPE,               "maxPE/F"); 
  fAnaTree->Branch("totPE",          &fTotPE,               "totPE/F");
  fAnaTree->Branch("nAbove",         &fNChanAbove,          "nAbove/b");
  fAnaTree->Branch("above",          fAbove,                "above[32]/O");
  fAnaTree->Branch("isNoise",        &fIsNoise,             "isNoise/O");
  fAnaTree->Branch("region",         &fRegion,              "region/I");
  fAnaTree->Branch("layer",          &fLayer,               "layer/I");
  fAnaTree->Branch("t0",             &fT0,                  "t0/l");
  fAnaTree->Branch("pollRate",       &fPollRate,            "pollRate/F");
  fAnaTree->Branch("instRate",       &fInstRate,            "instRate/F");

  //calTree
  if(pCalibrate){
     fCalTree->Branch("mac5",         &fMac5,        "mac5/b");
     fCalTree->Branch("active",       fActive,       "active[32]/O");
     fCalTree->Branch("gain",         fGain,         "gain[32]/F");
     fCalTree->Branch("gainErr",      fGainErr,      "gainErr[32]/F");
     fCalTree->Branch("gainXsqr",     fGainXsqr,     "gainXsqr[32]/F");
     fCalTree->Branch("gainNdf",      fGainNdf,      "gainNdf[32]/s");
     fCalTree->Branch("gainPed",      fGainPed,      "gainPed[32]/F");
     fCalTree->Branch("gainPedErr",   fGainPedErr,   "gainPedErr[32]/F");
     fCalTree->Branch("nPeak",        fNpeak,        "nPeak[32]/s");
     fCalTree->Branch("peakXsqr",     fPeakXsqr,     "peakXsqr[32][5]/F");
     fCalTree->Branch("peakNdf",      fPeakNdf,      "peakNdf[32][5]/s");
     fCalTree->Branch("peakMean",     fPeakMean,     "peakMean[32][5]/F");
     fCalTree->Branch("peakMeanErr",  fPeakMeanErr,  "peakMeanErr[32][5]/F");
     fCalTree->Branch("peakNorm",     fPeakNorm,     "peakNorm[32][5]/F");
     fCalTree->Branch("peakNormErr",  fPeakNormErr,  "peakNormErr[32][5]/F");
     fCalTree->Branch("peakSigma",    fPeakSigma,    "peakSigma[32][5]/F");
     fCalTree->Branch("peakSigmaErr", fPeakSigmaErr, "peakSigmaErr[32][5]/F");
     fCalTree->Branch("ped",          fPed,          "ped[32]/F");
     fCalTree->Branch("pedErr",       fPedErr,       "pedErr[32]/F");
     fCalTree->Branch("pedXsqr",      fPedXsqr,      "pedXsqr[32]/F");
     fCalTree->Branch("pedNdf",       fPedNdf,       "pedNdf[32]/s");
     fCalTree->Branch("pedSigma",     fPedSigma,     "pedSigma[32]/F");
     fCalTree->Branch("pedSigmaErr",  fPedSigmaErr,  "pedSigmaErr[32]/F");
     fCalTree->Branch("pedNorm",      fPedNorm,      "pedNorm[32]/F");
     fCalTree->Branch("pedNormErr",   fPedNormErr,   "pedNormErr[32]/F");
     fCalTree->Branch("threshAdc",    fThreshADC,    "threshAdc[32]/I");
     fCalTree->Branch("threshPe",     fThreshPE,     "threshPe[32]/F");
     fCalTree->Branch("nAbove",       fNabove,       "nAbove[32]/I");

  }//end if pCal


  //initialize vars
  fMac5 = 0;
  fFlags = 0;
  fLostcpu = 0;
  fLostfpga = 0;
  fTs0 = 0;
  fTs1 = 0;
  fAdc[32] = 0;
  fCoinc = 0;
  fRun_start_time = 0;
  fThis_poll_start = 0;
  fThis_poll_end = 0;
  fLast_poll_start = 0;
  fLast_poll_end = 0;
  fSystem_clock_deviation = 0;
  fFeb_per_poll = 0;
  fFeb_event_number = 0;
  fSequence_id = 0;
  fFragment_timestamp = 0;

  fIsNoise = false;
  fMaxChan = 33;
  fMaxPE = -1.;
  fTotPE = -1.;
  fNChanAbove = 0;
  fT0 = 0;
  fRegion = config().Region();
  fLayer = config().Layer();
  fPollRate = 0.;
  fInstRate = 0.;

  
  //generate histograms
  for(int i=0; i<(int)pMacs.size(); i++){
      fMacToHistos[pMacs[i]] = new vector<TH1F*>();
      fMacToPEHistos[pMacs[i]] = new vector<TH1F*>();
      for(int ch=0; ch<32; ch++){
          string hname = "hadc_"+to_string(pMacs[i])+"_"+to_string(ch);
          string htitle = "raw charge: mac5 "+to_string(pMacs[i])+", ch. "+to_string(ch);
          fMacToHistos[pMacs[i]]->push_back(tfs->make<TH1F>(hname.c_str(),htitle.c_str(),4100,0,4100));

	  string hnamePE = "hpe_"+to_string(pMacs[i])+"_"+to_string(ch);
          string htitlePE = "calibrated charge: mac5 "+to_string(pMacs[i])+", ch. "+to_string(ch);
          fMacToPEHistos[pMacs[i]]->push_back(tfs->make<TH1F>(hnamePE.c_str(),htitlePE.c_str(),2000,-5,95));
      }
  }

  for(int i=0; i<32; i++) {
      fPE[i] = -1.;
      fAbove[i] = false;
      fActive[i] = false;
      fGain[i] = FLT_MAX;
      fGainErr[i] = FLT_MAX;
      fGainXsqr[i] = FLT_MAX;
      fGainNdf[i] = SHRT_MAX;
      fGainPed[i] = FLT_MAX;
      fGainPedErr[i] = FLT_MAX;
      fNpeak[i] = SHRT_MAX;
      fPed[i] = FLT_MAX;
      fPedErr[i] = FLT_MAX;
      fPedXsqr[i] = FLT_MAX;
      fPedNdf[i] = SHRT_MAX;
      fPedSigma[i] = FLT_MAX;
      fPedSigmaErr[i] = FLT_MAX;
      fPedNorm[i] = FLT_MAX;
      fPedNormErr[i] = FLT_MAX;
      fThreshADC[i] = INT_MAX;
      fThreshPE[i] = FLT_MAX;
      fNabove[i] = INT_MAX;

      for(size_t j=0; j<5; j++){
          fPeakNorm[i][j] = FLT_MAX;
          fPeakNormErr[i][j] = FLT_MAX;
          fPeakSigma[i][j] = FLT_MAX;
          fPeakSigmaErr[i][j] = FLT_MAX;
          fPeakMean[i][j] = FLT_MAX;
          fPeakMeanErr[i][j] = FLT_MAX;
          fPeakXsqr[i][j] = FLT_MAX;
          fPeakNdf[i][j] = SHRT_MAX;
      }
  }

}

icarus::crt::AnaProducer::~AnaProducer()
{
}

void icarus::crt::AnaProducer::beginJob(){

}

void icarus::crt::AnaProducer::analyze_fragment(artdaq::Fragment & frag)
{
  sbndaq::BernCRTFragment bern_fragment(frag);

  fFragment_timestamp        = frag.timestamp();
  fSequence_id               = frag.sequenceID();

  //event data
  sbndaq::BernCRTEvent const* bevt = bern_fragment.eventdata();

  fMac5     = bevt->MAC5();
  fFlags    = bevt->flags;
  fLostcpu  = bevt->lostcpu;
  fLostfpga = bevt->lostfpga;
  fTs0      = bevt->Time_TS0();
  fTs1      = bevt->Time_TS1();
  fCoinc    = bevt->coinc;

  for(int ch=0; ch<32; ch++){
    fMacToHistos[fMac5]->at(ch)->Fill( bevt->ADC(ch));
    fAdc[ch] = bevt->ADC(ch);
  }

  //metadata
  const sbndaq::BernCRTFragmentMetadata* md = bern_fragment.metadata();

  fRun_start_time            = md->run_start_time();
  fThis_poll_start           = md->this_poll_start();
  fThis_poll_end             = md->this_poll_end();
  fLast_poll_start           = md->last_poll_start();
  fLast_poll_end             = md->last_poll_end();
  fSystem_clock_deviation    = md->system_clock_deviation();
  fFeb_per_poll              = md->feb_events_per_poll();
  fFeb_event_number          = md->feb_event_number();

  fRawTree->Fill();

} //analyze_fragment


void icarus::crt::AnaProducer::analyze(art::Event const & evt)
{ 
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  evt.getManyByType(fragmentHandles); 
  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;
    
    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      //Container fragment
      for (auto cont : *handle) { 
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() != sbndaq::detail::FragmentType::BERNCRT)
          continue; 
        for (size_t ii = 0; ii < contf.block_count(); ++ii)
          analyze_fragment(*contf[ii].get());
      }   
    }   
    else {
      //normal fragment
      if (handle->front().type() != sbndaq::detail::FragmentType::BERNCRT) continue;
      for (auto frag : *handle)
        analyze_fragment(frag);
    }   
  }
} //analyze



void icarus::crt::AnaProducer::endJob(){

        std::cout << "done filling histograms..." << std::endl;
        std::cout << "found " << fMacToHistos.size() << " FEBs" << std::endl;
        if(!fMacToHistos.begin()->second->empty()){
                std::cout << "first histo size: " << fMacToHistos.begin()->second->at(0)->Integral()
                          << " entries" << std::endl;
        }
        else {
                std::cout << "hist vect is empty!" << std::endl;
        }


	if(pCalibrate)
        for(auto const& macHist : fMacToHistos){

                uint8_t mac = macHist.first;
                vector<TH1F*>* hvec = macHist.second;

                std::cout << "construct instance of CrtCal for mac5 " << (short)mac << std::endl;
                CrtCal* cal = new CrtCal(hvec);
                cal->Cal();

                std::cout << "retreiving cal data..." << std::endl;
                fMac5         = macHist.first;
                bool*    ptrActive       = cal->GetActive();
                float*   ptrGain         = cal->GetGain();
                float*   ptrGainErr      = cal->GetGainErr();
                float*   ptrGainXsqr     = cal->GetGainXsqr();
                short*   ptrGainNdf      = cal->GetGainNdf();
                float*   ptrGainPed      = cal->GetGainPed();
                float*   ptrGainPedErr   = cal->GetGainPedErr();
                short*   ptrNpeak        = cal->GetNpeak();
                float**  ptrPeakNorm     = cal->GetPeakNorm();
                float**  ptrPeakNormErr  = cal->GetPeakNormErr();
                float**  ptrPeakSigma    = cal->GetPeakSigma();
                float**  ptrPeakSigmaErr = cal->GetPeakSigmaErr();
                float**  ptrPeakMean     = cal->GetPeakMean();
                float**  ptrPeakMeanErr  = cal->GetPeakMeanErr();
                float**  ptrPeakXsqr     = cal->GetPeakXsqr();
                short**  ptrPeakNdf      = cal->GetPeakNdf();
                float*   ptrPed          = cal->GetPed();
                float*   ptrPedErr       = cal->GetPedErr();
                float*   ptrPedXsqr      = cal->GetPedXsqr();
                short*   ptrPedNdf       = cal->GetPedNdf();
                float*   ptrPedSigma     = cal->GetPedSigma();
                float*   ptrPedSigmaErr  = cal->GetPedSigmaErr();
                float*   ptrPedNorm      = cal->GetPedNorm();
                float*   ptrPedNormErr   = cal->GetPedNormErr();
                int*  ptrThreshADC    = cal->GetThreshADC();
                float*   ptrThreshPE     = cal->GetThreshPE();
                int*  ptrNabove       = cal->GetNabove();

                //copy values to local arrays
                for(size_t i=0; i<32; i++){
                        if(ptrActive!=nullptr)      fActive[i]       = ptrActive[i];
                        if(ptrGain!=nullptr)        fGain[i]         = ptrGain[i];
                        if(ptrGainErr!=nullptr)     fGainErr[i]      = ptrGainErr[i];
                        if(ptrGainXsqr!=nullptr)    fGainXsqr[i]     = ptrGainXsqr[i];
                        if(ptrGainNdf!=nullptr)     fGainNdf[i]      = ptrGainNdf[i];
                        if(ptrGainPed!=nullptr)     fGainPed[i]      = ptrGainPed[i];
                        if(ptrGainPedErr!=nullptr)  fGainPedErr[i]   = ptrGainPedErr[i];
                        if(ptrNpeak!=nullptr)       fNpeak[i]        = ptrNpeak[i];
                        if(ptrPed!=nullptr)         fPed[i]          = ptrPed[i];
                        if(ptrPedErr!=nullptr)      fPedErr[i]       = ptrPedErr[i];
                        if(ptrPedXsqr!=nullptr)     fPedXsqr[i]      = ptrPedXsqr[i];
                        if(ptrPedNdf!=nullptr)      fPedNdf[i]       = ptrPedNdf[i];
                        if(ptrPedSigma!=nullptr)    fPedSigma[i]     = ptrPedSigma[i];
                        if(ptrPedSigmaErr!=nullptr) fPedSigmaErr[i]  = ptrPedSigmaErr[i];
                        if(ptrPedNorm!=nullptr)     fPedNorm[i]      = ptrPedNorm[i];
                        if(ptrPedNormErr!=nullptr)  fPedNormErr[i]   = ptrPedNormErr[i];
                        if(ptrThreshADC!=nullptr)   fThreshADC[i]    = ptrThreshADC[i];
                        if(ptrThreshPE!=nullptr)    fThreshPE[i]     = ptrThreshPE[i];
                        if(ptrNabove!=nullptr)      fNabove[i]       = ptrNabove[i];

                        for(size_t j=0; j<5; j++){
                                if(ptrPeakNorm!=nullptr)     fPeakNorm[i][j]     = ptrPeakNorm[i][j];
                                if(ptrPeakNormErr!=nullptr)  fPeakNormErr[i][j]  = ptrPeakNormErr[i][j];
                                if(ptrPeakSigma!=nullptr)    fPeakSigma[i][j]    = ptrPeakSigma[i][j];
                                if(ptrPeakSigmaErr!=nullptr) fPeakSigmaErr[i][j] = ptrPeakSigmaErr[i][j];
                                if(ptrPeakMean!=nullptr)     fPeakMean[i][j]     = ptrPeakMean[i][j];
                                if(ptrPeakMeanErr!=nullptr)  fPeakMeanErr[i][j]  = ptrPeakMeanErr[i][j];
                                if(ptrPeakXsqr!=nullptr)     fPeakXsqr[i][j]     = ptrPeakXsqr[i][j];
                                if(ptrPeakNdf!=nullptr)      fPeakNdf[i][j]      = ptrPeakNdf[i][j];
                        }
                }

                std::cout << "fill tree event" << std::endl;
                fCalTree->Fill();

                delete cal;

        }//end for mac5s if(pCalibrate)

	icarus::crt::CrtCalTree* ccl;
	if(pCalibrate) ccl = new CrtCalTree(fCalTree);
	else ccl = new CrtCalTree(pFile);

	//ccl->Dump();

	const size_t nRaw = fRawTree->GetEntriesFast();

	std::cout << "initiallizing time utility" << std::endl;
	icarus::crt::CRTRawTree traw(fRawTree);
	icarus::crt::CRTTiming time(traw);
	const map<size_t,size_t>* sortedToRaw = time.GetOrderedToRawMap();
	std::cout << "done. sorted through " << sortedToRaw->size() << " entries" << std::endl;
	if(sortedToRaw->size()!=nRaw)
		std::cout << "WARNING: sort map and rawTree are of different size!" << std::endl;

	//time.DumpRawTimes(0);
	//time.DumpSortedTimes(0);

	std::cout << "filling AnaTree..." << std::endl;
	for(size_t iraw=0; iraw<nRaw; iraw++) {
		//if(sortedToRaw->find(iraw)==sortedToRaw->end()) {
		//	std::cout << "raw tree index not found in time order map!" << std::endl;
		//}
		if(iraw%10000==0) std::cout << 100.0*iraw/nRaw << " %" << std::endl;
		fMaxPE = -1.;
		fTotPE = 0.;
		fMaxChan = 33;
		fIsNoise = false;
		fNChanAbove = 0;
		fMac5 = traw.GetMac(sortedToRaw->at(iraw));		
		fT0 = traw.GetAbsTime(sortedToRaw->at(iraw));
		//std::cout << "T0: " << fT0 << std::endl;
		fPollRate = traw.GetPollRate(sortedToRaw->at(iraw));
		if(iraw!=0 && iraw != nRaw-1) 
			fInstRate = traw.GetInstRate(sortedToRaw->at(iraw-1),sortedToRaw->at(iraw+1));
		else
			fInstRate = fPollRate; 

		for(size_t ch=0; ch<32; ch++) {
			fAbove[ch] = false;
			fActive[ch] = ccl->GetActive(fMac5,ch);
			fPE[ch] = 0.0;
			if(!fActive[ch]) continue;

			uint16_t adc = traw.GetADC(sortedToRaw->at(iraw),ch);
			fPE[ch] = (adc-ccl->GetPed(fMac5,ch))/ccl->GetGain(fMac5,ch);
			fMacToPEHistos[fMac5]->at(ch)->Fill(fPE[ch]);
			if(fPE[ch]>fMaxPE){
				fMaxPE = fPE[ch];
				fMaxChan = ch;
			}
			if(fPE[ch] > pPeThresh) {
				fTotPE+=fPE[ch];
				fAbove[ch] = true;
				fNChanAbove++;
			}
		}
		if(fMaxPE<pPeThresh) fIsNoise = true;

		fAnaTree->Fill();
	}//end 
	std::cout << "done" << std::endl;

}

DEFINE_ART_MODULE(icarus::crt::AnaProducer)
//this is where the name is specified
