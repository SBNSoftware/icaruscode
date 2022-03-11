//#undef CRT_CAL_CC
//#include "icaruscode/CRT/CRTDecoder/CrtCal.h"
#include "/icarus/app/users/tboone/decoder/icaruscode/srcs/icaruscode/icaruscode/CRT/CRTDecoder/edits/CrtCal.cc"
//#include "icaruscode/CRT/CRTDecoder/CrtCal.cc"

//ROOT includes
#include "TH1F.h"
#include "TNtuple.h"

//c++ includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <map>
#include <string>

void add_statistics_calib(){
  
  
  using std::map;
  using std::vector;
  using std::string;
  using std::to_string;
  using icarus::crt::CrtCal;
  
  TTree* calTree;

  //vars for calTree

  uint8_t fMac5;
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
  long    fChanStats[32];

  for(int i=0; i<32; i++) {
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
    fChanStats[i] = LONG_MAX;

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

  calTree = new TTree("calAnalyzerTree", "SiPM/FEB channel calibration data");

  calTree->Branch("mac5",         &fMac5,        "mac5/b");
  calTree->Branch("active",       fActive,       "active[32]/O");
  calTree->Branch("gain",         fGain,         "gain[32]/F");
  calTree->Branch("gainErr",      fGainErr,      "gainErr[32]/F");
  calTree->Branch("gainXsqr",     fGainXsqr,     "gainXsqr[32]/F");
  calTree->Branch("gainNdf",      fGainNdf,      "gainNdf[32]/s");
  calTree->Branch("gainPed",      fGainPed,      "gainPed[32]/F");
  calTree->Branch("gainPedErr",   fGainPedErr,   "gainPedErr[32]/F");
  calTree->Branch("nPeak",        fNpeak,        "nPeak[32]/s");
  calTree->Branch("peakXsqr",     fPeakXsqr,     "peakXsqr[32][5]/F");
  calTree->Branch("peakNdf",      fPeakNdf,      "peakNdf[32][5]/s");
  calTree->Branch("peakMean",     fPeakMean,     "peakMean[32][5]/F");
  calTree->Branch("peakMeanErr",  fPeakMeanErr,  "peakMeanErr[32][5]/F");
  calTree->Branch("peakNorm",     fPeakNorm,     "peakNorm[32][5]/F");
  calTree->Branch("peakNormErr",  fPeakNormErr,  "peakNormErr[32][5]/F");
  calTree->Branch("peakSigma",    fPeakSigma,    "peakSigma[32][5]/F");
  calTree->Branch("peakSigmaErr", fPeakSigmaErr, "peakSigmaErr[32][5]/F");
  calTree->Branch("ped",          fPed,          "ped[32]/F");
  calTree->Branch("pedErr",       fPedErr,       "pedErr[32]/F");
  calTree->Branch("pedXsqr",      fPedXsqr,      "pedXsqr[32]/F");
  calTree->Branch("pedNdf",       fPedNdf,       "pedNdf[32]/s");
  calTree->Branch("pedSigma",     fPedSigma,     "pedSigma[32]/F");
  calTree->Branch("pedSigmaErr",  fPedSigmaErr,  "pedSigmaErr[32]/F");
  calTree->Branch("pedNorm",      fPedNorm,      "pedNorm[32]/F");
  calTree->Branch("pedNormErr",   fPedNormErr,   "pedNormErr[32]/F");
  calTree->Branch("threshAdc",    fThreshADC,    "threshAdc[32]/I");
  calTree->Branch("threshPe",     fThreshPE,     "threshPe[32]/F");
  calTree->Branch("nAbove",       fNabove,       "nAbove[32]/I");
  calTree->Branch("stats",	  fChanStats,	 "stats[32]/L");

/*  std::vector<int> macadd = {1,3,6,7,
                             4,5,8,9,
                             10,11,13,14,15,16,    // west-center outer
                             17,19,21,22,23,24,    // west-center inner
                             25,26,27,28,29,30,    // east-center outer
                             31,32,33,34,35,36,    // east-center inner
                             37,38,39,40,41,42,    // east-north outer
                             43,44,45,46,47,48,    // east-north inner
                             49,50,51,52,53,54,    // east-south outer
                             55,56,57,58,59,60,    // east-south inner
                             61,62,63,64,65,66,    // west-north outer
                             67,68,69,70,71,72,    // west-north inner
                             73,74,75,76,77,78,    // west-south outer
                             79,80,81,82,83,84,    // west-south inner
			     85,86,87,88,89,90,91, // south outer
			     92,93,94,95,96,97};   // south inner
*/
std::vector<int> macadd = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93};


  TFile* fin = new TFile("/icarus/app/users/tboone/decoder/cal/adc_histos/run_7262_histos.root");
    
  if (fin->IsZombie()) {
    std::cout << "file " << fin << " not found" << std::endl;
  }
  
  map<uint8_t,vector<TH1F*>*> macToHistos;

  for (int mac=0; mac < macadd.size(); mac++){

    macToHistos[macadd[mac]] = new vector<TH1F*>();

    for(int ch=0; ch<32; ch++){
      string hname = "hadc_"+std::to_string(macadd[mac])+"_"+std::to_string(ch);
      string htitle = "raw charge: mac5 "+std::to_string(macadd[mac])+", ch. "+std::to_string(ch);
      TH1F* h = new TH1F( hname.c_str() , "htitle.c_str()" , 4100 , 0. , 4100.);
      h = (TH1F*)fin->Get(("CRTCalibrationAnalysis/hadc_"+std::to_string(macadd[mac])+"_"+std::to_string(ch)).c_str());
      if (!h) continue;
      macToHistos[macadd[mac]]->push_back(h);
    } // # channel 
  } // mac address

  std::cout << "done filling histograms..." << std::endl;
  std::cout << "found " << macToHistos.size() << " FEBs" << std::endl;
  if(!macToHistos.begin()->second->empty()){
    std::cout << "first histo size: " << macToHistos.begin()->second->at(0)->Integral() 
	      << " entries" << std::endl;
  }
  else {
    std::cout << "hist vect is empty!" << std::endl;
  }

vector<TH1F*>* hvec;
  for(auto const& macHist: macToHistos){
    uint8_t mac = macHist.first;
    hvec = macHist.second;
    if (hvec->size() == 0) continue;
    std::cout<< (short)mac <<"\t" <<hvec->size() << std::endl;
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
    long* ptrChanStats	     = cal->GetChanStats();

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
      if(ptrChanStats!=nullptr)	  fChanStats[i]	   = ptrChanStats[i];

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
    calTree->Fill();

    delete cal;
  }

  TFile fout("ftest.root", "RECREATE");
  calTree->Write();

}// void
