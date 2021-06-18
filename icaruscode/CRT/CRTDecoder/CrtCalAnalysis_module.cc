////////////////////////////////////////////////////////////////////////
// Class:       CrtCalAnalyzer
// Module Type: analyzer
// File:        CrtCalAnalyzer_module.cc
// Description: Makes a tree with waveform information.
////////////////////////////////////////////////////////////////////////

//framwork includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "canvas/Utilities/Exception.h"

//local includes
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
#include "icaruscode/CRT/CRTDecoder/CrtCal.h"

//ROOT includes
#include "art_root_io/TFileService.h"
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

namespace icarus {
 namespace crt {
  class CrtCalAnalyzer;
 }
}

using std::map;
using std::vector;
using std::string;
using std::to_string;
using icarus::crt::CrtCal;

class icarus::crt::CrtCalAnalyzer : public art::EDAnalyzer {

public:
  struct Config {

      // Save some typing:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      // One Atom for each parameter

      fhicl::Sequence<uint8_t> Mac {
        Name("mac5"),
        Comment("mac addresses for each FEB in the diasychain")
        };
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CrtCalAnalyzer(Parameters const& config); // explicit doesn't allow for copy initialization
  virtual ~CrtCalAnalyzer();

  virtual void beginJob(); 
  virtual void analyze(art::Event const & evt);
  virtual void endJob();


private:
  void analyze_fragment(artdaq::Fragment &);

  art::ServiceHandle<art::TFileService> tfs;

  vector<uint8_t> macs; 
  map<uint8_t,vector<TH1F*>*> macToHistos;
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
};

//Define the constructor
icarus::crt::CrtCalAnalyzer::CrtCalAnalyzer(Parameters const& config)
  : EDAnalyzer(config) , macs(config().Mac())
{
  fMac5 = 0;
  for(int i=0; i<(int)macs.size(); i++){
      macToHistos[macs[i]] = new vector<TH1F*>();
      for(int ch=0; ch<32; ch++){
          string hname = "hadc_"+to_string(macs[i])+"_"+to_string(ch);
          string htitle = "raw charge: mac5 "+to_string(macs[i])+", ch. "+to_string(ch);
          macToHistos[macs[i]]->push_back(tfs->make<TH1F>(hname.c_str(),htitle.c_str(),4100,0,4100));
      }
  }

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

icarus::crt::CrtCalAnalyzer::~CrtCalAnalyzer()
{
}

void icarus::crt::CrtCalAnalyzer::beginJob(){
  
  calTree = tfs->make<TTree>("calAnalyzerTree", "SiPM/FEB channel calibration data");

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

}

void icarus::crt::CrtCalAnalyzer::analyze_fragment(artdaq::Fragment & frag) {

  sbndaq::BernCRTFragment bern_fragment(frag);
  sbndaq::BernCRTEvent const* bevt = bern_fragment.eventdata();

  fMac5     = bevt->MAC5();
  
  if(macToHistos.find(fMac5) == macToHistos.end()) return;
  
  for(int ch=0; ch<32; ch++) {
    macToHistos[fMac5]->at(ch)->Fill( bevt->ADC(ch));
  }

}//end analyze

void icarus::crt::CrtCalAnalyzer::analyze(art::Event const & evt)
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


void icarus::crt::CrtCalAnalyzer::endJob(){

	std::cout << "done filling histograms..." << std::endl;
	std::cout << "found " << macToHistos.size() << " FEBs" << std::endl;
        if(!macToHistos.begin()->second->empty()){
		std::cout << "first histo size: " << macToHistos.begin()->second->at(0)->Integral() 
                          << " entries" << std::endl;
	}
	else {
		std::cout << "hist vect is empty!" << std::endl;
	}

	for(auto const& macHist	: macToHistos){

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
		calTree->Fill();

		delete cal;
	}

	//calTree->Write();

}

DEFINE_ART_MODULE(icarus::crt::CrtCalAnalyzer)
//this is where the name is specified
