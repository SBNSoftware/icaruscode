////////////////////////////////////////////////////////////////////////
// Class:       CrtNoiseMonTool
// Module Type: analyzer
// File:        CrtNoiseMonTool_module.cc
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

#include "sbndaq-artdaq-core/Overlays/Common/BernCRTTranslator.hh"

//#include "art/Framework/Services/Optional/TFileService.h"
#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TGraph.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>

namespace icarus {
 namespace crt {
  class CrtNoiseMonTool;
 }
}

class icarus::crt::CrtNoiseMonTool : public art::EDAnalyzer {

public:
  explicit CrtNoiseMonTool(fhicl::ParameterSet const & pset); // explicit doesn't allow for copy initialization
  virtual ~CrtNoiseMonTool();

  virtual void endJob();
  virtual void analyze(art::Event const & evt);
  

private:
  //  std::map<uint8_t,std::vector<float>> macToRate;
  //std::map<uint8_t,std::vector<uint64_t>> macToRateTime;
  bool     IsSideCRT(icarus::crt::BernCRTTranslator & hit);
  std::map<uint8_t,std::vector<float>> macToRateSide;
  std::map<uint8_t,std::vector<uint64_t>> macToRateSideTime;
  std::map<uint8_t,std::vector<float>> macToRateTop;
  std::map<uint8_t,std::vector<uint64_t>> macToRateTopTime;

  //AA: this list will need to be updated (and perhaps moved to a separate file)
  //  const uint8_t macsWestIn[6]   = {21,19,17,22,23,24};
  //const uint8_t macsWestOut[6]  = {13,11,10,14,15,16};
  //const uint8_t macsNorthIn[4]  = {7,6,3,1};
  //const uint8_t macsNorthOut[4] = {9,8,5,4};

  //TTree * tree;
  //std::map<std::pair<short,short>,TH1F*>> macChanToHist;
  //std::map<short,short> macToIndex;

//data payload
  /* uint8_t mac5; //last 8 bits of FEB mac5 address
   uint16_t flags;
   uint16_t lostcpu;
   uint16_t lostfpga;*/
   //uint32_t ts0;
   //uint32_t ts1;
   //uint16_t adc[32];
   //uint32_t coinc;

//metadata
   /*uint64_t  run_start_time;
   uint64_t  this_poll_start;
   uint64_t  this_poll_end;
   uint64_t  last_poll_start;
   uint64_t  last_poll_end;
   int32_t   system_clock_deviation;
   uint32_t  feb_tree_per_poll;
   uint32_t  feb_event_number;*/

   //information from fragment header
   //uint32_t  sequence_id;
   //uint64_t  fragment_timestamp;
 
};

//Define the constructor
icarus::crt::CrtNoiseMonTool::CrtNoiseMonTool(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{

  //this is how you setup/use the TFileService
  //I do this here in the constructor to setup the histogram just once
  //but you can call the TFileService and make new objects from anywhere
  //art::ServiceHandle<art::TFileService> tfs; //pointer to a file named tfs

  //tree = tfs->make<TTree>("tree", "FEB tree");

//event data
  /*tree->Branch("mac5",        &mac5,          "mac5/b");
  tree->Branch("flags",       &flags,         "flags/s");
  tree->Branch("lostcpu",     &lostcpu,       "lostcpu/s");
  tree->Branch("lostfpga",    &lostcpu,       "lostfpga/s");
  tree->Branch("ts0",         &ts0,           "ts0/i");
  tree->Branch("ts1",         &ts1,           "ts1/i");
  tree->Branch("adc",         &adc,           "adc[32]/s");
  tree->Branch("coinc",       &coinc,         "coinc/i");*/

//metadata
  /*tree->Branch("run_start_time",            &run_start_time,              "run_start_time/l");
  tree->Branch("this_poll_start",           &this_poll_start,             "this_poll_start/l");
  tree->Branch("this_poll_end",             &this_poll_end,               "this_poll_end/l");
  tree->Branch("last_poll_start",           &last_poll_start,             "last_poll_start/l");
  tree->Branch("last_poll_end",             &last_poll_end,               "last_poll_end/l");
  tree->Branch("system_clock_deviation",    &system_clock_deviation,      "system_clock_deviation/I");
  tree->Branch("feb_tree_per_poll",       &feb_tree_per_poll,         "feb_tree_per_poll/i");
  tree->Branch("feb_event_number",          &feb_event_number,            "feb_event_number/i");

  tree->Branch("sequence_id",               &sequence_id,                 "sequence_id/i");
  tree->Branch("fragment_timestamp",        &fragment_timestamp,          "fragment_timestamp/l");*/
}

icarus::crt::CrtNoiseMonTool::~CrtNoiseMonTool()
{
}

bool icarus::crt::CrtNoiseMonTool::IsSideCRT(icarus::crt::BernCRTTranslator & hit) {
  /**
   * Fragment ID described in SBN doc 16111
   **/
  return (hit.fragment_ID & 0x3100) == 0x3100;
}

void icarus::crt::CrtNoiseMonTool::analyze(art::Event const & evt) {

  //WK 09/02/21. Update to BernCRTTranslator in sbndaq_artdaq_core
  std::vector<icarus::crt::BernCRTTranslator> hit_vector;

  auto fragmentHandles = evt.getMany<artdaq::Fragments>();
  for (auto  handle : fragmentHandles) {
    if (!handle.isValid() || handle->size() == 0)
      continue;

    auto this_hit_vector = icarus::crt::BernCRTTranslator::getCRTData(*handle);

    hit_vector.insert(hit_vector.end(),this_hit_vector.begin(),this_hit_vector.end());

  }

  for(auto & hit : hit_vector) {

    uint8_t mac5     = hit.mac5;

    if(hit.flags==7) continue; //AA: 1. is it a correct check? 2. does it even work? :(

    const uint64_t runtime = (hit.this_poll_start-hit.run_start_time) *1.0e-9; //time since run start [s]
    const double   dt      = (hit.this_poll_start-hit.last_poll_start)*1.0e-9; //interval between consecutive polls [s]
    const float    rate    = 1.0e-3*hit.hits_in_poll/dt; //poll-averaged event rate [kHz]

    //avoid processing the same poll multiple times (they all have the same
    //    if(macToRateTime[mac5].size()!=0 && macToRateTime[mac5].back() == runtime) continue;
    if(IsSideCRT(hit)){
      //avoid processing the same poll multiple times (they all have the same
      if(macToRateSideTime[mac5].size()!=0 && macToRateSideTime[mac5].back() == runtime) continue;
      macToRateSideTime[mac5].push_back(runtime);
      macToRateSide[mac5].push_back(rate);
    }else{
      //avoid processing the same poll multiple times (they all have the same
      if(macToRateTopTime[mac5].size()!=0 && macToRateTopTime[mac5].back() == runtime) continue;
      macToRateTopTime[mac5].push_back(runtime);
      macToRateTop[mac5].push_back(rate);
    }
  }
}//end analyze


void icarus::crt::CrtNoiseMonTool::endJob(){

  art::ServiceHandle<art::TFileService> tfs; //pointer to a file named tfs
  /*
  std::string hviewName = "hview_";
  std::string hviewTitle = "noise rates: ";
  char reg ='e';//, lay='e';
  if(macToRate.find(macsWestIn[0])!=macToRate.end()) {
      hviewName+="w_i";
      hviewTitle+="West Inner";
      reg = 'w';
      //lay = 'i';
  }
  else if(macToRate.find(macsWestOut[0])!=macToRate.end()) {
      hviewName+="w_o";
      hviewTitle+="West Outer";
      reg = 'w';
      //lay = 'o';
  }
  else if(macToRate.find(macsNorthIn[0])!=macToRate.end()) {
      hviewName+="n_i";
      hviewTitle+="North Inner";
      reg = 'n';
      //lay = 'i';
  }
  else if(macToRate.find(macsNorthOut[0])!=macToRate.end()) {
      hviewName+="n_o";
      hviewTitle+="North Outer";
      reg = 'n';
      //lay = 'o';
  }
  else{
      std::cout << "no matching mac5s found. reg and lay not set!" << std::endl;
  }
 

  size_t yhigh=0;
  if(reg=='n') yhigh = 2;
  else yhigh = 3;
  TH2F *hview = tfs->make<TH2F>(hviewName.c_str(),hviewTitle.c_str(),2,0,2,yhigh,0,yhigh);
  if(reg=='w'){
      hview->GetXaxis()->SetBinLabel(1,"north");
      hview->GetXaxis()->SetBinLabel(2,"south");
      hview->GetYaxis()->SetBinLabel(1,"pit");
      hview->GetYaxis()->SetBinLabel(2,"mezzpit");
      hview->GetYaxis()->SetBinLabel(3,"mezz");
  }
  if(reg=='n'){
      hview->GetXaxis()->SetBinLabel(1,"east");
      hview->GetXaxis()->SetBinLabel(2,"west");
      hview->GetYaxis()->SetBinLabel(1,"pit");
      hview->GetYaxis()->SetBinLabel(2,"mezz");
  }
  */
  //size_t imac=0;
  for(auto it=macToRateSide.begin(); it!=macToRateSide.end(); it++){

      uint8_t mac = (*it).first;

      std::string hname = "hside_"+std::to_string(mac);
      std::string htitle = "noise rate: mac5 "+std::to_string(mac)+" [Side CRT]";
      TH1F *h = tfs->make<TH1F>(hname.c_str(),htitle.c_str(),4000,0.0,40.0);
      h->GetXaxis()->SetTitle("rate [kHz]");
      h->GetXaxis()->SetTitleSize(0.04);
      h->GetXaxis()->SetLabelSize(0.04);
      h->GetYaxis()->SetLabelSize(0.04);
      h->SetLineWidth(2);

      const size_t n = macToRateSide[mac].size();
      float* rate = new float[n];
      float* time = new float[n];
      float meanRate =0.0;
      for(size_t i=0; i<n; i++){
         rate[i] = macToRateSide[mac][i];
	 time[i] = macToRateSideTime[mac][i];
	 meanRate+=rate[i];

         h->Fill(rate[i]);
      }
      std::cout << "mean rate before division: " << meanRate << std::endl;
      meanRate*=1.0/n;
      std::cout << "mean rate averaged over " << n << " points: " << meanRate << std::endl;
      meanRate*=1.0/30;

      std::cout << "relative (to 30kHz) average "+htitle << " = " << meanRate << std::endl;
 
      std::string gname = "gside_"+std::to_string(mac);
      std::string gtitle = "noise rate: mac5 "+std::to_string(mac)+" [Side CRT]";
      TGraph *g = tfs->make<TGraph>(n,time,rate);
      g->SetName(gname.c_str());
      g->SetTitle(gtitle.c_str());
      g->GetXaxis()->SetTitle("run time [s]");
      g->GetYaxis()->SetTitle("rate [kHz]");
      g->GetXaxis()->SetTitleSize(0.04);
      g->GetYaxis()->SetTitleSize(0.04);
      g->GetXaxis()->SetLabelSize(0.04);
      g->GetYaxis()->SetLabelSize(0.04);

      g->SetMarkerStyle(8);

      h->Write();
      g->Write();

      delete[] rate;
      delete[] time;
  }

 for(auto it=macToRateTop.begin(); it!=macToRateTop.end(); it++){

      uint8_t mac = (*it).first;

      std::string hname = "htop_"+std::to_string(mac);
      std::string htitle = "noise rate: mac5 "+std::to_string(mac)+" [Top CRT]";
      TH1F *h = tfs->make<TH1F>(hname.c_str(),htitle.c_str(),4000,0.0,40.0);
      h->GetXaxis()->SetTitle("rate [kHz]");
      h->GetXaxis()->SetTitleSize(0.04);
      h->GetXaxis()->SetLabelSize(0.04);
      h->GetYaxis()->SetLabelSize(0.04);
      h->SetLineWidth(2);

      const size_t n = macToRateTop[mac].size();
      float* rate = new float[n];
      float* time = new float[n];
      float meanRate =0.0;
      for(size_t i=0; i<n; i++){
         rate[i] = macToRateTop[mac][i];
	 time[i] = macToRateTopTime[mac][i];
	 meanRate+=rate[i];

         h->Fill(rate[i]);
      }
      std::cout << "mean rate before division: " << meanRate << std::endl;
      meanRate*=1.0/n;
      std::cout << "mean rate averaged over " << n << " points: " << meanRate << std::endl;
      meanRate*=1.0/30;

      std::cout << "relative (to 30kHz) average "+htitle << " = " << meanRate << std::endl;
 
      std::string gname = "gtop_"+std::to_string(mac);
      std::string gtitle = "noise rate: mac5 "+std::to_string(mac)+" [Top CRT]";
      TGraph *g = tfs->make<TGraph>(n,time,rate);
      g->SetName(gname.c_str());
      g->SetTitle(gtitle.c_str());
      g->GetXaxis()->SetTitle("run time [s]");
      g->GetYaxis()->SetTitle("rate [kHz]");
      g->GetXaxis()->SetTitleSize(0.04);
      g->GetYaxis()->SetTitleSize(0.04);
      g->GetXaxis()->SetLabelSize(0.04);
      g->GetYaxis()->SetLabelSize(0.04);

      g->SetMarkerStyle(8);

      h->Write();
      g->Write();

      delete[] rate;
      delete[] time;
  }
  
  // hview->Write();

} //endJob

DEFINE_ART_MODULE(icarus::crt::CrtNoiseMonTool)
//this is where the name is specified
