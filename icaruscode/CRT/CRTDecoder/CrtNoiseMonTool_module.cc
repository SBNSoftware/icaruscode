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
  void analyze_fragment(artdaq::Fragment & frag);
  
  std::map<uint8_t,std::vector<float>> macToRate;
  std::map<uint8_t,std::vector<uint64_t>> macToRateTime;

  //AA: this list will need to be updated (and perhaps moved to a separate file)
  const uint8_t macsWestIn[6]   = {21,19,17,22,23,24};
  const uint8_t macsWestOut[6]  = {13,11,10,14,15,16};
  const uint8_t macsNorthIn[4]  = {7,6,3,1};
  const uint8_t macsNorthOut[4] = {9,8,5,4};

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

void icarus::crt::CrtNoiseMonTool::analyze_fragment(artdaq::Fragment & frag) {

  //this applies the 'overlay' to the fragment, to tell it to treat it like a BernCRT fragment
  sbndaq::BernCRTFragment bern_fragment(frag);

  //event data
  sbndaq::BernCRTEvent const* bevt = bern_fragment.eventdata();

  uint8_t mac5     = bevt->MAC5();

  if(bevt->flags==7) return; //AA: 1. is it a correct check? 2. does it even work? :(

  //metadata
  const sbndaq::BernCRTFragmentMetadata* md = bern_fragment.metadata();

  const uint64_t runtime = (md->this_poll_start()-md->run_start_time())*1.0e-9; //time since run start [s]
  const double dt = (md->this_poll_start() - md->last_poll_start())*1.0e-9; //interval between consecutive events [s]
  const float rate = 1.0e-3*md->feb_events_per_poll()/dt; //poll averaged event rate [kHz]

  if(macToRateTime[mac5].size()!=0 && macToRateTime[mac5].back() == runtime) return; //avoid processing the same poll multiple times (they all have the same 
  macToRateTime[mac5].push_back(runtime);
  macToRate[mac5].push_back(rate);


}//end analyze


void icarus::crt::CrtNoiseMonTool::endJob(){

  art::ServiceHandle<art::TFileService> tfs; //pointer to a file named tfs

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

  //size_t imac=0;
  for(auto it=macToRate.begin(); it!=macToRate.end(); it++){

      uint8_t mac = (*it).first;

      std::string hname = "h_"+std::to_string(mac);
      std::string htitle = "noise rate: mac5 "+std::to_string(mac);
      TH1F *h = tfs->make<TH1F>(hname.c_str(),htitle.c_str(),3000,0.0,30.0);
      h->GetXaxis()->SetTitle("rate [kHz]");
      h->GetXaxis()->SetTitleSize(0.04);
      h->GetXaxis()->SetLabelSize(0.04);
      h->GetYaxis()->SetLabelSize(0.04);
      h->SetLineWidth(2);

      const size_t n = macToRate[mac].size();
      float* rate = new float[n];
      float* time = new float[n];
      float meanRate =0.0;
      for(size_t i=0; i<n; i++){
         rate[i] = macToRate[mac][i];
	 time[i] = macToRateTime[mac][i];
	 meanRate+=rate[i];

         h->Fill(rate[i]);
      }
      std::cout << "mean rate before division: " << meanRate << std::endl;
      meanRate*=1.0/n;
      std::cout << "mean rate averaged over " << n << " points: " << meanRate << std::endl;
      meanRate*=1.0/30;

      std::cout << "relative (to 30kHz) average "+htitle << " = " << meanRate << std::endl;

      if(reg=='w'){
          for(size_t i=0; i<6; i++){
                  if(mac==macsWestIn[i]||mac==macsWestOut[i]) {
                      hview->SetBinContent(i/3+1,i%3+1,meanRate);
                      break;
                  }
          }
      }

      if(reg=='n'){
          for(size_t i=0; i<4; i++){
                  if(mac==macsNorthIn[i]||mac==macsNorthOut[i]) {
                      hview->SetBinContent(i/2+1,i%2+1,meanRate);
                      break;
                  }
          }
      }

      std::string gname = "g_"+std::to_string(mac);
      std::string gtitle = "noise rate: mac5 "+std::to_string(mac);
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
  
  hview->Write();

} //analyze_fragment

void icarus::crt::CrtNoiseMonTool::analyze(art::Event const & evt)
{
//  TLOG(TLVL_INFO)<<" Processing event "<<evt.event();

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

DEFINE_ART_MODULE(icarus::crt::CrtNoiseMonTool)
//this is where the name is specified
