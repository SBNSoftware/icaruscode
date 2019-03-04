////////////////////////////////////////////////////////////////////////
// Class:       ICARUSOpRecoAna
// Plugin Type: analyzer (art v3_01_01)
// File:        ICARUSOpRecoAna_module.cc
//
// Generated at Tue Feb 12 06:49:35 2019 by Kazuhiro Terao using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include <TTree.h>
#include <TFile.h>

class ICARUSOpRecoAna;

class ICARUSOpRecoAna : public art::EDAnalyzer {
public:
  explicit ICARUSOpRecoAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSOpRecoAna(ICARUSOpRecoAna const&) = delete;
  ICARUSOpRecoAna(ICARUSOpRecoAna&&) = delete;
  ICARUSOpRecoAna& operator=(ICARUSOpRecoAna const&) = delete;
  ICARUSOpRecoAna& operator=(ICARUSOpRecoAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
private:

  // Declare member data here.
  TFile* _f;
  std::vector<TTree*> _wftree_v;
  std::vector<TTree*> _hittree_v;
  std::vector<TTree*> _flashtree_v;
  std::vector<float> _wf;
  std::vector<double> _pe_v;
  int _ch;
  double _tstart;
  double _tpeak;
  double _tflash;
  double _amp;
  double _pe;
  double _area;
  std::string _output_fname;
  std::vector<std::string> _wf_label_v;
  std::vector<std::string> _hit_label_v;
  std::vector<std::string> _flash_label_v;
};


ICARUSOpRecoAna::ICARUSOpRecoAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  _output_fname = p.get<std::string>("OutputFileName");
  _wf_label_v  = p.get<std::vector<std::string> >("OpDetWaveformProducerList");
  _hit_label_v = p.get<std::vector<std::string> >("OpHitProducerList");
  _flash_label_v = p.get<std::vector<std::string> >("OpFlashProducerList");
}

void ICARUSOpRecoAna::beginJob()
{
  _f = TFile::Open(_output_fname.c_str(),"RECREATE");
  
  for(auto const& label : _wf_label_v) {
    std::string name = label + "_wftree";
    auto wftree = new TTree(name.c_str(),name.c_str());
    wftree->Branch("wf",&_wf);
    wftree->Branch("ch",&_ch,"ch/I");
    wftree->Branch("tstart",&_tstart,"tstart/D");
    _wftree_v.push_back(wftree);
  }
  for(auto const& label : _hit_label_v) {
    std::string name = label + "_hittree";
    auto hittree = new TTree(name.c_str(),name.c_str());
    hittree->Branch("ch",&_ch,"ch/I");
    hittree->Branch("tpeak",&_tpeak,"tpeak/D");
    hittree->Branch("amp",&_amp,"amp/D");
    hittree->Branch("area",&_area,"area/D");
    hittree->Branch("pe",&_pe,"pe/D");
    _hittree_v.push_back(hittree);
  }
  for(auto const& label : _flash_label_v) {
    std::string name = label + "_flashtree";
    auto flashtree = new TTree(name.c_str(),name.c_str());
    flashtree->Branch("pe",&_pe,"pe/D");
    flashtree->Branch("time",&_tflash,"time/D");
    flashtree->Branch("pe_v", &_pe_v);
    _flashtree_v.push_back(flashtree);
  }
}

void ICARUSOpRecoAna::endJob()
{
  for(auto& ptr : _wftree_v) { _f->cd(); ptr->Write(); }
  for(auto& ptr : _hittree_v) { _f->cd(); ptr->Write(); }
  for(auto& ptr : _flashtree_v) { _f->cd(); ptr->Write(); }
  if(_f) _f->Close();
}

void ICARUSOpRecoAna::analyze(art::Event const& e)
{
  for(size_t label_idx=0; label_idx<_wf_label_v.size(); ++label_idx) {
    auto const& label = _wf_label_v[label_idx];
    auto& wftree = _wftree_v[label_idx];
    art::Handle< std::vector< raw::OpDetWaveform > > wf_h;
    e.getByLabel(label,wf_h);
    if(!wf_h.isValid()){
      std::cerr << "Invalid producer for raw::OpDetWaveform: " << label << std::endl;
      throw std::exception();
    }
    for(auto const& wf : (*wf_h)) {
      _wf.resize(wf.size());
      for(size_t i=0; i<_wf.size(); ++i) { _wf.at(i) = wf.at(i); }
      _ch=wf.ChannelNumber();
      _tstart=wf.TimeStamp();
      wftree->Fill();
    }
  }

  for(size_t label_idx=0; label_idx<_hit_label_v.size(); ++label_idx) {
    auto const& label = _hit_label_v[label_idx];
    auto& hittree = _hittree_v[label_idx];
    art::Handle< std::vector< recob::OpHit > > hit_h;
    e.getByLabel(label,hit_h);
    if(!hit_h.isValid()){
      std::cerr << "Invalid producer for recob::OpHit: " << label << std::endl;
      throw std::exception();
    }

    for(auto const& hit : (*hit_h)) {
      _ch=hit.OpChannel();
      _tpeak = hit.PeakTime();
      _amp   = hit.Amplitude();
      _area  = hit.Area();
      _pe    = hit.PE();
      hittree->Fill();
    }
  }

  for(size_t label_idx=0; label_idx<_flash_label_v.size(); ++label_idx) {
    auto const& label = _flash_label_v[label_idx];
    auto& flashtree = _flashtree_v[label_idx];
    art::Handle< std::vector< recob::OpFlash > > flash_h;
    e.getByLabel(label,flash_h);
    if(!flash_h.isValid()){
      std::cerr << "Invalid producer for recob::OpFlash: " << label << std::endl;
      throw std::exception();
    }

    for(auto const& flash : (*flash_h)) {
      _tflash = flash.Time();
      _pe     = flash.TotalPE();
      _pe_v   = flash.PEs();
      flashtree->Fill();
    }
  }
  
}

DEFINE_ART_MODULE(ICARUSOpRecoAna)
