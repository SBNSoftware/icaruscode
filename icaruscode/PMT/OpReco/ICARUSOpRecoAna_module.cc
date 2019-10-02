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

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
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
  TFile *_f;
  TTree *_geotree, *_simedep_tree;
  std::vector<TTree*> _wftree_v;
  std::vector<TTree*> _hittree_v;
  std::vector<TTree*> _flashtree_v;
  std::vector<float> _wf;
  std::vector<double> _pe_v;
  int _run, _event;
  int _ch;
  double _tstart;
  double _tpeak;
  double _tflash;
  double _amp;
  double _pe;
  double _area;
  std::string _output_fname;
  std::string _simedep_label;
  std::vector<std::string> _wf_label_v;
  std::vector<std::string> _hit_label_v;
  std::vector<std::string> _flash_label_v;
  std::vector<double> _edep_x, _edep_y, _edep_z, _edep_e, _edep_t;
};


ICARUSOpRecoAna::ICARUSOpRecoAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, _geotree(nullptr), _simedep_tree(nullptr)
{
  _output_fname = p.get<std::string>("OutputFileName");
  _simedep_label = p.get<std::string>("SimEnergyDepositProducer","");
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
    wftree->Branch("run",&_run,"run/I");
    wftree->Branch("event",&_event,"event/I");
    wftree->Branch("wf",&_wf);
    wftree->Branch("ch",&_ch,"ch/I");
    wftree->Branch("tstart",&_tstart,"tstart/D");
    _wftree_v.push_back(wftree);
  }
  for(auto const& label : _hit_label_v) {
    std::string name = label + "_hittree";
    auto hittree = new TTree(name.c_str(),name.c_str());
    hittree->Branch("run",&_run,"run/I");
    hittree->Branch("event",&_event,"event/I");
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
    flashtree->Branch("run",&_run,"run/I");
    flashtree->Branch("event",&_event,"event/I");
    flashtree->Branch("pe",&_pe,"pe/D");
    flashtree->Branch("time",&_tflash,"time/D");
    flashtree->Branch("pe_v", &_pe_v);
    _flashtree_v.push_back(flashtree);
  }
  if(!_simedep_label.empty()) {
    _simedep_tree = new TTree("edep_tree","edep_tree");
    _simedep_tree->Branch("run",&_run,"run/I");
    _simedep_tree->Branch("event",&_event,"event/I");
    _simedep_tree->Branch("edep_x",&_edep_x);
    _simedep_tree->Branch("edep_y",&_edep_y);
    _simedep_tree->Branch("edep_z",&_edep_z);
    _simedep_tree->Branch("edep_e",&_edep_e);
    _simedep_tree->Branch("edep_t",&_edep_t);
  }

  _geotree = new TTree("geotree","geotree");
  std::vector<double> pmtX, pmtY, pmtZ;
  std::vector<double> minX, minY, minZ;
  std::vector<double> maxX, maxY, maxZ;
  auto const geop = lar::providerFrom<geo::Geometry>();
  double PMTxyz[3];
  for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
    geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
    pmtX.push_back(PMTxyz[0]);
    pmtY.push_back(PMTxyz[1]);
    pmtZ.push_back(PMTxyz[2]);
  }
  for(auto iter=geop->begin_TPC(); iter!=geop->end_TPC(); ++iter) {
    auto const& tpc = (*iter);
    minX.push_back(tpc.BoundingBox().MinX());
    minY.push_back(tpc.BoundingBox().MinY());
    minZ.push_back(tpc.BoundingBox().MinZ());
    maxX.push_back(tpc.BoundingBox().MaxX());
    maxY.push_back(tpc.BoundingBox().MaxY());
    maxZ.push_back(tpc.BoundingBox().MaxZ());
  }
  _geotree->Branch("pmtX",&pmtX);
  _geotree->Branch("pmtY",&pmtY);
  _geotree->Branch("pmtZ",&pmtZ);
  _geotree->Branch("minX",&minX);
  _geotree->Branch("minY",&minY);
  _geotree->Branch("minZ",&minZ);
  _geotree->Branch("maxX",&maxX);
  _geotree->Branch("maxY",&maxY);
  _geotree->Branch("maxZ",&maxZ);
  _geotree->Fill();
}

void ICARUSOpRecoAna::endJob()
{
  _f->cd();
  for(auto& ptr : _wftree_v) { _f->cd(); ptr->Write(); }
  for(auto& ptr : _hittree_v) { _f->cd(); ptr->Write(); }
  for(auto& ptr : _flashtree_v) { _f->cd(); ptr->Write(); }
  _f->cd(); _geotree->Write();
  if(_simedep_tree) {_f->cd(); _simedep_tree->Write();}
  if(_f) _f->Close();
}

void ICARUSOpRecoAna::analyze(art::Event const& e)
{

  _event = e.id().event();
  _run   = e.id().run();

  if(!_simedep_label.empty()) {
    art::Handle< std::vector<sim::SimEnergyDeposit> > edep_h;
    e.getByLabel(_simedep_label, edep_h);
    if(!edep_h.isValid()){
      std::cerr << "Invalid producer for sim::SimEnergyDeposit: " << _simedep_label << std::endl;
      throw std::exception();
    }
    _edep_x.resize(edep_h->size(),0.);
    _edep_y.resize(edep_h->size(),0.);
    _edep_z.resize(edep_h->size(),0.);
    _edep_e.resize(edep_h->size(),0.);
    _edep_t.resize(edep_h->size(),0.);
    for(size_t i=0; i<edep_h->size(); ++i) {
      auto const& edep = (*edep_h)[i];
      _edep_x[i] = edep.X();
      _edep_y[i] = edep.Y();
      _edep_z[i] = edep.Z();
      _edep_t[i] = edep.T();
      _edep_e[i] = edep.Energy();
    }
    _simedep_tree->Fill();
  }

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
