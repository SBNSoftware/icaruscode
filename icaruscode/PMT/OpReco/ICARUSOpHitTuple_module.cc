////////////////////////////////////////////////////////////////////////
// Class:       ICARUSOpHitAna
// Plugin Type: analyzer (art v3_01_01)
// File:        ICARUSOpHitAna_module.cc
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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include <TTree.h>
#include <TFile.h>

class ICARUSOpHitTuple;

class ICARUSOpHitTuple : public art::EDAnalyzer {
public:
  explicit ICARUSOpHitTuple(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSOpHitTuple(ICARUSOpHitTuple const&) = delete;
  ICARUSOpHitTuple(ICARUSOpHitTuple&&) = delete;
  ICARUSOpHitTuple& operator=(ICARUSOpHitTuple const&) = delete;
  ICARUSOpHitTuple& operator=(ICARUSOpHitTuple&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
private:

  // Declare member data here.
  TFile *_f;
  std::string _output_fname;
  // For data product labels
  std::string _wf_label;
  std::string _mct_label;
  std::vector<std::string> _hit_label_v;

  // Event wise info
  int _run, _event, _ch;
  double _event_time, _event_x, _event_y, _event_z;
  double _event_dr, _event_dx, _event_dy, _event_dz;
  int _tpc;
  
  // For waveform tree
  TTree* _wftree;
  double _tstart;
  std::vector<float> _wf;
  std::vector<std::string> _wf_label_v;

  // For hit trees
  std::vector<TTree*> _hittree_v;
  double _width;
  double _time;
  double _amp;
  double _area;
  double _pe;

  // For geometry info
  TTree *_geotree;

};


ICARUSOpHitTuple::ICARUSOpHitTuple(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, _wftree(nullptr), _geotree(nullptr)
{
  _output_fname = p.get<std::string>("OutputFileName"        );
  _wf_label     = p.get<std::string>("OpDetWaveformProducer" );
  _mct_label    = p.get<std::string>("MCTruthProducer","generator");
  _hit_label_v  = p.get<std::vector<std::string> >("OpHitProducerList");
}

void ICARUSOpHitTuple::beginJob()
{
  _f = TFile::Open(_output_fname.c_str(),"RECREATE");
  _wftree = nullptr;
  if(!_wf_label.empty()) {
    std::string name = _wf_label + "_wftree";
    _wftree = new TTree(name.c_str(),name.c_str());
    _wftree->Branch("run",&_run,"run/I");
    _wftree->Branch("event",&_event,"event/I");
    _wftree->Branch("event_time",&_event_time,"event_time/D");
    _wftree->Branch("event_x",&_event_x,"event_x/D");
    _wftree->Branch("event_y",&_event_y,"event_y/D");
    _wftree->Branch("event_z",&_event_z,"event_z/D");
    _wftree->Branch("event_dr",&_event_dr,"event_dr/D");
    _wftree->Branch("event_dx",&_event_dx,"event_dx/D");
    _wftree->Branch("event_dy",&_event_dy,"event_dy/D");
    _wftree->Branch("event_dz",&_event_dz,"event_dz/D");
    _wftree->Branch("tpc",&_tpc,"tpc/I");
    _wftree->Branch("ch",&_ch,"ch/I");
    _wftree->Branch("wf",&_wf);
    _wftree->Branch("tstart",&_tstart,"tstart/D");
  }

  for(auto const& label : _hit_label_v) {
    std::string name = label + "_hittree";
    auto hittree = new TTree(name.c_str(),name.c_str());
    hittree->Branch("run",&_run,"run/I");
    hittree->Branch("event",&_event,"event/I");
    hittree->Branch("event_time",&_event_time,"event_time/D");
    hittree->Branch("event_x",&_event_x,"event_x/D");
    hittree->Branch("event_y",&_event_y,"event_y/D");
    hittree->Branch("event_z",&_event_z,"event_z/D");
    hittree->Branch("event_dr",&_event_dr,"event_dr/D");
    hittree->Branch("event_dx",&_event_dx,"event_dx/D");
    hittree->Branch("event_dy",&_event_dy,"event_dy/D");
    hittree->Branch("event_dz",&_event_dz,"event_dz/D");
    hittree->Branch("tpc",&_tpc,"tpc/I");
    hittree->Branch("ch",&_ch,"ch/I");
    hittree->Branch("amp",&_amp,"amp/D");
    hittree->Branch("area",&_area,"area/D");
    hittree->Branch("width",&_width,"width/D");
    hittree->Branch("time",&_time,"time/D");
    hittree->Branch("pe",&_pe,"pe/D");
    _hittree_v.push_back(hittree);
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

void ICARUSOpHitTuple::endJob()
{
  _f->cd();
  if(_wftree) _wftree->Write();
  for(auto& ptr : _hittree_v) { _f->cd(); ptr->Write(); }
  _f->cd(); _geotree->Write();
  if(_f) _f->Close();
}

void ICARUSOpHitTuple::analyze(art::Event const& e)
{

  _event = e.id().event();
  _run   = e.id().run();

  // Fill waveforms if requested
  if(!_wf_label.empty()) {
    art::Handle< std::vector< raw::OpDetWaveform > > wf_h;
    e.getByLabel(_wf_label,wf_h);
    if(!wf_h.isValid()){
      std::cerr << "Invalid producer for raw::OpDetWaveform: " << _wf_label << std::endl;
      throw std::exception();
    }
    for(auto const& wf : (*wf_h)) {
      _wf.resize(wf.size());
      for(size_t i=0; i<_wf.size(); ++i) { _wf.at(i) = wf.at(i); }
      _ch=wf.ChannelNumber();
      _tstart=wf.TimeStamp();
      _wftree->Fill();
    }
  }

  // get MCTruth
  _event_time = std::numeric_limits<double>::max();
  _event_x = _event_y = _event_z = std::numeric_limits<double>::max();
  if(!_mct_label.empty()) {
    art::Handle< std::vector<simb::MCTruth> > mct_h;
    e.getByLabel(_mct_label,mct_h);
    if(!mct_h.isValid()) {
      std::cerr << "Invalid producer for simb::MCTruth: " << _mct_label << std::endl;
      throw std::exception();      
    }
    for(size_t mct_index=0; mct_index<mct_h->size(); ++mct_index) {
      
      auto const& mct = (*mct_h)[mct_index];
      
      for(int i=0; i<mct.NParticles(); ++i) {
	auto const& mcp = mct.GetParticle(i);
	if(mcp.StatusCode() != 1) continue;

	auto const& pos = mcp.Position(0);
	if(_event_time < pos.T()/1000.) continue;
	_event_time = pos.T()/1000.;
	_event_x = pos.X();
	_event_y = pos.Y();
	_event_z = pos.Z();
      }
    }
  }
  
  // set dr, dx, dy, dz
  _event_dr = _event_dx = _event_dy = _event_dz =std::numeric_limits<double>::max();
  _tpc = -1;
  if(_event_time != std::numeric_limits<double>::max()) {
    auto const geop = lar::providerFrom<geo::Geometry>();
    // measure smallest dr to any pmt
    double PMTxyz[3];
    for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
      geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
      double dx = PMTxyz[0] - _event_x;
      double dy = PMTxyz[1] - _event_y;
      double dz = PMTxyz[2] - _event_z;
      double dr = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
      if(_event_dr < dr) continue;
      _event_dr = dr;
      _event_dx = dx;
      _event_dy = dy;
      _event_dz = dz;
    }
  }

  // now fill ophit trees
  for(size_t label_idx=0; label_idx<_hit_label_v.size(); ++label_idx) {
    // Get data product handle
    auto const& label = _hit_label_v[label_idx];
    auto& hittree = _hittree_v[label_idx];
    art::Handle< std::vector< recob::OpHit > > hit_h;
    e.getByLabel(label,hit_h);
    if(!hit_h.isValid()){
      std::cerr << "Invalid producer for recob::OpHit: " << label << std::endl;
      throw std::exception();
    }

    // now loop over hit, identify mc hit, fill ttree
    for(auto const& hit : (*hit_h)) {
      // fill simple info
      _ch=hit.OpChannel();
      _time  = hit.PeakTime();
      _amp   = hit.Amplitude();
      _width = hit.Width();
      _area  = hit.Area();
      _pe    = hit.PE();
      hittree->Fill();
    }

  }
  
}

DEFINE_ART_MODULE(ICARUSOpHitTuple)
