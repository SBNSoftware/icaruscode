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

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include <TTree.h>
#include <TFile.h>

class ICARUSOpHitAna;

class ICARUSOpHitAna : public art::EDAnalyzer {
public:
  explicit ICARUSOpHitAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSOpHitAna(ICARUSOpHitAna const&) = delete;
  ICARUSOpHitAna(ICARUSOpHitAna&&) = delete;
  ICARUSOpHitAna& operator=(ICARUSOpHitAna const&) = delete;
  ICARUSOpHitAna& operator=(ICARUSOpHitAna&&) = delete;

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
  std::string _mchit_label;
  std::vector<std::string> _hit_label_v;
  
  // For waveform tree
  TTree* _wftree;
  int _run, _event, _ch;
  double _tstart;
  std::vector<float> _wf;
  std::vector<std::string> _wf_label_v;

  // For hit trees
  std::vector<TTree*> _hittree_v;
  double _time;
  double _amp;
  double _area;
  double _pe;
  double _time_true;
  double _pe_true;

  // Time period to match reco<=>MC (in micro-second)
  double _match_dt;

  // For geometry info
  TTree *_geotree;

};


ICARUSOpHitAna::ICARUSOpHitAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, _wftree(nullptr), _geotree(nullptr)
{
  _output_fname = p.get<std::string>("OutputFileName"        );
  _wf_label     = p.get<std::string>("OpDetWaveformProducer" );
  _mchit_label  = p.get<std::string>("MCOpHitProducer"       );
  _hit_label_v  = p.get<std::vector<std::string> >("OpHitProducerList");
  _match_dt     = p.get<double>("MatchDT",0.01); // in micro-seconds
}

void ICARUSOpHitAna::beginJob()
{
  _f = TFile::Open(_output_fname.c_str(),"RECREATE");
  _wftree = nullptr;
  if(!_wf_label.empty()) {
    std::string name = _wf_label + "_wftree";
    _wftree = new TTree(name.c_str(),name.c_str());
    _wftree->Branch("run",&_run,"run/I");
    _wftree->Branch("event",&_event,"event/I");
    _wftree->Branch("ch",&_ch,"ch/I");
    _wftree->Branch("wf",&_wf);
    _wftree->Branch("tstart",&_tstart,"tstart/D");
  }

  for(auto const& label : _hit_label_v) {
    std::string name = label + "_hittree";
    auto hittree = new TTree(name.c_str(),name.c_str());
    hittree->Branch("run",&_run,"run/I");
    hittree->Branch("event",&_event,"event/I");
    hittree->Branch("ch",&_ch,"ch/I");
    hittree->Branch("amp",&_amp,"amp/D");
    hittree->Branch("area",&_area,"area/D");
    hittree->Branch("time",&_time,"time/D");
    hittree->Branch("pe",&_pe,"pe/D");
    hittree->Branch("time_true",&_time_true,"time_true/D");
    hittree->Branch("pe_true",&_pe_true,"pe_true/D");
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

void ICARUSOpHitAna::endJob()
{
  _f->cd();
  if(_wftree) _wftree->Write();
  for(auto& ptr : _hittree_v) { _f->cd(); ptr->Write(); }
  _f->cd(); _geotree->Write();
  if(_f) _f->Close();
}

void ICARUSOpHitAna::analyze(art::Event const& e)
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

  // get MCOpHit
  art::Handle< std::vector< recob::OpHit > > mchit_h;
  e.getByLabel(_mchit_label, mchit_h);
  if(!mchit_h.isValid()) {
    std::cerr << "Invalid producer for truth recob::OpHit: " << _mchit_label << std::endl;
    throw std::exception();
  }
  // Create a "time-map" of MCOpHit
  // outer-array = op channel number
  // inner map ... key = mchit timing
  //               value = mchit location (i.e. array index number)
  std::vector<std::map<double,int> > mchit_db;
  // fill the map
  auto const geop = lar::providerFrom<geo::Geometry>();
  mchit_db.resize(geop->NOpChannels());
  for(size_t mchit_index=0; mchit_index < mchit_h->size(); ++mchit_index) {
    auto const& mchit = (*mchit_h)[mchit_index];
    auto& db = mchit_db.at(mchit.OpChannel());
    db[mchit.PeakTime()] = mchit_index;
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
    // keep the record of which mchit was used (to store un-tagged mchit at the end)
    std::vector<bool> mchit_used(mchit_h->size(),false);
    // now loop over hit, identify mc hit, fill ttree
    for(auto const& hit : (*hit_h)) {
      // fill simple info
      _ch=hit.OpChannel();
      _time = hit.PeakTime();
      _amp  = hit.Amplitude();
      _area = hit.Area();
      _pe   = hit.PE();
      // search for corresponding mchit
      auto const& db = mchit_db.at(_ch);
      auto low = db.lower_bound(_time);
      _pe_true   = -1;
      _time_true = std::numeric_limits<double>::max();
      if(low != db.begin()) {
	--low;
	// get mc ophit
	auto const& mchit = (*hit_h)[(*low).second];
	auto mctime = mchit.PeakTime();
	// Check if this is in the "match" range
	if( (_time - mctime) < _match_dt ) {
	  _pe_true   = mchit.PE();
	  _time_true = mchit.PeakTime();
	  mchit_used[(*low).second] = true;
	}
      }
      hittree->Fill();
    }
    // now fill mchit info that was not tagged
    for(size_t mchit_idx=0; mchit_idx < mchit_used.size(); ++mchit_idx) {
      if(mchit_used[mchit_idx]) continue;
      auto const& mchit = (*mchit_h)[mchit_idx];
      _ch = mchit.OpChannel();
      _pe_true = mchit.PE();
      _time_true = mchit.PeakTime();
      // fill the "reco hit" values with vogus values
      _time = std::numeric_limits<double>::max();
      _amp  = std::numeric_limits<double>::max();
      _area = std::numeric_limits<double>::max();
      _pe   = -1;
      hittree->Fill();
    }

  }
  
}

DEFINE_ART_MODULE(ICARUSOpHitAna)
