////////////////////////////////////////////////////////////////////////
// Class:       ICARUSOpFlashAna
// Plugin Type: analyzer (art v3_01_01)
// File:        ICARUSOpFlashAna_module.cc
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
#include "lardataobj/RecoBase/OpFlash.h"
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

class ICARUSOpFlashAna;

class ICARUSOpFlashAna : public art::EDAnalyzer {
public:
  explicit ICARUSOpFlashAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSOpFlashAna(ICARUSOpFlashAna const&) = delete;
  ICARUSOpFlashAna(ICARUSOpFlashAna&&) = delete;
  ICARUSOpFlashAna& operator=(ICARUSOpFlashAna const&) = delete;
  ICARUSOpFlashAna& operator=(ICARUSOpFlashAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
private:

  // Declare member data here.
  TFile *_f;
  std::string _output_fname;
  // For data product labels
  std::string _mcflash_label;
	std::string _mctruth_label;
  std::vector<std::string> _flash_label_v;

  // For waveform tree


  // For flash trees
  int _run, _event;
  std::vector<TTree*> _flashtree_v;
  double _time;
  double _pe_sum;
  std::vector<double> _pe_v;
  double _time_true;
  double _pe_sum_true;
  std::vector<double> _pe_true_v;
	double _x;
	double _y;
	double _z;
	double _nphotons;

  // Time period to match reco<=>MC (in micro-second)
  double _match_time_min;
  double _match_time_max;
  double _match_dt;

  // For geometry info
  TTree *_geotree;

};


ICARUSOpFlashAna::ICARUSOpFlashAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, _geotree(nullptr)
{
  _output_fname = p.get<std::string>("OutputFileName"        );
  _mcflash_label  = p.get<std::string>("MCOpFlashProducer"       );
	_mctruth_label = p.get<std::string>("MCTruthProducer");
  _flash_label_v  = p.get<std::vector<std::string> >("OpFlashProducerList");
  _match_time_min = p.get<double>("MatchTimeStart",0.105); // in micro-seconds
  _match_time_max = p.get<double>("MatchTimeEnd",0.120); // in micro-seconds
  _match_dt     = _match_time_max - _match_time_min;
  assert(_match_dt>0);
}

void ICARUSOpFlashAna::beginJob()
{
  _f = TFile::Open(_output_fname.c_str(),"RECREATE");

  for(auto const& label : _flash_label_v) {
    std::string name = label + "_flashtree";
    auto flashtree = new TTree(name.c_str(),name.c_str());
    flashtree->Branch("run",&_run,"run/I");
    flashtree->Branch("event",&_event,"event/I");
    flashtree->Branch("time",&_time,"time/D");
    flashtree->Branch("pe_v",&_pe_v);
    flashtree->Branch("time_true",&_time_true,"time_true/D");
    flashtree->Branch("pe_true_v",&_pe_true_v);
    flashtree->Branch("pe_sum",&_pe_sum,"pe_sum/D");
    flashtree->Branch("pe_sum_true",&_pe_sum_true,"pe_sum_true/D");
    flashtree->Branch("x",&_x,"x/D");
    flashtree->Branch("y",&_y,"y/D");
    flashtree->Branch("z",&_z,"z/D");
    flashtree->Branch("nphotons",&_nphotons,"nphotons/D");
    _flashtree_v.push_back(flashtree);
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
    minX.push_back(tpc.ActiveBoundingBox().MinX());
    minY.push_back(tpc.ActiveBoundingBox().MinY());
    minZ.push_back(tpc.ActiveBoundingBox().MinZ());
    maxX.push_back(tpc.ActiveBoundingBox().MaxX());
    maxY.push_back(tpc.ActiveBoundingBox().MaxY());
    maxZ.push_back(tpc.ActiveBoundingBox().MaxZ());
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

void ICARUSOpFlashAna::endJob()
{
  _f->cd();
  for(auto& ptr : _flashtree_v) { _f->cd(); ptr->Write(); }
  _f->cd(); _geotree->Write();
  if(_f) _f->Close();
}

void ICARUSOpFlashAna::analyze(art::Event const& e)
{

  _event = e.id().event();
  _run   = e.id().run();

	// get MCTruth
	art::Handle< std::vector< simb::MCTruth > > mctruth_h;
	e.getByLabel(_mctruth_label, mctruth_h);
	std::map<double, int> mctruth_db;
	for (size_t mctruth_index = 0; mctruth_index < mctruth_h->size(); ++mctruth_index) {
		auto const& mctruth = (*mctruth_h)[mctruth_index];
		for (int part_idx = 0; part_idx < mctruth.NParticles(); ++part_idx) {
			const simb::MCParticle & particle = mctruth.GetParticle(part_idx);
			//const TLorentzVector& pos = particle.Position();
			//const TLorentzVector& mom = particle.Momentum();
			mctruth_db[particle.T() + _match_time_min] = part_idx; // FIXME assumes mctruth_h->size() == 1 always?
		}
	}

  // get MCOpFlash
  art::Handle< std::vector< recob::OpFlash > > mcflash_h;
  e.getByLabel(_mcflash_label, mcflash_h);
  if(!mcflash_h.isValid()) {
    std::cerr << "Invalid producer for truth recob::OpFlash: " << _mcflash_label << std::endl;
    throw std::exception();
  }
  // Create a "time-map" of MCOpFlash
  // inner map ... key = mcflash timing
  //               value = mcflash location (i.e. array index number)
  std::map<double,int> mcflash_db;
  // fill the map
  //auto const geop = lar::providerFrom<geo::Geometry>();
  for(size_t mcflash_index=0; mcflash_index < mcflash_h->size(); ++mcflash_index) {
    auto const& mcflash = (*mcflash_h)[mcflash_index];
    mcflash_db[mcflash.Time() + _match_time_min] = mcflash_index;
  }
  // now fill opflash trees
  for(size_t label_idx=0; label_idx<_flash_label_v.size(); ++label_idx) {
    // Get data product handle
    auto const& label = _flash_label_v[label_idx];
    auto& flashtree = _flashtree_v[label_idx];
    art::Handle< std::vector< recob::OpFlash > > flash_h;
    e.getByLabel(label,flash_h);
    if(!flash_h.isValid()){
      std::cerr << "Invalid producer for recob::OpFlash: " << label << std::endl;
      throw std::exception();
    }

    // keep the record of which mcflash was used (to store un-tagged mcflash at the end)
    std::vector<bool> mcflash_used(mcflash_h->size(),false);
    // now loop over flash, identify mc flash, fill ttree
    for(auto const& flash : (*flash_h)) {
      // fill simple info
      _time = flash.Time();
      _pe_v = flash.PEs();
      _pe_sum = flash.TotalPE();//std::accumulate(_pe_v.begin(),_pe_v.end());
			// search for corresponding mctruth
			auto low_mct = mctruth_db.lower_bound(_time);
			if (low_mct != mctruth_db.begin()) {
				--low_mct;
				auto const& mctruth = (*mctruth_h).at(0);
				auto const& particle = mctruth.GetParticle((*low_mct).second);
				if ( (particle.T() - (*low_mct).first) < _match_dt) {
					_nphotons = particle.E();
					_x = particle.Vx();
					_y = particle.Vy();
					_z = particle.Vz();
				}
			}
      // search for corresponding mcflash
      auto low = mcflash_db.lower_bound(_time);
      _pe_true_v.resize(_pe_v.size());
      for(auto& pe : _pe_true_v) pe = 0.;
      _time_true = std::numeric_limits<double>::max();
      _pe_sum_true = -1;
      if(low != mcflash_db.begin()) {
	--low;
	// get mc opflash
	auto const& mcflash = (*mcflash_h).at((*low).second);
	// Check if this is in the "match" range
	if( (_time - (*low).first) < _match_dt ) {
	  _pe_true_v = mcflash.PEs();
	  _time_true = mcflash.Time();
	  _pe_sum_true = mcflash.TotalPE();
	  //_pe_sum_true = std::accumulate(_pe_true_v.begin(),_pe_true_v.end());
	  mcflash_used[(*low).second] = true;
		std::cout << mcflash.TotalPE() << " " << std::accumulate(_pe_true_v.begin(), _pe_true_v.end(), 0.) << std::endl;
	}
      }
      flashtree->Fill();
    }
    // now fill mcflash info that was not tagged
    for(size_t mcflash_idx=0; mcflash_idx < mcflash_used.size(); ++mcflash_idx) {
      if(mcflash_used[mcflash_idx]) continue;
      auto const& mcflash = (*mcflash_h)[mcflash_idx];
      _pe_true_v = mcflash.PEs();
      //_pe_sum_true = std::accumulate(_pe_true_v.begin(),_pe_true-v.end());
      _pe_sum_true = mcflash.TotalPE();
      _time_true = mcflash.Time();
      // fill the "reco flash" values with vogus values
      _time = std::numeric_limits<double>::max();
      _pe_v.clear();
      _pe_v.resize(_pe_true_v.size(),0.);
      _pe_sum = -1.;
      flashtree->Fill();
    }

  }

}

DEFINE_ART_MODULE(ICARUSOpFlashAna)
