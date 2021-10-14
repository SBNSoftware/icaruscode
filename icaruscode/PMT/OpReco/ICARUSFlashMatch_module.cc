////////////////////////////////////////////////////////////////////////
// Class:       ICARUSFlashMatch
// Plugin Type: analyzer (art v3_02_06)
// File:        ICARUSFlashMatch_module.cc
//
// Generated at Fri Oct 11 12:08:44 2019 by Laura Domine using cetskelgen
// from cetlib version v3_07_02.
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

#include <set>
#include <vector>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashMatchManager.h"
#include "sbncode/OpT0Finder/flashmatch/Algorithms/LightPath.h"
#include "sbncode/OpT0Finder/flashmatch/GeoAlgo/GeoAlgo.h"
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

class ICARUSFlashMatch;


class ICARUSFlashMatch : public art::EDAnalyzer {
public:
  explicit ICARUSFlashMatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSFlashMatch(ICARUSFlashMatch const&) = delete;
  ICARUSFlashMatch(ICARUSFlashMatch&&) = delete;
  ICARUSFlashMatch& operator=(ICARUSFlashMatch const&) = delete;
  ICARUSFlashMatch& operator=(ICARUSFlashMatch&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
	void beginJob() override;
	void endJob() override;

private:

	void
	CreateQClusterFromMC(
		flashmatch::QClusterArray_t& result,
		const std::vector<simb::MCParticle>& part_v,
		const std::vector<sim::SimEnergyDeposit>& edep_v) const;

	void
	QClusterFromMCParticle(
		flashmatch::QCluster_t& destination,
		const simb::MCParticle& part) const;

	flashmatch::QCluster_t
	QClusterFromMCParticle(const simb::MCParticle& part) const;

	void
	QClusterFromGeoTrajectory(
		flashmatch::QCluster_t& destination,
		const geoalgo::Trajectory& trj
		) const;


	void AppendFlashToTuple(const flashmatch::Flash_t& flash);
	void AppendQClusterToTuple(const flashmatch::QCluster_t& qcluster);
	void ClearTupleInfo();

  // Declare member data here.
	int _verbose;
	double _cluster_window;
	std::vector<size_t> _cryo_id_v;
	TFile *_f;
    std::string _output_fname;
	std::string _mcparticle_label;
	std::string _simedep_label;
	std::string _opflash_label;
	std::string _light_path_name;

	flashmatch::FlashMatchManager* _mgr;
	flashmatch::LightPath* _light_path;

	TTree* _match_tree;
	int _run, _event;
	int _store_tpc_id, _store_flash_id;
	int _store_touch_match;
	double _store_touch_score, _store_touch_point_x, _store_touch_point_y, _store_touch_point_z;
	double _store_score, _store_tpc_point_x, _store_tpc_point_y, _store_tpc_point_z;
	double _store_true_xmin;
	double _store_hypo_pe_sum, _store_reco_pe_sum;
	std::vector<double> _store_hypo_pe_v, _store_reco_pe_v;
	double _store_duration;
	int _store_num_steps;
	double _store_fit_range_xmin, _store_fit_range_xmax;
	double _store_flash_time, _store_flash_time_width, _store_flash_dt_prev, _store_flash_dt_next;
	double _store_part_time;


	bool _save_tuple_input, _save_tuple_match;

	TTree* _input_tree;
	int _n_tpc, _n_flash;
	std::vector<double> _ftime;
	std::vector<double> _ftime_width;
	std::vector<double> _ftime_dt_prev;
	std::vector<double> _ftime_dt_next;
	std::vector<double> _fpe_v;

	std::vector<double> _qx, _qy, _qz, _qv;
	std::vector<int> _qn;
	std::vector<double> _qtrue_t;
	std::vector<double> _qtrue_x;
};


ICARUSFlashMatch::ICARUSFlashMatch(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, _mgr(nullptr), _light_path(nullptr)
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
	_verbose = p.get<int>("Verbosity",0);
	_output_fname = p.get<std::string>("OutputFileName");
	_mcparticle_label = p.get<std::string>("MCParticleProducer","");
	_simedep_label = p.get<std::string>("SimEnergyDepositProducer","");
	_opflash_label = p.get<std::string>("OpFlashProducer");
	_light_path_name = p.get<std::string>("LightPathName","");
	_cluster_window = p.get<double>("ClusteringPeriod",8);
	_save_tuple_input = p.get<bool>("SaveTupleInput",false);
	_save_tuple_match = p.get<bool>("SaveTupleMatch",false);
	_cryo_id_v.clear();
	_cryo_id_v = p.get<std::vector<size_t> >("EnableCryostats",_cryo_id_v);
	auto& ds = flashmatch::DetectorSpecs::GetME(p.get<fhicl::ParameterSet>("DetectorSpecs"));
	ds.EnableCryostats(_cryo_id_v);
	ds.DumpInfo();

	if(_verbose)
		std::cout<<"Configuring FlashMatchManager..." << std::endl;
	_mgr = new flashmatch::FlashMatchManager();
	_mgr->Configure(p.get<fhicl::ParameterSet>("FlashMatchConfig"));
	if(!_light_path_name.empty())
		_light_path = (flashmatch::LightPath*)(_mgr->GetCustomAlgo(_light_path_name));
}

void ICARUSFlashMatch::beginJob()
{
	if(_save_tuple_input || _save_tuple_match)
		_f = TFile::Open(_output_fname.c_str(), "RECREATE");

	if(_save_tuple_match) {
		std::string name = "matchtree";
		_match_tree = new TTree(name.c_str(), name.c_str());
		_match_tree->Branch("run", &_run, "run/I");
		_match_tree->Branch("event", &_event, "event/I");
		_match_tree->Branch("tpc_id", &_store_tpc_id, "tpc_id/I");
		_match_tree->Branch("flash_id", &_store_flash_id, "flash_id/I");
		_match_tree->Branch("touch_score",&_store_touch_score,"touch_score/D");
		_match_tree->Branch("touch_point_x",&_store_touch_point_x,"touch_point_x/D");
		_match_tree->Branch("touch_point_y",&_store_touch_point_y,"touch_point_y/D");
		_match_tree->Branch("touch_point_z",&_store_touch_point_z,"touch_point_z/D");
		_match_tree->Branch("score",&_store_score,"score/D");
		_match_tree->Branch("tpc_point_x",&_store_tpc_point_x,"tpc_point_x/D");
		_match_tree->Branch("tpc_point_y",&_store_tpc_point_y,"tpc_point_y/D");
		_match_tree->Branch("tpc_point_z",&_store_tpc_point_z,"tpc_point_z/D");
		_match_tree->Branch("true_xmin",&_store_true_xmin,"true_xmin/D");
		_match_tree->Branch("flash_time",&_store_flash_time,"flash_time/D");
		_match_tree->Branch("flash_time_width",&_store_flash_time_width,"flash_time_width/D");
		_match_tree->Branch("flash_dt_prev",&_store_flash_dt_prev,"flash_dt_prev/D");
		_match_tree->Branch("flash_dt_next",&_store_flash_dt_next,"flash_dt_next/D");

		_match_tree->Branch("part_time",&_store_part_time,"part_time/D");
		_match_tree->Branch("hypo_pe_sum",&_store_hypo_pe_sum,"hypo_pe_sum/D");
		_match_tree->Branch("reco_pe_sum",&_store_reco_pe_sum,"reco_pe_sum/D");

		_store_hypo_pe_v.resize(180);
		_store_reco_pe_v.resize(180);
		_match_tree->Branch("hypo_pe_v",&_store_hypo_pe_v);
		_match_tree->Branch("reco_pe_v",&_store_reco_pe_v);

		_match_tree->Branch("duration",&_store_duration,"duration/D");
		_match_tree->Branch("num_steps",&_store_num_steps,"num_steps/I");
		_match_tree->Branch("fit_range_xmin",&_store_fit_range_xmin,"fit_range_xmin/D");
		_match_tree->Branch("fit_range_xmax",&_store_fit_range_xmax,"fit_range_xmax/D");

	}


	if(_save_tuple_input) {
		_input_tree = new TTree("input_tree","input_tree");
		_input_tree->Branch("n_tpc",&_n_tpc,"n_tpc/I");
		_input_tree->Branch("qn",&_qn);
		_input_tree->Branch("qtrue_t",&_qtrue_t);
		_input_tree->Branch("qtrue_x",&_qtrue_x);
		_input_tree->Branch("qx",&_qx);
		_input_tree->Branch("qy",&_qy);
		_input_tree->Branch("qz",&_qz);
		_input_tree->Branch("qv",&_qv);

		_input_tree->Branch("n_flash",&_n_flash,"n_flash/I");
		_input_tree->Branch("ftime",&_ftime);
		_input_tree->Branch("ftime_width",&_ftime_width);
		_input_tree->Branch("ftime_dt_next",&_ftime_dt_next);
		_input_tree->Branch("ftime_dt_prev",&_ftime_dt_prev);
		_input_tree->Branch("fpe",&_fpe_v);
	}
}

void ICARUSFlashMatch::endJob()
{
	if(_save_tuple_input || _save_tuple_match) {
		_f->cd();
		if(_save_tuple_match) _match_tree->Write();
		if(_save_tuple_input) _input_tree->Write();
		_f->Close();
	}
}


void ICARUSFlashMatch::AppendFlashToTuple(const flashmatch::Flash_t& flash)
{
	_n_flash++;
	_ftime.push_back(flash.time);
	_ftime_width.push_back(flash.time_width);
	_ftime_dt_next.push_back(flash.dt_next);
	_ftime_dt_prev.push_back(flash.dt_prev);
	_fpe_v.reserve(_fpe_v.size()+flash.pe_v.size());
	for(auto const& v : flash.pe_v) _fpe_v.push_back(v);
}


void ICARUSFlashMatch::AppendQClusterToTuple(const flashmatch::QCluster_t& qcluster)
{
	_n_tpc++;
	_qx.reserve(_qx.size()+qcluster.size());
	_qy.reserve(_qy.size()+qcluster.size());
	_qz.reserve(_qz.size()+qcluster.size());
	_qv.reserve(_qv.size()+qcluster.size());
	for(auto const& qpt : qcluster) {
		_qx.push_back(qpt.x);
		_qy.push_back(qpt.y);
		_qz.push_back(qpt.z);
		_qv.push_back(qpt.q);
	}
	_qn.push_back(qcluster.size());
	_qtrue_t.push_back(qcluster.time_true);
	_qtrue_x.push_back(qcluster.min_x_true);
}


void ICARUSFlashMatch::ClearTupleInfo()
{
	_n_tpc = 0;
	_n_flash = 0;
	_ftime.clear();
	_ftime_width.clear();
	_ftime_dt_next.clear();
	_ftime_dt_prev.clear();
	_fpe_v.clear();
	_qx.clear();
	_qy.clear();
	_qz.clear();
	_qv.clear();
	_qn.clear();
	_qtrue_t.clear();
	_qtrue_x.clear();
}

void
ICARUSFlashMatch::QClusterFromMCParticle(
	flashmatch::QCluster_t& destination,
	const simb::MCParticle& part) const
{

	const auto& trajectory = part.Trajectory();
	::geoalgo::Trajectory gtrj(trajectory.size(),3);

	for(size_t idx=0; idx<trajectory.size(); ++idx) {
		gtrj[idx][0] = trajectory.X(idx);
		gtrj[idx][1] = trajectory.Y(idx);
		gtrj[idx][2] = trajectory.Z(idx);
	}
	return this->QClusterFromGeoTrajectory(destination, gtrj);
}

flashmatch::QCluster_t
ICARUSFlashMatch::QClusterFromMCParticle(const simb::MCParticle& part) const
{
	flashmatch::QCluster_t result;
	this->QClusterFromMCParticle(result,part);
	return result;
}

void
ICARUSFlashMatch::QClusterFromGeoTrajectory(
	flashmatch::QCluster_t& destination,
	const geoalgo::Trajectory& trj
	) const
{
	if(!_light_path) {
		std::cerr << "LightPath algorithm is not configured!" << std::endl << std::endl << std::flush;
		throw std::exception();
	}

	auto qcls = _light_path->MakeQCluster(trj);
	destination.reserve(destination.size()+qcls.size());
	for(auto& qpt : qcls) destination.emplace_back(qpt);
}


void
ICARUSFlashMatch::CreateQClusterFromMC(
	flashmatch::QClusterArray_t& result,
	const std::vector<simb::MCParticle>& part_v,
	const std::vector<sim::SimEnergyDeposit>& edep_v) const
{

	// Create a list of all particle track IDs by looping over SimEnergyDeposit, estimate QCluster_t size per particle
	std::vector<flashmatch::QCluster_t> qcluster_v;
	std::vector<size_t> qcluster_size_v;
	for(auto const& edep : edep_v) {
		int track_id = std::abs(edep.TrackID());
		if(track_id >= (int)(qcluster_v.size())) {
			qcluster_v.resize(track_id+1);
			qcluster_size_v.resize(track_id+1,0);
		}
		qcluster_size_v[track_id]++;
		double edep_time = edep.T() * 1.e-3;
		if(edep_time < qcluster_v[track_id].time) {
			qcluster_v[track_id].time = edep_time;
			qcluster_v[track_id].time_true = edep_time;
		}
	}

	// Search for primaries + record timings
	std::map<double,int> start_time_m;
	std::vector<int> pdg_v(qcluster_v.size(),0);
	for(auto const& p : part_v) {
		int track_id = std::abs(p.TrackId());
		double p_time = p.T() * 1.e-3;
		if(track_id >= (int)(qcluster_v.size())) {
			qcluster_v.resize(track_id+1);
			qcluster_size_v.resize(track_id+1,0);
			pdg_v.resize(track_id+1,0);
		}
		if(p_time < qcluster_v[track_id].time) {
			qcluster_v[track_id].time = p_time;
			qcluster_v[track_id].time_true = p_time;
		}
		pdg_v[track_id] = p.PdgCode();
		if(p.Mother() != 0) continue;
		if(_verbose) std::cout << "Found a primary Track ID " << track_id
			<< " PDG " << pdg_v[track_id] << " at T " << p_time << std::endl;
		start_time_m.insert({p_time,track_id});
	}

	// Create a time-ordered array of events (QCluster_t)
	result.clear();
	result.reserve(start_time_m.size());
	std::vector<int> primary_track_id_v;
	for(auto const& keyval : start_time_m) {
		flashmatch::QCluster_t qc;
		qc.time = keyval.first;
		qc.time_true = keyval.first;
		result.emplace_back(qc);
		primary_track_id_v.push_back(keyval.second);
	}

	// Check which particle should be clustered into result array. 
	// The result array index is recorded per particle in cluster_id_v (-1 means should not be clustered)
	std::vector<int> cluster_id_v(qcluster_v.size(),-1);
	for(size_t track_id=0; track_id<qcluster_v.size(); ++track_id) {
		auto const& qcluster = qcluster_v[track_id];
		if(qcluster.time_true == flashmatch::kINVALID_DOUBLE) continue;

		auto iter = start_time_m.upper_bound(qcluster.time);
		if(iter == start_time_m.begin()) {
			std::cout<<"Ignoring a particle (too early at " << qcluster.time 
			<< " the first time " << start_time_m.begin()->first
			<< ") TrackID " << track_id << " PDG " << pdg_v[track_id] << std::endl;
			cluster_id_v[track_id]=-1;
			continue;
		}
		--iter;
		if( (qcluster.time - (iter->first)) > _cluster_window ) {
			std::cout<<"Ignoring a particle (too late at " << qcluster.time 
			<< " the last time " << std::prev(start_time_m.end())->first 
			<< ") TrackID " << track_id << " PDG " << pdg_v[track_id] << std::endl;
			cluster_id_v[track_id]=-1;
			continue;
		}
		// Identified a particle to be clustered
		auto qc_index =  std::distance(start_time_m.begin(), iter);
		cluster_id_v[track_id]=qc_index;
	}

	// Reserve qcluster size
	for(size_t i=0; i<qcluster_size_v.size(); ++i) {
		if(cluster_id_v[i]<0) continue;
		qcluster_v[i].reserve(qcluster_size_v[i]);
	}

	// Fill QCluster
	auto light_yield = flashmatch::DetectorSpecs::GetME().LightYield();
	for(auto const& edep : edep_v) {
		int track_id = std::abs(edep.TrackID());
		int cluster_id = cluster_id_v[track_id];
		if(cluster_id<0) continue;
		flashmatch::QPoint_t qpt(edep.X(),edep.Y(),edep.Z(),edep.Energy() * light_yield);
		qcluster_v[track_id].emplace_back(qpt);
	}

	// If LightPath is enabled, loop over particles and create QCluster from Trajectory and replace with the existing one
	if(_light_path) {
		for(auto const& p : part_v) {
			int track_id = std::abs(p.TrackId());
			if(cluster_id_v[track_id]<0) continue;
			if(std::abs(p.PdgCode()) == 11 || std::abs(p.PdgCode()) == 22 || p.Trajectory().size()<2)
				continue;

			qcluster_v[track_id].clear();
			const auto& trajectory = p.Trajectory();
			geoalgo::Trajectory geo_trj(trajectory.size(),3);
			for (size_t pt_idx = 0; pt_idx < trajectory.size(); ++pt_idx) {
				geo_trj[pt_idx][0] = trajectory.X(pt_idx);
				geo_trj[pt_idx][1] = trajectory.Y(pt_idx);
				geo_trj[pt_idx][2] = trajectory.Z(pt_idx);
			}
			//std::cout << "Track ID " << track_id << " ... " << qcluster_v[track_id].size() << " pts " << std::endl;
			this->QClusterFromGeoTrajectory(qcluster_v[track_id],geo_trj);
			/*
			std::cout << " added MCTrajectory ... "
				  << "(" << geo_trj[0][0] << "," << geo_trj[0][1] << "," << geo_trj[0][2] 
				  << ") => "
				  << "(" << geo_trj.back()[0] << "," << geo_trj.back()[1] << "," << geo_trj.back()[2]
				  << ") ... " << qcluster_v[track_id].size() << " pts " << std::endl;
			*/
		}
	}

	// Merge QCluster from each particle
	for(size_t track_id=0; track_id<cluster_id_v.size(); ++track_id){
		auto const& cluster_id = cluster_id_v[track_id];
		if(cluster_id<0) continue;

		//std::cout << "Merging QCluster: track ID " << track_id << " => " << primary_track_id_v[cluster_id] << std::endl;

		auto& qcluster = qcluster_v[track_id];
		result[cluster_id].reserve(result[cluster_id].size()+qcluster.size());
		for(auto& qpt : qcluster)
			result[cluster_id].emplace_back(qpt);
	}

	// Set true min X per QCluster, shift time unit from ns to us
	// Only use points that are within the active volume.
	//auto const& bbox = flashmatch::DetectorSpecs::GetME().ActiveVolume();
	auto const& bbox = flashmatch::DetectorSpecs::GetME().ActiveVolume();
	auto const& min_pt = bbox.Min();
	auto const& max_pt = bbox.Max();
	for(size_t idx=0; idx<result.size(); ++idx) {
		auto& qcluster = result[idx];
		qcluster.time = qcluster.time;
		qcluster.time_true = qcluster.time_true;

		if(_verbose) {
			std::cout << std::endl 
				  << "QCluster " << idx << " Track ID " << primary_track_id_v[idx] 
				  << " number of points " << qcluster.size() 
				  << " time " << qcluster.time_true
				  << " qsum " << qcluster.sum() 
				  << std::endl
				  << "  Initial range X    : " << qcluster.min_x() << " => " << qcluster.max_x() << std::endl; 
		}

		// Drop points outside the active volume
		std::cout << std::endl << "Dropping points for QCluster " << idx << " Track ID " << primary_track_id_v[idx]
			  << " ... before " << qcluster.size() << std::flush;
		qcluster.drop(min_pt[0],min_pt[1],min_pt[2],max_pt[0],max_pt[1],max_pt[2]);
		std::cout << " after " << qcluster.size() << std::endl;
		qcluster.min_x_true = qcluster.min_x();

		// Add an artificial x position shift w.r.t. beam time (=0 in case of icarus)
		double min_pt_x = qcluster.min_x();
		double max_pt_x = qcluster.max_x();
		if(_verbose)
			std::cout << "  After truncation X : " << min_pt_x << " => " << max_pt_x << std::endl; 
		double shift_x = qcluster.time_true * flashmatch::DetectorSpecs::GetME().DriftVelocity();
		// Shift is "+" if drifted toward the tpc plane at lower x position
		if( (min_pt_x - bbox.Min()[0]) > (bbox.Max()[0] - max_pt_x) ) shift_x *= -1.;
		qcluster += shift_x;

		if(_verbose)
			std::cout << "  After time-shift X : " << qcluster.min_x() << " => " << qcluster.max_x() << std::endl; 
	}
}


void ICARUSFlashMatch::analyze(art::Event const& e)
{
    // Implementation of required member function here.
	_event = e.id().event();
	_run   = e.id().run();

	// QCluster array to be used for flash matching
	flashmatch::QClusterArray_t tpc_event_v;
	flashmatch::FlashArray_t pmt_event_v;

	//
	// Create TPC objects
	//
	if(!_mcparticle_label.empty()) {
		// Create QClusterArray_t (TPC event list) from MC objects (MCParticle and SimEnergyDeposit)
		if(_simedep_label.empty()) {
			std::cerr << "MCParticle and SimEnergyDeposit producer labels must be both provided or both empty!"
			<< std::endl << std::flush;
			throw std::exception();
		}
		art::Handle<std::vector<simb::MCParticle> > mcp_h;
		art::Handle<std::vector<sim::SimEnergyDeposit> > edep_h;
		e.getByLabel(_mcparticle_label,mcp_h);
		e.getByLabel(_simedep_label,edep_h);
		this->CreateQClusterFromMC(tpc_event_v, *mcp_h, *edep_h);
	}else{
		// Create QClusterArray_t (TPC event list) from reconstructed objects
		std::cerr << "Reconstruction not implemented! " << std::endl << std::flush;
		throw std::exception();
	}

	//
	// Create PMT objects
	//
	art::Handle<std::vector<recob::OpFlash> > flash_h;
	e.getByLabel(_opflash_label,flash_h);
	pmt_event_v.resize(flash_h->size());
	std::set<double> flash_time_s;
	for(size_t idx=0; idx<flash_h->size(); ++idx) {
		auto const& flash = (*flash_h)[idx];
		auto& pmt_event = pmt_event_v[idx];
		//pmt_event.x = flash.XCenter();
		pmt_event.y = flash.YCenter();
		pmt_event.z = flash.ZCenter();
		//pmt_event.x_err = flash.XWidth();
		pmt_event.y_err = flash.YWidth();
		pmt_event.z_err = flash.ZWidth();
		pmt_event.idx = idx;
		pmt_event.time = flash.Time();
		pmt_event.time_width = flash.TimeWidth();
		pmt_event.pe_v = flash.PEs();
		pmt_event.pe_err_v.reserve(pmt_event.pe_v.size());
		for(auto const& v : pmt_event.pe_v)
			pmt_event.pe_err_v.push_back(sqrt(v));
		flash_time_s.insert(pmt_event.time);

		if(_verbose) {
			std::cout << std::endl 
			<< "Flash " << idx << " at " << pmt_event.time << " width " << pmt_event.time_width 
			<< " PE min " << pmt_event.min_pe() << " max " << pmt_event.max_pe() << " sum " << pmt_event.TotalPE()
			<< std::endl;
		}
		//std::cout<<idx<<" ... "<<pmt_event.pe_v.size()<< " v.s. " << flashmatch::DetectorSpecs::GetME().NOpDets()<<std::endl;
	}
	for(auto& pmt_event : pmt_event_v) {
		auto iter = flash_time_s.find(pmt_event.time);
		if(iter == flash_time_s.end()) continue;

		if(iter != flash_time_s.begin())
			pmt_event.dt_prev = (*iter) - (*(std::prev(iter)));

		if(std::next(iter) != flash_time_s.end())
			pmt_event.dt_next = (*(std::next(iter))) - (*iter);
	}

	// Prepare flash match manager
	_mgr->Reset();
	this->ClearTupleInfo();
	for(size_t tpc_idx=0; tpc_idx < tpc_event_v.size(); ++tpc_idx)
	{ 
		if(_save_tuple_input) this->AppendQClusterToTuple(tpc_event_v[tpc_idx]);
		_mgr->Emplace(std::move(tpc_event_v[tpc_idx]));
	}
	for(size_t pmt_idx=0; pmt_idx < pmt_event_v.size(); ++pmt_idx)
	{ 
		if(_save_tuple_input) this->AppendFlashToTuple(pmt_event_v[pmt_idx]);
		_mgr->Emplace(std::move(pmt_event_v[pmt_idx]));
	}

	// Run flash matching
	if(_verbose) {
		std::cout << "Running a flash match with " << _mgr->QClusterArray().size() << " TPC objects v.s. " 
		<< _mgr->FlashArray().size() << " PMT objects" << std::endl;
	}
	auto result = _mgr->Match();

	if(_verbose)
		std::cout << "Found " << result.size() << " matches..." << std::endl;

	// Record results
	auto const& qcluster_v = _mgr->QClusterArray();
	auto const& flash_v    = _mgr->FlashArray();

	for(size_t idx=0; idx<result.size(); ++idx) {
		auto const& match = result[idx];
		auto const& qcluster = qcluster_v[match.tpc_id];
		auto const& flash = flash_v[match.flash_id];

		_store_hypo_pe_sum = std::accumulate(match.hypothesis.begin(),match.hypothesis.end(),0.);
		_store_reco_pe_sum = std::accumulate(flash.pe_v.begin(),flash.pe_v.end(),0.);

		if(_verbose) {
			std::cout << "Match " << idx << " ... TPC " << match.tpc_id << " with PMT " << match.flash_id << std::endl
			<< "  Time comparison   : " << qcluster.time_true << " v.s. flash " << flash.time << std::endl
			<< "  Touch match " << match.touch_match << " score " << match.touch_score << std::endl 
			<< "  Flash match score : " << match.score << std::endl 
			<< "  QCluster   X hypo : " << match.tpc_point.x << " v.s. true " << qcluster.min_x_true << std::endl 
			<< "  Flash PE sum hypo : " << _store_hypo_pe_sum << " v.s. reco " << _store_reco_pe_sum << std::endl << std::endl;
		}

		if(_save_tuple_match) {

			_store_tpc_id = match.tpc_id;
			_store_flash_id = match.flash_id;

			_store_touch_match = (int)(match.touch_match);
			_store_touch_score = match.touch_score;
			_store_touch_point_x = match.touch_point.x;
			_store_touch_point_y = match.touch_point.y;
			_store_touch_point_z = match.touch_point.z;

			_store_score = match.score;
			_store_tpc_point_x = match.tpc_point.x;
			_store_tpc_point_y = match.tpc_point.y;
			_store_tpc_point_z = match.tpc_point.z;
			_store_true_xmin = qcluster.min_x_true;
			_store_flash_time = flash.time;
			_store_flash_time_width = flash.time_width;
			_store_flash_dt_prev = flash.dt_prev;
			_store_flash_dt_next = flash.dt_next;
			_store_part_time = qcluster.time_true;
			_store_hypo_pe_v = match.hypothesis;
			_store_reco_pe_v = flash.pe_v;

			_store_duration  = match.duration;
			_store_num_steps = match.num_steps;
			_store_fit_range_xmin = match.minimizer_min_x;
			_store_fit_range_xmax = match.minimizer_max_x;

			_match_tree->Fill();
		}
	}
	if(_save_tuple_input) _input_tree->Fill();
}

DEFINE_ART_MODULE(ICARUSFlashMatch)
