////////////////////////////////////////////////////////////////////////
// Class:       ICARUSParticleAna
// Plugin Type: analyzer (art v3_02_06)
// File:        ICARUSParticleAna_module.cc
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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
//#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

class ICARUSParticleAna;


class ICARUSParticleAna : public art::EDAnalyzer {
public:
  explicit ICARUSParticleAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSParticleAna(ICARUSParticleAna const&) = delete;
  ICARUSParticleAna(ICARUSParticleAna&&) = delete;
  ICARUSParticleAna& operator=(ICARUSParticleAna const&) = delete;
  ICARUSParticleAna& operator=(ICARUSParticleAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
	void beginJob() override;
	void endJob() override;

private:

  // Declare member data here.
	TFile *_f;
  std::string _output_fname;
	std::string _particle_label;
	std::string _trajectory_label;

	TTree* _particletree;
	int _run, _event;
	int _pdg_code;
	double _x, _y, _z;
	double _time;
	double _energy;
	double _momentum;
	int _track_id;
	std::vector<double> _x_v, _y_v, _z_v, _time_v, _energy_v;
};


ICARUSParticleAna::ICARUSParticleAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
	_output_fname = p.get<std::string>("OutputFileName");
	_particle_label = p.get<std::string>("ParticleProducer");
	_trajectory_label = p.get<std::string>("TrajectoryProducer");
}

void ICARUSParticleAna::beginJob()
{
	_f = TFile::Open(_output_fname.c_str(), "RECREATE");
	std::string name = _particle_label + "_particletree";
	_particletree = new TTree(name.c_str(), name.c_str());
  _particletree->Branch("run", &_run, "run/I");
	_particletree->Branch("event", &_event, "event/I");
	_particletree->Branch("pdg_code", &_pdg_code, "pdg_code/I");
	_particletree->Branch("x", &_x, "x/D");
	_particletree->Branch("y", &_y, "y/D");
	_particletree->Branch("z", &_z, "z/D");
	_particletree->Branch("time", &_time, "time/D");
	_particletree->Branch("energy", &_energy, "energy/D");
	_particletree->Branch("momentum", &_momentum, "momentum/D");
	_particletree->Branch("track_id", &_track_id, "track_id/I");
	_particletree->Branch("x_v", &_x_v);
	_particletree->Branch("y_v", &_y_v);
	_particletree->Branch("z_v", &_z_v);
	_particletree->Branch("time_v", &_time_v);
	_particletree->Branch("energy_v", &_energy_v);
}

void ICARUSParticleAna::endJob()
{
	_f->cd();
	_particletree->Write();
	if (_f) _f->Close();
}

void ICARUSParticleAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.
	_event = e.id().event();
	_run   = e.id().run();

	art::Handle<std::vector< simb::MCTruth > > mctruth_h;
	e.getByLabel(_particle_label, mctruth_h);
  /*for (size_t idx = 0; idx < mctruth_h->size(); ++idx) {
		auto const& mctruth = (*mctruth_h)[idx];
		for (int part_idx = 0; part_idx < mctruth.NParticles(); ++part_idx) {
			const auto& particle = mctruth.GetParticle(part_idx);
			_pdg_code = particle.PdgCode();
			_x = particle.Vx();
			_y = particle.Vy();
			_z = particle.Vz();
			_time = particle.T();
			_energy = particle.E();
			_momentum = particle.P();
			_track_id = particle.TrackId();

			const auto& trajectory = particle.Trajectory();
			_x_v.resize(trajectory.size(), 0.);
			_y_v.resize(trajectory.size(), 0.);
			_z_v.resize(trajectory.size(), 0.);
			for (size_t pt_idx = 0; pt_idx < trajectory.size(); ++pt_idx) {
				_x_v[pt_idx] = trajectory.X(pt_idx);
				_y_v[pt_idx] = trajectory.Y(pt_idx);
				_z_v[pt_idx] = trajectory.Z(pt_idx);
				_time_v[pt_idx] = trajectory.T(pt_idx);
				_energy_v[pt_idx] = trajectory.E(pt_idx);
			}
			for (auto const & pair : trajectory) {
				TLorentzVector const& pos = pair.first;
				_x_v.push_back(pos.X());
				_y_v.push_back(pos.Y());
				_z_v.push_back(pos.Z());
				_time_v.push_back(pos.T());
				_energy_v.push_back(pos.E());
				std::cout << pos << std::endl;
			}
			_particletree->Fill();
		}
	}*/

	art::Handle<std::vector< simb::MCParticle > > mcparticle_h;
	e.getByLabel(_trajectory_label, mcparticle_h);
	for (size_t idx = 0; idx < mcparticle_h->size(); ++idx) {
		auto const& particle = (*mcparticle_h)[idx];
		//for (auto const& particle : mcparticle_v) {
			_pdg_code = particle.PdgCode();
			_x = particle.Vx();
			_y = particle.Vy();
			_z = particle.Vz();
			_time = particle.T();
			_energy = particle.E();
			_momentum = particle.P();
			_track_id = particle.TrackId();

			const auto& trajectory = particle.Trajectory();
			_x_v.resize(trajectory.size(), 0.);
			_y_v.resize(trajectory.size(), 0.);
			_z_v.resize(trajectory.size(), 0.);
			_time_v.resize(trajectory.size(), 0.);
			_energy_v.resize(trajectory.size(), 0.);
			for (size_t pt_idx = 0; pt_idx < trajectory.size(); ++pt_idx) {
				_x_v[pt_idx] = trajectory.X(pt_idx);
				_y_v[pt_idx] = trajectory.Y(pt_idx);
				_z_v[pt_idx] = trajectory.Z(pt_idx);
				_time_v[pt_idx] = trajectory.T(pt_idx);
				_energy_v[pt_idx] = trajectory.E(pt_idx);
			}
			_particletree->Fill();
		//}
	}

	// Also record sim::SimEnergyDeposit from largeant?
	//art::Handle<std::vector< sim::SimEnergyDeposit > > simenergy_h;
	//e.getByLabel(_simenergy_label, simenergy_h);
	// Plus Kazu's algorithm to sparsify?
}

DEFINE_ART_MODULE(ICARUSParticleAna)
