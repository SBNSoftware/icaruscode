////////////////////////////////////////////////////////////////////////
// Class:       FakeFlash
// Plugin Type: producer (art v3_01_01)
// File:        FakeFlash_module.cc
//
// Generated at Tue Feb 12 13:46:25 2019 by Kazuhiro Terao using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandExponential.h"
#include "TRandom.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <memory>
#include <TLorentzVector.h>

class FakeFlash;


class FakeFlash : public art::EDProducer {
public:
  explicit FakeFlash(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FakeFlash(FakeFlash const&) = delete;
  FakeFlash(FakeFlash&&) = delete;
  FakeFlash& operator=(FakeFlash const&) = delete;
  FakeFlash& operator=(FakeFlash&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginRun(art::Run& run) override;
private:

  void GenPosition(double& x, double& y, double& z);
  void FillSimPhotons(std::vector<sim::SimPhotons>& simph_v,
		      int nphotons,
		      size_t mother_trackid,
		      const TLorentzVector& pos);
  std::vector<double> GenerateTime(size_t numphotons);
  // Declare member data here.
  bool _verbose; ///< verbosity for debugging
  double _frequency; ///< [MHz]
  double _duration;  ///< [us]
  double _tstart;    ///< [ns]
  size_t _min_photons;  ///< [photons]
  size_t _max_photons;  ///< [photons]
  std::vector<size_t> _tpc_v; ///< List of TPC ID to be used
  double _fast_frac; ///< fraction of prompt light
  double _fast_tau;  ///< scintillation emission time constant for fast component
  double _slow_tau;  ///< scintillation emission time constant for slow component
  size_t _ch_min;    ///< channel range min to produce SimPhotons
  size_t _ch_max;    ///< channel range max to produce SimPhotons
  double _xmax;      ///< 0.0-1.0 the x-position range in fraction of a TPC volume
  double _ymax;      ///< 0.0-1.0 the y-position range in fraction of a TPC volume
  double _zmax;      ///< 0.0-1.0 the z-position range in fraction of a TPC volume
  double _xmin;      ///< 0.0-1.0 the x-position range in fraction of a TPC volume
  double _ymin;      ///< 0.0-1.0 the y-position range in fraction of a TPC volume
  double _zmin;      ///< 0.0-1.0 the z-position range in fraction of a TPC volume
  CLHEP::HepRandomEngine& fFlatEngine;
  CLHEP::RandFlat *fFlatRandom;
  CLHEP::RandExponential* fExpoRandom;
  CLHEP::RandPoisson* fPoisRandom;
};


FakeFlash::FakeFlash(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fFlatEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "Gen", p, "Seed"))
  // More initializers here.
{
  _verbose = p.get<bool>("Verbose",false); // If you want someone to talk to you
  auto min_photons = p.get<int>("MinPhotons",24000);     // Min of the range of photons to be injected in one shot
  auto max_photons = p.get<int>("MaxPhotons",2400000);   // Max of the range of photons to be injected in one shot
  assert(min_photons < max_photons && min_photons>0 && max_photons>0); 
  _min_photons = min_photons;
  _max_photons = max_photons;

  // scintillation params (yeah fine you can replace these w/ service)
  _fast_frac = p.get<double>("PromptLightFraction",0.23);
  _fast_tau  = p.get<double>("FastTimeConstant",0.006);
  _slow_tau  = p.get<double>("SlowTimeConstant",1.5);

  // channel range to simulate
  auto geop = lar::providerFrom<geo::Geometry>();
  _ch_min = 0;
  _ch_max = geop->NOpChannels() - 1;
  _ch_min = p.get<size_t>("ChannelMin",_ch_min);
  _ch_max = p.get<size_t>("ChannelMax",_ch_max);
  assert(_ch_min<_ch_max);
  // Given in micro-seconds as larsoft default time unit, but then we convert to ns to record as photon time
  _fast_tau *= 1.e3;
  _slow_tau *= 1.e3;
  // TPC list + range
  _tpc_v = p.get<std::vector<size_t> >("TPCList"); // The list of TPC in which flash can happen
  _xmax  = p.get<double>("XMax",1.0);
  _xmin  = p.get<double>("XMin",0.0);
  _ymax  = p.get<double>("YMax",1.0);
  _ymin  = p.get<double>("YMin",0.0);
  _zmax  = p.get<double>("ZMax",1.0);
  _zmin  = p.get<double>("ZMin",0.0);
  assert(_xmax>=_xmin && _ymax>=_ymin && _zmax>=_zmin &&
	 _xmax<=1.0 && _ymax<=1.0 && _zmax<=1.0 &&
	 _xmin>=0.0 && _ymin>=0.0 && _zmin>=0.0 );
  // generation timing info
  _frequency = p.get<double>("Frequency"); // The frequency of photon(s) injection
  _duration  = p.get<double>("Duration");  // The duration of photon(s) injection period
  _tstart    = p.get<double>("G4TStart");  // The start time of photon injection period
  produces<std::vector<sim::SimPhotons> >();
  produces<std::vector<simb::MCTruth> >();
  produces< sumdata::RunData, art::InRun >();

  fFlatRandom = new CLHEP::RandFlat(fFlatEngine,0,1);
  fExpoRandom = new CLHEP::RandExponential(fFlatEngine);
  fPoisRandom = new CLHEP::RandPoisson(fFlatEngine);
}

void FakeFlash::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using                                                                                                  
  art::ServiceHandle<geo::Geometry> geo;

  std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runData));

  return;
}

void FakeFlash::GenPosition(double& x, double& y, double& z) {
    
    size_t tpc_id = (size_t)(fFlatRandom->fire(0,_tpc_v.size()));
    bool found = false;
    // Implementation of required member function here.
    auto geop = lar::providerFrom<geo::Geometry>();
    for(size_t c=0; c<geop->Ncryostats(); ++c) {
        auto const& cryostat = geop->Cryostat(c);
        if(!cryostat.HasTPC(tpc_id)) continue;
        auto const& tpc = cryostat.TPC(tpc_id);
        auto const& tpcabox = tpc.ActiveBoundingBox();
        double xmin = tpcabox.MinX() + (tpcabox.MaxX() - tpcabox.MinX()) * _xmin;
        double xmax = tpcabox.MinX() + (tpcabox.MaxX() - tpcabox.MinX()) * _xmax;
        double ymin = tpcabox.MinY() + (tpcabox.MaxY() - tpcabox.MinY()) * _ymin;
        double ymax = tpcabox.MinY() + (tpcabox.MaxY() - tpcabox.MinY()) * _ymax;
        double zmin = tpcabox.MinZ() + (tpcabox.MaxZ() - tpcabox.MinZ()) * _zmin;
        double zmax = tpcabox.MinZ() + (tpcabox.MaxZ() - tpcabox.MinZ()) * _zmax;
        x = fFlatRandom->fire(xmin,xmax);
        y = fFlatRandom->fire(ymin,ymax);
        z = fFlatRandom->fire(zmin,zmax);
        found = true;
        break;
    }
    if(!found) std::cerr<< "\033[93mTPC " << tpc_id << " not found...\033[00m" << std::endl;
}

void FakeFlash::FillSimPhotons(std::vector<sim::SimPhotons>& simph_v, 
			       int mother_trackid,
			       size_t nphotons, 
			       const TLorentzVector& pos)
{
  art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
  double xyz[3];
  xyz[0] = pos.X();
  xyz[1] = pos.Y();
  xyz[2] = pos.Z();
  auto const& Visibilities = pvs->GetAllVisibilities(xyz);
  // Now loop over optical channels and fill SimPhotons array
  if(simph_v.empty()) {
    for(size_t opch=_ch_min; opch<=_ch_max; ++opch) 
      simph_v.push_back(sim::SimPhotons(opch));
  }
  for(size_t opch=_ch_min; opch <= _ch_max; ++opch) {
    // get visibility
    size_t detected = ((double)(nphotons)) * Visibilities[opch];
    std::cout<<opch<<","<<detected<<std::endl;
    auto time_array = this->GenerateTime(detected);
    assert(time_array.size() == detected);
    // record
    auto& simph = simph_v[opch];
    for(size_t idx=0; idx<detected; ++idx) {
      	sim::OnePhoton phot;
	phot.Time = pos.T() + time_array[idx];
	phot.InitialPosition.SetX(pos.X());
	phot.InitialPosition.SetY(pos.Y());
	phot.InitialPosition.SetZ(pos.Z());
	phot.MotherTrackID = mother_trackid;
	simph.emplace_back(phot);
    }
  }
}

std::vector<double> FakeFlash::GenerateTime(size_t numphotons) {

  double fast_expected = _fast_frac * numphotons;
  // draw poison
  int fast_count = std::max((long int)0,std::min(fPoisRandom->fire(fast_expected),(long int)numphotons));
  std::vector<double> res(numphotons);
  for(int i=0; i<((int)numphotons); ++i) {
    if(i<fast_count)
      res[i] = fExpoRandom->fire(_fast_tau);
    else res[i] = fExpoRandom->fire(_slow_tau);
  }
  return res;
}

void FakeFlash::produce(art::Event& e)
{

  auto simph_v = std::unique_ptr<std::vector<sim::SimPhotons> >(new std::vector<sim::SimPhotons>());
  auto mct_v   = std::unique_ptr<std::vector<simb::MCTruth> >(new std::vector<simb::MCTruth>());
  std::vector<simb::MCParticle> mcp_v;
  double clock = 0.;
  while(clock <= _duration) {
    // Determine photon count
    int nphotons = _min_photons + (_max_photons - _min_photons) * (fFlatRandom->fire(0,1));
    
    // Generate position 4 vector
    double time = _tstart + clock * 1.e3;
    double x,y,z;
    this->GenPosition(x,y,z);
    TLorentzVector pos(x,y,z,time);

    // Fill SimPhotons
    this->FillSimPhotons(*simph_v, (int)(mcp_v.size()), nphotons, pos);

    // Record a fake photon particle
    simb::MCParticle part(mcp_v.size(), 22, "primary", 0, 0., 1);
    TLorentzVector mom(0.,0.,0.,(double)nphotons);
    part.AddTrajectoryPoint(pos,mom);
    mcp_v.emplace_back(std::move(part));

    // count time
    clock += (1./_frequency);
  }

  simb::MCTruth mct;
  for(auto& part : mcp_v)
    mct.Add(part);
  mct_v->emplace_back(std::move(mct));

  e.put(std::move(simph_v));
  e.put(std::move(mct_v));
}

DEFINE_ART_MODULE(FakeFlash)
