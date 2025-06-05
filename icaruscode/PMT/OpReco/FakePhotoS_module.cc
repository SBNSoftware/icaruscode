////////////////////////////////////////////////////////////////////////
// Class:       FakePhotoS
// Plugin Type: producer (art v3_01_01)
// File:        FakePhotoS_module.cc
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
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandExponential.h"
#include "TRandom.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include <memory>

class FakePhotoS;


class FakePhotoS : public art::EDProducer {
public:
  explicit FakePhotoS(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FakePhotoS(FakePhotoS const&) = delete;
  FakePhotoS(FakePhotoS&&) = delete;
  FakePhotoS& operator=(FakePhotoS const&) = delete;
  FakePhotoS& operator=(FakePhotoS&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginRun(art::Run& run) override;

private:
  std::vector<double> GenerateTime(size_t numphotons);
  // Declare member data here.
  bool _verbose;
  double _frequency; ///< [MHz]
  double _duration;  ///< [us]
  double _tstart;    ///< [ns]
  size_t _min_pe;    ///< [p.e.]
  size_t _max_pe;    ///< [p.e.]
  bool   _simulate_time; ///< whether or not to simulate scintillation fall time
  double _fast_frac; ///< fraction of prompt light
  double _fast_tau;  ///< scintillation emission time constant for fast component
  double _slow_tau;  ///< scintillation emission time constant for slow component
  std::vector<unsigned int> _ch_v; // opchannels to create photons for

  CLHEP::HepRandomEngine &fFlatEngine;
  CLHEP::RandFlat        *fFlatRandom;
  CLHEP::RandExponential *fExpoRandom;
};


FakePhotoS::FakePhotoS(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fFlatEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(createEngine(0, "HepJamesRandom", "Gen"), "HepJamesRandom", "Gen", p, "Seed"))
  // More initializers here.
{
  _verbose = p.get<bool>("Verbose",false); // If you want someone to talk to you
  auto min_pe = p.get<int>("MinPE",1);     // Min of the range of PE count to be injected in one shot
  auto max_pe = p.get<int>("MaxPE",1);     // Max of the range of PE count to be injected in one shot
  assert(min_pe <= max_pe && min_pe>0 && max_pe>0); 
  _min_pe = min_pe;
  _max_pe = max_pe;
  _ch_v.clear();
  _ch_v = p.get<std::vector<unsigned int> >("Channels",_ch_v); // Specify if you want to use only select channels
  _frequency = p.get<double>("Frequency"); // The frequency of photon(s) injection
  _duration  = p.get<double>("Duration");  // The duration of photon(s) injection period
  _tstart    = p.get<double>("G4TStart");  // The start time of photon injection period
  produces<std::vector<sim::SimPhotons> >();
  produces< sumdata::RunData, art::InRun >();
  // scintillation params (yeah fine you can replace these w/ service)
  _simulate_time = p.get<bool>("SimulateTime",false); // scintillation fall timing simulation
  _fast_frac = p.get<double>("PromptLightFraction",0.23);
  _fast_tau  = p.get<double>("FastTimeConstant",0.006);
  _slow_tau  = p.get<double>("SlowTimeConstant",1.5);

  fFlatRandom = new CLHEP::RandFlat(fFlatEngine,_min_pe,_max_pe);
  fExpoRandom = new CLHEP::RandExponential(fFlatEngine);
}

void FakePhotoS::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using                                                                                                  
  art::ServiceHandle<geo::Geometry> geo;

  std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runData), art::fullRun());

  return;
}


std::vector<double> FakePhotoS::GenerateTime(size_t numphotons) {

  std::vector<double> res(numphotons);
  for(int i=0; i<((int)numphotons); ++i) {
    if(fFlatRandom->fire(0.,1.) < _fast_frac)
      res[i] = fExpoRandom->fire(_fast_tau);
    else res[i] = fExpoRandom->fire(_slow_tau);
  }
  return res;
}

void FakePhotoS::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto simph_v = std::unique_ptr<std::vector<sim::SimPhotons> >(new std::vector<sim::SimPhotons>());

  if(_ch_v.empty()) {
    auto const nOpChannels = art::ServiceHandle<geo::WireReadout const>()->Get().NOpChannels();
    _ch_v.reserve(nOpChannels);
    for(size_t opch=0; opch<nOpChannels; ++opch)
      _ch_v.push_back(opch);
  }
  simph_v->reserve(_ch_v.size());
  for(auto const& opch : _ch_v) {
    if(_verbose) std::cout << "OpChannel " << opch << std::endl;
    sim::SimPhotons sphot(opch);
    sphot.reserve(int(_duration * _frequency));
    double clock = 0.;
    while(clock <= _duration) {
      
      size_t npe = fFlatRandom->fire(_min_pe,_max_pe+0.999999);
      if(_verbose) std::cout << npe << "@" << (int)(_tstart + clock * 1.e3) << "[ns] " << std::flush; 
      std::vector<double> time_array(npe,0.);
      if(_simulate_time) time_array = this->GenerateTime(npe);
      for(size_t ctr=0;ctr<npe;++ctr) {
	sim::OnePhoton phot;
	phot.Time = _tstart + clock * 1.e3 + time_array[ctr];
	sphot.emplace_back(std::move(phot));
      }
      clock += (1./_frequency);
    }
    if(_verbose) std::cout << std::endl;
    simph_v->emplace_back(std::move(sphot));
  }

  e.put(std::move(simph_v));
}

DEFINE_ART_MODULE(FakePhotoS)
