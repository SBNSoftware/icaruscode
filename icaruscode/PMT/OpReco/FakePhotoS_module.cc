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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

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

  // Declare member data here.
  double _frequency; // [MHz]
  double _duration;  // [us]
  double _tstart;    // [ns]
  std::map<unsigned int,unsigned short> _pemap; // opch <=> #pe to be injected
};


FakePhotoS::FakePhotoS(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  auto const ch_v = p.get<std::vector<unsigned int> >("Channels");
  auto const pe_v = p.get<std::vector<unsigned int> >("PEs");
  assert(ch_v.size() == pe_v.size());
  for(size_t i=0; i<ch_v.size(); ++i) _pemap[ch_v[i]] = pe_v[i];
  _frequency = p.get<double>("Frequency");
  _duration  = p.get<double>("Duration");
  _tstart    = p.get<double>("G4TStart");
  produces<std::vector<sim::SimPhotons> >();
  produces< sumdata::RunData, art::InRun >();

}

void FakePhotoS::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using                                                                                                  
  art::ServiceHandle<geo::Geometry> geo;

  std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runData));

  return;
}

void FakePhotoS::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto simph_v = std::unique_ptr<std::vector<sim::SimPhotons> >(new std::vector<sim::SimPhotons>());
  simph_v->reserve(_pemap.size());
  for(auto const& key_val : _pemap) {
    
    auto const& ch  = key_val.first;
    auto const& npe = key_val.second;
    
    sim::SimPhotons sphot(ch);
    sphot.reserve(int(_duration * _frequency));
    double clock = 0.;
    while(clock <= _duration) {
      for(size_t ctr=0;ctr<npe;++ctr) {
	sim::OnePhoton phot;
	phot.Time = _tstart + clock * 1.e3;
	sphot.emplace_back(std::move(phot));
      }
      clock += (1./_frequency);
    }

    simph_v->emplace_back(std::move(sphot));
  }

  e.put(std::move(simph_v));
}

DEFINE_ART_MODULE(FakePhotoS)
