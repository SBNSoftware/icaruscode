////////////////////////////////////////////////////////////////////////
// Class:       ICARUSMCOpHit
// Plugin Type: producer (art v3_01_02)
// File:        ICARUSMCOpHit_module.cc
//
// Generated at Sun Mar  3 16:52:53 2019 by Kazuhiro Terao using cetskelgen
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

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpHit.h"

#include <memory>

class ICARUSMCOpHit;


class ICARUSMCOpHit : public art::EDProducer {
public:
  explicit ICARUSMCOpHit(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSMCOpHit(ICARUSMCOpHit const&) = delete;
  ICARUSMCOpHit(ICARUSMCOpHit&&) = delete;
  ICARUSMCOpHit& operator=(ICARUSMCOpHit const&) = delete;
  ICARUSMCOpHit& operator=(ICARUSMCOpHit&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  double _merge_period;
  std::string _simph_producer;
};


ICARUSMCOpHit::ICARUSMCOpHit(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  _merge_period = p.get<double>("MergePeriod");
  _simph_producer = p.get<std::string>("SimPhotonsProducer");
  produces<std::vector<recob::OpHit> >();
}

void ICARUSMCOpHit::produce(art::Event& e)
{
  auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
  auto oph_v = std::unique_ptr<std::vector<recob::OpHit> >(new std::vector<recob::OpHit>());

  art::Handle< std::vector< sim::SimPhotons > > simph_h;
  e.getByLabel(_simph_producer,simph_h);
  if(!simph_h.isValid()) {
    std::cerr << "Could not retrieve sim::SimPhotons from producer label: " << _simph_producer << std::endl;
    throw std::exception();
  }

  std::vector<bool> processed_v;
  for(auto const& simph : *simph_h) {
    // Make sure channel number is unique (e.g. one sim::SimPhotons per op channel)
    size_t opch = simph.OpChannel();

    if(opch >= processed_v.size()) processed_v.resize(opch+1,false);
    if(processed_v[opch]) {
      std::cerr << "Found duplicate channels in std::vector<sim::SimPhotons>! not expected (logic will fail).."<<std::endl;
      throw std::exception();
    }

    processed_v[opch] = true;
    bool in_window  = false;
    double oph_time = -1.e9;
    double pe = 0.;
    // Insert photon times into a sorted set
    std::map<double,size_t> time_m;
    for(auto const& oneph : simph) {
      double this_time = ts->G4ToElecTime(oneph.Time) - ts->TriggerTime();
      time_m[this_time] += 1;
    }

    // Loop over the time vector, emplace photons
    for(auto const& time_photon_pair : time_m) {

      auto const& this_time = time_photon_pair.first;
      if(this_time > (oph_time + _merge_period) && in_window) {
	recob::OpHit oph(opch, 
			 oph_time,
			 oph_time + ts->TriggerTime(),
			 0, // frame
			 1., // width
			 pe, // area,
			 pe, // peakheight,
			 pe, // pe
			 0.);
	oph_v->emplace_back(std::move(oph));
	in_window = false;
	pe = 0;
      }

      if(!in_window) oph_time = this_time;
      in_window = true;
      pe += time_photon_pair.second;
    }
    if(in_window) {
      recob::OpHit oph(opch, 
		       oph_time,
		       oph_time + ts->TriggerTime(),
		       0, // frame
		       1., // width
		       pe, // area,
		       pe, // peakheight,
		       pe, // pe
		       0.);
      oph_v->emplace_back(std::move(oph));
    }
  }

  e.put(std::move(oph_v));
}

DEFINE_ART_MODULE(ICARUSMCOpHit)
