////////////////////////////////////////////////////////////////////////
// Class:       FilterOpFlash
// Plugin Type: filter (art v3_06_03)
// File:        FilterOpFlash_module.cc
//
// Generated at Wed Mar 24 13:20:35 2021 by Kazuhiro Terao using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include <memory>

class FilterOpFlash;


class FilterOpFlash : public art::EDFilter {
public:
  explicit FilterOpFlash(fhicl::ParameterSet const& p);

  FilterOpFlash(FilterOpFlash const&) = delete;
  FilterOpFlash(FilterOpFlash&&) = delete;
  FilterOpFlash& operator=(FilterOpFlash const&) = delete;
  FilterOpFlash& operator=(FilterOpFlash&&) = delete;
  bool filter(art::Event& e) override;

private:
  std::vector<std::string> _flash_producer_v;
  double _time_start;
  double _time_end;
  double _pe_threshold;
};


FilterOpFlash::FilterOpFlash(fhicl::ParameterSet const& p)
  : EDFilter{p}
{
  _flash_producer_v = p.get<std::vector<std::string> >("OpFlashProducerList");
  _time_start = p.get<double>("WindowStartTime");
  _time_end   = p.get<double>("WindowEndTime");
  _pe_threshold = p.get<double>("FlashPEThreshold",-1);
  assert(_time_end>_time_start);
}

bool FilterOpFlash::filter(art::Event& e)
{

  bool pass = false;
  for(auto const& producer : _flash_producer_v) {
    
    art::Handle<std::vector<recob::OpFlash> > flash_handle;
    e.getByLabel(producer, flash_handle);
    if(!flash_handle.isValid()) {
	std::cerr << "Invalid producer for truth recob::OpFlash: " << producer<< std::endl;
	throw std::exception();
    }

    for(auto const& flash : *flash_handle) {
      if(flash.TotalPE() < _pe_threshold) continue;
      if(flash.Time() < _time_start || flash.Time() > _time_end) continue;
      pass = true;
      break;
    }

    if(pass) break;
  }
  return pass;
}

DEFINE_ART_MODULE(FilterOpFlash)
