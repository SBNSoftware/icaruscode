////////////////////////////////////////////////////////////////////////
// Class:       SubRunPOTInEvent
// Plugin Type: producer (art v3_00_00)
// File:        SubRunPOTInEvent_module.cc
//
// Generated at Mon Jan 21 09:09:38 2019 by Wesley Ketchum using cetskelgen
// from cetlib version v3_04_00.
//
// This module directly copies the subrun POT info into every event.
// **IT IS NOT A POT PER EVENT IT IS NOT A POT PER EVENT**
// **IT IS NOT A POT PER EVENT IT IS NOT A POT PER EVENT**
// **IT IS NOT A POT PER EVENT IT IS NOT A POT PER EVENT**
// **IT IS NOT A POT PER EVENT IT IS NOT A POT PER EVENT**
// **IT IS NOT A POT PER EVENT IT IS NOT A POT PER EVENT**
// **IT IS NOT A POT PER EVENT IT IS NOT A POT PER EVENT**
//
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

#include <memory>

#include "larcoreobj/SummaryData/POTSummary.h"

namespace mix {
  class SubRunPOTInEvent;
}


class mix::SubRunPOTInEvent : public art::EDProducer {
public:
  explicit SubRunPOTInEvent(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SubRunPOTInEvent(SubRunPOTInEvent const&) = delete;
  SubRunPOTInEvent(SubRunPOTInEvent&&) = delete;
  SubRunPOTInEvent& operator=(SubRunPOTInEvent const&) = delete;
  SubRunPOTInEvent& operator=(SubRunPOTInEvent&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginSubRun(art::SubRun& sr) override;

private:

  art::InputTag fInputTag;
  sumdata::POTSummary fSubRunPOT;

};


mix::SubRunPOTInEvent::SubRunPOTInEvent(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fInputTag = p.get<art::InputTag>("InputTag");

  produces<sumdata::POTSummary>("SubRunPOT");
}

void mix::SubRunPOTInEvent::produce(art::Event& e)
{
  std::unique_ptr<sumdata::POTSummary> srpot_ptr(new sumdata::POTSummary(fSubRunPOT));
  e.put(std::move(srpot_ptr),"SubRunPOT");
}

void mix::SubRunPOTInEvent::beginSubRun(art::SubRun& sr)
{
  art::Handle<sumdata::POTSummary> srpot_handle;
  sr.getByLabel(fInputTag,srpot_handle);
  if(srpot_handle.isValid()) fSubRunPOT = *srpot_handle;
  else std::cout << "BAD SUBRUNPOT\n";

  std::cout << "GOT POT: " << fSubRunPOT.totgoodpot << std::endl;
}

DEFINE_ART_MODULE(mix::SubRunPOTInEvent)

