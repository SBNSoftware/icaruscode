////////////////////////////////////////////////////////////////////////
// Class:       CopyDaqICARUSTPC
// Plugin Type: producer (Unknown Unknown)
// File:        CopyDaqICARUSTPC_module.cc
//
// Generated at Tue Nov  7 13:07:45 2023 by Giuseppe Cerati using cetskelgen
// from  version .
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

#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
// #include "sbnobj/Common/Trigger/BeamBits.h"
// #include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration
// #include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"
// #include "sbndaq-artdaq-core/Overlays/FragmentType.hh" // sbndaq::FragmentType
#include "artdaq-core/Data/Fragment.hh"

#include <memory>

namespace {

  /// Moves the content of `data` into a `std::unique_ptr`.
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& data)
    { return std::make_unique<T>(std::move(data)); }

} // local namespace

class CopyDaqICARUSTPC;

class CopyDaqICARUSTPC : public art::EDProducer {
public:
  explicit CopyDaqICARUSTPC(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CopyDaqICARUSTPC(CopyDaqICARUSTPC const&) = delete;
  CopyDaqICARUSTPC(CopyDaqICARUSTPC&&) = delete;
  CopyDaqICARUSTPC& operator=(CopyDaqICARUSTPC const&) = delete;
  CopyDaqICARUSTPC& operator=(CopyDaqICARUSTPC&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::vector<art::InputTag> const fInputTags;
};


CopyDaqICARUSTPC::CopyDaqICARUSTPC(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fInputTags(p.get<std::vector<art::InputTag>>("FragmentsLabelVec"))
{
  for (art::InputTag const& inputTag: fInputTags) {
    consumes<artdaq::Fragments>(inputTag);
    produces<artdaq::Fragments>(inputTag.instance());
  }
}

void CopyDaqICARUSTPC::produce(art::Event& e)
{
  for (art::InputTag const& inputTag: fInputTags) {
     auto const& thisHandle = e.getHandle<artdaq::Fragments>(inputTag);
     artdaq::Fragments frags = *thisHandle;
     e.put(moveToUniquePtr(frags), inputTag.instance());
  }  
}

DEFINE_ART_MODULE(CopyDaqICARUSTPC)
