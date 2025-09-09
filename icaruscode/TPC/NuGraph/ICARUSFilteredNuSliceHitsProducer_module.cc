#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <algorithm> 
#include <cmath> 
#include <memory>
#include <unordered_map>
#include <utility> 
#include <vector>

#include "lardataobj/RecoBase/OpFlash.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"

#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

class ICARUSFilteredNuSliceHitsProducer : public art::EDProducer {

public:
  explicit ICARUSFilteredNuSliceHitsProducer(fhicl::ParameterSet const& p);

  ICARUSFilteredNuSliceHitsProducer(ICARUSFilteredNuSliceHitsProducer const&) = delete;
  ICARUSFilteredNuSliceHitsProducer(ICARUSFilteredNuSliceHitsProducer&&) = delete;
  ICARUSFilteredNuSliceHitsProducer& operator=(ICARUSFilteredNuSliceHitsProducer const&) = delete;
  ICARUSFilteredNuSliceHitsProducer& operator=(ICARUSFilteredNuSliceHitsProducer&&) = delete;

private:

  art::InputTag const fHitLabel;
  art::InputTag const fHitScoreLabel;
  float fScoreCut;
  bool fWeightFilterByPFP;

  void produce(art::Event& e) override;

};

ICARUSFilteredNuSliceHitsProducer::ICARUSFilteredNuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitLabel(p.get<art::InputTag>("HitLabel", "nuslhitsCryoE")),
  fHitScoreLabel(p.get<art::InputTag>("HitScoreLabel", "NuGraphCryoE:filter")),
  fScoreCut(p.get<float>("ScoreCut", 0))
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ICARUSFilteredNuSliceHitsProducer::produce(art::Event& e)
{

  auto outputHits = std::make_unique<std::vector<recob::Hit>>();

  // get hits from the tagged slice
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);

  // get hit filter values
  art::Handle<std::vector<anab::FeatureVector<1>>> filterOutputHandle;
  e.getByLabel(fHitScoreLabel, filterOutputHandle);

  for (size_t ihit = 0; ihit < hitListHandle->size(); ihit++) {
    art::Ptr<recob::Hit> hit(hitListHandle, ihit);
    if (fScoreCut >= 0 && filterOutputHandle->at(ihit).at(0) >= fScoreCut) 
      outputHits->emplace_back(*hit);
  }

  std::cout << "Number of hits after ICARUSFilteredNuSliceHitsProducer: " << outputHits->size() << std::endl;
  e.put(std::move(outputHits));

}

DEFINE_ART_MODULE(ICARUSFilteredNuSliceHitsProducer)