/**
 * @file    icaruscode/TPC/NuGraph/ICARUSFilteredNuSliceHitsProducer_module.cc
 * @brief   Implementation of `ICARUSFilteredNuSliceHitsProducer` _art_ module.
 * @author  Riccardo Triozzi
 * @date    September 9, 2025
 */

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

  // input tags
  art::InputTag const fSliceLabel;
  art::InputTag const fPandoraLabel;
  art::InputTag const fNGLabel;

  // module parameters
  float fScoreCut;

  void produce(art::Event& e) override;

};

ICARUSFilteredNuSliceHitsProducer::ICARUSFilteredNuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fSliceLabel(p.get<art::InputTag>("SliceLabel", "NCCSlicesCryoE")),
  fPandoraLabel{p.get<std::string>("PandoraLabel", "pandoraGausCryoE")},
  fNGLabel{p.get<std::string>("NuGraphLabel", "NGMultiSliceCryoE")},
  fScoreCut(p.get<float>("ScoreCut", 0))
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
  produces<std::vector<anab::FeatureVector<1>>>();
  produces<std::vector<anab::FeatureVector<5>>>();
  produces<art::Assns<recob::Hit, anab::FeatureVector<1>>>("filter");
  produces<art::Assns<recob::Hit, anab::FeatureVector<5>>>("semantic");

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ICARUSFilteredNuSliceHitsProducer::produce(art::Event& e)
{
  auto outputHits = std::make_unique<std::vector<recob::Hit>>();
  auto outputFilter = std::make_unique<std::vector<anab::FeatureVector<1>>>();
  auto outputSemantic = std::make_unique<std::vector<anab::FeatureVector<5>>>();
  auto outputHitFilterAssns = std::make_unique<art::Assns<recob::Hit, anab::FeatureVector<1>>>();
  auto outputHitSemanticAssns = std::make_unique<art::Assns<recob::Hit, anab::FeatureVector<5>>>();

  // get slices
  const std::vector<art::Ptr<recob::Slice>> slices = e.getProduct<std::vector<art::Ptr<recob::Slice>>>(fSliceLabel);
  art::FindManyP<recob::Hit> sliceToHitsAssoc(slices, e, fPandoraLabel);

  // get input hits
  std::vector<art::Ptr<recob::Hit>> inputHits;
  for (size_t islc = 0; islc < slices.size(); ++islc) {
    const std::vector<art::Ptr<recob::Hit>> hitsInSlice = sliceToHitsAssoc.at(islc);
    for (auto const& h : hitsInSlice) 
      inputHits.emplace_back(h);
  }
  std::cout << "Number of hits before ICARUSFilteredNuSliceHitsProducer: " << inputHits.size() << std::endl;

  // get NuGraph filter and semantic predictions
  art::FindOneP<anab::FeatureVector<1>> hitToNGFilterAssoc(inputHits, e, art::InputTag(fNGLabel.label(), "filter"));
  art::FindOneP<anab::FeatureVector<5>> hitToNGSemanticAssoc(inputHits, e, art::InputTag(fNGLabel.label(), "semantic"));
  art::PtrMaker<recob::Hit> hitPtrMaker{e};

  // filter input hits, re-create filter and sematic objects (needed for Pandora) and associations (needed for CAFs)
  for (size_t ihit = 0; ihit < inputHits.size(); ihit++) {
    art::Ptr<recob::Hit> hit = inputHits[ihit];
    if (hitToNGFilterAssoc.at(ihit).isNonnull()) {
      if (fScoreCut >= 0 && hitToNGFilterAssoc.at(ihit)->at(0) >= fScoreCut) {
        outputHits->emplace_back(*hit);
        outputFilter->emplace_back(*hitToNGFilterAssoc.at(ihit));
        outputSemantic->emplace_back(*hitToNGSemanticAssoc.at(ihit));
        outputHitFilterAssns->addSingle(hitPtrMaker(outputHits->size() -1), hitToNGFilterAssoc.at(ihit));
        outputHitSemanticAssns->addSingle(hitPtrMaker(outputHits->size() -1), hitToNGSemanticAssoc.at(ihit));
      }
    }
  }
  std::cout << "Number of hits after ICARUSFilteredNuSliceHitsProducer: " << outputHits->size() << std::endl;
  
  e.put(std::move(outputHits));
  e.put(std::move(outputFilter));
  e.put(std::move(outputSemantic));
  e.put(std::move(outputHitFilterAssns), "filter");
  e.put(std::move(outputHitSemanticAssns), "semantic");
}

DEFINE_ART_MODULE(ICARUSFilteredNuSliceHitsProducer)