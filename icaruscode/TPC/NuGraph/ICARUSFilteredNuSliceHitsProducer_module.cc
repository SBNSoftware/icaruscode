/**
 * @file    icaruscode/TPC/NuGraph/ICARUSFilteredNuSliceHitsProducer_module.cc
 * @brief   Producer to filter the input hits according to a configurable NuGraph2 filter score
 * @author  Riccardo Triozzi ( triozzi@pd.infn.it )
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

/**
 * @brief Filters the input hits according to a configurable NuGraph2 filter score.
 * 
 * This module gets the slices from a `SliceLabel` (e.g., `NCCSlicesCryoE`, 
 * corresponding to all the slices tagged as not-clear-cosmic by Pandora in 
 * the east cryostat) and the associated hit via the corresponding `PandoraLabel`
 * (e.g., `pandoraGausCryoE`).
 * For each hit, it grabs the NuGraph2 filter and semantic predictions
 * from the corresponding `NuGraphLabel` (e.g., `NGMultiSliceCryoE`). 
 * In the current infrastructure, all the hits are guaranteed to
 * have a corresponding NuGraph2 prediction.
 * 
 * The hits are filtered when their NuGraph2 filter score is above a configurable
 * `ScoreCut`, by default set to zero. A zero `ScoreCut` means that every hit is
 * propagated to the output, and no filtering happens.
 * 
 * Output
 * -------
 * 
 * Collections of filtered hits, along with their filter and semantic predictions.
 * Associations between the filtered hits and their filter and semantic predictions.
 * Both are needed to fit into the current Pandora and CAF infrastructures.
 * 
 * Configuration parameters
 * -------------------------
 * * `SliceLabel` (string, mandatory): tag of input slices (with cryostat tag), and
 *     their associations to hits, space points, particle flow objects etc.
 * * `PandoraLabel` (string, mandatory): tag of Pandora objects (with cryostat tag).
 * * `NuGraphLabel` (string, mandatory): tag of NuGraph2 predictions 
 *     (with cryostat tag).
 * * `ScoreCut` (real, default: 0): hits are filtered when their
 *     NuGraph2 filter score is <= `ScoreCut`.
 * 
 */

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
  produces<std::vector<recob::Hit>>();
  produces<std::vector<anab::FeatureVector<1>>>();
  produces<std::vector<anab::FeatureVector<5>>>();
  produces<art::Assns<recob::Hit, anab::FeatureVector<1>>>("filter");
  produces<art::Assns<recob::Hit, anab::FeatureVector<5>>>("semantic");
}

void ICARUSFilteredNuSliceHitsProducer::produce(art::Event& e)
{
  auto outputHits = std::make_unique<std::vector<recob::Hit>>();
  auto outputFilter = std::make_unique<std::vector<anab::FeatureVector<1>>>();
  auto outputSemantic = std::make_unique<std::vector<anab::FeatureVector<5>>>();
  auto outputHitFilterAssns = std::make_unique<art::Assns<recob::Hit, anab::FeatureVector<1>>>();
  auto outputHitSemanticAssns = std::make_unique<art::Assns<recob::Hit, anab::FeatureVector<5>>>();

  // Get slices.
  const std::vector<art::Ptr<recob::Slice>> slices = e.getProduct<std::vector<art::Ptr<recob::Slice>>>(fSliceLabel);
  art::FindManyP<recob::Hit> sliceToHitsAssoc(slices, e, fPandoraLabel);

  // Get input hits.
  std::vector<art::Ptr<recob::Hit>> inputHits;
  for (size_t islc = 0; islc < slices.size(); ++islc) {
    const std::vector<art::Ptr<recob::Hit>> hitsInSlice = sliceToHitsAssoc.at(islc);
    for (auto const& h : hitsInSlice) 
      inputHits.emplace_back(h);
  }
  std::cout << "Number of hits before ICARUSFilteredNuSliceHitsProducer: " << inputHits.size() << std::endl;

  // Get NuGraph filter and semantic predictions.
  art::FindOneP<anab::FeatureVector<1>> hitToNGFilterAssoc(inputHits, e, art::InputTag(fNGLabel.label(), "filter"));
  art::FindOneP<anab::FeatureVector<5>> hitToNGSemanticAssoc(inputHits, e, art::InputTag(fNGLabel.label(), "semantic"));
  art::PtrMaker<recob::Hit> hitPtrMaker{e};

  // Filter input hits, re-create filter and sematic objects (needed for Pandora) and associations (needed for CAFs).
  // By setting `fScoreCut` to zero, not filtering happens and all the hits are propagated.
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