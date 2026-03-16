/**
 * @file icaruscode/TPC/NuGraph/ICARUSNCCSlicesProducer_module.cc
 * @brief Implementation of `ICARUSNCCSlicesProducer` _art_ module.
 * @author Leonardo Lena (https://github.com/leonardo-lena)
 * @date September 9, 2025
 */

#include <memory>
#include <map>
#include <utility> // std::move()
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" // for DEFINE_ART_MODULE
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"

/**
 * @class ICARUSNCCSlicesProducer
 *
 * Produces an std::vector<art::Ptr<recob::Slice>> whose IsPrimary PFParticle is NOT tagged as IsClearCosmic by Pandora.
 * Also produces an Association between all recob::Slices and all recob::SpacePoints associated with atleast 3 recob::Hits in the same recob::Slice.
 */
class ICARUSNCCSlicesProducer : public art::EDProducer {

public:
  /**
   * @brief Constructor
   * @param params
   */
  explicit ICARUSNCCSlicesProducer(fhicl::ParameterSet const &params);

  // Plugins should not be copied or assigned.
  ICARUSNCCSlicesProducer(ICARUSNCCSlicesProducer const &) = delete;
  ICARUSNCCSlicesProducer(ICARUSNCCSlicesProducer &&) = delete;
  ICARUSNCCSlicesProducer &operator=(ICARUSNCCSlicesProducer const &) = delete;
  ICARUSNCCSlicesProducer &operator=(ICARUSNCCSlicesProducer &&) = delete;

private:
  art::InputTag const fPandoraLabel;
  art::InputTag const fClusterLabel;

  void produce(art::Event &event) override;
};

ICARUSNCCSlicesProducer::ICARUSNCCSlicesProducer(fhicl::ParameterSet const &params) : EDProducer{params},
  fPandoraLabel(params.get<art::InputTag>("PandoraLabel", "pandoraGausCryoE")),
  fClusterLabel(params.get<art::InputTag>("ClusterLabel", "cluster3DCryoE")) {

  produces<std::vector<art::Ptr<recob::Slice>>>();
  produces<art::Assns<recob::Slice, recob::SpacePoint>>();
}

void ICARUSNCCSlicesProducer::produce(art::Event &event) {
  auto outputSlices = std::make_unique<std::vector<art::Ptr<recob::Slice>>>();
  auto outputSlcSpsAssns = std::make_unique<art::Assns<recob::Slice, recob::SpacePoint>>();
  auto outputHits = std::make_unique<std::vector<art::Ptr<recob::Hit>>>();


  // begin slice selection
  art::ValidHandle<std::vector<recob::Slice>> slicesHandle = event.getValidHandle<std::vector<recob::Slice>>(fPandoraLabel);
  art::ValidHandle<std::vector<recob::SpacePoint>> spacePointsHandle = event.getValidHandle<std::vector<recob::SpacePoint>>(fClusterLabel);
  art::ValidHandle<std::vector<recob::PFParticle>> PFPHandle = event.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);
  art::ValidHandle<std::vector<recob::Hit>> hitsHandle = event.getValidHandle<std::vector<recob::Hit>>(fClusterLabel);

  art::FindManyP<recob::PFParticle> findMPfpFromSlice(slicesHandle, event, fPandoraLabel);
  art::FindManyP<recob::Hit> findMHitFromSP(spacePointsHandle, event, fClusterLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> find1MetaFromPfp(PFPHandle, event, fPandoraLabel);
  art::FindOneP<recob::Slice> find1SliceFromHit(hitsHandle, event, fPandoraLabel);

  for (size_t sliceIdx = 0; sliceIdx < slicesHandle->size(); ++sliceIdx) {
    const art::Ptr<recob::Slice> singleSlice(slicesHandle, sliceIdx);
    std::vector<art::Ptr<recob::PFParticle>> PFParticlesInSlice = findMPfpFromSlice.at(sliceIdx);

    for (size_t pfpIdx = 0; pfpIdx < PFParticlesInSlice.size(); pfpIdx++) {
      art::Ptr<recob::PFParticle> singlePFParticle = PFParticlesInSlice[pfpIdx];
      if (singlePFParticle->IsPrimary()) {
        const art::Ptr<larpandoraobj::PFParticleMetadata> PFPMetadata = find1MetaFromPfp.at(singlePFParticle.key());
        bool IsClearCosmic = PFPMetadata->GetPropertiesMap().count("IsClearCosmic");

        if (IsClearCosmic) {continue;}
        else { // we want only those who are NOT clear cosmic; cannot use ! operator on int
          const art::Ptr<recob::Slice> NCCSlice(slicesHandle, sliceIdx);
          outputSlices->emplace_back(NCCSlice);
          break;
        }
      }
    }
  }

  for (size_t spIdx = 0; spIdx < spacePointsHandle->size(); ++spIdx) {
    const art::Ptr<recob::SpacePoint> singleSpacePoint(spacePointsHandle, spIdx);
    const std::vector<art::Ptr<recob::Hit>> hitsAtSpacePoint = findMHitFromSP.at(spIdx);
    if (hitsAtSpacePoint.size() < 3) {continue;}
    auto slicekey0 = find1SliceFromHit.at(hitsAtSpacePoint[0].key()).key();
    auto slicekey1 = find1SliceFromHit.at(hitsAtSpacePoint[1].key()).key();
    auto slicekey2 = find1SliceFromHit.at(hitsAtSpacePoint[2].key()).key();
    if (slicekey0 != slicekey1 || slicekey0 != slicekey2) continue;
    outputSlcSpsAssns->addSingle(find1SliceFromHit.at(hitsAtSpacePoint[0].key()), singleSpacePoint);
  }

  event.put(std::move(outputSlices));
  event.put(std::move(outputSlcSpsAssns));

}

DEFINE_ART_MODULE(ICARUSNCCSlicesProducer)
