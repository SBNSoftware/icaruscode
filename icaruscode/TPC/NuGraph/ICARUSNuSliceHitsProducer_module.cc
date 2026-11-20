/**
 * @file   icaruscode/TPC/NuGraph/ICARUSNuSliceHitsProducer_module.cc
 * @brief  Implementation of `ICARUSNuSliceHitsProducer` _art_ module.
 * @author Giuseppe Cerati (cerati@fnal.gov)
 * @date   May 25, 2021
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

#include <algorithm> // std::find()
#include <cmath> // std::abs()
#include <memory>
#include <unordered_map>
#include <utility> // std::move()
#include <vector>

#include "lardataobj/RecoBase/OpFlash.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

class ICARUSNuSliceHitsProducer : public art::EDProducer {

/**
 * @brief Produce separate collections of hits and space point from the slice identified as most likely from the neutrino interaction.
 *
 * Produce separate collections of hits and space point from the slice identified as most likely from the neutrino interaction.
 * The slice is the one identified as the one producing the trigger flash. Only space points made out of 3 hits are saved.
 * It also produces a vector of ints, that maps the new hit collection to the original one.
 * Optionally, it can store the association to the trigger flash used to pick the slice.
 *
 */

public:
  explicit ICARUSNuSliceHitsProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSNuSliceHitsProducer(ICARUSNuSliceHitsProducer const&) = delete;
  ICARUSNuSliceHitsProducer(ICARUSNuSliceHitsProducer&&) = delete;
  ICARUSNuSliceHitsProducer& operator=(ICARUSNuSliceHitsProducer const&) = delete;
  ICARUSNuSliceHitsProducer& operator=(ICARUSNuSliceHitsProducer&&) = delete;


private:

  // Declare member data here.
  art::InputTag const fBaryMatchLabel;
  art::InputTag const fSliceLabel;
  art::InputTag const fSpsLabel;
  bool const doFlashAssn;

  // Required functions.
  void produce(art::Event& e) override;

};


ICARUSNuSliceHitsProducer::ICARUSNuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fBaryMatchLabel( p.get<art::InputTag>("BaryMatchLabel","tpcpmtbarycentermatchCryoE")),
  fSliceLabel(p.get<art::InputTag>("SliceLabel","pandoraGausCryoE")),
  fSpsLabel(p.get<art::InputTag>("SpsLabel","cluster3DCryoE")),
  doFlashAssn(p.get<bool>("DoFlashAssn",false))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
  produces<std::vector<unsigned int>>();
  produces<std::vector<recob::SpacePoint>>();
  produces<art::Assns<recob::Hit, recob::SpacePoint>>();
  if (doFlashAssn) produces<art::Assns<recob::Slice, recob::OpFlash>>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ICARUSNuSliceHitsProducer::produce(art::Event& e)
{

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  auto assocSliceHitKeys = std::make_unique<std::vector<unsigned int> >();

  art::PtrMaker<recob::Hit> hitPtrMaker(e);
  art::PtrMaker<recob::SpacePoint> spsPtrMaker(e);

  auto outputSpacepoints = std::make_unique<std::vector<recob::SpacePoint> >();
  auto outputSpsHitAssns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
  auto outputSlcFlhAssns = std::make_unique<art::Assns<recob::Slice, recob::OpFlash>>();

  art::ValidHandle< std::vector< sbn::TPCPMTBarycenterMatch > > baryMatchListHandle = e.getValidHandle<std::vector<sbn::TPCPMTBarycenterMatch> >(fBaryMatchLabel);
  art::FindOneP<recob::Slice> msl(baryMatchListHandle, e, fBaryMatchLabel);
  art::FindOneP<recob::OpFlash> mfl(baryMatchListHandle, e, fBaryMatchLabel);
  float minRad = 99999.;
  art::Ptr<recob::Slice> triggeringSlice;
  art::Ptr<recob::OpFlash> triggeringFlash;
  for (size_t ibm=0; ibm<baryMatchListHandle->size(); ++ibm) {
    sbn::TPCPMTBarycenterMatch const& bm = (*baryMatchListHandle)[ibm];
    // require the radius to be set (i.e. that there is a best match) and to match the one of the trigger within an arbitrarily small enough margin
    if (bm.radius_Trigger<=0 || std::abs(bm.radius-bm.radius_Trigger)>=0.00001) continue;
    art::Ptr<recob::Slice> const& slice = msl.at(ibm);
    art::Ptr<recob::OpFlash> const& flash = mfl.at(ibm);
    if (!slice || !flash) throw std::logic_error("Unexpected missing flash or slice in barycenter matching associations");
    mf::LogTrace{ "ICARUSNuSliceHitsProducer" }
      << "ibm=" << ibm << " radius=" << bm.radius << " radius_Trigger=" << bm.radius_Trigger << " flashTime=" << bm.flashTime
      << " slkey=" << msl.at(ibm).key() << " center=" << bm.chargeCenter
      << " flkey=" << mfl.at(ibm).key() << " center=" << bm.flashCenter;
    if (bm.radius_Trigger<minRad) {
      minRad = bm.radius_Trigger;
      triggeringSlice = slice;
      triggeringFlash = flash;
    }
  }
  mf::LogTrace{ "ICARUSNuSliceHitsProducer" } << "keys: slice=" << triggeringSlice.key() << " flash=" << triggeringFlash.key();

  if (doFlashAssn) {
    if (triggeringSlice && triggeringFlash) {
      outputSlcFlhAssns->addSingle(triggeringSlice,triggeringFlash);
    }
  }

  art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(fSliceLabel);
  art::FindManyP<recob::Hit> assocSliceHit(inputSlice, e, fSliceLabel);

  if (triggeringSlice) {
    mf::LogTrace{ "ICARUSNuSliceHitsProducer" } << "slice hits=" << assocSliceHit.at(triggeringSlice.key()).size();
    auto const& hits = assocSliceHit.at(triggeringSlice.key());
    std::unordered_map<size_t,size_t> keyIdxMap;
    for (size_t ih=0; ih<hits.size(); ih++) {
      assocSliceHitKeys->push_back(hits[ih].key());
      keyIdxMap.insert({hits[ih].key(), ih});
    }
    std::sort(begin(*assocSliceHitKeys), end(*assocSliceHitKeys));
    outputHits->reserve(assocSliceHitKeys->size());
    for (auto const& key: *assocSliceHitKeys) outputHits->push_back(*hits[keyIdxMap.at(key)]);
  }

  // Get spacepoints from the event record
  auto splist = e.getValidHandle<std::vector<recob::SpacePoint>>(fSpsLabel);
  size_t countsps = 0;
  // Get assocations from spacepoints to hits
  art::FindManyP<recob::Hit> spacePointToHits(splist, e, fSpsLabel);
  for (size_t spIdx = 0; spIdx < splist->size(); ++spIdx) {
    auto const& assochits = spacePointToHits.at(spIdx);
    // Consider only space points with hits associated on all planes. This is enough for NuGraph.
    if (assochits.size()<3) continue;
    //if (splist[spIdx]->Chisq()>0.5) continue;
    //
    std::vector<size_t> hitpos;
    for (size_t j = 0; j < assochits.size(); ++j) {
      auto pos = std::lower_bound(assocSliceHitKeys->begin(),assocSliceHitKeys->end(),assochits[j].key());
      if ( (pos == assocSliceHitKeys->end()) || (*pos != assochits[j].key()) ) break;
      hitpos.push_back(pos-assocSliceHitKeys->begin());
    }
    if (hitpos.size()<3) continue;
    outputSpacepoints->emplace_back((*splist)[spIdx]);
    //
    const art::Ptr<recob::SpacePoint> sps = spsPtrMaker(outputSpacepoints->size()-1);
    for (size_t j = 0; j < hitpos.size(); ++j) {
      const art::Ptr<recob::Hit> ahp = hitPtrMaker(hitpos[j]);
      outputSpsHitAssns->addSingle(ahp,sps);
    }
    countsps++;
  } // for spacepoint
  mf::LogTrace{ "ICARUSNuSliceHitsProducer" } << "sps count=" << countsps;

  e.put(std::move(outputHits));
  e.put(std::move(assocSliceHitKeys));
  e.put(std::move(outputSpacepoints));
  e.put(std::move(outputSpsHitAssns));
  if (doFlashAssn) e.put(std::move(outputSlcFlhAssns));

}

DEFINE_ART_MODULE(ICARUSNuSliceHitsProducer)
