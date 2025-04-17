////////////////////////////////////////////////////////////////////////
// Class:       IcarusTrueNuSliceHitsProducer
// Plugin Type: producer (art v3_06_03)
// File:        IcarusTrueNuSliceHitsProducer_module.cc
//
// Generated at Tue May 25 10:39:19 2021 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_11_01.
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
#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

class IcarusTrueNuSliceHitsProducer;
class IcarusTrueNuSliceHitsProducer : public art::EDProducer {
public:
  explicit IcarusTrueNuSliceHitsProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  IcarusTrueNuSliceHitsProducer(IcarusTrueNuSliceHitsProducer const&) = delete;
  IcarusTrueNuSliceHitsProducer(IcarusTrueNuSliceHitsProducer&&) = delete;
  IcarusTrueNuSliceHitsProducer& operator=(IcarusTrueNuSliceHitsProducer const&) = delete;
  IcarusTrueNuSliceHitsProducer& operator=(IcarusTrueNuSliceHitsProducer&&) = delete;


private:

  // Declare member data here.
  std::string fBaryMatchLabel;
  std::string fSliceLabel;
  std::string fSpsLabel;
  bool doFlashAssn;

  std::vector<int> nearwires;

  // Required functions.
  void produce(art::Event& e) override;

};


IcarusTrueNuSliceHitsProducer::IcarusTrueNuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fBaryMatchLabel( p.get<std::string>("BaryMatchLabel","tpcpmtbarycentermatchCryoE")),
  fSliceLabel(p.get<std::string>("SliceLabel","pandoraGausCryoE")),
  fSpsLabel(p.get<std::string>("SpsLabel","cluster3DCryoE")),
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

void IcarusTrueNuSliceHitsProducer::produce(art::Event& e)
{

  auto outputHits = std::make_unique<std::vector<recob::Hit> >();
  auto assocSliceHitKeys = std::make_unique<std::vector<unsigned int> >();

  art::PtrMaker<recob::Hit> hitPtrMaker(e);
  art::PtrMaker<recob::SpacePoint> spsPtrMaker(e);

  auto outputSpacepoints = std::make_unique<std::vector<recob::SpacePoint> >();
  auto outputSpsHitAssns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
  auto outputSlcFlhAssns = std::make_unique<art::Assns<recob::Slice, recob::OpFlash>>();

  art::ServiceHandle<cheat::ParticleInventoryService> pi; 
  art::ServiceHandle<cheat::BackTrackerService> bt_h;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(fSliceLabel);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputSlice, e, fSliceLabel));
  int slkey = -1;
  int maxCnt = 0;
  for (size_t sk=0; sk<inputSlice->size(); sk++) {
    int slcnt = 0;
    for (auto hit : assocSliceHit->at(sk)) {
      std::vector<sim::TrackIDE> h_ides = bt_h->ChannelToTrackIDEs(clockData, hit->Channel(), hit->StartTick(), hit->EndTick());
      for (auto& tide : h_ides) {
      	int tid = tide.trackID;
	auto mct = pi->TrackIdToMCTruth_P(abs(tid));
	int fromNu = (mct.isAvailable() ? mct->NeutrinoSet() : 0);
	if (fromNu==1) {
	  slcnt++;
	  break;
	}
      }
      if (slcnt>maxCnt) {
	slkey = sk;
	maxCnt = slcnt;
      }
    }
  }
  std::cout << "keys: slice=" << slkey << std::endl;

  if (slkey>=0) {
    std::cout << "slice hits=" << assocSliceHit->at(slkey).size() << std::endl;
    for (auto h : assocSliceHit->at(slkey)) {
      assocSliceHitKeys->push_back(h.key());
      outputHits->emplace_back(*h);  
    }
  }

  // Get spacepoints from the event record
  art::Handle<std::vector<recob::SpacePoint>> spListHandle;
  std::vector<art::Ptr<recob::SpacePoint>> splist;
  if (e.getByLabel(fSpsLabel, spListHandle)) art::fill_ptr_vector(splist, spListHandle);
  size_t countsps = 0;
  // Get assocations from spacepoints to hits
  art::FindManyP<recob::Hit> fmp(spListHandle, e, fSpsLabel);
  for (size_t spIdx = 0; spIdx < splist.size(); ++spIdx) {
    // if (splist[spIdx]->Chisq()>0.5) continue;
    auto assochits = fmp.at(spIdx);
    if (assochits.size()<3) continue;
    //
    std::vector<size_t> hitpos;
    for (size_t j = 0; j < assochits.size(); ++j) {
      auto pos = std::find(assocSliceHitKeys->begin(),assocSliceHitKeys->end(),assochits[j].key());
      if (pos==assocSliceHitKeys->end()) break;
      hitpos.push_back(pos-assocSliceHitKeys->begin());
    }
    if (hitpos.size()<3) continue;
    outputSpacepoints->emplace_back(*splist[spIdx]);
    //
    const art::Ptr<recob::SpacePoint> sps = spsPtrMaker(outputSpacepoints->size()-1);
    for (size_t j = 0; j < hitpos.size(); ++j) {
      const art::Ptr<recob::Hit> ahp = hitPtrMaker(hitpos[j]);
      outputSpsHitAssns->addSingle(ahp,sps);
    }
    countsps++;
  } // for spacepoint
  std::cout << "sps count=" << countsps << std::endl;

  art::ValidHandle< std::vector< sbn::TPCPMTBarycenterMatch > > baryMatchListHandle = e.getValidHandle<std::vector<sbn::TPCPMTBarycenterMatch> >(fBaryMatchLabel);
  art::FindManyP<recob::Slice> msl(baryMatchListHandle, e, fBaryMatchLabel);
  art::FindManyP<recob::OpFlash> mfl(baryMatchListHandle, e, fBaryMatchLabel);
  float minRad = 99999.;
  int mbkey = -1;
  for (size_t ibm=0; ibm<baryMatchListHandle->size(); ++ibm) {
    art::Ptr<sbn::TPCPMTBarycenterMatch> bm(baryMatchListHandle,ibm);
    if (int(msl.at(ibm).at(0).key())!=slkey) continue;
    if (mfl.at(ibm).size()>0) {
      std::cout << "ibm=" << ibm << " radius=" << bm->radius << " radius_Trigger=" << bm->radius_Trigger << " flashTime=" << bm->flashTime 
      		<< " slkey=" << msl.at(ibm).at(0).key() << " center=" << bm->chargeCenter
      		<< " flkey=" << mfl.at(ibm).at(0).key() << " center=" << bm->flashCenter
      		<< std::endl;
      if (bm->radius<minRad) {
	//this is the triggering flash
	mbkey = ibm;
	minRad = bm->radius;
      }
    }
  }
  std::cout << "mbkey=" << mbkey << std::endl;
  if (doFlashAssn){
    if (slkey>=0 && mbkey>=0) {
      art::Ptr<recob::Slice> slp(inputSlice,slkey);
      outputSlcFlhAssns->addSingle(slp,mfl.at(mbkey).at(0));
    }
  }

  e.put(std::move(outputHits));
  e.put(std::move(assocSliceHitKeys));
  e.put(std::move(outputSpacepoints));
  e.put(std::move(outputSpsHitAssns));
  if (doFlashAssn) e.put(std::move(outputSlcFlhAssns));
}

DEFINE_ART_MODULE(IcarusTrueNuSliceHitsProducer)
