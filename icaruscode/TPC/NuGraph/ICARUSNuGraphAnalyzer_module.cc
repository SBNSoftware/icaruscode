
////////////////////////////////////////////////////////////////////////
// Class:       ICARUSNuGraphAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        ICARUSNuGraphAnalyzer_module.cc
//
// Riccardo Triozzi, based on the corresponding
// larrecodnn module by Giuseppe Cerati
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// saving output
#include "TTree.h"
#include "art_root_io/TFileService.h"

#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

class ICARUSNuGraphAnalyzer;

using std::vector;

class ICARUSNuGraphAnalyzer : public art::EDAnalyzer {
public:
  explicit ICARUSNuGraphAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSNuGraphAnalyzer(ICARUSNuGraphAnalyzer const&) = delete;
  ICARUSNuGraphAnalyzer(ICARUSNuGraphAnalyzer&&) = delete;
  ICARUSNuGraphAnalyzer& operator=(ICARUSNuGraphAnalyzer const&) = delete;
  ICARUSNuGraphAnalyzer& operator=(ICARUSNuGraphAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // Declare member data here.
  TTree *_treeHit, *_treeEvt;
  int _run, _subrun, _event, _id, _wire, _plane, _tpc, _cryo;
  float _x_filter, _MIP, _HIP, _shower, _michel, _diffuse, _time;
  int _islc, _icluster, _ipfp;
  float _vtx_x, _vtx_y, _vtx_z;
  std::string fNGLabel, fSliceLabel, fPandoraLabel;
};

ICARUSNuGraphAnalyzer::ICARUSNuGraphAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fNGLabel{p.get<std::string>("NuGraphLabel", "NGMultiSliceCryoE")},
  fSliceLabel{p.get<std::string>("SliceLabel", "NCCSlicesCryoE")},
  fPandoraLabel{p.get<std::string>("PandoraLabel", "pandoraGausCryoE")}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  _treeHit = tfs->make<TTree>("NuGraphHitOutput", "NuGraphHitOutput");
  _treeHit->Branch("run", &_run, "run/I");
  _treeHit->Branch("subrun", &_subrun, "subrun/I");
  _treeHit->Branch("event", &_event, "event/I");
  _treeHit->Branch("id", &_id, "id/I");
  _treeHit->Branch("wire", &_wire, "wire/I");
  _treeHit->Branch("plane", &_plane, "plane/I");
  _treeHit->Branch("tpc", &_tpc, "tpc/I");
  _treeHit->Branch("cryo", &_cryo, "cryo/I");
  _treeHit->Branch("x_filter", &_x_filter, "x_filter/F");
  _treeHit->Branch("MIP", &_MIP, "MIP/F");
  _treeHit->Branch("HIP", &_HIP, "HIP/F");
  _treeHit->Branch("shower", &_shower, "shower/F");
  _treeHit->Branch("michel", &_michel, "michel/F");
  _treeHit->Branch("diffuse", &_diffuse, "diffuse/F");
  _treeHit->Branch("time", &_time, "time/F");
  _treeHit->Branch("islc", &_islc, "islc/I");
  _treeHit->Branch("icluster", &_icluster, "icluster/I");
  _treeHit->Branch("ipfp", &_ipfp, "ipfp/I");
  _treeEvt = tfs->make<TTree>("NuGraphEventOutput", "NuGraphEventOutput");
  _treeEvt->Branch("run", &_run, "run/I");
  _treeEvt->Branch("subrun", &_subrun, "subrun/I");
  _treeEvt->Branch("event", &_event, "event/I");
  _treeEvt->Branch("vtx_x", &_vtx_x, "vtx_x/F");
  _treeEvt->Branch("vtx_y", &_vtx_y, "vtx_y/F");
  _treeEvt->Branch("vtx_z", &_vtx_z, "vtx_z/F");
}

void ICARUSNuGraphAnalyzer::analyze(art::Event const& e)
{

  // art::ValidHandle<std::vector<anab::FeatureVector<1>>> filterHandle = e.getValidHandle<std::vector<anab::FeatureVector<1>>>(fNGLabel);
  // art::ValidHandle<std::vector<anab::FeatureVector<5>>> semanticHandle = e.getValidHandle<std::vector<anab::FeatureVector<5>>>(fNGLabel);

  // get slices
  const std::vector<art::Ptr<recob::Slice>> slices = e.getProduct<vector<art::Ptr<recob::Slice>>>(fSliceLabel);
  art::FindManyP<recob::Hit> findMHitsFromSlice(slices, e, fPandoraLabel);

  // map hits to slices
  std::vector<art::Ptr<recob::Hit>> allHits;
  std::map<unsigned int, int> hitToSliceID;
  for (size_t islc = 0; islc < slices.size(); ++islc) {
    art::Ptr<recob::Slice> slice = slices[islc];
    const std::vector<art::Ptr<recob::Hit>> hitsInSlice = findMHitsFromSlice.at(islc);
    for (auto const& h : hitsInSlice) {
      hitToSliceID[h.key()] = slice.key();
      allHits.emplace_back(h);
    }
  }

  std::cout << "allHits size: " << allHits.size() << '\n';

  art::InputTag fNGFilterLabel{fNGLabel, "filter"};
  art::InputTag fNGSemanticLabel{fNGLabel, "semantic"};

  art::FindOneP<anab::FeatureVector<1>> find1FilterFromHit(allHits, e, fNGFilterLabel);
  art::FindOneP<anab::FeatureVector<5>> find1SemanticFromHit(allHits, e, fNGSemanticLabel);

  std::cout << "allHit keys " << '\n';
  for (auto hit : allHits) {
    std::cout << ' ' << hit.key();
  }

  for (size_t hitIdx = 0; hitIdx < allHits.size(); hitIdx++) {
    art::Ptr<recob::Hit> hit = allHits[hitIdx];
    std::cout << "Begin hit #" << hitIdx+1 << " with key #" << hit.key() << ' ';
    // event information
    _event  = e.event();
    _subrun = e.subRun();
    _run    = e.run();
    _id     = hit.key();

    // NuGraph predictions
    const art::Ptr<anab::FeatureVector<1>> filter = find1FilterFromHit.at(hitIdx);
    const art::Ptr<anab::FeatureVector<5>> semantic = find1SemanticFromHit.at(hitIdx);
    _x_filter = filter->at(0);
    _MIP      = semantic->at(0);
    _HIP      = semantic->at(1);
    _shower   = semantic->at(2);
    _michel   = semantic->at(3);
    _diffuse  = semantic->at(4);

    // hit description
    _wire  = hit->WireID().Wire;
    _plane = hit->WireID().Plane;
    _tpc   = hit->WireID().TPC;
    _cryo  = hit->WireID().Cryostat;
    _time  = hit->PeakTime();

    auto itSlice = hitToSliceID.find(hit.key());
    _islc = (itSlice != hitToSliceID.end()) ? itSlice->second : -1;

    _treeHit->Fill();
    std::cout << "end" << '\n';
  }
}

DEFINE_ART_MODULE(ICARUSNuGraphAnalyzer)
