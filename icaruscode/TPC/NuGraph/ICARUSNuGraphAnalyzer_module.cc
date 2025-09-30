
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
  std::string fNGLabel, fHitLabel;
};

ICARUSNuGraphAnalyzer::ICARUSNuGraphAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, 
  fNGLabel{p.get<std::string>("NuGraphLabel", "NuGraph")},
  fHitLabel{p.get<std::string>("HitLabel", "nuslhitsCryoE")}
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

  // get hits from the tagged slice
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  std::cout << hitListHandle->size() << std::endl;

  // get NuGraph predictions
  art::Handle<std::vector<anab::FeatureVector<1>>> filterHandle;
  e.getByLabel(art::InputTag(fNGLabel, "filter"), filterHandle);

  art::Handle<std::vector<anab::FeatureVector<5>>> semanticHandle;
  e.getByLabel(art::InputTag(fNGLabel, "semantic"), semanticHandle);

  for (size_t ihit = 0; ihit < hitListHandle->size(); ihit++) {
    art::Ptr<recob::Hit> hit(hitListHandle, ihit);

    // event information
    _event  = e.event();
    _subrun = e.subRun();
    _run    = e.run();
    _id     = hit.key();

    // NuGraph predictions
    _x_filter = filterHandle->at(ihit).at(0);
    _MIP      = semanticHandle->at(ihit).at(0);
    _HIP      = semanticHandle->at(ihit).at(1);
    _shower   = semanticHandle->at(ihit).at(2);
    _michel   = semanticHandle->at(ihit).at(3);
    _diffuse  = semanticHandle->at(ihit).at(4);

    // hit description
    _wire  = hit->WireID().Wire;
    _plane = hit->WireID().Plane;
    _tpc   = hit->WireID().TPC;
    _cryo  = hit->WireID().Cryostat;
    _time  = hit->PeakTime();

    _treeHit->Fill();
  }

   auto PredVertexColl = e.getHandle<std::vector<recob::Vertex>>(art::InputTag(fNGLabel, "vertex"));
   if (PredVertexColl.isValid() && PredVertexColl->size() > 0) { //there should be only one
      _vtx_x = PredVertexColl->at(0).position().X();
      _vtx_y = PredVertexColl->at(0).position().Y();
      _vtx_z = PredVertexColl->at(0).position().Z();
      _treeEvt->Fill();
   }
}

DEFINE_ART_MODULE(ICARUSNuGraphAnalyzer)
