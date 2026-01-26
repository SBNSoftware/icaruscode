
////////////////////////////////////////////////////////////////////////
// Class:       ICARUSPandoraHitDumpAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        ICARUSPandoraHitDumpAnalyzer_module.cc
//
// Riccardo Triozzi
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
#include "lardataobj/RecoBase/PFParticleMetadata.h"

class ICARUSPandoraHitDumpAnalyzer;

using std::vector;

class ICARUSPandoraHitDumpAnalyzer : public art::EDAnalyzer {
public:
  explicit ICARUSPandoraHitDumpAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSPandoraHitDumpAnalyzer(ICARUSPandoraHitDumpAnalyzer const&) = delete;
  ICARUSPandoraHitDumpAnalyzer(ICARUSPandoraHitDumpAnalyzer&&) = delete;
  ICARUSPandoraHitDumpAnalyzer& operator=(ICARUSPandoraHitDumpAnalyzer const&) = delete;
  ICARUSPandoraHitDumpAnalyzer& operator=(ICARUSPandoraHitDumpAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // Declare member data here.
  TTree *_treeHit, *_treeEvt;

  // hit information
  int _run, _subrun, _event, _id, _wire, _plane, _tpc, _cryo;

  // NuGraph2 information
  float _x_filter, _MIP, _HIP, _shower, _michel, _diffuse, _time;

  // Pandora information
  int _islc, _icluster, _ipfp, _pdg_hit;
  float _ipfpslc, _vtx_x, _vtx_y, _vtx_z, _pdg;
  std::string fHitLabel, fPandoraLabel;
};

ICARUSPandoraHitDumpAnalyzer::ICARUSPandoraHitDumpAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, 
  fHitLabel{p.get<std::string>("HitLabel", "cluster3DCryoE")},
  fPandoraLabel{p.get<std::string>("PandoraLabel", "pandoraGausCryoE")}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  _treeHit = tfs->make<TTree>("PandoraHitOutput", "PandoraHitOutput");
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
  _treeHit->Branch("pfppdg", &_pdg_hit, "pfppdg/I");    
  _treeEvt = tfs->make<TTree>("PandoraPFPOutput", "PandoraPFPOutput");
  _treeEvt->Branch("run", &_run, "run/I");
  _treeEvt->Branch("subrun", &_subrun, "subrun/I");
  _treeEvt->Branch("event", &_event, "event/I");
  _treeEvt->Branch("ipfpslc", &_ipfpslc, "ipfpslc/I");
  _treeEvt->Branch("vtx_x", &_vtx_x, "vtx_x/F");
  _treeEvt->Branch("vtx_y", &_vtx_y, "vtx_y/F");
  _treeEvt->Branch("vtx_z", &_vtx_z, "vtx_z/F");
  _treeEvt->Branch("pdg", &_pdg, "pdg/I");
}

void ICARUSPandoraHitDumpAnalyzer::analyze(art::Event const& e)
{
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  std::cout << hitListHandle->size() << std::endl;

  // map hits to Pandora slices
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fPandoraLabel, sliceHandle);
  art::FindManyP<recob::Hit> sliceAssoc(sliceHandle, e, fPandoraLabel);

  std::map<unsigned int, int> hitToSliceID;
  for (size_t islc = 0; islc < sliceHandle->size(); ++islc) {
    art::Ptr<recob::Slice> slice(sliceHandle, islc);
    auto const& hitsInSlice = sliceAssoc.at(islc);
    for (auto const& h : hitsInSlice) 
      hitToSliceID[h.key()] = slice->ID();
  }

  // map hits to Pandora clusters
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  e.getByLabel(fPandoraLabel, clusterHandle);
  art::FindManyP<recob::Hit> clusterAssoc(clusterHandle, e, fPandoraLabel);

  std::map<unsigned int, int> hitToClusterID;
  for (size_t iclus = 0; iclus < clusterHandle->size(); ++iclus) {
    art::Ptr<recob::Cluster> cluster(clusterHandle, iclus);
    auto const& hitsInCluster = clusterAssoc.at(iclus);
    for (auto const& h : hitsInCluster) 
      // hitToClusterID[h.key()] = cluster->ID();
      hitToClusterID[h.key()] = cluster.key();
  }

  // Pandora PFPs
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPandoraLabel, pfpHandle);
  art::FindManyP<recob::PFParticle> pfpAssoc(clusterHandle, e, fPandoraLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssoc(pfpHandle, e, fPandoraLabel);

  std::map<unsigned int, int> pfpIDMap;
  for (size_t ipfp = 0; ipfp < pfpHandle->size(); ++ipfp) {
    art::Ptr<recob::PFParticle> pfp(pfpHandle, ipfp);
    pfpIDMap[pfp.key()] = ipfp; 
  }

  // Pandora clusters
  std::map<unsigned int, int> clusterToPFPID;

  for (size_t iclus = 0; iclus < clusterHandle->size(); ++iclus) {
    art::Ptr<recob::Cluster> cluster(clusterHandle, iclus);
    auto const& pfps = pfpAssoc.at(iclus);

    auto it = pfpIDMap.find(pfps.front().key());
    if (it != pfpIDMap.end()) 
      clusterToPFPID[cluster.key()] = it->second;
  }

  // fill PFP-level tree
  for (size_t ihit = 0; ihit < hitListHandle->size(); ihit++) {
    art::Ptr<recob::Hit> hit(hitListHandle, ihit);

    // event information
    _event  = e.event();
    _subrun = e.subRun();
    _run    = e.run();
    _id     = hit.key();

    // hit description
    _wire  = hit->WireID().Wire;
    _plane = hit->WireID().Plane;
    _tpc   = hit->WireID().TPC;
    _cryo  = hit->WireID().Cryostat;
    _time  = hit->PeakTime();

    // map to Pandora information
    auto itSlice = hitToSliceID.find(hit.key());
    _islc = (itSlice != hitToSliceID.end()) ? itSlice->second : -1;

    auto itCluster = hitToClusterID.find(hit.key());
    _icluster = (itCluster != hitToClusterID.end()) ? itCluster->second : -1;

    _ipfp = -1;
    _pdg_hit = -1;

    if (_icluster > -1) {
      auto itPFP = clusterToPFPID.find(_icluster);
      if (itPFP != clusterToPFPID.end()) {
        _ipfp = itPFP->second; 
        if (_ipfp > -1 && _ipfp < (int)pfpHandle->size()) {
          art::Ptr<recob::PFParticle> pfp(pfpHandle, _ipfp);
          _pdg_hit = pfp->PdgCode();
        }
      }
    }

    _treeHit->Fill();
  }

  // fill PFP-level tree
  art::FindManyP<recob::Slice> slcToPFPAssoc(pfpHandle, e, fPandoraLabel);
  art::FindManyP<recob::Vertex> vtxToPFPAssoc(pfpHandle, e, fPandoraLabel);
  for (size_t ipfp = 0; ipfp < pfpHandle->size(); ipfp++) {
    art::Ptr<recob::PFParticle> pfp(pfpHandle, ipfp);

    // slice index
    std::vector<art::Ptr<recob::Slice>> slcList = slcToPFPAssoc.at(pfp.key());
    if (slcList.size()) {
      art::Ptr<recob::Slice> slc = slcList[0];
      _ipfpslc = slc->ID();
    }
    
    // vertex
    std::vector<art::Ptr<recob::Vertex>> vtxList = vtxToPFPAssoc.at(pfp.key());
    if (vtxList.size()) {
      art::Ptr<recob::Vertex> vtx = vtxList[0];
      _vtx_x = vtx->position().X();
      _vtx_y = vtx->position().Y();
      _vtx_z = vtx->position().Z();
    }
    else {
      _vtx_x = _vtx_y = _vtx_z = -1.;
    }

    // truth information
    _pdg = pfp->PdgCode();

    _treeEvt->Fill();
    
  }
}

DEFINE_ART_MODULE(ICARUSPandoraHitDumpAnalyzer)