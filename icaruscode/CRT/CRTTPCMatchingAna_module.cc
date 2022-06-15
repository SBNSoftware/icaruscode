////////////////////////////////////////////////////////////////////////
// Class:       CRTTPCMatchingAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTTPCMatchingAna_module.cc
//
// Generated at Wed Feb 16 22:07:45 2022 by Biswaranjan Behera using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"


// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcorealg/CoreUtils/NumericUtils.h"

// icaruscode includes
#include "icaruscode/CRT/CRTUtils/RecoUtils.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"
//#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

// C++ includes
#include <map>
#include <vector>
#include <string>

namespace icarus {
  class CRTTPCMatchingAna;
}


class icarus::CRTTPCMatchingAna : public art::EDAnalyzer {
public:
  explicit CRTTPCMatchingAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTPCMatchingAna(CRTTPCMatchingAna const&) = delete;
  CRTTPCMatchingAna(CRTTPCMatchingAna&&) = delete;
  CRTTPCMatchingAna& operator=(CRTTPCMatchingAna const&) = delete;
  CRTTPCMatchingAna& operator=(CRTTPCMatchingAna&&) = delete;

  // Called once, at start of the job
  virtual void beginJob() override;

  // Called once per event
  virtual void analyze(const art::Event& event) override;

  // Called once, at end of the job
  virtual void endJob() override;

  // Calculate the distance from the track crossing point to CRT overlap coordinates
  double DistToCrtHit(TVector3 trackPos, sbn::crt::CRTHit crtHit);

  void reconfigure(fhicl::ParameterSet const & p);
  // Required functions.
  //  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag fCRTHitLabel;   ///< name of CRT producer
  art::InputTag fTPCTrackLabel; ///< name of CRT producer
  bool          fVerbose;             ///< print information about what's going on

  CRTT0MatchAlg t0Alg;

  //  CRTGeoAlg fCrtGeo;
  //TPCGeoAlg fTpcGeo;

  geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
  icarus::crt::CRTCommonUtils* fCrtutils;
  //  CRTCommonUtils* fCrtutils;

  // Histograms
  std::map<std::string, TH1D*> hDCA;
  std::map<std::string, TH1D*> hMatchDCA;
  std::map<std::string, TH1D*> hNoMatchDCA;

  std::map<std::string, TH1D*> hDoL;
  std::map<std::string, TH1D*> hMatchDoL;
  std::map<std::string, TH1D*> hNoMatchDoL;

  std::map<std::string, TH1D*> hT0;
  std::map<std::string, TH1D*> hMatchT0;
  std::map<std::string, TH1D*> hNoMatchT0;
};


icarus::CRTTPCMatchingAna::CRTTPCMatchingAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  , t0Alg(p.get<fhicl::ParameterSet>("t0Alg"))
  , fCrtutils(new icarus::crt::CRTCommonUtils())
  
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  reconfigure(p);
}

void icarus::CRTTPCMatchingAna::reconfigure(fhicl::ParameterSet const & p)
{
  fCRTHitLabel = (p.get<art::InputTag> ("CRTHitLabel"));
  fTPCTrackLabel= (p.get<art::InputTag> ("TPCTrackLabel"));
  fVerbose = (p.get<bool>("Verbose"));
  //fTpcTrackModuleLabel = (p.get<art::InputTag> ("TpcTrackModuleLabel"));
  // fTpcTrackModuleLabel = p.get< std::vector<art::InputTag>>("TpcTrackModuleLabel",std::vector<art::InputTag>() = {"pandoraTrackGausCryoE"});
  //fCrtHitModuleLabel   = (p.get<art::InputTag> ("CrtHitModuleLabel"));

} // CRTT0Matching::reconfigure()

void icarus::CRTTPCMatchingAna::beginJob()
{

  // Access tfileservice to handle creating and writing histograms
  art::ServiceHandle<art::TFileService> tfs;
  for(int i = 30; i < 50 + 1; i++){
    std::string tagger = "All";
    if (i < 40) continue;
    if (i==48 || i==49) continue;
    // if(i < ){
    tagger = fCrtutils->GetRegionNameFromNum(i);//fCrtGeo.GetTagger(i).name;
    std::cout << "tagger: " << tagger.c_str() << std::endl;
    hDCA[tagger]        = tfs->make<TH1D>(Form("DCA_%s", tagger.c_str()),        "", 50, 0, 100);
    hDoL[tagger]        = tfs->make<TH1D>(Form("DoL_%s", tagger.c_str()),        "", 100, 0, 0.25);
    hT0[tagger]        = tfs->make<TH1D>(Form("T0_%s", tagger.c_str()),        "", 600, -3000, 3000);
  }
}

void icarus::CRTTPCMatchingAna::analyze(const art::Event& event)
{

  // Fetch basic event info
  if(fVerbose){
    std::cout<<"============================================"<<std::endl
	     <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
	     <<"============================================"<<std::endl;
  }

  //----------------------------------------------------------------------------------------------------------
  //                                          GETTING PRODUCTS
  //----------------------------------------------------------------------------------------------------------

  // Get CRT hits from the event
  art::Handle< std::vector<sbn::crt::CRTHit>> crtHitHandle;
  std::vector<art::Ptr<sbn::crt::CRTHit> > crtHitList;
  if (event.getByLabel(fCRTHitLabel, crtHitHandle))
    art::fill_ptr_vector(crtHitList, crtHitHandle);

  // Get reconstructed tracks from the event
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);

  std::vector<sbn::crt::CRTHit> crtHits;
  int hit_i = 0;
  double minHitTime = 99999;
  double maxHitTime = -99999;

  for(auto const& hit : (*crtHitHandle)){
    double hitTime = (double)(int)hit.ts0_ns * 1e-3;
    if(hitTime < minHitTime) minHitTime = hitTime;
    if(hitTime > maxHitTime) maxHitTime = hitTime;

    crtHits.push_back(hit);
    hit_i++;
  }

  //----------------------------------------------------------------------------------------------------------
  //                                DISTANCE OF CLOSEST APPROACH ANALYSIS
  //----------------------------------------------------------------------------------------------------------

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

  // if(fVerbose) std::cout<<"----------------- DCA Analysis -------------------"<<std::endl;
  // Loop over reconstructed tracks
  for (auto const& tpcTrack : (*tpcTrackHandle)){
    //      std::cout<<"----------------- why only 4 times -------------------" <<std::endl;
    // if(fVerbose) std::cout<<"----------------- # track -------------------" << (*tpcTrackHandle).size()<<std::endl;
    // Get the associated hits
    std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
    matchCand closest = t0Alg.GetClosestCRTHit(detProp, tpcTrack, crtHits, event);
    double sin_angle = -99999;
    if(closest.dca != -99999){
      hDCA[closest.thishit.tagger]->Fill(closest.dca);
      //hDCA["All"]->Fill(closest.dca);
      sin_angle = closest.dca/closest.extrapLen;
      hDoL[closest.thishit.tagger]->Fill(sin_angle);
      // hDoL["All"]->Fill(sin_angle);
      hT0[closest.thishit.tagger]->Fill(closest.t0);
    }
  }
}

void icarus::CRTTPCMatchingAna::endJob()
{

} // CRTT0Matching::endJob()


DEFINE_ART_MODULE(icarus::CRTTPCMatchingAna)
