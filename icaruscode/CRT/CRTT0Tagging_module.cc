/////////////////////////////////////////////////////////////////////////////
/// Class:       CRTT0Tagging
/// Module Type: producer
/// File:        CRTT0Tagging_module.cc
///
/// Author:         Francesco Poppi
/// E-mail address: poppi@bo.infn.it 
/// October 2024
///
/////////////////////////////////////////////////////////////////////////////

#include "sbnobj/Common/CRT/CRTHit.hh"
//#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "icaruscode/CRT/CRTUtils/CRTMatchingUtils.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
// ROOT
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TTree.h"

namespace icarus {
namespace crt {

class CRTT0Tagging;

}  // namespace crt
}  // namespace icarus

using namespace icarus::crt;

class icarus::crt::CRTT0Tagging : public art::EDProducer {
public:
  
  using CRTHit = sbn::crt::CRTHit;
  using CRTPMTMatching = sbn::crt::CRTPMTMatching;

  explicit CRTT0Tagging(fhicl::ParameterSet const& p);
  
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTT0Tagging(CRTT0Tagging const&) = delete;
  CRTT0Tagging(CRTT0Tagging&&) = delete;
  CRTT0Tagging& operator=(CRTT0Tagging const&) = delete;
  CRTT0Tagging& operator=(CRTT0Tagging&&) = delete;

  // Required functions.
  void beginRun(art::Run& r) override;
  void produce(art::Event& e) override;
  void endJob() override;

private:

  // Declare member data here.

  art::InputTag fCrtHitModuleLabel;
  art::InputTag fTriggerLabel;
  art::InputTag fTriggerConfigurationLabel;
  
  std::vector<art::InputTag> fTPCTrackLabel; ///< labels for source of tracks
  std::vector<art::InputTag> fPFParticleLabel; ///< labels for source of PFParticle
  std::vector<art::InputTag> fHitLabel; ///< labels for source of hits
  art::InputTag fTRKHMproducer; ///< labels for hit metadata
  
  std::optional<icarus::TriggerConfiguration> fTriggerConfiguration;

  geo::GeometryCore const* fGeometryService;  ///< pointer to Geometry provider
  CRTCommonUtils* fCrtUtils; 
  icarus::crt::TopCRTCentersMap fTopCRTCenterMap;
  icarus::crt::TopCrtTransformations fTopCrtTransformations;
  CRTMatchingAlg fMatchingAlg;
  double fMinimalTrackLength;
  int fMinimumGoodHits;
  double fMaximalCRTDistance;
  double fGoodCandidateDistance;

};

icarus::crt::CRTT0Tagging::CRTT0Tagging(fhicl::ParameterSet const& p)
    : EDProducer{p},
      fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit")),
      fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")),
      fTriggerConfigurationLabel(
          p.get<art::InputTag>("TriggerConfiguration", "triggerconfig")),
      fTPCTrackLabel(p.get< std::vector<art::InputTag> >("TPCTrackLabel",             {""})),
      fPFParticleLabel(p.get< std::vector<art::InputTag> >("PFParticleLabel",             {""})),  
      fHitLabel(p.get< std::vector<art::InputTag> >("HitLabel",             {""})),
      fTRKHMproducer(p.get< art::InputTag   > ("TRKHMproducer", "")),
      fCrtUtils(new CRTCommonUtils()),
      fMatchingAlg(p.get<fhicl::ParameterSet> ("MatchingAlg"))
{
  
  produces< std::vector<anab::T0>                   >();
  produces< art::Assns<recob::Track , anab::T0>     >();
  produces< art::Assns<sbn::crt::CRTHit, anab::T0>  >();

  // Get a pointer to the geometry service provider.
  fGeometryService = lar::providerFrom<geo::Geometry>();
  fMinimalTrackLength = p.get<double>("MinimalTrackLength", 40.0);
  fMinimumGoodHits = p.get<double>("MinimumGoodHits", 5);
  fMaximalCRTDistance = p.get<double>("MaximalCRTDistance", 300.);
  fGoodCandidateDistance = p.get<double>("GoodCandidateDistance", 100.);

  art::ServiceHandle<art::TFileService> tfs;

}

void icarus::crt::CRTT0Tagging::beginRun(art::Run& r)
{
  // we don't know if this is data or not; if not, there will be no trigger config
  auto const& trigConfHandle = 
    r.getHandle<icarus::TriggerConfiguration>(fTriggerConfigurationLabel);
  
  fTriggerConfiguration
    = trigConfHandle.isValid()? std::make_optional(*trigConfHandle): std::nullopt;

  fTopCRTCenterMap=icarus::crt::LoadTopCRTCenters();
  fTopCrtTransformations=icarus::crt::LoadTopCrtTransformations();
}

void icarus::crt::CRTT0Tagging::produce(art::Event& e)
{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);
  if (!fTriggerConfiguration) {
    mf::LogDebug("CRTT0Tagging")
      << "Skipping because no data (or at least no trigger configuration).";
    return;
  }

  mf::LogDebug("CRTT0Tagging: ") << "beginning production" << '\n';

  std::unique_ptr< std::vector<anab::T0> > t0col( new std::vector<anab::T0>);
  std::unique_ptr< art::Assns<recob::Track, anab::T0> > trackAssn( new art::Assns<recob::Track, anab::T0>);
  std::unique_ptr< art::Assns <sbn::crt::CRTHit, anab::T0> > t0CrtHitAssn( new art::Assns<sbn::crt::CRTHit, anab::T0> );

  // CRTHits
  std::vector<art::Ptr<CRTHit>> CRTHitList;
  art::ValidHandle<std::vector<CRTHit>> crthits = e.getValidHandle<std::vector<CRTHit>>(fCrtHitModuleLabel);
  art::fill_ptr_vector(CRTHitList, crthits);

  int cryo=-1;
  for(const auto& PFPLabel : fPFParticleLabel) {
    cryo++;
    auto it = &PFPLabel - fPFParticleLabel.data();
    std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
    art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(PFPLabel);
    art::fill_ptr_vector(PFParticleList, pfparticles);

    // Pandora MetaData
    art::FindManyP<anab::T0> fmt0pandora(pfparticles, e, PFPLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> PFPMetaDataAssoc(pfparticles, e, PFPLabel);

    // Tracks
    art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel[it]);

    // Track - associated data
    art::FindManyP<recob::Track> fmTracks(PFParticleList, e, fTPCTrackLabel[it]);

    // Collect all hits
    art::ValidHandle<std::vector<recob::Hit>> allhit_handle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel[it]);
    std::vector<art::Ptr<recob::Hit>> allHits;
    art::fill_ptr_vector(allHits, allhit_handle);

    // Start looping on the particles
    for (art::Ptr<recob::PFParticle> p_pfp: PFParticleList) {
      const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(p_pfp.key());
      if (thisTrack.size() != 1) continue;    
      art::Ptr<recob::Track> trkPtr = thisTrack.at(0);
      const recob::Track &track = *trkPtr;
      art::InputTag thm_label = fTRKHMproducer.empty() ? fTPCTrackLabel[it] : fTRKHMproducer;
      art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e, thm_label);
      std::vector<art::Ptr<recob::Hit>> emptyHitVector;
      const std::vector<art::Ptr<recob::Hit>> &trkHits  = fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;
      std::vector<const recob::TrackHitMeta*> emptyTHMVector;
      const std::vector<const recob::TrackHitMeta*> &trkHitMetas = fmtrkHits.isValid() ? fmtrkHits.data(trkPtr.key()) : emptyTHMVector;
      
      if(track.Length()<fMinimalTrackLength) continue;

      // T0
      float t0 = std::numeric_limits<float>::signaling_NaN();
      auto t0s = fmt0pandora.at(p_pfp.key());
      if (!t0s.empty()){  
        t0 = t0s[0]->Time();   //Get T0  
	    }
      if(!isnan(t0)) std::cout<<"This Track T0 "<<t0<<std::endl;
      int goodHits=0;
      int countE=0, countW=0;
      // These counters are used to determine if track is CC-E, EE, EW, CC-W, WE, WW
      // depending on the track type, the Top CRT uses the appropriate position corretions
      double trackRes=0, deltaP=0;

      std::vector<float> hx, hy, hz, ht;
      for(size_t i=0; i<trkHits.size(); i++){
        bool badhit = (trkHitMetas[i]->Index() == std::numeric_limits<unsigned int>::max()) ||
                    (!track.HasValidPoint(trkHitMetas[i]->Index()));
        geo::Point_t loc = track.LocationAtPoint(trkHitMetas[i]->Index());
        if(loc.X()==-999) continue;
        if(badhit) continue;        
        hx.push_back(loc.X()); hy.push_back(loc.Y()); hz.push_back(loc.Z());
        ht.push_back(trkHits[i]->PeakTime());
        goodHits++;
        if(trkHits[i]->WireID().TPC==0 || trkHits[i]->WireID().TPC==1) countE++;
        else countW++;

      }
      int trackType=-1;
      if(countW!=0 && countE!=0 && cryo==0) trackType=0; //CCEast
      else if(countW!=0 && countE==0 && cryo==0) trackType=2; //East-West
      else if(countW==0 && countE!=0 && cryo==0) trackType=1; //East-East      
      else if(countW!=0 && countE!=0 && cryo==1) trackType=3; //CCWest
      else if(countW!=0 && countE==0 && cryo==1) trackType=5; //West-West
      else if(countW==0 && countE!=0 && cryo==1) trackType=4; //West-East      
      icarus::crt::TopCRTCorrectionMap TopCrtCorrection;
      if(trackType==0) TopCrtCorrection=fTopCrtTransformations.EastCC;
      else if(trackType==1) TopCrtCorrection=fTopCrtTransformations.EE;
      else if(trackType==2) TopCrtCorrection=fTopCrtTransformations.EW;
      else if(trackType==3) TopCrtCorrection=fTopCrtTransformations.WestCC;
      else if(trackType==4) TopCrtCorrection=fTopCrtTransformations.WE;
      else if(trackType==5) TopCrtCorrection=fTopCrtTransformations.WW;
      for(size_t j=1; j<hx.size();j++) deltaP+= sqrt(pow(hx[j]-hx[j-1],2) + pow(hy[j]-hy[j-1],2) + pow(hz[j]-hz[j-1],2));
      
      trackRes= track.Length()/goodHits;
      deltaP = deltaP/(goodHits-1);
      if(trackRes>10) std::cout<<"Track Res very bad "<<trackRes<<" "<<deltaP<<std::endl;
      icarus::crt::Direction trackPCADir={-5,-5,-5,0,0,0};

      if(goodHits<fMinimumGoodHits) {
        mf::LogDebug("CRTT0Tagging:")<<"Track does not have the minimal requirements of good hits: "<<goodHits;
        continue;
      }
      trackPCADir=fMatchingAlg.PCAfit(hx, hy, hz);
      std::vector<icarus::crt::CandCRT> crtCands;
      for(art::Ptr<CRTHit> p_crthit: CRTHitList){
        const CRTHit &crtHit = *p_crthit;
        double crtTime=crtHit.ts1_ns/1e3;
        icarus::crt::DriftedTrack thisDriftedTrack = fMatchingAlg.DriftTrack(trkHits, trkHitMetas, fGeometryService, detProp, crtTime, cryo, track);    
        if(thisDriftedTrack.outbound>0) continue;
        icarus::crt::Direction driftedTrackDir=fMatchingAlg.PCAfit(thisDriftedTrack.spx, thisDriftedTrack.spy, thisDriftedTrack.spz);
        int crtSys=-1;
        if(crtHit.plane<=34) crtSys=0;
        else if (crtHit.plane==50) crtSys=2;
        else crtSys=1;
        if(crtSys==2) continue; // lets discard Bottom CRT Hits for the moment

        double deltaX=std::numeric_limits<float>::signaling_NaN();
        double deltaZ=std::numeric_limits<float>::signaling_NaN();
        double crtDistance=std::numeric_limits<float>::signaling_NaN();

        icarus::crt::CrtPlane thisCrtPlane = fMatchingAlg.DeterminePlane(crtHit);
        icarus::crt::CrossPoint crossPoint = fMatchingAlg.DetermineProjection(driftedTrackDir, thisCrtPlane);
        double crtX=crtHit.x_pos;
        double crtY=crtHit.y_pos;
        double crtZ=crtHit.z_pos;

        if(crtSys==0){
          double centerDX, centerDY, centerDZ;
          centerDX=crtX-fTopCRTCenterMap[(int)crtHit.feb_id[0]].X;
          centerDY=crtY-fTopCRTCenterMap[(int)crtHit.feb_id[0]].Y;
          centerDZ=crtZ-fTopCRTCenterMap[(int)crtHit.feb_id[0]].Z;
          icarus::crt::AffineTrans thisAffine=TopCrtCorrection[(int)crtHit.feb_id[0]];
          std::pair<double,double> transCrt;
          if(crtHit.plane==30) {
            transCrt=icarus::crt::AffineTransformation(centerDX, centerDZ, thisAffine);
            crtX=transCrt.first;
            crtZ=transCrt.second;
          } else if(crtHit.plane==31 ||crtHit.plane==32) {
            transCrt=icarus::crt::AffineTransformation(centerDY, centerDZ, thisAffine);
            crtY=transCrt.first;
            crtZ=transCrt.second;            
          } else if(crtHit.plane==33 ||crtHit.plane==34){
            transCrt=icarus::crt::AffineTransformation(centerDX, centerDY, thisAffine);
            crtX=transCrt.first;
            crtY=transCrt.second;            
          }
          
        }
        if(thisCrtPlane.first==0){
          deltaX=crtX-crossPoint.X;
          deltaZ=crtZ-crossPoint.Z;
        } else if(thisCrtPlane.first==1){
          deltaX=crtY-crossPoint.Y;
          deltaZ=crtZ-crossPoint.Z;
        } else if(thisCrtPlane.first==2){
          deltaX=crtX-crossPoint.X;
          deltaZ=crtY-crossPoint.Y;          
        }
        crtDistance=sqrt(pow(deltaX,2)+pow(deltaZ,2));
        if(crtDistance>fMaximalCRTDistance) continue;
        icarus::crt::CandCRT thisCrtCand={crtHit,p_crthit, crtDistance, deltaX, deltaZ};
        crtCands.push_back(thisCrtCand);
      } // End of CRT Hit loop
      if(crtCands.empty()) {
        mf::LogDebug("CRTT0Tagging:")<<"No Good CRT match candidates for this track";
        continue;
      }
      auto minElementIt = std::min_element(crtCands.begin(), crtCands.end(), [](const icarus::crt::CandCRT& a, const icarus::crt::CandCRT& b) {
        return a.distance < b.distance;
      });
      icarus::crt::CandCRT& bestCrtCand=*minElementIt;
      if(bestCrtCand.distance<=fGoodCandidateDistance){ 
        int matchedSys=-1;
        if(bestCrtCand.CRThit.plane<=34) matchedSys=0;
        else if (bestCrtCand.CRThit.plane==50) matchedSys=2;
        else matchedSys=1;
        if(matchedSys==2) continue; // lets discard Bottom CRT Hits for the moment
        mf::LogInfo("CRTT0Tagging")
	      <<"Matched CRT time = "<<bestCrtCand.CRThit.ts1_ns/1e3<<" [us] to track "<<track.ID()<<" with projection-hit distance = "<<bestCrtCand.distance;        
        t0col->push_back(anab::T0(bestCrtCand.CRThit.ts1_ns, track.ID(), matchedSys, bestCrtCand.CRThit.plane,bestCrtCand.distance));
        util::CreateAssn(*this, e, *t0col, trkPtr, *trackAssn);
        util::CreateAssn(*this, e, *t0col, bestCrtCand.ptrCRThit, *t0CrtHitAssn);

      }
	 
	  } // End of Track Loop
	  
	} // End of Cryo Loop
  e.put(std::move(t0col));
  e.put(std::move(trackAssn));
  e.put(std::move(t0CrtHitAssn));
}

void CRTT0Tagging::endJob()
{

} // CRTT0Tagging::endJob()


DEFINE_ART_MODULE(CRTT0Tagging)
