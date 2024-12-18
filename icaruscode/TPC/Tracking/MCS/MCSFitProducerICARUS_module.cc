#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/T0.h"


//#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "icaruscode/TPC/Tracking/MCS/TrajectoryMCSFitterICARUS.h"
#include "lardata/RecoBaseProxy/Track.h" //needed only if you do use the proxies
#include <memory>

namespace trkf {
  /**
   * @file  MCSFitProducerICARUS_module.cc
   * @class trkf::MCSFitProducerICARUS
   *
   * @brief Producer for TrajectoryMCSFitterICARUS.
   *
   * Producer for TrajectoryMCSFitterICARUS, which performs a Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory. 
   * It reads a recob::Track collection and produces a collection of recob::MCSFitResult where the elements are in the same order as the input collection (no explicit association is written).
   *
   * For configuration options see MCSFitProducer#Inputs and MCSFsitProducer#Config
   *
   * @author  F. Varanini (Padova, ICARUS)
   * @date    2018
   * @version 1.0
   */
  class MCSFitProducerICARUS : public art::EDProducer {
  public:

    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputLabel {
	Name("inputLabel"),
	Comment("Label of recob::TrackTrajectory Collection to be fit")
      };
           fhicl::Atom<art::InputTag> thmLabel {
	Name("ThmLabel"),
	Comment("Label of TrackHitMetas")
      };   
             fhicl::Atom<art::InputTag> crtT0Label {
	Name("CRTT0Label"),
	Comment("Label of CRTT0")
      };
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<MCSFitProducerICARUS::Inputs> inputs {
	Name("inputs"),
      };
      fhicl::Table<TrajectoryMCSFitterICARUS::Config> fitter {
	Name("fitter")
      };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit MCSFitProducerICARUS(Parameters const & p);
    ~MCSFitProducerICARUS();

    // Plugins should not be copied or assigned.
    MCSFitProducerICARUS(MCSFitProducerICARUS const &) = delete;
    MCSFitProducerICARUS(MCSFitProducerICARUS &&) = delete;
    MCSFitProducerICARUS & operator = (MCSFitProducerICARUS const &) = delete;
    MCSFitProducerICARUS & operator = (MCSFitProducerICARUS &&) = delete;

    void produce(art::Event & e) override;
std::vector<recob::Hit> projectHitsOnPlane(art::Event & e,const recob::Track& traj,int idx,unsigned int p, std::vector<proxy::TrackPointData>& pdata) const;
//float CRTT0Shift(const std::vector<recob::Hit>& trkHits, const std::vector<const recob::TrackHitMeta*>& trkHitMetas, const geo::GeometryCore *GeometryService, detinfo::DetectorPropertiesData const& detProp, double time, const recob::Track& tpcTrack);

  private:
    Parameters p_;
    art::InputTag inputTag;
    art::InputTag thmTag;
     art::InputTag crtT0tag;
    TrajectoryMCSFitterICARUS mcsfitter;
  };
}

trkf::MCSFitProducerICARUS::MCSFitProducerICARUS(trkf::MCSFitProducerICARUS::Parameters const & p)
  : EDProducer{p}, p_(p), mcsfitter(p_().fitter)
{
  inputTag = art::InputTag(p_().inputs().inputLabel());
    thmTag = art::InputTag(p_().inputs().thmLabel());
    crtT0tag = art::InputTag(p_().inputs().crtT0Label());

  produces<std::vector<recob::MCSFitResult> >();
}

trkf::MCSFitProducerICARUS::~MCSFitProducerICARUS() {}

void trkf::MCSFitProducerICARUS::produce(art::Event & e)
{
  std::cout << " MCSFitProducerICARUS produce " << std::endl;
  //
  auto output  = std::make_unique<std::vector<recob::MCSFitResult> >();
  //
  art::Handle<std::vector<recob::Track> > inputH;
  bool ok = e.getByLabel(inputTag,inputH);
  if (!ok) throw cet::exception("MCSFitProducerICARUS") << "Cannot find input art::Handle with inputTag " << inputTag;
  const auto& inputVec = *(inputH.product());



// Track - associated data
//art::FindManyP<recob::Track> fmTracks(PFParticleList, e, fTPCTrackLabel[it]);
//
// Collect all hits
//art::ValidHandle<std::vector<recob::Hit>> allhit_handle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel[it]);
//std::vector<art::Ptr<recob::Hit>> allHits;
//art::fill_ptr_vector(allHits, allhit_handle);

//for (art::Ptr<recob::PFParticle> p_pfp: PFParticleList) {
//      const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(p_pfp.key());
//      if (thisTrack.size() != 1) continue;    
//      art::Ptr<recob::Track> trkPtr = thisTrack.at(0);
//      const recob::Track &track = *trkPtr;
  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(inputTag); 

 art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e,thmTag);
      

std::cout << " inputh size " << inputVec.size() << std::endl;

float minLen=40;
//float cutFinLen=40;
std::vector<int> count;
count.push_back(-1);count.push_back(-1);

 std::vector<art::Ptr<recob::Track>> allPtrs;
    art::fill_ptr_vector(allPtrs, inputH);
    std::cout << " allptrs size " << allPtrs.size() << std::endl;



    //t0
art::FindManyP<anab::T0> fmCRTTaggedT0(allPtrs,e,crtT0tag);
if(fmCRTTaggedT0.isValid()) {
std::cout << " fmCRTT0 size" << fmCRTTaggedT0.size() << std::endl;
std::cout << " fmCRTT0 0" << fmCRTTaggedT0.at(0).size() << std::endl;
}
for (unsigned int je=0;je<inputVec.size();je++) {
  auto element=inputVec.at(je);
 std::vector<float> dum;
  recob::MCSFitResult result=recob::MCSFitResult(0,
			    0,0,0,
0,0,0,
			     dum,dum);
 
//cut final length to simulate non-contained

    //fit
    std::vector<recob::Hit> hits2dC;
    hits2dC.clear();
      std::vector<recob::Hit> hits2dI2;
    hits2dI2.clear();
      std::vector<recob::Hit> hits2dI1;
    hits2dI1.clear();
    std::vector<proxy::TrackPointData> pdata;
    pdata.clear();
   
auto x=element.LocationAtPoint(0).X();
    int cryo=-1;
    if(x>0) cryo=1;
    else cryo=0;
    count[cryo]++;
     //if(element.Length()>minLen) std::cout << " track " << count << " length " << element.Length() << std::endl;

   // std::cout << " x " << x << " cryo " << cryo << std::endl;

   // if(count[cryo]==4&&cryo==1) { //0 EAST 1 WEST

       if(element.Length()>minLen) {
       hits2dC=projectHitsOnPlane(e,element,count[cryo],2,pdata);
       hits2dI2=projectHitsOnPlane(e,element,count[cryo],1,pdata);
       hits2dI1=projectHitsOnPlane(e,element,count[cryo],0,pdata);
       }
     std::vector<art::Ptr<recob::Hit>> emptyHitVector;
     art::Ptr<recob::Track> trkPtr=allPtrs.at(je);
   //   const std::vector<art::Ptr<recob::Hit>> &trkHits  = fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;
   //  std::vector<const recob::TrackHitMeta*> emptyTHMVector;
   //const std::vector<const recob::TrackHitMeta*> &trkHitMetas = fmtrkHits.isValid() ? fmtrkHits.data(trkPtr.key()) : emptyTHMVector;
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

 //geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);

float CRTT0=std::numeric_limits<float>::signaling_NaN();
if(fmCRTTaggedT0.isValid()) {
for(unsigned i_t0=0;i_t0<fmCRTTaggedT0.size();i_t0++) {

//std::cout << " key "<< trkPtr.key() << " fmcrttaggedt0 size " << fmCRTTaggedT0.size() << " fmcrttaggedt0[0] size " << fmCRTTaggedT0[i_t0].size() << std::endl;

	if(fmCRTTaggedT0.at(trkPtr.key()).size()) {
         std::cout << " triggerbits " << fmCRTTaggedT0.at(trkPtr.key()).at(0)->TriggerBits() << " triggerconfidence " << fmCRTTaggedT0.at(trkPtr.key()).at(0)->TriggerConfidence() << std::endl;
            if(fmCRTTaggedT0.at(trkPtr.key()).at(0)->TriggerBits() != 0) continue;
        if(fmCRTTaggedT0.at(trkPtr.key()).at(0)->TriggerConfidence() > 70.) continue;
           CRTT0=fmCRTTaggedT0.at(trkPtr.key()).at(0)->Time();
  }}}
  float CRTshift=0;
   double vDrift=detProp.DriftVelocity();
   float fTickAtAnode=850.;
   float fTickPeriod=0.4;
   std::cout << " crt ratio " << CRTT0/fTickPeriod/1000. << std::endl;
   std::cout << " crt tickshift " << CRTT0/1000./fTickPeriod-fTickAtAnode << std::endl;

  if(!std::isnan(CRTT0)) 
CRTshift=(CRTT0/1000./fTickPeriod-fTickAtAnode)*fTickPeriod*vDrift;
    
  std::cout << " CRTT0 " << CRTT0 << std::endl;
    std::cout << " CRTshift " << CRTshift << std::endl;
   // float shift=CRTT0Shift(hits2dC, trkHitMetas, fGeometryService, detProp, time, element);

    mcsfitter.set2DHitsC(hits2dC);
    mcsfitter.set2DHitsI2(hits2dI2);
    mcsfitter.set2DHitsI1(hits2dI1);
    mcsfitter.setPointData(pdata);

    mcsfitter.ComputeD3P(0);
    mcsfitter.ComputeD3P(1);
    mcsfitter.ComputeD3P(2);

    mcsfitter.setCRTShift(CRTshift);
    
 
    std::cout << " fitting icarus trackIdx " << count[cryo] << " cryo " << cryo << " length " << element.Length() << std::endl;
    std::cout << " 3dpoints " << element.NPoints() << " coll hits " << hits2dC.size() << " length " << element.Length() << std::endl;
    std::cout << " average 3d pitch " << element.Length()/element.NPoints() << " average coll pitch " << element.Length()/hits2dC.size() << std::endl;

    if(element.Length()>minLen) result = mcsfitter.fitMcs(element);
    if(result.fwdMomentum()>0.) std::cout << " fitting icarus trackIdx " << count[cryo] << " cryo " << cryo << " mcs momentum " << result.fwdMomentum() << std::endl;

    
output->emplace_back(std::move(result));
 // }

  }

  e.put(std::move(output));

}

std::vector<recob::Hit> trkf::MCSFitProducerICARUS::projectHitsOnPlane(art::Event & e,const recob::Track& traj, int idx, unsigned int p,std::vector<proxy::TrackPointData>& pdata) const
{
std::vector<recob::Hit> v;
  // Get track collection proxy and parallel mcs fit data (associated hits loaded by default)
  // Note: if tracks were produced from a TrackTrajectory collection you could access the original trajectories adding ',proxy::withOriginalTrajectory()' to the list of arguments

auto const& tracks   = proxy::getCollection<proxy::Tracks>(e,inputTag);
std::cout << " tracks size in projecting " << tracks.size() << std::endl;
const auto& track = tracks[idx];
//auto hap=track.hitAtPoint(0);
        for (const art::Ptr<recob::Hit>& h : track.hits()) {
         if(h->WireID().Plane == p) {
          std::cout << v.size() << "th hit wire=" << h->WireID() << " peak time=" << h->PeakTime() << std::endl;
          v.emplace_back(*h);
        }
        }
        std::cout << " making trackpointdata npoints " << traj.NPoints() << std::endl;
               for (unsigned int jp; jp<traj.NPoints();jp++) {
auto pd=proxy::makeTrackPointData(track,jp);
         //std::cout << " emplacing pdata " << jp << std::endl;
          pdata.emplace_back(pd);
          
        }

return v;
}

DEFINE_ART_MODULE(trkf::MCSFitProducerICARUS)
