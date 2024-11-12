#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"

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

  private:
    Parameters p_;
    art::InputTag inputTag;
    TrajectoryMCSFitterICARUS mcsfitter;
  };
}

trkf::MCSFitProducerICARUS::MCSFitProducerICARUS(trkf::MCSFitProducerICARUS::Parameters const & p)
  : EDProducer{p}, p_(p), mcsfitter(p_().fitter)
{
  inputTag = art::InputTag(p_().inputs().inputLabel());
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

std::cout << " inputh size " << inputVec.size() << std::endl;

float minLen=40;
//float cutFinLen=40;
std::vector<int> count;
count.push_back(-1);count.push_back(-1);
for (const auto& element : inputVec) {
 std::vector<float> dum;
  recob::MCSFitResult result=recob::MCSFitResult(0,
			    0,0,0,
0,0,0,
			     dum,dum);
 
//cut final length to simulate non-contained

    //fit
    std::vector<recob::Hit> hits2d;
    hits2d.clear();
    std::vector<proxy::TrackPointData> pdata;
    pdata.clear();
   
 
  
auto x=element.LocationAtPoint(0).X();
    int cryo=-1;
    if(x>0) cryo=1;
    else cryo=0;
    count[cryo]++;
     //if(element.Length()>minLen) std::cout << " track " << count << " length " << element.Length() << std::endl;

   // std::cout << " x " << x << " cryo " << cryo << std::endl;

    //if(count[cryo]==1&&cryo==0) { //0 EAST 1 WEST

       if(element.Length()>minLen) hits2d=projectHitsOnPlane(e,element,count[cryo],2,pdata);

    mcsfitter.set2DHits(hits2d);
    mcsfitter.setPointData(pdata);

    mcsfitter.ComputeD3P();
 
    std::cout << " fitting icarus trackIdx " << count[cryo] << " cryo " << cryo << " length " << element.Length() << std::endl;
    std::cout << " 3dpoints " << element.NPoints() << " coll hits " << hits2d.size() << " length " << element.Length() << std::endl;
    std::cout << " average 3d pitch " << element.Length()/element.NPoints() << " average coll pitch " << element.Length()/hits2d.size() << std::endl;

    if(element.Length()>minLen) result = mcsfitter.fitMcs(element);
    if(result.fwdMomentum()>0.) std::cout << " fitting icarus trackIdx " << count[cryo] << " cryo " << cryo << " mcs momentum " << result.fwdMomentum() << std::endl;

    
output->emplace_back(std::move(result));
  //}

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
