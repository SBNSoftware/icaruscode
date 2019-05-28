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
std::vector<recob::Hit> projectHitsOnPlane(art::Event & e,const recob::Track& traj,unsigned int p) const;

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

for (const auto& element : inputVec) {
    //fit
    std::vector<recob::Hit> hits2d=projectHitsOnPlane(e,element,2);
    mcsfitter.set2DHits(hits2d);
    mcsfitter.ComputeD3P();
    recob::MCSFitResult result = mcsfitter.fitMcs(element);
    output->emplace_back(std::move(result));
  }

  e.put(std::move(output));
}

std::vector<recob::Hit> trkf::MCSFitProducerICARUS::projectHitsOnPlane(art::Event & e,const recob::Track& traj,unsigned int p) const
{
std::vector<recob::Hit> v;
  // Get track collection proxy and parallel mcs fit data (associated hits loaded by default)
  // Note: if tracks were produced from a TrackTrajectory collection you could access the original trajectories adding ',proxy::withOriginalTrajectory()' to the list of arguments
std::cout << " before calling proxy " << std::endl;
auto const& tracks   = proxy::getCollection<proxy::Tracks>(e,inputTag);
std::cout << " after calling proxy " << std::endl;
const auto& track = tracks[0];
std::cout << " proxy nhits " << track.nHits() << std::endl;
 
        for (const art::Ptr<recob::Hit>& h : track.hits())
         if(h->WireID().Plane == p) {
          std::cout << "collection hit wire=" << h->WireID() << " peak time=" << h->PeakTime() << std::endl;
          v.emplace_back(*h);
        }
return v;
}


DEFINE_ART_MODULE(trkf::MCSFitProducerICARUS)
