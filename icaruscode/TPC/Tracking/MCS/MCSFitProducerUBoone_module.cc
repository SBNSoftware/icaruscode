#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/ProductRetriever.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "TrajectoryMCSFitterUBoone.h"
#include "lardata/RecoBaseProxy/Track.h" //needed only if you do use the proxies

#include <memory>

namespace trkf {
  /**
   * @file  MCSFitProducer_module.cc
   * @class trkf::MCSFitProducer
   *
   * @brief Producer for TrajectoryMCSFitter.
   *
   * Producer for TrajectoryMCSFitter, which performs a Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory.
   * It reads a recob::Track collection and produces a collection of recob::MCSFitResult where the elements are in the same order as the input collection (no explicit association is written).
   *
   * For configuration options see MCSFitProducer#Inputs and MCSFitProducer#Config
   *
   * @author  G. Cerati (FNAL, MicroBooNE), based on code from L. Kalousis and D. Kaleko
   * @date    2017
   * @version 1.0
   */

  class MCSFitProducer : public art::EDProducer {
  public:
    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputLabel {
        Name("inputLabel"),
        Comment("Label of recob::TrackTrajectory Collection to be fit")};
    };
    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<MCSFitProducer::Inputs> inputs{
        Name("inputs"),
      };
      fhicl::Table<TrajectoryMCSFitterUBoone::Config> fitter{Name("fitter")};
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit MCSFitProducer(Parameters const& p);
    ~MCSFitProducer();

    // Plugins should not be copied or assigned.
    MCSFitProducer(MCSFitProducer const&) = delete;
    MCSFitProducer(MCSFitProducer&&) = delete;
    MCSFitProducer& operator=(MCSFitProducer const&) = delete;
    MCSFitProducer& operator=(MCSFitProducer&&) = delete;
    std::vector<recob::Hit> projectHitsOnPlane(art::Event & e, const recob::Track& traj, int idx, unsigned int p, std::vector<proxy::TrackPointData>& pdata) const;

  private:
    void produce(art::Event& e) override;

    Parameters p_;
    art::InputTag inputTag;
    art::InputTag PFPTag;
    TrajectoryMCSFitterUBoone mcsfitter;
  };
}

trkf::MCSFitProducer::MCSFitProducer(trkf::MCSFitProducer::Parameters const& p) : EDProducer{p}, p_(p), mcsfitter(p_().fitter) {
  inputTag = art::InputTag(p_().inputs().inputLabel());
  produces<std::vector<recob::MCSFitResult>>();
}

trkf::MCSFitProducer::~MCSFitProducer() {}

void trkf::MCSFitProducer::produce(art::Event& e) {
  auto output = std::make_unique<std::vector<recob::MCSFitResult>>();
  art::Handle<std::vector<recob::Track>> inputH;
  bool ok = e.getByLabel(inputTag, inputH);
  if (!ok) throw cet::exception("MCSFitProducer") << "Cannot find input art::Handle with inputTag " << inputTag;
  auto inputVec = *(inputH.product());
  float minLen = 100;
  std::vector<int> count;
  count.push_back(-1); count.push_back(-1);
  for (size_t trackIdx = 0; trackIdx < inputVec.size(); trackIdx++) {
    std::vector<float> dum;
    recob::MCSFitResult result = recob::MCSFitResult(0, 0, 0, 0, 0, 0, 0, dum, dum);;
    auto element = inputVec.at(trackIdx);
    std::vector<recob::Hit> hits2dC; hits2dC.clear();
    std::vector<recob::Hit> hits2dI2; hits2dI2.clear();
    std::vector<recob::Hit> hits2dI1; hits2dI1.clear();
    std::vector<proxy::TrackPointData> pdata; pdata.clear();
    auto x = element.LocationAtPoint(0).X();
    int cryo = -1;
    if (x > 0) cryo = 1;
    else cryo = 0;
    std::cout << " " << std::endl;
    if (cryo == 0) std::cout << "starting uboone algorithm, track with id = " << trackIdx << " in cryostat EAST" << std::endl;
    if (cryo == 1) std::cout << "starting uboone algorithm, track with id = " << trackIdx << " in cryostat WEST" << std::endl;
    std::cout << " " << std::endl;
    count[cryo]++;
    if (element.Length() > minLen) {
      hits2dI1 = projectHitsOnPlane(e, element, count[cryo], 0, pdata);
      hits2dI2 = projectHitsOnPlane(e, element, count[cryo], 1, pdata);
      hits2dC = projectHitsOnPlane(e, element, count[cryo], 2, pdata);
    }
    mcsfitter.set2DHitsI1(hits2dI1);
    mcsfitter.set2DHitsI2(hits2dI2);
    mcsfitter.set2DHitsC(hits2dC);
    mcsfitter.setPointData(pdata);
    if (element.Length() > minLen) result = mcsfitter.fitMcs(element);
    output->emplace_back(std::move(result));
  }
  e.put(std::move(output));
}

std::vector<recob::Hit> trkf::MCSFitProducer::projectHitsOnPlane(art::Event & e, const recob::Track& traj, int idx, unsigned int p, std::vector<proxy::TrackPointData>& pdata) const {
  std::vector<recob::Hit> v;
  auto const& tracks = proxy::getCollection<proxy::Tracks>(e, inputTag);
  const auto& track = tracks[idx];
  for (const art::Ptr<recob::Hit>& h : track.hits()) {
    if (h->WireID().Plane == p) v.emplace_back(*h);
  }
  for (unsigned int jp; jp < traj.NPoints(); jp++) {
    auto pd = proxy::makeTrackPointData(track, jp);
    pdata.emplace_back(pd);
  }
  return v;
}

DEFINE_ART_MODULE(trkf::MCSFitProducer)