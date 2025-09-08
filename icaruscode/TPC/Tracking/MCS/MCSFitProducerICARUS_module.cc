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

#include "icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "sbnobj/Common/Reco/RangeP.h"

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
   * Producer for TrajectoryMCSFitterICARUS, which performs a Maximum Likelihood fit of Multiple Coulomb Scattering angles between segments within a Track or Trajectory. It reads a recob::Track collection and produces a collection of recob::MCSFitResultGS where the elements are in the same order as the input collection (no explicit association is written).
   *
   * @author  G. Chiello (Pisa, ICARUS) and F. Varanini (Padova, ICARUS) based on code from G. Cerati 
   * @date    2025
   * @version 2.0
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
      fhicl::Atom<art::InputTag> rangeLabel {
        Name("rangeLabel"),
        Comment("Label of range")
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

    void produce(
      art::Event & e) override;

    std::vector<recob::Hit> projectHitsOnPlane(
      art::Event & e,
      const recob::Track& traj,
      int idx,
      unsigned int p, 
      std::vector<proxy::TrackPointData>& pdata) const;
    
    bool GeoStopCheck(
      const recob::Track& traj) const;

    void SelectEvent(
      int, int);

  private:
    Parameters p_;
    art::InputTag inputTag;
    art::InputTag thmTag;
    art::InputTag crtT0tag;
    art::InputTag rangeTag;
    TrajectoryMCSFitterICARUS mcsfitter;
    int selectCount = -1;
    int selectCryo = -1;
  };
}

trkf::MCSFitProducerICARUS::MCSFitProducerICARUS(
  trkf::MCSFitProducerICARUS::Parameters const & p) 
  : EDProducer{p}, 
  p_(p), 
  mcsfitter(p_().fitter) 
{
  inputTag = art::InputTag(p_().inputs().inputLabel());
  thmTag = art::InputTag(p_().inputs().thmLabel());
  crtT0tag = art::InputTag(p_().inputs().crtT0Label());
  rangeTag = art::InputTag(p_().inputs().rangeLabel());
  produces<std::vector<recob::MCSFitResultGS> >();
}

trkf::MCSFitProducerICARUS::~MCSFitProducerICARUS() {}

void trkf::MCSFitProducerICARUS::produce(art::Event & e) {
  auto output = std::make_unique<std::vector<recob::MCSFitResultGS> >();
  art::Handle<std::vector<recob::Track> > inputH;
  bool ok = e.getByLabel(inputTag, inputH);
  if (!ok) throw cet::exception("MCSFitProducerICARUS") << "Cannot find input art::Handle with inputTag " << inputTag;
  const auto& inputVec = *(inputH.product());
  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(inputTag); 
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e, thmTag);
  float minLen = 100;
  std::vector<int> count;
  count.push_back(-1); count.push_back(-1);

  std::vector<art::Ptr<recob::Track>> allPtrs;
  art::fill_ptr_vector(allPtrs, inputH);
  art::FindManyP<anab::T0> fmCRTTaggedT0(allPtrs, e, crtT0tag);

  auto const& range_handle = e.getValidHandle<std::vector<sbn::RangeP>>(rangeTag);
  std::vector<art::Ptr<sbn::RangeP>> rangePs;
  art::fill_ptr_vector(rangePs, range_handle);

  for (unsigned int je = 0; je < inputVec.size(); je++) {
    std::cout << " " << std::endl;
    auto element = inputVec.at(je);
    auto x = element.LocationAtPoint(0).X();
    int cryo = -1;
    if (x > 0) cryo = 1;
    else cryo = 0;
    if (cryo == 0) std::cout << "starting gsasso algorithm, event " << e.event() << " track id = " << je << " in cryostat EAST" << std::endl;
    if (cryo == 1) std::cout << "starting gsasso algorithm, event " << e.event() << " track id = " << je << " in cryostat WEST" << std::endl;
    std::cout << " " << std::endl;

    count[cryo]++;
    if ((count[cryo] == selectCount && cryo == selectCryo) || selectCount < 0) {//0 EAST 1 WEST
      auto const& start = element.Start(); 
      auto const& end = element.End(); 
      auto const& length = element.Length(); 
      std::cout << "first point [cm] = "
            << start.X() << ", "
            << start.Y() << ", "
            << start.Z() << std::endl;
      std::cout << "last point [cm] = "
            << end.X() << ", "
            << end.Y() << ", "
            << end.Z() << std::endl;
      
      std::vector<recob::Hit> hits2dC; hits2dC.clear();
      std::vector<recob::Hit> hits2dI2; hits2dI2.clear();
      std::vector<recob::Hit> hits2dI1; hits2dI1.clear();
      std::vector<proxy::TrackPointData> pdata; pdata.clear();
      
      bool geocheck = GeoStopCheck(element);
      if (!geocheck) {
        std::cout << "track not stopping, go to next track" << std::endl; 
        recob::MCSFitResultGS result = recob::MCSFitResultGS();
        output->emplace_back(std::move(result)); 
        continue; }
      if (length < minLen) {
        std::cout << "track too short, go to next track" << std::endl; 
        recob::MCSFitResultGS result = recob::MCSFitResultGS();
        output->emplace_back(std::move(result)); 
        continue; }

      hits2dI1 = projectHitsOnPlane(e, element, count[cryo], 0, pdata);
      hits2dI2 = projectHitsOnPlane(e, element, count[cryo], 1, pdata);
      hits2dC = projectHitsOnPlane(e, element, count[cryo], 2, pdata);
      std::vector<art::Ptr<recob::Hit>> emptyHitVector;
      art::Ptr<recob::Track> trkPtr = allPtrs.at(je);
      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);

      float CRTT0 = std::numeric_limits<float>::signaling_NaN();
      if (fmCRTTaggedT0.isValid()) {
        for (unsigned i_t0 = 0; i_t0 < fmCRTTaggedT0.size(); i_t0++) {
          if (fmCRTTaggedT0.at(trkPtr.key()).size()) {
            if (fmCRTTaggedT0.at(trkPtr.key()).at(0)->TriggerBits() != 0) continue;
            if (fmCRTTaggedT0.at(trkPtr.key()).at(0)->TriggerConfidence() > 100.) continue;
            CRTT0 = fmCRTTaggedT0.at(trkPtr.key()).at(0)->Time();
          }
        }
      }

      art::Ptr<sbn::RangeP> rangeP = rangePs.at(trkPtr.key());
      std::cout << "track length [cm] = " << length << std::endl;
      std::cout << "range momentum [GeV/c] = " << rangeP->range_p << std::endl;

      float CRTshift = 0;
      double vDrift = detProp.DriftVelocity();
      float fTickAtAnode = 850.;
      float fTickPeriod = 0.4;
      if (!std::isnan(CRTT0)) CRTshift = (CRTT0 / 1000. / fTickPeriod - fTickAtAnode) * fTickPeriod * vDrift;
      std::cout << "CRTT0 = " << CRTT0 << std::endl; std::cout << " " << std::endl;

      mcsfitter.set2DHitsI1(hits2dI1);
      mcsfitter.set2DHitsI2(hits2dI2);
      mcsfitter.set2DHitsC(hits2dC);
      mcsfitter.setPointData(pdata);
      mcsfitter.ComputeD3P(0);
      mcsfitter.ComputeD3P(1);
      mcsfitter.ComputeD3P(2);
      mcsfitter.setCRTShift(CRTshift);
      mcsfitter.setRangeP(rangeP->range_p);
      
      recob::MCSFitResultGS result = mcsfitter.fitMcs(element);
      output->emplace_back(std::move(result)); 
    }
  }
  e.put(std::move(output));
}

std::vector<recob::Hit> trkf::MCSFitProducerICARUS::projectHitsOnPlane(
  art::Event & e, 
  const recob::Track& traj, 
  int idx, 
  unsigned int p, 
  std::vector<proxy::TrackPointData>& pdata
) const {
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

bool trkf::MCSFitProducerICARUS::GeoStopCheck(const recob::Track& traj) const {
  size_t lastIndex = traj.LastValidPoint();
  geo::Point_t lastPoint = traj.LocationAtPoint(lastIndex); 
  double step = 20;

  double low_x; double high_x; 
  double last_x = lastPoint.X();
  if (last_x > 0) { low_x = 61.7; high_x = 358.73; }
  if (last_x < 0) { low_x = -358.73; high_x = -61.7; }
  //std::cout << "last point x = " << last_x << " low x = " << low_x << " high x = " << high_x << std::endl;
  bool check_x = true; 
  if ((abs(last_x - low_x) < step) || (abs(last_x - high_x) < step)) {
    check_x = false;
    std::cout << "WARNING! too near to borders in x direction!" << std::endl; }

  double low_y = -181.86; double high_y = 134.36; 
  double last_y = lastPoint.Y();
  //std::cout << "last point y = " << last_y << " low y = " << low_y << " high y = " << high_y << std::endl;
  bool check_y = true; 
  if ((abs(last_y - low_y) < step) || (abs(last_y - high_y) < step)) {
    check_y = false;
    std::cout << "WARNING! too near to borders in y direction!" << std::endl; }
  
  double low_z = -894.951; double high_z = 894.951; 
  double last_z = lastPoint.Z();
  //std::cout << "last point z = " << last_z << " low z = " << low_z << " high z = " << high_z << std::endl;
  bool check_z = true; 
  if ((abs(last_z - low_z) < step) || (abs(last_z - high_z) < step)) {
    check_z = false;
    std::cout << "WARNING! too near to borders in z direction!" << std::endl; }

  if (check_x && check_y && check_z) {
    std::cout << "stopping track found inside MCSFitProducer!" << std::endl; std::cout << " " << std::endl;
    return true; }
  else {
    std::cout << "crossing track found inside MCSFitProducer!" << std::endl; std::cout << " " << std::endl;
    return false; }
}

void trkf::MCSFitProducerICARUS::SelectEvent(int count, int cryo) {
  // Function to select event based on count and cryostat
  // Implementation here...
  selectCount = count;
  selectCryo = cryo;
  //std::cout << "Selected event with count: " << count << " and cryostat: " << cryo << std::endl;
}

DEFINE_ART_MODULE(trkf::MCSFitProducerICARUS)
