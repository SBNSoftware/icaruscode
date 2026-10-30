#include "TrajectoryMCSFitterICARUS.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "icaruscode/TPC/Tracking/MCS/MCSFitResultGS.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "sbnobj/Common/CRT/CRTHitT0TaggingInfo.hh"
#include "sbnobj/Common/CRT/CRTHitT0TaggingTruthInfo.hh"

#include "sbnobj/Common/Reco/RangeP.h"

#include "sbncode/Calibration/ITCSSelectionTool.h"
#include "art/Utilities/make_tool.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "icaruscode/TPC/Tracking/MCS/TrajectoryMCSFitterICARUS.h"
#include "lardata/RecoBaseProxy/Track.h" 
#include <memory>

namespace trkf {
  class MCSFitProducerICARUS : public art::EDProducer {
  public:
    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> TRKLabel   { Name("TRKLabel") };
      fhicl::Atom<art::InputTag> TRKHMLabel { Name("TRKHMLabel") };
      fhicl::Atom<art::InputTag> PFPLabel   { Name("PFPLabel") };
      fhicl::Atom<art::InputTag> PFPT0Label { Name("PFPT0Label") };
      fhicl::Atom<art::InputTag> CRTT0Label { Name("CRTT0Label") };
      fhicl::Atom<art::InputTag> rangeLabel { Name("rangeLabel") };
      fhicl::Atom<art::InputTag> caloLabel  { Name("caloLabel") }; };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<MCSFitProducerICARUS::Inputs> inputs      { Name("inputs") };
      fhicl::Table<TrajectoryMCSFitterICARUS::Config> fitter { Name("fitter") }; };
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

  private:
    Parameters p_;
    art::InputTag TRKTag;
    art::InputTag TRKHMTag;
    art::InputTag PFPTag;
    art::InputTag PFPT0Tag;
    art::InputTag CRTT0Tag;
    art::InputTag rangeTag;
    art::InputTag caloTag;
    TrajectoryMCSFitterICARUS mcsfitter; }; }

trkf::MCSFitProducerICARUS::MCSFitProducerICARUS(trkf::MCSFitProducerICARUS::Parameters const & p) 
  : EDProducer{p}, 
    p_(p), 
    mcsfitter(p_().fitter) {
      TRKTag   = art::InputTag(p_().inputs().TRKLabel());
      TRKHMTag = art::InputTag(p_().inputs().TRKHMLabel());
      PFPTag   = art::InputTag(p_().inputs().PFPLabel());
      PFPT0Tag = art::InputTag(p_().inputs().PFPT0Label());
      CRTT0Tag = art::InputTag(p_().inputs().CRTT0Label());
      rangeTag = art::InputTag(p_().inputs().rangeLabel());
      caloTag  = art::InputTag(p_().inputs().caloLabel());

      produces<std::vector<recob::MCSFitResultGS>>(); }

trkf::MCSFitProducerICARUS::~MCSFitProducerICARUS() {}

void trkf::MCSFitProducerICARUS::produce(art::Event & e) {
  mcsfitter.selectionMode();
  mcsfitter.verbose();

  auto const& wr = art::ServiceHandle<geo::WireReadout const>()->Get();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);
  auto const* geom = lar::providerFrom<geo::Geometry>();

  mcsfitter.setWR(wr);
  mcsfitter.setDetProp(detProp);
  mcsfitter.setGeom(geom);

  // build geometry, TPC and Active Volumes (AVs)
  std::vector<std::vector<geo::BoxBoundedGeo>> TPCVols;
  std::vector<geo::BoxBoundedGeo> AVs;
  for (auto const &cryo: geom->Iterate<geo::CryostatGeo>()) {
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    for (auto const& TPC : geom->Iterate<geo::TPCGeo>(cryo.ID())) {
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox()); }
    TPCVols.push_back(std::move(this_tpc_volumes)); }
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: TPCVols) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();
    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();
    AVs.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax); }

  // define tracks and MCS results vector
  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(TRKTag); 
  std::vector<recob::MCSFitResultGS> results(tracks->size());

  // range momentum module
  auto const& range_handle = e.getValidHandle<std::vector<sbn::RangeP>>(rangeTag);
  std::vector<art::Ptr<sbn::RangeP>> rangePs;
  art::fill_ptr_vector(rangePs, range_handle);

  // PFP information and T0 module
  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(PFPTag);
  art::fill_ptr_vector(PFParticleList, pfparticles);
  art::FindManyP<anab::T0> fPFPT0(PFParticleList, e, PFPTag);
  art::FindManyP<recob::SpacePoint> PFParticleSPs(PFParticleList, e, PFPTag);

  // CRT based T0 module
  art::FindManyP<anab::T0> fCRTT0(tracks, e, CRTT0Tag);
  art::FindManyP<sbn::crt::CRTHitT0TaggingInfo> fmCRTHitT0TaggingInfo(tracks, e, CRTT0Tag);

  // track associated data
  art::FindManyP<recob::Track> fmTracks(pfparticles, e, TRKTag);
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e, TRKHMTag);
  art::FindManyP<anab::Calorimetry> fmCalo(tracks, e, caloTag);
  auto const& proxyTracks = proxy::getCollection<proxy::Tracks>(e, TRKTag);

  // loop over PFP just like skimmer
  for (art::Ptr<recob::PFParticle> pfp : PFParticleList) {
    if (pfp->PdgCode() == 11) continue;
    auto const& thisTrack = fmTracks.at(pfp.key());
    if (thisTrack.size() != 1) continue;
    art::Ptr<recob::Track> trkPtr = thisTrack.at(0);
    std::vector<art::Ptr<recob::Hit>> trkHits = fmtrkHits.at(trkPtr.key());
    std::vector<const recob::TrackHitMeta*> trkHitMetas = fmtrkHits.data(trkPtr.key());
    std::vector<art::Ptr<anab::Calorimetry>> calo = fmCalo.at(trkPtr.key());
    if (trkHits.size() != trkHitMetas.size()) continue;

    // identify cryostat and print track information
    geo::Point_t const& pos = trkPtr->LocationAtPoint(0);
    geo::Point_t point{pos.X(), pos.Y(), pos.Z()};
    std::size_t cryo = geom->PositionToCryostatID(point).Cryostat;
    const char* cryoName = (cryo == 0) ? "EAST" : "WEST";
    if (verbose_ > 0) std::cout << "producer MCS: event ID = " << e.event() 
              << " track ID = " << trkPtr->ID() 
              << " cryostat = " << cryoName
              << std::endl;

    // recall range momentum and print it with track length
    art::Ptr<sbn::RangeP> rangeP = rangePs.at(trkPtr.key());
    mcsfitter.setRangeP(rangeP->range_p);
    if (verbose_ > 0) std::cout << "range momentum [GeV/c] = " << rangeP->range_p 
              << " track length [cm] = " << trkPtr->Length() 
              << std::endl;

    // T0 from PFP associated data
    bool hasPFPT0 = false;
    double PFPT0 = std::numeric_limits<float>::signaling_NaN();
    if (fPFPT0.isValid()) {
      auto const& PFPT0s = fPFPT0.at(pfp.key());
      if (!PFPT0s.empty()) {
        PFPT0 = PFPT0s.at(0)->Time();
        hasPFPT0 = true;
        if (verbose_ > 0) std::cout << "track " << trkPtr->ID() << " has PFPT0 = " << PFPT0 << std::endl; } }
    mcsfitter.sethasPFPT0(hasPFPT0);
    mcsfitter.setPFPT0(PFPT0);

    // T0 from CRT Hit
    bool hasCRTT0 = false;
    double CRTT0 = std::numeric_limits<float>::signaling_NaN();
    if (fCRTT0.isValid()) {
      auto const& CRTT0s = fCRTT0.at(trkPtr.key());
      if (!CRTT0s.empty()) {
        const sbn::crt::CRTHitT0TaggingInfo &tag = *fmCRTHitT0TaggingInfo.at(trkPtr.key()).at(0);
        double time = CRTT0s.at(0)->Time();
        bool fIncludeTopCRT = true;
        bool fIncludeSideCRT = false;
        bool crtHitSysRejected = (tag.Sys == 0 && !fIncludeTopCRT) || (tag.Sys==1 && !fIncludeSideCRT); 
        geo::Point_t end {trkPtr->Start().X(), trkPtr->Start().Y(), trkPtr->Start().Z()};
        bool trackIsStopping = false;
        for (auto const &v: AVs) {
          if (v.ContainsPosition(end)) trackIsStopping = true; }
        bool crtHitDistanceRejected = trackIsStopping ? 
            ((tag.Sys == 0) ? (tag.Distance > 100.) : (tag.Distance > 100.)):
            ((tag.Sys == 0) ? (tag.Distance > 100.) : (tag.Distance > 100.));
        if (!crtHitSysRejected && !crtHitDistanceRejected) {
          CRTT0 = time;
          hasCRTT0 = true;
          if (verbose_ > 0) std::cout << "track " << trkPtr->ID() << " has CRTT0 = " << CRTT0 << std::endl; } } }
    mcsfitter.sethasCRTT0(hasCRTT0);
    mcsfitter.setCRTT0(CRTT0);

    // convert T0 from ns to µs
    bool hasT0 = false;
    double T0 = std::numeric_limits<float>::signaling_NaN();
    if (hasPFPT0) T0 = PFPT0 * 1e-3; 
    else if (hasCRTT0) T0 = CRTT0 * 1e-3;
    hasT0 = hasPFPT0 || hasCRTT0;
    mcsfitter.sethasT0(hasT0);
    mcsfitter.setT0(T0);

    // "whicht0" should reflect the T0 used for the reconstruction of the drift coordinate
    int whicht0 = -1;
    if (hasPFPT0) whicht0 = 0;
    else if (hasCRTT0) whicht0 = 2;
    if (verbose_ > 0) std::cout << "whichT0 = " << whicht0 << std::endl;
    
    // build track object for selection
    sbn::TrackInfo SBNTrack;
    if (hasT0 && !trkHits.empty()) {
      const double shift = geom->TPC(trkHits.back()->WireID()).DriftDir().X() * detProp.DriftVelocity() * T0;
      SBNTrack.end.x = trkPtr->End().X() + shift; }
    SBNTrack.end.y = trkPtr->End().Y();
    SBNTrack.end.z = trkPtr->End().Z();
    SBNTrack.dir.x = trkPtr->StartDirection().X();
    SBNTrack.dir.y = trkPtr->StartDirection().Y();
    SBNTrack.dir.z = trkPtr->StartDirection().Z(); 
    SBNTrack.whicht0 = whicht0;

    // fill hit information with calorimetry
    for (size_t i = 0; i < trkHits.size(); ++i) {
      const art::Ptr<recob::Hit>& hit = trkHits[i];
      const recob::TrackHitMeta* meta = trkHitMetas[i];
      if (!hit || !meta) continue;
      unsigned hindex = hit.key();
      sbn::TrackHitInfo hinfo; 
      hinfo.h.time  = hit->PeakTime();
      hinfo.h.tpc   = hit->WireID().TPC;
      hinfo.h.plane = hit->WireID().Plane;
      hinfo.h.id = (int)hindex;
      bool badhit = (meta->Index() == std::numeric_limits<unsigned int>::max()) ||
                    (!trkPtr->HasValidPoint(meta->Index()));
      if (!badhit) {
        for (const art::Ptr<anab::Calorimetry> &c: calo) {
          if (c->PlaneID().Plane != hinfo.h.plane) continue;
          for (unsigned i_calo = 0; i_calo < c->dQdx().size(); i_calo++) {
            if (c->TpIndices()[i_calo] == hindex) { 
              hinfo.oncalo = true;
              hinfo.pitch = c->TrkPitchVec()[i_calo];
              hinfo.dqdx = c->dQdx()[i_calo];
              hinfo.rr = c->ResidualRange()[i_calo];
              break; } }
          if (hinfo.oncalo) break; } }
      if (hinfo.h.plane == 2) {
        SBNTrack.hits2.push_back(hinfo); } }

    // compute min and max hit time
    double minTPCE = -1;
    double maxTPCE = -1;
    double minTPCW = -1;
    double maxTPCW = -1;
    bool inTPCE = false;
    for (const sbn::TrackHitInfo &h: SBNTrack.hits2) {
      inTPCE = (h.h.tpc <= 1);
      if (h.oncalo && inTPCE == true) {
        if (minTPCE < 0. || h.h.time < minTPCE) minTPCE = h.h.time;
        if (maxTPCE < 0. || h.h.time > maxTPCE) maxTPCE = h.h.time; } 
      if (h.oncalo && inTPCE == false) {
        if (minTPCW < 0. || h.h.time < minTPCW) minTPCW = h.h.time;
        if (maxTPCW < 0. || h.h.time > maxTPCW) maxTPCW = h.h.time; } }
    SBNTrack.hit_min_time_p2_tpcE = minTPCE;
    SBNTrack.hit_max_time_p2_tpcE = maxTPCE;
    SBNTrack.hit_min_time_p2_tpcW = minTPCW;
    SBNTrack.hit_max_time_p2_tpcW = maxTPCW;

    // check if track is downgoing
    bool RequireDownwards = true;
    bool downwards = (SBNTrack.dir.y < 0.) || !RequireDownwards;

    // check if track is within fiducial volume
    std::vector<geo::BoxBoundedGeo> FVs;
    double FVInsetMinX = 15; double FVInsetMaxX = 15;
    double FVInsetMinY = 25; double FVInsetMaxY = 15;
    double FVInsetMinZ = 15; double FVInsetMaxZ = 15;
    for (const geo::BoxBoundedGeo &AV: AVs) {
      FVs.emplace_back(AV.MinX() + FVInsetMinX, AV.MaxX() - FVInsetMaxX, 
                      AV.MinY() + FVInsetMinY, AV.MaxY() - FVInsetMaxY, 
                      AV.MinZ() + FVInsetMinZ, AV.MaxZ() - FVInsetMaxZ); }
    bool CheckFiducialX = true;
    bool end_is_fid = false;
    for (const geo::BoxBoundedGeo &g: FVs) {
      geo::Point_t end {SBNTrack.end.x, SBNTrack.end.y, SBNTrack.end.z};
      bool is_contained = CheckFiducialX ? g.ContainsPosition(end) : g.ContainsYZ(end.Y(), end.Z());
      if (is_contained) {
        end_is_fid = true;
        break; } }

    // check if collection plane times are fiducial
    unsigned NumberTimeSamples = 4096;
    double TickMin = 0;
    double TickMax = TickMin + NumberTimeSamples;
    double MinTimeTickInset = 100;
    double MaxTimeTickInset = 100; 
    double FidTickMin = TickMin + MinTimeTickInset;
    double FidTickMax = TickMax - MaxTimeTickInset;
    bool time_is_fid = (SBNTrack.hit_min_time_p2_tpcE < 0. || SBNTrack.hit_min_time_p2_tpcE > FidTickMin) &&
                      (SBNTrack.hit_max_time_p2_tpcE < 0. || SBNTrack.hit_max_time_p2_tpcE < FidTickMax) &&
                      (SBNTrack.hit_min_time_p2_tpcW < 0. || SBNTrack.hit_min_time_p2_tpcW > FidTickMin) &&
                      (SBNTrack.hit_max_time_p2_tpcW < 0. || SBNTrack.hit_max_time_p2_tpcW < FidTickMax);

    // check if median dQdx of the last 5 cm is greater than 1000 ADC/cm
    double MediandQdxRRMax = 5;
    std::vector<double> endp_dqdx;
    for (const sbn::TrackHitInfo &h: SBNTrack.hits2) {
      if (h.oncalo && h.rr < MediandQdxRRMax) endp_dqdx.push_back(h.dqdx); }
    double med_dqdx = -1;
    if (endp_dqdx.size()) {
      unsigned middle = endp_dqdx.size() / 2;
      std::nth_element(endp_dqdx.begin(), endp_dqdx.begin() + middle, endp_dqdx.end());
      med_dqdx = endp_dqdx[middle];
      if (endp_dqdx.size() % 2 == 0) {
        unsigned other_middle = middle - 1;
        std::nth_element(endp_dqdx.begin(), endp_dqdx.begin() + other_middle, endp_dqdx.end());
        med_dqdx = (med_dqdx + endp_dqdx[other_middle]) / 2.; } }
    double EndMediandQdxCut = 1000;
    bool valid_med_dqdx = ((med_dqdx > 0.) && (med_dqdx > EndMediandQdxCut)) || (EndMediandQdxCut < 0.);

    // selection of cosmic stopping muons
    if (verbose_ > 0) std::cout << "hasT0 = " << hasT0
                  << " downwards = " << downwards
                  << " end_is_fid = " << end_is_fid
                  << " time_is_fid = " << time_is_fid
                  << " valid_med_dqdx = " << valid_med_dqdx
                  << std::endl;
    if (selectionMode_ == 0) {
      bool stopping = hasT0 && downwards && end_is_fid && time_is_fid && valid_med_dqdx;
      if (stopping) {
        if (verbose_ > 0) std::cout << "stopping track found!" << std::endl; }
      if (!stopping) continue; }

    // build proxy and 2D hits for MCS fitter
    std::vector<proxy::TrackPointData> pdata; 
    std::vector<recob::Hit> hits2dI1, hits2dI2, hits2dC;
    for (const art::Ptr<recob::Hit>& h : trkHits) {
      const auto& wid = h->WireID();
      if (wid.Plane == 0) hits2dI1.emplace_back(*h);
      if (wid.Plane == 1) hits2dI2.emplace_back(*h);
      if (wid.Plane == 2) hits2dC.emplace_back(*h); }
    auto const& ptrack = proxyTracks[trkPtr.key()];
    for (std::size_t i = 0; i < trkPtr->NPoints(); ++i) {
      pdata.emplace_back(proxy::makeTrackPointData(ptrack, i)); }
    if (trkPtr.key() >= proxyTracks.size()) continue;
    mcsfitter.setPointData(pdata);
    mcsfitter.set2DHitsI1(hits2dI1);
    mcsfitter.set2DHitsI2(hits2dI2);
    mcsfitter.set2DHitsC(hits2dC);

    // compute D3P
    mcsfitter.ComputeD3P(0);
    mcsfitter.ComputeD3P(1);
    mcsfitter.ComputeD3P(2);
    
    // return MCS fit result for each track
    recob::MCSFitResultGS result = mcsfitter.fitMcs(*trkPtr);
    results[trkPtr.key()] = std::move(result); }
  auto output = std::make_unique<std::vector<recob::MCSFitResultGS>>(std::move(results));
  e.put(std::move(output)); }

DEFINE_ART_MODULE(trkf::MCSFitProducerICARUS)
