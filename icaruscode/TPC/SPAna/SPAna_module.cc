////////////////////////////////////////////////////////////////////////
// Class:       SPAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        SPAna_module.cc
//
// Generated at Sun Jun 23 05:48:04 2024 by Gray Putnam using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TTree.h"
#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbnobj/Common/Calibration/TrackCaloSkimmerObj.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

class SPAna;


class SPAna : public art::EDAnalyzer {
public:
  explicit SPAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SPAna(SPAna const&) = delete;
  SPAna(SPAna&&) = delete;
  SPAna& operator=(SPAna const&) = delete;
  SPAna& operator=(SPAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void Clear();
  void setBranches();
  void setMetaBranches(TTree *);

  sbn::HitInfo MakeHit(const recob::Hit &hit,
    unsigned hkey,
    const art::Ptr<recob::SpacePoint> &sp,
    const geo::GeometryCore *geo,
    const detinfo::DetectorClocksData &dclock,
    const cheat::BackTrackerService *bt_serv);

  std::pair<std::vector<sbn::TrueHit>, std::vector<sbn::Vector3D>> ParticleTrueHits(
    const simb::MCParticle &particle,
    const std::map<int, std::vector<std::tuple<geo::WireID, short, const sim::IDE*>>> id_to_ide_map,
    const detinfo::DetectorPropertiesData &dprop,
    const geo::GeometryCore *geo);

  void BuildHitRanges(const std::vector<art::Ptr<recob::Hit>> &hits);
  void BuildWireRanges(const std::vector<art::Ptr<recob::Wire>> &wires);

  bool FindHit(unsigned channel, unsigned tick);
  bool FindWire(unsigned channel, unsigned tick);

  void respondToOpenInputFile(const art::FileBlock& fb) override {
    (void) fb;
    _fileno ++;
  }

private:
  std::vector<art::InputTag> fHitProducers;
  std::vector<art::InputTag> fWireProducers;
  art::InputTag fG4Producer;
  art::InputTag fSimChannelProducer;

  TTree *tRecoHits;
  TTree *tRecoWire;
  TTree *tTrueHits;

  int _fileno;
  int _run;
  int _subrun;
  int _evt;
  int _cryo;

  std::vector<sbn::HitInfo> _reco_hits;
  std::vector<sbn::TrueHit> _true_hits;

  std::vector<float> _true_hit_dir_x;
  std::vector<float> _true_hit_dir_y;
  std::vector<float> _true_hit_dir_z;
  std::vector<int> _true_hits_has_hit;
  std::vector<int> _true_hits_has_wire;

  std::vector<int> _reco_wire_start;
  std::vector<int> _reco_wire_end;
  std::vector<int> _reco_wire_wire;
  std::vector<int> _reco_wire_plane;
  std::vector<int> _reco_wire_channel;
  std::vector<float> _reco_wire_charge;
  std::vector<float> _reco_wire_true_charge;
  std::vector<float> _reco_wire_true_energy;

  std::map<unsigned, std::vector<std::pair<unsigned, unsigned>>> fHitRanges;
  std::map<unsigned, std::vector<std::pair<unsigned, unsigned>>> fWireRanges;
};


SPAna::SPAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  tRecoHits = tfs->make<TTree>("RecoHits", "Reconstructed hits");
  tRecoWire = tfs->make<TTree>("RecoWire", "Reconstructed wires");
  tTrueHits = tfs->make<TTree>("TrueHits", "True hits");

  setBranches();

  _fileno = 0;
  _cryo = p.get<int>("Cryostat");

  fHitProducers = p.get<std::vector<art::InputTag>>("HitProducers");
  fWireProducers = p.get<std::vector<art::InputTag>>("WireProducers");
  fG4Producer = p.get<art::InputTag>("G4Producer");
  fSimChannelProducer = p.get<art::InputTag>("SimChannelProducer");
}

void SPAna::BuildHitRanges(const std::vector<art::Ptr<recob::Hit>> &hits) {
  fHitRanges.clear();
  for (const art::Ptr<recob::Hit> &h: hits) fHitRanges[h->Channel()].push_back({h->StartTick(), h->EndTick()});
}

void SPAna::BuildWireRanges(const std::vector<art::Ptr<recob::Wire>> &wires) {
  fWireRanges.clear();
  for (const art::Ptr<recob::Wire> &w: wires) {
    unsigned channel = w->Channel();
    for (auto const &range: w->SignalROI().get_ranges()) {
      fWireRanges[channel].push_back({range.begin_index(), range.end_index()});
    }
  }
}
  
bool FindRange(const std::map<unsigned, std::vector<std::pair<unsigned, unsigned>>> &ranges, unsigned channel, unsigned tick) {
  if (!ranges.count(channel)) return false;

  const std::vector<std::pair<unsigned, unsigned>> &thisrange = ranges.at(channel);
  for (const std::pair<unsigned, unsigned> &r: thisrange) {
    if (tick >= r.first && tick <= r.second) return true;
  }

  return false;
}

bool SPAna::FindHit(unsigned channel, unsigned tick) {
  return FindRange(fHitRanges, channel, tick);
}

bool SPAna::FindWire(unsigned channel, unsigned tick) {
  return FindRange(fWireRanges, channel, tick);
}

void SPAna::setBranches() {
  setMetaBranches(tRecoHits);
  setMetaBranches(tRecoWire);
  setMetaBranches(tTrueHits);

  tRecoHits->Branch("h", &_reco_hits);

  tTrueHits->Branch("t", &_true_hits);
  tTrueHits->Branch("t.dir_x", &_true_hit_dir_x);
  tTrueHits->Branch("t.dir_y", &_true_hit_dir_y);
  tTrueHits->Branch("t.dir_z", &_true_hit_dir_z);
  tTrueHits->Branch("t.has_hit", &_true_hits_has_hit);
  tTrueHits->Branch("t.has_wire", &_true_hits_has_wire);

  tRecoWire->Branch("w.start", &_reco_wire_start);
  tRecoWire->Branch("w.end", &_reco_wire_end);
  tRecoWire->Branch("w.wire", &_reco_wire_wire);
  tRecoWire->Branch("w.plane", &_reco_wire_plane);
  tRecoWire->Branch("w.channel", &_reco_wire_channel);
  tRecoWire->Branch("w.charge", &_reco_wire_charge);
  tRecoWire->Branch("w.true_charge", &_reco_wire_true_charge);
  tRecoWire->Branch("w.true_energy", &_reco_wire_true_energy);

}

void SPAna::setMetaBranches(TTree *t) {
  t->Branch("fileno", &_fileno, "fileno/i");
  t->Branch("run", &_run, "run/i");
  t->Branch("sub", &_subrun, "sub/i");
  t->Branch("evt", &_evt, "evt/i");
  t->Branch("cryo", &_cryo, "cryo/i");
}

void SPAna::Clear() {
  _reco_hits.clear();
  _true_hits.clear();

  _true_hit_dir_x.clear();
  _true_hit_dir_y.clear();
  _true_hit_dir_z.clear();
  _true_hits_has_hit.clear();
  _true_hits_has_wire.clear();

  _reco_wire_start.clear();
  _reco_wire_end.clear();
  _reco_wire_wire.clear();
  _reco_wire_plane.clear();
  _reco_wire_channel.clear();
  _reco_wire_charge.clear();
  _reco_wire_true_charge.clear();
  _reco_wire_true_energy.clear();
}

std::map<int, std::vector<std::tuple<geo::WireID, short, const sim::IDE*>>> PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels, const geo::GeometryCore &geo) {
    std::map<int, std::vector<std::tuple<geo::WireID, short, const sim::IDE*>>> ret;

    for (const art::Ptr<sim::SimChannel> sc : simchannels) {
      // Lookup the wire of this channel
      raw::ChannelID_t channel = sc->Channel();
      std::vector<geo::WireID> maybewire = geo.ChannelToWire(channel);
      geo::WireID thisWire; // Default constructor makes invalid wire
      if (maybewire.size()) thisWire = maybewire[0];

      for (const auto &item : sc->TDCIDEMap()) {
        for (const sim::IDE &ide: item.second) {
          // indexing initializes empty vector
          ret[abs(ide.trackID)].push_back({thisWire, item.first, &ide});
        }
      }
    }
    return ret;
}

// Turn a particle position to a space-charge induced position
geo::Point_t TrajectoryToWirePosition(const geo::Point_t &loc, const geo::TPCID &tpc) {
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  art::ServiceHandle<geo::Geometry const> geom;

  geo::Point_t ret = loc;

  // Returned X is the drift -- multiply by the drift direction to undo this
  int corr = geom->TPC(tpc).DriftDir().X();
  
  if (sce && sce->EnableSimSpatialSCE()) {
    geo::Vector_t offset = sce->GetPosOffsets(ret);
  
    ret.SetX(ret.X() + corr * offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
}

// Turn a space-charge induced position to a trajectory Position
geo::Point_t WireToTrajectoryPosition(const geo::Point_t &loc, const geo::TPCID &tpc) {
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t ret = loc;

  if (sce && sce->EnableSimSpatialSCE()) {
    geo::Vector_t offset = sce->GetCalPosOffsets(ret, tpc.TPC);

    ret.SetX(ret.X() + offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
  
}

std::pair<std::vector<sbn::TrueHit>, std::vector<sbn::Vector3D>> SPAna::ParticleTrueHits(
    const simb::MCParticle &particle,
    const std::map<int, std::vector<std::tuple<geo::WireID, short, const sim::IDE*>>> id_to_ide_map,
    const detinfo::DetectorPropertiesData &dprop,
    const geo::GeometryCore *geo) {

  // Organize deposition info into per-wire true "Hits" -- key is the Channel Number
  std::map<unsigned, sbn::TrueHit> truehits; 

  const std::vector<std::tuple<geo::WireID, short, const sim::IDE *>> empty;
  const std::vector<std::tuple<geo::WireID, short, const sim::IDE *>> &particle_ides = id_to_ide_map.count(particle.TrackId()) ? id_to_ide_map.at(particle.TrackId()) : empty;

  for (auto const &ide_tup: particle_ides) {
    const geo::WireID &w = std::get<0>(ide_tup);
    short tick = std::get<1>(ide_tup);

    // cut on cryostat
    if ((int)w.Cryostat != _cryo) continue;

    unsigned c = geo->PlaneWireToChannel(w);
    const sim::IDE *ide = std::get<2>(ide_tup);

    // Set stuff
    truehits[c].cryo = w.Cryostat;
    truehits[c].tpc = w.TPC;
    truehits[c].plane = w.Plane;
    truehits[c].wire = w.Wire;
    truehits[c].channel = c;

    // Average stuff using charge-weighting
    float old_elec = truehits[c].nelec;
    float new_elec = old_elec + ide->numElectrons;
    truehits[c].p.x = (truehits[c].p.x*old_elec + ide->x*ide->numElectrons) / new_elec;
    truehits[c].p.y = (truehits[c].p.y*old_elec + ide->y*ide->numElectrons) / new_elec;
    truehits[c].p.z = (truehits[c].p.z*old_elec + ide->z*ide->numElectrons) / new_elec;
    truehits[c].time = (truehits[c].time*old_elec + tick*ide->numElectrons) / new_elec;

    // Also get the position with space charge un-done
    geo::Point_t ide_p(ide->x, ide->y, ide->z);
    geo::Point_t ide_p_scecorr = WireToTrajectoryPosition(ide_p, w);

    truehits[c].p_scecorr.x = (truehits[c].p_scecorr.x*old_elec + ide_p_scecorr.x()*ide->numElectrons) / new_elec; 
    truehits[c].p_scecorr.y = (truehits[c].p_scecorr.y*old_elec + ide_p_scecorr.y()*ide->numElectrons) / new_elec; 
    truehits[c].p_scecorr.z = (truehits[c].p_scecorr.z*old_elec + ide_p_scecorr.z()*ide->numElectrons) / new_elec; 
    
    // Sum stuff
    truehits[c].nelec += ide->numElectrons;
    truehits[c].e += ide->energy;
    truehits[c].ndep += 1;
  }

  // Compute widths
  for (auto const &ide_tup: particle_ides) {
    const geo::WireID &w = std::get<0>(ide_tup);

    // cut on cryostat
    if ((int)w.Cryostat != _cryo) continue;

    unsigned c = geo->PlaneWireToChannel(w);
    const sim::IDE *ide = std::get<2>(ide_tup);

    geo::Point_t ide_p(ide->x, ide->y, ide->z);
    geo::Point_t ide_p_scecorr = WireToTrajectoryPosition(ide_p, w);

    // Average stuff using charge-weighting
    float this_elec = ide->numElectrons;

    truehits[c].p_width.x += (ide_p.x() - truehits[c].p.x) * (ide_p.x() - truehits[c].p.x) * this_elec / truehits[c].nelec;
    truehits[c].p_width.y += (ide_p.y() - truehits[c].p.y) * (ide_p.y() - truehits[c].p.y) * this_elec / truehits[c].nelec;
    truehits[c].p_width.z += (ide_p.z() - truehits[c].p.z) * (ide_p.z() - truehits[c].p.z) * this_elec / truehits[c].nelec;

    truehits[c].p_scecorr_width.x += (ide_p_scecorr.x() - truehits[c].p_scecorr.x) * (ide_p_scecorr.x() - truehits[c].p_scecorr.x) * this_elec / truehits[c].nelec;
    truehits[c].p_scecorr_width.y += (ide_p_scecorr.y() - truehits[c].p_scecorr.y) * (ide_p_scecorr.y() - truehits[c].p_scecorr.y) * this_elec / truehits[c].nelec;
    truehits[c].p_scecorr_width.z += (ide_p_scecorr.z() - truehits[c].p_scecorr.z) * (ide_p_scecorr.z() - truehits[c].p_scecorr.z) * this_elec / truehits[c].nelec;
  }

  // Convert to vector
  std::vector<sbn::TrueHit> truehits_v;
  for (auto const &p: truehits) {
    truehits_v.push_back(p.second);
  }

  // Save true directions
  std::vector<sbn::Vector3D> truehitdirs;

  // Compute the time of each hit
  for (sbn::TrueHit &h: truehits_v) {
    // TODO: fix magic number
    h.time -= 2900; // == (G4RefTime - TriggerOffsetTPC)/TickPeriod = (1500 - 340)/0.4
    double xdrift = abs(h.p.x - geo->Plane(geo::PlaneID(h.cryo, h.tpc, 0)).GetCenter().X());
    h.tdrift = xdrift / dprop.DriftVelocity(); 
  }

  // Compute the pitch of each hit and order it in the trajectory
  for (sbn::TrueHit &h: truehits_v) {
    // Use the SCE-undone hit since this matches to the Trajectory
    TVector3 h_p(h.p_scecorr.x, h.p_scecorr.y, h.p_scecorr.z);

    TVector3 direction;
    float closest_dist = -1.;
    int traj_index = -1;
    for (unsigned i_traj = 0; i_traj < particle.NumberTrajectoryPoints(); i_traj++) {
      if (closest_dist < 0. || (particle.Position(i_traj).Vect() - h_p).Mag() < closest_dist) {
        direction = particle.Momentum(i_traj).Vect().Unit();
        closest_dist = (particle.Position(i_traj).Vect() - h_p).Mag();
        traj_index = i_traj;
      }
    }

    sbn::Vector3D dir_v; 
    dir_v.x = direction.x();
    dir_v.y = direction.y();
    dir_v.z = direction.z();
    truehitdirs.push_back(dir_v);

    // If we got a direction, get the pitch
    if (closest_dist >= 0. && direction.Mag() > 1e-4) {
      geo::PlaneID plane(h.cryo, h.tpc, h.plane);
      float angletovert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();
      float cosgamma = abs(cos(angletovert) * direction.Z() + sin(angletovert) * direction.Y());
      float pitch = geo->WirePitch(plane) / cosgamma;
      h.pitch = pitch;
    }
    else {
      h.pitch = -1.;
    }
    // And the pitch induced by SCE
    if (closest_dist >= 0. && direction.Mag() > 1e-4) {
      geo::PlaneID plane(h.cryo, h.tpc, h.plane);
      float angletovert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();

      TVector3 loc_mdx_v = h_p - direction * (geo->WirePitch(geo->View(plane)) / 2.);
      TVector3 loc_pdx_v = h_p + direction * (geo->WirePitch(geo->View(plane)) / 2.);

      // Convert types for helper functions
      geo::Point_t loc_mdx(loc_mdx_v.X(), loc_mdx_v.Y(), loc_mdx_v.Z());
      geo::Point_t loc_pdx(loc_pdx_v.X(), loc_pdx_v.Y(), loc_pdx_v.Z());
      geo::Point_t h_p_point(h_p.X(), h_p.Y(), h_p.Z());

      loc_mdx = TrajectoryToWirePosition(loc_mdx, plane);
      loc_pdx = TrajectoryToWirePosition(loc_pdx, plane);
      
      // Direction at wires
      geo::Vector_t dir = (loc_pdx - loc_mdx) /  (loc_mdx - loc_pdx).r(); 

      // Pitch at wires
      double cosgamma = std::abs(std::sin(angletovert)*dir.Y() + std::cos(angletovert)*dir.Z());
      double pitch;
      if (cosgamma) {
        pitch = geo->WirePitch(geo->View(plane))/cosgamma;
      }
      else {
        pitch = 0.;
      }

      // Now bring that back to the particle trajectory
      geo::Point_t loc_w = TrajectoryToWirePosition(h_p_point, plane);
      
      geo::Point_t locw_pdx_traj = WireToTrajectoryPosition(loc_w + pitch*dir, plane);
      geo::Point_t loc = WireToTrajectoryPosition(loc_w, plane);
      
      h.pitch_sce = (locw_pdx_traj - loc).R();
    }
    else {
      h.pitch_sce = -1.;
    }

    // And the trajectory location
    h.itraj = traj_index;

    // And the residual range of the hit
    h.rr = 0.;
    if (traj_index >= 0) {
      for (int i_traj = traj_index+1; i_traj < (int)particle.NumberTrajectoryPoints(); i_traj++) {
        h.rr += (particle.Position(i_traj).Vect() - particle.Position(i_traj-1).Vect()).Mag();
      }

      // Also account for the distance from the Hit point to the matched trajectory point
      double hit_distance_along_particle = (h_p - particle.Position(traj_index).Vect()).Dot(particle.Momentum(traj_index).Vect().Unit());
      h.rr += -hit_distance_along_particle;
    }
  }

  return {truehits_v, truehitdirs};

}

sbn::HitInfo SPAna::MakeHit(const recob::Hit &hit,
    unsigned hkey,
    const art::Ptr<recob::SpacePoint> &sp,
    const geo::GeometryCore *geo,
    const detinfo::DetectorClocksData &dclock,
    const cheat::BackTrackerService *bt_serv) {

  // TrackHitInfo to save
  sbn::HitInfo h;

  // information from the hit object
  h.integral = hit.Integral();
  h.sumadc = hit.ROISummedADC();
  h.width = hit.RMS();
  h.time = hit.PeakTime();
  h.mult = hit.Multiplicity();
  h.wire = hit.WireID().Wire;
  h.plane = hit.WireID().Plane;
  h.channel = geo->PlaneWireToChannel(hit.WireID());
  h.tpc = hit.WireID().TPC;
  h.end = hit.EndTick();
  h.start = hit.StartTick();
  h.id = (int)hkey;

  // Do back-tracking on each hit
  if (bt_serv) {
    // The default BackTracking function goes from (peak - width, peak + width).
    //
    // This time range does not match well hits with a non-Gaussian shape where
    // the Gaussian-fit-width does not replicate the width of the pulse. 
    //
    // Instead, we use the Hit (start, end) time range. This is also more relevant
    // for (e.g.) the SummedADC charge extraction method.
    //
    // Don't use this:
    // std::vector<sim::TrackIDE> ides = bt_serv->HitToTrackIDEs(dclock, hit);
    //
    // Use this:
    std::vector<sim::TrackIDE> ides = bt_serv->ChannelToTrackIDEs(dclock, hit.Channel(), hit.StartTick(), hit.EndTick());

    h.truth.e = 0.;
    h.truth.nelec = 0.;

    for (const sim::TrackIDE &ide: ides) {
      h.truth.e += ide.energy;
      h.truth.nelec += ide.numElectrons;
    }
  }
  else {
    h.truth.e = -1.;
    h.truth.nelec = -1.;
  }


  // Save SpacePoint information
  if (sp) {
    h.sp.x = sp->position().x();
    h.sp.y = sp->position().y();
    h.sp.z = sp->position().z();

    h.hasSP = true;
  }
  else {
    h.hasSP = false;
  }

  return h;
}

void SPAna::analyze(art::Event const& e)
{
  unsigned evt = e.event();
  unsigned sub = e.subRun();
  unsigned run = e.run();
  _evt = evt;
  _subrun = sub;
  _run = run;

  std::cout << "[SPAna::analyze] Run: " << run << ", SubRun: " << sub << ", Event: "<< evt << ", Is Data: " << e.isRealData() << std::endl;

  Clear();

  // Load services
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;


  // Collect products

  // Hits
  std::vector<art::Ptr<recob::Hit>> hitList;
  for (const art::InputTag &t: fHitProducers) {
    art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(t);
    art::fill_ptr_vector(hitList, hitHandle);
  }
  // Add spacepoints later?
  // art::FindManyP<recob::SpacePoint> allHitSPs(hitList, e, fPFPproducer);

  // Wires
  std::vector<art::Ptr<recob::Wire>> wireList;
  for (const art::InputTag &t: fWireProducers) {
    art::ValidHandle<std::vector<recob::Wire>> wireHandle = e.getValidHandle<std::vector<recob::Wire>>(t);
    art::fill_ptr_vector(wireList, wireHandle);
  }

  // Truth
  std::vector<art::Ptr<simb::MCParticle>> mcparticles;
  art::ValidHandle<std::vector<simb::MCParticle>> mcparticle_handle = e.getValidHandle<std::vector<simb::MCParticle>>(fG4Producer);
  art::fill_ptr_vector(mcparticles, mcparticle_handle);

  std::vector<art::Ptr<sim::SimChannel>> simchannels;
  art::ValidHandle<std::vector<sim::SimChannel>> simchannel_handle = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelProducer);
  art::fill_ptr_vector(simchannels, simchannel_handle);

  // Prep matching info
  std::map<int, std::vector<std::tuple<geo::WireID, short, const sim::IDE*>>> id_to_ide_map = PrepSimChannels(simchannels, *geometry);
  BuildHitRanges(hitList);
  BuildWireRanges(wireList);

  // Save hits
  for (unsigned i = 0; i < hitList.size(); i++) {
    _reco_hits.push_back(MakeHit(*hitList[i], hitList[i].key(), {} /*allHitSps.at(hitList[i].key())*/, geometry, clock_data, bt_serv.get()));
  }

  // Save wires
  for (unsigned i = 0; i < wireList.size(); i++) {
    art::Ptr<recob::Wire> wire = wireList[i];

    unsigned channel = wire->Channel();
    unsigned plane_id = geometry->ChannelToWire(wire->Channel()).at(0).Plane;
    unsigned wire_id = geometry->ChannelToWire(wire->Channel()).at(0).Wire;
    for (auto const &range: wire->SignalROI().get_ranges()) {

      _reco_wire_start.push_back(range.begin_index());
      _reco_wire_end.push_back(range.end_index());
      _reco_wire_wire.push_back(wire_id);
      _reco_wire_plane.push_back(plane_id);
      _reco_wire_channel.push_back(channel);

      float charge = 0;
      for (float val: range) charge += val;
      _reco_wire_charge.push_back(charge);

      float true_charge = -1;
      float true_energy = -1;
      // Do back-tracking on each wire
      if (bt_serv.get()) {
        std::vector<sim::TrackIDE> ides = bt_serv->ChannelToTrackIDEs(clock_data, channel, range.begin_index(), range.end_index());

        true_charge = 0;
        true_energy = 0;

        for (const sim::TrackIDE &ide: ides) {
          true_energy += ide.energy;
          true_charge += ide.numElectrons;
        }
      }  
      _reco_wire_true_charge.push_back(true_charge);
      _reco_wire_true_energy.push_back(true_energy);
    }
  }

  // Save truth info
  for (unsigned i = 0; i < mcparticles.size(); i++) {
    const simb::MCParticle &p = *mcparticles[i];
    auto ret = ParticleTrueHits(p, id_to_ide_map, dprop, geometry);

    for (const sbn::TrueHit &th: ret.first) _true_hits.push_back(th);
    for (const sbn::Vector3D &d: ret.second) {
      _true_hit_dir_x.push_back(d.x);
      _true_hit_dir_y.push_back(d.y);
      _true_hit_dir_z.push_back(d.z);
    }
  }

  // true to reco matching
  for (const sbn::TrueHit &th: _true_hits) {
    _true_hits_has_hit.push_back((int)FindHit(th.channel, th.time));
    _true_hits_has_wire.push_back((int)FindWire(th.channel, th.time));
  }

  tRecoHits->Fill();
  tTrueHits->Fill();
  tRecoWire->Fill();

}

DEFINE_ART_MODULE(SPAna)
