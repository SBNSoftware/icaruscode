////////////////////////////////////////////////////////////////////////
//// Class:       IcarusHDF5Maker
//// Plugin Type: analyzer (art v3_06_03)
//// File:        IcarusHDF5Maker_module.cc
////
//// Generated at Wed May  5 08:23:31 2021 by V Hewes using cetskelgen
/// from cetlib version v3_11_01.
//////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcore/Geometry/WireReadout.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "hep_hpc/hdf5/make_ntuple.hpp"
#include "art/Framework/Services/Registry/ServiceHandle.h"

class IcarusHDF5Maker : public art::EDAnalyzer {
public:
  explicit IcarusHDF5Maker(fhicl::ParameterSet const& p);
  ~IcarusHDF5Maker() noexcept {}; // bare pointers are cleaned up by endSubRun

  IcarusHDF5Maker(IcarusHDF5Maker const&) = delete;
  IcarusHDF5Maker(IcarusHDF5Maker&&) = delete;
  IcarusHDF5Maker& operator=(IcarusHDF5Maker const&) = delete;
  IcarusHDF5Maker& operator=(IcarusHDF5Maker&&) = delete;

  void beginSubRun(art::SubRun const& sr) override;
  void endSubRun(art::SubRun const& sr) override;
  void analyze(art::Event const& e) override;

private:

  std::string fTruthLabel;
  std::string fHitLabel;
  std::string fHitTruthLabel;
  std::string fSPLabel;
  std::string fOpHitLabel;
  std::string fOpFlashLabel;
  std::string fPFParticleLabel;

  bool fUseMap;
  std::string fEventInfo;
  std::string fOutputName;

  hep_hpc::hdf5::File fFile;  ///< output HDF5 file

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>     // event id (run, subrun, event)
  >* fEventNtuple; ///< event ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // is cc
                        hep_hpc::hdf5::Column<int, 1>, // nu pdg
                        hep_hpc::hdf5::Column<float, 1>,  // nu energy
                        hep_hpc::hdf5::Column<float, 1>,  // lep energy
                        hep_hpc::hdf5::Column<float, 1>   // nu dir (x, y, z)
  >* fEventNtupleNu; ///< event ntuple with neutrino information

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // spacepoint id
                        hep_hpc::hdf5::Column<float, 1>,  // 3d position (x, y, z)
                        hep_hpc::hdf5::Column<int, 1>,     // 2d hit (u, v, y)
                        hep_hpc::hdf5::Column<float, 1>     //ChiSquared of the hit 
  >* fSpacePointNtuple; ///< spacepoint ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // hit id
                        hep_hpc::hdf5::Column<float, 1>,  // hit integral
                        hep_hpc::hdf5::Column<float, 1>,  // hit rms
                        hep_hpc::hdf5::Column<int, 1>,    // tpc id
                        hep_hpc::hdf5::Column<int, 1>,    // global plane
                        hep_hpc::hdf5::Column<float, 1>,  // global wire
                        hep_hpc::hdf5::Column<float, 1>,  // global time
                        hep_hpc::hdf5::Column<int, 1>,    // raw plane
                        hep_hpc::hdf5::Column<float, 1>,  // raw wire
                        hep_hpc::hdf5::Column<float, 1>,  // raw time
                        hep_hpc::hdf5::Column<int, 1>    //cryostat
  >* fHitNtuple; ///< hit ntuple


  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // g4 id
                        hep_hpc::hdf5::Column<int, 1>,    // particle type
                        hep_hpc::hdf5::Column<int, 1>,    // parent g4 id
                        hep_hpc::hdf5::Column<int, 1>,    // is from nu
                        hep_hpc::hdf5::Column<float, 1>,  // momentum
                        hep_hpc::hdf5::Column<float, 1>,  // start position (x, y, z)
                        hep_hpc::hdf5::Column<float, 1>,  // end position (x, y, z)
                        hep_hpc::hdf5::Column<std::string, 1>, // start process
                        hep_hpc::hdf5::Column<std::string, 1>  // end process
  >* fParticleNtuple; ///< particle ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // hit id
                        hep_hpc::hdf5::Column<int, 1>,    // g4 id
                        hep_hpc::hdf5::Column<float, 1>,  // deposited energy [ MeV ]
                        hep_hpc::hdf5::Column<float, 1>,  // x position
                        hep_hpc::hdf5::Column<float, 1>,  // y position
                        hep_hpc::hdf5::Column<float, 1>   // z position
  >* fEnergyDepNtuple; ///< energy deposition ntuple

  hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // hit id
			hep_hpc::hdf5::Column<int, 1>,    // hit_channel
			hep_hpc::hdf5::Column<int, 1>,    // wire pos
			hep_hpc::hdf5::Column<float, 1>,  // peaktime
			hep_hpc::hdf5::Column<float, 1>,  // width
			hep_hpc::hdf5::Column<float, 1>,  // area
			hep_hpc::hdf5::Column<float, 1>,  // amplitude
			hep_hpc::hdf5::Column<float, 1>,  // pe
			hep_hpc::hdf5::Column<int, 1>     // usedInFlash
  >* fOpHitNtuple; ///< PMT hit ntuple

 hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // flash id
			hep_hpc::hdf5::Column<int, 1>,    // wire pos
			hep_hpc::hdf5::Column<float, 1>,  // time
			hep_hpc::hdf5::Column<float, 1>,  // time width
			hep_hpc::hdf5::Column<float, 1>,  // Y center
			hep_hpc::hdf5::Column<float, 1>,  // Y width
			hep_hpc::hdf5::Column<float, 1>,  // Z center
			hep_hpc::hdf5::Column<float, 1>,  // Z width
			hep_hpc::hdf5::Column<float, 1>   // totalpe
  >* fOpFlashNtuple; ///< Flash ntuple

hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                        hep_hpc::hdf5::Column<int, 1>,    // sumpe id
                        hep_hpc::hdf5::Column<int, 1>,    // flash id
                        hep_hpc::hdf5::Column<int, 1>,    // PMT channel
  			hep_hpc::hdf5::Column<float, 1>   // pe
  >* fOpFlashSumPENtuple; ///< Flash SumPE ntuple

       



  using ProxyPfpColl_t = decltype(proxy::getCollection<std::vector<recob::PFParticle> >(
											std::declval<art::Event>(),std::declval<art::InputTag>(),
											proxy::withAssociated<larpandoraobj::PFParticleMetadata>(std::declval<art::InputTag>()),
											proxy::withAssociated<recob::Slice>(std::declval<art::InputTag>()),
											proxy::withAssociated<recob::Cluster>(std::declval<art::InputTag>()),
											proxy::withAssociated<recob::Vertex>(std::declval<art::InputTag>())  ));
  using ProxyPfpElem_t = ProxyPfpColl_t::element_proxy_t;

  // proxy to connect cluster to hit
  using ProxyClusColl_t = decltype(proxy::getCollection<std::vector<recob::Cluster>>(
										     std::declval<art::Event>(), std::declval<art::InputTag>(),
										     proxy::withAssociated<recob::Hit>(std::declval<art::InputTag>())));
  using ProxyClusElem_t = ProxyClusColl_t::element_proxy_t;

  std::vector<int> nearwires;
  int NearWire(const geo::WireReadoutGeom& geo, const geo::PlaneID &ip, const float x, const float y, const float z);

  int stitchedPlane(const geo::WireID wid) const {
    int plane = wid.Plane;
    if(wid.TPC==2 || wid.TPC==3) {
      if(wid.Plane==1) plane=2;
      else if(wid.Plane==2) plane=1;
    }
    return plane;
  }
  float stitchedTime(const geo::WireID wid, float timein) const {
    float time = timein;
    if(wid.TPC==2 || wid.TPC==3) {
      //correction = 2*(tpcgeo.DriftDistance()/detProp.DriftVelocity()-clockData.TriggerOffsetTPC())/clockData.TPCClock().TickPeriod() = 6442.15
      time = 6442.15 - timein;
    }
    return time;
  }
  size_t stitchedWire(const geo::WireID wid) const {
    size_t wire = wid.Wire;
    int plane = stitchedPlane(wid);
    if(wid.TPC==1 || wid.TPC==3) {
      if(plane==1 || plane == 2) {
	wire = wid.Wire + 2535; //2535 is the last part of the wires in cryos 0 an 2 before the cut in z=0
      }
    }
    return wire;
  }

};

IcarusHDF5Maker::IcarusHDF5Maker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTruthLabel(p.get<std::string>("TruthLabel")),
    fHitLabel(  p.get<std::string>("HitLabel")),
    fHitTruthLabel(  p.get<std::string>("HitTruthLabel","")),
    fSPLabel(   p.get<std::string>("SPLabel")),
    fOpHitLabel(  p.get<std::string>("OpHitLabel")),
    fOpFlashLabel(  p.get<std::string>("OpFlashLabel")),
    fPFParticleLabel(  p.get<std::string>("PFParticleLabel")),
    fUseMap(    p.get<bool>("UseMap", false)),
    fEventInfo( p.get<std::string>("EventInfo")),
    fOutputName(p.get<std::string>("OutputName"))
{
  if (fEventInfo != "none" && fEventInfo != "nu")
    throw art::Exception(art::errors::Configuration)
      << "EventInfo must be \"none\" or \"nu\", not " << fEventInfo;
}
std::vector<double> tpc_ids_checked = {}; // add the TPC ids here so we only check them once  

void IcarusHDF5Maker::analyze(art::Event const& e) {

  const cheat::BackTrackerService* bt = 0;
  if (!fUseMap) {
    art::ServiceHandle<cheat::BackTrackerService> bt_h;
    bt = bt_h.get();
  }
  geo::WireReadoutGeom const& geo = art::ServiceHandle<geo::WireReadout const>()->Get();

  std::set<int> g4id;
  art::ServiceHandle<cheat::ParticleInventoryService> pi; 
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);
  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  // auto tpcgeo = art::ServiceHandle<geo::Geometry const>()->TPC();
  // std::cout << "tpc drift distance=" << tpcgeo.DriftDistance() << std::endl;
  // std::cout << "drift velocity=" << detProp.DriftVelocity() << std::endl;
  // std::cout << "drift time=" << tpcgeo.DriftDistance()/detProp.DriftVelocity() << " 2*dt/0.4=" << 2*tpcgeo.DriftDistance()/detProp.DriftVelocity()/clockData.TPCClock().TickPeriod() << std::endl;
  // std::cout << "tick period=" << clockData.TPCClock().TickPeriod() << std::endl;
  // std::cout << "tpcTime=" << clockData.TPCTime() << std::endl;
  // std::cout << "triggerOffsetTPC=" << clockData.TriggerOffsetTPC() << std::endl;
  // std::cout << "triggerTime=" << clockData.TriggerTime() << std::endl;
  // std::cout << "correction=" << 2*(tpcgeo.DriftDistance()/detProp.DriftVelocity()-clockData.TriggerOffsetTPC())/clockData.TPCClock().TickPeriod() << std::endl;

  std::cout << "NEW EVENT: " << run <<"  "<< subrun << " "<< event<<std::endl;

  std::array<int, 3> evtID { run, subrun, event };
  // Fill event table
  if (fEventInfo == "none") {
    fEventNtuple->insert( evtID.data() );
    mf::LogInfo("IcarusHDF5Maker") << "Filling event table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2];
  }

  if (fEventInfo == "nu") {
    // Get MC truth
    art::Handle<std::vector<simb::MCTruth>> truthHandle;
    e.getByLabel(fTruthLabel, truthHandle);
    if (!truthHandle.isValid() || truthHandle->size() == 0) {
      throw art::Exception(art::errors::LogicError)
        << "Expected to find exactly one MC truth object!";
    }
    simb::MCNeutrino nutruth = truthHandle->at(0).GetNeutrino();

    auto up = nutruth.Nu().Momentum().Vect().Unit();
    std::array<float, 3> nuMomentum {(float)up.X(),(float)up.Y(),(float)up.Z()};

    fEventNtupleNu->insert( evtID.data(),
			    nutruth.CCNC() == simb::kCC,
			    nutruth.Nu().PdgCode(),
			    nutruth.Nu().E(),
			    nutruth.Lepton().E(),
			    nuMomentum.data()
			    );

    // for (int ip=0;ip<truthHandle->at(0).NParticles();ip++) {
    //   std::cout << "mcp tkid=" << truthHandle->at(0).GetParticle(ip).TrackId() << " pdg=" << truthHandle->at(0).GetParticle(ip).PdgCode() 
    // 		<< " mother=" << truthHandle->at(0).GetParticle(ip).Mother()
    // 		<< " vtx=" << truthHandle->at(0).GetParticle(ip).Vx() << " " << truthHandle->at(0).GetParticle(ip).Vy() << " " << truthHandle->at(0).GetParticle(ip).Vz()
    // 		<< std::endl;
    // }

    mf::LogInfo("IcarusHDF5Maker") << "Filling event table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\nis cc? " << (nutruth.CCNC() == simb::kCC)
                             << ", nu energy " << nutruth.Nu().E()
                             << ", lepton energy " << nutruth.Lepton().E()
                             << "\nnu momentum x " << nuMomentum[0] << ", y "
                             << nuMomentum[1] << ", z " << nuMomentum[2];
  } // if nu event info

  // Get spacepoints from the event record
  art::Handle<std::vector<recob::SpacePoint>> spListHandle;
  std::vector<art::Ptr<recob::SpacePoint>> splist;
  if (e.getByLabel(fSPLabel, spListHandle))  art::fill_ptr_vector(splist, spListHandle);

  // Get hits from the event record
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (e.getByLabel(fHitLabel, hitListHandle)) art::fill_ptr_vector(hitlist, hitListHandle);
  
  // Get assocations from spacepoints to hits
  art::FindManyP<recob::Hit> fmp(spListHandle, e, fSPLabel);

  // Fill spacepoint table
  for (size_t i = 0; i < spListHandle->size(); ++i) {

    art::Ptr<recob::SpacePoint> spacePointPtr(spListHandle,i); 
    std::vector<art::Ptr<recob::Hit>> associatedHits(fmp.at(spacePointPtr.key()));        
    if (associatedHits.size() < 3) {
      std::cout << "I am certain this cannot happen... but here you go, space point with " << associatedHits.size() << " hits" << std::endl;
      exit(1);
      continue;
    }

    std::array<float, 3> pos {(float)splist[i]->XYZ()[0],(float)splist[i]->XYZ()[1],(float)splist[i]->XYZ()[2]};

    std::array<int, 3> hitID { -1, -1, -1 };
    for (size_t j = 0; j < associatedHits.size(); ++j) { 
      int plane = stitchedPlane(associatedHits[j]->WireID());
      //std::cout << "j=" << j << " tpc=" << associatedHits[j]->WireID().TPC << " plane=" << associatedHits[j]->WireID().Plane 
      //          << " plane2=" << plane << " key=" <<  associatedHits[j].key() << std::endl;      
      hitID[plane] = associatedHits[j].key();
    }
    if (hitID[0]==0 && hitID[1]==-1 && hitID[2]==-1) {
      std::cout << "THIS SHOULD NOT HAPPEN -- sp i=" << i << " hit ids=" << hitID[0] << " " << hitID[1] << " " << hitID[2] << std::endl;
      exit(1);
    }
    fSpacePointNtuple->insert(evtID.data(),splist[i]->ID(), pos.data(), hitID.data(),splist[i]->Chisq() );

    mf::LogInfo("IcarusHDF5Maker") << "Filling spacepoint table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\nspacepoint id " << splist[i]->ID()
                             << "\nposition x " << pos[0] << ", y " << pos[1]
                             << ", z " << pos[2]
                             << "\nhit ids " << hitID[0] << ", " << hitID[1]
                             << ", " << hitID[2]
                             << "Chi_squared "  << splist[i]->Chisq();
  } // for spacepoint

  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (fUseMap) {
    hittruth = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, e, fHitTruthLabel));
  }

  // Loop over hits
  for (art::Ptr<recob::Hit> hit : hitlist) {

    // Fill hit table
    geo::WireID wireid = hit->WireID();
    size_t plane = stitchedPlane(wireid);
    double time = stitchedTime(wireid,hit->PeakTime());
    size_t wire = stitchedWire(wireid);
    fHitNtuple->insert(evtID.data(),
		       hit.key(), hit->Integral(), hit->RMS(), wireid.TPC,
		       plane, wire, time,
		       wireid.Plane, wireid.Wire, hit->PeakTime(),wireid.Cryostat
		       );

    mf::LogInfo("IcarusHDF5Maker") << "Filling hit table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\nhit id " << hit.key() << ", integral "
                             << hit->Integral() << ", RMS " << hit->RMS()
                             << ", TPC " << wireid.TPC
                             << "\nglobal plane " << plane << ", global wire "
                             << wire << ", global time " << time
                             << "\nlocal plane " << wireid.Plane
                             << ", local wire " << wireid.Wire
                             << ", local time " << hit->PeakTime()
                             << ", Cryostat " << wireid.Cryostat;

    // Fill energy deposit table
    if (fUseMap) {
      //not supported in icarus at the moment
      exit(1);
      /*
      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
      //loop over particles
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
        g4id.insert(particle_vec[i_p]->TrackId());
        fEnergyDepNtuple->insert(evtID.data(), hit.key(),
                                 particle_vec[i_p]->TrackId(),
                                 match_vec[i_p]->ideFraction,
                                 std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN());
        mf::LogInfo("IcarusHDF5Maker") << "Filling energy deposit table"
				 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
				 << ", event " << evtID[2]
				 << "\nhit id " << hit.key() << ", g4 id "
				 << particle_vec[i_p]->TrackId() << ", energy fraction "
				 << match_vec[i_p]->ideFraction << ", position ("
				 << std::nan << ", " << std::nan << ", " << std::nan << ")";;
      } // for particle 
      */
    }  else {

      std::vector<sim::TrackIDE> h_ides = bt->ChannelToTrackIDEs(clockData, hit->Channel(), hit->StartTick(), hit->EndTick());
      for (auto& tide : h_ides) {
      	int tid = tide.trackID;
      	if (pi->TrackIdToParticle_P(abs(tid))) tid = abs(tid);
        g4id.insert(tid);
        fEnergyDepNtuple->insert(evtID.data(),hit.key(), tid, tide.energyFrac, -1., -1., -1.);
      }

    } // if using microboone map method or not
  } // for hit

  std::set<int> allIDs = g4id; // Copy original so we can safely modify it

  // Add invisible particles to hierarchy
  for (int id : g4id) {
    const simb::MCParticle* p = pi->TrackIdToParticle_P(abs(id));
    while (p!= NULL && p->Mother() != 0 ) {
      allIDs.insert(abs(p->Mother()));
      p = pi->TrackIdToParticle_P(abs(p->Mother()));
    }
  }
  // Loop over true particles and fill table
  for (int id : allIDs) {
    const simb::MCParticle* p = pi->TrackIdToParticle_P(abs(id));
    if(p==NULL) std::cout<< "we have a problem p is null"<<std::endl;
    if (p==NULL) continue;
    auto mct = pi->TrackIdToMCTruth_P(abs(id));
    int fromNu = (mct.isAvailable() ? mct->NeutrinoSet() : 0);
    // std::cout << "all mcp tkid=" << p->TrackId() << " pdg=" << p->PdgCode()
    // 	      << " mother=" << p->Mother()
    // 	      << " vtx=" << p->Vx() << " " << p->Vy() << " " << p->Vz()
    // 	      << " mctruth=" << mct->NeutrinoSet()
    // 	      << std::endl;
    std::array<float, 3> particleStart { (float)p->Vx(), (float)p->Vy(), (float)p->Vz() };
    std::array<float, 3> particleEnd { (float)p->EndX(), (float)p->EndY(), (float)p->EndZ() };
    fParticleNtuple->insert(evtID.data(),
			    abs(id), p->PdgCode(), p->Mother(), fromNu, (float)p->P(),
			    particleStart.data(), particleEnd.data(),
			    p->Process(), p->EndProcess()
			    );
    mf::LogInfo("IcarusHDF5Maker") << "Filling particle table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
			     << "\ng4 id " << abs(id) << ", pdg code "
                             << p->PdgCode() << ", parent " << p->Mother()
                             << ", momentum " << p->P()
                             << "\nparticle start x " << particleStart[0]
                             << ", y " << particleStart[1]
                             << ", z " << particleStart[2]
                             << "\nparticle end x " << particleEnd[0] << ", y "
                             << particleEnd[1] << ", z " << particleEnd[2]
                             << "\nstart process " << p->Process()
                             << ", end process " << p->EndProcess();
  }

  // Get OpFlashs from the event record
  art::Handle< std::vector< recob::OpFlash > > opFlashListHandle;
  std::unique_ptr<art::FindManyP<recob::OpHit> > assocFlashHit;
  std::vector< art::Ptr< recob::OpFlash > > opflashlist;
  if (fOpFlashLabel != "") {
    if (e.getByLabel(fOpFlashLabel, opFlashListHandle)) art::fill_ptr_vector(opflashlist, opFlashListHandle);
    assocFlashHit = std::unique_ptr<art::FindManyP<recob::OpHit> >(new art::FindManyP<recob::OpHit>(opFlashListHandle, e, fOpFlashLabel));
  }

  // get the flash matching the slice we selected
  art::Handle< art::Assns<recob::Slice, recob::OpFlash> > sliceOpFlashAssnsHandle;
  e.getByLabel(fHitLabel, sliceOpFlashAssnsHandle);
  int flkey = -1;
  if (sliceOpFlashAssnsHandle->size()) {
    art::Ptr<recob::OpFlash> fp = sliceOpFlashAssnsHandle->at(0).second;
    flkey = fp.key();
  }

  std::vector<std::vector<int> > flashsumpepmtmap;//this has at most one element, due to the check on flkey

  // Loop over opflashs, pick only the one we care about
  for (art::Ptr< recob::OpFlash > opflash : opflashlist) {
    if (int(opflash.key())!=flkey) continue;
    std::vector<float> pes;
    for (auto pe : opflash->PEs()) pes.push_back(pe);
    std::vector<int> nearwires;
    double xyz[3] = {0.,opflash->YCenter(),opflash->ZCenter()};
    for (auto const& p : geo.Iterate<geo::PlaneID>()) nearwires.push_back(NearWire(geo,p,xyz[0],xyz[1],xyz[2]));

    // Fill opflash table
    fOpFlashNtuple->insert(evtID.data(),
			   opflash.key(),
			   nearwires.data(),
			   opflash->Time(),opflash->TimeWidth(),
			   opflash->YCenter(),opflash->YWidth(),opflash->ZCenter(),opflash->ZWidth(),
			   opflash->TotalPE()
			   );
    size_t count = 0;
    std::vector<int> sumpepmtmap(art::ServiceHandle<geo::Geometry>()->NOpDets(),-1);
    for (size_t ipmt=0;ipmt<pes.size();ipmt++) {
      if (pes[ipmt]<=0.) continue;
      fOpFlashSumPENtuple->insert(evtID.data(),count,opflash.key(),ipmt,pes[ipmt]);
      sumpepmtmap[ipmt] = count;
      count++;
    }
    flashsumpepmtmap.push_back(sumpepmtmap);
    mf::LogInfo("IcarusHDF5Maker") << "Filling opflash table"
			     << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			     << ", event " << evtID[2]
			     << "\nflash id " << opflash.key() << ", Time " << opflash->Time()
			     << ", TotalPE " << opflash->TotalPE()//;
			     << "\nWireCenters size0 " << opflash->WireCenters().size()//;
			     << "\nYCenter " << opflash->YCenter()<< " ZCenter " << opflash->ZCenter()
			     << "\nYWidth " << opflash->YWidth()<< " ZWidth " << opflash->ZWidth()
			     << "\nInBeamFrame " << opflash->InBeamFrame()<< " OnBeamTime " << opflash->OnBeamTime()
			     << "\nAbsTime " << opflash->AbsTime() << " TimeWidth " << opflash->TimeWidth() << " FastToTotal " << opflash->FastToTotal();
  }

  // Get OpHits from the event record
  art::Handle< std::vector< recob::OpHit > > opHitListHandle;
  std::vector< art::Ptr< recob::OpHit > > ophitlist;
  if (fOpHitLabel != "") {
    if (e.getByLabel(fOpHitLabel, opHitListHandle))
      art::fill_ptr_vector(ophitlist, opHitListHandle);
  }

  // Loop over ophits
  for (art::Ptr< recob::OpHit > ophit : ophitlist) {
    std::vector<int> nearwires;
    auto xyz = geo.OpDetGeoFromOpChannel(ophit->OpChannel()).GetCenter();
    for (auto const& p : geo.Iterate<geo::PlaneID>()) nearwires.push_back(NearWire(geo,p,xyz.X(),xyz.Y(),xyz.Z()));

    bool isInFlash = false;
    if (flkey>=0) {
      auto ophv = assocFlashHit->at(flkey);
      isInFlash = (std::find(ophv.begin(),ophv.end(),ophit)!=ophv.end());
    }
    // Fill ophit table
    fOpHitNtuple->insert(evtID.data(),ophit.key(),
			 ophit->OpChannel(),nearwires.data(),
			 ophit->PeakTime(),ophit->Width(),
			 ophit->Area(),ophit->Amplitude(),ophit->PE(),
			 (isInFlash ? flashsumpepmtmap[0][ophit->OpChannel()] : -1)
			 );
    mf::LogInfo("IcarusHDF5Maker") << "\nFilling ophit table"
			     << "\nrun " << evtID[0] << ", subrun " << evtID[1]
			     << ", event " << evtID[2]
			     << "\nhit id " << ophit.key() << ", channel "
			     << ophit->OpChannel() << ", PeakTime " << ophit->PeakTime()
			     << ", Width " << ophit->Width()
			     << "\narea " << ophit->Area() << ", amplitude "
			     << ophit->Amplitude() << ", PE " << ophit->PE();
  }
  //End optical analyzer

} // function IcarusHDF5Maker::analyze

void IcarusHDF5Maker::beginSubRun(art::SubRun const& sr) {

  struct timeval now;
  gettimeofday(&now, NULL);
  // Open HDF5 output
  std::ostringstream fileName;
  fileName << fOutputName << "_r" << std::setfill('0') << std::setw(5) << sr.run()
	   << "_s" << std::setfill('0') << std::setw(5) << sr.subRun() << "_ts" << std::setw(6) << now.tv_usec << ".h5";

  fFile = hep_hpc::hdf5::File(fileName.str(), H5F_ACC_TRUNC);

  if (fEventInfo == "none")
    fEventNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "event_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3)
    ));
  if (fEventInfo == "nu")
    fEventNtupleNu = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "event_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("is_cc"),
      hep_hpc::hdf5::make_scalar_column<int>("nu_pdg"),    
      hep_hpc::hdf5::make_scalar_column<float>("nu_energy"),
      hep_hpc::hdf5::make_scalar_column<float>("lep_energy"),
      hep_hpc::hdf5::make_column<float>("nu_dir", 3)
    ));

  fSpacePointNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "spacepoint_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("spacepoint_id"),
      hep_hpc::hdf5::make_column<float>("position", 3),
      hep_hpc::hdf5::make_column<int>("hit_id", 3),
      hep_hpc::hdf5::make_scalar_column<float>("Chi_squared")

  ));

  fHitNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "hit_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
      hep_hpc::hdf5::make_scalar_column<float>("integral"),
      hep_hpc::hdf5::make_scalar_column<float>("rms"),
      hep_hpc::hdf5::make_scalar_column<int>("tpc"),
      hep_hpc::hdf5::make_scalar_column<int>("global_plane"),
      hep_hpc::hdf5::make_scalar_column<float>("global_wire"),
      hep_hpc::hdf5::make_scalar_column<float>("global_time"),
      hep_hpc::hdf5::make_scalar_column<int>("local_plane"),
      hep_hpc::hdf5::make_scalar_column<float>("local_wire"),
      hep_hpc::hdf5::make_scalar_column<float>("local_time"),
      hep_hpc::hdf5::make_scalar_column<int>("Cryostat")

  ));

  fParticleNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "particle_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("g4_id"),
      hep_hpc::hdf5::make_scalar_column<int>("type"),
      hep_hpc::hdf5::make_scalar_column<int>("parent_id"),
      hep_hpc::hdf5::make_scalar_column<int>("from_nu"),
      hep_hpc::hdf5::make_scalar_column<float>("momentum"),
      hep_hpc::hdf5::make_column<float>("start_position", 3),
      hep_hpc::hdf5::make_column<float>("end_position", 3),
      hep_hpc::hdf5::make_scalar_column<std::string>("start_process"),
      hep_hpc::hdf5::make_scalar_column<std::string>("end_process")
  ));

  fEnergyDepNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "edep_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
      hep_hpc::hdf5::make_scalar_column<int>("g4_id"),
      hep_hpc::hdf5::make_scalar_column<float>("energy"),
      hep_hpc::hdf5::make_scalar_column<float>("x_position"),
      hep_hpc::hdf5::make_scalar_column<float>("y_position"),
      hep_hpc::hdf5::make_scalar_column<float>("z_position")
  ));

  fOpHitNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "ophit_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
      hep_hpc::hdf5::make_scalar_column<int>("hit_channel"),
      hep_hpc::hdf5::make_column<int>("wire_pos", 3),//3 views
      hep_hpc::hdf5::make_scalar_column<float>("peaktime"),
      hep_hpc::hdf5::make_scalar_column<float>("width"),
      hep_hpc::hdf5::make_scalar_column<float>("area"),
      hep_hpc::hdf5:: make_scalar_column<float>("amplitude"),
      hep_hpc::hdf5::make_scalar_column<float>("pe"),
      hep_hpc::hdf5::make_scalar_column<int>("sumpe_id")
  ));

  fOpFlashNtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "opflash_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("flash_id"),
      hep_hpc::hdf5::make_column<int>("wire_pos", 3),//3 views
      hep_hpc::hdf5::make_scalar_column<float>("time"),
      hep_hpc::hdf5::make_scalar_column<float>("time_width"),
      hep_hpc::hdf5::make_scalar_column<float>("y_center"),
      hep_hpc::hdf5::make_scalar_column<float>("y_width"),
      hep_hpc::hdf5::make_scalar_column<float>("z_center"),
      hep_hpc::hdf5::make_scalar_column<float>("z_width"),
      hep_hpc::hdf5::make_scalar_column<float>("totalpe")
  ));

  fOpFlashSumPENtuple = new hep_hpc::hdf5::Ntuple(
      hep_hpc::hdf5::make_ntuple({fFile, "opflashsumpe_table", 1000},
      hep_hpc::hdf5::make_column<int>("event_id", 3),
      hep_hpc::hdf5::make_scalar_column<int>("sumpe_id"),
      hep_hpc::hdf5::make_scalar_column<int>("flash_id"),
      hep_hpc::hdf5::make_scalar_column<int>("pmt_channel"),
      hep_hpc::hdf5::make_scalar_column<float>("sumpe")
  ));
}

void IcarusHDF5Maker::endSubRun(art::SubRun const& sr) {
  if (fEventInfo == "none") delete fEventNtuple;
  if (fEventInfo == "nu") delete fEventNtupleNu;
  delete fSpacePointNtuple;
  delete fHitNtuple;
  delete fParticleNtuple;
  delete fEnergyDepNtuple;
  delete fOpHitNtuple;
  delete fOpFlashNtuple;
  delete fOpFlashSumPENtuple;
  fFile.close();
}

int IcarusHDF5Maker::NearWire(const geo::WireReadoutGeom& geo, const geo::PlaneID &ip, const float x, const float y, const float z)
//art::ServiceHandle<geo::Geometry> geo;
{
  geo::PlaneGeo const& plane = geo.Plane(ip);
  geo::WireID wireID;
  try {
    wireID = plane.NearestWireID(geo::Point_t(x,y,z));
  }
  catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane.ClosestWireID(e.suggestedWireID());
  }
  return wireID.Wire;
}

DEFINE_ART_MODULE(IcarusHDF5Maker)
