// std inlcudes
#include <string>
#include <utility>
#include <vector>

// ROOT includes
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"

// wairemlod
#include "sbncode/WireMod/AIML/WireModInfer.hh"

//namespace
namespace wairemlod
{
  // ─────────────────────────────────────────────────────────────────────────
  // TTree row schema (one row per hit)
  //
  // Reco-only feature columns are computed identically for data and MC, so
  // the closure analysis can bin both samples in the same feature space.
  // Direction-dependent features (track angle, dQ/dx with path-length
  // correction) are not available from a single hit and are therefore not
  // stored — the analyzer uses raw `dqdx_uncorr = old_integral / pitch_nom`
  // as a proxy when binning by charge density.
  // ─────────────────────────────────────────────────────────────────────────
  struct HitRow
  {
    // Sample classification
    Bool_t  is_mc;       // backtracker hit ⇒ true; otherwise data/overlay
    Bool_t  ischanged;   // true if model produced a different (I, w)

    // Reco-only features
    Float_t channel;
    UChar_t view;        // 0=U, 1=V, 2=Y (matches geo::View_t)
    UChar_t tpc;         // hit.WireID().TPC
    UChar_t plane;       // hit.WireID().Plane
    UChar_t cryo;        // hit.WireID().Cryostat
    Float_t x_drift;     // [cm] from peak time
    Float_t wire_y;      // [cm] wire-center y
    Float_t wire_z;      // [cm] wire-center z
    Float_t peak_time;   // [ticks]
    Float_t pitch_nom;   // [cm] nominal wire pitch (no track-angle correction)

    // Hit (I, w) before and after model
    Float_t old_integral;
    Float_t old_width;
    Float_t new_integral;
    Float_t new_width;
  };

  class WaireMLodifier : public art::EDProducer
  {
    public:
      explicit WaireMLodifier(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

    private:
      const geo::WireReadoutGeom* fWireReadout = &(art::ServiceHandle<geo::WireReadout const>()->Get());
      const art::ServiceHandle<cheat::BackTrackerService> fBT;
      const art::ServiceHandle<cheat::ParticleInventoryService> fPI;
      std::string fWeightFileName; // there is where we try to grab the weights (full path)
      art::InputTag fHitLabel;  // which hits are we pulling in?

      // Per-hit closure tree. One row per hit, written to the module's
      // TFileDirectory (so it lands inside `WaireMLodTPC??` /
      // `WaireMLodPTTPC??` when the module is configured under those labels).
      TTree*  fHitTree;
      HitRow  fRow;

  }; // end WaireMLodifier class

  //------------------------------------------------
  WaireMLodifier::WaireMLodifier(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void WaireMLodifier::reconfigure(fhicl::ParameterSet const& pset)
  {
    fHitLabel  = pset.get<art::InputTag>("HitLabel");
    fWeightFileName = pset.get<std::string>("WeightFileName");

    // we make these things
    produces<std::vector<recob::Hit>>();

    // --- Per-hit closure tree -------------------
    art::ServiceHandle<art::TFileService> tfs;
    fHitTree = tfs->make<TTree>("hits",
      "WaireMLodifier closure tree (one row per hit)");

    fHitTree->Branch("is_mc",        &fRow.is_mc,        "is_mc/O");
    fHitTree->Branch("ischanged",    &fRow.ischanged,    "ischanged/O");
    fHitTree->Branch("channel",      &fRow.channel,      "channel/F");
    fHitTree->Branch("view",         &fRow.view,         "view/b");
    fHitTree->Branch("tpc",          &fRow.tpc,          "tpc/b");
    fHitTree->Branch("plane",        &fRow.plane,        "plane/b");
    fHitTree->Branch("cryo",         &fRow.cryo,         "cryo/b");
    fHitTree->Branch("x_drift",      &fRow.x_drift,      "x_drift/F");
    fHitTree->Branch("wire_y",       &fRow.wire_y,       "wire_y/F");
    fHitTree->Branch("wire_z",       &fRow.wire_z,       "wire_z/F");
    fHitTree->Branch("peak_time",    &fRow.peak_time,    "peak_time/F");
    fHitTree->Branch("pitch_nom",    &fRow.pitch_nom,    "pitch_nom/F");
    fHitTree->Branch("old_integral", &fRow.old_integral, "old_integral/F");
    fHitTree->Branch("old_width",    &fRow.old_width,    "old_width/F");
    fHitTree->Branch("new_integral", &fRow.new_integral, "new_integral/F");
    fHitTree->Branch("new_width",    &fRow.new_width,    "new_width/F");
  }

  //------------------------------------------------
  void WaireMLodifier::produce(art::Event& evt)
  {
    // get a clock for the event
    const detinfo::DetectorClocksData detClock =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    const detinfo::DetectorPropertiesData detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, detClock);

    // get the old hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    bool valid_handle = evt.getByLabel(fHitLabel, hitHandle);
    if (not valid_handle)
    {
      mf::LogVerbatim("WaireMLodifier")
        << "WaireMLodifier --- Could not grab handle at " << fHitLabel << ", skip.";
      return;
    }
    std::vector<art::Ptr<recob::Hit>> hitPtrVec;
    art::fill_ptr_vector(hitPtrVec, hitHandle);

    // if there are no hits, just don't bother
    if (hitPtrVec.size() == 0)
    {
      mf::LogVerbatim("WaireMLodifier")
        << "WaireMLodifier --- No hits found at " << fHitLabel << ", skip.";
      return;
    }

    // WaireMLodifier
    sys::WaireMLod WM(fWeightFileName);
    std::unique_ptr<std::vector<recob::Hit>> new_hits =
      std::make_unique<std::vector<recob::Hit>>(
        WM.produceNew(hitPtrVec, fBT.get(), fPI.get(), &detClock, fWireReadout)
      );

    if (hitPtrVec.size() != new_hits->size())
      throw std::runtime_error("WaireMLodifier --- Unequal number of input and output hits");

    for (size_t hit_idx = 0; hit_idx < hitPtrVec.size(); ++hit_idx)
    {
      recob::Hit*          newHitPtr = &(new_hits->at(hit_idx));
      art::Ptr<recob::Hit> oldHitPtr =  hitPtrVec.at(hit_idx);

      const float old_integral = oldHitPtr->Integral();
      const float old_width    = oldHitPtr->RMS();
      const float new_integral = newHitPtr->Integral();
      const float new_width    = newHitPtr->RMS();
      const bool  ischanged    = (old_integral != new_integral) || (old_width != new_width);

      // ── Sample classification: backtracker presence ⇒ MC ──────────────────
      bool is_mc = true;
      try
      {
        std::vector<double> hitXYZ = fBT->HitToXYZ(detClock, *oldHitPtr);
        (void)hitXYZ;
      }
      catch (...)
      {
        is_mc = false;
      }

      // ── Reco-only feature values (identical method for data and MC) ───────
      const geo::WireID  wireID  = oldHitPtr->WireID();
      const geo::PlaneID planeID = wireID.asPlaneID();

      const auto wcenter   = fWireReadout->Wire(wireID).GetCenter();
      const float pitch_nom = fWireReadout->Plane(planeID).WirePitch();

      // ── Fill the row ──────────────────────────────────────────────────────
      fRow.is_mc        = is_mc;
      fRow.ischanged    = ischanged;
      fRow.channel      = static_cast<float>(oldHitPtr->Channel());
      fRow.view         = static_cast<unsigned char>(oldHitPtr->View());
      fRow.tpc          = static_cast<unsigned char>(wireID.TPC);
      fRow.plane        = static_cast<unsigned char>(wireID.Plane);
      fRow.cryo         = static_cast<unsigned char>(wireID.Cryostat);
      fRow.x_drift      = static_cast<float>(detProp.ConvertTicksToX(oldHitPtr->PeakTime(), planeID));
      fRow.wire_y       = static_cast<float>(wcenter.Y());
      fRow.wire_z       = static_cast<float>(wcenter.Z());
      fRow.peak_time    = oldHitPtr->PeakTime();
      fRow.pitch_nom    = pitch_nom;
      fRow.old_integral = old_integral;
      fRow.old_width    = old_width;
      fRow.new_integral = new_integral;
      fRow.new_width    = new_width;
      fHitTree->Fill();
    }
    evt.put(std::move(new_hits));
  }
  DEFINE_ART_MODULE(WaireMLodifier)
} // end namespace
