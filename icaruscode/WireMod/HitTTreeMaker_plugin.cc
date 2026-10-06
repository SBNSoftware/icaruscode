//std includes
#include <cmath>
#include <vector>

//ROOT includes
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"

//Framework includes
#include "art/Framework/Core/ResultsProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Results.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

namespace WireMod {
  class HitTTreeMaker : public art::ResultsProducer {
  public:
    explicit HitTTreeMaker(fhicl::ParameterSet const& pset);
    ~HitTTreeMaker() override = default;
    void event(art::Event const& evt) override;
    void reconfigure(fhicl::ParameterSet const& p);
    void reset_values();
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
    const double fPi = std::acos(-1); // get pi so I don't have to care later

    art::InputTag fTag; // how the hits/wires are labeled in the input file
    bool fUseTrackOnly; // do we want to only use hits from tacks?
    TTree* fHitTree;    // where we store the info

    // TTree Branches
    geo::WireID fWireID = geo::WireID(geo::CryostatID::InvalidID, geo::TPCID::InvalidID, geo::PlaneID::InvalidID, geo::WireID::InvalidID); // the ID of the wire the hit was on
    double fHitAmplitude = std::numeric_limits<double>::quiet_NaN(); // how tall was the hit?
    double fHitWidth = std::numeric_limits<double>::quiet_NaN();     // how wide was the hit?
    double fHitGoF = std::numeric_limits<double>::quiet_NaN();       // how good was the hit fit?
    bool fFromTrack = false;                                         // was the hit from a track?
    double fTrackDirX = std::numeric_limits<double>::quiet_NaN();    // track direction (in GLOBAL coordinates)
    double fTrackDirY = std::numeric_limits<double>::quiet_NaN();    // track direction (in GLOBAL coordinates)
    double fTrackDirZ = std::numeric_limits<double>::quiet_NaN();    // track direction (in GLOBAL coordinates)

  };// end HitTTreeMaker class

  //-------------------------------------------------------------
  // Define the constructor
  // The way ART works is it will call the construct on the fhicl parameter set
  // Basically we want it to pull out the variables we defined there
  // We're doing this in the reconfigure function (because that makes it clearer what we are doing imo)
  // so just call that function
  HitTTreeMaker::HitTTreeMaker(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------------------
  // this function pulls out the stuff we defined in the fhicl file
  void HitTTreeMaker::reconfigure(fhicl::ParameterSet const& pset)
  {
    // the first arguement is where in the fhicl to look, the second is the default value if that isn't found
    std::string label   = pset.get<std::string>("Label"  , "gaushitTPCEE");
    std::string process = pset.get<std::string>("Process", "MCstage0Var" );
    fUseTrackOnly = pset.get<bool>("UseTrackOnly", false);
    mf::LogVerbatim("HitTTreeMaker")
      << "Getting hits with label " << label << " from process " << process;
    std::string tpcStr = (fUseTrackOnly) ? label.substr(label.size() - 1) : label.substr(label.size() - 2);
    mf::LogVerbatim("HitTTreeMaker")
      << "  These are in the " << tpcStr << " TPC";
    fTag = art::InputTag(label, {}, process);

    // set up the output TTree
    art::ServiceHandle<art::TFileService> tfs;
    std::string treeName  = label + "_" + process;
    std::string treeTitle = (fUseTrackOnly) ? "Hits From Tracks" : "All Hits";
    fHitTree = tfs->make<TTree>(treeName.c_str(), treeTitle.c_str());
    fHitTree->Branch("WireID"       , &fWireID);
    fHitTree->Branch("Amplitude"    , &fHitAmplitude);
    fHitTree->Branch("Width"        , &fHitWidth);
    fHitTree->Branch("GoodnessOfFit", &fHitGoF);
    fHitTree->Branch("IsTrackHit"   , &fFromTrack);
    fHitTree->Branch("TrackDirX"    , &fTrackDirX);
    fHitTree->Branch("TrackDirY"    , &fTrackDirY);
    fHitTree->Branch("TrackDirZ"    , &fTrackDirZ);
  }

  //-------------------------------------------------------------
  // Reset the values for our TTree
  void HitTTreeMaker::reset_values()
  {
    fWireID = geo::WireID(geo::CryostatID::InvalidID, geo::TPCID::InvalidID, geo::PlaneID::InvalidID, geo::WireID::InvalidID);
    fHitAmplitude = std::numeric_limits<double>::quiet_NaN();
    fHitWidth = std::numeric_limits<double>::quiet_NaN();
    fHitGoF = std::numeric_limits<double>::quiet_NaN();
    fFromTrack = false;
    fTrackDirX = std::numeric_limits<double>::quiet_NaN();
    fTrackDirY = std::numeric_limits<double>::quiet_NaN();
    fTrackDirZ = std::numeric_limits<double>::quiet_NaN();
  }

  //-------------------------------------------------------------
  // this function is run on every event in the art file
  // the event stores the information we want to analyze
  void HitTTreeMaker::event(art::Event const& evt)
  {

    // with or without tracks?
    if (fUseTrackOnly)
    {
      // get the tracks and the associations for the hits
      art::Handle<std::vector<recob::Track>> trackHandle;
      evt.getByLabel(fTag, trackHandle);
      art::FindManyP<recob::Hit, recob::TrackHitMeta> hitsFromTracks(trackHandle, evt, fTag);

      // loop over tracks, get the hits and the metadata
      for (size_t trackIdx = 0; trackIdx < trackHandle->size(); ++trackIdx)
      {
        art::Ptr<recob::Track> track(trackHandle, trackIdx);
        std::vector<art::Ptr<recob::Hit>> trackHits = hitsFromTracks.at(track.key());
        std::vector<const recob::TrackHitMeta*> trackMetas = hitsFromTracks.data(track.key());

        for (size_t hitIdx = 0; hitIdx < trackHits.size(); ++hitIdx)
        {
          art::Ptr<recob::Hit> hit = trackHits.at(hitIdx);
          const recob::TrackHitMeta* meta = trackMetas.at(hitIdx);

          // get the info
          size_t hitInTrack = meta->Index();
          if (hitInTrack == std::numeric_limits<unsigned int>::max() || not track->HasValidPoint(hitInTrack))
            continue;

          fWireID       = hit->WireID();
          fHitAmplitude = hit->PeakAmplitude();
          fHitWidth     = hit->RMS();
          fHitGoF       = hit->GoodnessOfFit();
          fFromTrack    = true;
          fTrackDirX    = track->DirectionAtPoint(hitInTrack).X();
          fTrackDirY    = track->DirectionAtPoint(hitInTrack).Y();
          fTrackDirZ    = track->DirectionAtPoint(hitInTrack).Z();

          fHitTree->Fill();

          this->reset_values();
        }
      }
    } else {
      // get the recob::Hits
      art::Handle<std::vector<recob::Hit>> hitHandle;
      evt.getByLabel(fTag, hitHandle);

      // loop over the wires
      for (auto const& hit : *hitHandle)
      {
        fWireID       = hit.WireID();
        fHitAmplitude = hit.PeakAmplitude();
        fHitWidth     = hit.RMS();
        fHitGoF       = hit.GoodnessOfFit();

        fHitTree->Fill();

        this->reset_values();
      }
    }
  }

  //-------------------------------------------------------------
  // writeResults is currently unused, but we still need to define it
  // just have it do nothing for now
  void HitTTreeMaker::writeResults(art::Results& r)
  {
  }
  
  //-------------------------------------------------------------
  // lear is currently unused, but we still need to define it
  // just have it do nothing for now
  void HitTTreeMaker::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(HitTTreeMaker)
}// end WireMod namespace
