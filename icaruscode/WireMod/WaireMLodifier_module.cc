// std inlcudes
#include <string>
#include <vector>

// ROOT includes
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2F.h"

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "sbnobj/ICARUS/TPC/ChannelROI.h"
#include "icaruscode/TPC/Utilities/ChannelROICreator.h"

// wairemlod
#include "sbncode/WireMod/AIML/WireModInfer.hh"

//namespace
namespace wairemlod
{

  class WaireMLodifier : public art::EDProducer
  {
    public:
      explicit WaireMLodifier(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

    private:
      const geo::WireReadoutGeom* fWireReadout = &(art::ServiceHandle<geo::WireReadout const>()->Get());
      const art::ServiceHandle<cheat::BackTrackerService> fBT;
      std::string fWeightFileName; // there is where we try to grab the weights (full path)
      art::InputTag fHitLabel;  // which hits are we pulling in?
      TH2F* integrals; // store the new/old integrals for the modified hits
      TH2F* widths; // store the new/old widths for the modified hits

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

    // setup output histograms
    art::ServiceHandle<art::TFileService> tfs;
    integrals = tfs->make<TH2F>("integrals", ";Original Integral;Modified Integral",
                                1000, 0, 1000, 1000, 0, 1000);
    widths    = tfs->make<TH2F>("widths",    ";Original Width;Modified Width",
                                1000, 0, 1000, 1000, 0, 1000);
  }

  //------------------------------------------------
  void WaireMLodifier::produce(art::Event& evt)
  {
    // get a clock for the event
    const detinfo::DetectorClocksData detClock =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // get the old hits
    art::Handle< std::vector<recob::Hit> > hitHandle;
    bool valid_handle = evt.getByLabel(fHitLabel, hitHandle);
    if (not valid_handle)
    {
      mf::LogVerbatim("WaireMLodifier")
        << "WaireMLodifier --- Could not grab handle at " << fHitLabel << ", skip.";
      return;
    }
    auto const& hitVec(*hitHandle);

    // if there are no hits, just don't bother
    if (hitVec.size() == 0)
    {
      mf::LogVerbatim("WaireMLodifier")
        << "WaireMLodifier --- No hits found at " << fHitLabel << ", skip.";
      return;
    }

    // WaireMLodifier
    sys::WaireMLod WM(fWeightFileName);
    std::unique_ptr<std::vector<recob::Hit>> new_hits =
      std::make_unique<std::vector<recob::Hit>>(WM.produceNew(hitVec, fBT.get(), &detClock, fWireReadout));

    if (hitVec.size() != new_hits->size())
      throw std::runtime_error("WaireMLodifier --- Unequal number of input and output hits");

    for (size_t hit_idx = 0; hit_idx < hitVec.size(); ++hit_idx)
    {
      float old_integral = hitVec.at(hit_idx).Integral();
      float old_width    = hitVec.at(hit_idx).RMS();
      float new_integral = new_hits->at(hit_idx).Integral();
      float new_width    = new_hits->at(hit_idx).RMS();
      if (old_integral == new_integral && old_width == new_width)
        continue;
      integrals->Fill(old_integral, new_integral);
      widths->Fill(old_width, new_width);
    }

    evt.put(std::move(new_hits));
  }
  DEFINE_ART_MODULE(WaireMLodifier)
} // end namespace
