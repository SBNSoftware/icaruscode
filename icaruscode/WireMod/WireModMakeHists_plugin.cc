//std includes
#include <vector>

//ROOT includes
#include "TFile.h"

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
#include "lardata/Utilities/AssociationUtil.h"

namespace WireMod {
  class MakeHistograms : public art::ResultsProducer {
  public:
    explicit MakeHistograms(fhicl::ParameterSet const& pset);
    ~MakeHistograms() override = default;
    void event(art::Event const& evt) override;
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    art::InputTag fTracksLabel;                // how the tracks are labeled in the input file
    std::vecttor<art::InputTag> fHitsLabelVec; // how the hits are labeled in the input file
  }// end MakeHistograms class

  //-------------------------------------------------------------
  // Define the constructor
  // The way ART works is it will call the construct on the fhicl parameter set
  // Basically we want it to pull out the variables we defined there
  // We're doing this in the reconfigure function (because that makes it clearer what we are doing imo)
  // so just call that function
  MakeHistograms::MakeHistograms(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------------------
  // this function pulls out the stuff we defined in the fhicl file
  void MakeHistograms::reconfigure(fhicl::ParameterSet const& pset)
  {
    // the first arguement is where in the fhicl to look, the second is the default value if that isn't found
    fTracksLabel = pset.get<art::InputTag>("TracksLabel", "pandoraTrackGausCryoE"); // do CryoE and CryoW in fhicl
    fHitsLabelVec = pset.get<std::vector<art::InputTag>>("HitsLabelVec"); // EE & EW or WE & WW
  }

  //-------------------------------------------------------------
  // this function is run on every event in the art file
  // the event stores the information we want to analyze
  void MakeHistograms::event(art::Event const& evt)
  {
    for (auto const& fHitsLabel : fHitsLabelVec)
    {
      // get the trackss out of the event and put them in a handle
      // the handle is needed to get the associated hits
      // we can use the handle to get a vector of tracks
      art::Handle<std::vector<recob::Track>> trackHandle;
      evt.getByLabel(fTracksLabel, trackHandle);
      if (not trackHandle.valid())
      {
        MF_LOG_VERBATIM("WireModMakeHistograms")
          << "Track hand is not valid! abort";
        return;
      }
      std::vector<art::Ptr<recob::Track>> trackPtrVec;
      art::fill_ptr_vector(trackPtrVec, trackHandle);

      // do a similar thing for the hits
      art::Handle<std::vector<recob::Hit>> hitHandle;
      evt.getByLabel(fHitsLabel, hitHandle);
      if (not hitHandle.valid())
      {
        MF_LOG_VERBATIM("WireModMakeHistograms"
          << "Hit hand is not valid! abort";
        return;
      }
      std::vector<art::Ptr<recob::Hit>> hitPtrVec;
      art::fill_ptr_vector(hitPtrVec, hitHandle)

      // this tool will let us get the wire for each hit
      art::FindOneP<recob::Wire> hitToWireAssns(hitHandle, evt, fHitsLabel);

      // this one the hits for each track
      art::FindManyP<recod::Hit> trackToHitsAssns(trackHandle, evt, fTracksLabel);

      // loop over the Tracks
      // these are art::Ptr<recob::Track>, which in most respects can be treated as a pointer to a recob::Track
      for (auto const& trackPtr : trackPtrVec)
      {
        // do checks on the track here, eg is it anode/cathode crossing?

        // get hits for the track
        // the ID identifies the art::Ptr<recob::Track> to get the associated hits
        std::vector<art::Ptr<recob::Hit>> hitsInTrack = trackToHitsAssns.at(trackPtr.ID());

        // loop over the hits
        for (auto const& hitPtr : hitsInTrack)
        {
          // get info from hit here

          // get the associated wire
          art::Ptr<recob::Wire> wirePtr = hitToWireAssns.at(hitPtr.ID());

          // put the waveform in the histogram
        }
      }
    }
  }
 
}// end WireMod namespace
