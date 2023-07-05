//std includes
#include <vector>

//ROOT includes
#include "TFile.h"
#include "TH1.h"

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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

namespace WireMod {
  class WireModMakeHists : public art::ResultsProducer {
  public:
    explicit WireModMakeHists(fhicl::ParameterSet const& pset);
    ~WireModMakeHists() override = default;
    void event(art::Event const& evt) override;
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    art::InputTag fHitsLabel; // how the hits are labeled in the input file
  };// end WireModMakeHists class

  //-------------------------------------------------------------
  // Define the constructor
  // The way ART works is it will call the construct on the fhicl parameter set
  // Basically we want it to pull out the variables we defined there
  // We're doing this in the reconfigure function (because that makes it clearer what we are doing imo)
  // so just call that function
  WireModMakeHists::WireModMakeHists(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------------------
  // this function pulls out the stuff we defined in the fhicl file
  void WireModMakeHists::reconfigure(fhicl::ParameterSet const& pset)
  {
    // the first arguement is where in the fhicl to look, the second is the default value if that isn't found
    fHitsLabel= pset.get<art::InputTag>("HitsLabel", "decon1droi"); // EE & EW or WE & WW
  }

  //-------------------------------------------------------------
  // this function is run on every event in the art file
  // the event stores the information we want to analyze
  void WireModMakeHists::event(art::Event const& evt)
  {
    // this is what will make our histograms for us
    art::ServiceHandle<art::TFileService>  tfs;

    // get the hits out of the event and put them in a handle
    // the handle is needed to get the associated wires
    // we can use the handle to get a vector of hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    evt.getByLabel(fHitsLabel, hitHandle);
    if (not hitHandle.isValid())
    {
      MF_LOG_VERBATIM("WireModWireModMakeHists")
        << "Hit handle is not valid!" << '\n'
        << "Tried " << fHitsLabel << '\n'
        << "abort";
      return;
    }
    std::vector<art::Ptr<recob::Hit>> hitPtrVec;
    art::fill_ptr_vector(hitPtrVec, hitHandle);

    // this tool will let us get the wire for each hit
    art::FindOneP<recob::Wire> hitToWireAssns(hitHandle, evt, fHitsLabel);

    // loop over the hits
    for (auto const& hitPtr : hitPtrVec)
    {
      // get info from hit here

      // get the associated wire
      art::Ptr<recob::Wire> wirePtr = hitToWireAssns.at(hitPtr.key());

      // put the waveform in the histogram
      // tfs will make a whatever is in the <>, (in this case a TH1F)
      // the agruements past to it should be the same as for the constructor for the <object>
      TH1F* waveformHist = tfs->make<TH1F>(("adc_"+fHitsLabel.label()+"_"+std::to_string(wirePtr.key())).c_str(), //> name of the object
                                           ";Sample;Arbitrary Units",                                             //> axes labels
                                           wirePtr->NSignal(),                                                    //> numbeer of bins
                                           0,                                                                     //> lowest edge
                                           wirePtr->NSignal());                                                   //> upper edge

      // ROOT counts from 1, everyone else counts from 0
      for (size_t bin = 1; bin < wirePtr->NSignal() + 1; ++bin)
      {
        waveformHist->SetBinContent(bin, (wirePtr->Signal())[bin-1]);
      }

      // In testing this I just want one
      break;
    }
  }

  //-------------------------------------------------------------
  // writeResults is currently unused, but we still need to define it
  // just have it do nothing for now
  void WireModMakeHists::writeResults(art::Results& r)
  {
  }
  
  //-------------------------------------------------------------
  // lear is currently unused, but we still need to define it
  // just have it do nothing for now
  void WireModMakeHists::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(WireModMakeHists)
}// end WireMod namespace
