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
    art::InputTag fLabel; // how the hits/wires are labeled in the input file
    bool fGetHits;        // are we getting hits? if false the label is for the wires
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
    fLabel = pset.get<art::InputTag>("Label", "decon1droi");
    fGetHits = pset.get<bool>("GetHits", false);
  }

  //-------------------------------------------------------------
  // this function is run on every event in the art file
  // the event stores the information we want to analyze
  void WireModMakeHists::event(art::Event const& evt)
  {
    // this is what will make our histograms for us
    art::ServiceHandle<art::TFileService>  tfs;

    // get a unique string for this event
    std::string evtStr = std::to_string(evt.id().run()) + "_"
                       + std::to_string(evt.id().subRun()) + "_"
                       + std::to_string(evt.id().event()) + "_";

    if (fGetHits)
    {
      // get the hits out of the event and put them in a handle
      // the handle is needed to get the associated wires
      // we can use the handle to get a vector of hits
      art::Handle<std::vector<recob::Hit>> hitHandle;
      evt.getByLabel(fLabel, hitHandle);
      if (not hitHandle.isValid())
      {
        MF_LOG_VERBATIM("WireModWireModMakeHists")
          << "Hit handle is not valid!" << '\n'
          << "Tried " << fLabel << '\n'
          << "abort";
        return;
      }
      std::vector<art::Ptr<recob::Hit>> hitPtrVec;
      art::fill_ptr_vector(hitPtrVec, hitHandle);

      // this tool will let us get the wire for each hit
      art::FindOneP<recob::Wire> hitToWireAssns(hitHandle, evt, fLabel);

      // loop over the hits
      for (auto const& hitPtr : hitPtrVec)
      {
        // get info from hit here
        // like the start and stop times of the hit
        size_t fstTick = hitPtr->StartTick();
        size_t endTick = hitPtr->EndTick();
        size_t hitWidth = endTick - fstTick;

        // get the associated wire
        art::Ptr<recob::Wire> wirePtr = hitToWireAssns.at(hitPtr.key());

        // the start and end ticks were aquired assuming the hits are gaussian
        // we want to plot a bit of buffer around the region of interest to get the full shape
        // default to a hit width on either side, but if the wire doesn't have enough ticks then use what we can
        size_t nTicks = wirePtr->NSignal();
        size_t fstBuffer = (fstTick > hitWidth)          ? hitWidth : fstTick;
        size_t endBuffer = (nTicks > endTick + hitWidth) ? hitWidth : nTicks - endTick;

        // put the waveform in the histogram
        // tfs will make a whatever is in the <>, (in this case a TH1F)
        // the agruements past to it should be the same as for the constructor for the <object>
        TH1F* waveformHist = tfs->make<TH1F>(("adc_"+evtStr+fLabel.label()+"_"+std::to_string(wirePtr.key())).c_str(), //> name of the object
                                             ";Sample;Arbitrary Units",                                                //> axes labels
                                             hitWidth + fstBuffer + endBuffer,                                         //> number of bins
                                             fstTick - fstBuffer,                                                      //> lowest edge
                                             endTick + endBuffer);                                                     //> upper edge

        // ROOT counts from 1, everyone else counts from 0
        //for (size_t bin = 1; bin < endTick - fstTick + 1; ++bin)
        for (size_t bin = 1; bin < hitWidth + fstBuffer + endBuffer + 1; ++bin)
        {
          float wvfmVal = (wirePtr->Signal())[fstTick - fstBuffer + bin - 1]; 
          waveformHist->SetBinContent(bin, wvfmVal);
        }

        // In testing this I just want one
        break;
      }
    } else {
      // get the wires directly since we aren't getting the hits
      art::Handle<std::vector<recob::Wire>> wireHandle;
      evt.getByLabel(fLabel, wireHandle);
      if (not wireHandle.isValid())
      {
        MF_LOG_VERBATIM("WireModWireModMakeHists")
          << "Wire handle is not valid!" << '\n'
          << "Tried " << fLabel << '\n'
          << "abort";
        return;
      }
      std::vector<art::Ptr<recob::Wire>> wirePtrVec;
      art::fill_ptr_vector(wirePtrVec, wireHandle);

      // loop over the wires
      for (auto const& wirePtr : wirePtrVec)
      {
        // put the waveform in the histogram
        // tfs will make a whatever is in the <>, (in this case a TH1F)
        // the agruements past to it should be the same as for the constructor for the <object>
        TH1F* waveformHist = tfs->make<TH1F>(("adc_"+evtStr+fLabel.label()+":"+fLabel.instance()+"_"+std::to_string(wirePtr.key())).c_str(), //> name of the object
                                             ";Sample;Arbitrary Units",                                                                      //> axes labels
                                             wirePtr->NSignal(),                                                                             //> numbeer of bins
                                             0,                                                                                              //> lowest edge
                                             wirePtr->NSignal());                                                                            //> upper edge

        // ROOT counts from 1, everyone else counts from 0
        for (size_t bin = 1; bin < wirePtr->NSignal() + 1; ++bin)
        {
          waveformHist->SetBinContent(bin, (wirePtr->Signal())[bin-1]);
        }

        // In testing this I just want one
        break;
      }
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
