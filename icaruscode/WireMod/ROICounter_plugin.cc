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
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

namespace WireMod {
  class ROICounter : public art::ResultsProducer {
  public:
    explicit ROICounter(fhicl::ParameterSet const& pset);
    ~ROICounter() override = default;
    void event(art::Event const& evt) override;
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    art::InputTag fLabel;      // how the hits/wires are labeled in the input file
    TH1F* fHist;               // store N ROI per channel
    TH2F* fHeight;             // store a hist for the heights
    TH2F* fWidth;              // store a hist for the widths
    TProfile* fHeightProf;     // store a profile for the heights
    TProfile* fWidthProf;      // store a profile for the widths
  };// end ROICounter class

  //-------------------------------------------------------------
  // Define the constructor
  // The way ART works is it will call the construct on the fhicl parameter set
  // Basically we want it to pull out the variables we defined there
  // We're doing this in the reconfigure function (because that makes it clearer what we are doing imo)
  // so just call that function
  ROICounter::ROICounter(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------------------
  // this function pulls out the stuff we defined in the fhicl file
  void ROICounter::reconfigure(fhicl::ParameterSet const& pset)
  {
    // the first arguement is where in the fhicl to look, the second is the default value if that isn't found
    fLabel = pset.get<art::InputTag>("Label", "decon1droi");

    art::ServiceHandle<art::TFileService> tfs;
    fHist       = tfs->make<TH1F>    ("ROI_vs_Channel"        , ";Channel;Number of ROIs", 55296, -0.5, 55295.5);
    fHeight     = tfs->make<TH2F>    ("Height_vs_Channel"     , ";Channel;ROI Height"    , 55296, -0.5, 55295.5, 1000, 0, 200);
    fWidth      = tfs->make<TH2F>    ("Width_vs_Channel"      , ";Channel;ROI Width"     , 55296, -0.5, 55295.5, 1000, 0, 200);
    fHeightProf = tfs->make<TProfile>("Height_vs_Channel_prof", ";Channel;ROI Height"    , 55296, -0.5, 55295.5);
    fWidthProf  = tfs->make<TProfile>("Width_vs_Channel_prof" , ";Channel;ROI Width"     , 55296, -0.5, 55295.5);
  }

  //-------------------------------------------------------------
  // this function is run on every event in the art file
  // the event stores the information we want to analyze
  void ROICounter::event(art::Event const& evt)
  {
    // loop for the TPCs
    std::string instances[4] = {"PHYSCRATEDATATPCEE",
                                "PHYSCRATEDATATPCEW",
                                "PHYSCRATEDATATPCWE",
                                "PHYSCRATEDATATPCWW"};

    for (auto const& inst : instances)
    {
      // get the tag
      art::InputTag fullTag(fLabel.label(), inst);

      // get the Wires
      art::Handle< std::vector<recob::Wire> > wireHandle;
      evt.getByLabel(fullTag, wireHandle);

      // loop over the wires
      for (auto const& wire : *wireHandle)
      {
        // get the channel the wire lives on
        raw::ChannelID_t channel = wire.Channel();

        // loop over ROIs on wire
        for (size_t range_idx = 0; range_idx < wire.SignalROI().n_ranges(); ++range_idx)
        {
          // couting ROIs...
          fHist->Fill(channel);

          // get ROI
          std::vector<float> roi = wire.SignalROI().range(range_idx).data();
          if (roi.size() == 0)
            continue;
          
          // make a TH1 and fit it with a Gaussian
          TH1F roiHist("roiHist", ";Tick;Abritrary Units", roi.size(), 0, roi.size());
          TF1* gaus = new TF1("gausFunc", "gaus", 0, roi.size());
          
          // loop to get min and max
          float  roi_min = std::numeric_limits<float>::max();
          float  roi_max = std::numeric_limits<float>::min();
          float  roi_int = 0;
          size_t roi_bin = 0;
          for (auto const& roi_val : roi)
          {
            roi_min = (roi_min > roi_val) ? roi_val : roi_min;
            roi_max = (roi_max < roi_val) ? roi_val : roi_max;
            roi_int += roi_val;
            roiHist.SetBinContent(roi_bin, roi_val);
            ++roi_bin;
          }
          
          // seed and fit
          float heightSeed = roi_max - roi_min;
          float widthSeed  = roi_int / (std::sqrt(2*std::acos(-1)) * heightSeed);
          float meanSeed   = roi.size() / 2.0;
          gaus->SetParameters(heightSeed, widthSeed, meanSeed);
          TFitResultPtr gausFit = roiHist.Fit(gaus, "S");

          // fill profiles
          // proxy the height with the difference between min and max
          // and the width as the number of ticks in the ROI
          //fHeightProf->Fill(channel, roi_max - roi_min);
          //fWidthProf ->Fill(channel, roi.size());
          if ((gausFit.Get()     != nullptr) &&
              (gausFit->Status() ==       0)   )
          {
            mf::LogVerbatim("ROICounter")
              << "Filling Hist";
            fHeight    ->Fill(channel, gausFit->Parameter(0));
            fWidth     ->Fill(channel, gausFit->Parameter(2));
            fHeightProf->Fill(channel, gausFit->Parameter(0));
            fWidthProf ->Fill(channel, gausFit->Parameter(2));
          }
        }
      }
    }
  }

  //-------------------------------------------------------------
  // writeResults is currently unused, but we still need to define it
  // just have it do nothing for now
  void ROICounter::writeResults(art::Results& r)
  {
  }
  
  //-------------------------------------------------------------
  // lear is currently unused, but we still need to define it
  // just have it do nothing for now
  void ROICounter::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(ROICounter)
}// end WireMod namespace
