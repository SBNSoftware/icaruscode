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
  class HitComp : public art::ResultsProducer {
  public:
    explicit HitComp(fhicl::ParameterSet const& pset);
    ~HitComp() override = default;
    void event(art::Event const& evt) override;
    void reconfigure(fhicl::ParameterSet const& p);
    void writeResults(art::Results& r) override;
    void clear() override;

  private:
    art::InputTag fTag;         // how the hits/wires are labeled in the input file
    TH1F* fHitHeightPlane0;     // how tall do the hits look? (1st Ind)
    TH1F* fHitWidthPlane0;      // how wide do the hits look? (1st Ind)
    TH1F* fHitChi2perDOFPlane0; // how well are the hits fit? (1st Ind)
    TH1F* fHitHeightPlane1;     // how tall do the hits look? (2nd Ind)
    TH1F* fHitWidthPlane1;      // how wide do the hits look? (2nd Ind)
    TH1F* fHitChi2perDOFPlane1; // how well are the hits fit? (2nd Ind)
    TH1F* fHitHeightPlane2;     // how tall do the hits look? (Collection)
    TH1F* fHitWidthPlane2;      // how wide do the hits look? (Collection)
    TH1F* fHitChi2perDOFPlane2; // how well are the hits fit? (Collection)
  };// end HitComp class

  //-------------------------------------------------------------
  // Define the constructor
  // The way ART works is it will call the construct on the fhicl parameter set
  // Basically we want it to pull out the variables we defined there
  // We're doing this in the reconfigure function (because that makes it clearer what we are doing imo)
  // so just call that function
  HitComp::HitComp(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------------------
  // this function pulls out the stuff we defined in the fhicl file
  void HitComp::reconfigure(fhicl::ParameterSet const& pset)
  {
    // the first arguement is where in the fhicl to look, the second is the default value if that isn't found
    std::string label   = pset.get<std::string>("Label"  , "gaushitTPCEE");
    std::string process = pset.get<std::string>("Process", "pMCstage0Var");
    mf::LogVerbatim("HitComp")
      << "Getting hits with label " << label << " from process " << process;
    std::string tpcStr = label.substr(label.size() - 2);
    mf::LogVerbatim("HitComp")
      << "  These are in the " << tpcStr << " TPC";
    fTag = art::InputTag(label, {}, process);

    art::ServiceHandle<art::TFileService> tfs;
    fHitHeightPlane0     = tfs->make<TH1F>(("HitHeight_"    +process+"_"+label+"_1stInd").c_str(), (tpcStr+" 1^{st}_{} Induction;Hit Peak Amplitude (ADC Units);Number of Hits").c_str(), 1000, 0, 100);
    fHitWidthPlane0      = tfs->make<TH1F>(("HitWidth_"     +process+"_"+label+"_1stInd").c_str(), (tpcStr+" 1^{st}_{} Induction;Hit Width (Ticks);Number of Hits"             ).c_str(), 1000, 0, 100);
    fHitChi2perDOFPlane0 = tfs->make<TH1F>(("HitChi2perDOF_"+process+"_"+label+"_1stInd").c_str(), (tpcStr+" 1^{st}_{} Induction;Hit Height #chi^{2}_{}/DOF;Number of Hits"    ).c_str(), 1000, 0, 100);
    fHitHeightPlane1     = tfs->make<TH1F>(("HitHeight_"    +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" 2^{nd}_{} Induction;Hit Peak Amplitude (ADC Units);Number of Hits").c_str(), 1000, 0, 100);
    fHitWidthPlane1      = tfs->make<TH1F>(("HitWidth_"     +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" 2^{nd}_{} Induction;Hit Width (Ticks);Number of Hits"             ).c_str(), 1000, 0, 100);
    fHitChi2perDOFPlane1 = tfs->make<TH1F>(("HitChi2perDOF_"+process+"_"+label+"_2ndInd").c_str(), (tpcStr+" 2^{nd}_{} Induction;Hit Height #chi^{2}_{}/DOF;Number of Hits"    ).c_str(), 1000, 0, 100);
    fHitHeightPlane2     = tfs->make<TH1F>(("HitHeight_"    +process+"_"+label+"_Collec").c_str(), (tpcStr+" Collection;Hit Peak Amplitude (ADC Units);Number of Hits"         ).c_str(), 1000, 0, 100);
    fHitWidthPlane2      = tfs->make<TH1F>(("HitWidth_"     +process+"_"+label+"_Collec").c_str(), (tpcStr+" Collection;Hit Width (Ticks);Number of Hits"                      ).c_str(), 1000, 0, 100);
    fHitChi2perDOFPlane2 = tfs->make<TH1F>(("HitChi2perDOF_"+process+"_"+label+"_Collec").c_str(), (tpcStr+" Collection;Hit Height #chi^{2}_{}/DOF;Number of Hits"             ).c_str(), 1000, 0, 100);
  }

  //-------------------------------------------------------------
  // this function is run on every event in the art file
  // the event stores the information we want to analyze
  void HitComp::event(art::Event const& evt)
  {
    // DEBUGGING
    mf::LogVerbatim("HitComp")
      << "List InputTags for std::vector<recob::Hit>" << '\n' << "looking for" << fTag;
    for (auto const& tag : evt.getInputTags<std::vector<recob::Hit>>())
    {
      mf::LogVerbatim("HitComp")
        << "  " << tag;
      if (tag == fTag)
      {
        mf::LogVerbatim("HitComp")
          << "    MATCH!";
      } else {
        mf::LogVerbatim("HitComp")
          << "    NOT A MATCH";
      }
    }
    // DEBUGGING

    // get the recob::Hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    evt.getByLabel(fTag, hitHandle);

    // loop over the wires
    for (auto const& hit : *hitHandle)
    {
      geo::PlaneID planeID = hit.WireID().planeID();
      if        (planeID.Plane == 0)
      {
        fHitHeightPlane0    ->Fill(hit.PeakAmplitude());
        fHitWidthPlane0     ->Fill(hit.RMS());
        fHitChi2perDOFPlane0->Fill(hit.GoodnessOfFit());
      } else if (planeID.Plane == 1) {
        fHitHeightPlane1    ->Fill(hit.PeakAmplitude());
        fHitWidthPlane1     ->Fill(hit.RMS());
        fHitChi2perDOFPlane1->Fill(hit.GoodnessOfFit());
      } else if (planeID.Plane == 2) {
        fHitHeightPlane2    ->Fill(hit.PeakAmplitude());
        fHitWidthPlane2     ->Fill(hit.RMS());
        fHitChi2perDOFPlane2->Fill(hit.GoodnessOfFit());
      }
    }
  }

  //-------------------------------------------------------------
  // writeResults is currently unused, but we still need to define it
  // just have it do nothing for now
  void HitComp::writeResults(art::Results& r)
  {
  }
  
  //-------------------------------------------------------------
  // lear is currently unused, but we still need to define it
  // just have it do nothing for now
  void HitComp::clear()
  {
  }

  DEFINE_ART_RESULTS_PLUGIN(HitComp)
}// end WireMod namespace
