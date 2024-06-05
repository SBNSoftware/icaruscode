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
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
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
    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
    const double fPi = std::acos(-1); // get pi so I don't have to care later

    art::InputTag fTag;            // how the hits/wires are labeled in the input file
    TH1D* fHitHeightPlane0;        // how tall do the hits look? (1st Ind)
    TH1D* fHitWidthPlane0;         // how wide do the hits look? (1st Ind)
    TH1D* fHitChi2perDOFPlane0;    // how well are the hits fit? (1st Ind)
    TH2D* fHitHeightVsWidthPlane0; // 2D height vs width distributions (1st Ind)
    TH1D* fHitHeightPlane1;        // how tall do the hits look? (2nd Ind)
    TH1D* fHitWidthPlane1;         // how wide do the hits look? (2nd Ind)
    TH1D* fHitChi2perDOFPlane1;    // how well are the hits fit? (2nd Ind)
    TH2D* fHitHeightVsWidthPlane1; // 2D height vs width distributions (2nd Ind)
    TH1D* fHitHeightPlane2;        // how tall do the hits look? (Collection)
    TH1D* fHitWidthPlane2;         // how wide do the hits look? (Collection)
    TH1D* fHitChi2perDOFPlane2;    // how well are the hits fit? (Collection)
    TH2D* fHitHeightVsWidthPlane2; // 2D height vs width distributions (Collection)
    bool fUseTrackOnly;            // do we want to only use hits from tacks?
    TH2D* fHitHeightVsThXYPlane0;  // 2D height vs ThXY distributions (1st Ind)
    TH2D* fHitHeightVsThYZPlane0;  // 2D height vs ThYZ distributions (1st Ind)
    TH2D* fHitWidthVsThXYPlane0;   // 2D width vs ThXY distributions (1st Ind)
    TH2D* fHitWidthVsThYZPlane0;   // 2D width vs ThYZ distributions (1st Ind)
    TH2D* fHitHeightVsThXYPlane1;  // 2D height vs ThXY distributions (2nd Ind)
    TH2D* fHitHeightVsThYZPlane1;  // 2D height vs ThYZ distributions (2nd Ind)
    TH2D* fHitWidthVsThXYPlane1;   // 2D width vs ThXY distributions (2nd Ind)
    TH2D* fHitWidthVsThYZPlane1;   // 2D width vs ThYZ distributions (2nd Ind)
    TH2D* fHitHeightVsThXYPlane2;  // 2D height vs ThXY distributions (Collection)
    TH2D* fHitHeightVsThYZPlane2;  // 2D height vs ThYZ distributions (Collection)
    TH2D* fHitWidthVsThXYPlane2;   // 2D width vs ThXY distributions (Collection)
    TH2D* fHitWidthVsThYZPlane2;   // 2D width vs ThYZ distributions (Collection)

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
    fUseTrackOnly = pset.get<bool>("UseTrackOnly", false);
    mf::LogVerbatim("HitComp")
      << "Getting hits with label " << label << " from process " << process;
    std::string tpcStr = (fUseTrackOnly) ? label.substr(label.size() - 1) : label.substr(label.size() - 2);
    mf::LogVerbatim("HitComp")
      << "  These are in the " << tpcStr << " TPC";
    fTag = art::InputTag(label, {}, process);

    art::ServiceHandle<art::TFileService> tfs;
    if (fUseTrackOnly)
    {
      fHitHeightPlane0        = tfs->make<TH1D>(("HitHeight_"       +process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;Hit Peak Amplitude (ADC Units);Number of Hits"          ).c_str(), 1000, 0, 100);
      fHitWidthPlane0         = tfs->make<TH1D>(("HitWidth_"        +process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;Hit Width (Ticks);Number of Hits"                       ).c_str(), 1000, 0, 100);
      fHitChi2perDOFPlane0    = tfs->make<TH1D>(("HitChi2perDOF_"   +process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;Hit Height #chi^{2}_{}/DOF;Number of Hits"              ).c_str(), 1000, 0, 100);
      fHitHeightVsWidthPlane0 = tfs->make<TH2D>(("HitHeightVsWidth_"+process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;Hit Width (Ticks);Hit Peak Amplitude (ADC Units)"       ).c_str(), 1000, 0, 100, 1000, 0, 100);
      fHitHeightPlane1        = tfs->make<TH1D>(("HitHeight_"       +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;Hit Peak Amplitude (ADC Units);Number of Hits"          ).c_str(), 1000, 0, 100);
      fHitWidthPlane1         = tfs->make<TH1D>(("HitWidth_"        +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;Hit Width (Ticks);Number of Hits"                       ).c_str(), 1000, 0, 100);
      fHitChi2perDOFPlane1    = tfs->make<TH1D>(("HitChi2perDOF_"   +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;Hit Height #chi^{2}_{}/DOF;Number of Hits"              ).c_str(), 1000, 0, 100);
      fHitHeightVsWidthPlane1 = tfs->make<TH2D>(("HitHeightVsWidth_"+process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;Hit Width (Ticks);Hit Peak Amplitude (ADC Units)"       ).c_str(), 1000, 0, 100, 1000, 0, 100);
      fHitHeightPlane2        = tfs->make<TH1D>(("HitHeight_"       +process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;Hit Peak Amplitude (ADC Units);Number of Hits"                   ).c_str(), 1000, 0, 100);
      fHitWidthPlane2         = tfs->make<TH1D>(("HitWidth_"        +process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;Hit Width (Ticks);Number of Hits"                                ).c_str(), 1000, 0, 100);
      fHitChi2perDOFPlane2    = tfs->make<TH1D>(("HitChi2perDOF_"   +process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;Hit Height #chi^{2}_{}/DOF;Number of Hits"                       ).c_str(), 1000, 0, 100);
      fHitHeightVsWidthPlane2 = tfs->make<TH2D>(("HitHeightVsWidth_"+process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;Hit Width (Ticks);Hit Peak Amplitude (ADC Units)"                ).c_str(), 1000, 0, 100, 1000, 0, 100);
      fHitHeightVsThXYPlane0  = tfs->make<TH2D>(("HitHeightVsThXY_" +process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;#theta^{}_{XY} (radians);Hit Peak Amplitude (ADC Units)").c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitHeightVsThYZPlane0  = tfs->make<TH2D>(("HitHeightVsThYZ_" +process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;#theta^{}_{YZ} (radians);Hit Peak Amplitude (ADC Units)").c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitWidthVsThXYPlane0   = tfs->make<TH2D>(("HitWidthVsThXY_"  +process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;#theta^{}_{XY} (radians);Hit Width (Ticks)"             ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitWidthVsThYZPlane0   = tfs->make<TH2D>(("HitWidthVsThYZ_"  +process+"_"+label+"_1stInd").c_str(), (tpcStr+" Cryostat 1^{st}_{} Induction;#theta^{}_{YZ} (radians);Hit Width (Ticks)"             ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitHeightVsThXYPlane1  = tfs->make<TH2D>(("HitHeightVsThXY_" +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;#theta^{}_{XY} (radians);Hit Peak Amplitude (ADC Units)").c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitHeightVsThYZPlane1  = tfs->make<TH2D>(("HitHeightVsThYZ_" +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;#theta^{}_{YZ} (radians);Hit Peak Amplitude (ADC Units)").c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitWidthVsThXYPlane1   = tfs->make<TH2D>(("HitWidthVsThXY_"  +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;#theta^{}_{XY} (radians);Hit Width (Ticks)"             ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitWidthVsThYZPlane1   = tfs->make<TH2D>(("HitWidthVsThYZ_"  +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" Cryostat 2^{nd}_{} Induction;#theta^{}_{YZ} (radians);Hit Width (Ticks)"             ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitHeightVsThXYPlane2  = tfs->make<TH2D>(("HitHeightVsThXY_" +process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;#theta^{}_{XY} (radians);Hit Peak Amplitude (ADC Units)"         ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitHeightVsThYZPlane2  = tfs->make<TH2D>(("HitHeightVsThYZ_" +process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;#theta^{}_{YZ} (radians);Hit Peak Amplitude (ADC Units)"         ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitWidthVsThXYPlane2   = tfs->make<TH2D>(("HitWidthVsThXY_"  +process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;#theta^{}_{XY} (radians);Hit Width (Ticks)"                      ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
      fHitWidthVsThYZPlane2   = tfs->make<TH2D>(("HitWidthVsThYZ_"  +process+"_"+label+"_Collec").c_str(), (tpcStr+" Cryostat Collection;#theta^{}_{YZ} (radians);Hit Width (Ticks)"                      ).c_str(), 1000, -1*fPi, fPi, 1000, 0, 100);
    } else {
      fHitHeightPlane0        = tfs->make<TH1D>(("HitHeight_"       +process+"_"+label+"_1stInd").c_str(), (tpcStr+" TPC 1^{st}_{} Induction;Hit Peak Amplitude (ADC Units);Number of Hits"   ).c_str(), 1000, 0, 100);
      fHitWidthPlane0         = tfs->make<TH1D>(("HitWidth_"        +process+"_"+label+"_1stInd").c_str(), (tpcStr+" TPC 1^{st}_{} Induction;Hit Width (Ticks);Number of Hits"                ).c_str(), 1000, 0, 100);
      fHitChi2perDOFPlane0    = tfs->make<TH1D>(("HitChi2perDOF_"   +process+"_"+label+"_1stInd").c_str(), (tpcStr+" TPC 1^{st}_{} Induction;Hit Height #chi^{2}_{}/DOF;Number of Hits"       ).c_str(), 1000, 0, 100);
      fHitHeightVsWidthPlane0 = tfs->make<TH2D>(("HitHeightVsWidth_"+process+"_"+label+"_1stInd").c_str(), (tpcStr+" TPC 1^{st}_{} Induction;Hit Width (Ticks);Hit Peak Amplitude (ADC Units)").c_str(), 1000, 0, 100, 1000, 0, 100);
      fHitHeightPlane1        = tfs->make<TH1D>(("HitHeight_"       +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" TPC 2^{nd}_{} Induction;Hit Peak Amplitude (ADC Units);Number of Hits"   ).c_str(), 1000, 0, 100);
      fHitWidthPlane1         = tfs->make<TH1D>(("HitWidth_"        +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" TPC 2^{nd}_{} Induction;Hit Width (Ticks);Number of Hits"                ).c_str(), 1000, 0, 100);
      fHitChi2perDOFPlane1    = tfs->make<TH1D>(("HitChi2perDOF_"   +process+"_"+label+"_2ndInd").c_str(), (tpcStr+" TPC 2^{nd}_{} Induction;Hit Height #chi^{2}_{}/DOF;Number of Hits"       ).c_str(), 1000, 0, 100);
      fHitHeightVsWidthPlane1 = tfs->make<TH2D>(("HitHeightVsWidth_"+process+"_"+label+"_2ndInd").c_str(), (tpcStr+" TPC 2^{nd}_{} Induction;Hit Width (Ticks);Hit Peak Amplitude (ADC Units)").c_str(), 1000, 0, 100, 1000, 0, 100);
      fHitHeightPlane2        = tfs->make<TH1D>(("HitHeight_"       +process+"_"+label+"_Collec").c_str(), (tpcStr+" TPC Collection;Hit Peak Amplitude (ADC Units);Number of Hits"            ).c_str(), 1000, 0, 100);
      fHitWidthPlane2         = tfs->make<TH1D>(("HitWidth_"        +process+"_"+label+"_Collec").c_str(), (tpcStr+" TPC Collection;Hit Width (Ticks);Number of Hits"                         ).c_str(), 1000, 0, 100);
      fHitChi2perDOFPlane2    = tfs->make<TH1D>(("HitChi2perDOF_"   +process+"_"+label+"_Collec").c_str(), (tpcStr+" TPC Collection;Hit Height #chi^{2}_{}/DOF;Number of Hits"                ).c_str(), 1000, 0, 100);
      fHitHeightVsWidthPlane2 = tfs->make<TH2D>(("HitHeightVsWidth_"+process+"_"+label+"_Collec").c_str(), (tpcStr+" TPC Collection;Hit Width (Ticks);Hit Peak Amplitude (ADC Units)"         ).c_str(), 1000, 0, 100, 1000, 0, 100);
    }
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
          geo::PlaneID planeID = hit->WireID().planeID();
          geo::PlaneGeo planeGeo = fGeometry->Plane(planeID);
          double trackDirX = track->DirectionAtPoint(hitInTrack).X();
          double trackDirY = track->DirectionAtPoint(hitInTrack).Y();
          double trackDirZ = track->DirectionAtPoint(hitInTrack).Z();
          double planeCosTh = std::cos(planeGeo.ThetaZ());
          double planeSinTh = std::sin(planeGeo.ThetaZ());
          double relDirX = trackDirX;
          double relDirY = trackDirY * planeCosTh - trackDirZ * planeSinTh;
          double relDirZ = trackDirZ * planeCosTh + trackDirY * planeSinTh;
          double planeThXY = std::atan2(relDirY, relDirX);
          double planeThYZ = std::atan2(relDirY, relDirZ);

          // fill hists
          if        (planeID.Plane == 0)
          {
            fHitHeightPlane0       ->Fill(hit->PeakAmplitude());
            fHitWidthPlane0        ->Fill(hit->RMS());
            fHitChi2perDOFPlane0   ->Fill(hit->GoodnessOfFit());
            fHitHeightVsWidthPlane0->Fill(hit->RMS(), hit->PeakAmplitude());
            fHitHeightVsThXYPlane0 ->Fill(planeThXY, hit->PeakAmplitude());
            fHitHeightVsThYZPlane0 ->Fill(planeThYZ, hit->PeakAmplitude());
            fHitWidthVsThXYPlane0  ->Fill(planeThXY, hit->RMS());
            fHitWidthVsThYZPlane0  ->Fill(planeThYZ, hit->RMS());
          } else if (planeID.Plane == 1) {
            fHitHeightPlane1       ->Fill(hit->PeakAmplitude());
            fHitWidthPlane1        ->Fill(hit->RMS());
            fHitChi2perDOFPlane1   ->Fill(hit->GoodnessOfFit());
            fHitHeightVsWidthPlane1->Fill(hit->RMS(), hit->PeakAmplitude());
            fHitHeightVsThXYPlane1 ->Fill(planeThXY, hit->PeakAmplitude());
            fHitHeightVsThYZPlane1 ->Fill(planeThYZ, hit->PeakAmplitude());
            fHitWidthVsThXYPlane1  ->Fill(planeThXY, hit->RMS());
            fHitWidthVsThYZPlane1  ->Fill(planeThYZ, hit->RMS());
          } else if (planeID.Plane == 2) {
            fHitHeightPlane2       ->Fill(hit->PeakAmplitude());
            fHitWidthPlane2        ->Fill(hit->RMS());
            fHitChi2perDOFPlane2   ->Fill(hit->GoodnessOfFit());
            fHitHeightVsWidthPlane2->Fill(hit->RMS(), hit->PeakAmplitude());
            fHitHeightVsThXYPlane2 ->Fill(planeThXY, hit->PeakAmplitude());
            fHitHeightVsThYZPlane2 ->Fill(planeThYZ, hit->PeakAmplitude());
            fHitWidthVsThXYPlane2  ->Fill(planeThXY, hit->RMS());
            fHitWidthVsThYZPlane2  ->Fill(planeThYZ, hit->RMS());
          }
        }
      }
    } else {
      // get the recob::Hits
      art::Handle<std::vector<recob::Hit>> hitHandle;
      evt.getByLabel(fTag, hitHandle);

      // loop over the wires
      for (auto const& hit : *hitHandle)
      {
        geo::PlaneID planeID = hit.WireID().planeID();
        if        (planeID.Plane == 0)
        {
          fHitHeightPlane0       ->Fill(hit.PeakAmplitude());
          fHitWidthPlane0        ->Fill(hit.RMS());
          fHitChi2perDOFPlane0   ->Fill(hit.GoodnessOfFit());
          fHitHeightVsWidthPlane0->Fill(hit.RMS(), hit.PeakAmplitude());
        } else if (planeID.Plane == 1) {
          fHitHeightPlane1        ->Fill(hit.PeakAmplitude());
          fHitWidthPlane1         ->Fill(hit.RMS());
          fHitChi2perDOFPlane1    ->Fill(hit.GoodnessOfFit());
          fHitHeightVsWidthPlane1->Fill(hit.RMS(), hit.PeakAmplitude());
        } else if (planeID.Plane == 2) {
          fHitHeightPlane2        ->Fill(hit.PeakAmplitude());
          fHitWidthPlane2         ->Fill(hit.RMS());
          fHitChi2perDOFPlane2    ->Fill(hit.GoodnessOfFit());
          fHitHeightVsWidthPlane2->Fill(hit.RMS(), hit.PeakAmplitude());
        }
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
