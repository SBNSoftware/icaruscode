// std inlcudes
#include <cstring>
#include <string>
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

// wiremod
#include "sbncode/WireMod/Utility/WireModUtility.hh"

//namespace
namespace wiremod
{

  class WireModifierXXW : public art::EDProducer
  {
    public:
      explicit WireModifierXXW(fhicl::ParameterSet const& pset);
      void shapeGraphAngle(TGraph2D& graph);
      void shapeGraphPos(TGraph2D& graph, const readout::TPCsetID& tpcsetid);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

    private:
      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
      const geo::WireReadoutGeom* fWireReadout = &(art::ServiceHandle<geo::WireReadout const>()->Get());
      const detinfo::DetectorClocksData fDetClocksData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
      const art::ServiceHandle<cheat::BackTrackerService> fBT;
      std::string fRatioFileName_XXW; // there is where we try to grab the splines/graphs (if they exist)
      std::vector<TGraph2D*> fGraph_charge_XXW;
      std::vector<TGraph2D*> fGraph_sigma_XXW;
      std::string fRatioFileName_YZ;
      std::vector<TGraph2D*> fGraph_charge_YZ;
      std::vector<TGraph2D*> fGraph_sigma_YZ;
      art::InputTag fWireLabel; // which wires are we pulling in?
      art::InputTag fHitLabel;  // which hits are we pulling in?
      art::InputTag fEDepLabel; // which are the EDeps?
      art::InputTag fEDepShiftedLabel; // which are the EDeps?
      unsigned int fCryo; // which Cryo are we in?
      unsigned int fTPCset; // which TPC are we in?
      bool fLocalRatios;           // is the ratio file local?
      bool fSaveHistsByChannel;    // save modified signals by channel?
      bool fSaveHistsByWire;       // save modified signals by wire?
      bool fSaveROITree;           // fill ROITree TTree (before/after waveforms)?
      bool fSaveTickHistograms;    // fill EDep/Hit tick-vs-channel TH2Fs?
      bool fSavePerHitData;        // write per-hit scaleQ/scaleSigma/truthX/thetaXW/isMC art products?
      bool fInRads;                // is the TGraph2D angle axis in radians?
      bool fXAbs;                  // is the TGraph2D x an absolute value?
      bool fAdditiveModification;  // additive (true) vs multiplicative (false) ROI modification
      bool fSetNullScaleIntegral;  // if true, force r_Q=1 for all sub-ROIs (disable integral scaling)
      bool fSetNullScaleWidth;     // if true, force r_sigma=1 for all sub-ROIs (disable width scaling)
      bool fUseChannelROIMode;     // true = read ChannelROI directly (float precision), false = Wire (nominal)
      art::InputTag fChannelLabel; // ChannelROI input label, used only when fUseChannelROIMode is true
      double fOffset;           // ad hoc offset
      TH2F* fEdepTickInd1      = nullptr;
      TH2F* fEdepTickInd2      = nullptr;
      TH2F* fEdepTickColl      = nullptr;
      TH2F* fHitTickInd1       = nullptr;
      TH2F* fHitTickInd2       = nullptr;
      TH2F* fHitTickColl       = nullptr;
      TH2F* fEdepXvsChan       = nullptr;
      TH2F* fHitXvsChan        = nullptr;
      TH2F* fHitBTvsDetProp    = nullptr;

      // diagnostic tree for ROI modification studies
      TTree* fROITree           = nullptr;
      Int_t   tb_run;
      Int_t   tb_subrun;
      Int_t   tb_event;
      Int_t   tb_channel;
      Int_t   tb_wire;
      Int_t   tb_tpc;
      Int_t   tb_plane;
      Int_t   tb_roi_idx;
      Float_t tb_roi_begin;
      Float_t tb_roi_before_total_q;
      Float_t tb_roi_before_center;
      Float_t tb_roi_before_sigma;
      Float_t tb_roi_after_total_q;
      Float_t tb_roi_after_center;
      Float_t tb_roi_after_sigma;
      std::vector<float> tb_roi_before;
      std::vector<float> tb_roi_after;
      std::vector<float> tb_hit_before_center;
      std::vector<float> tb_hit_before_integral;
      std::vector<float> tb_hit_before_sigma;
      std::vector<float> tb_hit_after_center;
      std::vector<float> tb_hit_after_integral;
      std::vector<float> tb_hit_after_sigma;
      std::vector<float> tb_scale_Q;
      std::vector<float> tb_scale_sigma;
      std::vector<float> tb_theta_xw;
      std::vector<float> tb_truth_x;        // truth X coordinate used for XXW graph lookup
      std::vector<int>   tb_track_id;
      // 0=no EDep matched (overlay data hit)
      // 1=energy veto (truth_E<0.3 && Q>80)
      // 2=scaling applied (scale!=1)
      // 3=EDep matched but map evaluates to 1
      std::vector<int>   tb_scale_reason;
      char               tb_xxw_graph_name[128]; // name of the sigma TGraph2D used for this ROI's plane
      std::vector<std::string> fNameVec_sigma_XXW; // kept from reconfigure for tree output
      Int_t tb_roi_n_hits;             // total hits matched to this ROI
      Int_t tb_roi_n_pulse_train_hits; // hits with GoodnessOfFit==-1 && DoF==1
      Int_t tb_roi_is_pulse_train;     // 1 if any hit is a pulse train

  }; // end WireModifierXXW class

  //------------------------------------------------
  WireModifierXXW::WireModifierXXW(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void WireModifierXXW::shapeGraphAngle(TGraph2D& graph)
  {
    // The graphs can be passed in degrees, but we want radians
    // convert the y-axis to radians
    double* xPtr = graph.GetX();
    double* yPtr = graph.GetY();
    double* zPtr = graph.GetZ();
    for (int idx = 0; idx < graph.GetN(); ++idx)
    {
      double x = *(xPtr+idx);
      double y = *(yPtr+idx);
      double z = *(zPtr+idx);
      graph.SetPoint(idx, x, TMath::DegToRad() * y, z);
    }
  }

  //------------------------------------------------
  void WireModifierXXW::shapeGraphPos(TGraph2D& graph, const readout::TPCsetID& tpcsetid)
  {
    // The graphs can be passed in abs(x), but we want global x
    // convert the x-axis to detector x,y,z coordinates
    mf::LogDebug("WireModifierXXW")
      << "Starting graph x range: [" << graph.GetXmin() 
      << ", " << graph.GetXmax() << "]";
    
    // get the cathode possition, the drift direction, and the drift distance from geom
    // get the first TPC Geom since the cathode should be the same in all of them
    geo::TPCID tpcid = fWireReadout->TPCsetToTPCs(tpcsetid).front();
    const geo::TPCGeo* tpcGeom =  fGeometry->TPCPtr(tpcid);
    double cathode = tpcGeom->GetCathodeCenter().X();

    double* xPtr = graph.GetX();
    double* yPtr = graph.GetY();
    double* zPtr = graph.GetZ();
    for (int idx = 0; idx < graph.GetN(); ++idx)
    {
      double x = *(xPtr+idx);
      double y = *(yPtr+idx);
      double z = *(zPtr+idx);
      if (cathode > 0)
      {
        graph.SetPoint(idx, x, y, z);
      } else {
        graph.SetPoint(idx, -x, y, z);
      }
    }

    mf::LogDebug("WireModifierXXW")
      << "Ending graph x range: [" << graph.GetXmin() 
      << ", " << graph.GetXmax() << "]";
  }

  //------------------------------------------------
  void WireModifierXXW::reconfigure(fhicl::ParameterSet const& pset)
  {
    fWireLabel = pset.get<art::InputTag>("WireLabel");
    fHitLabel  = pset.get<art::InputTag>("HitLabel");
    fEDepLabel = pset.get<art::InputTag>("EDepLabel");
    fEDepShiftedLabel = pset.get<art::InputTag>("EDepShiftedLabel");
    fCryo = pset.get<unsigned int>("Cryo");
    fTPCset = pset.get<unsigned int>("TPCset");

    // what, if anything, are we putting in the histogram files
    fSaveHistsByChannel = pset.get<bool>("SaveByChannel",        false);
    fSaveHistsByWire    = pset.get<bool>("SaveByWire",           false);
    fSaveROITree        = pset.get<bool>("SaveROITree",          false);
    fSaveTickHistograms = pset.get<bool>("SaveTickHistograms",   false);
    fSavePerHitData     = pset.get<bool>("SavePerHitData",       false);

    // do we need to reshape the TGraph2D?
    fInRads = pset.get<bool>("InRadians", true);
    fXAbs = pset.get<bool>("XAbs", false);

    // additive vs multiplicative ROI modification
    fAdditiveModification  = pset.get<bool>("AdditiveModification",  false);
    fSetNullScaleIntegral  = pset.get<bool>("SetNullScaleIntegral",  false);
    fSetNullScaleWidth     = pset.get<bool>("SetNullScaleWidth",     false);

    // ChannelROI direct mode (float precision, bypasses Wire intermediate)
    fUseChannelROIMode = pset.get<bool>("UseChannelROIMode", false);
    if (fUseChannelROIMode)
      fChannelLabel = pset.get<art::InputTag>("ChannelLabel");

    // try to read in the graphs/splines from a file
    // if that file does not exist then fake them
    fRatioFileName_XXW = pset.get<std::string>("RatioFileName_XXW", "NOFILE");
    fLocalRatios = pset.get<bool>("LocalRatios", false);
    if (fRatioFileName_XXW == "NOFILE")
    {
      mf::LogVerbatim("WireModifierXXW")
        << "WireModifierXXW::reconfigure - No ratio file given. No scaling is applied...";
    } else if (fRatioFileName_XXW == "DUMMYSCALE")
    {
      mf::LogVerbatim("WireModifierXXW")
        << "WireModifierXXW::reconfigure - Fake universal 2x scaling is applied...";

      // make a graph spanning -400 to 400 in X and -4 to 4 in ThXW (for completeness)
      // all point have a value of 2
      TGraph2D* temp = new TGraph2D();
      for (size_t idx = 0; idx < 4; ++idx)
      {
        double x = (static_cast<double>(idx & 2) - 1)*400;
        double y = (static_cast<double>(2*(idx & 1)) - 1)*4;
        temp->AddPoint(x, y, 2);
      }
      for (size_t idx = 0; idx < 3; ++idx)
      {
        fGraph_charge_XXW.push_back(temp);
        fGraph_sigma_XXW.push_back(temp);
      }
    } else {
      std::string dir_path;
      if (not fLocalRatios)
      {
        char* icaruscode_dir = std::getenv("ICARUSCODE_DIR");
        assert(icaruscode_dir
          && "WireModifierXXW::reconfigure - ICARUSCODE_DIR environment variable must be set!");
        dir_path = std::string(icaruscode_dir) + "/root/WireMod";
      } else
      {
        dir_path = ".";
      }
      mf::LogDebug("WireModifierXXW")
        << "WireModifierXXW::reconfigure - Get file " << fRatioFileName_XXW
        << " from directory " << dir_path;
      TFile* ratioFile = new TFile((dir_path + "/" + fRatioFileName_XXW).c_str(), "READ"); // read only
      assert(ratioFile && "WireModifierXXW::reconfigure - Could not open ratio file");
      // the file exists! pull the ratios
      assert(ratioFile && "WireModifierXXW::reconfigure - WireMod Ratio File Must Exist!");
      assert(!ratioFile->IsZombie()
        && "WireModifierXXW::reconfigure - WireMod Ratio File Must Not be a Zombie!");
      mf::LogVerbatim("WireModifierXXW")
        << "WireModifierXXW::reconfigure - Getting XXW Scales...";
      std::vector<std::string> nameVec_charge_XXW = pset.get<std::vector<std::string>>("XXWScaleIntegral");
      for (auto const& name : nameVec_charge_XXW)
      {
        mf::LogDebug("WireModifierXXW")
          << "WireModifierXXW::reconfigure - Looking for " << name << " in TFile " << fRatioFileName_XXW << "...";
        TGraph2D* temp = static_cast<TGraph2D*>(ratioFile->Get<TGraph2D>(name.c_str())->Clone());
        if (temp != nullptr)
        {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...found";
          if (not fInRads)
            shapeGraphAngle(*temp);
          if (fXAbs)
            shapeGraphPos(*temp, readout::TPCsetID(fCryo, fTPCset));
          fGraph_charge_XXW.push_back(temp);
        } else {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...not found. Don't scale for this plane.";
          // make a graph spanning -400 to 400 in X and -4 to 4 in ThXW (for completeness)
          // all point have a value of 1
          temp = new TGraph2D();
          for (size_t idx = 0; idx < 4; ++idx)
          {
            double x = (static_cast<double>(idx & 2) - 1)*400;
            double y = (static_cast<double>(2*(idx & 1)) - 1)*4;
            temp->AddPoint(x, y, 1);
          }
          fGraph_charge_XXW.push_back(temp);
        }
      }
      std::vector<std::string> nameVec_sigma_XXW = pset.get<std::vector<std::string>>("XXWScaleWidth");
      fNameVec_sigma_XXW = nameVec_sigma_XXW;
      for (auto const& name : nameVec_sigma_XXW)
      {
        mf::LogDebug("WireModifierXXW")
          << "WireModifierXXW::reconfigure - Looking for " << name << " in TFile " << fRatioFileName_XXW << "...";
        TGraph2D* temp = static_cast<TGraph2D*>(ratioFile->Get<TGraph2D>(name.c_str())->Clone());
        if (temp != nullptr)
        {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...found";
          if (not fInRads)
            shapeGraphAngle(*temp);
          if (fXAbs)
            shapeGraphPos(*temp, readout::TPCsetID(fCryo, fTPCset));
          fGraph_sigma_XXW.push_back(temp);
        } else {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...not found. Don't scale for this plane.";
          // make a graph spanning -400 to 400 in X and -4 to 4 in ThXW (for completeness)
          // all point have a value of 1
          temp = new TGraph2D();
          for (size_t idx = 0; idx < 4; ++idx)
          {
            double x = (static_cast<double>(idx & 2) - 1)*400;
            double y = (static_cast<double>(2*(idx & 1)) - 1)*4;
            temp->AddPoint(x, y, 1);
          }
          fGraph_sigma_XXW.push_back(temp);
        }
      }
    }
    fRatioFileName_YZ = pset.get<std::string>("RatioFileName_YZ", "NOFILE");
    fLocalRatios = pset.get<bool>("LocalRatios", false);
    if (fRatioFileName_YZ == "NOFILE")
    {
      mf::LogVerbatim("WireModifierXXW")
        << "WireModifierXXW::reconfigure - No ratio file given. No scaling is applied...";
    } else {
      std::string dir_path;
      if (not fLocalRatios)
      {
        char* icaruscode_dir = std::getenv("ICARUSCODE_DIR");
        assert(icaruscode_dir
          && "WireModifierXXW::reconfigure - ICARUSCODE_DIR environment variable must be set!");
        dir_path = std::string(icaruscode_dir) + "/root/WireMod";
      } else
      {
        dir_path = ".";
      }
      mf::LogDebug("WireModifierXXW")
        << "WireModifierXXW::reconfigure - Get file " << fRatioFileName_YZ
        << " from directory " << dir_path;
      TFile* ratioFile = new TFile((dir_path + "/" + fRatioFileName_YZ).c_str(), "READ"); // read only
      assert(ratioFile && "WireModifierXXW::reconfigure - Could not open ratio file");
      // the file exists! pull the ratios
      assert(ratioFile && "WireModifierXXW::reconfigure - WireMod Ratio File Must Exist!");
      assert(!ratioFile->IsZombie()
        && "WireModifierXXW::reconfigure - WireMod Ratio File Must Not be a Zombie!");
      mf::LogVerbatim("WireModifierXXW")
        << "WireModifierXXW::reconfigure - Getting XXW Scales...";
      std::vector<std::string> nameVec_charge_YZ = pset.get<std::vector<std::string>>("YZScaleIntegral");
      for (auto const& name : nameVec_charge_YZ)
      {
        mf::LogDebug("WireModifierXXW")
          << "WireModifierXXW::reconfigure - Looking for " << name << " in TFile " << fRatioFileName_YZ << "...";
        TGraph2D* temp = static_cast<TGraph2D*>(ratioFile->Get<TGraph2D>(name.c_str())->Clone());
        if (temp != nullptr)
        {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...found";
          fGraph_charge_YZ.push_back(temp);
        } else {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...not found. Don't scale for this plane.";
          // make a graph spanning -400 to 400 in X and -4 to 4 in ThXW (for completeness)
          // all point have a value of 1
          temp = new TGraph2D();
          for (size_t idx = 0; idx < 4; ++idx)
          {
            double x = (static_cast<double>(idx & 2) - 1)*400;
            double y = (static_cast<double>(2*(idx & 1)) - 1)*4;
            temp->AddPoint(x, y, 1);
          }
          fGraph_charge_YZ.push_back(temp);
        }
      }
      std::vector<std::string> nameVec_sigma_YZ = pset.get<std::vector<std::string>>("YZScaleWidth");
      for (auto const& name : nameVec_sigma_YZ)
      {
        mf::LogDebug("WireModifierXXW")
          << "WireModifierXXW::reconfigure - Looking for " << name << " in TFile " << fRatioFileName_YZ << "...";
        TGraph2D* temp = static_cast<TGraph2D*>(ratioFile->Get<TGraph2D>(name.c_str())->Clone());
        if (temp != nullptr)
        {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...found";
          fGraph_sigma_YZ.push_back(temp);
        } else {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...not found. Don't scale for this plane.";
          // make a graph spanning -400 to 400 in X and -4 to 4 in ThXW (for completeness)
          // all point have a value of 1
          temp = new TGraph2D();
          for (size_t idx = 0; idx < 4; ++idx)
          {
            double x = (static_cast<double>(idx & 2) - 1)*400;
            double y = (static_cast<double>(2*(idx & 1)) - 1)*4;
            temp->AddPoint(x, y, 1);
          }
          fGraph_sigma_YZ.push_back(temp);
        }
      }
    }

    // we make these things (mode-dependent)
    if (!fUseChannelROIMode)
      produces<std::vector<recob::Wire>>();
    else
      produces<std::vector<recob::ChannelROI>>();
    // per-hit scale info: one entry per hit in fHitLabel, sentinel -999 if no EDep match
    produces<std::vector<float>>("scaleQ");
    produces<std::vector<float>>("scaleSigma");
    produces<std::vector<float>>("truthX");
    produces<std::vector<float>>("thetaXW");
    // per-hit MC flag: 1 = MC hit (BackTracker matched), 0 = overlay data hit
    produces<std::vector<int>>("isMC");
    // per-hit truth direction and derived geometry (filled only when SavePerHitData: true)
    produces<std::vector<float>>("dirX");   // truth dxdr
    produces<std::vector<float>>("dirY");   // truth dydr
    produces<std::vector<float>>("dirZ");   // truth dzdr
    produces<std::vector<float>>("pitch");  // wirePitch / cosG [cm]
    produces<std::vector<float>>("dQdX");   // integral_before / pitch [ADC/cm], no lifetime correction

    art::ServiceHandle<art::TFileService> tfs;

    if (fSaveTickHistograms)
    {
      fEdepTickInd1 = tfs->make<TH2F>("EdepTickInd1", "1st Induction;Projected Tick;Projected Channel",
                                       4096, -128,  4096,
                                       4603,    0, 55296);
      fEdepTickInd2 = tfs->make<TH2F>("EdepTickInd2", "2nd Induction;Projected Tick;Projected Channel",
                                       4096, -128,  4096,
                                       4603,    0, 55296);
      fEdepTickColl = tfs->make<TH2F>("EdepTickColl", "Collection;Projected Tick;Projected Channel",
                                       4096, -128,  4096,
                                       4603,    0, 55296);
      fHitTickInd1  = tfs->make<TH2F>("HitTickInd1",  "1st Induction;Tick;Channel",
                                       4096, -128,  4096,
                                       4603,    0, 55296);
      fHitTickInd2  = tfs->make<TH2F>("HitTickInd2",  "2nd Induction;Tick;Channel",
                                       4096, -128,  4096,
                                       4603,    0, 55296);
      fHitTickColl  = tfs->make<TH2F>("HitTickColl",  "Collection;Tick;Channel",
                                       4096, -128,  4096,
                                       4603,    0, 55296);
      fEdepXvsChan  = tfs->make<TH2F>("EdepXvsChan", ";X (cm);Channel",
                                        800, -400, 400,
                                       4603,    0, 55296);
      fHitXvsChan   = tfs->make<TH2F>("HitXvsChan", ";X (cm);Channel",
                                        800, -400, 400,
                                       4603,    0, 55296);
      fHitBTvsDetProp = tfs->make<TH2F>("HitBTvsDetProp",
                                        ";X (cm) From BackTracker;X (cm) From DetectorProperties",
                                        800, -400, 400,
                                        800, -400, 400);
      fEdepTickInd1->SetMarkerColor(kRed);
      fEdepTickInd2->SetMarkerColor(kRed);
      fEdepTickColl->SetMarkerColor(kRed);
      fEdepXvsChan->SetMarkerColor(kRed);
      fHitTickInd1->SetMarkerColor(kBlack);
      fHitTickInd2->SetMarkerColor(kBlack);
      fHitTickColl->SetMarkerColor(kBlack);
      fHitXvsChan->SetMarkerColor(kBlack);
    }

    fROITree = tfs->make<TTree>("ROITree", "ROI modification diagnostic");
    if (fSaveROITree)
    {
      fROITree->Branch("run",    &tb_run,    "run/I");
      fROITree->Branch("subrun", &tb_subrun, "subrun/I");
      fROITree->Branch("event",  &tb_event,  "event/I");
      fROITree->Branch("channel",            &tb_channel,            "channel/I");
      fROITree->Branch("wire",               &tb_wire,               "wire/I");
      fROITree->Branch("tpc",                &tb_tpc,                "tpc/I");
      fROITree->Branch("plane",              &tb_plane,              "plane/I");
      fROITree->Branch("roi_idx",            &tb_roi_idx,            "roi_idx/I");
      fROITree->Branch("roi_begin",          &tb_roi_begin,          "roi_begin/F");
      fROITree->Branch("roi_before_total_q", &tb_roi_before_total_q, "roi_before_total_q/F");
      fROITree->Branch("roi_before_center",  &tb_roi_before_center,  "roi_before_center/F");
      fROITree->Branch("roi_before_sigma",   &tb_roi_before_sigma,   "roi_before_sigma/F");
      fROITree->Branch("roi_after_total_q",  &tb_roi_after_total_q,  "roi_after_total_q/F");
      fROITree->Branch("roi_after_center",   &tb_roi_after_center,   "roi_after_center/F");
      fROITree->Branch("roi_after_sigma",    &tb_roi_after_sigma,    "roi_after_sigma/F");
      fROITree->Branch("roi_before",         &tb_roi_before);
      fROITree->Branch("roi_after",          &tb_roi_after);
      fROITree->Branch("hit_before_center",   &tb_hit_before_center);
      fROITree->Branch("hit_before_integral", &tb_hit_before_integral);
      fROITree->Branch("hit_before_sigma",    &tb_hit_before_sigma);
      fROITree->Branch("hit_after_center",    &tb_hit_after_center);
      fROITree->Branch("hit_after_integral",  &tb_hit_after_integral);
      fROITree->Branch("hit_after_sigma",     &tb_hit_after_sigma);
      fROITree->Branch("scale_Q",            &tb_scale_Q);
      fROITree->Branch("scale_sigma",        &tb_scale_sigma);
      fROITree->Branch("theta_xw",           &tb_theta_xw);
      fROITree->Branch("truth_x",            &tb_truth_x);
      fROITree->Branch("track_id",           &tb_track_id);
      fROITree->Branch("scale_reason",       &tb_scale_reason);
      fROITree->Branch("xxw_graph_name",     tb_xxw_graph_name, "xxw_graph_name/C");
      fROITree->Branch("roi_n_hits",             &tb_roi_n_hits,             "roi_n_hits/I");
      fROITree->Branch("roi_n_pulse_train_hits", &tb_roi_n_pulse_train_hits, "roi_n_pulse_train_hits/I");
      fROITree->Branch("roi_is_pulse_train",     &tb_roi_is_pulse_train,     "roi_is_pulse_train/I");
    }
  }

  //------------------------------------------------
  void WireModifierXXW::produce(art::Event& evt)
  {
    // here's where the "magic" happens
    art::ServiceHandle<art::TFileService> tfs;

    // get a clock and det props for the event
    const detinfo::DetectorPropertiesData detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    // handles, read conditionally based on mode
    art::Handle< std::vector<recob::Wire> >       wireHandle;
    art::Handle< std::vector<recob::ChannelROI> > chanROIHandle;
    if (!fUseChannelROIMode)
      evt.getByLabel(fWireLabel, wireHandle);
    else
      evt.getByLabel(fChannelLabel, chanROIHandle);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepHandle;
    evt.getByLabel(fEDepLabel, edepHandle);
    auto const& edepVec(*edepHandle);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepShiftedHandle;
    evt.getByLabel(fEDepShiftedLabel, edepShiftedHandle);
    auto const& edepShiftedVec(*edepShiftedHandle);

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitLabel, hitHandle);
    auto const& hitVec(*hitHandle);

    // per-hit scale products, only sized when SavePerHitData is on
    const size_t perHitN = fSavePerHitData ? hitVec.size() : 0;
    auto vec_scale_q     = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_scale_sigma = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_truth_x     = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_theta_xw    = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_is_mc       = std::make_unique<std::vector<int>>  (perHitN, 0);
    auto vec_dirx        = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_diry        = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_dirz        = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_pitch       = std::make_unique<std::vector<float>>(perHitN, -999.f);
    auto vec_dqdx        = std::make_unique<std::vector<float>>(perHitN, -999.f);

    // Get the tick offset by comparing the hit ticks with the tick gotten by
    // projecting the backtracked hit X to a tick. Yes this is silly thanks for asking.
    // Also fills vec_is_mc: 1 if BackTracker finds a match (MC hit), 0 otherwise (overlay).
    double BT_Offset = 0;
    unsigned int BT_Hits = 0;
    for (size_t i_h = 0; i_h < hitVec.size(); ++i_h)
    {
      auto const& hit = hitVec[i_h];
      double hitX = std::numeric_limits<double>::max();
      try
      {
        hitX = fBT->HitToXYZ(fDetClocksData, hit).at(0);
      }
      catch (...)
      {
        // Overlay hit
        continue;
      }
      if (fSavePerHitData) (*vec_is_mc)[i_h] = 1;
      double hitTick  = hit.PeakTime();
      double projTick = detProp.ConvertXToTicks(hitX, hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat);
      BT_Offset += (hitTick - projTick);
      ++BT_Hits;
    }
    if (BT_Hits > 0) BT_Offset /= BT_Hits;

    sys::WireModUtility wmUtil(fGeometry, fWireReadout, detProp,
                               false, // Channel Scale
                               false, // X
                               (fRatioFileName_YZ != "NOFILE"), // YZ
                               false, // XZ-Angle
                               false, // YZ-Angle
                               false, // dE/dx
                               (fRatioFileName_XXW != "NOFILE"),  // X-ThXW
                               BT_Offset); // Tick Offset
    if (fRatioFileName_YZ != "NOFILE")
    {
      wmUtil.graph2Ds_Charge_YZ = fGraph_charge_YZ;
      wmUtil.graph2Ds_Sigma_YZ  = fGraph_sigma_YZ;
    }
    if (fRatioFileName_XXW != "NOFILE")
    {
      wmUtil.graph2Ds_Charge_XXW = fGraph_charge_XXW;
      wmUtil.graph2Ds_Sigma_XXW  = fGraph_sigma_XXW;
    }

    wmUtil.additiveModification = fAdditiveModification;

    // add some debugging here
    mf::LogVerbatim("WireModifierXXW")
      << "DUMP CONFIG:" << '\n'
      << "---------------------------------------------------" << '\n'
      << "  applyChannelScale:     " << wmUtil.applyChannelScale     << '\n'
      << "  applyXScale:           " << wmUtil.applyXScale           << '\n'
      << "  applyYZScale:          " << wmUtil.applyYZScale          << '\n'
      << "  applyXXWScale:         " << wmUtil.applyXXWScale         << '\n'
      << "  applyXZAngleScale:     " << wmUtil.applyXZAngleScale     << '\n'
      << "  applyYZAngleScale:     " << wmUtil.applyYZAngleScale     << '\n'
      << "  applydEdXScale:        " << wmUtil.applydEdXScale        << '\n'
      << "  additiveModification:  " << wmUtil.additiveModification  << '\n'
      << "  setNullScaleIntegral:  " << fSetNullScaleIntegral        << '\n'
      << "  setNullScaleWidth:     " << fSetNullScaleWidth           << '\n'
      << "  readoutWindowTicks:    " << wmUtil.readoutWindowTicks    << '\n'
      << "  tickOffset:            " << wmUtil.tickOffset            << '\n'
      << "---------------------------------------------------";

    // tick check
    mf::LogDebug("WireModifierXXW")
      << "edepVec.size() = " << edepVec.size() << '\n'
      << "edepShiftedVec.size() = " << edepShiftedVec.size() << '\n'
      << "hitVec.size() = " << hitVec.size();
    if (fSaveTickHistograms)
    {
      for (auto const& edep : edepVec)
      {
        if (fWireReadout->FindTPCsetAtPosition(edep.MidPoint()) != readout::TPCsetID(fCryo, fTPCset))
          continue;

        geo::TPCGeo const* curTPCGeomPtr = fGeometry->PositionToTPCptr(edep.MidPoint());
        for (auto const& plane : fWireReadout->Iterate<geo::PlaneGeo>(curTPCGeomPtr->ID()))
        {
          TH2F* targetHist = (plane.View() == geo::kY) ? fEdepTickInd1
                           : (plane.View() == geo::kV) ? fEdepTickInd2
                           : (plane.View() == geo::kU) ? fEdepTickColl
                           :                             nullptr;
          if (targetHist == nullptr) continue;

          geo::WireID wireID;
          try { wireID = plane.NearestWireID(edep.MidPoint()); }
          catch (...) { continue; }
          float projTick = detProp.ConvertXToTicks(edep.X(), wireID.Plane, wireID.TPC, wireID.Cryostat) + wmUtil.tickOffset;
          float projChan = fWireReadout->PlaneWireToChannel(wireID);
          targetHist->Fill(projTick, projChan);
          fEdepXvsChan->Fill(edep.X(), projChan);
        }
      }
      for (auto const& hit : hitVec)
      {
        TH2F* targetHist = (hit.View() == geo::kY) ? fHitTickInd1
                         : (hit.View() == geo::kV) ? fHitTickInd2
                         : (hit.View() == geo::kU) ? fHitTickColl
                         :                           nullptr;
        if (targetHist == nullptr) continue;

        float tick = hit.PeakTime();
        float chan = hit.Channel();
        targetHist->Fill(tick, chan);
        float hitX = detProp.ConvertTicksToX(tick,
                                             hit.WireID().Plane, hit.WireID().TPC, hit.WireID().Cryostat);
        try
        {
          float hitXBT = fBT->HitToXYZ(fDetClocksData, hit).at(0);
          fHitBTvsDetProp->Fill(hitXBT, hitX);
          hitX = hitXBT;
        }
        catch (...) { }
        fHitXvsChan->Fill(hitX, chan);
      }
    }

    double offset_ADC = 0;

    if (!fUseChannelROIMode) {
    // Wire (nominal) mode
    auto const& wireVec(*wireHandle);
    std::unique_ptr<std::vector<recob::Wire>> new_wires(new std::vector<recob::Wire>());

    mf::LogVerbatim("WireModifierXXW") << "Get Edep Map";
    wmUtil.FillROIMatchedEdepMap(edepShiftedVec, wireVec, offset_ADC);
    mf::LogVerbatim("WireModifierXXW") << "Got Edep Map." << '\n' << "Get Hit Map";
    wmUtil.FillROIMatchedHitMap(hitVec, wireVec);
    mf::LogVerbatim("WireModifierXXW") << "Got Hit Map.";

    tb_run    = static_cast<Int_t>(evt.run());
    tb_subrun = static_cast<Int_t>(evt.subRun());
    tb_event  = static_cast<Int_t>(evt.event());

    // loop-de-loop
    for(size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    {
      mf::LogDebug("WireModifierXXW")
        << "Checking wire " << i_w;

      auto const& wire = wireVec.at(i_w);
      recob::Wire::RegionsOfInterest_t new_rois;
      bool isModified = false;

      if (wire.NSignal() == 0)
        continue;

      new_rois.resize(wire.NSignal());

      unsigned int my_plane = geo::kUnknown;
      if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 0)).View())
      {
        mf::LogDebug("WireModifierXXW")
          << "Wire is on plane 0, view " << wire.View();
        my_plane = 0;
      } else if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 1)).View()) {
        mf::LogDebug("WireModifierXXW")
          << "Wire is on plane 1, view " << wire.View();
        my_plane = 1;
      } else if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 2)).View()) {
        mf::LogDebug("WireModifierXXW")
          << "Wire is on plane 2, view " << wire.View();
        my_plane = 2;
      }

      if (my_plane == geo::kUnknown)
      {
        mf::LogDebug("WireModifierXXW")
          << "Wire is on unsupported plane. Skip.";
      }

      for(size_t i_r = 0; i_r < wire.SignalROI().get_ranges().size(); ++i_r)
      {
        mf::LogDebug("WireModifierXXW")
          << "  Checking ROI " << i_r;

        auto const& range = wire.SignalROI().get_ranges()[i_r];

        sys::WireModUtility::ROI_Key_t roi_key(wire.Channel(), i_r);

        std::vector<float> modified_data(range.data());

        auto it_map = wmUtil.ROIMatchedEdepMap.find(roi_key);
        if(it_map==wmUtil.ROIMatchedEdepMap.end()){
          new_rois.add_range(range.begin_index(), modified_data);
          mf::LogDebug("WireModifierXXW")
            << "    Could not find matching Edep. Skip";
          continue;
        }
        std::vector<size_t> matchedEdepIdxVec = it_map->second;
        if(matchedEdepIdxVec.size() == 0)
        {
          new_rois.add_range(range.begin_index(), modified_data);
          mf::LogDebug("WireModifierXXW")
            << "    No indices for Edep. Skip";
          continue;
        }
        std::vector<const sim::SimEnergyDeposit*> matchedEdepPtrVec;
        std::vector<const sim::SimEnergyDeposit*> matchedEdepShiftedPtrVec;
        for(auto i_e : matchedEdepIdxVec)
        {
          matchedEdepPtrVec.push_back(&edepVec[i_e]);
          matchedEdepShiftedPtrVec.push_back(&edepShiftedVec[i_e]);
        }
        mf::LogDebug("WireModifierXXW")
          << "  Found " << matchedEdepPtrVec.size() << " Edeps";

        std::vector<const recob::Hit*> matchedHitPtrVec;
        auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ) {
          for( auto i_h : it_hit_map->second )
            matchedHitPtrVec.push_back(&hitVec[i_h]);
        }

        mf::LogDebug("WireModifierXXW")
          << "    Found " << matchedHitPtrVec.size() << " matching hits";

        auto roi_properties = wmUtil.CalcROIProperties(wire, i_r);
        mf::LogDebug("WireModifierXXW")
          << "    ROI Properties:" << '\n'
          << "                    key:     (" << roi_properties.key.first
                                      << ", " << roi_properties.key.second << ")" << '\n'
          << "                    view:    " << roi_properties.view << '\n'
          << "                    begin:   " << roi_properties.begin << '\n'
          << "                    end:     " << roi_properties.end << '\n'
          << "                    total_q: " << roi_properties.total_q << '\n'
          << "                    center:  " << roi_properties.center << '\n'
          << "                    sigma:   " << roi_properties.sigma;

        auto subROIPropVec = wmUtil.CalcSubROIProperties(roi_properties, matchedHitPtrVec);
        
        mf::LogDebug("WireModifierXXW")
          << "    have " << subROIPropVec.size() << " SubROI";

        auto SubROIMatchedEdepShiftedMap =
          wmUtil.MatchEdepsToSubROIs(subROIPropVec, matchedEdepShiftedPtrVec, offset_ADC);
        std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>>
          SubROIMatchedEdepMap;
        for (auto const& key_edepPtrVec_pair : SubROIMatchedEdepShiftedMap)
        {
          auto key = key_edepPtrVec_pair.first;
          for (auto const& shifted_edep_ptr : key_edepPtrVec_pair.second)
          {
            for (size_t i_e = 0; i_e < matchedEdepShiftedPtrVec.size(); ++i_e)
            {
              if (shifted_edep_ptr == matchedEdepShiftedPtrVec[i_e])
              {
                SubROIMatchedEdepMap[key].push_back(matchedEdepPtrVec[i_e]);
                break;
              }
            }
          }
        }

        std::map<sys::WireModUtility::SubROI_Key_t,
                 sys::WireModUtility::ScaleValues_t> SubROIMatchedScalesMap;
        std::map<sys::WireModUtility::SubROI_Key_t,
                 sys::WireModUtility::TruthProperties_t> SubROIMatchedTruthMap;
        std::map<sys::WireModUtility::SubROI_Key_t, int> SubROIScaleReasonMap;
        for (auto const& subroi_prop : subROIPropVec)
        {
          sys::WireModUtility::ScaleValues_t scale_vals;
          auto key = subroi_prop.key;
          auto key_it =  SubROIMatchedEdepMap.find(key);
          int scale_reason = 0;

          if (key_it != SubROIMatchedEdepMap.end() && key_it->second.size() > 0)
          {
            auto truth_vals = wmUtil.CalcPropertiesFromEdeps(key_it->second, offset_ADC);
            SubROIMatchedTruthMap[key] = truth_vals;

            if (truth_vals.total_energy < 0.3 && subroi_prop.total_q > 80)
            {
              scale_vals.r_Q     = 1.;
              scale_vals.r_sigma = 1.;
              scale_reason = 1; // energy veto
            } else
            {
              scale_vals = wmUtil.GetScaleValues(truth_vals, roi_properties);
              if (scale_vals.r_Q != 1. || scale_vals.r_sigma != 1.)
              {
                mf::LogDebug("WireModifierXXW")
                  << "Scaling! Q scale: " << scale_vals.r_Q
                  << "     sigma scale: " << scale_vals.r_sigma;
                isModified = true;
                scale_reason = 2; // scaling applied
              } else {
                scale_reason = 3; // EDep matched but map evaluates to 1
              }
            }
          } else
          {
            scale_vals.r_Q     = 1.;
            scale_vals.r_sigma = 1.;
            scale_reason = 0; // no EDep matched (likely overlay data hit)
          }

          SubROIMatchedScalesMap[key]  = scale_vals;
          SubROIScaleReasonMap[key]    = scale_reason;
        }
        
        // null-scale overrides: force r_Q and/or r_sigma to 1 for all sub-ROIs.
        // applied in-place so the diagnostic tree records what was actually used.
        if (fSetNullScaleIntegral || fSetNullScaleWidth) {
          for (auto& kv : SubROIMatchedScalesMap) {
            if (fSetNullScaleIntegral) kv.second.r_Q     = 1.0;
            if (fSetNullScaleWidth)    kv.second.r_sigma  = 1.0;
          }
        }

        wmUtil.ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap, 1.5);

        // fill diagnostic tree
        if (fSaveROITree)
        {
          // wire/geometry info
          std::vector<geo::WireID> wireIDs = fWireReadout->ChannelToWire(wire.Channel());
          geo::WireID firstWireID = wireIDs.empty() ? geo::WireID{} : wireIDs.front();
          tb_channel = static_cast<Int_t>(wire.Channel());
          tb_wire    = static_cast<Int_t>(firstWireID.Wire);
          tb_tpc     = static_cast<Int_t>(firstWireID.TPC);
          tb_plane   = static_cast<Int_t>(my_plane);
          tb_roi_idx = static_cast<Int_t>(i_r);

          // ROI before
          tb_roi_begin         = roi_properties.begin;
          tb_roi_before_total_q = roi_properties.total_q;
          tb_roi_before_center  = roi_properties.center;
          tb_roi_before_sigma   = roi_properties.sigma;
          tb_roi_before.assign(range.data().begin(), range.data().end());

          // ROI after: compute charge-weighted properties from modified waveform
          float mod_total_q = 0.f, mod_center = 0.f, mod_sigma = 0.f;
          for (size_t i_t = 0; i_t < modified_data.size(); ++i_t)
          {
            mod_total_q += modified_data[i_t];
            mod_center  += modified_data[i_t] * (i_t + roi_properties.begin);
          }
          if (mod_total_q > 0.f) mod_center /= mod_total_q;
          for (size_t i_t = 0; i_t < modified_data.size(); ++i_t)
          {
            float t = static_cast<float>(i_t) + roi_properties.begin;
            mod_sigma += modified_data[i_t] * (t - mod_center) * (t - mod_center);
          }
          if (mod_total_q > 0.f) mod_sigma = std::sqrt(mod_sigma / mod_total_q);
          tb_roi_after_total_q = mod_total_q;
          tb_roi_after_center  = mod_center;
          tb_roi_after_sigma   = mod_sigma;
          tb_roi_after = modified_data;

          // pulse train detection (GoodnessOfFit==-1 && DoF==1 is gaushit's sentinel)
          {
            int n_pt = 0;
            for (auto const& hptr : matchedHitPtrVec)
              if (hptr->GoodnessOfFit() == -1.0f && hptr->DegreesOfFreedom() == 1)
                ++n_pt;
            tb_roi_n_hits             = static_cast<Int_t>(matchedHitPtrVec.size());
            tb_roi_n_pulse_train_hits = static_cast<Int_t>(n_pt);
            tb_roi_is_pulse_train     = (n_pt > 0) ? 1 : 0;
          }

          // per-subROI hit candidate info and scale factors
          tb_hit_before_center.clear();
          tb_hit_before_integral.clear();
          tb_hit_before_sigma.clear();
          tb_hit_after_center.clear();
          tb_hit_after_integral.clear();
          tb_hit_after_sigma.clear();
          tb_scale_Q.clear();
          tb_scale_sigma.clear();
          tb_theta_xw.clear();
          tb_truth_x.clear();
          tb_track_id.clear();
          tb_scale_reason.clear();
          {
            std::string gname = (tb_plane < (int)fNameVec_sigma_XXW.size())
                                ? fNameVec_sigma_XXW[tb_plane] : "unknown";
            std::strncpy(tb_xxw_graph_name, gname.c_str(), 127);
            tb_xxw_graph_name[127] = '\0';
          }

          for (auto const& subroi_prop : subROIPropVec)
          {
            auto key        = subroi_prop.key;
            auto scale_vals = SubROIMatchedScalesMap.at(key);

            tb_hit_before_center.push_back(subroi_prop.center);
            tb_hit_before_integral.push_back(subroi_prop.total_q);
            tb_hit_before_sigma.push_back(subroi_prop.sigma);
            tb_hit_after_center.push_back(subroi_prop.center);
            tb_hit_after_integral.push_back(subroi_prop.total_q * static_cast<float>(scale_vals.r_Q));
            tb_hit_after_sigma.push_back(subroi_prop.sigma    * static_cast<float>(scale_vals.r_sigma));
            tb_scale_Q.push_back(static_cast<float>(scale_vals.r_Q));
            tb_scale_sigma.push_back(static_cast<float>(scale_vals.r_sigma));

            // dominant track ID: pick the TrackID with highest summed energy
            int dominant_tid = -999;
            auto edep_it = SubROIMatchedEdepMap.find(key);
            if (edep_it != SubROIMatchedEdepMap.end() && !edep_it->second.empty())
            {
              std::map<int, double> tidEnergy;
              for (auto const& ep : edep_it->second)
                tidEnergy[ep->TrackID()] += ep->E();
              double maxE = 0.;
              for (auto const& te : tidEnergy)
                if (te.second > maxE) { maxE = te.second; dominant_tid = te.first; }
            }
            tb_track_id.push_back(dominant_tid);
            tb_scale_reason.push_back(SubROIScaleReasonMap.at(key));

            // angle in degrees (ThetaXW returns radians) + truth X for graph lookup
            float theta_deg = -999.f;
            float truth_x_val = -9999.f;
            auto truth_it = SubROIMatchedTruthMap.find(key);
            if (truth_it != SubROIMatchedTruthMap.end())
            {
              auto const& tv = truth_it->second;
              truth_x_val = static_cast<float>(tv.x);
              geo::TPCGeo const* tpcPtr = fGeometry->PositionToTPCptr(
                  geo::Point_t{(double)tv.x, (double)tv.y, (double)tv.z});
              if (tpcPtr != nullptr)
              {
                try
                {
                  auto const& plane_obj = fWireReadout->Plane(tpcPtr->ID(), wire.View());
                  theta_deg = static_cast<float>(wmUtil.ThetaXW(
                      tv.dxdr, tv.dydr, tv.dzdr, plane_obj.ThetaZ()) * TMath::RadToDeg());
                }
                catch (...) {}
              }
            }
            tb_theta_xw.push_back(theta_deg);
            tb_truth_x.push_back(truth_x_val);
          }

          fROITree->Fill();
        }

        // Fill per-hit scale product vectors indexed by position in hitVec
        if (fSavePerHitData && it_hit_map != wmUtil.ROIMatchedHitMap.end())
        {
          for (auto const& subroi_prop : subROIPropVec)
          {
            size_t i_subroi = subroi_prop.key.second;
            if (i_subroi >= it_hit_map->second.size()) continue;
            size_t hitVec_idx = it_hit_map->second[i_subroi];

            auto const& scale_vals = SubROIMatchedScalesMap.at(subroi_prop.key);
            (*vec_scale_q    )[hitVec_idx] = static_cast<float>(scale_vals.r_Q);
            (*vec_scale_sigma)[hitVec_idx] = static_cast<float>(scale_vals.r_sigma);

            auto truth_it = SubROIMatchedTruthMap.find(subroi_prop.key);
            if (truth_it != SubROIMatchedTruthMap.end())
            {
              auto const& tv = truth_it->second;
              (*vec_truth_x)[hitVec_idx] = static_cast<float>(tv.x);
              geo::TPCGeo const* tpcPtr = fGeometry->PositionToTPCptr(
                  geo::Point_t{tv.x, tv.y, tv.z});
              if (tpcPtr)
              {
                try {
                  auto const& plane_obj = fWireReadout->Plane(
                      hitVec[hitVec_idx].WireID().planeID());
                  double thetaZ = plane_obj.ThetaZ();
                  (*vec_theta_xw)[hitVec_idx] = static_cast<float>(
                      wmUtil.ThetaXW(tv.dxdr, tv.dydr, tv.dzdr, thetaZ) * TMath::RadToDeg());
                  (*vec_dirx)[hitVec_idx] = static_cast<float>(tv.dxdr);
                  (*vec_diry)[hitVec_idx] = static_cast<float>(tv.dydr);
                  (*vec_dirz)[hitVec_idx] = static_cast<float>(tv.dzdr);
                  double cosG = std::abs(tv.dydr * std::sin(thetaZ) + tv.dzdr * std::cos(thetaZ));
                  float wp    = static_cast<float>(plane_obj.WirePitch());
                  float pval  = (cosG > 1e-6) ? wp / static_cast<float>(cosG) : -999.f;
                  (*vec_pitch)[hitVec_idx] = pval;
                  (*vec_dqdx)[hitVec_idx]  = (pval > 0.f)
                                             ? hitVec[hitVec_idx].Integral() / pval : -999.f;
                }
                catch (...) {}
              }
            }
          }
        }

        new_rois.add_range(roi_properties.begin, modified_data);
      }

      new_wires->emplace_back(new_rois, wire.Channel(), wire.View());

      if (fSaveHistsByChannel && isModified)
      {
        readout::ROPID ropID = fWireReadout->ChannelToROP(wire.Channel());  
        std::string titleStr =  "Cryo-"         + std::to_string(ropID.Cryostat)
                             + "_TPCset-"       + std::to_string(ropID.TPCset)
                             + "_ReadOutPlane-" + std::to_string(ropID.ROP)
                             + "_Channel-"      + std::to_string(wire.Channel());
        TH1F* oldChannelHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                        wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        TH1F* newChannelHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                        wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
        {
          float oldSample =
            (tick <     wire         .Signal().size() ) ?     wire         .Signal().at(tick) : 0;
          float newSample =
            (tick < new_wires->back().Signal().size() ) ? new_wires->back().Signal().at(tick) : 0;
          oldChannelHist->SetBinContent(tick + 1, oldSample);
          newChannelHist->SetBinContent(tick + 1, newSample);
        }
          
        TH1F* oldChannelHist_toSave = tfs->make<TH1F>(*oldChannelHist);
        TH1F* newChannelHist_toSave = tfs->make<TH1F>(*newChannelHist);
        mf::LogDebug("WireModifierXXW")
          << "Saved histograms " << oldChannelHist_toSave->GetName() << '\n'
          << "             and " << newChannelHist_toSave->GetName();
      }

      if (fSaveHistsByWire && isModified)
      {
        std::vector<geo::WireID> wireIDs = fWireReadout->ChannelToWire(wire.Channel());
        mf::LogDebug("WireModifierXXW")
          << "Channel " << wire.Channel() << " has " << wireIDs.size() << " wire(s)";
        for (auto const& wireID : wireIDs)
        {
          std::string titleStr =  "Cryo-"  + std::to_string(wireID.Cryostat)
                               + "_TPC-"   + std::to_string(wireID.TPC)
                               + "_Plane-" + std::to_string(wireID.Plane)
                               + "_Wire-"  + std::to_string(wireID.Wire);
          TH1F* oldWireHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                       wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          TH1F* newWireHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                       wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
          {
            float oldSample =
              (tick <     wire         .Signal().size() ) ?     wire         .Signal().at(tick) : 0;
            float newSample =
              (tick < new_wires->back().Signal().size() ) ? new_wires->back().Signal().at(tick) : 0;
            oldWireHist->SetBinContent(tick + 1, oldSample);
            newWireHist->SetBinContent(tick + 1, newSample);
          }

          TH1F* oldWireHist_toSave = tfs->make<TH1F>(*oldWireHist);
          TH1F* newWireHist_toSave = tfs->make<TH1F>(*newWireHist);
          mf::LogDebug("WireModifierXXW")
            << "Saved histograms " << oldWireHist_toSave->GetName() << '\n'
            << "             and " << newWireHist_toSave->GetName();
        }
      }

    } // end loop over wires

    evt.put(std::move(new_wires));

    } else {
    // ChannelROI (float precision) mode
    auto const& chanROIVec(*chanROIHandle);
    std::unique_ptr<std::vector<recob::ChannelROI>> new_chanROIs(new std::vector<recob::ChannelROI>());

    mf::LogVerbatim("WireModifierXXW") << "Get Edep Map";
    wmUtil.FillROIMatchedEdepMap(edepShiftedVec, chanROIVec, offset_ADC);
    mf::LogVerbatim("WireModifierXXW") << "Got Edep Map." << '\n' << "Get Hit Map";
    wmUtil.FillROIMatchedHitMap(hitVec, chanROIVec);
    mf::LogVerbatim("WireModifierXXW") << "Got Hit Map.";

    tb_run    = static_cast<Int_t>(evt.run());
    tb_subrun = static_cast<Int_t>(evt.subRun());
    tb_event  = static_cast<Int_t>(evt.event());

    for(size_t i_c = 0; i_c < chanROIVec.size(); ++i_c)
    {
      mf::LogDebug("WireModifierXXW") << "Checking chanROI " << i_c;

      auto const& chanROI = chanROIVec.at(i_c);
      recob::ChannelROI::RegionsOfInterest_t new_rois;
      bool isModified = false;

      if (chanROI.NSignal() == 0)
        continue;

      new_rois.resize(chanROI.NSignal());

      geo::View_t chanView = fWireReadout->View(chanROI.Channel());
      unsigned int my_plane = geo::kUnknown;
      if      (chanView == fWireReadout->Plane(geo::PlaneID(0, 0, 0)).View()) my_plane = 0;
      else if (chanView == fWireReadout->Plane(geo::PlaneID(0, 0, 1)).View()) my_plane = 1;
      else if (chanView == fWireReadout->Plane(geo::PlaneID(0, 0, 2)).View()) my_plane = 2;

      if (my_plane == geo::kUnknown)
        mf::LogDebug("WireModifierXXW") << "ChannelROI on unsupported plane. Skip.";

      // SignalROIF() returns by value, so a const& to its sub-objects would dangle.
      // SignalROI() is a stable const ref instead.
      auto const& roiShort = chanROI.SignalROI();
      const float adcScale = 1.0f / static_cast<float>(chanROI.ADCScaleFactor());

      for(size_t i_r = 0; i_r < roiShort.get_ranges().size(); ++i_r)
      {
        mf::LogDebug("WireModifierXXW") << "  Checking ROI " << i_r;

        auto const& range_short = roiShort.get_ranges()[i_r];
        sys::WireModUtility::ROI_Key_t roi_key(chanROI.Channel(), i_r);

        // convert short -> float (same as SignalROIF() but without the temporary)
        std::vector<float> modified_data;
        modified_data.reserve(range_short.data().size());
        for (auto s : range_short.data())
          modified_data.push_back(static_cast<float>(s) * adcScale);

        auto toShort = [&chanROI](float v) {
          return static_cast<short>(std::round(v * chanROI.ADCScaleFactor()));
        };

        auto it_map = wmUtil.ROIMatchedEdepMap.find(roi_key);
        if(it_map==wmUtil.ROIMatchedEdepMap.end()){
          std::vector<short> pass_short(modified_data.size());
          std::transform(modified_data.begin(), modified_data.end(), pass_short.begin(), toShort);
          new_rois.add_range(range_short.begin_index(), pass_short);
          mf::LogDebug("WireModifierXXW") << "    Could not find matching Edep. Skip";
          continue;
        }
        std::vector<size_t> matchedEdepIdxVec = it_map->second;
        if(matchedEdepIdxVec.size() == 0)
        {
          std::vector<short> pass_short(modified_data.size());
          std::transform(modified_data.begin(), modified_data.end(), pass_short.begin(), toShort);
          new_rois.add_range(range_short.begin_index(), pass_short);
          mf::LogDebug("WireModifierXXW") << "    No indices for Edep. Skip";
          continue;
        }
        std::vector<const sim::SimEnergyDeposit*> matchedEdepPtrVec;
        std::vector<const sim::SimEnergyDeposit*> matchedEdepShiftedPtrVec;
        for(auto i_e : matchedEdepIdxVec)
        {
          matchedEdepPtrVec.push_back(&edepVec[i_e]);
          matchedEdepShiftedPtrVec.push_back(&edepShiftedVec[i_e]);
        }
        mf::LogDebug("WireModifierXXW") << "  Found " << matchedEdepPtrVec.size() << " Edeps";

        std::vector<const recob::Hit*> matchedHitPtrVec;
        auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ) {
          for( auto i_h : it_hit_map->second )
            matchedHitPtrVec.push_back(&hitVec[i_h]);
        }
        mf::LogDebug("WireModifierXXW") << "    Found " << matchedHitPtrVec.size() << " matching hits";

        auto roi_properties = wmUtil.CalcROIProperties(chanROI, i_r);
        mf::LogDebug("WireModifierXXW")
          << "    ROI Properties:" << '\n'
          << "                    key:     (" << roi_properties.key.first
                                      << ", " << roi_properties.key.second << ")" << '\n'
          << "                    view:    " << roi_properties.view << '\n'
          << "                    begin:   " << roi_properties.begin << '\n'
          << "                    end:     " << roi_properties.end << '\n'
          << "                    total_q: " << roi_properties.total_q << '\n'
          << "                    center:  " << roi_properties.center << '\n'
          << "                    sigma:   " << roi_properties.sigma;

        auto subROIPropVec = wmUtil.CalcSubROIProperties(roi_properties, matchedHitPtrVec);
        mf::LogDebug("WireModifierXXW") << "    have " << subROIPropVec.size() << " SubROI";

        auto SubROIMatchedEdepShiftedMap =
          wmUtil.MatchEdepsToSubROIs(subROIPropVec, matchedEdepShiftedPtrVec, offset_ADC);
        std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>>
          SubROIMatchedEdepMap;
        for (auto const& key_edepPtrVec_pair : SubROIMatchedEdepShiftedMap)
        {
          auto key = key_edepPtrVec_pair.first;
          for (auto const& shifted_edep_ptr : key_edepPtrVec_pair.second)
          {
            for (size_t i_e = 0; i_e < matchedEdepShiftedPtrVec.size(); ++i_e)
            {
              if (shifted_edep_ptr == matchedEdepShiftedPtrVec[i_e])
              {
                SubROIMatchedEdepMap[key].push_back(matchedEdepPtrVec[i_e]);
                break;
              }
            }
          }
        }

        std::map<sys::WireModUtility::SubROI_Key_t,
                 sys::WireModUtility::ScaleValues_t> SubROIMatchedScalesMap;
        std::map<sys::WireModUtility::SubROI_Key_t,
                 sys::WireModUtility::TruthProperties_t> SubROIMatchedTruthMap;
        std::map<sys::WireModUtility::SubROI_Key_t, int> SubROIScaleReasonMap;
        for (auto const& subroi_prop : subROIPropVec)
        {
          sys::WireModUtility::ScaleValues_t scale_vals;
          auto key = subroi_prop.key;
          auto key_it = SubROIMatchedEdepMap.find(key);
          int scale_reason = 0;

          if (key_it != SubROIMatchedEdepMap.end() && key_it->second.size() > 0)
          {
            auto truth_vals = wmUtil.CalcPropertiesFromEdeps(key_it->second, offset_ADC);
            SubROIMatchedTruthMap[key] = truth_vals;

            if (truth_vals.total_energy < 0.3 && subroi_prop.total_q > 80)
            {
              scale_vals.r_Q     = 1.;
              scale_vals.r_sigma = 1.;
              scale_reason = 1;
            } else
            {
              scale_vals = wmUtil.GetScaleValues(truth_vals, roi_properties);
              if (scale_vals.r_Q != 1. || scale_vals.r_sigma != 1.)
              {
                mf::LogDebug("WireModifierXXW")
                  << "Scaling! Q scale: " << scale_vals.r_Q
                  << "     sigma scale: " << scale_vals.r_sigma;
                isModified = true;
                scale_reason = 2;
              } else {
                scale_reason = 3;
              }
            }
          } else
          {
            scale_vals.r_Q     = 1.;
            scale_vals.r_sigma = 1.;
            scale_reason = 0;
          }

          SubROIMatchedScalesMap[key]  = scale_vals;
          SubROIScaleReasonMap[key]    = scale_reason;
        }

        if (fSetNullScaleIntegral || fSetNullScaleWidth) {
          for (auto& kv : SubROIMatchedScalesMap) {
            if (fSetNullScaleIntegral) kv.second.r_Q     = 1.0;
            if (fSetNullScaleWidth)    kv.second.r_sigma  = 1.0;
          }
        }

        // snapshot of float data before in-place modification (used by ROITree "before" branch)
        std::vector<float> modified_data_before = fSaveROITree ? modified_data : std::vector<float>{};

        wmUtil.ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap, 1.5);

        if (fSaveROITree)
        {
          std::vector<geo::WireID> wireIDs = fWireReadout->ChannelToWire(chanROI.Channel());
          geo::WireID firstWireID = wireIDs.empty() ? geo::WireID{} : wireIDs.front();
          tb_channel = static_cast<Int_t>(chanROI.Channel());
          tb_wire    = static_cast<Int_t>(firstWireID.Wire);
          tb_tpc     = static_cast<Int_t>(firstWireID.TPC);
          tb_plane   = static_cast<Int_t>(my_plane);
          tb_roi_idx = static_cast<Int_t>(i_r);

          tb_roi_begin          = roi_properties.begin;
          tb_roi_before_total_q = roi_properties.total_q;
          tb_roi_before_center  = roi_properties.center;
          tb_roi_before_sigma   = roi_properties.sigma;
          // modified_data holds the float-converted original (before ModifyROI runs)
          tb_roi_before = modified_data_before;

          float mod_total_q = 0.f, mod_center = 0.f, mod_sigma = 0.f;
          for (size_t i_t = 0; i_t < modified_data.size(); ++i_t)
          {
            mod_total_q += modified_data[i_t];
            mod_center  += modified_data[i_t] * (i_t + roi_properties.begin);
          }
          if (mod_total_q > 0.f) mod_center /= mod_total_q;
          for (size_t i_t = 0; i_t < modified_data.size(); ++i_t)
          {
            float t = static_cast<float>(i_t) + roi_properties.begin;
            mod_sigma += modified_data[i_t] * (t - mod_center) * (t - mod_center);
          }
          if (mod_total_q > 0.f) mod_sigma = std::sqrt(mod_sigma / mod_total_q);
          tb_roi_after_total_q = mod_total_q;
          tb_roi_after_center  = mod_center;
          tb_roi_after_sigma   = mod_sigma;
          tb_roi_after = modified_data;

          {
            int n_pt = 0;
            for (auto const& hptr : matchedHitPtrVec)
              if (hptr->GoodnessOfFit() == -1.0f && hptr->DegreesOfFreedom() == 1)
                ++n_pt;
            tb_roi_n_hits             = static_cast<Int_t>(matchedHitPtrVec.size());
            tb_roi_n_pulse_train_hits = static_cast<Int_t>(n_pt);
            tb_roi_is_pulse_train     = (n_pt > 0) ? 1 : 0;
          }

          tb_hit_before_center.clear();
          tb_hit_before_integral.clear();
          tb_hit_before_sigma.clear();
          tb_hit_after_center.clear();
          tb_hit_after_integral.clear();
          tb_hit_after_sigma.clear();
          tb_scale_Q.clear();
          tb_scale_sigma.clear();
          tb_theta_xw.clear();
          tb_truth_x.clear();
          tb_track_id.clear();
          tb_scale_reason.clear();
          {
            std::string gname = (tb_plane < (int)fNameVec_sigma_XXW.size())
                                ? fNameVec_sigma_XXW[tb_plane] : "unknown";
            std::strncpy(tb_xxw_graph_name, gname.c_str(), 127);
            tb_xxw_graph_name[127] = '\0';
          }

          for (auto const& subroi_prop : subROIPropVec)
          {
            auto key        = subroi_prop.key;
            auto scale_vals = SubROIMatchedScalesMap.at(key);

            tb_hit_before_center.push_back(subroi_prop.center);
            tb_hit_before_integral.push_back(subroi_prop.total_q);
            tb_hit_before_sigma.push_back(subroi_prop.sigma);
            tb_hit_after_center.push_back(subroi_prop.center);
            tb_hit_after_integral.push_back(subroi_prop.total_q * static_cast<float>(scale_vals.r_Q));
            tb_hit_after_sigma.push_back(subroi_prop.sigma    * static_cast<float>(scale_vals.r_sigma));
            tb_scale_Q.push_back(static_cast<float>(scale_vals.r_Q));
            tb_scale_sigma.push_back(static_cast<float>(scale_vals.r_sigma));

            int dominant_tid = -999;
            auto edep_it = SubROIMatchedEdepMap.find(key);
            if (edep_it != SubROIMatchedEdepMap.end() && !edep_it->second.empty())
            {
              std::map<int, double> tidEnergy;
              for (auto const& ep : edep_it->second)
                tidEnergy[ep->TrackID()] += ep->E();
              double maxE = 0.;
              for (auto const& te : tidEnergy)
                if (te.second > maxE) { maxE = te.second; dominant_tid = te.first; }
            }
            tb_track_id.push_back(dominant_tid);
            tb_scale_reason.push_back(SubROIScaleReasonMap.at(key));

            float theta_deg = -999.f;
            float truth_x_val = -9999.f;
            auto truth_it = SubROIMatchedTruthMap.find(key);
            if (truth_it != SubROIMatchedTruthMap.end())
            {
              auto const& tv = truth_it->second;
              truth_x_val = static_cast<float>(tv.x);
              geo::TPCGeo const* tpcPtr = fGeometry->PositionToTPCptr(
                  geo::Point_t{(double)tv.x, (double)tv.y, (double)tv.z});
              if (tpcPtr != nullptr)
              {
                try
                {
                  auto const& plane_obj = fWireReadout->Plane(tpcPtr->ID(), roi_properties.view);
                  theta_deg = static_cast<float>(wmUtil.ThetaXW(
                      tv.dxdr, tv.dydr, tv.dzdr, plane_obj.ThetaZ()) * TMath::RadToDeg());
                }
                catch (...) {}
              }
            }
            tb_theta_xw.push_back(theta_deg);
            tb_truth_x.push_back(truth_x_val);
          }

          fROITree->Fill();
        }

        if (fSavePerHitData && it_hit_map != wmUtil.ROIMatchedHitMap.end())
        {
          for (auto const& subroi_prop : subROIPropVec)
          {
            size_t i_subroi = subroi_prop.key.second;
            if (i_subroi >= it_hit_map->second.size()) continue;
            size_t hitVec_idx = it_hit_map->second[i_subroi];

            auto const& scale_vals = SubROIMatchedScalesMap.at(subroi_prop.key);
            (*vec_scale_q    )[hitVec_idx] = static_cast<float>(scale_vals.r_Q);
            (*vec_scale_sigma)[hitVec_idx] = static_cast<float>(scale_vals.r_sigma);

            auto truth_it = SubROIMatchedTruthMap.find(subroi_prop.key);
            if (truth_it != SubROIMatchedTruthMap.end())
            {
              auto const& tv = truth_it->second;
              (*vec_truth_x)[hitVec_idx] = static_cast<float>(tv.x);
              geo::TPCGeo const* tpcPtr = fGeometry->PositionToTPCptr(
                  geo::Point_t{tv.x, tv.y, tv.z});
              if (tpcPtr)
              {
                try {
                  auto const& plane_obj = fWireReadout->Plane(
                      hitVec[hitVec_idx].WireID().planeID());
                  double thetaZ = plane_obj.ThetaZ();
                  (*vec_theta_xw)[hitVec_idx] = static_cast<float>(
                      wmUtil.ThetaXW(tv.dxdr, tv.dydr, tv.dzdr, thetaZ) * TMath::RadToDeg());
                  (*vec_dirx)[hitVec_idx] = static_cast<float>(tv.dxdr);
                  (*vec_diry)[hitVec_idx] = static_cast<float>(tv.dydr);
                  (*vec_dirz)[hitVec_idx] = static_cast<float>(tv.dzdr);
                  double cosG = std::abs(tv.dydr * std::sin(thetaZ) + tv.dzdr * std::cos(thetaZ));
                  float wp    = static_cast<float>(plane_obj.WirePitch());
                  float pval  = (cosG > 1e-6) ? wp / static_cast<float>(cosG) : -999.f;
                  (*vec_pitch)[hitVec_idx] = pval;
                  (*vec_dqdx)[hitVec_idx]  = (pval > 0.f)
                                             ? hitVec[hitVec_idx].Integral() / pval : -999.f;
                }
                catch (...) {}
              }
            }
          }
        }

        std::vector<short> mod_short(modified_data.size());
        std::transform(modified_data.begin(), modified_data.end(), mod_short.begin(), toShort);
        new_rois.add_range(roi_properties.begin, mod_short);
      } // end ROI loop

      new_chanROIs->emplace_back(new_rois, chanROI.Channel(), chanROI.ADCScaleFactor());

      if (fSaveHistsByChannel && isModified)
      {
        readout::ROPID ropID = fWireReadout->ChannelToROP(chanROI.Channel());
        std::string titleStr =  "Cryo-"         + std::to_string(ropID.Cryostat)
                             + "_TPCset-"       + std::to_string(ropID.TPCset)
                             + "_ReadOutPlane-" + std::to_string(ropID.ROP)
                             + "_Channel-"      + std::to_string(chanROI.Channel());
        TH1F* oldChannelHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                        wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        TH1F* newChannelHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                        wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
        {
          float oldSample = (tick < chanROI.Signal().size())              ? chanROI.Signal().at(tick)              : 0;
          float newSample = (tick < new_chanROIs->back().Signal().size()) ? new_chanROIs->back().Signal().at(tick) : 0;
          oldChannelHist->SetBinContent(tick + 1, oldSample);
          newChannelHist->SetBinContent(tick + 1, newSample);
        }
        TH1F* oldChannelHist_toSave = tfs->make<TH1F>(*oldChannelHist);
        TH1F* newChannelHist_toSave = tfs->make<TH1F>(*newChannelHist);
        mf::LogDebug("WireModifierXXW")
          << "Saved histograms " << oldChannelHist_toSave->GetName() << '\n'
          << "             and " << newChannelHist_toSave->GetName();
      }

      if (fSaveHistsByWire && isModified)
      {
        std::vector<geo::WireID> wireIDs = fWireReadout->ChannelToWire(chanROI.Channel());
        mf::LogDebug("WireModifierXXW")
          << "Channel " << chanROI.Channel() << " has " << wireIDs.size() << " wire(s)";
        for (auto const& wireID : wireIDs)
        {
          std::string titleStr =  "Cryo-"  + std::to_string(wireID.Cryostat)
                               + "_TPC-"   + std::to_string(wireID.TPC)
                               + "_Plane-" + std::to_string(wireID.Plane)
                               + "_Wire-"  + std::to_string(wireID.Wire);
          TH1F* oldWireHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                       wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          TH1F* newWireHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units",
                                       wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
          {
            float oldSample = (tick < chanROI.Signal().size())              ? chanROI.Signal().at(tick)              : 0;
            float newSample = (tick < new_chanROIs->back().Signal().size()) ? new_chanROIs->back().Signal().at(tick) : 0;
            oldWireHist->SetBinContent(tick + 1, oldSample);
            newWireHist->SetBinContent(tick + 1, newSample);
          }
          TH1F* oldWireHist_toSave = tfs->make<TH1F>(*oldWireHist);
          TH1F* newWireHist_toSave = tfs->make<TH1F>(*newWireHist);
          mf::LogDebug("WireModifierXXW")
            << "Saved histograms " << oldWireHist_toSave->GetName() << '\n'
            << "             and " << newWireHist_toSave->GetName();
        }
      }

    } // end loop over chanROIs

    evt.put(std::move(new_chanROIs));

    } // end else (ChannelROI mode)

    evt.put(std::move(vec_scale_q),     "scaleQ");
    evt.put(std::move(vec_scale_sigma), "scaleSigma");
    evt.put(std::move(vec_truth_x),     "truthX");
    evt.put(std::move(vec_theta_xw),    "thetaXW");
    evt.put(std::move(vec_is_mc),       "isMC");
    evt.put(std::move(vec_dirx),        "dirX");
    evt.put(std::move(vec_diry),        "dirY");
    evt.put(std::move(vec_dirz),        "dirZ");
    evt.put(std::move(vec_pitch),       "pitch");
    evt.put(std::move(vec_dqdx),        "dQdX");
  }
  DEFINE_ART_MODULE(WireModifierXXW)
} // end namespace
