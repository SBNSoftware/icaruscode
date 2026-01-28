// std inlcudes
#include <string>
#include <vector>

// ROOT includes
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"

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
      void shapeGraph(TGraph2D& graph);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

    private:
      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
      const geo::WireReadoutGeom* fWireReadout = &(art::ServiceHandle<geo::WireReadout const>()->Get());
      std::string fRatioFileName; // there is where we try to grab the splines/graphs (if they exist)
      std::vector<TGraph2D*> fGraph_charge_XXW;
      std::vector<TGraph2D*> fGraph_sigma_XXW;
      art::InputTag fWireLabel; // which wires are we pulling in?
      art::InputTag fHitLabel;  // which hits are we pulling in?
      art::InputTag fEDepLabel; // which are the EDeps?
      art::InputTag fEDepShiftedLabel; // which are the EDeps?
      bool fSaveHistsByChannel; // save modified signals by channel?
      bool fSaveHistsByWire;    // save modified signals by wire?
      bool fInRads;             // is the TGraph2D angle axis in radians?

  }; // end WireModifierXXW class

  //------------------------------------------------
  WireModifierXXW::WireModifierXXW(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void WireModifierXXW::shapeGraph(TGraph2D& graph)
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
  void WireModifierXXW::reconfigure(fhicl::ParameterSet const& pset)
  {
    fWireLabel = pset.get<art::InputTag>("WireLabel");
    fHitLabel  = pset.get<art::InputTag>("HitLabel");
    fEDepLabel = pset.get<art::InputTag>("EDepLabel");
    fEDepShiftedLabel = pset.get<art::InputTag>("EDepLabel");

    // what, if anything, are we putting in the histogram files
    fSaveHistsByChannel = pset.get<bool>("SaveByChannel", false);
    fSaveHistsByWire    = pset.get<bool>("SaveByWire"   , false);

    // do we need to reshape the TGraph2D?
    fInRads = pset.get<bool>("InRadians", true);

    // try to read in the graphs/splines from a file
    // if that file does not exist then fake them
    fRatioFileName = pset.get<std::string>("RatioFileName", "NOFILE");
    if (fRatioFileName == "NOFILE")
    {
      mf::LogVerbatim("WireModifierXXW")
        << "WireModifierXXW::reconfigure - No ratio file given. No scaling is applied...";
    } else {
      char* icaruscode_dir = std::getenv("ICARUSCODE_DIR");
      assert(icaruscode_dir
        && "WireModifierXXW::reconfigure - ICARUSCODE_DIR environment variable must be set!");
      std::string dir_path = std::string(icaruscode_dir) + "/root/WireMod";
      mf::LogDebug("WireModifierXXW")
        << "WireModifierXXW::reconfigure - Get file " << fRatioFileName
        << " from directory " << dir_path;
      TFile* ratioFile = new TFile((dir_path + "/" + fRatioFileName).c_str(), "READ"); // read only
      // the file exists! pull the ratios
      assert(ratioFile && "WireModifierXXW::reconfigure - WireMod Ratio File Must Exist!");
      assert(!ratioFile->IsZombie()
        && "WireModifierXXW::reconfigure - WireMod Ratio File Must Not be a Zombie!");
      mf::LogVerbatim("WireModifierXXW")
        << "WireModifierXXW::reconfigure - Getting XXW Scales...";
      std::vector<std::string> nameVec_charge_XXW = pset.get<std::vector<std::string>>("XXWScaleHeight");
      for (auto const& name : nameVec_charge_XXW)
      {
        mf::LogDebug("WireModifierXXW")
          << "WireModifierXXW::reconfigure - Looking for " << name << " in TFile " << fRatioFileName << "...";
        TGraph2D* temp = static_cast<TGraph2D*>(ratioFile->Get<TGraph2D>(name.c_str())->Clone());
        if (temp != nullptr)
        {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...found";
          if (not fInRads)
            shapeGraph(*temp);
          fGraph_charge_XXW.push_back(temp);
        } else {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...not found";
        }
      }
      std::vector<std::string> nameVec_sigma_XXW = pset.get<std::vector<std::string>>("XXWScaleWidth");
      for (auto const& name : nameVec_sigma_XXW)
      {
        mf::LogDebug("WireModifierXXW")
          << "WireModifierXXW::reconfigure - Looking for " << name << " in TFile " << fRatioFileName << "...";
        TGraph2D* temp = static_cast<TGraph2D*>(ratioFile->Get<TGraph2D>(name.c_str())->Clone());
        if (temp != nullptr)
        {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...found";
          if (not fInRads)
            shapeGraph(*temp);
          fGraph_sigma_XXW.push_back(temp);
        } else {
          mf::LogDebug("WireModifierXXW")
            << "WireModifierXXW::reconfigure -  ...not found";
        }
      }
    }

    // we make these things
    produces<std::vector<recob::Wire>>();
  }

  //------------------------------------------------
  void WireModifierXXW::produce(art::Event& evt)
  {
    // here's where the "magic" happens
    art::ServiceHandle<art::TFileService> tfs;

    // get a clock and det props for the event
    const detinfo::DetectorPropertiesData detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    // get the things to do the things on
    art::Handle< std::vector<recob::Wire> > wireHandle;
    evt.getByLabel(fWireLabel, wireHandle);
    auto const& wireVec(*wireHandle);

    art::FindManyP<raw::RawDigit> digit_assn(wireHandle, evt, fWireLabel);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepHandle;
    evt.getByLabel(fEDepLabel, edepHandle);
    auto const& edepVec(*edepHandle);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepShiftedHandle;
    evt.getByLabel(fEDepShiftedLabel, edepShiftedHandle);
    auto const& edepShiftedVec(*edepShiftedHandle);

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitLabel, hitHandle);
    auto const& hitVec(*hitHandle);

    // put the new stuff somewhere
    std::unique_ptr<std::vector<recob::Wire>> new_wires(new std::vector<recob::Wire>());

    sys::WireModUtility wmUtil(fGeometry, fWireReadout, detProp); // detector geometry & properties
    wmUtil.applyChannelScale = false;
    wmUtil.applyXScale       = false;
    wmUtil.applyYZScale      = false;
    wmUtil.applyXZAngleScale = false;
    wmUtil.applyYZAngleScale = false;
    wmUtil.applydEdXScale    = false;
    wmUtil.graph2Ds_Charge_XXW = fGraph_charge_XXW;
    wmUtil.graph2Ds_Sigma_XXW  = fGraph_sigma_XXW;

    // add some debugging here
    MF_LOG_VERBATIM("WireModifierXXW")
      << "DUMP CONFIG:" << '\n'
      << "---------------------------------------------------" << '\n'
      << "  applyChannelScale:  " << wmUtil.applyChannelScale  << '\n'
      << "  applyXScale:        " << wmUtil.applyXScale        << '\n'
      << "  applyYZScale:       " << wmUtil.applyYZScale       << '\n'
      << "  applyXXWScale:      " << wmUtil.applyXXWScale      << '\n'
      << "  applyXZAngleScale:  " << wmUtil.applyXZAngleScale  << '\n'
      << "  applyYZAngleScale:  " << wmUtil.applyYZAngleScale  << '\n'
      << "  applydEdXScale:     " << wmUtil.applydEdXScale     << '\n'
      << "  readoutWindowTicks: " << wmUtil.readoutWindowTicks << '\n'
      << "  tickOffset:         " << wmUtil.tickOffset         << '\n'
      << "---------------------------------------------------";

    // do the things
    double offset_ADC = 0; // don't use an offset atm
    MF_LOG_VERBATIM("WireModifierXXW")
      << "Get Edep Map";
    wmUtil.FillROIMatchedEdepMap(edepShiftedVec, wireVec, offset_ADC);
    MF_LOG_VERBATIM("WireModifierXXW")
      << "Got Edep Map." << '\n'
      << "Get Hit Map";
    wmUtil.FillROIMatchedHitMap(hitVec, wireVec);
    MF_LOG_VERBATIM("WireModifierXXW")
      << "Got Hit Map.";

    // loop-de-loop
    for(size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    {
      MF_LOG_DEBUG("WireModifierXXW")
        << "Checking wire " << i_w;

      auto const& wire = wireVec.at(i_w);
      if (wire.NSignal() == 0)
        continue;

      recob::Wire::RegionsOfInterest_t new_rois;
      new_rois.resize(wire.SignalROI().size());

      unsigned int my_plane = geo::kUnknown;
      if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 0)).View())
      {
        MF_LOG_DEBUG("WireModifierXXW")
          << "Wire is on plane 0, view " << wire.View();
        my_plane = 0;
      } else if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 1)).View()) {
        MF_LOG_DEBUG("WireModifierXXW")
          << "Wire is on plane 1, view " << wire.View();
        my_plane = 1;
      } else if (wire.View() == fWireReadout->Plane(geo::PlaneID(0, 0, 2)).View()) {
        MF_LOG_DEBUG("WireModifierXXW")
          << "Wire is on plane 2, view " << wire.View();
        my_plane = 2;
      }

      if (my_plane == geo::kUnknown)
      {
        MF_LOG_DEBUG("WireModifierXXW")
          << "Wire is on unsupported plane. Skip.";
      }

      // keep track of if this wire is modified
      bool isModified = false;

      for(size_t i_r = 0; i_r < wire.SignalROI().get_ranges().size(); ++i_r)
      {
        MF_LOG_DEBUG("WireModifierXXW")
          << "  Checking ROI " << i_r;

        auto const& range = wire.SignalROI().get_ranges()[i_r];
        sys::WireModUtility::ROI_Key_t roi_key(wire.Channel(), i_r);

        std::vector<float> modified_data(range.data());

        auto it_map = wmUtil.ROIMatchedEdepMap.find(roi_key);
        if(it_map==wmUtil.ROIMatchedEdepMap.end()){
          new_rois.add_range(range.begin_index(), modified_data);
          MF_LOG_DEBUG("WireModifierXXW")
            << "    Could not find matching Edep. Skip";
          continue;
        }
        std::vector<size_t> matchedEdepIdxVec = it_map->second;
        if(matchedEdepIdxVec.size() == 0)
        {
          new_rois.add_range(range.begin_index(), modified_data);
          MF_LOG_DEBUG("WireModifierXXW")
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
        MF_LOG_DEBUG("WireModifierXXW")
          << "  Found " << matchedEdepPtrVec.size() << " Edeps";

        std::vector<const recob::Hit*> matchedHitPtrVec;
        auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ) {
          for( auto i_h : it_hit_map->second )
            matchedHitPtrVec.push_back(&hitVec[i_h]);
        }

        MF_LOG_DEBUG("WireModifierXXW")
          << "    Found " << matchedHitPtrVec.size() << " matching hits";

        auto roi_properties = wmUtil.CalcROIProperties(wire, i_r);
        MF_LOG_DEBUG("WireModifierXXW")
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
        
        MF_LOG_DEBUG("WireModifierXXW")
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
        for (auto const& subroi_prop : subROIPropVec)
        {
          sys::WireModUtility::ScaleValues_t scale_vals;
          auto key = subroi_prop.key;
          auto key_it =  SubROIMatchedEdepMap.find(key);

          if (key_it != SubROIMatchedEdepMap.end() && key_it->second.size() > 0)
          {
            auto truth_vals = wmUtil.CalcPropertiesFromEdeps(key_it->second, offset_ADC);

            if (truth_vals.total_energy < 0.3 && subroi_prop.total_q > 80)
            {
              scale_vals.r_Q     = 1.;
              scale_vals.r_sigma = 1.;
            } else
            {
              scale_vals = wmUtil.GetScaleValues(truth_vals, roi_properties);
              mf::LogDebug("WireModifierXXW")
                << "Scaling! Q scale: " << scale_vals.r_Q
                << "     sigma scale: " << scale_vals.r_sigma;
              isModified = true;
            }
          } else
          {
            scale_vals.r_Q     = 1.;
            scale_vals.r_sigma = 1.;
          }

          SubROIMatchedScalesMap[key] = scale_vals;
        }
        
        wmUtil.ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap);
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
  }
  DEFINE_ART_MODULE(WireModifierXXW)
} // end namespace
