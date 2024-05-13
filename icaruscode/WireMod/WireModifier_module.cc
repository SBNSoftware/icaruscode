// std inlcudes
#include <string>
#include <vector>

// ROOT includes
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
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "sbnobj/ICARUS/TPC/ChannelROI.h"
#include "icaruscode/TPC/Utilities/ChannelROICreator.h"

// wiremod
#include "sbncode/WireMod/Utility/WireModUtility.hh"
#include "sbncode/WireMod/Utility/WireModUtility.cc"

//namespace
namespace wiremod
{

  class WireModifier : public art::EDProducer
  {
    public:
      explicit WireModifier(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

    private:
      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
      std::string fRatioFileName; // there is where we try to grab the splines/graphs (if they exist)
      std::vector<TSpline3*> fSpline_charge;
      std::vector<TSpline3*> fSpline_sigma;
      std::vector<TGraph2D*> fGraph_charge; 
      std::vector<TGraph2D*> fGraph_sigma;
      art::InputTag fWireLabel; // which wires are we pulling in?
      art::InputTag fHitLabel; // which hits are we pulling in?
      bool fSaveHistsByChannel; // save modified signals by channel?
      bool fSaveHistsByWire; // save modified signals by wire?

  }; // end WireModifier class

  //------------------------------------------------
  WireModifier::WireModifier(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void WireModifier::reconfigure(fhicl::ParameterSet const& pset)
  {
    // For now we aren't configuring anything because we're just getting it running
    fWireLabel  = pset.get<art::InputTag>("WireLabel", "roifinder:PHYSCRATEDATATPCEE");
    fHitLabel   = pset.get<art::InputTag>("HitLabel", "gaushitTPCEE");

    // what, if anything, are we putting in the histogram files
    fSaveHistsByChannel = pset.get<bool>("SaveByChannel", false);
    fSaveHistsByWire    = pset.get<bool>("SaveByWire"   , false);

    // try to read in the graphs/splines from a file
    // if that file does not exist then fake them
    fRatioFileName = pset.get<std::string>("RatioFileName", "NOFILE");
    if (fRatioFileName == "NOFILE")
    {
      // if we can't find the file, take the ratios
      double pnts[3] = {-1000, 0, 1000}; // start and end "points" for the splines/graphs. just make sure this covers the detector
      double vals[3] = {0.9, 0.9, 0.9};
      double vals_inv[3] = {1.0, 1.0, 1.0};
      for (size_t plane = 0; plane < 3; ++plane)
      {
        fSpline_charge.push_back(new TSpline3("dummySpline", pnts[0], pnts[2], vals, 3)    );
        fSpline_sigma .push_back(new TSpline3("dummySpline", pnts[0], pnts[2], vals_inv, 3));
        fGraph_charge .push_back(new TGraph2D(3, pnts, pnts, vals)                         );
        fGraph_sigma  .push_back(new TGraph2D(3, pnts, pnts, vals_inv)                     );
      }
    } else {
      TFile* ratioFile = new TFile(fRatioFileName.c_str(), "READ"); // read only
      // the file exists! pull the ratios
      // for now hard code what the names should be
      for (size_t plane = 0; plane < 3; ++plane)
      {
        fSpline_charge.push_back(ratioFile->Get<TSpline3>(("Spline_charge_" + std::to_string(plane)).c_str()));
        fSpline_sigma .push_back(ratioFile->Get<TSpline3>(("Spline_sigma_"  + std::to_string(plane)).c_str()));
        fGraph_charge .push_back(ratioFile->Get<TGraph2D>(("Graph_charge_"  + std::to_string(plane)).c_str()));
        fGraph_sigma  .push_back(ratioFile->Get<TGraph2D>(("Graph_sigma_"   + std::to_string(plane)).c_str()));
      }
    }

    // we make these things
    produces<std::vector<recob::Wire      >>();
    produces<std::vector<recob::ChannelROI>>();
    //produces<art::Assns<raw::RawDigit, recob::Wire>>();
  }

  //------------------------------------------------
  void WireModifier::produce(art::Event& evt)
  {
    // here's where the "magic" happens
    art::ServiceHandle<art::TFileService> tfs;

    // get a clock and det props for the event
    const detinfo::DetectorPropertiesData detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    // get the things to do the things on
    art::Handle< std::vector<recob::Wire> > wireHandle;
    evt.getByLabel(fWireLabel, wireHandle);
    auto const& wireVec(*wireHandle);

    art::FindManyP<raw::RawDigit> digit_assn(wireHandle, evt, fWireLabel);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepShiftedHandle;
    evt.getByLabel("largeant:TPCActive", edepShiftedHandle);
    auto const& edepShiftedVec(*edepShiftedHandle);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepOrigHandle;
    evt.getByLabel("ionization", edepOrigHandle);
    auto const& edepOrigVec(*edepOrigHandle);
      
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitLabel, hitHandle);
    auto const& hitVec(*hitHandle);

    // put the new stuff somewhere
    std::unique_ptr<std::vector<recob::Wire      >> new_wires(new std::vector<recob::Wire      >());
    std::unique_ptr<std::vector<recob::ChannelROI>> new_crois(new std::vector<recob::ChannelROI>());
    //std::unique_ptr<art::Assns<raw::RawDigit, recob::Wire>> new_digit_assn(new art::Assns<raw::RawDigit, recob::Wire>());

    // let's just try making the damn thing
    sys::WireModUtility wmUtil(fGeometry, detProp); // detector geometry & properties

    // set the ratios
    wmUtil.splines_Charge_X       = fSpline_charge;
    wmUtil.splines_Sigma_X        = fSpline_sigma;
    wmUtil.splines_Charge_XZAngle = fSpline_charge;
    wmUtil.splines_Sigma_XZAngle  = fSpline_sigma;
    wmUtil.splines_Charge_YZAngle = fSpline_charge;
    wmUtil.splines_Sigma_YZAngle  = fSpline_sigma;
    wmUtil.splines_Charge_dEdX    = fSpline_charge;
    wmUtil.splines_Sigma_dEdX     = fSpline_sigma;
    wmUtil.graph2Ds_Charge_YZ     = fGraph_charge;
    wmUtil.graph2Ds_Sigma_YZ      = fGraph_sigma;

    // add some debugging here
    MF_LOG_VERBATIM("WireModifier")
      << "DUMP CONFIG:" << '\n'
      << "---------------------------------------------------" << '\n'
      << "  applyChannelScale:  " << wmUtil.applyChannelScale  << '\n'
      << "  applyXScale:        " << wmUtil.applyXScale        << '\n'
      << "  applyYZScale:       " << wmUtil.applyYZScale       << '\n'
      << "  applyXZAngleScale:  " << wmUtil.applyXZAngleScale  << '\n'
      << "  applyYZAngleScale:  " << wmUtil.applyYZAngleScale  << '\n'
      << "  applydEdXScale:     " << wmUtil.applydEdXScale     << '\n'
      << "  readoutWindowTicks: " << wmUtil.readoutWindowTicks << '\n'
      << "  tickOffset:         " << wmUtil.tickOffset         << '\n'
      << "---------------------------------------------------";

    // do the things
    double offset_ADC = 0; // don't use an offset atm
    MF_LOG_VERBATIM("WireModifier")
      << "Get Edep Map";
    wmUtil.FillROIMatchedEdepMap(edepShiftedVec, wireVec, offset_ADC);
    MF_LOG_VERBATIM("WireModifier")
      << "Got Edep Map." << '\n'
      << "Get Hit Map";
    wmUtil.FillROIMatchedHitMap(hitVec, wireVec);
    MF_LOG_VERBATIM("WireModifier")
      << "Got Hit Map.";

    // loop-de-loop
    for(size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    {
      MF_LOG_DEBUG("WireModifier")
        << "Checking wire " << i_w;

      auto const& wire = wireVec[i_w];

      recob::Wire::RegionsOfInterest_t new_rois;
      recob::ChannelROI::RegionsOfInterest_t new_rois_ints;
      new_rois     .resize(wire.SignalROI().size());
      new_rois_ints.resize(wire.SignalROI().size());

      unsigned int my_plane = geo::kUnknown;
      if        (wire.View() == fGeometry->View(geo::PlaneID(0, 0, 0)))
      {
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on plane 0, view " << wire.View();
        my_plane = 0;
      } else if (wire.View() == fGeometry->View(geo::PlaneID(0, 0, 1))) {
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on plane 1, view " << wire.View();
        my_plane = 1;
      } else if (wire.View() == fGeometry->View(geo::PlaneID(0, 0, 2))) {
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on plane 2, view " << wire.View();
        my_plane = 2;
      }

      if (my_plane == geo::kUnknown)
      {
        MF_LOG_DEBUG("WireModifier")
          << "Wire is on unsupported plane. Skip.";
      }

      // keep track of if this wire is modified
      bool isModified = false;

      for(size_t i_r = 0; i_r < wire.SignalROI().get_ranges().size(); ++i_r)
      {
        MF_LOG_DEBUG("WireModifier")
          << "  Checking ROI " << i_r;

        auto const& range = wire.SignalROI().get_ranges()[i_r];
        sys::WireModUtility::ROI_Key_t roi_key(wire.Channel(), i_r);

        std::vector<float> modified_data(range.data());

        auto it_map = wmUtil.ROIMatchedEdepMap.find(roi_key);
        if(it_map==wmUtil.ROIMatchedEdepMap.end()){
          new_rois     .add_range(range.begin_index(), modified_data);
          new_rois_ints.add_range(range.begin_index(), modified_data);
          MF_LOG_DEBUG("WireModifier")
            << "    Could not find matching Edep. Skip";
          continue;
        }
        std::vector<size_t> matchedEdepIdxVec = it_map->second;
        if(matchedEdepIdxVec.size() == 0)
        {
          new_rois     .add_range(range.begin_index(), modified_data);
          new_rois_ints.add_range(range.begin_index(), modified_data);
          MF_LOG_DEBUG("WireModifier")
            << "    No indices for Edep. Skip";
          continue;
        }
        std::vector<const sim::SimEnergyDeposit*> matchedEdepPtrVec;
        std::vector<const sim::SimEnergyDeposit*> matchedShiftedEdepPtrVec;
        for(auto i_e : matchedEdepIdxVec)
        {
          matchedEdepPtrVec.push_back(&edepOrigVec[i_e]);
          matchedShiftedEdepPtrVec.push_back(&edepShiftedVec[i_e]);
        }
        MF_LOG_DEBUG("WireModifier")
          << "  Found " << matchedShiftedEdepPtrVec.size() << " shifted Edeps";

        std::vector<const recob::Hit*> matchedHitPtrVec;
        auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ) {
          for( auto i_h : it_hit_map->second )
            matchedHitPtrVec.push_back(&hitVec[i_h]);
        }

        MF_LOG_DEBUG("WireModifier")
          << "    Found " << matchedHitPtrVec.size() << " matching hits";

        auto roi_properties = wmUtil.CalcROIProperties(wire, i_r);
        MF_LOG_DEBUG("WireModifier")
          << "    ROI Properties:" << '\n'
          << "                    key:     (" << roi_properties.key.first << ", " << roi_properties.key.second << ")" << '\n'
          << "                    view:    " << roi_properties.view << '\n'
          << "                    begin:   " << roi_properties.begin << '\n'
          << "                    end:     " << roi_properties.end << '\n'
          << "                    total_q: " << roi_properties.total_q << '\n'
          << "                    center:  " << roi_properties.center << '\n'
          << "                    sigma:   " << roi_properties.sigma;

        auto subROIPropVec = wmUtil.CalcSubROIProperties(roi_properties, matchedHitPtrVec);
        
        MF_LOG_DEBUG("WireModifier")
          << "    have " << subROIPropVec.size() << " SubROT";

        auto SubROIMatchedShiftedEdepMap = wmUtil.MatchEdepsToSubROIs(subROIPropVec, matchedShiftedEdepPtrVec, offset_ADC);
        MF_LOG_DEBUG("WireModifier")
          << "    size of SubROIMatchedShiftedEdepMap: " << SubROIMatchedShiftedEdepMap.size();
        std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> SubROIMatchedEdepMap;
        for ( auto const& key_edepPtrVec_pair : SubROIMatchedShiftedEdepMap ) {
          auto key = key_edepPtrVec_pair.first;
          for ( auto const& shifted_edep_ptr : key_edepPtrVec_pair.second ) {
            for ( unsigned int i_e=0; i_e < matchedShiftedEdepPtrVec.size(); i_e++ ) {
              if ( shifted_edep_ptr == matchedShiftedEdepPtrVec[i_e] ) {
                MF_LOG_DEBUG("WireModifier")
                  << "    found matching shifted Edep!";
                SubROIMatchedEdepMap[key].push_back(matchedEdepPtrVec[i_e]);
                break;
              }
            }
          }
        }
        MF_LOG_DEBUG("WireModifier")
          << "    size of SubROIMatchedEdepMap: " << SubROIMatchedEdepMap.size();

        std::map<sys::WireModUtility::SubROI_Key_t, sys::WireModUtility::ScaleValues_t> SubROIMatchedScalesMap;
        for ( auto const& subroi_prop : subROIPropVec ) {
          sys::WireModUtility::ScaleValues_t scale_vals;
          auto key = subroi_prop.key;
          auto key_it =  SubROIMatchedEdepMap.find(key);

          if ( key_it != SubROIMatchedEdepMap.end() && key_it->second.size() > 0 ) {
            auto truth_vals = wmUtil.CalcPropertiesFromEdeps(key_it->second, offset_ADC);

            if ( truth_vals.total_energy < 0.3 && subroi_prop.total_q > 80 ) {
              scale_vals.r_Q     = 1.;
              scale_vals.r_sigma = 1.;
            } 
            else {
              scale_vals = wmUtil.GetScaleValues(truth_vals, roi_properties);
              mf::LogDebug("WireModifier")
                << "Scaling! Q scale: " << scale_vals.r_Q
                << "     sigma sclae: " << scale_vals.r_sigma;
              isModified = true;
            }
          }
          else {
            scale_vals.r_Q     = 1.;
            scale_vals.r_sigma = 1.;
          }

          SubROIMatchedScalesMap[key] = scale_vals;
        }
        
        wmUtil.ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap);
        new_rois     .add_range(roi_properties.begin, modified_data);
        new_rois_ints.add_range(roi_properties.begin, modified_data);
      }

        

      new_wires->emplace_back(new_rois,      wire.Channel(), wire.View());
      new_crois->emplace_back(new_rois_ints, wire.Channel()             );

      if (fSaveHistsByChannel && isModified)
      {
        readout::ROPID ropID = fGeometry->ChannelToROP(wire.Channel());  
        std::string titleStr =  "Cryo-"         + std::to_string(ropID.Cryostat)
                             + "_TPCset-"       + std::to_string(ropID.TPCset)
                             + "_ReadOutPlane-" + std::to_string(ropID.ROP)
                             + "_Channel-"      + std::to_string(wire.Channel());
        TH1F* oldChannelHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        TH1F* newChannelHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
        for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
        {
          float oldSample = (tick <     wire         .Signal().size() ) ?     wire         .Signal().at(tick) : 0;
          float newSample = (tick < new_wires->back().Signal().size() ) ? new_wires->back().Signal().at(tick) : 0;
          oldChannelHist->SetBinContent(tick + 1, oldSample);
          newChannelHist->SetBinContent(tick + 1, newSample);
        }
          
        TH1F* oldChannelHist_toSave = tfs->make<TH1F>(*oldChannelHist);
        TH1F* newChannelHist_toSave = tfs->make<TH1F>(*newChannelHist);
        mf::LogDebug("WireModifier")
          << "Saved histograms " << oldChannelHist_toSave->GetName() << '\n'
          << "             and " << newChannelHist_toSave->GetName();
      }

      if (fSaveHistsByWire && isModified)
      {
        std::vector<geo::WireID> wireIDs = fGeometry->ChannelToWire(wire.Channel());
        mf::LogDebug("WireModifier")
          << "Channel " << wire.Channel() << " has " << wireIDs.size() << " wire(s)";
        for (auto const& wireID : wireIDs)
        {
          std::string titleStr =  "Cryo-"  + std::to_string(wireID.Cryostat)
                               + "_TPC-"   + std::to_string(wireID.TPC)
                               + "_Plane-" + std::to_string(wireID.Plane)
                               + "_Wire-"  + std::to_string(wireID.Wire);
          TH1F* oldWireHist = new TH1F(("Old_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          TH1F* newWireHist = new TH1F(("New_" + titleStr).c_str(), ";Sample;Arbitrary Units", wmUtil.readoutWindowTicks, 0, wmUtil.readoutWindowTicks);
          for (size_t tick = 0; tick < wmUtil.readoutWindowTicks; ++tick)
          {
            float oldSample = (tick <     wire         .Signal().size() ) ?     wire         .Signal().at(tick) : 0;
            float newSample = (tick < new_wires->back().Signal().size() ) ? new_wires->back().Signal().at(tick) : 0;
            oldWireHist->SetBinContent(tick + 1, oldSample);
            newWireHist->SetBinContent(tick + 1, newSample);
          }

          TH1F* oldWireHist_toSave = tfs->make<TH1F>(*oldWireHist);
          TH1F* newWireHist_toSave = tfs->make<TH1F>(*newWireHist);
          mf::LogDebug("WireModifier")
            << "Saved histograms " << oldWireHist_toSave->GetName() << '\n'
            << "             and " << newWireHist_toSave->GetName();
        }
      }

      //auto const& rd_ptrs = digit_assn.at(i_w);
      //for(auto const& rd_ptr : rd_ptrs)
      //  util::CreateAssn(*this, evt, *new_wires, rd_ptr, *new_digit_assn, new_wires->size() - 1);
    } // end loop over wires

    evt.put(std::move(new_wires));
    evt.put(std::move(new_crois));
    //evt.put(std::move(new_digit_assn));
  }
  DEFINE_ART_MODULE(WireModifier)
} // end namespace
