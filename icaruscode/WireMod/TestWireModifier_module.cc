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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"

// wiremod testing
#include "sbncode/WireMod/Utility/WireModUtility.hh"
#include "sbncode/WireMod/Utility/WireModUtility.cc"

//namespace
namespace wiremod
{

  class TestWireModifier : public art::EDProducer
  {
    public:
      explicit TestWireModifier(fhicl::ParameterSet const& pset);
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt) override;

    private:
      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry
      TSpline3 fDummySpline_charge; // we're faking it!
      TSpline3 fDummySpline_sigma;
      TGraph2DErrors fDummyGraph_charge; // yup, still fake 
      TGraph2DErrors fDummyGraph_sigma;
      art::InputTag fWireLabel; // which wires are we pulling in?
      art::InputTag fHitLabel; // which hits are we pulling in?

  }; // end TestWireModifier class

  //------------------------------------------------
  TestWireModifier::TestWireModifier(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    // get the fhicl parameters
    this->reconfigure(pset);
  }

  //------------------------------------------------
  void TestWireModifier::reconfigure(fhicl::ParameterSet const& pset)
  {
    // For now we aren't configuring anything because we're just getting it running
    fWireLabel  = pset.get<art::InputTag>("WireLabel", "roifinder:PHYSCRATEDATATPCEE");
    fHitLabel   = pset.get<art::InputTag>("HitLabel", "gaushitTPCEE");

    // we'll make the splines and the like here since they are FAKE atm
    // later you'd read them in from a file
    double pnts[3] = {-1000, 0, 1000}; // start and end "points" for the splines/graphs. just make sure this covers the detector
    double vals[3] = {1.2, 1.2, 1.2};
    double vals_inv[3] = {1.0 / 1.2, 1.0 / 1.2, 1.0 / 1.2};
    fDummySpline_charge = TSpline3("dummySpline", pnts[0], pnts[2], vals, 3);
    fDummySpline_sigma  = TSpline3("dummySpline", pnts[0], pnts[2], vals_inv, 3);
    fDummyGraph_charge  = TGraph2DErrors(3, pnts, pnts, vals);
    fDummyGraph_sigma   = TGraph2DErrors(3, pnts, pnts, vals_inv);

    // we make these things
    //produces<std::vector<recob::Wire>>();
    //produces<art::Assns<raw::RawDigit, recob::Wire>>();
  }

  //------------------------------------------------
  void TestWireModifier::produce(art::Event& evt)
  {
    // here's where the "magic" happens
    art::ServiceHandle<art::TFileService> tfs;

    // get a clock and det props for the event
    const detinfo::DetectorClocksData    detClock = art::ServiceHandle<detinfo::DetectorClocksService const    >()->DataFor(evt);
    const detinfo::DetectorPropertiesData detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    // get the things to do the things on
    art::Handle< std::vector<recob::Wire> > wireHandle;
    evt.getByLabel(fWireLabel, wireHandle);
    auto const& wireVec(*wireHandle);

    art::FindManyP<raw::RawDigit> digit_assn(wireHandle, evt, fWireLabel);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepShiftedHandle;
    evt.getByLabel("ionization", edepShiftedHandle);
    auto const& edepShiftedVec(*edepShiftedHandle);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepOrigHandle;
    evt.getByLabel("largeant:TPCActive", edepOrigHandle);
    auto const& edepOrigVec(*edepOrigHandle);
      
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitLabel, hitHandle);
    auto const& hitVec(*hitHandle);

    // put the new stuff somewhere
    std::unique_ptr<std::vector<recob::Wire>> new_wires(new std::vector<recob::Wire>());
    std::unique_ptr<art::Assns<raw::RawDigit, recob::Wire>> new_digit_assn(new art::Assns<raw::RawDigit, recob::Wire>());

    // let's just try making the damn thing
    sys::WireModUtility wmUtil(fGeometry, detClock, detProp); // geometry, clock, properties

    // it's all fake, so each plane is the same
    for (size_t i = 0; i < 3; ++i)
    {
      wmUtil.splines_Charge_X      .push_back(&fDummySpline_charge);
      wmUtil.splines_Sigma_X       .push_back(&fDummySpline_sigma);
      wmUtil.splines_Charge_XZAngle.push_back(&fDummySpline_charge);
      wmUtil.splines_Sigma_XZAngle .push_back(&fDummySpline_sigma);
      wmUtil.splines_Charge_YZAngle.push_back(&fDummySpline_charge);
      wmUtil.splines_Sigma_YZAngle .push_back(&fDummySpline_sigma);
      wmUtil.splines_Charge_dEdX   .push_back(&fDummySpline_charge);
      wmUtil.splines_Sigma_dEdX    .push_back(&fDummySpline_sigma);

      wmUtil.graph2Ds_Charge_YZ.push_back(&fDummyGraph_charge);
      wmUtil.graph2Ds_Sigma_YZ .push_back(&fDummyGraph_sigma);
    } // end loop over planes

    // add some debugging here
    MF_LOG_VERBATIM("TestWireModifier")
      << "DUMP CONFIG:" << '\n'
      << "---------------------------------------------------" << '\n'
      << "  applyXScale:        " << wmUtil.applyXScale        << '\n'
      << "  applyYZScale:       " << wmUtil.applyYZScale       << '\n'
      << "  applyXZAngleScale:  " << wmUtil.applyXZAngleScale  << '\n'
      << "  applyYZAngleScale:  " << wmUtil.applyYZAngleScale  << '\n'
      << "  applydEdXScale:     " << wmUtil.applydEdXScale     << '\n'
      << "  readoutWindowTicks: " << wmUtil.readoutWindowTicks << '\n'
      << "  tickOffset:         " << wmUtil.tickOffset         << '\n'
      //<< "  wiresPercm:         " << wmUtil.wiresPercm         << '\n'
      //<< "  originWire_plane0:  " << wmUtil.originWire_plane0  << '\n'
      //<< "  originWire_plane1:  " << wmUtil.originWire_plane1  << '\n'
      //<< "  originWire_plane2:  " << wmUtil.originWire_plane2  << '\n'
      //<< "  ticksPercm:         " << wmUtil.ticksPercm         << '\n'
      //<< "  zeroTick:           " << wmUtil.zeroTick           << '\n'
      //<< "  angle_plane0:       " << wmUtil.angle_plane0       << '\n'
      //<< "  angle_plane1:       " << wmUtil.angle_plane1       << '\n'
      //<< "  angle_plane2:       " << wmUtil.angle_plane2       << '\n'
      //<< "  plane0_boundary:    " << wmUtil.plane0_boundary    << '\n'
      //<< "  plane1_boundary:    " << wmUtil.plane1_boundary    << '\n'
      //<< "  plane2_boundary:    " << wmUtil.plane2_boundary    << '\n'
      << "---------------------------------------------------";

    // add some debugging here

    // timing things
    double nu_t = 0.;
    //art::Handle<std::vector<simb::MCTruth> > mctruth_h;
    //evt.getByLabel("generator", mctruth_h);
    //if(not mctruth_h.isValid())
    //{
    //  MF_LOG_VERBATIM("TestWireModifier")
    //    << "Invalid truth info. Bail.";
    //    return;
    //}
    //std::vector<art::Ptr<simb::MCTruth> > MCTruthVec;
    //art::fill_ptr_vector(MCTruthVec, mctruth_h);
    //for (auto& mctruth : MCTruthVec ) {
    //  const simb::MCParticle& truth_mcparticle = mctruth->GetParticle(0);
    //  nu_t = truth_mcparticle.T( 0 );
    //}
    double offset_ADC = (nu_t + detClock.TriggerOffsetTPC()) / detClock.TPCClock().TickPeriod();
    //double offset_ADC = 0;
    MF_LOG_VERBATIM("TestWireModifier")
      << "nu time is " << nu_t << ", so offset_ADC is " << offset_ADC << " = (" <<nu_t << " + " << detClock.TriggerOffsetTPC() << ") / " << detClock.TPCClock().TickPeriod();

    // do the things
    MF_LOG_VERBATIM("TestWireModifier")
      << "Get Edep Map";
    wmUtil.FillROIMatchedEdepMap(edepShiftedVec, wireVec, offset_ADC);
    MF_LOG_VERBATIM("TestWireModifier")
      << "Got Edep Map." << '\n'
      << "Get Hit Map";
    wmUtil.FillROIMatchedHitMap(hitVec, wireVec);
    MF_LOG_VERBATIM("TestWireModifier")
      << "Got Hit Map.";

    // match old and new wires
    std::vector<std::pair<size_t, size_t>> oldNew_wireIdx;

    // loop-de-loop
    for(size_t i_w = 0; i_w < wireVec.size(); ++i_w)
    {
      bool isModified = false;

      MF_LOG_DEBUG("TestWireModifier")
        << "Checking wire " << i_w;

      auto const& wire = wireVec[i_w];

      recob::Wire::RegionsOfInterest_t new_rois;
      new_rois.resize(wire.SignalROI().size());

      unsigned int my_plane = geo::kUnknown;
      if        (wire.View() == fGeometry->View(geo::PlaneID(0, 0, 0)))
      {
        MF_LOG_DEBUG("TestWireModifier")
          << "Wire is on plane 0, view " << wire.View();
        my_plane = 0;
      } else if (wire.View() == fGeometry->View(geo::PlaneID(0, 0, 1))) {
        MF_LOG_DEBUG("TestWireModifier")
          << "Wire is on plane 1, view " << wire.View();
        my_plane = 1;
      } else if (wire.View() == fGeometry->View(geo::PlaneID(0, 0, 2))) {
        MF_LOG_DEBUG("TestWireModifier")
          << "Wire is on plane 2, view " << wire.View();
        my_plane = 2;
      }

      if (my_plane == geo::kUnknown)
      {
        MF_LOG_DEBUG("TestWireModifier")
          << "Wire is on unsupported plane. Skip.";
      }

      for(size_t i_r = 0; i_r < wire.SignalROI().get_ranges().size(); ++i_r)
      {
        MF_LOG_DEBUG("TestWireModifier")
          << "  Checking ROI " << i_r;

        auto const& range = wire.SignalROI().get_ranges()[i_r];
        sys::WireModUtility::ROI_Key_t roi_key(wire.Channel(),i_r);

        std::vector<float> modified_data(range.data());

        auto it_map = wmUtil.ROIMatchedEdepMap.find(roi_key);
        if(it_map==wmUtil.ROIMatchedEdepMap.end()){
          new_rois.add_range(range.begin_index(),modified_data);
          MF_LOG_DEBUG("TestWireModifier")
            << "    Could not find matching Edep. Skip";
          continue;
        }
        std::vector<size_t> matchedEdepIdxVec = it_map->second;
        if(matchedEdepIdxVec.size() == 0)
        {
          new_rois.add_range(range.begin_index(),modified_data);
          MF_LOG_DEBUG("TestWireModifier")
            << "    No idices for Edep. Skip";
          continue;
        }
        std::vector<const sim::SimEnergyDeposit*> matchedEdepPtrVec;
        std::vector<const sim::SimEnergyDeposit*> matchedShiftedEdepPtrVec;
        for(auto i_e : matchedEdepIdxVec)
        {
          matchedEdepPtrVec.push_back(&edepOrigVec[i_e]);
          matchedShiftedEdepPtrVec.push_back(&edepShiftedVec[i_e]);
        }
        MF_LOG_DEBUG("TestWireModifier")
          << "  Found " << matchedShiftedEdepPtrVec.size() << " shifted Edeps";

        std::vector<const recob::Hit*> matchedHitPtrVec;
        auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ) {
          for( auto i_h : it_hit_map->second )
            matchedHitPtrVec.push_back(&hitVec[i_h]);
        }

        MF_LOG_DEBUG("TestWireModifier")
          << "    Found " << matchedHitPtrVec.size() << " matching hits";

        auto roi_properties = wmUtil.CalcROIProperties(range);
        roi_properties.key   = roi_key;
        roi_properties.plane = my_plane;
        //roi_properties.plane = wire.View();
        MF_LOG_DEBUG("TestWireModifier")
          << "    ROI Properties:" << '\n'
          << "                    key:     (" << roi_properties.key.first << ", " << roi_properties.key.second << ")" << '\n'
          << "                    plane:   " << roi_properties.plane << '\n'
          << "                    begin:   " << roi_properties.begin << '\n'
          << "                    end:     " << roi_properties.end << '\n'
          << "                    total_q: " << roi_properties.total_q << '\n'
          << "                    center:  " << roi_properties.center << '\n'
          << "                    sigma:   " << roi_properties.sigma;

        auto subROIPropVec = wmUtil.CalcSubROIProperties(roi_properties, matchedHitPtrVec);
        
        MF_LOG_DEBUG("TestWireModifier")
          << "    have " << subROIPropVec.size() << " SubROT";

        auto SubROIMatchedShiftedEdepMap = wmUtil.MatchEdepsToSubROIs(subROIPropVec, matchedShiftedEdepPtrVec, offset_ADC);
        MF_LOG_DEBUG("TestWireModifier")
          << "    size of SubROIMatchedShiftedEdepMap: " << SubROIMatchedShiftedEdepMap.size();
        std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> SubROIMatchedEdepMap;
        for ( auto const& key_edepPtrVec_pair : SubROIMatchedShiftedEdepMap ) {
          auto key = key_edepPtrVec_pair.first;
          for ( auto const& shifted_edep_ptr : key_edepPtrVec_pair.second ) {
            for ( unsigned int i_e=0; i_e < matchedShiftedEdepPtrVec.size(); i_e++ ) {
              if ( shifted_edep_ptr == matchedShiftedEdepPtrVec[i_e] ) {
                MF_LOG_DEBUG("TestWireModifier")
                  << "    found matching shifted Edep!";
                SubROIMatchedEdepMap[key].push_back(matchedEdepPtrVec[i_e]);
                break;
              }
            }
          }
        }
        MF_LOG_DEBUG("TestWireModifier")
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
              isModified = false;
            } 
            else {
              scale_vals = wmUtil.GetScaleValues(truth_vals, roi_properties);
              mf::LogVerbatim("TestWireModifier")
                << "Scaling! Q scale: " << scale_vals.r_Q
                << "     sigma sclae: " << scale_vals.r_sigma;
              isModified = true;
            }
          }
          else {
            scale_vals.r_Q     = 1.;
            scale_vals.r_sigma = 1.;
            isModified = false;
          }

          SubROIMatchedScalesMap[key] = scale_vals;
        }

        // store old ROI
        TH1F* oldROI = new TH1F(("oldROI__Wire_"+std::to_string(i_w)+"_ROI_"+std::to_string(i_r)).c_str(), ";Sample;Arbitrary Units", modified_data.size(), 0, modified_data.size());
        for (size_t bin = 1; bin < modified_data.size() + 1; ++bin)
          oldROI->SetBinContent(bin, modified_data[bin - 1]);

        wmUtil.ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap);

        // store new ROI
        TH1F* newROI = new TH1F(("newROI_Wire_"+std::to_string(i_w)+"_ROI_"+std::to_string(i_r)).c_str(), ";Sample;Arbitrary Units", modified_data.size(), 0, modified_data.size());
        for (size_t bin = 1; bin < modified_data.size() + 1; ++bin)
          newROI->SetBinContent(bin, modified_data[bin - 1]);

        if (isModified)
        {
          TH1F* oldROI_saved = tfs->make<TH1F>(*oldROI);
          TH1F* newROI_saved = tfs->make<TH1F>(*newROI);
          mf::LogVerbatim("TestWireModifier")
            << "Wrote histograms " << oldROI_saved->GetName() << " and " << newROI_saved->GetName();
        }

        new_rois.add_range(roi_properties.begin, modified_data);
      }

      new_wires->emplace_back(new_rois, wire.Channel(), wire.View());

      //auto const& rd_ptrs = digit_assn.at(i_w);
      //for(auto const& rd_ptr : rd_ptrs)
      //  util::CreateAssn(*this, evt, *new_wires, rd_ptr, *new_digit_assn, new_wires->size() - 1);
    } // end loop over wires

    //evt.put(std::move(new_wires));
    //evt.put(std::move(new_digit_assn));
  }
  DEFINE_ART_MODULE(TestWireModifier)
} // end namespace
