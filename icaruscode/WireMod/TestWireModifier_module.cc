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
      TSpline3 fDummySpline; // we're faking it!
      TGraph2DErrors fDummyGraph; // yup, still fake 

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
    MF_LOG_VERBATIM("TestWireModifier")
      << "Here's where you'd read in the fhicls and such";

    // we'll make the splines and the like here since they are FAKE atm
    // later you'd read them in from a file
    double pnts[3] = {-1000, 0, 1000}; // start and end "points" for the splines/graphs. just make sure this covers the detector
    double vals[3] = {0.288675135, 0.288675135, 0.288675135}; // just double everything, it doesn't matter right now (1/sqrt(12) and 12 scale factors are applied)
    fDummySpline = TSpline3("dummySpline", pnts[0], pnts[2], vals, 3);
    fDummyGraph = TGraph2DErrors(3, pnts, pnts, vals);
  }

  //------------------------------------------------
  void TestWireModifier::produce(art::Event& evt)
  {
    // here's where the "magic" happens

    // get a clock and det props for the event
    const detinfo::DetectorClocksData    detClock = art::ServiceHandle<detinfo::DetectorClocksService const    >()->DataFor(evt);
    const detinfo::DetectorPropertiesData detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    // let's just try making the damn thing
    // take cryo 0 tpc 0 becasue who cares?
    sys::WireModUtility wmUtil(fGeometry, detClock, detProp, geo::CryostatID(0), geo::TPCID(0, 0)); // geometry, clock, properties, cry, tpc

    // it's all fake, so each plane is the same
    size_t nPlanes = fGeometry->Cryostat(geo::CryostatID(0)).TPC(geo::TPCID(0, 0)).Nplanes();
    MF_LOG_VERBATIM("TestWireModifier")
      << "TPC 0 Cryo 0 has " << nPlanes << " planes";
    for (size_t i = 0; i < nPlanes; ++i)
    {
      wmUtil.splines_Charge_X      .push_back(&fDummySpline);
      wmUtil.splines_Sigma_X       .push_back(&fDummySpline);
      wmUtil.splines_Charge_XZAngle.push_back(&fDummySpline);
      wmUtil.splines_Sigma_XZAngle .push_back(&fDummySpline);
      wmUtil.splines_Charge_YZAngle.push_back(&fDummySpline);
      wmUtil.splines_Sigma_YZAngle .push_back(&fDummySpline);
      wmUtil.splines_Charge_dEdX   .push_back(&fDummySpline);
      wmUtil.splines_Sigma_dEdX    .push_back(&fDummySpline);

      wmUtil.graph2Ds_Charge_YZ.push_back(&fDummyGraph);
      wmUtil.graph2Ds_Sigma_YZ .push_back(&fDummyGraph);
    }

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
      << "  wiresPercm:         " << wmUtil.wiresPercm         << '\n'
      << "  originWire_plane0:  " << wmUtil.originWire_plane0  << '\n'
      << "  originWire_plane1:  " << wmUtil.originWire_plane1  << '\n'
      << "  originWire_plane2:  " << wmUtil.originWire_plane2  << '\n'
      << "  ticksPercm:         " << wmUtil.ticksPercm         << '\n'
      << "  zeroTick:           " << wmUtil.zeroTick           << '\n'
      << "  angle_plane0:       " << wmUtil.angle_plane0       << '\n'
      << "  angle_plane1:       " << wmUtil.angle_plane1       << '\n'
      << "  angle_plane2:       " << wmUtil.angle_plane2       << '\n'
      << "  plane0_boundary:    " << wmUtil.plane0_boundary    << '\n'
      << "  plane1_boundary:    " << wmUtil.plane1_boundary    << '\n'
      << "  plane2_boundary:    " << wmUtil.plane2_boundary    << '\n'
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
    MF_LOG_VERBATIM("TestWireModifier")
      << "nu time is " << nu_t << ", so offset_ADC is " << offset_ADC;

    // get the things to do the things on
    art::Handle< std::vector<recob::Wire> > wireHandle;
    evt.getByLabel("decon1droi:PHYSCRATEDATATPCEE", wireHandle);
    auto const& wireVec(*wireHandle);

    art::FindManyP<raw::RawDigit> digit_assn(wireHandle, evt, "decon1droi:PHYSCRATEDATATPCEE");

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepShiftedHandle;
    evt.getByLabel("ionization", edepShiftedHandle);
    auto const& edepShiftedVec(*edepShiftedHandle);

    art::Handle< std::vector<sim::SimEnergyDeposit> > edepOrigHandle;
    evt.getByLabel("largeant:TPCActive", edepOrigHandle);
    auto const& edepOrigVec(*edepOrigHandle);
    
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel("gaushitTPCEE", hitHandle);
    auto const& hitVec(*hitHandle);

    // put the new stuff somewhere
    std::unique_ptr<std::vector<recob::Wire>> new_wires(new std::vector<recob::Wire>());
    std::unique_ptr<art::Assns<raw::RawDigit, recob::Wire>> new_digit_assn(new art::Assns<raw::RawDigit, recob::Wire>());

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

      MF_LOG_VERBATIM("TestWireModifier")
        << "Checking wire " << i_w;

      auto const& wire = wireVec[i_w];

      recob::Wire::RegionsOfInterest_t new_rois;
      new_rois.resize(wire.SignalROI().size());

      unsigned int my_plane = wire.View();
      if (my_plane >= nPlanes)
        continue;

      for(size_t i_r = 0; i_r < wire.SignalROI().get_ranges().size(); ++i_r)
      {

        auto const& range = wire.SignalROI().get_ranges()[i_r];
        sys::WireModUtility::ROI_Key_t roi_key(wire.Channel(),i_r);

        std::vector<float> modified_data(range.data());

        auto it_map = wmUtil.ROIMatchedEdepMap.find(roi_key);
        if(it_map==wmUtil.ROIMatchedEdepMap.end()){
          new_rois.add_range(range.begin_index(),modified_data);
          continue;
        }
        std::vector<size_t> matchedEdepIdxVec = it_map->second;
        if(matchedEdepIdxVec.size() == 0)
        {
          new_rois.add_range(range.begin_index(),modified_data);
          continue;
        }
        std::vector<const sim::SimEnergyDeposit*> matchedEdepPtrVec;
        std::vector<const sim::SimEnergyDeposit*> matchedShiftedEdepPtrVec;
        for(auto i_e : matchedEdepIdxVec)
        {
          matchedEdepPtrVec.push_back(&edepOrigVec[i_e]);
          matchedShiftedEdepPtrVec.push_back(&edepShiftedVec[i_e]);
        }

        std::vector<const recob::Hit*> matchedHitPtrVec;
        auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
        if( it_hit_map != wmUtil.ROIMatchedHitMap.end() ) {
          for( auto i_h : it_hit_map->second )
            matchedHitPtrVec.push_back(&hitVec[i_h]);
        }

        auto roi_properties = wmUtil.CalcROIProperties(range);
        roi_properties.key   = roi_key;
        roi_properties.plane = my_plane;

        auto subROIPropVec = wmUtil.CalcSubROIProperties(roi_properties, matchedHitPtrVec);

        auto SubROIMatchedShiftedEdepMap = wmUtil.MatchEdepsToSubROIs(subROIPropVec, matchedShiftedEdepPtrVec, offset_ADC);
        std::map<sys::WireModUtility::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> SubROIMatchedEdepMap;
        for ( auto const& key_edepPtrVec_pair : SubROIMatchedShiftedEdepMap ) {
          auto key = key_edepPtrVec_pair.first;
          for ( auto const& shifted_edep_ptr : key_edepPtrVec_pair.second ) {
            for ( unsigned int i_e=0; i_e < matchedShiftedEdepPtrVec.size(); i_e++ ) {
              if ( shifted_edep_ptr == matchedShiftedEdepPtrVec[i_e] ) {
                SubROIMatchedEdepMap[key].push_back(matchedEdepPtrVec[i_e]);
                break;
              }
            }
          }
        }

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
        TH1F* oldROI = new TH1F(("oldROI_Wire_"+std::to_string(i_w)+"_ROI_"+std::to_string(i_r)).c_str(), ";Sample;Arbitrary Units", modified_data.size(), 0, modified_data.size());
        for (size_t bin = 1; bin < modified_data.size() + 1; ++bin)
          oldROI->SetBinContent(bin, modified_data[bin - 1]);

        wmUtil.ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap);

        // store new ROI
        TH1F* newROI = new TH1F(("newROI_Wire_"+std::to_string(i_w)+"_ROI_"+std::to_string(i_r)).c_str(), ";Sample;Arbitrary Units", modified_data.size(), 0, modified_data.size());
        for (size_t bin = 1; bin < modified_data.size() + 1; ++bin)
          newROI->SetBinContent(bin, modified_data[bin - 1]);

        if (isModified)
        {
          art::ServiceHandle<art::TFileService> tfs;
          TH1F* oldROI_saved = tfs->make<TH1F>(*oldROI);
          TH1F* newROI_saved = tfs->make<TH1F>(*newROI);
          mf::LogVerbatim("TestWireModifier")
            << "Wrote histograms " << oldROI_saved->GetName() << " and " << newROI_saved->GetName();
        }

        new_rois.add_range(roi_properties.begin, modified_data);
      }

      new_wires->emplace_back(new_rois, wire.Channel(), wire.View());

      auto const& rd_ptrs = digit_assn.at(i_w);
      for(auto const& rd_ptr : rd_ptrs)
        util::CreateAssn(*this, evt, *new_wires, rd_ptr, *new_digit_assn, new_wires->size() - 1);
    }

    evt.put(std::move(new_wires));
    evt.put(std::move(new_digit_assn));
  }
  DEFINE_ART_MODULE(TestWireModifier)
} // end namespace
