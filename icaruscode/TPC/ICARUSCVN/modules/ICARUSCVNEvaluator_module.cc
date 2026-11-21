////////////////////////////////////////////////////////////////////////
// \file    LArCVNEvaluator_module.cc
// \brief   Producer module creating CVN neural net results
// \author  Alexander Radovic - a.radovic@gmail.com
//          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larrecodnn/CVN/interfaces/ITFNetHandler.h"
#include "icaruscode/TPC/ICARUSCVN/module_helpers/ICARUSPixelMapProducer.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/Result.h"

namespace lcvn {

class ICARUSCVNEvaluator : public art::EDProducer {
public:
  explicit ICARUSCVNEvaluator(fhicl::ParameterSet const &pset);

  void produce(art::Event &evt);

private:
  // /// Module label for input pixel maps
  std::string fModule;
  std::string fPixelMapInput;
  std::string fPixelMapModuleLabel;
  art::InputTag fSliceLabel;
  art::InputTag fPFParticleModuleLabel;
  art::InputTag fT0Label;
  std::string fCVNType;
  bool fVerbose;
  std::unique_ptr<lcvn::ITFNetHandler> fTFHandler;
  ICARUSPixelMapHitProducer fPMProducer;
};

ICARUSCVNEvaluator::ICARUSCVNEvaluator(fhicl::ParameterSet const &pset)
    : EDProducer{pset},
      fModule(pset.get<std::string>("Module")),
      fPixelMapInput(pset.get<std::string>("PixelMapInput")),
      fPixelMapModuleLabel(pset.get<std::string>("PixelMapModuleLabel")),
      fSliceLabel(pset.get<art::InputTag>("SliceLabel", "")),
      fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel", "")),
      fT0Label(pset.get<art::InputTag>("T0Label", "")),
      fCVNType(pset.get<std::string>("CVNType")),
      fVerbose(pset.get<bool>("verbose", false)),
      fTFHandler{art::make_tool<lcvn::ITFNetHandler>(pset.get<fhicl::ParameterSet>("ICARUSTFHandler"))},
      fPMProducer(pset.get<fhicl::ParameterSet>("PixelMapProducer"))
  {
  produces<std::vector<lcvn::Result>>();
  if (fSliceLabel != "") {
    produces<art::Assns<recob::Slice, lcvn::Result>>();
  } else {
    produces<art::Assns<lcvn::PixelMap, lcvn::Result>>();
  }
}
//......................................................................
void ICARUSCVNEvaluator::produce(art::Event &evt) {

  /// Define containers for the things we're going to produce
  std::unique_ptr<std::vector<Result>> resultCol(new std::vector<Result>);
  auto assn1 = std::make_unique<art::Assns<recob::Slice, lcvn::Result>>();
  auto assn2 = std::make_unique<art::Assns<lcvn::PixelMap, lcvn::Result>>();
  bool isSlice = false;
  // for (auto &cfg : fEvaluators) {
    if (fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow") {
      if (fSliceLabel != "") { // by default use slice information
        if (fVerbose) std::cout << "Using slice information" << std::endl;
        isSlice = true;
        std::vector<art::Ptr<recob::Slice>> slcList;
        auto slcHandle = evt.getHandle<std::vector<recob::Slice>>(fSliceLabel);
        if (!slcHandle.isValid()) {
          throw cet::exception("ICARUSCVNEvaluator")
              << "Unable to get slices using label " << fSliceLabel;
        } else {
          if (fVerbose) std::cout << "Filling slice" << std::endl;
          art::fill_ptr_vector(slcList, slcHandle);
        }
        art::Handle<std::vector<recob::PFParticle>> PFPListHandle;
        std::vector<art::Ptr<recob::PFParticle>> PFPList;
        if (evt.getByLabel(fPFParticleModuleLabel, PFPListHandle)) {
          if (fVerbose) std::cout << "Filling pfp" << std::endl;
          art::fill_ptr_vector(PFPList, PFPListHandle);
        }

        art::FindManyP<recob::Hit> findManyHits(slcHandle, evt,
                                                fSliceLabel);
        art::FindManyP<recob::PFParticle> findManyPFPs(
            slcHandle, evt, fPFParticleModuleLabel);
        art::FindManyP<anab::T0> findManyT0s(PFPListHandle, evt, fT0Label);
        auto const detProp =
            art::ServiceHandle<detinfo::DetectorPropertiesService const>()
                ->DataFor(evt);
        for (auto const &slice : slcList) {
          if (fVerbose)
            std::cout << "********* " << fModule << " ********* "
                      << evt.run() << "  " << evt.subRun() << "  "
                      << evt.id().event() << "  " << slice->ID()
                      << "  **************\n";
          std::vector<float> pfp_T0_vec;
          if (findManyPFPs.isValid()) {
            std::vector<art::Ptr<recob::PFParticle>> slicePFPs = findManyPFPs.at(slice.key());
            if (slicePFPs.size()) {
              for (auto const &pfp : slicePFPs) {
                if (findManyT0s.isValid()) {
                  std::vector<art::Ptr<anab::T0>> T0_vec =
                      findManyT0s.at(pfp.key());
                  if (T0_vec.size()) {
                    for (auto const &T0 : T0_vec) {
                      pfp_T0_vec.push_back(T0->Time());
                    }
                  }
                }
              }
            }
          }

          float min_T0 = 0.;
          if (pfp_T0_vec.size()) {
            min_T0 = *min_element(pfp_T0_vec.begin(), pfp_T0_vec.end());
            if (fVerbose) std::cout << "min t0: " << min_T0 << std::endl;
          }

          if (findManyHits.isValid()) {
            std::vector<art::Ptr<recob::Hit>> slicehits = findManyHits.at(slice.key());
            if (fVerbose) std::cout << "slice key: " << slice.key() << " slice id: " << slice.id() << std::endl;
            fPMProducer.Set_fT0_value(min_T0);
            PixelMap pm = fPMProducer.ICARUSCreateMap(detProp, slicehits);
            auto nhits = fPMProducer.TotHits();
            pm.SetTotHits(nhits);
            // pm.fSliceID = slice->ID();
            std::string moduleSlcId = std::string(fModule) + " - slice id: " + std::to_string(slice.id()) + " - slice key: " + std::to_string(slice.key() + 1);
            std::vector<std::vector<float>> output = fTFHandler->Predict(pm);
	    for (auto &vec : output) {
	      if (vec.size() >= 2) vec.insert(vec.begin() + 2, 0.0f);
	    }
            resultCol->emplace_back(output);
            util::CreateAssn(*this, evt, *resultCol, slice, *assn1);
	  } else {
            if (fVerbose) std::cout << "findManyHits: " << findManyHits.isValid() << std::endl;
          }
        }
      } else {
	mf::LogError("ICARUSCVNEvaluator::produce")
          << "No Slice Label. Exiting without processing events" << std::endl;
	return;
      }
    } else {
      mf::LogError("ICARUSCVNEvaluator::produce")
          << "CVN Type not in the allowed list: Tensorflow" << std::endl;
      mf::LogError("ICARUSCVNEvaluator::produce")
          << "Exiting without processing events" << std::endl;
      return;
    }
  // }
    if (fVerbose) {
      for (const auto &r : *resultCol) {
	for (const auto &row : r.fOutput) {
	  for (const auto &val : row) {
	    std::cout << val << " ";
	  }
	  std::cout << std::endl;
        }
      }
    }

    evt.put(std::move(resultCol));
    if (isSlice) {
      evt.put(std::move(assn1));
    } else {
      evt.put(std::move(assn2));
    }
}

DEFINE_ART_MODULE(lcvn::ICARUSCVNEvaluator)
} // namespace lcvn
////////////////////////////////////////////////////////////////////////
