////////////////////////////////////////////////////////////////////////
// FilterMCTruthVolume_module.cc
//
// art::EDFilter that selects events containing at least one MCTruth
// neutrino interaction with a vertex __outside__ a user-specified box.
// Produces filtered MCTruth, GTruth, and MCFlux collections with
// their associations. Passes through POT summary at the SubRun level.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include <array>
#include <memory>
#include <vector>

class FilterMCTruthVolume : public art::EDFilter {
public:
  explicit FilterMCTruthVolume(fhicl::ParameterSet const& pset);

  bool filter(art::Event& evt) override;
  bool endSubRun(art::SubRun& sr) override;

private:
  art::InputTag fMCTruthLabel;
  art::InputTag fPOTLabel;
  std::array<double, 3> fVolumeLow;
  std::array<double, 3> fVolumeHigh;

  bool inVolume(double x, double y, double z) const;
};

// ---------------------------------------------------------------------------
FilterMCTruthVolume::FilterMCTruthVolume(fhicl::ParameterSet const& pset)
  : EDFilter{pset}
  , fMCTruthLabel{pset.get<art::InputTag>("MCTruthLabel", "generator")}
  , fPOTLabel{pset.get<art::InputTag>("POTLabel", "generator")}
{
  auto lo = pset.get<std::vector<double>>("VolumeLow");
  auto hi = pset.get<std::vector<double>>("VolumeHigh");
  if (lo.size() != 3 || hi.size() != 3)
    throw art::Exception(art::errors::Configuration)
      << "VolumeLow and VolumeHigh must each have exactly 3 elements (x, y, z).";
  std::copy(lo.begin(), lo.end(), fVolumeLow.begin());
  std::copy(hi.begin(), hi.end(), fVolumeHigh.begin());

  produces<std::vector<simb::MCTruth>>();
  produces<std::vector<simb::GTruth>>();
  produces<std::vector<simb::MCFlux>>();
  produces<art::Assns<simb::MCTruth, simb::GTruth>>();
  produces<art::Assns<simb::MCTruth, simb::MCFlux>>();
  produces<sumdata::POTSummary, art::InSubRun>();
}

// ---------------------------------------------------------------------------
bool FilterMCTruthVolume::inVolume(double x, double y, double z) const
{
  return x >= fVolumeLow[0] && x <= fVolumeHigh[0] &&
         y >= fVolumeLow[1] && y <= fVolumeHigh[1] &&
         z >= fVolumeLow[2] && z <= fVolumeHigh[2];
}

// ---------------------------------------------------------------------------
bool FilterMCTruthVolume::filter(art::Event& evt)
{
  auto outMCTruth = std::make_unique<std::vector<simb::MCTruth>>();
  auto outGTruth  = std::make_unique<std::vector<simb::GTruth>>();
  auto outMCFlux  = std::make_unique<std::vector<simb::MCFlux>>();
  auto assnsTG    = std::make_unique<art::Assns<simb::MCTruth, simb::GTruth>>();
  auto assnsTF    = std::make_unique<art::Assns<simb::MCTruth, simb::MCFlux>>();

  auto const& mctruthHandle = evt.getValidHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);

  art::FindOneP<simb::GTruth> findGTruth(mctruthHandle, evt, fMCTruthLabel);
  art::FindOneP<simb::MCFlux> findMCFlux(mctruthHandle, evt, fMCTruthLabel);

  art::PtrMaker<simb::MCTruth> makeTruthPtr(evt);
  art::PtrMaker<simb::GTruth>  makeGTruthPtr(evt);
  art::PtrMaker<simb::MCFlux>  makeMCFluxPtr(evt);

  for (size_t i = 0; i < mctruthHandle->size(); ++i) {
    auto const& mct = mctruthHandle->at(i);

    if (!mct.NeutrinoSet()) continue;

    double vx = mct.GetNeutrino().Nu().Vx();
    double vy = mct.GetNeutrino().Nu().Vy();
    double vz = mct.GetNeutrino().Nu().Vz();

    if (inVolume(vx, vy, vz)) continue;

    outMCTruth->push_back(mct);
    size_t idx = outMCTruth->size() - 1;

    if (findGTruth.isValid()) {
      auto const& gtp = findGTruth.at(i);
      if (gtp.isNonnull()) {
        outGTruth->push_back(*gtp);
        assnsTG->addSingle(makeTruthPtr(idx), makeGTruthPtr(outGTruth->size() - 1));
      }
    }

    if (findMCFlux.isValid()) {
      auto const& mfp = findMCFlux.at(i);
      if (mfp.isNonnull()) {
        outMCFlux->push_back(*mfp);
        assnsTF->addSingle(makeTruthPtr(idx), makeMCFluxPtr(outMCFlux->size() - 1));
      }
    }
  }

  bool pass = !outMCTruth->empty();

  evt.put(std::move(outMCTruth));
  evt.put(std::move(outGTruth));
  evt.put(std::move(outMCFlux));
  evt.put(std::move(assnsTG));
  evt.put(std::move(assnsTF));

  return pass;
}

// ---------------------------------------------------------------------------
bool FilterMCTruthVolume::endSubRun(art::SubRun& sr)
{
  auto const& potHandle = sr.getValidHandle<sumdata::POTSummary>(fPOTLabel);
  sr.put(std::make_unique<sumdata::POTSummary>(*potHandle), art::subRunFragment());
  return true;
}

DEFINE_ART_MODULE(FilterMCTruthVolume)
