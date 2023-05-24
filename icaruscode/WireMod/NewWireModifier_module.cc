// std inlcudes
#include <string>
#include <vector>

// ROOT includes
#include "TH1D.h"

// art includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataobj/RecoBase/Wire.h"

//namespace
namespace wiremod
{

  class NewWireMod : public art::ReplicatedProducer
  {
    public:
      explicit NewWireMod(fhicl::ParameterSet const& pseti, art::ProcessingFrame const& frame);
      ~NewWireMod() override = default;
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt, art::ProcessingFrame const&) override;

    private:
      std::vector<art::InputTag> fWireLabelVec; // what wires are we modifying

      const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>(); // get the geometry

  }; // end NewWireMod class

  DEFINE_ART_MODULE(NewWireMod)

  //------------------------------------------------
  NewWireMod::NewWireMod(fhicl::ParameterSet const& pset, art::ProcessingFrame const& frame) : art::ReplicatedProducer(pset, frame)
  {
    // get the fhicl parameters
    this->reconfigure(pset);

    // each input has an output
    for (const auto& wireLabel : fWireLabelVec)
      produces<std::vector<recob::Wire>>(wireLabel.instance());

    mf::LogVerbatim("WireMod")
      << "looks like you made it this far...";
  }

  //------------------------------------------------
  void NewWireMod::reconfigure(fhicl::ParameterSet const& pset)
  {
    fWireLabelVec = pset.get<std::vector<art::InputTag>>("WireLabelVec");
  }

  //------------------------------------------------
  void NewWireMod::produce(art::Event& evt, art::ProcessingFrame const& frame)
  {
  }

} // end namespace
