////////////////////////////////////////////////////////////////////////
// Class:       testCheated
// Plugin Type: analyzer (Unknown Unknown)
// File:        testCheated_module.cc
//
// Generated at Thu Feb 20 02:28:49 2025 by Mattia Sotgia using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include <TTree.h>

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "icaruscode/RecoUtils/RecoUtils.h"

namespace icaruscode {
  namespace pandoraCheating {
    class testCheated;
  }
}


class icaruscode::pandoraCheating::testCheated : public art::EDAnalyzer {
public:
  explicit testCheated(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  testCheated(testCheated const&) = delete;
  testCheated(testCheated&&) = delete;
  testCheated& operator=(testCheated const&) = delete;
  testCheated& operator=(testCheated&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Output TTree + branches
  // TTree* fTree;


  // FHiCL settings;
  std::string fPFParticleLabel, fPFParticleMetadataLabel;
  bool fVerbose;
  std::vector<std::string> fPandoraSuffixes;

};


icaruscode::pandoraCheating::testCheated::testCheated(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fPFParticleMetadataLabel(p.get<std::string>("PFParticleMetadataLabel")),
  fVerbose(p.get<bool>("VerboseOutput")),
  fPandoraSuffixes(p.get<std::vector<std::string>>("PandoraSuffixes"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void icaruscode::pandoraCheating::testCheated::analyze(art::Event const& e)
{
  
  // loops over two cryostats FUCKING ICARUS
  for (auto const &pandora_suffix: fPandoraSuffixes) {

    if (fVerbose)
      std::cout << "[LOOKing at " << pandora_suffix << "]" << std::endl;

    art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = 
      e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel + pandora_suffix);
    std::vector<art::Ptr<recob::PFParticle>> Pfps;
  
    if (pfpHandle.isValid())
      art::fill_ptr_vector(Pfps, pfpHandle);
  
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpHandle, e, fPFParticleMetadataLabel + pandora_suffix);

    if (Pfps.size() == 0) {
      if (fVerbose)
        std::cout << "No slice(s) found for event on " << pandora_suffix << std::endl;
      continue;
    }
  
    // int iPfp = 0;
    int iUnambiguousCosmics = 0;
    for (const art::Ptr<recob::PFParticle> &pfp: Pfps) {
      const larpandoraobj::PFParticleMetadata *primary_meta = pfpMetadataAssns.at(pfp.key()).at(0).get();
      // const larpandoraobj::PFParticleMetadata *primary_meta = pfpMetadataAssns.at(iPfp).at(0).get();
      auto const &properties = primary_meta->GetPropertiesMap();
      if (properties.count("IsClearCosmic")) {
        if (fVerbose)
          std::cout << "Unambiguous cosmic with PDG" << pfp->PdgCode() << std::endl;
        iUnambiguousCosmics++;
      }
      // iPfp++;
    }

    if (fVerbose && iUnambiguousCosmics > 0)
      std::cout << "LOOKing at event " << e.id().event() << " found " << iUnambiguousCosmics << " unambiguous cosmics" << std::endl;
    
    if (fVerbose && iUnambiguousCosmics == 0)
      std::cout << "LOOKing at event " << e.id().event() << " found no unambiguous cosmics" << std::endl;
  }
  
}

void icaruscode::pandoraCheating::testCheated::beginJob()
{
  // Implementation of optional member function here.
}

void icaruscode::pandoraCheating::testCheated::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(icaruscode::pandoraCheating::testCheated)
