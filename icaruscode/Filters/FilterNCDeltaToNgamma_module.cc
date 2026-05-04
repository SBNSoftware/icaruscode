////////////////////////////////////////////////////////////////////////
/// \file FilterNCDeltaToNgamma_module.cc
/// \brief Simple filter to select only NC Delta to Ng events
/// 
/// \author hhausner@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef FILTER_FILTERNCDELTATONGAMMA_H
#define FILTER_FILTERNCDELTATONGAMMA_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"

// LArSoft Includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// C++ Includes
#include <string>

namespace ncres
{
  constexpr uint64_t pdgHash(const int& pdg1, const int& pdg2)
  {
    auto [low, high] = std::minmax({pdg1, pdg2});
    return (static_cast<uint64_t>(static_cast<uint32_t>(low)) << 32) | static_cast<uint32_t>(high);
  }
  constexpr uint64_t hash_1g0p = pdgHash(22, 2112);
  constexpr uint64_t hash_1g1p = pdgHash(22, 2212);

  class FilterNCDeltaToNgamma : public art::EDFilter
  {
    public:
      explicit FilterNCDeltaToNgamma(const fhicl::ParameterSet& pset);
      virtual void reconfigure(const fhicl::ParameterSet&);
      bool isNeutrino(const simb::MCParticle&);
      bool isDelta(const simb::MCParticle&);
      bool isNg(const std::vector<int>&);
      std::string PDGToName(const int&);
      bool filter(art::Event&);
    private:
      std::string fGenieModuleLabel;
  };

  //-----------------------------------------------------------------------
  // Constructor
  FilterNCDeltaToNgamma::FilterNCDeltaToNgamma(const fhicl::ParameterSet& pset) : EDFilter{pset}
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Reconfigure
  void FilterNCDeltaToNgamma::reconfigure(const fhicl::ParameterSet& pset)
  {
    fGenieModuleLabel = pset.get<std::string>("GenieModuleLabel");
  }

  //-----------------------------------------------------------------------
  // Is Neutrino?
  bool FilterNCDeltaToNgamma::isNeutrino(const simb::MCParticle& particle)
  {
    const int pdg = particle.PdgCode();
    return std::abs(pdg) == 12 || // nu-e
           std::abs(pdg) == 14 || // nu-mu
           std::abs(pdg) == 16  ; // nu-tau
  }

  //-----------------------------------------------------------------------
  // Is Delta?
  bool FilterNCDeltaToNgamma::isDelta(const simb::MCParticle& particle)
  {
    const int pdg = particle.PdgCode();
    return pdg == 1114 || // Delta-
           pdg == 2114 || // Delta0
           pdg == 2214 || // Delta+
           pdg == 2224  ; // Delta++
  }

  //-----------------------------------------------------------------------
  // Is Ng?
  bool FilterNCDeltaToNgamma::isNg(const std::vector<int>& daughters)
  {
    if (daughters.size() != 2)
      return false;
    uint64_t hash = pdgHash(daughters[0], daughters[1]);
    return (hash == hash_1g1p) || (hash == hash_1g0p);

  }

  //-----------------------------------------------------------------------
  // PDG to Name
  std::string FilterNCDeltaToNgamma::PDGToName(const int& pdg)
  {
    std::string pdgStr;
    switch (pdg)
    {
      case 1114:
        pdgStr = "Delta-";
        break;
      case 2114:
        pdgStr = "Delta0";
        break;
      case 2214:
        pdgStr = "Delta+";
        break;
      case 2224:
        pdgStr = "Delta++";
        break;
      case 2112:
        pdgStr = "neutron";
        break;
      case 2212:
        pdgStr = "proton";
        break;
      case 111:
        pdgStr = "pi0";
        break;
      case 211:
        pdgStr = "pi+";
        break;
      case -211:
        pdgStr = "pi-";
        break;
      case 22:
        pdgStr = "photon";
        break;
      default:
        pdgStr = std::to_string(pdg);
        break;
    }
    return pdgStr;
  }

  //-----------------------------------------------------------------------
  // Filter
  bool FilterNCDeltaToNgamma::filter(art::Event& evt)
  {
    art::Handle<std::vector<simb::MCTruth>> mcHandle;
    evt.getByLabel(fGenieModuleLabel, mcHandle);
    art::Ptr<simb::MCTruth> mc(mcHandle, 0);

    int delta_idx = std::numeric_limits<int>::min();
    int deltaPDG = std::numeric_limits<int>::min();
    std::vector<int> daughterPDG;
    for (int particle_idx = 0; particle_idx < mc->NParticles(); ++particle_idx)
    {
      const simb::MCParticle& particle = mc->GetParticle(particle_idx);
      if (delta_idx == std::numeric_limits<int>::min())
      {
        if (particle.Mother() == 0 && not isNeutrino(particle))
          return false;
        if (isDelta(particle))
        {
          delta_idx = particle_idx;
          deltaPDG = particle.PdgCode();
        }
      } else {
        if (particle.Mother() == delta_idx)
          daughterPDG.push_back(particle.PdgCode());
      }
    }

    mf::LogVerbatim("FilterNCDeltaToNgamma")
      << "╔════════════";
    if (daughterPDG.size() != 0)
    {
      mf::LogVerbatim("FilterNCDeltaToNgamma")
        << "║ " << PDGToName(deltaPDG);
      for (auto const& pdg : daughterPDG)
        mf::LogVerbatim("FilterNCDeltaToNgamma")
          << "║  ↳ " << PDGToName(pdg);
    } else
    {
      mf::LogVerbatim("FilterNCDeltaToNgamma")
        << "║ No Delta";
    }
    mf::LogVerbatim("FilterNCDeltaToNgamma")
      << "╚════════════";

    return (daughterPDG.size() != 0) && isNg(daughterPDG);
  }

  DEFINE_ART_MODULE(FilterNCDeltaToNgamma)

} // namespace ncres

#endif // FILTER_FILTERNCDELTATONGAMMA_H
