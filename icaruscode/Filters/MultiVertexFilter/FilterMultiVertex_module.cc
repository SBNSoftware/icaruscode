////////////////////////////////////////////////////////////////////////
// Class:       FilterMultiVertex
// Plugin Type: filter (Unknown Unknown)
// File:        FilterMultiVertex_module.cc
//
// Generated at Mon Jan  5 15:27:14 2026 by Mattia Sotgia using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Misc includes (?)
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft Includes
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"

#include <memory>
#include <map>
#include <vector>
#include <string>
#include <sstream>

namespace simfilter {
  class FilterMultiVertex;
}


class simfilter::FilterMultiVertex : public art::EDFilter {
public:
  explicit FilterMultiVertex(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterMultiVertex(FilterMultiVertex const&) = delete;
  FilterMultiVertex(FilterMultiVertex&&) = delete;
  FilterMultiVertex& operator=(FilterMultiVertex const&) = delete;
  FilterMultiVertex& operator=(FilterMultiVertex&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Parameters
  bool m_keepAll;

  // Data members
  int m_numberOfVertices;
  std::map<int, std::vector<int>> m_mapVertexNumberToEvents;
};

simfilter::FilterMultiVertex::FilterMultiVertex(fhicl::ParameterSet const& p)
  : EDFilter{p},
  m_keepAll{p.get<bool>("keep_all")},
  m_numberOfVertices{0}
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool simfilter::FilterMultiVertex::filter(art::Event& e)
{
  bool isOneInteraction(true);
  m_numberOfVertices = 0;
  auto allMCLists = e.getMany<std::vector<simb::MCTruth>>();

  for (size_t i = 0; i < allMCLists.size(); i++)
  {
    art::Handle<std::vector<simb::MCTruth>> mcListHandle = allMCLists.at(i);
    for (size_t j = 0; j < mcListHandle->size(); j++)
    {
      art::Ptr<simb::MCTruth> mcTruth(mcListHandle, j);
      for(int k = 0; k < mcTruth->NParticles(); k++)
      {
        auto p = mcTruth->GetParticle(k);
        if (p.Mother() < 0 && (std::abs(p.PdgCode()) == 12 || std::abs(p.PdgCode()) == 14))
        {
          m_numberOfVertices++;
        }
      }
    }
  }

  isOneInteraction = (m_numberOfVertices == 1) || m_keepAll;
  mf::LogDebug("FilterMultiVertex") << "This has number of vertices = " << m_numberOfVertices 
    << " and keep all is " << (m_keepAll ? "true":"false") << " and isOneInteraction is " 
    << (isOneInteraction ? "true":"false") << std::endl;
  m_mapVertexNumberToEvents[m_numberOfVertices].push_back(e.event());

  return isOneInteraction;
}//

void simfilter::FilterMultiVertex::beginJob()
{
}

void simfilter::FilterMultiVertex::endJob()
{
  for (auto const& [vertices, events]: m_mapVertexNumberToEvents)
  {
    std::cout << "Counted " << events.size() << " events with "
              << vertices << " vertices" << std::endl
              << "--> event(s): ";
    for (auto const& e: events)
    {
      std::cout << e << ", ";
    }
    std::cout << std::endl;
  }
}

DEFINE_ART_MODULE(simfilter::FilterMultiVertex)
