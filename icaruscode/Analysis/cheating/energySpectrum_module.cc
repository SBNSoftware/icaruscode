////////////////////////////////////////////////////////////////////////
// Class:       energySpectrum
// Plugin Type: analyzer (Unknown Unknown)
// File:        energySpectrum_module.cc
//
// Generated at Mon Feb 24 10:09:58 2025 by Mattia Sotgia using cetskelgen
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

namespace icaruscode {
  namespace pandoraCheating {
    class energySpectrum;
  }
}


class icaruscode::pandoraCheating::energySpectrum : public art::EDAnalyzer {
public:
  explicit energySpectrum(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  energySpectrum(energySpectrum const&) = delete;
  energySpectrum(energySpectrum&&) = delete;
  energySpectrum& operator=(energySpectrum const&) = delete;
  energySpectrum& operator=(energySpectrum&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


icaruscode::pandoraCheating::energySpectrum::energySpectrum(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void icaruscode::pandoraCheating::energySpectrum::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void icaruscode::pandoraCheating::energySpectrum::beginJob()
{
  // Implementation of optional member function here.
}

void icaruscode::pandoraCheating::energySpectrum::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(icaruscode::pandoraCheating::energySpectrum)
