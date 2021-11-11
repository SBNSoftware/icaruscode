////////////////////////////////////////////////////////////////////////
// Class:       DemoAna
// Plugin Type: analyzer
// File:        DemoAna_module.cc
//
// Generated at Thu Nov 11 13:28:52 2021 by Gray Putnam using cetskelgen
// from  version .
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

namespace icarus {
  class DemoAna;
}


class icarus::DemoAna : public art::EDAnalyzer {
public:
  explicit DemoAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DemoAna(DemoAna const&) = delete;
  DemoAna(DemoAna&&) = delete;
  DemoAna& operator=(DemoAna const&) = delete;
  DemoAna& operator=(DemoAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


icarus::DemoAna::DemoAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
}

void icarus::DemoAna::analyze(art::Event const& e)
{

  std::cout << "Processing Event: " << e.event() << " in run: " << e.run();

}

DEFINE_ART_MODULE(icarus::DemoAna)
