////////////////////////////////////////////////////////////////////////
// Class:       MyFirstAnalyzer
// Plugin Type: analyzer (art v3_05_00)
// File:        MyFirstAnalyzer_module.cc
//
// Generated at Thu Apr  9 12:46:35 2020 by Bruce Howard using cetskelgen
// from cetlib version v3_10_00.
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

// Added for this analyzer
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "TH1D.h"

class MyFirstAnalyzer;


class MyFirstAnalyzer : public art::EDAnalyzer {
public:
  explicit MyFirstAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MyFirstAnalyzer(MyFirstAnalyzer const&) = delete;
  MyFirstAnalyzer(MyFirstAnalyzer&&) = delete;
  MyFirstAnalyzer& operator=(MyFirstAnalyzer const&) = delete;
  MyFirstAnalyzer& operator=(MyFirstAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

  art::InputTag fHitsLabel; // art label for hits info
  TH1D* histHitSummedADC;
};


MyFirstAnalyzer::MyFirstAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fHitsLabel (p.get<art::InputTag>("HitsLabel","gaushitTPC0")) // default gaushit in tpc0
{
  // Following example as in Andy Mastbaum's SBN LArSoft tutorial, we'll set up the
  // TFileService to draw the plot. Are other ways of storing data for later, for
  // example writing a tree with branches for later analysis
  art::ServiceHandle<art::TFileService> tfs;
  histHitSummedADC = tfs->make<TH1D>("histHitSummedADC",";Summed ADC;",100,0,4000);
}

void MyFirstAnalyzer::analyze(art::Event const& e)
{
  // Follow similar example to Andy Mastbaum's SBN LArSoft tutorial, let's just use 
  // GetByLabel. Are numerous ways to access the data in principle...
  art::Handle< std::vector<recob::Hit> > hits;
  e.getByLabel(fHitsLabel,hits);

  // Loop through the hits in the event and fill the histogram with SummedADC
  for( auto const& hit : *hits ){
    histHitSummedADC->Fill(hit.SummedADC());
  }
}

DEFINE_ART_MODULE(MyFirstAnalyzer)
