////////////////////////////////////////////////////////////////////////
// Class:       InputHitAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        InputHitAnalysis_module.cc
//
// Generated at Tue Feb 25 12:46:04 2025 by Bruce Howard using cetskelgen
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

#include "lardataobj/RecoBase/Hit.h"

// Trying to follow the idea of the HitEfficiencyAna_module to an extent and write a tree
#include "art_root_io/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

class InputHitAnalysis;


class InputHitAnalysis : public art::EDAnalyzer {
public:
  explicit InputHitAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  InputHitAnalysis(InputHitAnalysis const&) = delete;
  InputHitAnalysis(InputHitAnalysis&&) = delete;
  InputHitAnalysis& operator=(InputHitAnalysis const&) = delete;
  InputHitAnalysis& operator=(InputHitAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  TTree* fTree;
  int run;
  int subrun;
  int event;
  int eventIterator;

  float chargeSum;
  float chargeInt;
  float chargeAmp;
  float chargeTime;
  float chargeWidth;
  float chargeMult;
  int chargeCryo;
  int chargeTPC;
  int chargePlane;
  int chargeWire;
  float chargeGoodness;
  float chargeDOF;
};


InputHitAnalysis::InputHitAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree","HitsTree"); // as in hit efficiency analyzer

  fTree->Branch("run",&run);
  fTree->Branch("subrun",&subrun);
  fTree->Branch("event",&event);
  fTree->Branch("eventIterator",&eventIterator);
  fTree->Branch("chargeInt",&chargeInt);
  fTree->Branch("chargeAmp",&chargeAmp);
  fTree->Branch("chargeSum",&chargeSum);
  fTree->Branch("chargeTime",&chargeTime);
  fTree->Branch("chargeWidth",&chargeWidth);
  fTree->Branch("chargeMult",&chargeMult);
  fTree->Branch("chargeCryo",&chargeCryo);
  fTree->Branch("chargeTPC",&chargeTPC);
  fTree->Branch("chargePlane",&chargePlane);
  fTree->Branch("chargeWire",&chargeWire);
  fTree->Branch("chargeGoodness",&chargeGoodness);
  fTree->Branch("chargeDOF",&chargeDOF);

  eventIterator=-1;
}

void InputHitAnalysis::analyze(art::Event const& e)
{
  // --> Let's make a flattened tree of input hit characteristics:
  // - R/S/E >> the count of these per event is then the input number of hits we have
  // ---> need some other identifier for MC here...
  // - hit charge
  // - hit peak time
  // - hit width
  // - multiplicity
  // - Cryo/TPC/Plane
  // - wireID / ChannelID
  // - GoodnessOfFit / DegreesOfFreedom

  eventIterator+=1;

  // Get the event info
  run = e.id().run();
  subrun = e.id().subRun();
  event = e.id().event();

  // Get the input hit list for the event --> as elsewhere
  art::Handle< std::vector<recob::Hit> > hitsHandle;
  std::vector< art::Ptr<recob::Hit> > hits;
  if ( e.getByLabel("harps",hitsHandle) ) {
    art::fill_ptr_vector(hits,hitsHandle);
  }

  // Loop the hits
  for ( auto const& hit : hits ) {
    // Save the hit products
    chargeInt = hit->Integral();
    chargeSum = hit->SummedADC();
    chargeAmp = hit->PeakAmplitude();
    chargeTime = hit->PeakTime();
    chargeWidth = hit->RMS();
    chargeMult = hit->Multiplicity();
    chargeCryo = hit->WireID().Cryostat;
    chargeTPC = hit->WireID().TPC;
    chargePlane = hit->WireID().Plane;
    chargeWire = hit->WireID().Wire;
    chargeGoodness = hit->GoodnessOfFit();
    chargeDOF = hit->DegreesOfFreedom();
    // And write the output
    fTree->Fill();
  }

}

DEFINE_ART_MODULE(InputHitAnalysis)
