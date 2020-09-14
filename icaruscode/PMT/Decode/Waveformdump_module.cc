////////////////////////////////////////////////////////////////////////
// Class:       Waveformdump
// Module Type: analyzer
// File:        Waveformdump_module.cc
// Description: Prints out information about each event.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"

#include "art_root_io/TFileService.h"
#include "TH1F.h"
#include "TNtuple.h"

#include <stdio.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <bitset>
#include <iostream>


namespace opdet {
  class Waveformdump;
}

/**************************************************************************************************/

class opdet::Waveformdump : public art::EDAnalyzer {

public:

  explicit Waveformdump(const fhicl::ParameterSet & pset);
  virtual ~Waveformdump();

  void analyze(const art::Event& evt) override;
  void beginJob() override;
  void endJob() override;

private:

  TTree *fEventTree;
  std::string fInputModule;

  pmtana::PulseRecoManager  fPulseRecoMgr;
  pmtana::PMTPulseRecoBase* fThreshAlg;
  pmtana::PMTPedestalBase*  fPedAlg;

  int fRun, fEvent;
  std::vector<int> fChannelVec;
  std::vector<std::vector<double>> fWvfmsVec;
  std::vector<std::vector<double>> fBaseline;
  std::vector<std::vector<double>> fBaselineRMS;

};


opdet::Waveformdump::Waveformdump(const fhicl::ParameterSet & pset): art::EDAnalyzer(pset)
{

  fInputModule = pset.get< std::string >("InputModule");

  // Initialize the hit finder algorithm
  auto const hit_alg_pset = pset.get< fhicl::ParameterSet >("HitAlgoPset");
  std::string threshAlgName = hit_alg_pset.get< std::string >("Name");
  if      (threshAlgName == "Threshold")
    fThreshAlg = new pmtana::AlgoThreshold(hit_alg_pset);
  else if (threshAlgName == "SiPM")
    fThreshAlg = new pmtana::AlgoSiPM(hit_alg_pset);
  else if (threshAlgName == "SlidingWindow")
    fThreshAlg = new pmtana::AlgoSlidingWindow(hit_alg_pset);
  else if (threshAlgName == "FixedWindow")
    fThreshAlg = new pmtana::AlgoFixedWindow(hit_alg_pset);
  else if (threshAlgName == "CFD" )
    fThreshAlg = new pmtana::AlgoCFD(hit_alg_pset);
  else throw art::Exception(art::errors::UnimplementedFeature)
                  << "Cannot find implementation for "
                  << threshAlgName << " algorithm.\n";

  auto const ped_alg_pset = pset.get< fhicl::ParameterSet >("PedAlgoPset");
  std::string pedAlgName = ped_alg_pset.get< std::string >("Name");
  if      (pedAlgName == "Edges")
    fPedAlg = new pmtana::PedAlgoEdges(ped_alg_pset);
  else if (pedAlgName == "RollingMean")
    fPedAlg = new pmtana::PedAlgoRollingMean(ped_alg_pset);
  else if (pedAlgName == "UB"   )
    fPedAlg = new pmtana::PedAlgoUB(ped_alg_pset);
  else throw art::Exception(art::errors::UnimplementedFeature)
                  << "Cannot find implementation for "
                  << pedAlgName << " algorithm.\n";

  fPulseRecoMgr.AddRecoAlgo(fThreshAlg);
  fPulseRecoMgr.SetDefaultPedAlgo(fPedAlg);

}

void opdet::Waveformdump::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;

  fEventTree = tfs->make<TTree>("events","waveform tree");
  fEventTree->Branch("fRun",&fRun,"fRun/I");
  fEventTree->Branch("fEvent",&fEvent,"fEvent/I");
  fEventTree->Branch("fChannelVec",&fChannelVec);
  fEventTree->Branch("fWvfmsVec",&fWvfmsVec);
  fEventTree->Branch("fBaseline",&fBaseline);
  fEventTree->Branch("fBaselineRMS",&fBaselineRMS);

}

void opdet::Waveformdump::endJob()
{
  delete fThreshAlg;
  delete fPedAlg;

  fWvfmsVec.clear();
  fChannelVec.clear();

  std::cout << "Ending Waveformdump...\n";
}


opdet::Waveformdump::~Waveformdump()
{
}


void opdet::Waveformdump::analyze(const art::Event& evt)
{

  fRun = evt.run();
  fEvent = evt.event();

  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  if ( !evt.getByLabel(fInputModule, wfHandle) ) {
    std::cerr << "Requested fragments with label : "
                               << fInputModule << "but none exist" << std::endl;
    return;
  }

  if (!wfHandle.isValid()) {
    std::cerr << "NOT VALD HANDLE!" << std::endl; return;
  }

  if (wfHandle.isValid()) {

    // Loop over the waveforms
    for(auto const& wf : *wfHandle) {

      // GetChannell
      fChannelVec.push_back( wf.ChannelNumber() );

      // Run pedestal subtracton
      fPulseRecoMgr.Reconstruct(wf);

      auto const  _ped_mean_v = fPedAlg->Mean();
      auto const  _ped_sigma_v = fPedAlg->Sigma();

      fBaseline.push_back(_ped_mean_v);
      fBaselineRMS.push_back(_ped_sigma_v);

      // Save waveform after pedestal removal
      std::vector<double> tmpwfm;
      for( size_t i=0; i< wf.size(); i++ )
        tmpwfm.push_back( wf[i] - _ped_mean_v[i] ) ;

      fWvfmsVec.push_back(tmpwfm);

    }

    fEventTree->Fill();

  }
}

DEFINE_ART_MODULE(opdet::Waveformdump)
