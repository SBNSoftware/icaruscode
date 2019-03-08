////////////////////////////////////////////////////////////////////////
// Class:       OpDetWFDump
// Plugin Type: analyzer (art v3_01_01)
// File:        OpDetWFDump_module.cc
//
// Generated at Tue Feb 12 06:49:35 2019 by Kazuhiro Terao using cetskelgen
// from cetlib version v3_05_01.
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

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpHit.h"
#include <TTree.h>
#include <TFile.h>

class OpDetWFDump;

class OpDetWFDump : public art::EDAnalyzer {
public:
  explicit OpDetWFDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpDetWFDump(OpDetWFDump const&) = delete;
  OpDetWFDump(OpDetWFDump&&) = delete;
  OpDetWFDump& operator=(OpDetWFDump const&) = delete;
  OpDetWFDump& operator=(OpDetWFDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
private:

  // Declare member data here.
  TFile* _f;
  TTree* _wftree;
  TTree* _hittree;
  std::vector<float> _wf;
  int _ch;
  double _tstart;
  double _tend;
  double _tpeak;
  double _amp;
  double _area;
};


OpDetWFDump::OpDetWFDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, _wftree{nullptr}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _wftree = nullptr;
}

void OpDetWFDump::beginJob()
{
  if(!_wftree) {
    _f = TFile::Open("opdetwf.root","RECREATE");
    
    _wftree = new TTree("anatree","anatree");
    //_wftree->Branch("wf","std::vector<float>",&_wf);
    _wftree->Branch("wf",&_wf);
    _wftree->Branch("ch",&_ch,"ch/I");
    _wftree->Branch("tstart",&_tstart,"tstart/D");

    _hittree = new TTree("hittree","hittree");
    _hittree->Branch("ch",&_ch,"ch/I");
    _hittree->Branch("tpeak",&_tpeak,"tpeak/D");
    _hittree->Branch("amp",&_amp,"amp/D");
    _hittree->Branch("area",&_area,"area/D");
  }
  std::cout<<_wftree<<std::endl;
}

void OpDetWFDump::endJob()
{
  if(_wftree) {
    _f->cd();
    _wftree->Write();
    _hittree->Write();
    _f->Close();
  }
}

void OpDetWFDump::analyze(art::Event const& e)
{
  auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
  // Implementation of required member function here.
  art::Handle< std::vector< raw::OpDetWaveform > > hwf;
  e.getByLabel("opdaq",hwf);

  if(!hwf.isValid()){
    std::cout << "Invalid producer..." << std::endl;
  }

  std::cout << "Number of photon channel: " <<hwf->size() << std::endl;
  //  unsigned int nChannels = pmtHandle->size();
  auto geop = lar::providerFrom<geo::Geometry>(); 

  std::cout << geop->Ncryostats() << " cryostats" << std::endl;

  for(auto const& opwf : (*hwf)) {
    std::cout << opwf.ChannelNumber() << " ... " << opwf.size() << std::endl;
    _wf.resize(opwf.size());
    //std::cout<<_wftree<<std::endl;
    for(size_t i=0; i<_wf.size(); ++i) { _wf.at(i) = opwf.at(i); }
    _ch=opwf.ChannelNumber();
    _tstart=opwf.TimeStamp();
    _wftree->Fill();
  }

  std::cout << (*hwf)[0].TimeStamp() << " ... " << ts->TriggerTime() << " ... " << ts->BeamGateTime() << " ... " << ts->G4ToElecTime(0.) << std::endl;

  art::Handle< std::vector< recob::OpHit> > hhit;
  e.getByLabel("opreco",hhit);

  if(!hhit.isValid())
    std::cout << "aho" << std::endl;

  for(auto const& oph : (*hhit)) {
    _ch=oph.OpChannel();
    _tpeak = oph.PeakTime();
    _amp   = oph.Amplitude();
    _area  = oph.Area();
    _hittree->Fill();
  }
  
}

DEFINE_ART_MODULE(OpDetWFDump)
