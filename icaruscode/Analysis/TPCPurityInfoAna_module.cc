////////////////////////////////////////////////////////////////////////
// Class:       TPCPurityInfoAna
// Plugin Type: analyzer (art v3_04_00)
// File:        TPCPurityInfoAna_module.cc
//
// Generated at Sun Jan 26 22:13:22 2020 by Wesley Ketchum using cetskelgen
// from cetlib version v3_09_00.
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

//purity info class
#include "icaruscode/IcarusObj/TPCPurityInfo.hh"

//output to ntuple
#include "art_root_io/TFileService.h"
#include "TNtuple.h"

namespace ana {
  class TPCPurityInfoAna;
}


class ana::TPCPurityInfoAna : public art::EDAnalyzer {
public:
  explicit TPCPurityInfoAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCPurityInfoAna(TPCPurityInfoAna const&) = delete;
  TPCPurityInfoAna(TPCPurityInfoAna&&) = delete;
  TPCPurityInfoAna& operator=(TPCPurityInfoAna const&) = delete;
  TPCPurityInfoAna& operator=(TPCPurityInfoAna&&) = delete;

  void analyze(art::Event const& e) override;
  void beginJob() override;

private:

  // Declare member data here.
  art::InputTag const fPurityInfoLabel;
  TNtuple* fPurityTuple;

  bool fPrintInfo;

};


ana::TPCPurityInfoAna::TPCPurityInfoAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  , fPurityInfoLabel(p.get<art::InputTag>("PurityInfoLabel"))
  , fPrintInfo(p.get<bool>("PrintInfo",true))
{
  consumes< std::vector<anab::TPCPurityInfo> >(fPurityInfoLabel);
}

void ana::TPCPurityInfoAna::beginJob()
{
    
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  fPurityTuple = tfs->make<TNtuple>("purityTuple","Purity Tuple","run:subrun:ev:cryo:tpc:wires:ticks:attNear:attFar:att:err");

}
void ana::TPCPurityInfoAna::analyze(art::Event const& e)
{

  if(fPrintInfo)
    std::cout << "Processing Run " << e.run() 
	      << ", Subrun " << e.subRun()
	      << ", Event " << e.event() 
	      << ":" << std::endl;

  art::Handle< std::vector<anab::TPCPurityInfo> > purityInfoHandle;
  e.getByLabel(fPurityInfoLabel,purityInfoHandle);
  auto const& purityInfoVec(*purityInfoHandle);

  if(fPrintInfo)
    std::cout << "\tThere are " << purityInfoVec.size() << " purity info objects in the event."
	      << std::endl;
  
  for(auto const& pinfo : purityInfoVec){
    if(fPrintInfo) pinfo.Print();
    fPurityTuple->Fill(e.run(),e.subRun(),e.event(),pinfo.Cryostat,pinfo.TPC,pinfo.Wires,pinfo.Ticks,pinfo.AttenuationNEAR,pinfo.AttenuationFAR,pinfo.Attenuation,pinfo.FracError);
  }


}

DEFINE_ART_MODULE(ana::TPCPurityInfoAna)
