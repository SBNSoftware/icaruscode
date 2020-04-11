////////////////////////////////////////////////////////////////////////
// Class:       PMTOpFlash
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTOpFlash_module.cc
//
// OpFlash Demo
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

// For optical detector
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

class PMTOpFlash;


class PMTOpFlash : public art::EDAnalyzer {
public:
  explicit PMTOpFlash(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTOpFlash(PMTOpFlash const&) = delete;
  PMTOpFlash(PMTOpFlash&&) = delete;
  PMTOpFlash& operator=(PMTOpFlash const&) = delete;
  PMTOpFlash& operator=(PMTOpFlash&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  art::InputTag fOpFlashModuleLabel;

  TTree *flashtree;

  int NFlash;
  std::vector<double> flashTime;
  std::vector<double> flashTimeWidth;
  std::vector<double> flashYCenter;
  std::vector<double> flashYWidth;
  std::vector<double> flashZCenter;
  std::vector<double> flashZWidth;
  std::vector<double> flashTotalPE;
};


PMTOpFlash::PMTOpFlash(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fOpFlashModuleLabel (p.get<art::InputTag>("OpFlashModuleLabel","opflashTPC0"))
{
  art::ServiceHandle<art::TFileService> tfs;

  flashtree = tfs->make<TTree>("opflashtree","tree for optical flashes");
  flashtree->Branch("NFlash",&NFlash,"NFlash/I");
  flashtree->Branch("time","std::vector<double>",&flashTime);
  flashtree->Branch("timewidth","std::vector<double>",&flashTimeWidth);
  flashtree->Branch("ycenter","std::vector<double>",&flashYCenter);
  flashtree->Branch("zcenter","std::vector<double>",&flashZCenter);
  flashtree->Branch("ywidth","std::vector<double>",&flashYWidth);
  flashtree->Branch("zwidth","std::vector<double>",&flashZWidth);
  flashtree->Branch("totalpe","std::vector<double>",&flashTotalPE);
}

void PMTOpFlash::analyze(art::Event const& e)
{
  NFlash = 0;

  // Similar to uboone anaTree code for flashes
  art::Handle< std::vector<recob::OpFlash> > flashListHandle;
  std::vector< art::Ptr<recob::OpFlash> >    flashlist;
  if( e.getByLabel(fOpFlashModuleLabel,flashListHandle) )
      art::fill_ptr_vector(flashlist, flashListHandle);

  // Loop through the flashes
  NFlash = flashlist.size();
  for (unsigned int iFlash=0; iFlash < flashlist.size(); ++iFlash){
    flashTime.push_back( flashlist.at(iFlash)->Time() );
    flashTimeWidth.push_back( flashlist.at(iFlash)->TimeWidth() );
    flashYCenter.push_back( flashlist.at(iFlash)->YCenter() );
    flashYWidth.push_back( flashlist.at(iFlash)->YWidth() );
    flashZCenter.push_back( flashlist.at(iFlash)->ZCenter() );
    flashZWidth.push_back( flashlist.at(iFlash)->ZWidth() );
    flashTotalPE.push_back( flashlist.at(iFlash)->TotalPE() );
  }

  flashtree->Fill();

  flashTime.clear();
  flashTimeWidth.clear();
  flashYCenter.clear();
  flashYWidth.clear();
  flashZCenter.clear();
  flashZWidth.clear();
  flashTotalPE.clear();
}

DEFINE_ART_MODULE(PMTOpFlash)
