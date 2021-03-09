////////////////////////////////////////////////////////////////////////
// Class:       PMTBackgroundphotonsCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTBackgroundphotonsCalibration_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
//
//  mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "TTree.h"

namespace pmtcalo {
  class PMTBackgroundphotonsCalibration;
}


class pmtcalo::PMTBackgroundphotonsCalibration : public art::EDAnalyzer {

public:

  explicit PMTBackgroundphotonsCalibration(fhicl::ParameterSet const& pset);

  PMTBackgroundphotonsCalibration(PMTBackgroundphotonsCalibration const&) = delete;
  PMTBackgroundphotonsCalibration(PMTBackgroundphotonsCalibration&&) = delete;
  PMTBackgroundphotonsCalibration& operator=(PMTBackgroundphotonsCalibration const&) = delete;
  PMTBackgroundphotonsCalibration& operator=(PMTBackgroundphotonsCalibration&&) = delete;

  virtual void beginJob() override;

  void analyze(art::Event const& event) override;

  void clean();

private:

  art::InputTag m_ophit_label;


  TTree *m_ophit_ttree;

  int m_run;
  int m_subrun;
  int m_event;
  
  std::vector<int>   *m_channel_id = NULL;
  std::vector<float> *m_tstart = NULL;
  std::vector<float> *m_tmax = NULL;
  std::vector<float> *m_amplitude = NULL;
  std::vector<float> *m_integral = NULL;

  art::ServiceHandle<art::TFileService> tfs;

};


//------------------------------------------------------------------------------


pmtcalo::PMTBackgroundphotonsCalibration::PMTBackgroundphotonsCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{
   m_ophit_label = pset.get<art::InputTag>("OpHitModule", "ophit");
}


//------------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::beginJob()
{

  m_ophit_ttree = tfs->make<TTree>("pulsetree","OpHit TTree");
  
  m_ophit_ttree->Branch("run", &m_run, "run/I" );
  m_ophit_ttree->Branch("subrun", &m_subrun, "subrun/I" );
  m_ophit_ttree->Branch("event", &m_event, "event/I" );
  m_ophit_ttree->Branch("channel_id", &m_channel_id);
  m_ophit_ttree->Branch("tstart", &m_tstart );
  m_ophit_ttree->Branch("tmax", &m_tmax );
  m_ophit_ttree->Branch("amplitude", &m_amplitude );
  m_ophit_ttree->Branch("integral", &m_integral );
  
}


//-----------------------------------------------------------------------------



void pmtcalo::PMTBackgroundphotonsCalibration::analyze(art::Event const& event)
{

   m_run = event.id().run();
   m_subrun = event.id().subRun();
   m_event = event.id().event();

   // First thing we sort the ophit in their respective channels
   art::Handle< std::vector< recob::OpHit > > ophitHandle;
   event.getByLabel(m_ophit_label, ophitHandle);

   std::map<int, std::vector<recob::OpHit>> ophitch;
   for( auto const& ophit : (*ophitHandle) ){
      
      int ch = ophit.OpChannel();
      m_channel_id->push_back( ch );
      m_tstart->push_back( ophit.PeakTimeAbs() );
      m_tmax->push_back( ophit.PeakTime() ); 
      m_integral->push_back( ophit.Area() );
      m_amplitude->push_back( ophit.Amplitude() );

   }

   m_ophit_ttree->Fill();
   
   clean();

   ophitch.clear();
  
} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTBackgroundphotonsCalibration::clean(){

  m_channel_id->clear();
  m_tstart->clear();
  m_tmax->clear();
  m_amplitude->clear();
  m_integral->clear();

}


DEFINE_ART_MODULE(pmtcalo::PMTBackgroundphotonsCalibration)
