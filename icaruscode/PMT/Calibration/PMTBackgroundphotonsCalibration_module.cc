////////////////////////////////////////////////////////////////////////
// Class:       PMTBackgroundphotonsCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTBackgroundphotonsCalibration_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
//
//  mailto:ascarpell@bnl.gov
// 
// 
//  Fist step of the data preparation for the PMT equalization using 
//  the backgorund photons light. 
// 
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
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

#include "icaruscode/PMT/Calibration/CaloTools/Utils.h"

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"



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

  bool inTime( double const & time, std::vector<double> & bounds );

  void clean();

private:

  art::InputTag m_ophit_label;
  art::InputTag fTriggerLabel;


  TTree *m_ophit_ttree;

  std::map<unsigned int, TH1D*> hintegral;
  std::map<unsigned int, TH1D*> hamplitude;
  std::map<unsigned int, TH1D*> hpulses;

  int m_run;
  int m_subrun;
  int m_event;
  uint32_t m_timestamp;

  double adc_to_mV;
  double adc_to_pC;
  double m_threshold;
  std::vector< unsigned int > m_channel_mask;
  std::vector<double> m_filter_intime;

  art::ServiceHandle<art::TFileService> tfs;

  std::map<unsigned int, unsigned int> m_pulses_count;

  double echarge = 1.602176634; // In units of 10^-7 pC

};


//------------------------------------------------------------------------------


pmtcalo::PMTBackgroundphotonsCalibration::PMTBackgroundphotonsCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{

   m_ophit_label = pset.get<art::InputTag>("OpHitModule", "ophit");
   fTriggerLabel = pset.get<art::InputTag>("TriggerModule");

   m_threshold = pset.get<double>("AmplitudeThreshold");

   adc_to_mV = pset.get<double>("ADCmV");  // Conversion from ADC to mV
   adc_to_pC = pset.get<double>("ADCpC");  // Conversion from ADC to pC 

   m_channel_mask = pset.get< std::vector< unsigned int > >
                            ("ChannelMasks", std::vector< unsigned int >());

   m_filter_intime = pset.get< std::vector<double> >("FilterInTime", std::vector< double >()); // in us 


  // Add histogram bins and range in form of a triplet (hist low, hist high, binsize)
    
}


//------------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::beginJob()
{

  m_ophit_ttree = tfs->make<TTree>("event","Event metadata TTree");
  
  m_ophit_ttree->Branch("run", &m_run, "run/I" );
  m_ophit_ttree->Branch("subrun", &m_subrun, "subrun/I" );
  m_ophit_ttree->Branch("event", &m_event, "event/I" );
  m_ophit_ttree->Branch("timestamp", &m_timestamp, "timestamp/I" );

  auto const geop = lar::providerFrom<geo::Geometry>();
  const unsigned int nPMTs = geop->NOpChannels();

  char histname[100]; 
  char histtitle[100];
  
  for(unsigned int opch=0; opch<nPMTs; ++opch)
  {

    // Create a histogram for any valid PMT 
    if( hasChannel(opch, m_channel_mask) )
      continue;

    sprintf(histname, "hintegral%u", opch);
    sprintf(histtitle, "PMT: %u;Pulse charge [10^7 electrons]", opch);
    hintegral[opch] = tfs->make<TH1D>(histname, histtitle, 200, 0, 8);

    sprintf(histname, "hamplitude%u", opch);
    sprintf(histtitle, "PMT: %u;Pulse amplitude [mV]", opch);
    hamplitude[opch] = tfs->make<TH1D>( histname, histtitle, 200, 0, 120 );

    sprintf(histname, "hpulses%u", opch);
    sprintf(histtitle, "PMT: %u;", opch);
    hpulses[opch] = tfs->make<TH1D>( histname, histtitle, 20, 0, 20 );

    m_pulses_count[opch] = 0;

  }
  
}


//-----------------------------------------------------------------------------

bool pmtcalo::PMTBackgroundphotonsCalibration::inTime( 
                        double const & time, std::vector<double> & bounds ){
  return (time >= bounds[0]) & ( time <= bounds[1] );
}


//-----------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::analyze(art::Event const& event)
{

   // Timestamp of the first event 

   m_run = event.id().run();
   m_subrun = event.id().subRun();
   m_event = event.id().event();
   m_timestamp = event.time().timeHigh(); // We just need precision at the s level


   // Here we filter out triggers that are not MinBias
   art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
   event.getByLabel( fTriggerLabel, trigger_handle );

   bool _is_minbias=true;

   if( trigger_handle.isValid() ) {

        _is_minbias = value( trigger_handle->triggerType ) ==1;
    
    }
      else{
         mf::LogError("pmtcalo::PMTBackgroundphotonsCalibration") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ; 
    }

   // First thing we sort the ophit in their respective channels
   art::Handle< std::vector< recob::OpHit > > ophitHandle;
   event.getByLabel(m_ophit_label, ophitHandle);

   
   for( auto const& ophit : (*ophitHandle) )
   {
      
      unsigned int opch = ophit.OpChannel();

      // Remove OpHits from known channels that are not working + OpHits that are not
      // In the trigger window
      if( hasChannel(opch, m_channel_mask) | 
            !inTime( ophit.StartTime(), m_filter_intime ) | !_is_minbias )
        continue;

      hintegral[opch]->Fill( ophit.Area()*adc_to_pC/echarge );
      hamplitude[opch]->Fill( ophit.Amplitude()*adc_to_mV );

      if( ophit.Amplitude()*adc_to_mV >= m_threshold )
        m_pulses_count[opch]++;

   }

   for( auto counts : m_pulses_count )
   {
    unsigned int opch = counts.first;
    hpulses[opch]->Fill( counts.second );
   }

   m_ophit_ttree->Fill();
   
   clean();
  
} // end analyze


//-----------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::clean()
{

  m_pulses_count.clear();

}


DEFINE_ART_MODULE(pmtcalo::PMTBackgroundphotonsCalibration)
