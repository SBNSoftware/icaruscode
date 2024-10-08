////////////////////////////////////////////////////////////////////////
// Class:       PMTBackgroundphotonsCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTBackgroundphotonsCalibration_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
// Updated at Fri Mar 08 10:32:21 2024 by Matteo Vicenzi  
// 
// mailto: ascarpell@bnl.gov, mvicenzi@bnl.gov
// 
//  Fist step of the data preparation for the PMT equalization using 
//  the background photons light (offbeamminbias streams).
//  This modules relies on the OpHit objects from the standard optical
//  reconstruction flow.
// 
////////////////////////////////////////////////////////////////////////

#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"

#include "art_root_io/TFileService.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/OpHit.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/PMT/Calibration/CaloTools/Utils.h"

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

#include <vector>
#include <algorithm>

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
  bool inGate(double const & time, std::vector<double> & bounds);

  void analyze(art::Event const& event) override;

private:

  // input tags: ophits + trigger
  art::InputTag m_ophit_label;
  art::InputTag m_trigger_label;

  // output histograms
  TTree *m_ophit_tree; // all ophits
  std::map<unsigned int, TH1D*> hamplitude;  //ophit amplitude
  std::map<unsigned int, TH1D*> hintegral;   //ophit integral

  int m_run;
  int m_event;
  uint32_t m_timestamp;

  unsigned int m_channel;
  double m_time;
  double m_amplitude;
  double m_integral;

  // conversions
  double adc_to_mV;
  double adc_to_pC;
  double echarge = 1.602176634; // In units of 10^-7 pC

  // dead/broken channels to be skipped
  std::vector<unsigned int> m_channel_mask;
  // limit to select in-gate ophits
  std::vector<double> m_filter_ingate;

  // select isolated ophits  
  bool m_filter_intime;
  double m_time_window;

  // histogram binnings: low, high binsize
  std::vector<double> m_amp_binning;
  std::vector<double> m_integral_binning;

  art::ServiceHandle<art::TFileService> tfs;
  
  static bool sortByChannelTime(const recob::OpHit &a, const recob::OpHit &b){
    if (a.OpChannel() == b.OpChannel()) return a.StartTime() < b.StartTime(); 
    return a.OpChannel() < b.OpChannel();
  };

};


//------------------------------------------------------------------------------


pmtcalo::PMTBackgroundphotonsCalibration::PMTBackgroundphotonsCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{

   m_ophit_label = pset.get<art::InputTag>("OpHitModule", "ophit");
   m_trigger_label = pset.get<art::InputTag>("TriggerModule");\

   adc_to_mV = pset.get<double>("ADCmV");  // Conversion from ADC to mV
   adc_to_pC = pset.get<double>("ADCpC");  // Conversion from ADC to pC 

   m_channel_mask = pset.get< std::vector< unsigned int > >
                            ("ChannelMask", std::vector< unsigned int >());
   m_filter_ingate = pset.get< std::vector<double> >("FilterInGate", std::vector< double >()); // in us
   m_filter_intime = pset.get<bool>("FilterInTime", true); 
   m_time_window = pset.get<double>("TimeWindow", 0.2 ); //in us

   m_amp_binning = pset.get< std::vector<double> >("AmplitudeBinning");
   m_integral_binning = pset.get< std::vector<double> >("IntegralBinning");
    
}


//------------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::beginJob()
{

  // create full ophits tree for possible amplitude-cut selection
  // include event/timestamp data for later use
  m_ophit_tree = tfs->make<TTree>("ophits","ophits tree");
  
  m_ophit_tree->Branch("run", &m_run, "run/I" );
  m_ophit_tree->Branch("event", &m_event, "event/I" );
  m_ophit_tree->Branch("timestamp", &m_timestamp, "timestamp/I" );
  m_ophit_tree->Branch("channel", &m_channel, "channel/I" );
  m_ophit_tree->Branch("time", &m_time, "time/D" );
  m_ophit_tree->Branch("amplitude", &m_amplitude, "amplitude/D" );
  m_ophit_tree->Branch("integral", &m_integral, "integral/D" );

  const unsigned int nPMTs = art::ServiceHandle<geo::WireReadout const>()->Get().NOpChannels();

  char histname[100]; 
  char histtitle[100];

  double alow = m_amp_binning.at(0);
  double ahigh = m_amp_binning.at(1);
  double absize = m_amp_binning.at(2);
  int abins = int((ahigh-alow)/absize);
  
  double ilow = m_integral_binning.at(0);
  double ihigh = m_integral_binning.at(1);
  double ibsize = m_integral_binning.at(2);
  int ibins = int((ihigh-ilow)/ibsize);

  for(unsigned int opch=0; opch<nPMTs; ++opch)
  {

    // Create a histogram for any valid PMT 
    if( hasChannel(opch, m_channel_mask) ) continue;

    sprintf(histname, "hintegral%u", opch);
    sprintf(histtitle, "PMT: %u;Pulse charge [10^7 electrons]", opch);
    hintegral[opch] = tfs->make<TH1D>(histname, histtitle, ibins, ilow, ihigh);

    sprintf(histname, "hamplitude%u", opch);
    sprintf(histtitle, "PMT: %u;Pulse amplitude [mV]", opch);
    hamplitude[opch] = tfs->make<TH1D>( histname, histtitle, abins, alow, ahigh);
  
  }

}

//-----------------------------------------------------------------------------

bool pmtcalo::PMTBackgroundphotonsCalibration::inGate(double const & time, std::vector<double> & bounds ){
  return (time >= bounds[0]) && ( time <= bounds[1] );
}


//-----------------------------------------------------------------------------

void pmtcalo::PMTBackgroundphotonsCalibration::analyze(art::Event const& event)
{

   // event info: this is mainly to save the timestamps
   // the timestamp of the first event is used to sort files in the calibration db 
  
   m_run = event.id().run();
   m_event = event.id().event();
   m_timestamp = event.time().timeHigh(); // We just need precision at the s level

   // here we filter out triggers that are not MinBias
   // this modules is supposed to be run on (offbeam)(bnb/numi)minbias streams
   // or laser runs without laser (equivalent to minbias)
   
   art::Handle trigger_handle = event.getHandle<sbn::ExtraTriggerInfo>(m_trigger_label);
   
   if( trigger_handle.isValid() ) {
    
     bool _is_minbias = value( trigger_handle->triggerType ) == 1;
     if( !_is_minbias ) return;    

   } else{
     mf::LogError("pmtcalo::PMTBackgroundphotonsCalibration") << "No raw::Trigger associated to label: " << m_trigger_label.label();
   }

   // first thing we sort the ophit in their respective channels
   // creating the amplitude and integral histograms (w/o cuts)

   auto const & ophitCollection = event.getProduct<std::vector<recob::OpHit>>(m_ophit_label);

   // make sure they are sorted
   auto sortableOpHits = ophitCollection;
   std::sort( sortableOpHits.begin(), sortableOpHits.end(), sortByChannelTime );
 
   for (auto it = sortableOpHits.begin(); it != sortableOpHits.end(); ++it) {
    
     auto const &ophit = *it;  // Current element
     unsigned int opch = ophit.OpChannel();
      
     // remove OpHits from known channels that are not working
     if( hasChannel(opch, m_channel_mask) ) continue;
     // remove OpHits that are not within a subset of the beam gate (thus unbiased)
     if( !inGate( ophit.StartTime(), m_filter_ingate ) ) continue; 
   
     auto nextIt = std::next(it);
     auto prevIt = std::prev(it); 

     std::cout << "--" << std::endl;
     std::cout << "current ophit: " << opch << " " << ophit.StartTime();

     if( m_filter_intime ){
       
       // skip if too close in time with prev ophit 
       // check same channel + within time window
       if ( prevIt != sortableOpHits.begin() ){
         
         auto const &prevOphit = *prevIt;
         unsigned int prevopch = prevOphit.OpChannel();
         if( prevopch == opch && ophit.StartTime()-prevOphit.StartTime() < m_time_window ){
           continue;
         }
       }
     
       // skip if too close in time with next ophit 
       // check same channel + within time window
       if ( nextIt != sortableOpHits.end() ) {

         auto const &nextOphit = *nextIt;
         unsigned int nextopch = nextOphit.OpChannel();
         if( nextopch == opch && nextOphit.StartTime()-ophit.StartTime() < m_time_window ){
           continue;
         }
       }
  
     }
     
     m_channel = opch;
     m_time = ophit.StartTime();
     m_amplitude = ophit.Amplitude()*adc_to_mV;
     m_integral = ophit.Area()*adc_to_pC/echarge; 
   
     // fill histograms  
     hintegral[opch]->Fill( m_integral );
     hamplitude[opch]->Fill( m_amplitude );

     // fill ophit tree
     m_ophit_tree->Fill();   

   } 
   
} // end analyze


DEFINE_ART_MODULE(pmtcalo::PMTBackgroundphotonsCalibration)
