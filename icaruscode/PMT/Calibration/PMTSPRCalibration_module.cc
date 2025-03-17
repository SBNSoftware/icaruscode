////////////////////////////////////////////////////////////////////////
// Class:       PMTSPRCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTSPRCalibration_module.cc
//
// Generated at Mon Mar 10 16:10:37 2025 by Matteo Vicenzi
// 
// mailto: mvicenzi@bnl.gov
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
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/PMT/Calibration/CaloTools/Utils.h"
#include "larana/OpticalDetector/IPedAlgoMakerTool.h"

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

#include <vector>
#include <algorithm>

namespace pmtcalo {
  class PMTSPRCalibration;
}


class pmtcalo::PMTSPRCalibration : public art::EDAnalyzer {


public:

  explicit PMTSPRCalibration(fhicl::ParameterSet const& pset);

  PMTSPRCalibration(PMTSPRCalibration const&) = delete;
  PMTSPRCalibration(PMTSPRCalibration&&) = delete;
  PMTSPRCalibration& operator=(PMTSPRCalibration const&) = delete;
  PMTSPRCalibration& operator=(PMTSPRCalibration&&) = delete;

  virtual void beginJob() override;
  bool inGate(double const & time, std::vector<double> & bounds) const;
  double getMedian(std::vector<short> v) const;

  void analyze(art::Event const& event) override;

private:

  // input tags: ophits + trigger
  art::InputTag m_ophit_label;
  art::InputTag m_trigger_label;
  art::InputTag m_wf_label;

  // output histograms
  TTree *m_wf_tree; // extracted waveforms
  int m_run;
  int m_event;
  uint32_t m_timestamp;
  unsigned int m_channel;
  double m_start_time;
  double m_rise_time;
  double m_peak_time;
  double m_time_width;
  double m_amplitude;
  double m_integral;
  std::size_t m_sample;
  unsigned int m_nsize;
  double m_wfstart;
  std::size_t m_hitpeak;
  std::size_t m_hitstart;
  std::size_t m_hitend;
  double m_median;
  std::vector<double> m_baselines;
  std::vector<short> m_wf;

  // dead/broken channels to be skipped
  std::vector<unsigned int> m_channel_mask;
  // limit to select in-gate ophits
  std::vector<double> m_filter_ingate;

  // select isolated ophits  
  bool m_filter_intime;
  double m_time_window;

  // waveform length
  std::size_t m_prepulseSamples;
  std::size_t m_afterpulseSamples;

  std::unique_ptr<pmtana::PMTPedestalBase> fPedAlgo;

  art::ServiceHandle<art::TFileService> tfs;

  // constants
  const double fTriggerTime = 1500.; // us
  const double fOpticalTick = 0.002; // 2ns

  static bool sortByChannelTime(const recob::OpHit &a, const recob::OpHit &b){
    if (a.OpChannel() == b.OpChannel()) return a.StartTime() < b.StartTime(); 
    return a.OpChannel() < b.OpChannel();
  };

};


//------------------------------------------------------------------------------


pmtcalo::PMTSPRCalibration::PMTSPRCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{

   m_ophit_label = pset.get<art::InputTag>("OpHitModule", "ophit");
   m_trigger_label = pset.get<art::InputTag>("TriggerModule");
   m_wf_label = pset.get<art::InputTag>("WaveformModule", "daqPMT");

   m_channel_mask = pset.get< std::vector< unsigned int > >
                            ("ChannelMask", std::vector< unsigned int >());
   m_filter_ingate = pset.get< std::vector<double> >("FilterInGate", std::vector< double >()); // in us
   m_filter_intime = pset.get<bool>("FilterInTime", true); 
   m_time_window = pset.get<double>("TimeWindow", 0.2 ); //in us
  
   m_prepulseSamples = pset.get<std::size_t>("prePulseSamples",100);
   m_afterpulseSamples = pset.get<std::size_t>("afterPulseSamples",400);

   fPedAlgo = art::make_tool<opdet::IPedAlgoMakerTool>(pset.get<fhicl::ParameterSet>("PedAlgoRollingMeanMaker"))->makeAlgo();
}


//------------------------------------------------------------------------------


void pmtcalo::PMTSPRCalibration::beginJob()
{

  // create full ophits tree for possible amplitude-cut selection
  // include event/timestamp data for later use
  m_wf_tree = tfs->make<TTree>("wfs","waveforms tree");

  m_wf_tree->Branch("run", &m_run, "run/I" );
  m_wf_tree->Branch("event", &m_event, "event/I" );
  m_wf_tree->Branch("timestamp", &m_timestamp, "timestamp/I" );
  m_wf_tree->Branch("channel", &m_channel, "channel/I" );
  m_wf_tree->Branch("start_time", &m_start_time, "start_time/D" );
  m_wf_tree->Branch("peak_time", &m_peak_time, "peak_time/D" );
  m_wf_tree->Branch("rise_time", &m_rise_time, "rise_time/D" );
  m_wf_tree->Branch("time_width", &m_time_width, "time_width/D" );
  m_wf_tree->Branch("amplitude", &m_amplitude, "amplitude/D" );
  m_wf_tree->Branch("integral", &m_integral, "integral/D" );
  m_wf_tree->Branch("sample", &m_sample);
  m_wf_tree->Branch("nsize", &m_nsize);
  m_wf_tree->Branch("wfstart", &m_wfstart, "wfstart/D" );
  m_wf_tree->Branch("hitpeak", &m_hitpeak);
  m_wf_tree->Branch("hitstart", &m_hitstart);
  m_wf_tree->Branch("hitend", &m_hitend);
  m_wf_tree->Branch("median", &m_median, "median/D");
  m_wf_tree->Branch("baselines", &m_baselines);
  m_wf_tree->Branch("wf", &m_wf);

}

//-----------------------------------------------------------------------------

bool pmtcalo::PMTSPRCalibration::inGate(double const & time, std::vector<double> & bounds ) const
{
  return (time >= bounds[0]) && ( time <= bounds[1] );
}

double pmtcalo::PMTSPRCalibration::getMedian(std::vector<short> v) const
{
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  double vn = v[v.size() / 2];
  
  // if odd
  if( v.size() % 2 != 0) 
    return vn;

  // if even
  std::nth_element(v.begin(), v.begin() + v.size()/2 - 1, v.end());
  double vnn = v[v.size() / 2 - 1];
  return (vn+vnn)/2;
}

//-----------------------------------------------------------------------------

void pmtcalo::PMTSPRCalibration::analyze(art::Event const& event)
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
     mf::LogError("pmtcalo::PMTSPRCalibration") << "No raw::Trigger associated to label: " << m_trigger_label.label();
   }

   // first thing we sort the ophit in their respective channels
   // creating the amplitude and integral histograms (w/o cuts)

   auto const & ophitCollection = event.getProduct<std::vector<recob::OpHit>>(m_ophit_label);
   auto const & wfCollection = event.getProduct<std::vector<raw::OpDetWaveform>>(m_wf_label);

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

     if( m_filter_intime ){
       
       // skip if too close in time with prev ophit 
       // check same channel + within time window
       if ( it != sortableOpHits.begin() ){
         
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
     m_start_time = ophit.StartTime();
     m_peak_time = ophit.PeakTime();
     m_rise_time = ophit.StartTime()+ophit.RiseTime();
     m_time_width = ophit.Width();
     m_amplitude = ophit.Amplitude(); // ADC 
     m_integral = ophit.Area(); // ADC x tick = ADC x 2ns


     for (auto jt = wfCollection.begin(); jt != wfCollection.end(); ++jt) {

      std::size_t nsize = jt->Waveform().size();

      // waveform must be from the same channel
      if( jt->ChannelNumber() != m_channel ) continue; 

      m_wfstart = jt->TimeStamp() - fTriggerTime;
      double wfend = m_wfstart + nsize*fOpticalTick;

      // ophit must be within the time window covered by the waveform
      if ( m_start_time < m_wfstart || m_start_time > wfend ) continue;
      
      // find hit peak sample 
      double diff_peak = m_peak_time - m_wfstart;
      m_sample = static_cast<std::size_t>(std::round(diff_peak/fOpticalTick));

      // find hit start sample
      double diff_start = m_start_time - m_wfstart;
      std::size_t start_sample = static_cast<std::size_t>(std::round(diff_start/fOpticalTick));
      std::size_t start_to_peak = m_sample - start_sample;

      //find hit end sample
      double diff_end = m_start_time + m_time_width - m_wfstart;
      std::size_t end_sample = static_cast<std::size_t>(std::round(diff_end/fOpticalTick));
      std::size_t peak_to_end = end_sample - m_sample;

      // find baseline
      fPedAlgo->Evaluate(jt->Waveform());
      std::vector<double> baselines = fPedAlgo->Mean();

      // select samples around the peak for saving
      std::size_t tick_start = (m_sample - m_prepulseSamples > 0) ? m_sample -  m_prepulseSamples : 0;
      std::size_t tick_end = (m_sample + m_afterpulseSamples < nsize) ? m_sample + m_afterpulseSamples : nsize;
      
      // size of waveform snippet
      m_nsize = tick_end - tick_start;
      // position of hit start/peak/end sample in waveform snippet
      m_hitpeak = m_prepulseSamples;
      m_hitstart = m_hitpeak - start_to_peak;
      m_hitend = m_hitpeak + peak_to_end;

      // cut the waveform
      m_wf = std::vector<short>(jt->Waveform().begin() + tick_start, jt->Waveform().begin() + tick_end);
      m_baselines = std::vector<double>(baselines.begin() + tick_start, baselines.begin() + tick_end);

      // get median from pre-pulse sample (up to hit start)
      std::vector<short> prewf = std::vector<short>(jt->Waveform().begin() + tick_start, jt->Waveform().begin() + start_sample);
      m_median = getMedian(prewf);

     }

    // fill tree
    m_wf_tree->Fill(); 

   } 
   
} // end analyze


DEFINE_ART_MODULE(pmtcalo::PMTSPRCalibration)
