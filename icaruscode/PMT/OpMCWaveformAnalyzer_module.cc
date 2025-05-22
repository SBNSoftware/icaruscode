////////////////////////////////////////////////////////////////////////
// Class:       OpMCWaveformAnalyzer
// Plugin Type: analyzer (art v3_05_00)
// File:        OpMCWaveformAnalyzer_module.cc
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
#include "lardataobj/Simulation/SimPhotons.h"

#include "larana/OpticalDetector/IPedAlgoMakerTool.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"

#include "TTree.h"
#include <vector>
#include <algorithm>

namespace icarus {
  class OpMCWaveformAnalyzer;
}


class icarus::OpMCWaveformAnalyzer : public art::EDAnalyzer {


public:

  explicit OpMCWaveformAnalyzer(fhicl::ParameterSet const& pset);

  OpMCWaveformAnalyzer(OpMCWaveformAnalyzer const&) = delete;
  OpMCWaveformAnalyzer(OpMCWaveformAnalyzer&&) = delete;
  OpMCWaveformAnalyzer& operator=(OpMCWaveformAnalyzer const&) = delete;
  OpMCWaveformAnalyzer& operator=(OpMCWaveformAnalyzer&&) = delete;

  virtual void beginJob() override;
  void analyze(art::Event const& event) override;

private:

  // input tags: ophits, wfs, simphotons
  art::InputTag m_ophit_label;
  art::InputTag m_simphotons_label;
  art::InputTag m_wf_label;
  // input offset
  double m_simphoffset_ns;
  bool m_onlysimphotons;	

  // output ttree
  TTree *m_wf_tree; // extracted waveforms
  int m_run;
  int m_event;
  uint32_t m_timestamp;
  unsigned int m_channel;
  unsigned int m_nsize;
  double m_wfstart;
  std::vector<double> m_baselines;
  std::vector<short> m_wf;

  // optical hits
  unsigned int m_n_ophits;
  std::vector<double> m_start_time;
  std::vector<double> m_rise_time;
  std::vector<double> m_peak_time;
  std::vector<double> m_time_width;
  std::vector<double> m_amplitude;
  std::vector<double> m_integral;
  std::vector<std::size_t> m_hitstart;
  std::vector<std::size_t> m_hitpeak;
  std::vector<std::size_t> m_hitend;

  // simphotons
  unsigned int m_n_simphotons;
  std::vector<double> m_sim_time;
  std::vector<double> m_sim_energy;
  std::vector<std::size_t> m_sim_sample;
  std::vector<double> m_sim_start_x;
  std::vector<double> m_sim_start_y;
  std::vector<double> m_sim_start_z;

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


icarus::OpMCWaveformAnalyzer::OpMCWaveformAnalyzer(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{

   m_ophit_label = pset.get<art::InputTag>("OpHitModule", "ophit");
   m_simphotons_label = pset.get<art::InputTag>("SimPhotonsModule");
   m_wf_label = pset.get<art::InputTag>("WaveformModule", "opdaq");
   m_simphoffset_ns = pset.get<double>("SimPhotonsOffset", 55.1); //ns
   m_onlysimphotons = pset.get<bool>("SaveOnlyWithSimPhotons", false);

   fPedAlgo = art::make_tool<opdet::IPedAlgoMakerTool>(pset.get<fhicl::ParameterSet>("PedAlgoRollingMeanMaker"))->makeAlgo();
}


//------------------------------------------------------------------------------


void icarus::OpMCWaveformAnalyzer::beginJob()
{

  // create full tree
  m_wf_tree = tfs->make<TTree>("mcwfs","mc waveforms tree");

  m_wf_tree->Branch("run", &m_run, "run/I" );
  m_wf_tree->Branch("event", &m_event, "event/I" );
  m_wf_tree->Branch("timestamp", &m_timestamp, "timestamp/I" );
  m_wf_tree->Branch("channel", &m_channel, "channel/I" );
  m_wf_tree->Branch("nsize", &m_nsize);
  m_wf_tree->Branch("wfstart", &m_wfstart, "wfstart/D" );
  m_wf_tree->Branch("baselines", &m_baselines);
  m_wf_tree->Branch("wf", &m_wf);

  m_wf_tree->Branch("n_ophits", &m_n_ophits);
  m_wf_tree->Branch("start_time", &m_start_time );
  m_wf_tree->Branch("peak_time", &m_peak_time );
  m_wf_tree->Branch("rise_time", &m_rise_time);
  m_wf_tree->Branch("time_width", &m_time_width );
  m_wf_tree->Branch("amplitude", &m_amplitude);
  m_wf_tree->Branch("integral", &m_integral );
  m_wf_tree->Branch("hitpeak", &m_hitpeak);
  m_wf_tree->Branch("hitstart", &m_hitstart);
  m_wf_tree->Branch("hitend", &m_hitend);

  m_wf_tree->Branch("n_simphotons", &m_n_simphotons);
  m_wf_tree->Branch("sim_time", &m_sim_time );
  m_wf_tree->Branch("sim_energy", &m_sim_energy );
  m_wf_tree->Branch("sim_sample", &m_sim_sample );
  m_wf_tree->Branch("sim_start_x", &m_sim_start_x );
  m_wf_tree->Branch("sim_start_y", &m_sim_start_y );
  m_wf_tree->Branch("sim_start_z", &m_sim_start_z );

}

//-----------------------------------------------------------------------------

void icarus::OpMCWaveformAnalyzer::analyze(art::Event const& event)
{

  // event info: this is mainly to save the timestamps
  // the timestamp of the first event is used to sort files in the calibration db 
  detinfo::DetectorTimings const& timings = detinfo::makeDetectorTimings(
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob());

  m_run = event.id().run();
  m_event = event.id().event();
  m_timestamp = event.time().timeHigh(); // We just need precision at the s level

  // get the data products: ophits, waveforms, and simphotons
  auto const & ophitCollection = event.getProduct<std::vector<recob::OpHit>>(m_ophit_label);
  auto const & wfCollection = event.getProduct<std::vector<raw::OpDetWaveform>>(m_wf_label);

  static const std::vector<sim::SimPhotons> emptySimPhotons;
  auto const & simphotonsCollection = (!m_simphotons_label.empty())
                                      ? event.getProduct<std::vector<sim::SimPhotons>>(m_simphotons_label)
                                      : emptySimPhotons;

  // make sure they are sorted
  auto sortableOpHits = ophitCollection;
  std::sort( sortableOpHits.begin(), sortableOpHits.end(), sortByChannelTime );
 
  // LOOP OVER ALL WAVEFORMS 
  for (auto jt = wfCollection.begin(); jt != wfCollection.end(); ++jt) {

    //clean up vectors at each loop
    m_baselines.clear();
    m_wf.clear();
    m_start_time.clear();
    m_rise_time.clear();
    m_peak_time.clear();
    m_time_width.clear();
    m_amplitude.clear();
    m_integral.clear();
    m_hitstart.clear();
    m_hitpeak.clear();
    m_hitend.clear();
    m_sim_time.clear();
    m_sim_energy.clear();
    m_sim_sample.clear();
    m_sim_start_x.clear();
    m_sim_start_y.clear();
    m_sim_start_z.clear();

    // pick this particular waveform
    m_nsize = jt->Waveform().size();
    m_channel = jt->ChannelNumber();
    m_wfstart = jt->TimeStamp() - fTriggerTime;
    m_wf = jt->Waveform();

    fPedAlgo->Evaluate(jt->Waveform());
    m_baselines = fPedAlgo->Mean();

    double wfend = m_wfstart + m_nsize*fOpticalTick;

    //std::cout << "Opch " << m_channel << " nsize " << m_nsize << " start " << m_wfstart << " end " << wfend << std::endl;

    // LOOP OVER ALL OPHITS: hits are saved only if:
    // 1) coming from the same channel
    // 2)  start time is within wf boundaries
    for (auto it = sortableOpHits.begin(); it != sortableOpHits.end(); ++it) {
    
      auto const &ophit = *it;  // current element
      unsigned int opch = ophit.OpChannel();
      
      if( opch != m_channel ) continue; //skip if different channel
      if( ophit.StartTime() < m_wfstart || ophit.StartTime() > wfend ) continue; //skip if out of current wf

      // find hit peak sample 
      double diff_peak = ophit.PeakTime() - m_wfstart;
      m_peak_time.push_back(ophit.PeakTime());
      m_hitpeak.push_back(static_cast<std::size_t>(std::round(diff_peak/fOpticalTick)));
    
      // find hit start sample
      double diff_start = ophit.StartTime() - m_wfstart;
      m_start_time.push_back(ophit.StartTime());
      m_hitstart.push_back(static_cast<std::size_t>(std::round(diff_start/fOpticalTick)));
          
      //find hit end sample
      double diff_end = ophit.StartTime() + ophit.Width() - m_wfstart;
      m_hitend.push_back(static_cast<std::size_t>(std::round(diff_end/fOpticalTick)));
      
      m_time_width.push_back(ophit.Width());
      m_rise_time.push_back(ophit.StartTime()+ophit.RiseTime());
      m_amplitude.push_back(ophit.Amplitude()); // ADC 
      m_integral.push_back(ophit.Area()); // ADC x tick = ADC x 2ns

    }
          
    m_n_ophits = m_peak_time.size(); // count how many ophits
    //std::cout << "nophits " << m_n_ophits << std::endl;
    
    // LOOP OVER ALL SIMPHOTONS: photons are saved only if:
    // 1) coming from the same channel
    // 2) start time is within wf boundaries
    for (auto it = simphotonsCollection.begin(); it != simphotonsCollection.end(); ++it) {
    
      auto const &simphotons = *it;  // current element
      unsigned int opch = simphotons.OpChannel();
      
      if( opch != m_channel ) continue; //skip if different channel

      // looping through each photon in this channel
      for(auto const& ph : simphotons){
        
        detinfo::timescales::simulation_time const photonTime { ph.Time };
        detinfo::timescales::trigger_time const mytime = timings.toTriggerTime(photonTime);
        double t = mytime.value() + m_simphoffset_ns*0.001; // in us

        if( t < m_wfstart || t > wfend ) continue; //skip if out of current wf

        m_sim_time.push_back(t);
        m_sim_energy.push_back(ph.Energy);
	m_sim_start_x.push_back(ph.InitialPosition.X());
	m_sim_start_y.push_back(ph.InitialPosition.Y());
	m_sim_start_z.push_back(ph.InitialPosition.Z());

        double diff = t - m_wfstart;
        m_sim_sample.push_back(static_cast<std::size_t>(std::round(diff/fOpticalTick)));

      }
    }

    m_n_simphotons = m_sim_sample.size(); // count how many ophits
    //std::cout << "nsimphotons " << m_n_simphotons << std::endl;

    if( m_n_ophits == 0 && m_n_simphotons == 0){
      std::cout << "No Ophits or SimPhotons!" << std::endl;
      continue; //don't save!
    }

    if(m_onlysimphotons && m_n_simphotons == 0) continue;

    // fill tree
    m_wf_tree->Fill(); 

   } 
   
} // end analyze


DEFINE_ART_MODULE(icarus::OpMCWaveformAnalyzer)
