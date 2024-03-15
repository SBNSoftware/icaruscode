////////////////////////////////////////////////////////////////////////
// Class:       PMTBeamSignalsExtractor
// Plugin Type: producer (Unknown Unknown)
// File:        PMTBeamSignalsExtractor_module.cc
//
// Generated at Sun Feb 11 11:37:14 2024 by Matteo Vicenzi using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "icaruscode/IcarusObj/PMTWaveformTimeCorrection.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataobj/RawData/OpDetWaveform.h"

#include "TTree.h"

#include <memory>
#include <vector>
#include <string>
#include <algorithm>

namespace icarus {
  namespace timing {
    class PMTBeamSignalsExtractor;
  }
}


class icarus::timing::PMTBeamSignalsExtractor : public art::EDProducer {
public:

  explicit PMTBeamSignalsExtractor(fhicl::ParameterSet const& pset);
		
  //process waveforms
  template<typename T> T Median(std::vector<T> data) const;
  template<typename T> static size_t getMaxBin(std::vector<T> const& vv, size_t startElement, size_t endElement);
  template<typename T> static size_t getMinBin(std::vector<T> const& vv, size_t startElement, size_t endElement);
  template<typename T> static size_t getStartSample( std::vector<T> const& vv, T thres );

  //trigger-hardware timing correction
  double getTriggerCorrection(int channel, std::vector<icarus::timing::PMTWaveformTimeCorrection> const& corrections);

  // Plugins should not be copied or assigned.
  PMTBeamSignalsExtractor(PMTBeamSignalsExtractor const&) = delete;
  PMTBeamSignalsExtractor(PMTBeamSignalsExtractor&&) = delete;
  PMTBeamSignalsExtractor& operator=(PMTBeamSignalsExtractor const&) = delete;
  PMTBeamSignalsExtractor& operator=(PMTBeamSignalsExtractor&&) = delete;

  void beginJob();
  void produce(art::Event& e) override;

private:

  art::ServiceHandle<art::TFileService> tfs;
  bool fDebugTrees;
  bool fSaveWaveforms;
  art::InputTag fTriggerCorrectionLabel;
  short int fADCThreshold;
  std::vector<fhicl::ParameterSet> fBoardSetup;

  std::map<int, std::string> fBoardBySpecialChannel;

  TTree* fRWMTree;
  TTree* fEWTree;

  int m_run;
  int m_event;
  int m_timestamp; 
  double m_trigger_time;
		
  // rwm signal
  int m_n_rwm;
  int m_rwm_channel;
  double m_rwm_wfstart;
  double m_rwm_utime;
  double m_rwm_time;
  std::vector<short> m_rwm_wf;
 
  // ew signal 
  int m_n_ew;
  int m_ew_channel;
  double m_ew_wfstart;
  double m_ew_utime;
  double m_ew_time;
  std::vector<short> m_ew_wf;

  // these channels have been chosen because they 
  // were not affected by the mapping change 
  // (they stayed on the same board)
  // FIXME: get it from the mapping db
  
  std::map< std::string, int> fSingleChannelPerBoard = 
  {
    { "EE-BOT-C", 4 },
    { "EE-BOT-B", 24 },
    { "EE-TOP-C", 54 },
    { "EE-TOP-B", 64 },
    { "EW-BOT-C", 94 },
    { "EW-BOT-B", 114 },
    { "EW-TOP-C", 144 },
    { "EW-TOP-B", 154 },
    { "WE-BOT-C", 184 },
    { "WE-BOT-B", 204 },
    { "WE-TOP-C", 234 },
    { "WE-TOP-B", 244 },
    { "WW-BOT-C", 274 },
    { "WW-BOT-B", 294 },
    { "WW-TOP-C", 320 },
    { "WW-TOP-B", 339 },
  };

};


icarus::timing::PMTBeamSignalsExtractor::PMTBeamSignalsExtractor(fhicl::ParameterSet const& pset)
  : EDProducer{pset}
  , fDebugTrees( pset.get<bool>("DebugTrees") )
  , fSaveWaveforms( pset.get<bool>("SaveWaveforms") )
  , fTriggerCorrectionLabel( pset.get<art::InputTag>("TriggerCorrectionLabel") )
  , fADCThreshold( pset.get<short int>("ADCThreshold") )
  , fBoardSetup( pset.get<std::vector<fhicl::ParameterSet>>("BoardSetup") )
{
  
  for (fhicl::ParameterSet const& setup : fBoardSetup ) {
    auto innerSet = setup.get<std::vector<fhicl::ParameterSet>>("SpecialChannels");
    fBoardBySpecialChannel[ innerSet[0].get<int>("Channel") ] = setup.get<std::string>("Name");
  }

  // Call appropriate produces<>() functions here.

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::beginJob()
{
  if( !fDebugTrees ) return;

  fRWMTree = tfs->make<TTree>("rwmtree","RWM info");
  fRWMTree->Branch("run",&m_run);
  fRWMTree->Branch("event",&m_event);
  fRWMTree->Branch("timestamp",&m_timestamp);
  fRWMTree->Branch("trigger_time",&m_trigger_time);
  fRWMTree->Branch("n_rwm",&m_n_rwm);
  fRWMTree->Branch("rwm_channel",&m_rwm_channel);
  fRWMTree->Branch("rwm_wfstart", &m_rwm_wfstart);
  fRWMTree->Branch("rwm_utime", &m_rwm_utime);
  fRWMTree->Branch("rwm_time", &m_rwm_time);
  
  fEWTree = tfs->make<TTree>("ewtree","EW info");
  fEWTree->Branch("run",&m_run);
  fEWTree->Branch("event",&m_event);
  fEWTree->Branch("timestamp",&m_timestamp);
  fEWTree->Branch("trigger_time",&m_trigger_time);
  fEWTree->Branch("n_ew",&m_n_ew);
  fEWTree->Branch("ew_channel",&m_ew_channel);
  fEWTree->Branch("ew_wfstart", &m_ew_wfstart);
  fEWTree->Branch("ew_utime", &m_ew_utime);
  fEWTree->Branch("ew_time", &m_ew_time);

  if(fSaveWaveforms){
    fRWMTree->Branch("rwm_wf",&m_rwm_wf);
    fEWTree->Branch("ew_wf",&m_ew_wf);
  }

}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::produce(art::Event& e)
{
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second    
  
  detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));
  m_trigger_time = detTimings.TriggerTime().value();

  m_rwm_wf.clear();
  m_ew_wf.clear();

  auto const& wfCorrections = e.getProduct<std::vector<icarus::timing::PMTWaveformTimeCorrection>>(fTriggerCorrectionLabel);
  int ntrig = wfCorrections.size();
			
  if( ntrig < 1 ) 
    mf::LogError("PMTBeamSignalsExtractor") << "Not found PMTWaveformTimeCorrections with label '" 
                                            << fTriggerCorrectionLabel.label() << "'"; 
  else if ( ntrig < 360 )
    mf::LogError("PMTBeamSignalsExtractor") << "Missing some PMTWaveformTimeCorrections with label '" 
                                            << fTriggerCorrectionLabel.label() << "'"; 
  
  auto const& ewWaveforms = e.getProduct<std::vector<raw::OpDetWaveform>>("daqPMT:EW");
  auto const& rwmWaveforms = e.getProduct<std::vector<raw::OpDetWaveform>>("daqPMT:RWM");
  m_n_ew = ewWaveforms.size();
  m_n_rwm = rwmWaveforms.size();
  if( m_n_ew< 1 ) 
    mf::LogError("PMTBeamSignalsExtractor") << "Not found raw::OpDetWaveform with label 'daqPMT:EW'"; 
  if( m_n_rwm< 1 ) 
    mf::LogError("PMTBeamSignalsExtractor") << "Not found raw::OpDetWaveform with label 'daqPMT:RWM'"; 

  for( auto const & wave : ewWaveforms ){
    
    detinfo::timescales::electronics_time tstart = util::quantities::points::microsecond{wave.TimeStamp()};
    double sample_rise = getStartSample( wave.Waveform(), fADCThreshold );
    
    m_ew_channel = wave.ChannelNumber();
    m_ew_wfstart = tstart.value() ;
    m_ew_utime = (sample_rise > 0) ? tstart.value() + 0.002*sample_rise : -1;
    m_ew_time = (sample_rise > 0) ? m_ew_utime + getTriggerCorrection(m_ew_channel, wfCorrections) : -1;
    if(fSaveWaveforms) m_ew_wf = wave.Waveform();
  
    fEWTree->Fill();
  }
  
  for( auto const & wave : rwmWaveforms ){
    
    detinfo::timescales::electronics_time tstart = util::quantities::points::microsecond{wave.TimeStamp()};
    double sample_rise = getStartSample( wave.Waveform(), fADCThreshold );
    
    m_rwm_channel = wave.ChannelNumber();
    m_rwm_wfstart = tstart.value() ;
    m_rwm_utime = (sample_rise > 0) ? tstart.value() + 0.002*sample_rise : -1;
    m_rwm_time = (sample_rise > 0) ? m_rwm_utime + getTriggerCorrection(m_rwm_channel, wfCorrections) : -1;
    if(fSaveWaveforms) m_rwm_wf = wave.Waveform();
  
    fRWMTree->Fill();
  }

}

// ---------------------------------------------------------------------------

template<typename T> T icarus::timing::PMTBeamSignalsExtractor::Median( std::vector<T> data ) const {

    std::nth_element( data.begin(), data.begin() + data.size()/2, data.end() );
    return data[ data.size()/2 ];

}

// -----------------------------------------------------------------------------

template<typename T>
  size_t icarus::timing::PMTBeamSignalsExtractor::getMinBin( 
        std::vector<T> const& vv, size_t startElement, size_t endElement ){

    auto minel = 
        std::min_element( vv.begin()+startElement, vv.begin()+endElement );
    size_t minsample = std::distance( vv.begin()+startElement, minel );

    return minsample;
}

// -----------------------------------------------------------------------------

template<typename T>
  size_t icarus::timing::PMTBeamSignalsExtractor::getMaxBin( 
            std::vector<T> const& vv, size_t startElement, size_t endElement){

    auto maxel = 
        std::max_element( vv.begin()+startElement, vv.begin()+endElement );
    
    size_t maxsample = std::distance( vv.begin()+startElement, maxel );

    return maxsample;
} 

// -----------------------------------------------------------------------------

template<typename T>
  size_t icarus::timing::PMTBeamSignalsExtractor::getStartSample( std::vector<T> const& vv, T thres ){
    
    // We are thinking in inverted polarity
    size_t minbin = getMinBin( vv, 0, vv.size() );

    //Search only a cropped region of the waveform backward from the min
    size_t maxbin =  minbin-20;

    // Now we crawl betweem maxbin and minbin and we stop when:
    // bin value  > (maxbin value - bin value )*0.2
    size_t startbin = maxbin;
    auto delta = vv[maxbin]-vv[minbin];

    if( delta < thres ) //just noise
      return 0; //return first bin   

    for( size_t bin=maxbin; bin<minbin; bin++ ){
      auto val = vv[maxbin]-vv[bin];
      if( val >= 0.2*delta ){ //20%
        startbin = bin - 1;
        break;
      }
    }

    if( startbin < maxbin ){
      startbin=maxbin;
    }

    return startbin;
}

// -----------------------------------------------------------------------------
  
double icarus::timing::PMTBeamSignalsExtractor::getTriggerCorrection(int channel, 
       std::vector<icarus::timing::PMTWaveformTimeCorrection> const& corrections){

  // get the board name, convert to digitizer_label
  std::string board = fBoardBySpecialChannel[channel];
  
  std::string head = "icaruspmt";
  std::string dash = "-";
  std::string letter = (board.substr(board.size()-2, board.size())=="02") ? "B" : "C";
  
  board.erase( board.find(head), head.size() );
  std::transform(board.begin(), board.end(), board.begin(), ::toupper);
  board.insert( 2, dash);
  board.insert( 6, dash);
  std::string digitizer_label = board.substr(0,board.size()-2) + letter;

  // trigger-hardware corrections are shared by all channels on the same board
  // mapping currently does not expose channel<->board relationship
  // using ad-hoc configuration... FIXME!
  int pmtch = fSingleChannelPerBoard[digitizer_label];

  // trigger-hardware correction are in order
  // index of vector is pmtch
  return corrections.at(pmtch).startTime;

}

DEFINE_ART_MODULE(icarus::timing::PMTBeamSignalsExtractor)
