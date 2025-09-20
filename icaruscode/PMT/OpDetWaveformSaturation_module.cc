////////////////////////////////////////////////////////////////////////
// Class:       OpDetWaveformSaturation
// Plugin Type: analyzer (Unknown Unknown)
// File:        OpDetWaveformSaturation_module.cc
//
// Generated at Tue Aug 22 17:27:03 2023 by Matteo Vicenzi using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

// framework libraries
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art_root_io/TFileService.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// ROOT libraries
#include "TTree.h"
#include "TFile.h"

// C/C++ standard libraries
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
namespace icarus {   class OpDetWaveformSaturation; }

class icarus::OpDetWaveformSaturation : public art::EDAnalyzer {
public:
  /// configuration of the module.
  struct Config {

    fhicl::Atom<art::InputTag> InputLabel {
        fhicl::Name("InputLabel"), 
        fhicl::Comment("Raw waveform input label to be used")
    };

    fhicl::Atom<art::InputTag> TriggerLabel {
        fhicl::Name("TriggerLabel"),
        fhicl::Comment("Label for the Trigger fragment label")
    };

    fhicl::Atom<int> SaturationLevel {
	fhicl::Name("SaturationLevel"),
	fhicl::Comment("Absolute ADC level for waveform saturation")
    };
	    
  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  
  /// constructor: just reads the configuration
  explicit OpDetWaveformSaturation(Parameters const& config);

  /// process the event
  void analyze(art::Event const& e) override;

  short getMedian(std::vector<short> v);
  short getMin(std::vector<short> v);

  void beginJob();
  void endJob();

private:

  /// inputs
  art::InputTag const fInputLabel;
  art::InputTag const fTriggerLabel;
  int const fSaturationLevel;

  /// data members
  art::ServiceHandle<art::TFileService> tfs;
  TTree *fTree;
  int m_run;
  int m_event;
  int m_timestamp;
  uint64_t m_trigger_timestamp;
  int m_nsaturated; //channels with saturation in event
  std::vector<int> m_channels; //all channel ids
  std::vector<short> m_baselines; //all baselines
  std::vector<short> m_wf_mins;   //all wf minimum values

};

// --------------------------------------------------------------------------
icarus::OpDetWaveformSaturation::OpDetWaveformSaturation(Parameters const& config)
  : art::EDAnalyzer(config)
  , fInputLabel{ config().InputLabel() }
  , fTriggerLabel{ config().TriggerLabel() }
  , fSaturationLevel{ config().SaturationLevel() }
{
    // configuration checks
    if (fInputLabel.empty()) {
        throw art::Exception{ art::errors::Configuration }
          << "The input waveform data product ('"
          << config().InputLabel.name() << "') is empty.\n";
    }
    if (fTriggerLabel.empty()) {
        throw art::Exception{ art::errors::Configuration }
          << "The input trigger data product ('"
          << config().TriggerLabel.name() << "') is empty.\n";
    }
    
    // consumes
    consumes<std::vector<raw::OpDetWaveform>>(fInputLabel);
}

// ---------------------------------------------------------------------------
void icarus::OpDetWaveformSaturation::beginJob() {

  fTree = tfs->make<TTree>("sattree", "saturation info" );
  fTree->Branch("run",&m_run);
  fTree->Branch("event",&m_event);
  fTree->Branch("timestamp",&m_timestamp);
  fTree->Branch("trigger_timestamp",&m_trigger_timestamp);
  fTree->Branch("n_saturated",&m_nsaturated);
  fTree->Branch("channels",&m_channels);
  fTree->Branch("baselines",&m_baselines);
  fTree->Branch("min_wf",&m_wf_mins);

}

// ---------------------------------------------------------------------------
short icarus::OpDetWaveformSaturation::getMedian(std::vector<short> v)
{
	std::sort(v.begin(),v.end());
	if( v.size() % 2 == 0 ) //even
		return 0.5*(v.at(v.size()/2-1)+v.at(v.size()/2)); //average the two
	else
		return v.at(v.size()/2);
}

short icarus::OpDetWaveformSaturation::getMin(std::vector<short> v)
{
	std::sort(v.begin(),v.end());
	return v.at(0);
}

void icarus::OpDetWaveformSaturation::analyze(art::Event const& e)
{
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second    

  // get trigger timestamp
  art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
  e.getByLabel( fTriggerLabel, trigger_handle );
  m_trigger_timestamp = trigger_handle->triggerTimestamp;

  // get waveforms 
  auto const& opDetWfs = e.getProduct<std::vector<raw::OpDetWaveform>>(fInputLabel);
  m_nsaturated = 0;

  for(  auto const & wf : opDetWfs  ){

    m_channels.push_back(wf.ChannelNumber());
    m_baselines.push_back( getMedian(wf.Waveform()) );  
   
    short min = getMin(wf.Waveform());
    m_wf_mins.push_back( min );
	
    if (min < fSaturationLevel)
	m_nsaturated++;  
  }

  fTree->Fill();

  m_channels.clear();
  m_baselines.clear();
  m_wf_mins.clear();
}

// ----------------------------------------------------------------------------
void icarus::OpDetWaveformSaturation::endJob()
{}

// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(icarus::OpDetWaveformSaturation)
