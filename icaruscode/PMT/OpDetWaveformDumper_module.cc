////////////////////////////////////////////////////////////////////////
// Class:       OpDetWaveformDumper
// Plugin Type: analyzer (Unknown Unknown)
// File:        OpDetWaveformDumper_module.cc
//
// Generated at Fri May 19 15:59:03 2023 by Matteo Vicenzi using cetskelgen
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
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// ROOT libraries
#include "TTree.h"

// C/C++ standard libraries
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
namespace icarus {   class OpDetWaveformDumper; }
/**
 * @brief Dump PMT optical waveforms to a ROOT tree
 *
 * This module reads the raw optical waveforms and saves them in
 * a simple ROOT tree for quick access and visualization.
 * If a baseline data product is specified, it dumps it in the tree as well.
 * 
 * 
 * Input
 * ------
 * 
 * * `std::vector<raw::OpDetWaveform>` data products (as for `InputLabel`)
 * * `std::vector<icarus::WaveformBaseline>` (optional) data product (as for `BaselineLabel`)
 * 
 * Output
 * -------
 * 
 * * A ROOT tree, storing optical waveforms per event and channel. 
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * * `InputLabel` (input tag, mandatory): list of optical waveforms data
 *   products to dump. It must be non-empty.
 * * `BaselineLabel` (input tag, optional): baseline data product to be used. If not 
 *   specified or available, a simple median is used.
 *
 */

class icarus::OpDetWaveformDumper : public art::EDAnalyzer {
public:
  /// configuration of the module.
  struct Config {

    fhicl::Atom<art::InputTag> InputLabel {
        fhicl::Name("InputLabel"), 
        fhicl::Comment("Raw waveform input label to be used")
	    };
    
    fhicl::Atom<art::InputTag> BaselineLabel {
        fhicl::Name("BaselineLabel"), 
        fhicl::Comment("Waveform baseline label to be used"),
        ""
	    };
	    
  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  
  /// constructor: just reads the configuration
  explicit OpDetWaveformDumper(Parameters const& config);

  /// compute simple waveform median
  template<typename T> T Median(std::vector<T> data) const;

  /// process the event
  void analyze(art::Event const& e) override;
  void beginJob();

private:

  art::ServiceHandle<art::TFileService> tfs;
  
  /// inputs
  art::InputTag const fInputLabel;
  art::InputTag const fBaselineLabel;
  bool fBaselineOK = true;

  /// data members
  TTree *fTree;
  int m_run;
  int m_event;
  int m_timestamp;
  int m_channel_id;
  double m_tstart;
  double m_baseline;
  int m_nsize;
  std::vector<short> m_wf;

};

// --------------------------------------------------------------------------
icarus::OpDetWaveformDumper::OpDetWaveformDumper(Parameters const& config)
  : art::EDAnalyzer(config)
  , fInputLabel{ config().InputLabel() }
  , fBaselineLabel{ config().BaselineLabel() }
{
    // configuration checks
    if (fInputLabel.empty()) {
        throw art::Exception{ art::errors::Configuration }
          << "The input waveform data product ('"
          << config().InputLabel.name() << "') is empty.\n";
    }
    
    if (fBaselineLabel.empty()) {
        mf::LogInfo("OpDetWaveformDumper") << "No waveform baseline product selected, defaulting to simple median";
        fBaselineOK = false;
    }
    
    // consumes
    consumes<std::vector<raw::OpDetWaveform>>(fInputLabel);
}

// ---------------------------------------------------------------------------
void icarus::OpDetWaveformDumper::beginJob() {

  fTree = tfs->make<TTree>("wftree", "waveform info" );
  fTree->Branch("run",&m_run);
  fTree->Branch("event",&m_event);
  fTree->Branch("timestamp",&m_timestamp);
  fTree->Branch("channel_id",&m_channel_id);
  fTree->Branch("tstart",&m_tstart);
  fTree->Branch("baseline",&m_baseline);
  fTree->Branch("nsize",&m_nsize);
  fTree->Branch("wf",&m_wf);

}

// ---------------------------------------------------------------------------

void icarus::OpDetWaveformDumper::analyze(art::Event const& e)
{
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second    

  // get waveforms 
  auto const& opDetWfs = e.getProduct<std::vector<raw::OpDetWaveform>>(fInputLabel);

  // get baselines
  std::vector<icarus::WaveformBaseline> wfBaselines;
  if( fBaselineOK ) wfBaselines = e.getProduct<std::vector<icarus::WaveformBaseline>>(fBaselineLabel);

  for( size_t i=0; i<opDetWfs.size(); i++ ){

    auto &wf = opDetWfs[i];
    m_channel_id = wf.ChannelNumber();
    m_tstart = wf.TimeStamp();
    m_wf = wf.Waveform();
    m_nsize = m_wf.size();	

    if( fBaselineOK ) 
      m_baseline = wfBaselines[i].baseline();
    else
      m_baseline = Median( wf.Waveform() );

    fTree->Fill();
  }
}

// --------------------------------------------------------------------------
template<typename T> T icarus::OpDetWaveformDumper::Median(std::vector<T> data) const
{
  std::nth_element( data.begin(), data.begin() + data.size()/2, data.end() );
  return data[ data.size()/2 ];
}

// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(icarus::OpDetWaveformDumper)
