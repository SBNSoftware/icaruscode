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

// ROOT libraries
#include "TTree.h"
#include "TFile.h"

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
 * 
 * 
 * Input
 * ------
 * 
 * * `std::vector<raw::OpDetWaveform>` data products (as for `InputLabels`)
 * 
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
 * * `InputLabels` (list of input tags, mandatory): the list of optical waveforms data
 *   products to apply the filter on. It must be non-empty.
 * * `OutputFile` (string, default: `wf_output.root`): name of output ROOT file
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
	    
    fhicl::Atom<std::string> OutputFile {
      fhicl::Name("OutputFile"),
      fhicl::Comment("output ROOT tree filename"),
      "wf_output.root" // default
    };
    
  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  
  /// constructor: just reads the configuration
  explicit OpDetWaveformDumper(Parameters const& config);

  /// process the event
  void analyze(art::Event const& e) override;
  void beginJob();
  void endJob();

private:

  /// inputs
  art::InputTag const fInputLabel;
  std::string const fOutputFile;

  /// data members
  TFile *fFile;
  TTree *fTree;
  int m_run;
  int m_event;
  int m_timestamp;
  int m_channel_id;
  double m_tstart;
  int m_nsize;
  std::vector<short> m_wf;

};

// --------------------------------------------------------------------------
icarus::OpDetWaveformDumper::OpDetWaveformDumper(Parameters const& config)
  : art::EDAnalyzer(config)
  , fInputLabel{ config().InputLabel() }
  , fOutputFile{ config().OutputFile() }
{
    // configuration checks
    if (fInputLabel.empty()) {
        throw art::Exception{ art::errors::Configuration }
          << "The input waveform data product ('"
          << config().InputLabel.name() << "') is empty.\n";
    }
    
    // consumes
    consumes<std::vector<raw::OpDetWaveform>>(fInputLabel);
}

// ---------------------------------------------------------------------------
void icarus::OpDetWaveformDumper::beginJob() {

  fFile = new TFile(fOutputFile.c_str(),"RECREATE");
    
  fTree = new TTree("wftree", "waveform info" );
  fTree->Branch("run",&m_run);
  fTree->Branch("event",&m_event);
  fTree->Branch("timestamp",&m_timestamp);
  fTree->Branch("channel_id",&m_channel_id);
  fTree->Branch("tstart",&m_tstart);
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

  for(  auto const & wf : opDetWfs  ){

    m_channel_id = wf.ChannelNumber();
    m_tstart = wf.TimeStamp();
    m_wf = wf.Waveform();
    m_nsize = m_wf.size();	

    fTree->Fill();
  }
}

// ----------------------------------------------------------------------------
void icarus::OpDetWaveformDumper::endJob()
{
  fFile->cd();
  fTree->Write();
  fFile->Close();
}

// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(icarus::OpDetWaveformDumper)
