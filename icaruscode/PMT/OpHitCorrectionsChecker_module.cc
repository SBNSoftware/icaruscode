////////////////////////////////////////////////////////////////////////
// Class:       OpHitCorrectionsChecker
// Plugin Type: analyzer (Unknown Unknown)
// File:        OpHitCorrectionsChecker_module.cc
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
#include "lardataobj/RecoBase/OpHit.h"
#include "icaruscode/IcarusObj/PMTWaveformTimeCorrection.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

// ROOT libraries
#include "TTree.h"
#include "TFile.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

// -----------------------------------------------------------------------------
namespace icarus {   class OpHitCorrectionsChecker; }
/**
 * @brief Dump corrected and uncorrected OpHit times to a ROOT tree
 *
 * This module reads and saves the corrected/uncorrected OpHit times in
 * a simple ROOT tree for quick access and visualization.
 * 
 * 
 * Input
 * ------
 * 
 * * Two `std::vector<recob::OpHit>` data products
 * 
 * 
 * Output
 * -------
 * 
 * * A ROOT tree, storing ophits per event and channel. 
 * 
 * 
*/

class icarus::OpHitCorrectionsChecker : public art::EDAnalyzer {
public:
  
  struct Config {

  }; // struct Config 
 
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// constructor
  explicit OpHitCorrectionsChecker(Parameters const& config);

  /// process the event
  void analyze(art::Event const& e) override;
  void beginJob();
  void endJob();

private:

  art::ServiceHandle<art::TFileService> tfs;

  /// data members
  TTree *fTree;
  TTree *fTreeTrig;
  int m_run;
  int m_event;
  int m_timestamp;
  int m_channel_id;

  double m_start_time;
  double m_peak_time;
  double m_rise_time;
  double m_start_time_corr;
  double m_peak_time_corr;
  double m_rise_time_corr;
  double m_corr;

  double m_trg_channel;
  double m_wf_corr;

};

// --------------------------------------------------------------------------
icarus::OpHitCorrectionsChecker::OpHitCorrectionsChecker(Parameters const& config)
  : art::EDAnalyzer(config)
{   
    // consumes
    consumes<std::vector<recob::OpHit>>("ophituncorrected");
    consumes<std::vector<recob::OpHit>>("ophit");
}

// ---------------------------------------------------------------------------
void icarus::OpHitCorrectionsChecker::beginJob() {

  fTree = tfs->make<TTree>("ophittree", "ophit info" );
  fTree->Branch("run",&m_run);
  fTree->Branch("event",&m_event);
  fTree->Branch("timestamp",&m_timestamp);
  fTree->Branch("channel_id",&m_channel_id);
  fTree->Branch("tstart",&m_start_time);
  fTree->Branch("tpeak",&m_peak_time);
  fTree->Branch("trise",&m_rise_time);
  fTree->Branch("tstart_corr",&m_start_time_corr);
  fTree->Branch("tpeak_corr",&m_peak_time_corr);
  fTree->Branch("trise_corr",&m_rise_time_corr);
  fTree->Branch("corr",&m_corr);

  fTreeTrig = tfs->make<TTree>("trigtree","trig corr");
  fTreeTrig->Branch("run",&m_run);
  fTreeTrig->Branch("event",&m_event);
  fTreeTrig->Branch("timestamp",&m_timestamp);
  fTreeTrig->Branch("trg_channel",&m_trg_channel);
  fTreeTrig->Branch("wf_corr",&m_wf_corr);
 
}

// ---------------------------------------------------------------------------

void icarus::OpHitCorrectionsChecker::analyze(art::Event const& e)
{
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second    

  // get waveforms 
  auto const& opHitUncorrected = e.getProduct<std::vector<recob::OpHit>>("ophituncorrected");
  auto const& opHit = e.getProduct<std::vector<recob::OpHit>>("ophit");
  auto const& wfCorrections = e.getProduct<std::vector<icarus::timing::PMTWaveformTimeCorrection>>("daqPMT:globtrg");

  std::cout << "Uncorrected size: " << opHitUncorrected.size() << std::endl;
  std::cout << "Corrected size: " << opHit.size() << std::endl;
 
  for (unsigned int k=0; k<wfCorrections.size(); k++){
        auto wfc = wfCorrections.at(k);
	m_trg_channel = wfc.channelID;
        m_wf_corr = wfc.startTime;
        fTreeTrig->Fill();  
  }

  for( unsigned int i=0; i<opHit.size(); i++){

    auto h = opHit.at(i);
    auto uh = opHitUncorrected.at(i);

    if ( h.OpChannel() != uh.OpChannel() ) 
	std::cout << "ERROR: vector mismatch!!!!!" << std::endl;

    m_channel_id = h.OpChannel();

    m_start_time = uh.StartTime();
    m_peak_time = uh.PeakTime();
    m_rise_time = uh.RiseTime();

    m_start_time_corr = h.StartTime();
    m_peak_time_corr = h.PeakTime();
    m_rise_time_corr = h.RiseTime();

    m_corr = m_start_time - m_start_time_corr;

    fTree->Fill();
  }
}

// ----------------------------------------------------------------------------
void icarus::OpHitCorrectionsChecker::endJob()
{
}

// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(icarus::OpHitCorrectionsChecker)
