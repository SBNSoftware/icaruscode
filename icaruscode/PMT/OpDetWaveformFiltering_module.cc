/**
 * @file   icaruscode/PMT/OpDetWaveformFiltering_module.cc
 * @brief  Apply a low pass filter on PMT waveforms.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date   May 02, 2023
 */

// framework libraries
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time

// ROOT libraries
#include "TTree.h"
#include "TFile.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <string>
#include <utility> // std::move()
#include <vector>

// -----------------------------------------------------------------------------
namespace icarus { class OpDetWaveformFiltering; }
/**
 * @brief Creates a new collection of filtered optical waveforms.
 * 
 * This module reads the raw optical waveforms and applies a low pass
 * filter to smooth them out and reduce noise, helping the downstream hit finding.
 * 
 * A new collection of optical waveforms is produced containing a filtered
 * copy of all the waveforms from the input collections.
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
 * * a single `std::vector<raw::OpDetWaveform>` data product with the filtered
 *   waveforms from the input collections; the waveforms are in the same order of
 *   the data products specified in input.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * * `InputLabels` (list of input tags, mandatory): the list of optical waveforms data
 *   products to apply the filter on. It must be non-empty.
 * * `TimeConstant` (double, default: `20`): RC time constant for the low pass filter,
 *   which determines the smoothing factor. The expected unit is ns. 
 *   from the laser, so `CorrectLaser` must also be set.
 * * `LogCategory` (string, default: `OpHitWaveformFiltering`): name of the
 *   message stream for console output.
 * * `SaveDebug` (bool, default: `false`): create and save ROOT tree with debugging
 *   information, with raw and filtered waveforms.
 * * `OutputFile` (string, defaulf: `filter_debug.root`): name of the debug output
 *   ROOT file, created only if `SaveDebug` is set to true.
 */
class icarus::OpDetWaveformFiltering: public art::EDProducer {
  
public:
  
  /// configuration of the module.
  struct Config {

    fhicl::Atom<art::InputTag> InputLabel {
        fhicl::Name("InputLabel"), 
        fhicl::Comment("Raw waveform input label to be used")
    };

    fhicl::Atom<double> TimeConstant {
        fhicl::Name("TimeConstant"),
        fhicl::Comment("Low pass filter RC time constant"),
        20. //ns, default 
    };

    fhicl::Atom<std::string> LogCategory {
      fhicl::Name("LogCategory"),
      fhicl::Comment("category tag used for messages to message facility"),
      "OpDetWaveformFiltering" // default
    };
    
    fhicl::Atom<bool> SaveDebug {
      fhicl::Name("SaveDebug"),
      fhicl::Comment("save debugging info in ROOT tree"),
      false // default
    };
    
    fhicl::Atom<std::string> OutputFile {
      fhicl::Name("OutputFile"),
      fhicl::Comment("output ROOT tree filename"),
      "filter_debug.root" // default
    };
    
  }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;
  using nanoseconds = util::quantities::intervals::nanoseconds;

  /// constructor: just reads the configuration
  explicit OpDetWaveformFiltering(Parameters const& config);
    
  /// process the event
  void produce(art::Event& event);
  void beginJob();
  void endJob();

  /// low pass discrete filter
  template<typename T> std::vector<T> lowPassFilter(std::vector<T> data, double dt, double RC) const;

private:

  art::InputTag const fInputLabel;
  double const fTimeConstant; ///< RC time constant
  std::string const fLogCategory; ///< Category tag for messages.
  bool const fSaveDebug;
  std::string const fOutputFile;

  TFile *fFile;
  TTree *fTree;
  int m_run;
  int m_event;
  int m_timestamp;
  int m_channel_id;
  double m_tstart;
  int m_nsize;
  std::vector<short> m_wf;
  std::vector<short> m_filt_wf;

};


// -----------------------------------------------------------------------------
icarus::OpDetWaveformFiltering::OpDetWaveformFiltering
    ( Parameters const& config )
    : art::EDProducer(config)
    , fInputLabel{ config().InputLabel() }
    , fTimeConstant{ config().TimeConstant() }
    , fLogCategory{ config().LogCategory() }
    , fSaveDebug{ config().SaveDebug() }
    , fOutputFile{ config().OutputFile() }
{
    
    // configuration checks
    if (fInputLabel.empty()) {
        throw art::Exception{ art::errors::Configuration }
          << "The input waveform data product ('"
          << config().InputLabel.name() << "') is empty.\n";
    }
    if (fTimeConstant < 0) {
        throw art::Exception{ art::errors::Configuration }
          << "The RC time constant set with '"
          << config().TimeConstant.name()
          << "') must be a positive number. Fix the configuration!\n";
    }

    // consumes
    consumes<std::vector<raw::OpDetWaveform>>(fInputLabel);
	
    // produces
    produces<std::vector<raw::OpDetWaveform>>();

}

// -----------------------------------------------------------------------------

void icarus::OpDetWaveformFiltering::beginJob() {

  if(fSaveDebug) {
    
    fFile = new TFile(fOutputFile.c_str(),"RECREATE");
    
    fTree = new TTree("wftree", "waveform info" );
    fTree->Branch("run",&m_run);
    fTree->Branch("event",&m_event);
    fTree->Branch("timestamp",&m_timestamp);
    fTree->Branch("channel_id",&m_channel_id);
    fTree->Branch("tstart",&m_tstart);
    fTree->Branch("nsize",&m_nsize);
    fTree->Branch("wf",&m_wf);
    fTree->Branch("filt_wf",&m_filt_wf);
  }
}

// -----------------------------------------------------------------------------
void icarus::OpDetWaveformFiltering::produce( art::Event& event ) {

    m_run = event.id().run();
    m_event = event.id().event();
    m_timestamp = event.time().timeHigh(); // precision to the second    

    // get optical detector sampling time
    detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event));
    nanoseconds opticalTick = detTimings.OpticalClockPeriod();

    // create a copy of the waveforms 
    std::vector<raw::OpDetWaveform> filteredOpDetWfs;
     
    auto const& opDetWfs = event.getProduct<std::vector<raw::OpDetWaveform>>(fInputLabel);

    for(  auto const & wf : opDetWfs  ){

    	m_channel_id = wf.ChannelNumber();
        m_tstart = wf.TimeStamp();
        m_wf = wf.Waveform();

	// filtering waveform data
	m_filt_wf = lowPassFilter( m_wf, opticalTick.value(), fTimeConstant);	      
        std::vector<uint16_t> filtwf( m_filt_wf.begin(), m_filt_wf.end());

	// fill new collection
	filteredOpDetWfs.emplace_back(
           m_tstart,           // timestamp
           m_channel_id,       // channel
           filtwf              // waveform
        );
	
	if(fSaveDebug) fTree->Fill();

    }
 
    // the filtered waveforms are also saved in the event stream
    event.put(
      std::make_unique<std::vector<raw::OpDetWaveform>>(std::move(filteredOpDetWfs)) 
    );

}

// ----------------------------------------------------------------------------

template<typename T> std::vector<T> icarus::OpDetWaveformFiltering::lowPassFilter(std::vector<T> data, double dt, double RC) const
{
    std::vector<T> out(data.size());
    double alpha = dt/(RC+dt);
    
    out[0] = alpha * data[0];
    
    for( size_t i=1; i<data.size(); i++ )
	out[i] = alpha * data[i] + (1-alpha) * out[i-1];

    return out;
}


void icarus::OpDetWaveformFiltering::endJob()
{
    fFile->cd();
    fTree->Write();
    fFile->Close();
}


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::OpDetWaveformFiltering)


// -----------------------------------------------------------------------------
