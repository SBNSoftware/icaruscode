/**
 * @file   icaruscode/PMT/OpDetWaveformFiltering_module.cc
 * @brief  Apply a low pass filter on PMT waveforms.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date   May 02, 2023
 */

// framework libraries
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <string>
#include <utility> // std::move()


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
 * 
 */
class icarus::OpDetWaveformFiltering: public art::SharedProducer {
  
public:
  
  /// configuration of the module.
  struct Config {

    fhicl::Sequence<art::InputTag> InputLabels {
        fhicl::Name("InputLabels"), 
        fhicl::Comment("list of the input labels to be used")
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
    
  }; // struct Config

  using Parameters = art::SharedProducer::Table<Config>;
  using nanoseconds = util::quantities::intervals::nanoseconds;

  /// constructor: just reads the configuration
  explicit OpDetWaveformFiltering(Parameters const& config, art::ProcessingFrame const&);
    
  /// process the event
  void produce(art::Event& event, art::ProcessingFrame const&) override;

  /// low pass discrete filter
  template<typename T> std::vector<T> lowPassFilter(std::vector<T> data, double dt, double RC) const;

private:

  std::vector<art::InputTag> const fInputLabels;
  double const fTimeConstant; ///< RC time constant
  std::string const fLogCategory; ///< Category tag for messages.

};


// -----------------------------------------------------------------------------
icarus::OpDetWaveformFiltering::OpDetWaveformFiltering
    ( Parameters const& config, art::ProcessingFrame const& )
    : art::SharedProducer(config)
    , fInputLabels{ config().InputLabels() }
    , fTimeConstant{ config().TimeConstant() }
    , fLogCategory{ config().LogCategory() }
{
    async<art::InEvent>();
    
    // configuration checks
    if (fInputLabels.empty()) {
        throw art::Exception{ art::errors::Configuration }
          << "The list of input waveform data products ('"
          << config().InputLabels.name() << "') is empty.\n";
    }
    if (fTimeConstant < 0) {
        throw art::Exception{ art::errors::Configuration }
          << "The RC time constant set with '"
          << config().TimeConstant.name()
          << "') must be a positive number. Fix the configuration!\n";
    }

    /// Consumes
    for ( auto const & tag : fInputLabels )
        consumes<std::vector<raw::OpDetWaveform>>(tag);

    produces<std::vector<raw::OpDetWaveform>>();

}

// -----------------------------------------------------------------------------
void icarus::OpDetWaveformFiltering::produce( art::Event& event, art::ProcessingFrame const& ) {

    // get optical detector sampling time
    detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event));
    nanoseconds opticalTick = detTimings.OpticalClockPeriod();
    std::cout << "Tick [hopefully in ns] " << opticalTick.value() << std::endl;
    std::cout << "Time constant [ns] " << fTimeConstant << std::endl;

    // create a copy of the waveforms 
    std::vector<raw::OpDetWaveform> filteredOpDetWfs;

    for(art::InputTag const& label: fInputLabels) {
        
        auto const& opDetWfs = event.getProduct<std::vector<raw::OpDetWaveform>>(label);

        for(  auto const & wf : opDetWfs  ){

            double channel = wf.ChannelNumber();
            double time = wf.TimeStamp();

	    // filtering waveform data
	    std::vector<raw::ADC_Count_t> filt = lowPassFilter( wf.Waveform(), opticalTick.value(), fTimeConstant);	      
            std::vector<uint16_t> filtwf( filt.begin(), filt.end());

	    // fill new collection
	    filteredOpDetWfs.emplace_back(
                time,               // timestamp
                channel,            // channel
                filtwf              // waveform
            );
        }
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

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::OpDetWaveformFiltering)


// -----------------------------------------------------------------------------
