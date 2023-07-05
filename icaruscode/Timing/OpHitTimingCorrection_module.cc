/**
 * @file   icaruscode/Timing/OpHitTimingCorrection_module.cc
 * @brief  Extract timing correction and adjust the OpHit starting point.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   June 03, 2022
 */

// ICARUS libraries
#include "icaruscode/Timing/PMTTimingCorrections.h"
#include "icaruscode/Timing/IPMTTimingCorrectionService.h"

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
#include "lardataobj/RecoBase/OpHit.h"

// C/C++ standard libraries
#include <optional>
#include <memory> // std::unique_ptr<>
#include <string>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus { class OpHitTimingCorrection; }
/**
 * @brief Creates a new collection of optical hits with corrected time.
 * 
 * This module reads reconstructed optical detector hits and applies time
 * corrections from laser and cosmic ray analyses as found in the
 * `icarusDB::IPMTTimingCorrectionService` service. The corrections are based
 * solely on the PMT channel the hits were reconstructed from.
 * 
 * A new collection of hits is produced containing a corrected copy of all the
 * hits from the input collections.
 * 
 * 
 * Input
 * ------
 * 
 * * `std::vector<recob::OpHit>` data products (as for `InputLabels`)
 * 
 * 
 * Output
 * -------
 * 
 * * a single `std::vector<recob::OpHit>` data product with the hits from the
 *   input collections, all with corrected times; the hits are in the order of
 *   the data products specified in input.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * * `InputLabels` (list of input tags, mandatory): the list of optical hit data
 *   products to apply the time corrections on. It must be non-empty.
 * * `CorrectLaser` (flag, default: `true`): if set, applies the correction
 *   extracted from laser runs.
 * * `CorrectCosmics` (flag, default: `true`): if set, applies the correction
 *   extracted from cosmic ray analysis. This is a correction on top of the one
 *   from the laser, so `CorrectLaser` must also be set.
 * * `Verbose` (flag, default: `false`): prints on screen the corrections being
 *   applied.
 * * `LogCategory` (string, default: `OpHitTimingCorrection`): name of the
 *   message stream for console output.
 * 
 */
class icarus::OpHitTimingCorrection: public art::SharedProducer {
  
public:
  
  /// Configuration of the module.
  struct Config {

    fhicl::Sequence<art::InputTag> InputLabels {
        fhicl::Name("InputLabels"), 
        fhicl::Comment("list of the input lables to be used")
    };

    fhicl::Atom<bool> CorrectLaser {
        fhicl::Name("CorrectLaser"),
        fhicl::Comment("Equalize the PMT timing using the laser"),
        true //default 
    };

    fhicl::Atom<bool> CorrectCosmics {
        fhicl::Name("CorrectCosmics"),
        fhicl::Comment("Fine tune the equalization of the PMT adding the cosmics calibration"),
        true //default 
    };

    fhicl::Atom<bool> Verbose {
      fhicl::Name("Verbose"),
      fhicl::Comment("print the times read and the associated channels"),
      false // default
    };
    
    fhicl::Atom<std::string> LogCategory {
      fhicl::Name("LogCategory"),
      fhicl::Comment("category tag used for messages to message facility"),
      "OpHitTimingCorrection" // default
    };
    
  }; // struct Config

  using Parameters = art::SharedProducer::Table<Config>;

  /// Constructor: just reads the configuration.
  explicit OpHitTimingCorrection(Parameters const& config, art::ProcessingFrame const&);
    
  /// process the event
  void produce(art::Event& event, art::ProcessingFrame const&) override;

private:

  std::vector<art::InputTag> const fInputLabels;

  bool const fCorrectLaser;

  bool const fCorrectCosmics;

  bool const fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string const fLogCategory; ///< Category tag for messages.

  /// Pointer to the online pmt corrections service
  icarusDB::PMTTimingCorrections const& fPMTTimingCorrectionsService;

};


// -----------------------------------------------------------------------------
icarus::OpHitTimingCorrection::OpHitTimingCorrection
    ( Parameters const& config, art::ProcessingFrame const& )
    : art::SharedProducer(config)
    , fInputLabels{ config().InputLabels() }
    , fCorrectLaser{ config().CorrectLaser() }
    , fCorrectCosmics{ config().CorrectCosmics() }
    , fVerbose{ config().Verbose() }
    , fLogCategory{ config().LogCategory() }
    , fPMTTimingCorrectionsService
    { *(lar::providerFrom<icarusDB::IPMTTimingCorrectionService const>()) }
{
    async<art::InEvent>();
    
    // configuration checks
    if (fInputLabels.empty()) {
        throw art::Exception{ art::errors::Configuration }
          << "The list of input hit data products ('"
          << config().InputLabels.name() << "') is empty.\n";
    }
    if (fCorrectCosmics && !fCorrectLaser) {
        throw art::Exception{ art::errors::Configuration }
          << "The corrections from cosmic rays (enabled with '"
          << config().CorrectCosmics.name()
          << "') must be applied on top of the ones from laser runs ('"
          << config().CorrectLaser.name()
          << "'), but the latter are disabled. Fix the configuration.\n";
    }
    if (!fCorrectCosmics && !fCorrectLaser) {
        throw art::Exception{ art::errors::Configuration }
          << "No correction to be applied.\n";
    }

    /// Consumes
    for ( auto const & tag : fInputLabels )
        consumes<std::vector<recob::OpHit>>(tag);

    produces<std::vector<recob::OpHit>>();

}

// -----------------------------------------------------------------------------
void icarus::OpHitTimingCorrection::produce( art::Event& event, art::ProcessingFrame const& ) {

    // Create a copy of the OpHits 
    std::vector<recob::OpHit> correctedOpHits;
    
    auto log = fVerbose? std::make_optional<mf::LogInfo>(fLogCategory): std::nullopt;

    for(art::InputTag const& label: fInputLabels) {
        
        auto const& opHits = event.getProduct<std::vector<recob::OpHit>>(label);

        for(  auto const & opHit : opHits  ){

            double peakTime = opHit.PeakTime();
            double peakTimeAbs = opHit.PeakTimeAbs();
            double startTime = opHit.StartTime();
            double riseTime = opHit.RiseTime();

            double laserTimeCorrection=0;
            double cosmicsCorrection=0;

            if( fCorrectLaser ){
                laserTimeCorrection = 
                    fPMTTimingCorrectionsService.getLaserCorrections(opHit.OpChannel());
            }

            if( fCorrectLaser && fCorrectCosmics ){

                cosmicsCorrection = 
                    fPMTTimingCorrectionsService.getCosmicsCorrections(opHit.OpChannel());
            }

            double const totalCorrection = laserTimeCorrection + cosmicsCorrection;
            peakTime += totalCorrection;
            peakTimeAbs += totalCorrection;
            startTime += totalCorrection;
            // riseTime is currently relative to the start time, no correction needed
            
            if(log){
                *log << opHit.OpChannel() << ", " 
                     << opHit.PeakTime() << ", " 
                     << peakTime << ", " 
                     << cosmicsCorrection << ", " 
                     << laserTimeCorrection
                     << totalCorrection;
            }


            correctedOpHits.emplace_back(
                opHit.OpChannel(),   // channel
                peakTime,            // peaktime
                peakTimeAbs,         // peaktimeabs
                startTime,           // starttime
                riseTime,            // risetime
                opHit.Frame(),       // frame
                opHit.Width(),       // width
                opHit.Area(),        // area
                opHit.Amplitude(),   // peakheight
                opHit.PE(),          // pe
                opHit.FastToTotal()  // fasttototal
            );
        }

    }


    // The new OpHits collection is also saved in the event stream
    event.put(
      std::make_unique<std::vector<recob::OpHit>>(std::move(correctedOpHits)) 
    );

} //icarus::TimingCorrectionExtraction::produce


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::OpHitTimingCorrection)


// -----------------------------------------------------------------------------
