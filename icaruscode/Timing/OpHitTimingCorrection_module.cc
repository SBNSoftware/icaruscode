/**
 * @file   icarus::OpHitTimingCorrection_module.cc
 * @brief  Extract timing correction and adjust the OpHit starting point.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   June 03, 2022
 */

#include "icaruscode/Timing/PMTTimingCorrections.h"
#include "icaruscode/Timing/IPMTTimingCorrectionService.h"

// framework libraries
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"

// Database interface helpers
#include "wda.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>


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
class icarus::OpHitTimingCorrection: public art::EDProducer {
  
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

  using Parameters = art::EDProducer::Table<Config>;

  /// Constructor: just reads the configuration.
  explicit OpHitTimingCorrection(Parameters const& config);
    
  /// process the event
  void produce(art::Event& event ) override;

private:

  std::vector<art::InputTag> fInputLabels;

  bool fCorrectLaser;

  bool fCorrectCosmics;

  bool fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string fLogCategory; ///< Category tag for messages.

  /// Pointer to the online pmt corrections service
  icarusDB::PMTTimingCorrections const& fPMTTimingCorrectionsService;

};


// -----------------------------------------------------------------------------
icarus::OpHitTimingCorrection::OpHitTimingCorrection( Parameters const& config ) 
    : art::EDProducer(config)
    , fInputLabels{ config().InputLabels() }
    , fCorrectLaser{ config().CorrectLaser() } 
    , fCorrectCosmics{ config().CorrectCosmics() } 
    , fVerbose{ config().Verbose() }
    , fLogCategory{ config().LogCategory() }
    , fPMTTimingCorrectionsService
    { *(lar::providerFrom<icarusDB::IPMTTimingCorrectionService const>()) }
{

    /// Consumes
    for ( auto const & tag : fInputLabels )
        consumes<std::vector<recob::OpHit>>(tag);

    produces<std::vector<recob::OpHit>>();

}

// -----------------------------------------------------------------------------
void icarus::OpHitTimingCorrection::produce( art::Event& event ) {

    // Create a copy of the OpHits 
    std::vector<recob::OpHit> correctedOpHits;

    if( !fInputLabels.empty() ){
        
        for( size_t iTag=0; iTag<fInputLabels.size(); iTag++ ){

            art::InputTag label = fInputLabels[iTag];

            art::Handle<std::vector<recob::OpHit>> opHitHandle;
            event.getByLabel( label, opHitHandle );

            if( opHitHandle.isValid() && !opHitHandle->empty() ){

                for(  auto const & opHit : *opHitHandle  ){

                    double peakTime = opHit.PeakTime();
                    double peakTimeAbs = opHit.PeakTimeAbs();

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

                        
                    peakTime += (laserTimeCorrection+cosmicsCorrection);
                    peakTimeAbs += (laserTimeCorrection+cosmicsCorrection);
                    
                    if(fVerbose){
                        std::cout << opHit.OpChannel() << ", " 
                              << opHit.PeakTime() << ", " 
                              << peakTime << ", " 
                              << cosmicsCorrection << ", " 
                              << laserTimeCorrection << std::endl;
                    }


                    correctedOpHits.emplace_back(
                        opHit.OpChannel(),  //
                        peakTime, //
                        peakTimeAbs, //
                        opHit.Frame(), //
                        opHit.Width(), //
                        opHit.Area(), //
                        opHit.Amplitude(), //
                        opHit.PE(), //
                        0.0 //
                    );
                }

            } else {
                mf::LogError(fLogCategory) 
                    << " Not found recob::OpHit data product with label '" 
                    << label.encode() << "'" << std::endl;
                throw;
            }

        }

    } else {
        mf::LogError(fLogCategory) 
            << " InputLabels array should contain more than 1 valid entry " << std::endl;
        throw;
    }


    // The new OpHits collection is also saved in the event stream
    event.put(
      std::make_unique<std::vector<recob::OpHit>>(std::move(correctedOpHits)) 
    );

} //icarus::TimingCorrectionExtraction::produce


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::OpHitTimingCorrection)


// -----------------------------------------------------------------------------
