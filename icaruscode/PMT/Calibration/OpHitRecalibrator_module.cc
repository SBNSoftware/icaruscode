/**
 * @file   icaruscode/PMT/Calibration/OpHitRecalibrator_module.cc
 * @brief  Recalibrate PEs and times of optical hits.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date   December 03, 2025
 */

// ICARUS libraries
#include "icaruscode/Timing/PMTTimingCorrections.h"
#include "icaruscode/Timing/IPMTTimingCorrectionService.h"
#include "icaruscode/Timing/PMTTimingCorrectionsProvider.h"

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
namespace icarus
{
    class OpHitRecalibrator;
}
/**
 * @brief Creates a new collection of optical hits with re-calibrated PEs and times.
 *
 * This module reads reconstructed optical detector hits, removes previously-applied
 * gain and/or timing calibrations, and applies new ones.
 * A new collection of hits is produced containing a re-calibrated copy of all the
 * hits from the input collections.
 *
 * Input
 * ------
 * * `std::vector<recob::OpHit>` data products (as for `InputLabels`)
 *
 * Output
 * -------
 * * a single `std::vector<recob::OpHit>` data product with the recalibrated hits from the
 *   input collections; the hits are in the order of the data products specified in input.
 *
 * Configuration parameters
 * -------------------------
 *
 * * `InputLabels` (list of input tags, mandatory): the list of optical hit data
 *   products to recalibrated It must be non-empty.
 * * `RecalibratePE`   (flag, mandatory): if set, recalibrate hit PE values.
 * * `UseGainDatabase` (flag, default: true): if set, use gain values from database
 *    to re-calibrate hit PEs from hit area.
 * * `SPEArea` (double, default: -1): if not using the gain database, single-photoelectron
 *    area in ADC x tick to be used in the PE calibration.
 * * `RecalibrateTime` (flag, mandatory): if set, recalibrate hit times.
 * * `Verbose` (flag, default: `false`): verbose printing
 *
 */
class icarus::OpHitRecalibrator : public art::SharedProducer
{

public:
    /// Constructor: just reads the configuration.
    explicit OpHitRecalibrator(fhicl::ParameterSet const &config, art::ProcessingFrame const &);

    /// process the event
    void produce(art::Event &event, art::ProcessingFrame const &) override;

    /// configuration at every run change
    void beginRun(art::Run &run, art::ProcessingFrame const &) override;

private:
    art::InputTag const fInputLabel;
    bool const fRecalibratePE;
    bool const fRecalibrateTime;
    bool const fUseGainDatabase;
    double const fSPEArea;
    bool const fVerbose = false; ///< Whether to print the configuration we read.

    /// Pointer to the online pmt corrections service
    icarusDB::PMTTimingCorrections const &fPMTTimingCorrectionsService;

    /// Pointer to the provider for the old pmt corrections
    std::unique_ptr<icarusDB::PMTTimingCorrectionsProvider> fOldTimingProvider;
};

// -----------------------------------------------------------------------------
icarus::OpHitRecalibrator::OpHitRecalibrator(fhicl::ParameterSet const &config, art::ProcessingFrame const &)
    : art::SharedProducer(config),
      fInputLabel{config.get<art::InputTag>("InputLabel")},
      fRecalibratePE{config.get<bool>("RecalibratePE")},
      fRecalibrateTime{config.get<bool>("RecalibrateTime")},
      fUseGainDatabase{config.get<bool>("UseGainDatabase", false)},
      fSPEArea{config.get<double>("SPEArea", -1.)},
      fVerbose{config.get<bool>("Verbose", false)},
      fPMTTimingCorrectionsService{*(lar::providerFrom<icarusDB::IPMTTimingCorrectionService const>())},
      fOldTimingProvider{std::make_unique<icarusDB::PMTTimingCorrectionsProvider>(config.get<fhicl::ParameterSet>("OldTimingDBTags"))}
{
    async<art::InEvent>();

    // configuration checks
    if (fInputLabel.empty())
    {
        throw art::Exception{art::errors::Configuration}
            << "The list of input hit data products ('InputLabel') is empty.\n";
    }

    if (!fRecalibratePE && !fRecalibrateTime)
    {
        throw art::Exception{art::errors::Configuration}
            << "No re-calibration selected. Why are you running meeee!?!?! :/\n";
    }

    if (fRecalibratePE && !fUseGainDatabase && (fSPEArea < 0))
    {
        throw art::Exception{art::errors::Configuration}
            << "The gain database for PE recalibration has been disabled ('UseGainDatabase'),"
            << " but 'SPEArea' (" << fSPEArea
            << ") has not been set correctly. Fix the configuration!\n";
    }

    if (fRecalibrateTime)
    {
        mf::LogInfo("OpHitRecalibrator") << "Re-calibration of timing corrections enabled:\n"
                                         << "LaserTag: old " << fOldTimingProvider->getLaserDatabaseTag()
                                         << " -> new " << fPMTTimingCorrectionsService.getLaserDatabaseTag()
                                         << "\nCosmicsTag: old " << fOldTimingProvider->getCosmicsDatabaseTag()
                                         << " -> new " << fPMTTimingCorrectionsService.getCosmicsDatabaseTag();
    }

    // FIXME: temporary since no gain db exists yet...
    if (fUseGainDatabase)
    {
        throw art::Exception{art::errors::Configuration}
            << "Gain database interface doesn't exist yet. Try again later.\n";
    }

    // Consumes
    consumes<std::vector<recob::OpHit>>(fInputLabel);
    // Produces
    produces<std::vector<recob::OpHit>>();
}

// -----------------------------------------------------------------------------
void icarus::OpHitRecalibrator::beginRun(art::Run &run, art::ProcessingFrame const &)
{
    // make sure we are updating the locally-instantiated correction provider
    if (fOldTimingProvider)
    {
        if (fVerbose)
            mf::LogInfo("OpHitRecalibrator") << "Updating old timing corrections for " << run.id().run();
        fOldTimingProvider->readTimeCorrectionDatabase(run);
    }
}

// -----------------------------------------------------------------------------
void icarus::OpHitRecalibrator::produce(art::Event &event, art::ProcessingFrame const &)
{
    auto log = fVerbose ? std::make_optional<mf::LogInfo>("OpHitRecalibrator") : std::nullopt;

    // Create a copy of the OpHits
    std::vector<recob::OpHit> recalibratedOpHits;

    auto const &opHits = event.getProduct<std::vector<recob::OpHit>>(fInputLabel);

    for (auto const &opHit : opHits)
    {
        // read current times
        double peakTime = opHit.PeakTime();
        double peakTimeAbs = opHit.PeakTimeAbs();
        double startTime = opHit.StartTime();
        double riseTime = opHit.RiseTime();

        // read current PE
        double hitPE = opHit.PE();

        // First, recalibrate PE values (if enabled).
        // - if using gain database, fetch new SPEArea for this run/channel
        // - if not using gain database, use fSPEArea
        // - re-compute PE value with new SPE area
        if (fRecalibratePE)
        {
            double oldSPEArea = opHit.Area() / hitPE;
            double newSPEArea = fSPEArea;

            if (fUseGainDatabase)
            {
                // soon...
                // newSPEArea = get_from_db(opHit.OpChannel())
            }

            if (log)
            {
                *log << "Channel: " << opHit.OpChannel()
                     << ", Area: " << opHit.Area() << " [ADC x tick]"
                     << ", PE " << hitPE
                     << " (old SPEArea: " << oldSPEArea
                     << ") --> new PE " << opHit.Area() / newSPEArea
                     << " (new SPEArea: " << newSPEArea << ")\n";
            }

            hitPE = opHit.Area() / newSPEArea;
        }

        // Second, recalibrate PMT times (if enabled)
        // Previous corrections need to be subtracted
        // use locally instantiated provider class to fetch old corrections
        // use online pmt corrections servicer to fetch current/new corrections
        if (fRecalibrateTime)
        {
            // get the old timing corrections: these will need to be subtracted!
            double const oldLaserTimeCorrection = fOldTimingProvider->getLaserCorrections(opHit.OpChannel());
            double const oldCosmicsCorrection = fOldTimingProvider->getCosmicsCorrections(opHit.OpChannel());
            double const oldTotalCorrection = oldLaserTimeCorrection + oldCosmicsCorrection;

            // get new/current timing: these will need to be added!
            double const laserTimeCorrection = fPMTTimingCorrectionsService.getLaserCorrections(opHit.OpChannel());
            double const cosmicsCorrection = fPMTTimingCorrectionsService.getCosmicsCorrections(opHit.OpChannel());
            double const totalCorrection = laserTimeCorrection + cosmicsCorrection;

            double timeCorr = totalCorrection - oldTotalCorrection;

            if (log)
            {
                *log << "Channel: " << opHit.OpChannel()
                     << ", startTime " << startTime
                     << ", peakTime " << peakTime
                     << ", peakTimeAbs " << peakTimeAbs
                     << " (total correction: " << oldTotalCorrection
                     << ") --> new startTime " << startTime + timeCorr
                     << ", peakTime " << peakTime + timeCorr
                     << ", peakTimeAbs " << peakTimeAbs + timeCorr
                     << " (total correction: " << totalCorrection << ")\n";
            }

            peakTime = peakTime + timeCorr;
            peakTimeAbs = peakTimeAbs + timeCorr;
            startTime = startTime + timeCorr;
            // NOTE: riseTime is currently relative to the start time, so no correction needed
        }

        recalibratedOpHits.emplace_back(
            opHit.OpChannel(),  // channel
            peakTime,           // recalibrated peaktime
            peakTimeAbs,        // recalibrated peaktimeabs
            startTime,          // recalibrated starttime
            riseTime,           // recalibrated risetime
            opHit.Frame(),      // frame
            opHit.Width(),      // width
            opHit.Area(),       // area
            opHit.Amplitude(),  // peakheight
            hitPE,              // recalibrated pe
            opHit.FastToTotal() // fasttototal
        );
    }

    // The new OpHits collection is also saved in the event stream
    event.put(
        std::make_unique<std::vector<recob::OpHit>>(std::move(recalibratedOpHits)));

} // icarus::OpHitRecalibrator::produce

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::OpHitRecalibrator)

// -----------------------------------------------------------------------------
