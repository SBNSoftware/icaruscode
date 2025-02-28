/**
 * @file   icaruscode/Timing/PMTWaveformTimeCorrectionExtractor.h
 * @brief  Extract timing correction and adjust waveform starting point.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   June 03, 2022
 */

#ifndef ICARUSCODE_TIMING_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H
#define ICARUSCODE_TIMING_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H

// ICARUS/SBN libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Timing/PMTTimingCorrections.h"
#include "icaruscode/IcarusObj/PMTWaveformTimeCorrection.h"
#include "icaruscode/Timing/Tools/PulseStartExtractor.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>

namespace icarus::timing { class PMTWaveformTimeCorrectionExtractor; }

/**
 * @brief Extracts timing corrections from a waveform.
 * 
 * This algorithm extracts a time correction from a reference signal, and
 * associates that correction with all the channels in the same PMT readout
 * crate.
 * 
 * The correction extraction is performed by `findWaveformTimeCorrections()`.
 * That method analyzes a reference waveform searching for a reference signal
 * assumed to have been generated in time with the global trigger, and emits
 * a correction so that this signal would appear at exactly the time of the
 * global trigger (`detinfo::DetectorClocks::TriggerTime()`).
 * 
 * 
 * Signal timing extraction
 * -------------------------
 * 
 * The reference signal is expected to be a sharp square wave in negative
 * polarity. The time of the correction is based on the front side of that wave.
 * The extraction is perfomed by `icarus::timing::PulseStartExtractor` by
 * specifying one of the available extraction methods (constant-fraction discrimation,
 * logitstic function fit) and an ADC threshold for the signal identification.
 * 
 * @see `icarus::timing::PulseStartExtractor` for details.
 * 
 */
class icarus::timing::PMTWaveformTimeCorrectionExtractor {
	
	public: 

        // --- BEGIN -- Exceptions ---------------------------------------------

        /// Exception thrown when trying to overwrite a correction.
        struct Error: cet::exception { Error(std::string const& msg = ""); };

        /// Exception thrown when trying to overwrite a correction.
        struct MultipleCorrectionsForChannel: Error {
            MultipleCorrectionsForChannel(
                unsigned int channel,
                unsigned int existing, unsigned int additional
            );
        };

        /// Exception thrown when correction requested from a non-special channel.
        struct NotASpecialChannel: Error {
            NotASpecialChannel(unsigned int channel);
            private: static Error makeBaseException(unsigned int channel);
        };

        /// Exception thrown when PMT readout crate not recognised.
        struct UnknownCrate: Error {
            UnknownCrate(unsigned int channel);
            private: static Error makeBaseException(unsigned int channel);
        };

        /// Exception thrown when no signal is found above threshold
        struct NoSignalFound: Error {
            NoSignalFound(unsigned int channel, double threshold);
            private: static Error makeBaseException(unsigned int channel, double threshold);
        };

        // --- END ---- Exceptions ---------------------------------------------


        PMTWaveformTimeCorrectionExtractor(
            detinfo::DetectorClocksData const detTimingService,
            icarusDB::IICARUSChannelMap const & channelMapService,
            icarusDB::PMTTimingCorrections const* pmtTimingCorrectionsService,
            icarus::timing::ExtractionMethod const pulseStartExtractionMethod,
            double threshold,
            bool verbose );

        /**
         * @brief Extracts a correction from `wave` and assigns it to channels.
         * @param wave the reference waveform to extract the correction from
         * @param correctCableDelay whether to apply the correction for cable delays
         * @param[in,out] corrections where to add the newly extracted corrections
         * @throw NotASpecialChannel if `wave` is not a special channel
         * @throw UnknownCrate if `wave` is from an unexpected readout crate
         * @throw MultipleCorrectionsForChannel if `corrections` already contains
         *   a correction for any of the channels we are associating to the new correction
         * @throw NoSignalFound if `wave` does not contain any pulse above threshold
         * 
         * This function performs the analysis of the reference waveform `wave`,
         * extracts the correction out of it and associates it in `corrections`
         * for all the channels that are in the same readout crate as the
         * reference waveform itself.
         * 
         * The `corrections` list is extended as needed. Since `corrections`
         * is a dense structure with as index the channel ID, it may happen that
         * the extension produces correction information for channels not in
         * this readout crate: in this case, these correction values are invalid
         * (as in `PMTTimingCorrections::isValid()` returning `false`).
         * Any attempt to overwrite a valid correction already in `corrections`
         * will cause throwing an exception.
         * 
         * The channel ID in each stored correction is the one of the reference
         * waveform the correction was extracted from, and likewise the `sample`
         * data member is where the reference signal can be found within that
         * waveform. The `startTime` correction is an offset to be _added_ to
         * the waveform timestamps to correct them.
         * 
         */
        void findWaveformTimeCorrections(   
            raw::OpDetWaveform const & wave,
            bool correctCableDelay,
            std::vector<PMTWaveformTimeCorrection> & corrections ) const;


	private:

        detinfo::DetectorClocksData const fClocksData;

        icarusDB::IICARUSChannelMap const & fChannelMap;

        icarusDB::PMTTimingCorrections const* fPMTTimingCorrectionsService = nullptr;

        icarus::timing::ExtractionMethod const fPulseStartExtractionMethod;

        double const fPulseStartExtractionThreshold;

        bool const fVerbose;

        std::map<unsigned int, std::vector<unsigned int>> const
            fCrateFragmentMap {
                {0x0070 , { 0 , 1 , 2  }}, 
                {0x0060 , { 3 , 4 , 5  }}, 
                {0x0050 , { 6 , 7 , 8  }}, 
                {0x0040 , { 9 , 10, 11 }},
                {0x0030 , { 12, 13, 14 }}, 
                {0x0020 , { 15, 16, 17 }}, 
                {0x0010 , { 18, 19, 20 }}, 
                {0x0000 , { 21, 22, 23 }}, 
            }; 

};


#endif //ICARUSCODE_TIMING_PMTWAVEFORMTIMECORRECTIONEXTRACTOR_H
