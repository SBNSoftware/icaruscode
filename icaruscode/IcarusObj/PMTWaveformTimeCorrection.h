/**
 * @file   icaruscode/IcarusObj/PMTWaverformTimeCorrection.h
 * @brief  Holds the event-by-event waveform timing adjustment.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   March 23 2021
 * @see    icaruscode/Timing/DataProducts/PMTWaveformTimeCorrection.cxx
 */

#ifndef ICARUSCODE_ICARUSOBJ_PMTWAVEFORMTIMECORRECTION_H
#define ICARUSCODE_ICARUSOBJ_PMTWAVEFORMTIMECORRECTION_H

// C/C++ standard libraries
#include <limits>

namespace icarus::timing{

	/**
	 * @brief Corrections to the PMT waveform time.
	 * 
	 * An (uncorrected) waveform expects that the trigger signal happened at the
	 * nominal `TriggerTime()` from `detinfo::DetectorClocks`.
	 * 
	 * If that is not the case, the `startTime` correction can be used to shift
	 * the waveform start time so that that assumption holds.
	 * The new, corrected waveform is the old one _plus_ the correction.
	 * 
	 */
	struct PMTWaveformTimeCorrection {
		
		/// Special value to denote no channel ID information (and no correction).
		static constexpr auto NoChannelID = std::numeric_limits<unsigned int>::max();
		
		/// The channel this correction was extracted from.
		unsigned int channelID = NoChannelID;
		
		/// Time within the reference waveforms where the reference signal is found [ns]
		double sample = std::numeric_limits<double>::lowest();

		/// How earlier than the nominal trigger time the actual trigger signal
		/// seems to come [us]
		double startTime = 0.0;
		
		
		/// Returns whether the correction values are valid.
		bool isValid() const { return channelID != NoChannelID; }

	};

} // namespace icarus::timing

#endif //ICARUSCODE_ICARUSOBJ_PMTWAVEFORMTIMECORRECTION_H
