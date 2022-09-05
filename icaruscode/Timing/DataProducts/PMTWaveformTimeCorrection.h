/**
 * @file   icaruscode/Timing/DataProducts/PMTWaveformTimeCorrection.h
 * @brief  Holds the event-by-event waveform timing adjustment.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   March 23 2021
 * @see    icaruscode/Timing/DataProducts/PMTWaveformTimeCorrection.cxx
 */

#ifndef ICARUSCODE_TIMING_DATAPRODUCTS_PMTWAVEFORMTIMECORRECTION_H
#define ICARUSCODE_TIMING_DATAPRODUCTS_PMTWAVEFORMTIMECORRECTION_H

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <string>
#include <array>
#include <cassert>

namespace icarus::timing{

	struct PMTWaveformTimeCorrection {
		
		unsigned int channelID;
		
		double sample;

		double startTime;

	};

} // namespace icarus::timing

#endif //ICARUSCODE_TIMING_DATAPRODUCTS_PMTWAVEFORMTIMECORRECTION_H