/**
 * @file   icaruscode/IcarusObj/PMTBeamSignal.h
 * @brief  Holds the event-by-event RWM or EW time
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date   March 14 2024
 */

#ifndef ICARUSCODE_ICARUSOBJ_PMTBEAMSIGNAL_H
#define ICARUSCODE_ICARUSOBJ_PMTBEAMSIGNAL_H

// C/C++ standard libraries
#include <limits>

namespace icarus::timing{

	struct PMTBeamSignal {
		
		/// Special value to denote no special channel information
		static constexpr auto NoChannel = std::numeric_limits<unsigned int>::max();
		
		/// The special channel this time was extracted from
		unsigned int channel = NoChannel;
	 	/// Board on which the special channel is on
                std::string digitizerLabel = "";    	
	 	/// Crate this time applies to
                std::string crate = "";    
		/// Sample within the waveform where the reference signal is found
		size_t sample = std::numeric_limits<size_t>::lowest();
		/// Start time w.r.t. trigger time [us]
		double startTime = 0.0;
		
                PMTBeamSignal(unsigned int ch, std::string b, std::string c,
			      size_t s, double t):
			     channel(ch), digitizerLabel(b), crate(c), sample(s),
			     startTime(t) {};

		/// Returns whether the time is valid.
		bool isValid() const { return channel != NoChannel; }

	};

} // namespace icarus::timing

#endif //ICARUSCODE_ICARUSOBJ_PMTBEAMSIGNAL_H
