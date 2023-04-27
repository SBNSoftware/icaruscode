#include "icaruscode/IcarusObj/CRTPMTMatching.h"

namespace icarus::crt {

	CRTPMTMatching FillCRTPMT (FlashType thisFlash, int event, int run, int gate) {
		CRTPMTMatching crtpmt;
		crtpmt.event = event;
		crtpmt.run = run;
		crtpmt.gateType = gate;
		crtpmt.flashID; //
		crtpmt.flashTime_us = thisFlash.FlashTime_us;
		crtpmt.flashGateTime_ns = thisFlash.FlashGateTime_ns;
		crtpmt.firstOpHitPeakTime_us; //
		crtpmt.firstOpHitStartTime_us; //
		crtpmt.flashInGate = thisFlash.inGate;
		crtpmt.flashInBeam = thisFlash.inBeam;
		crtpmt.flashAmplitude_pe; //
		crtpmt.flashPosition = thisFlash.FlashPos;
		crtpmt.flashYWidth; //
		crtpmt.flashZWidth; //
		crtpmt.flashClassification = thisFlash.matchType;
		for (auto const& crts : thisFlash.CRTmatches){
			matchedCRT thismatchedCRT;
			thismatchedCRT.CRTHitModule = ; //
			thismatchedCRT.CRTRegion = crts.CRTRegion;
			thismatchedCRT.CRTSys = crts.CRTSys;
			thismatchedCRT.CRTHitPosition = crts.CRTHitPos;
			thismatchedCRT.CRTHitTime_us = crts.CRTTime_us;
			thismatchedCRT.CRTHitGateTime_ns = ;//
			thismatchedCRT.CRTHitAmplitude_pe = ;//
			thismatchedCRT.CRTPMTTimeDiff_ns = crts.CRTPMTTimeDiff_ns;
			thismatchedCRT.CRTHitFlashDistance = ; //
		}
		
		//crtpmt.matchedCRTHits = thisFlash.CRTmatches.size();
		crtpmt.topCRTBefore; //
		crtpmt.topCRTAfter; //
		crtpmt.sideCRTBefore; //
		crtpmt.sideCRTAfter; //
		
		return crtpmt;
	}
}
