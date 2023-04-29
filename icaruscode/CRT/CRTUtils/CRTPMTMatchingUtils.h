#include "icaruscode/IcarusObj/CRTPMTMatching.h"

namespace icarus::crt {

  struct CRTPMT {
    double tof;    ///< Time difference between CRT Hit and optical flash [ns]
    double distance;    ///< Distance between CRT Hit and optical flash centroid [cm]
    art::Ptr<sbn::crt::CRTHit> CRTHit;
  };

  struct CRTMatches {
    // vector of pairs where first is the matched Time of Flight and second is the
    // matched CRTHit
    std::vector<CRTPMT> entering;
    std::vector<CRTPMT> exiting;
    enum matchType FlashType;
  };

  struct FlashType {
    geo::Point_t FlashPos;
    double FlashTime_us;
    double FlashGateTime_ns;
    bool inBeam;
    bool inGate;
    enum matchType Classification;
    std::vector<MatchedCRT> CRTmatches;
  };

bool flashInTime(double flashTime, int gateType, double gateDiff,
                 double gateWidth) {
  // As a reminder, I will leave here the commented part of the vetoOffset, in
  // case something changes in the future
  /*double vetoOffset = 3.5;*/
  double activeGate = gateWidth /*- vetoOffset*/;

  double relFlashTime = flashTime + gateDiff / 1000. /*- vetoOffset*/;
  mf::LogDebug("CRTPMTMatchingProducer FlashInTime")
      << "Gate Diff " << gateDiff / 1000 << " Ftime+gateDiff "
      << flashTime + gateDiff / 1000. << " " << activeGate;

  return ((relFlashTime > 0) && (relFlashTime < activeGate));
}

CRTMatches CRTHitmatched(
    double flashTime, geo::Point_t const& flashpos,
    std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, double interval) {

  	std::vector<icarus::crt::CRTPMT> enteringCRTHits;
  	std::vector<icarus::crt::CRTPMT> exitingCRTHits;
  	enum matchType FlashType;
  	int topen = 0, topex = 0, sideen = 0, sideex = 0;
  	for (auto const& crtHit : crtHits) {
    		double tof = crtHit->ts1_ns - flashTime * 1e3;
    		double distance =
		  (flashpos - geo::Point_t{crtHit->x_pos, crtHit->y_pos, crtHit->z_pos})
            	    .R();
    		if (abs(tof) >= interval) continue;
    		if (tof < 0) {
      			if (crtHit->plane > 36)
        			sideen++;
      			else
        			topen++;
      			CRTPMT m_match = {tof, distance, crtHit};
      			enteringCRTHits.push_back(m_match);
    		} else if (tof >= 0) {
      			if (crtHit->plane > 36)
        			sideex++;
      			else
        			topex++;
      			CRTPMT m_match = {tof, distance, crtHit};
      			exitingCRTHits.push_back(m_match);
    		}
  	}
  	if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0)
    		FlashType = noMatch;
  	else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0)
    		FlashType = enTop;
  	else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0)
    		FlashType = enSide;
  	else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1)
    		FlashType = enTop_exSide;
  	else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0)
    		FlashType = exTop;
  	else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1)
    		FlashType = exSide;
  	else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0)
    		FlashType = enTop_mult;
  	else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1)
    		FlashType = enTop_exSide_mult;
  	else
    		FlashType = others;

  	return {std::move(enteringCRTHits), std::move(exitingCRTHits), FlashType};
}


  CRTPMTMatching FillCRTPMT (FlashType thisFlash, int event, int run, int gate){
    CRTPMTMatching crtpmt;
    crtpmt.event = event;
    crtpmt.run = run;
    crtpmt.gateType = gate;
    crtpmt.flashID = 0; //
    crtpmt.flashTime_us = thisFlash.FlashTime_us;
    crtpmt.flashGateTime_ns = thisFlash.FlashGateTime_ns;
    crtpmt.firstOpHitPeakTime_us = 0; //
    crtpmt.firstOpHitStartTime_us = 0; //
    crtpmt.flashInGate = thisFlash.inGate;
    crtpmt.flashInBeam = thisFlash.inBeam;
    crtpmt.flashAmplitude_pe = 0; //
    crtpmt.flashPosition = thisFlash.FlashPos;
    crtpmt.flashYWidth = 0; //
    crtpmt.flashZWidth = 0; //
    crtpmt.flashClassification = thisFlash.Classification;
    for (auto const& crts : thisFlash.CRTmatches){
      MatchedCRT thismatchedCRT;
      /*thismatchedCRT.CRTHitModule = 0; //
      thismatchedCRT.CRTRegion = crts.CRTRegion;
      thismatchedCRT.CRTSys = crts.CRTSys;
      thismatchedCRT.CRTHitPosition = crts.CRTHitPosition;
      thismatchedCRT.CRTHitTime_us = crts.CRTHitTime_us;
      thismatchedCRT.CRTHitGateTime_ns = 0;//
      thismatchedCRT.CRTHitAmplitude_pe = 0;//
      thismatchedCRT.CRTPMTTimeDiff_ns = crts.CRTPMTTimeDiff_ns;
      thismatchedCRT.CRTHitFlashDistance = 0; //*/
      // push stuff back to MatchedCRT (for now onlt CRTHitPos, CRTPMTTimeDiff, CRTTime, CRTSys, CRTReg)
      thismatchedCRT.CRTHitPos = crts.CRTHitPos;
      thismatchedCRT.CRTPMTTimeDiff_ns = crts.CRTPMTTimeDiff_ns; 
      thismatchedCRT.CRTTime_us = crts.CRTTime_us; 
      thismatchedCRT.CRTSys = crts.CRTSys;
      thismatchedCRT.CRTRegion = crts.CRTRegion; 
    }
    crtpmt.topCRTBefore = 0; //
    crtpmt.topCRTAfter = 0; //
    crtpmt.sideCRTBefore = 0; //
    crtpmt.sideCRTAfter = 0; //
    
    return crtpmt;
    }
}

