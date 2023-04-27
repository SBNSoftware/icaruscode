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
    matchType FlashType;
  };

  struct MatchedCRT {
    geo::Point_t CRTHitPos;
    double CRTPMTTimeDiff_ns;
    double CRTTime_us;
    int CRTSys;
    int CRTRegion;
  };

  struct FlashType {
    geo::Point_t FlashPos;
    double FlashTime_us;
    double FlashGateTime_ns;
    bool inBeam;
    bool inGate;
    matchType Classification;
    std::vector<MatchedCRT> CRTmatches;
  };



  CRTPMTMatching FillCRTPMT (FlashType thisFlash, int event, int run, int gate); /*{ // not yet sure how we expect to fill this
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
    }*/
	/*struct FlashType {
	  geo::Point_t FlashPos;
	  double FlashTime_us;
	  double FlashGateTime_ns;
	  bool inBeam;
	  bool inGate;
	  matchType Classification;
	  std::vector<MatchedCRT> CRTmatches;
	  };*/

}

