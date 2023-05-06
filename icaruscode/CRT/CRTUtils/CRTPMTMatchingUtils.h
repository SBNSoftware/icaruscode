/**
 * @file   icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h
 * @brief  Common functions for CRT/PMT matching.
 * @author Francesco Poppi ( poppi@bo.infn.it )  and Anna Heggestuen 
 * @date   April 26 2023.
 */

#include "icaruscode/IcarusObj/CRTPMTMatching.h"

namespace icarus::crt {

  struct CRTPMT {
    double tof;    ///< Time difference between CRT Hit and optical flash [ns]
    double distance;    ///< Distance between CRT Hit and optical flash centroid [cm]
    art::Ptr<sbn::crt::CRTHit> CRTHit; ///< Pointer to the matching CRT hit.
  };

  /// Contains vectors of records with matched Time of Flight and matched CRTHit.
  struct CRTMatches {
    std::vector<CRTPMT> entering; ///< Matches from outside inward.
    std::vector<CRTPMT> exiting; ///< Matches from inside outward.
    enum matchType FlashType; ///< Type of match.
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

/**
 * @brief Returns whether a flash is in time with the specified gate.
 * @param flashTime time of the flash, in trigger time scale [us]
 * @param gateType type of the gate [unused]
 * @param gateDiff time offset between PMT flashes and beam gate opening time [ns]
 * @param gateWidth the duration of the gate opening [us]
 * @return whether `flashTime` is included the gate of specified duration
 * 
 * The `gateDiff` offset represents the delay observed from the reconstructed
 * time for a flash that happens exactly at the opening of the beam gate,
 * and the opening time of the beam gate as reported by the trigger system.
 * Due to delays on the trigger decision after the digitization of the PMT
 * signals happened, this offset is typically larger than zero.
 * 
 * @note With Monte Carlo conventions, `gateDiff` needs to be `0` since the
 *       flash time is already in beam gate time scale rather than trigger's.
 */
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

/**
 * @brief Returns all the CRT hits matching the specified flash time.
 * @param flashTime the time of the flash to be matched [us]
 * @param flashpos nominal position of the flash source [cm]
 * @param crtHits list of the CRT hits to consider
 * @param interval time difference for matching flash and hit [ns]
 * @param isRealData `true` for detector data, `false` for simulation
 * @param fGlobalT0Offset CRT timing offset [ns]
 * @return a `CRTMatches` record with information about all the matched hits
 * 
 * Hits are separated between entering (before the flash) and exiting
 * (after the flash). The match is tagged according to how many entering and
 * exiting hits are found, and where.
 */
CRTMatches CRTHitmatched(
  double flashTime, geo::Point_t const& flashpos,
  std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, double interval, bool isRealData, double fGlobalT0Offset) {

  std::vector<icarus::crt::CRTPMT> enteringCRTHits;
  std::vector<icarus::crt::CRTPMT> exitingCRTHits;
  enum matchType FlashType;
  int topen = 0, topex = 0, sideen = 0, sideex = 0;
  for (auto const& crtHit : crtHits) {

    double CRTHitTime_ns = isRealData ? (crtHit->ts1_ns) : ((long long)(crtHit->ts0())-fGlobalT0Offset);
    double tof = CRTHitTime_ns - flashTime * 1e3;
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


/// Fills a `CRTPMTMatching` record out of the specified flash information.
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
  crtpmt.matchedCRTHits = thisFlash.CRTmatches;
  crtpmt.topCRTBefore = 0; //
  crtpmt.topCRTAfter = 0; //
  crtpmt.sideCRTBefore = 0; //
  crtpmt.sideCRTAfter = 0; //
  
  return crtpmt;
}

}

