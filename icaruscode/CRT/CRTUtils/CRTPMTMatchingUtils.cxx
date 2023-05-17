/**
 * @file   icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.cxx
 * @brief  Common functions for CRT/PMT matching.
 * @author Francesco Poppi ( poppi@bo.infn.it )  and Anna Heggestuen 
 * @date   April 26 2023.
 * @see    icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h
 */

// Library header
#include "icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <utility>


// -----------------------------------------------------------------------------
bool icarus::crt::flashInTime(double flashTime, int gateType, double gateDiff,
                 double gateWidth) {
  // As a reminder, I will leave here the commented part of the vetoOffset, in
  // case something changes in the future
  /*double vetoOffset = 3.5;*/
  double activeGate = gateWidth /*- vetoOffset*/;

  double relFlashTime = flashTime + gateDiff / 1000. /*- vetoOffset*/;
  mf::LogDebug("CRTPMTMatchingProducer_FlashInTime")
    << "Gate Diff " << gateDiff / 1000 << " Ftime+gateDiff "
    << flashTime + gateDiff / 1000. << " " << activeGate;

  return ((relFlashTime > 0) && (relFlashTime < activeGate));
}


// -----------------------------------------------------------------------------
double icarus::crt::CRTHitTime
  (sbn::crt::CRTHit const& hit, double globalT0Offset, bool isRealData)
{
  return isRealData? hit.ts1_ns
    : static_cast<std::int64_t>(hit.ts0()) - static_cast<std::int64_t>(globalT0Offset);
}


// -----------------------------------------------------------------------------
icarus::crt::CRTMatches icarus::crt::CRTHitmatched(
  double flashTime, geo::Point_t const& flashpos,
  std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, double interval, bool isRealData, double globalT0Offset) {

  std::vector<icarus::crt::CRTPMT> enteringCRTHits;
  std::vector<icarus::crt::CRTPMT> exitingCRTHits;
  MatchType flashType;
  int topen = 0, topex = 0, sideen = 0, sideex = 0, bottomen = 0, bottomex = 0;
  for (auto const& crtHit : crtHits) {
    // care with conversions: if either side of a subtraction is a `double`,
    // the other side is also converted to `double` just before the subtraction,
    // and in this conversion it may lose precision; subtraction itself must
    // operate between 64-bit integers, then the conversion of the result may happen
    double const CRTHitTime_ns = CRTHitTime(*crtHit, globalT0Offset, isRealData);
    double tof = CRTHitTime_ns - flashTime * 1e3;
    double distance =
      (flashpos - geo::Point_t{crtHit->x_pos, crtHit->y_pos, crtHit->z_pos})
      .R();
    if (abs(tof) >= interval) continue;
    if (tof < 0) {
      if (crtHit->plane > 49){
	std::cout << "found a bottom match, plane = " << crtHit->plane << ", crthit time = " << CRTHitTime_ns << " ns. \n";
	bottomen++;
      }
      else if (crtHit->plane > 36)
        sideen++;
      else
        topen++;
      CRTPMT m_match = {tof, distance, crtHit};
      enteringCRTHits.push_back(std::move(m_match));
    } else if (tof >= 0) {
      if (crtHit->plane > 49){
	std::cout << "found a bottom match, plane = " << crtHit->plane << ", crthit time = " << CRTHitTime_ns << " ns. \n";
	bottomex++;
      }
      if (crtHit->plane > 36)
        sideex++;
      else
        topex++;
      CRTPMT m_match = {tof, distance, crtHit};
      exitingCRTHits.push_back(std::move(m_match));
    }
  }
  if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::noMatch;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::enTop;
  else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::enSide;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::enTop_exSide;
  else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::exTop;
  else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::exSide;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::enTop_mult;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1 && bottomen == 0 && bottomex == 0)
    flashType = MatchType::enTop_exSide_mult;
  else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0 && bottomen == 1 && bottomex == 0)
    flashType = MatchType::enBottom;
  else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0 && bottomen == 0 && bottomex == 1)
    flashType = MatchType::exBottom;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0 && bottomen == 0 && bottomex == 1)
    flashType = MatchType::enTop_exBottom;
  else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0 && bottomen == 0 && bottomex == 1)
    flashType = MatchType::enSide_exBottom;
  else
    flashType = MatchType::others;

  return {std::move(enteringCRTHits), std::move(exitingCRTHits), flashType};
}


// -----------------------------------------------------------------------------
icarus::crt::CRTPMTMatching icarus::crt::FillCRTPMT(FlashType const& flash)
{
  // UPDATE CRTPMTMatchingProducer class documentation when some missing items
  //        are added
  CRTPMTMatching crtpmt;
  crtpmt.flashID = 0; //
  crtpmt.flashTime = flash.flashTime;
  crtpmt.flashGateTime = flash.flashGateTime;
  crtpmt.firstOpHitPeakTime = flash.firstOpHitPeakTime; 
  crtpmt.firstOpHitStartTime = flash.firstOpHitStartTime;
  crtpmt.flashInGate = flash.inGate;
  crtpmt.flashInBeam = flash.inBeam;
  crtpmt.flashPE = flash.flashPE; 
  crtpmt.flashPosition = flash.flashPos;
  crtpmt.flashYWidth = flash.flashYWidth;
  crtpmt.flashZWidth = flash.flashZWidth;
  crtpmt.flashClassification = flash.classification;
  crtpmt.matchedCRTHits = flash.CRTmatches;
  // crtpmt.nTopCRTHitsBefore = 0; //
  // crtpmt.nTopCRTHitsAfter = 0; //
  // crtpmt.nSideCRTHitsBefore = 0; //
  // crtpmt.nSideCRTHitsAfter = 0; //
  
  return crtpmt;
}


// -----------------------------------------------------------------------------
