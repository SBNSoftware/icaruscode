/**
 * @file   icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h
 * @brief  Common functions for CRT/PMT matching.
 * @author Francesco Poppi ( poppi@bo.infn.it )  and Anna Heggestuen (aheggest@colostate.edu)
 * @date   April 26 2023.
 * @see    icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.cxx
 */

#ifndef ICARUSCODE_CRT_CRTUTILS_CRTPMTMATCHINGUTILS_H
#define ICARUSCODE_CRT_CRTUTILS_CRTPMTMATCHINGUTILS_H


//#include "icaruscode/IcarusObj/CRTPMTMatching.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "canvas/Persistency/Common/Ptr.h" 

#include <vector>

namespace recob { class OpFlash; } // no need to know the details

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
    sbn::crt::MatchType flashType; ///< Type of match.
  };



  /**
   * @brief Returns whether a flash is in time with the specified gate.
   * @param flashTime time of the flash, in trigger time scale [us]
   *  " fill in for new vars "
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
                   double gateWidth);

  /**
   * @brief Returns the relative time of a CRT hit [ns]
   * @param hit the CRT hit
   * @param globalT0Offset for simulation, the global T0 offset in nanoseconds
   * @param isRealData whether the hit is from data (as opposed to simulation)
   * @return the time of the hit [ns]
   * 
   * For data hits, the time is directly taken from the timestamp 1 (`ts1_ns`).
   * For simulation, the time is taken from the "absolute" timestamp 0 (`ts0()`)
   * after subtracting a known offset (`globalT0Offset`).
   */
  double CRTHitTime
    (sbn::crt::CRTHit const& hit, double globalT0Offset, bool isRealData);

  /**
   * @brief Returns all the CRT hits matching the specified flash time.
   * @param flashTime the time of the flash to be matched [us]
   * @param flashpos nominal position of the flash source [cm]
   * @param crtHits list of the CRT hits to consider
   * @param interval time difference for matching flash and hit [ns]
   * @param isRealData `true` for detector data, `false` for simulation
   * @param globalT0Offset CRT timing offset [ns]
   * @return a `CRTMatches` record with information about all the matched hits
   * 
   * Hits are separated between entering (before the flash) and exiting
   * (after the flash). The match is tagged according to how many entering and
   * exiting hits are found, and where.
   */
  CRTMatches CRTHitmatched(
    double flashTime, geo::Point_t const& flashpos,
    std::vector<art::Ptr<sbn::crt::CRTHit>>& crtHits, double interval, bool isRealData, double globalT0Offset,  bool matchBottomCRT);


  //@{
  /**
   * @brief Returns an unique ID for the flash in `ptr`.
   * @param flash a flash or its pointer
   * @param version (default: latest) the ID creation algorithm version to use
   * @return an unique ID for the flash in `ptr`
   * 
   * An ID is made up out of the flash.
   * 
   * Supported versions:
   *  * `10`: ID is the flash time in picoseconds, rounded;
   *          `icarus::crt::CRTPMTMatching::NoID` if flash is invalid
   *          (e.g. null pointer)
   */
  int makeFlashID
    (recob::OpFlash const& flash, int version = std::numeric_limits<int>::max());

  int makeFlashID
    (recob::OpFlash const* flash, int version = std::numeric_limits<int>::max());
  //@}
} // end icarus::crt namespace 

namespace sbn::crt { 
  struct FlashType {
    //static constexpr auto NoID = icarus::crt::CRTPMTMatching::NoID;
    static constexpr int NoID = std::numeric_limits<int>::min();
    
    int flashID = NoID;
    geo::Point_t flashPos;
    double flashTime;
    double flashGateTime;
    double firstOpHitPeakTime;
    double firstOpHitStartTime;
    //bool matchBottomCRT;
    bool inBeam;
    bool inGate;
    double flashPE;
    double flashYWidth;
    double flashZWidth;
    MatchType classification;
    std::vector<MatchedCRT> CRTmatches;
  };   
    /// Fills a `CRTPMTMatching` record out of the specified `flash` information.
  CRTPMTMatching FillCRTPMT (FlashType const& flash);

} // end sbn::crt namespace 





#endif // ICARUSCODE_CRT_CRTUTILS_CRTPMTMATCHINGUTILS_H
