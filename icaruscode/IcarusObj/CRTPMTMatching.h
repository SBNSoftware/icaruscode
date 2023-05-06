/**
 * @file   icaruscode/IcarusObj/CRTPMTMatching.h
 * @brief  Holds the flashes matched with CRT Hits.
 * @author Francesco Poppi ( poppi@bo.infn.it )  and Anna Heggestuen 
 * @date   April 26 2023.
 */

#ifndef ICARUSCODE_ICARUSOBJ_CRTPMTMATCHING_H
#define ICARUSCODE_ICARUSOBJ_CRTPMTMATCHING_H

// C++ includes
#include <vector>
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
//#include "icaruscode/CRT/CRTPMTMatchingUtils.h"
//#include "icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h"
//#include "Math/GenVector/Cartesian3D.h"
//#include "Math/GenVector/PositionVector3D.h"
//#include "Math/GenVector/DisplacementVector3D.h"
namespace icarus::crt {

  /// Type of match between a TPC object (e.g. PMT flash) and CRT system.
  enum matchType {

    noMatch = 0,            ///< No CRT match.
    enTop = 1,              ///< Matched with Top CRT hit before optical flash.
    enSide = 2,             ///< Matched with Side CRT hit before optical flash.
    enTop_exSide = 3,       ///< Matched with one Top CRT hit before the optical flash and matched with one Side CRT hit after the optical flash.
    exTop = 4,              ///< Matched with a Top CRT after the optical flash.
    exSide = 5,             ///< Matched with a Side CRT after the optical flash.
    enTop_mult = 6,         ///< Matched with multiple Top CRT hits before the optical flash.
    enTop_exSide_mult = 7,  ///< Matched with multiple Top CRT hits before the optical flash and more then 1 side CRT hits after the optical flash.
    others = 9              ///< All the other cases.

  };
  
  /// Information about a CRT hit matched with a PMT flash.
  struct MatchedCRT {
    geo::Point_t CRTHitPos;   ///< Hit location [cm]
    double CRTPMTTimeDiff_ns; ///< CRT hit time minus PMT flash time [ns]
    double CRTTime_us;        ///< CRT hit time [us]
    int CRTSys;               ///< CRT subdetector the hit fell into.
    int CRTRegion;            ///< Region the matched CRT hit fell into.
  };

  /// Information about the match between a PMT flash and CRT hits at detector entrance and exit.
  struct CRTPMTMatching{

    int          event;         ///< Event number.
    int          run;           ///< Run number.
    unsigned int gateType;      ///< Beam gate type.
   
    int          flashID;       ///< ID of the optical flash.
    double       flashTime_us;  ///< Time of the optical flash w.r.t. the global trigger in us.
    double       flashGateTime_ns;  ///< Time of the optical flash w.r.t. the beam gate opening in ns.
    double       firstOpHitPeakTime_us;  ///< Time of the first optical hit peak time w.r.t. the global trigger in us.
    double       firstOpHitStartTime_us; ///< Time of the first optical hit start time w.r.t. the global trigger in us.
    bool         flashInGate;   ///< Flash within gate or not.
    bool         flashInBeam;   ///< Flash within the beam window of the gate or not.
    double       flashAmplitude_pe;  ///< Flash amplitude in PEs.
    geo::Point_t flashPosition; ///< Flash baricenter coordinates evaluated using ADCs as weights.
    double       flashYWidth;   ///< Flash spread along Y.
    double       flashZWidth;   ///< Flash spread along Z. 
    enum matchType flashClassification;  ///< Classification of the optical flash.
    std::vector<MatchedCRT> matchedCRTHits;    ///< Matched CRT Hits with the optical flash.
    int          topCRTBefore;  ///< Number of Top CRT Hits before the optical flash.
    int          topCRTAfter;   ///< Number of Top CRT Hits after the optical flash.
    int          sideCRTBefore; ///< Number of Side CRT Hits before the optical flash.
    int          sideCRTAfter;  ///< Number of Side CRT Hits after the optical flash.
    //std::vector<recob::OpHit>  opHits;      ///< Optical hits of the flash.
  };


}

#endif
