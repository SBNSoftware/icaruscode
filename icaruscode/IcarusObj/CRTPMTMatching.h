
/*
 * @file   icaruscode/IcarusObj/CRTPMTMatching.h
 * @brief  Holds the flashes matched with CRT Hits..
 * @author Francesco Poppi ( poppi@bo.infn.it )  and Anna Heggestuen 
 * @date   April 26 2023.
 */

#ifndef ICARUSCODE_ICARUSOBJ_CRTPMTMATCHING_H
#define ICARUSCODE_ICARUSOBJ_CRTPMTMATCHING_H


namespace icarus::crt {

enum matchType {

	noMatch = 0,       	///< No CRT match
	enTop = 1,         	///< Matched with Top CRT hit before optical flash.
	enSide = 2,        	///< Matched with Side CRT hit before optical flash.
 	enTop_exSide = 3,  	///< Matched with one Top CRT hit before the optical flash and matched with one Side CRT hit after the optical flash.
	exTop = 4,         	///< Matched with a Top CRT after the optical flash.
	exSide = 5,        	///< Matched with a Side CRT after the optical flash.
	enTop_mult = 6,  	///< Matched with multiple Top CRT hits before the optical flash.
	enTop_exSide_mult = 7,  ///< Matched with multiple Top CRT hits before the optical flash and more then 1 side CRT hits after the optical flash.
	others = 9		///< All the other cases.

};

struct matchedCRT{

	int		CRTHitModule;		///< Module ID of the matched CRT.
	int		CRTRegion;		///< Region identifier of the matched CRT.
	int		CRTSys;			///< CRT subsystem identifier: 0 Top CRT, 1 Side CRT, 2 Bottom CRT.
	Point_t		CRTHitPosition;		///< Coordinated of the matched CRT.
	double		CRTHitTime_us;		///< Time of the CRT Hit w.r.t. the global trigger in us.
	double		CRTHitGateTime_ns;	///< Time of the CRT Hit w.r.t. the beam gate opening in ns.
	double		CRTHitAmplitude_pe;	///< CRTHit amplitude in PEs.
	double 		CRTPMTTimeDiff_ns;	///< Time difference between the CRT Hit and the optical flash in ns.
	double		CRTHitFlashDistance;	///< Distance between the CRT Hit and the optical flash barycenter.

}


struct CRTPMTMatching{

    int				event;			///< Event number.
    int				run;			///< Run number.
    unsigned int		gateType;		///< Beam gate type.
   
    int				flashID;		///< ID of the optical flash.
    double			flashTime_us;		///< Time of the optical flash w.r.t. the global trigger in us.
    double			flashGateTime_ns;	///< Time of the optical flash w.r.t. the beam gate opening in ns.
    double			firstOpHitPeakTime_us;  ///< Time of the first optical hit peak time w.r.t. the global trigger in us.
    double			firstOpHitStartTime_us; ///< Time of the first optical hit start time w.r.t. the global trigger in us.
    bool			flashInGate;		///< Flash within gate or not.
    bool			flashInBeam;		///< Flash within the beam window of the gate or not.
    double			flashAmplitude_pe	///< Flash amplitude in PEs.
    Point_t			flashPosition;		///< Flash barycenter coordinates evaluated using ADCs as weights.
    double			flashYWidth;		///< Flash spread along Y.
    double			flashZWidth;		///< Flash spread along Z.
   
    matchType			flashClassification;	///< Classication of the optical flash.	
    std::vector<matchedCRT>	matchedCRTHits;		///< Matched CRT Hits with the optical flash.
    int				topCRTBefore;		///< Number of Top CRT Hits before the optical flash.
    int				topCRTAfter;		///< Number of Top CRT Hits after the optical flash.
    int				sideCRTBefore;		///< Number of Side CRT Hits before the optical flash.
    int				sideCRTAfter;		///< Number of Side CRT Hits after the optical flash.
    std::vector<recob::OpHit>	opHits;			///< Optical hits of the flash.
}


}


