/**
 * @file   icaruscode/NuGraph/StitchingUtils.h
 * @brief  Utilities to stitch together the different logical TPCs in each cryostat.
 * @author Giuseppe Cerati (cerati@fnal.gov)
 * @date   June 05, 2025
 *
 * This library provides utilities to stitch together the different logical TPCs in each cryostat.
 * As a result, a single time vs wire coordinate system can be used in each cryostat.
 *
 * This is a header-only library.
 */

#ifndef ICARUSCODE_NUGRAPH_STITCHINGUTILS_H
#define ICARUSCODE_NUGRAPH_STITCHINGUTILS_H

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace util {

  int stitchedPlane(const geo::WireID wid) {
    int plane = wid.Plane;
    if(wid.TPC==2 || wid.TPC==3) {
      if(wid.Plane==1) plane=2;
      else if(wid.Plane==2) plane=1;
    }
    return plane;
  }

  float stitchedTime(const geo::WireID wid, float timein) {
    float time = timein;
    if(wid.TPC==2 || wid.TPC==3) {
      //correction = 2*(tpcgeo.DriftDistance()/detProp.DriftVelocity()-clockData.TriggerOffsetTPC())/clockData.TPCClock().TickPeriod() = 6442.15 us
      time = 6442.15 - timein;
    }
    return time;
  }

  size_t stitchedWire(const geo::WireID wid) {
    size_t wire = wid.Wire;
    int plane = stitchedPlane(wid);
    if(wid.TPC==1 || wid.TPC==3) {
      if(plane==1 || plane == 2) {
	wire = wid.Wire + 2535; //2535 is the last part of the wires in cryos 0 an 2 before the cut in z=0
      }
    }
    return wire;
  }

} // namespace util
#endif
