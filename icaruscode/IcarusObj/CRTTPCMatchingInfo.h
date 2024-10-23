/**
 * @file icaruscode/IcarusObj/CRTTPCMatchingInfo.h
 * @brief Additional information on the matching between CRT and TPC.
 * @authors Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date August 29, 2022
 * 
 */

#ifndef ICARUSCODE_ICARUSOBJ_CRTTPCMATCHINGINFO_H
#define ICARUSCODE_ICARUSOBJ_CRTTPCMATCHINGINFO_H


// C/C++ standard libraries
#include <limits>


// -----------------------------------------------------------------------------
namespace icarus { struct CRTTPCMatchingInfo; }
/**
 * @brief Additional information on the matching between CRT and TPC tracks.
 * 
 */
struct icarus::CRTTPCMatchingInfo {
  
  /// Magic value denoting the absence of DCA information.
  static constexpr double NoDistance = std::numeric_limits<double>::lowest();
  static constexpr double noTime = std::numeric_limits<double>::lowest();

  /// Magic value denoting the absence of CRTs.
  static constexpr int noCRT = std::numeric_limits<int>::lowest();
  
  // --- BEGIN -- Data members -------------------------------------------------
  
  /// Distance of closest approach between track extension and CRT hit [cm]
  double trackCrtHitDistance = NoDistance;

  /// Matched CRT Hit sub detector: 0 Top, 1 Side, 2 Bottom
  int crtSys = noCRT;

  /// Matched CRT Hit region: 30-34 Top CRT, 40-48 Side CRT, 50 Bottom CRT
  int crtRegion = noCRT;

  /// Matched CRT Hit time w.r.t. trigger [ns]
  double crtTime = noTime;

  /// Distance distinguished into its components DX, DY, DZ [cm]
  double deltaX = NoDistance;
  double deltaY = NoDistance;
  double deltaZ = NoDistance;

  /// Extrapolated Track Projection Crossing Point onto the CRT Plane
  double crossX = NoDistance;
  double crossY = NoDistance;
  double crossZ = NoDistance;

  /// Fix Coordinate for the CRT Plane: 
  /// For Top CRT region 30 Y coordinate is constant: plane=0
  /// For Top CRT region 31/32 and Side CRT 40/41/42/43/44/45 X coordinate is constant: plane=1
  /// For Top CRT region 33/34 and Side CRT 46/47 Z coordinate is constant: plane=2
  int plane=-1;  

  // --- END ---- Data members -------------------------------------------------
    
}; // icarus::CRTTPCMatchingInfo


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_ICARUSOBJ_CRTTPCMATCHINGINFO_H