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
  
  
  // --- BEGIN -- Data members -------------------------------------------------
  
  /// Distance of closest approach between track extension and CRT hit [cm]
  double DCA = NoDistance;
  
  /// How far the track extension meets the relevant CRT plane [cm]
  double extrapLength = NoDistance;
  
  // --- END ---- Data members -------------------------------------------------
  
  
  
}; // icarus::CRTTPCMatchingInfo


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_ICARUSOBJ_CRTTPCMATCHINGINFO_H