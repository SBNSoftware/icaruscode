#ifndef SBN_TRACKTREESTOREOBJ_H
#define SBN_TRACKTREESTOREOBJ_H

#include "cstdint" // std::uint64_t

namespace sbn {

  struct selTrackInfo {
    int trackID = -1;
    float t0 = -1.0;
    float start_x = -1.0;
    float start_y = -1.0;
    float start_z = -1.0;
    float end_x = -1.0;
    float end_y = -1.0;
    float end_z = -1.0;
    float length = -1.0;
    float dir_x = -1.0;
    float dir_y = -1.0;
    float dir_z = -1.0;
  };  // selTrackInfo

  struct selBeamInfo {
    std::uint64_t beamGateSimStart = 0;
    float beamGateDuration = -1.0;
    unsigned int beamGateType = 999;
  };
  
  struct selTriggerInfo {
    unsigned int beamType = 0;
    std::uint64_t triggerTime = 0;
    std::uint64_t beamGateTime = 0;
    std::uint64_t triggerID = 0;
    std::uint64_t gateID = 0;
  };
  
} // namespace sbn

#endif // SBN_TRACKTREESTOREOBJ_H
