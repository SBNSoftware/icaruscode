#ifndef SBN_TRACKTREESTOREOBJ_H
#define SBN_TRACKTREESTOREOBJ_H

namespace sbn {

  struct selTrackInfo {
    int trackID;
    float t0;
    float start_x;
    float start_y;
    float start_z;
    float end_x;
    float end_y;
    float end_z;
    float length;
    float dir_x;
    float dir_y;
    float dir_z;
  selTrackInfo():
    trackID(-1),
      t0(-1),
      start_x(-1),
      start_y(-1),
      start_z(-1),
      end_x(-1),
      end_y(-1),
      end_z(-1),
      length(-1),
      dir_x(-1),
      dir_y(-1),
      dir_z(-1)
    {}
  };

  struct selBeamInfo {
    uint64_t beamGateSimStart;
    float beamGateDuration;
    unsigned int beamGateType;
  selBeamInfo():
    beamGateSimStart(0),
      beamGateDuration(-1),
    beamGateType(999)
    {}
  };
  
  struct selTriggerInfo {
    unsigned int beamType;
    uint64_t triggerTime;
    uint64_t beamGateTime;
    uint64_t triggerID;
    uint64_t gateID;
  selTriggerInfo():
    beamType(0),
    triggerTime(0),
      beamGateTime(0),
      triggerID(0),
      gateID(0)
    {}
  };
  
} // namespace sbn

#endif // SBN_TRACKTREESTOREOBJ_H
