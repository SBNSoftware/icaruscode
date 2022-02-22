#ifndef SBN_TRACKTREESTOREOBJ_H
#define SBN_TRACKTREESTOREOBJ_H

#include "cstdint" // std::uint64_t
#include <limits>

namespace sbn {

  struct selTrackInfo {
    static constexpr float NoPosition = -999999.0;
    
    int trackID = -1;
    float t0 = NoPosition;
    float start_x = NoPosition;
    float start_y = NoPosition;
    float start_z = NoPosition;
    float end_x = NoPosition;
    float end_y = NoPosition;
    float end_z = NoPosition;
    float length = -1.0;
    float energy = NoPosition;
    float charge_int = NoPosition;
    float charge_dqdx = NoPosition;
    float dir_x = NoPosition;
    float dir_y = NoPosition;
    float dir_z = NoPosition;
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
    unsigned int triggerID = 0;
    unsigned int gateID = 0;
  };
  
  struct selSimTriggerInfo {
    
    static constexpr double NoTime = -999999.0;
    
    double time = NoTime; ///< Time of the trigger in electronics time scale.
    
  };

  struct selLightInfo {
    static constexpr float NoPosition = -999999.0;
    int flash_id = -1;
    float sum_pe = NoPosition;
    float flash_time = NoPosition;
    float flash_x = NoPosition;
    float flash_y = NoPosition;
    float flash_z = NoPosition;
    float diff_flash_t0 = NoPosition;
  };

  struct selHitInfo {
    static constexpr float NoValue = -999999.0;
    static constexpr uint16_t NoVal = 65535;
    float integral = NoValue;
    float sumadc = NoValue;
    float width = NoValue; ///< RMS
    float pk_time = NoVal;
    uint16_t mult = NoVal;
    uint16_t wire = NoVal;
    uint16_t plane = NoVal;
    int channel = NoVal;
    uint16_t tpc = NoVal;
    int16_t end = -32767;
    int16_t start = -32767;
    int id = NoValue;
    bool oncalo = false;
    float pitch = NoValue;
    float dqdx = NoValue;
    float dEdx = NoValue;
    float rr = NoValue;
    
  };
  
} // namespace sbn

#endif // SBN_TRACKTREESTOREOBJ_H
