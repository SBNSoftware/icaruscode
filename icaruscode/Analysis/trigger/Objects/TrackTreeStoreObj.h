#ifndef SBN_TRACKTREESTOREOBJ_H
#define SBN_TRACKTREESTOREOBJ_H

#include "cstdint" // std::uint64_t
#include <limits>

namespace sbn {

  /**
   * 
   * Middle point is cutting the track in two chunks with the same length.
   * 
   * The part "before" the cathode is always the one at lower _x_ coordinates
   * than the cathode itself.
   * 
   */
  struct selTrackInfo {
    static constexpr float NoPosition = -999999.0;
    static constexpr float NoEnergy = -999999.0;
    static constexpr double NoTime = -99999.0;
    
    int trackID = -1;
    double t0 = NoTime;
    double t0_TPC = NoTime;
    double t0_CRT = NoPosition;
    float start_x = NoPosition;
    float start_y = NoPosition;
    float start_z = NoPosition;
    float end_x = NoPosition;
    float end_y = NoPosition;
    float end_z = NoPosition;
    float middle_x = NoPosition; ///< Half-way distance along the trajectory.
    float middle_y = NoPosition; ///< Half-way distance along the trajectory.
    float middle_z = NoPosition; ///< Half-way distance along the trajectory.
    float atcathode_x = NoPosition; ///< Cathode crossing point of trajectory.
    float atcathode_y = NoPosition; ///< Cathode crossing point of trajectory.
    float atcathode_z = NoPosition; ///< Cathode crossing point of trajectory.
    float midbeforecathode_x = NoPosition; ///< Midpoint of subpath before cathode.
    float midbeforecathode_y = NoPosition; ///< Midpoint of subpath before cathode.
    float midbeforecathode_z = NoPosition; ///< Midpoint of subpath before cathode.
    float midaftercathode_x = NoPosition; ///< Midpoint of subpath after cathode.
    float midaftercathode_y = NoPosition; ///< Midpoint of subpath after cathode.
    float midaftercathode_z = NoPosition; ///< Midpoint of subpath after cathode.
    float beforecathode = -1; ///< Track path length before cathode (lower _x_).
    float aftercathode = -1; ///< Track path length before cathode (higher _x_).
    float length = -1.0;
    float energy = NoEnergy;
    float energy_int = NoEnergy;
    float charge_int = NoEnergy;
    float energy_range = NoEnergy;
    // float charge_dqdx = NoPosition;
    float dir_x = NoPosition;
    float dir_y = NoPosition;
    float dir_z = NoPosition;
  };  // selTrackInfo

  struct selBeamInfo {
    static constexpr float NoTime = -99999.0;
    
    float beamGateSimStart = NoTime;
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
    
    static constexpr double NoTime = -99999.0;
    
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
    float px = NoValue;
    float py = NoValue;
    float pz = NoValue;
    float dirx = NoValue;
    float diry = NoValue;
    float dirz = NoValue;
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
  
  struct selCRTpoint {
    static constexpr double Nowhere = -999999.0;
    static constexpr double NoTime = -99999.0;
    static constexpr int NoRegion = -1;
    double time = NoTime; ///< Hit time from trigger [us]
    double x = Nowhere;
    double y = Nowhere;
    double z = Nowhere;
    float resolution = -1.0;
    float DCA = Nowhere;
    int region = NoRegion;
  };
  
  struct selCRTInfo {
    static constexpr double NoTime = -99999.0;
    selCRTpoint entry;
    selCRTpoint exit;
    double time = NoTime; ///< Matching time from trigger [us]
  };
  
} // namespace sbn

#endif // SBN_TRACKTREESTOREOBJ_H
