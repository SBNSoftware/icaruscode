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
    
    int trackID = -1; ///< `recob::Track::ID()` from the reconstructed track.
    double t0 = NoTime; ///< Best time of the track (used for selection) [us]
    double t0_TPC = NoTime; ///< Time of (cathode crossing) track from TPC reconstruction [us]
    double t0_CRT = NoPosition; ///< Time of track from association with CRT hit [us]
    float t0_TPC_min = NoTime; ///< Minimum track time allowed by its hits [us]
    float t0_TPC_max = NoTime; ///< Maximum track time allowed by its hits [us]
    float t0_diff = NoTime; ///< Distance of `t0` time from the one allowed by hits [us]
    float t0_CRT_diff = NoTime; ///< Distance of CRT time from the one allowed by hits [us]
    float hitTick_min = NoTime; ///< Lowest tick time among hits of this track.
    float hitTick_max = NoTime; ///< Highest tick time among hits of this track.
    float start_x = NoPosition; ///< Start of the track (assumed downgoing) [cm]
    float start_y = NoPosition; ///< Start of the track (assumed downgoing) [cm]
    float start_z = NoPosition; ///< Start of the track (assumed downgoing) [cm]
    float end_x = NoPosition; ///< End of the track (assumed downgoing) [cm]
    float end_y = NoPosition; ///< End of the track (assumed downgoing) [cm]
    float end_z = NoPosition; ///< End of the track (assumed downgoing) [cm]
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
    float length = -1.0; ///< Full track length [cm]
    float energy = NoEnergy; ///< Total energy as reported by calorimetry [MeV]
    float energy_int = NoEnergy; ///< Integral of all (dE/dx) x dx [MeV]
    float charge_int = NoEnergy; ///< Integral of all (dQ/dx) x dx [electrons?]
    float energy_range = NoEnergy; ///< Energy from range (assuming stopping) [MeV]
    float dir_x = NoPosition; ///< Direction vector at the start if the track.
    float dir_y = NoPosition; ///< Direction vector at the start if the track.
    float dir_z = NoPosition; ///< Direction vector at the start if the track.
    float driftCorrX = 0.0; ///< Shift applied to track position for time [cm]
    bool  flipped = false; ///< Whether the track has been flipped backward.
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
    float flash_time = NoPosition; ///< Nominal flash time [us]
    float flash_x = NoPosition;
    float flash_y = NoPosition;
    float flash_z = NoPosition;
    float diff_flash_pos = NoPosition; ///< Distance from track middle point [cm]
    float diff_flash_t0 = NoPosition; ///< Time from track time `t0` [us]
    float diff_flash_TPCt0 = NoPosition; ///< Time from TPC time `t0_TPC` [us]
    float diff_flash_CRTt0 = NoPosition; ///< Time from CRT time `t0_CRT` [us]
    bool flash_closest_to_track = false; ///< Is this the smallest `diff_flash_t0`?
    bool flash_nearest_to_track = false; ///< Is this the smallest `diff_flash_pos`?
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
