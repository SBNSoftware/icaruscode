/**
 * @file   icaruscode/PMT/Trigger/Algorithms/FixedTriggerGateBuilder.h
 * @brief  Fixed-length gate builder.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Algorithms/FixedTriggerGateBuilder.cxx`
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FIXEDTRIGGERGATEBUILDER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FIXEDTRIGGERGATEBUILDER_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/ManagedTriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/Utilities/quantities/electronics.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>


namespace icarus::trigger {
  
  // ---------------------------------------------------------------------------
  //
  // declarations
  //
  
  class FixedTriggerGateBuilder;
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Fixed-length gate builder.
 * @see `icarus::trigger::DynamicTriggerGateBuilder`
 * 
 * Gates are opened when the waveform goes above threshold, and closed a fixed
 * time after that.
 * 
 * Note that gates at higher thresholds are usually not contained in the ones
 * with lower thresholds: the former ones start later _and_ end later.
 * 
 */
class icarus::trigger::FixedTriggerGateBuilder
  : public icarus::trigger::ManagedTriggerGateBuilder
{
  using Base_t = icarus::trigger::ManagedTriggerGateBuilder;
  
  class FixedGateManager: private GateManager {
    
    struct FixedGateInfo: public GateInfoBase {
      /// Activity ignored up to this tick excluded.
      optical_tick openUntil = TriggerGate_t::MinTick;
      
      ///Ticks to keep the gate open.
      optical_time_ticks gateDuration;
      
      FixedGateInfo(
        TriggerGate_t& gate, optical_time_ticks gateDuration,
        optical_tick openUntil = TriggerGate_t::MinTick
        )
        : GateInfoBase(gate), gateDuration(gateDuration) {}
      
      void belowThresholdAt(optical_tick /* tick */) {}
      void aboveThresholdAt(optical_tick tick)
        {
          using namespace util::quantities::electronics_literals;
          MF_LOG_TRACE(TriggerGateDebugLog)
            << "Declared above threshold at: " << tick;
          if (tick < openUntil) {
            MF_LOG_TRACE(TriggerGateDebugLog)
              << "  we are in dead time until " << openUntil
              << ", come back later.";
            return; // open only if not in "dead time"
          }
          gate().openFor(tick, gateDuration); // open the gate for this long
          openUntil = tick + gateDuration; // set some dead time
          MF_LOG_TRACE(TriggerGateDebugLog) << "  gate (re)opened ("
            << gate().openingCount(tick - 1_tick)
            << " => " << gate().openingCount(tick)
            << ") for " << gateDuration << " until " << openUntil
            << " (" << gate().openingCount(openUntil - 1_tick)
            << " => " << gate().openingCount(openUntil) << ")";
        }
      
    }; // struct FixedGateInfo
    
    
    optical_time_ticks gateDuration;
    
    
      public:
    using GateInfo_t = FixedGateInfo;
    
    FixedGateManager(optical_time_ticks gateDuration)
      : gateDuration(gateDuration) {}
    
    GateInfo_t create(GateInfo_t::TriggerGate_t& gate) const
      { return { gate, gateDuration }; }
    
  }; // struct FixedGateManager
  

    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config: public Base_t::Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
    fhicl::Atom<microsecond::value_t> GateDuration {
      Name("GateDuration"),
      Comment("duration of a trigger gate [us]")
      };
    
  }; // struct Config
  // --- END Configuration -----------------------------------------------------
  
  
  /// Constructor: sets the configuration.
  FixedTriggerGateBuilder(Config const& config);
  
  
  /// Algorithm setup.
  virtual void setup(detinfo::DetectorTimings const& timings) override;
  
  /// Returns a collection of `TriggerGates` objects sorted by threshold.
  virtual std::vector<TriggerGates> build
    (std::vector<raw::OpDetWaveform> const& waveforms) const override
    { return unifiedBuild(FixedGateManager(fGateTicks), waveforms); }
  
    private:
  
  // --- BEGIN Configuration parameters ----------------------------------------
  
  microsecond fGateDuration; ///< Duration of a channel gate [&micro;s]
  
  optical_time_ticks fGateTicks; ///< Duration of a channel gate [ticks]
  
  // --- END Configuration parameters ------------------------------------------
  
  
  // --- BEGIN Setup -----------------------------------------------------------
  
  // --- END Setup -------------------------------------------------------------
  
  
}; // class icarus::trigger::FixedTriggerGateBuilder


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FIXEDTRIGGERGATEBUILDER_H
