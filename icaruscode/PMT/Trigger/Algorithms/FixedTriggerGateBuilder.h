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
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/Utilities/quantities_fhicl.h" // for microsecond parameter
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
 * Gates are opened when the waveform transits above threshold, and closed a
 * fixed time after that.
 * Optionally, if a transition happens when the gate is already open, the
 * opening can be extended by another fixed time; otherwise, the opening time
 * is considered dead time and transitions occurring in that while will be
 * ignored.
 * 
 * Note that gates at higher thresholds are usually not contained in the ones
 * with lower thresholds: the former ones start later _and_ end later.
 * Also, transitions on higher thresholds do not affect the lower threshold
 * gates. For example, if the thresholds are at 50 and 100 counts, and the
 * waveforms raises to 110, then starts oscillating between 90 and 110 counts,
 * the gate for threshold 100 will be reopened or kept open, but the one for
 * threshold 50 will close after the fixed time has elapsed and stay closed
 * as long as no transition from under 50 to 50 or more happens.
 * 
 */
class icarus::trigger::FixedTriggerGateBuilder
  : public icarus::trigger::ManagedTriggerGateBuilder
{
  using Base_t = icarus::trigger::ManagedTriggerGateBuilder;
  
  class FixedGateManager: private GateManager {
    
    struct FixedGateInfo: public GateInfoBase {
      /// Activity ignored up to this tick excluded.
      optical_tick openUntil { TriggerGateData_t::MinTick };
      
      /// Ticks to keep the gate open.
      optical_time_ticks gateDuration;
      
      /// Whether new crossings extend the duration of opening or are ignored.
      bool const extendGate { false };
      
      FixedGateInfo(
        TriggerGate_t& gate, optical_time_ticks gateDuration,
        bool extendGate = false
        )
        : GateInfoBase(gate), gateDuration(gateDuration), extendGate(extendGate)
        {}
      
      void belowThresholdAt(optical_tick /* tick */) {}
      void aboveThresholdAt(optical_tick tick)
        {
          using namespace util::quantities::electronics_literals;
          MF_LOG_TRACE(details::TriggerGateDebugLog)
            << "Declared above threshold at: " << tick;
          if ((tick < openUntil) && !extendGate) {
            MF_LOG_TRACE(details::TriggerGateDebugLog)
              << "  we are in dead time until " << openUntil
              << ", come back later.";
            return; // open only if not in "dead time"
          }
          
          // open or extend the gate:
          auto const reopenAt = extendGate? openUntil: tick;
          openUntil = tick + gateDuration; // set some dead time
          gate().openBetween(reopenAt.value(), openUntil.value());
          
          MF_LOG_TRACE(details::TriggerGateDebugLog) << "  gate (re)opened ("
            << gate().openingCount((tick - 1_tick).value())
            << " => " << gate().openingCount(tick.value())
            << ") for " << gateDuration << " until " << openUntil
            << " (" << gate().openingCount((openUntil - 1_tick).value())
            << " => " << gate().openingCount(openUntil.value()) << ")";
        } // aboveThresholdAt()
      
    }; // struct FixedGateInfo
    
    
    optical_time_ticks gateDuration;
    
    bool const extendGate;
    
    
      public:
    using GateInfo_t = FixedGateInfo;
    
    FixedGateManager(optical_time_ticks gateDuration, bool extendGate = false)
      : gateDuration(gateDuration), extendGate(extendGate) {}
    
    GateInfo_t create(GateInfo_t::TriggerGate_t& gate) const
      { return { gate, gateDuration, extendGate }; }
    
  }; // struct FixedGateManager
  

    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config: public Base_t::Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
    fhicl::Atom<microsecond> GateDuration {
      Name("GateDuration"),
      Comment("duration of a trigger gate")
      };
    
    fhicl::Atom<bool> ExtendGate {
      Name("ExtendGate"),
      Comment("whether gate is extended by crossing activity inside it"),
      false
      };
    
  }; // struct Config
  // --- END Configuration -----------------------------------------------------
  
  
  /// Constructor: sets the configuration.
  FixedTriggerGateBuilder(Config const& config);
  
  
  /// Algorithm setup.
  virtual void setup(detinfo::DetectorTimings const& timings) override;
  
  /// Returns a collection of `TriggerGates` objects sorted by threshold.
  virtual std::vector<TriggerGates> build
    (std::vector<WaveformWithBaseline> const& waveforms) const override
    {
      return unifiedBuild(FixedGateManager(fGateTicks, fExtendGate), waveforms);
    }
  
    private:
  
  // --- BEGIN Configuration parameters ----------------------------------------
  
  microsecond fGateDuration; ///< Duration of a channel gate [&micro;s]
  
  optical_time_ticks fGateTicks; ///< Duration of a channel gate [ticks]
  
  bool const fExtendGate; ///< Whether gate opening time can be extended.
  
  // --- END Configuration parameters ------------------------------------------
  
  
  // --- BEGIN Setup -----------------------------------------------------------
  
  // --- END Setup -------------------------------------------------------------
  
  
  /// Prints the class configuration.
  /// 
  /// Assumes the start of a new line and does not break the last one.
  virtual void doDumpConfiguration(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const override;
  
  /// Prints the algorithm-specific configuration.
  ///
  /// Assumes the start of a new line and does not break the last one.
  void dumpLocalConfiguration(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const;

  
}; // class icarus::trigger::FixedTriggerGateBuilder


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FIXEDTRIGGERGATEBUILDER_H
