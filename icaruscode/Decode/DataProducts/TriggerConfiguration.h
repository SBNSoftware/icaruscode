
/**
 * @file   icaruscode/Decode/DataProducts/TriggerConfiguration.h (todo move to sbnobj)
 * @brief  Information from the configuration of the ICARUS trigger readout.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   March 23 2021
 * @see    icaruscode/Decode/DataProducts/TriggerConfiguration.cxx
 */

#ifndef ICARUSCODE_DECODE_DATAPRODUCT_TRIGGERCONFIGURATION_H
#define ICARUSCODE_DECODE_DATAPRODUCT_TRIGGERCONFIGURATION_H

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <string>
#include <array>
#include <cassert>

#include "sbnobj/Common/Trigger/BeamBits.h"

namespace icarus {

  struct TriggerConfiguration;

  /// Prints the configuration into a stream with default verbosity.
  std::ostream& operator<<(std::ostream& out, icarus::TriggerConfiguration const& config);

}


struct icarus::TriggerConfiguration {

  struct CryoConfig {

    /// Majority Level for in-time activity
    unsigned int majLevelInTime = 0U;

    /// Majority Level for out-of-time activity
    unsigned int majLevelDrift = 0U;

    /// Window selection "Fixed" (0) or "Overlapping" (1)
    unsigned int slidingWindow = 0U;

    #if __cplusplus < 202004L
      //@{
      /// Comparison: all fields need to have the same values.
      bool operator== ( CryoConfig const& other) const noexcept;
      bool operator!= ( CryoConfig const& other) const noexcept
        { return ! this->operator== (other); }
      //@}
    #else
    # error "With C++20 support, enable the default comparison operators"
    // probably the compiler will be generating these anyway, so don't bother
    // bool operator== (CryoConfig const& other) const = default;
    // bool operator!= (CryoConfig const& other) const = default;
    #endif

  };

  struct GateConfig {

    /// Return gate activation
    bool hasGate = false;

    /// Return drift gate activation status (for out-of-time light)
    bool hasDriftGate = false; 

    /// Return MinBias triggers activation status
    bool hasMinBiasGate = false;

    /// Return MinBias drift gate activation status (for out-of-time light)
    bool hasMinBiasDriftGate = false;

    /// Duration of the gate for the in-time activity in ns 
    unsigned int gateWidth = 0U;

    /// Duration of the drift gate for the out-of-time activity in ns
    unsigned int driftGateWidth = 0U;

    /// Prescale for the MinBias triggers (calculated with respect to the number of gates opened) 
    unsigned long prescaleMinBias=1U;

    /// Rate of gates opened outside the extraction (calculated with respect to the number of gates opened) 
    unsigned long offBeamGateRate = 1U;

    /// Early warning offset for the BNB (NuMI) GatedBES ($MIBS74) in ns; used for the beam gate.
    unsigned long earlyWarningOffset = 0U; 

    /// Early Early warning offset for the BNB (NuMI) $1D ($AE) in ns; used for the drift gate.
    unsigned long earlyEarlyWarningOffset = 0U; 

    /// Period of two consecutive pulses from the internal pulse generator (valid for calibration gate) in ns
    unsigned int period = 0U;

    #if __cplusplus < 202004L
      //@{
      /// Comparison: all fields need to have the same values.
      bool operator== ( GateConfig const& other) const noexcept;
      bool operator!= ( GateConfig const& other) const noexcept
        { return ! this->operator== (other); }
      //@}
    #else
    # error "With C++20 support, enable the default comparison operators"
    // probably the compiler will be generating these anyway, so don't bother
    // bool operator== (GateConfig const& other) const = default;
    // bool operator!= (GateConfig const& other) const = default;
    #endif

  };

  // --- BEGIN -- Data members -------------------------------------------------
 
  // NOTE when adding data members, remember to add an element to the comparison

  /// Use the WR time reference
  bool useWrTime = false;

  /// Add an offset between the npt and tai time as used in the wr reference (normally it is 1 or 2 leap seconds) in ns
  unsigned int wrTimeOffset = 1'000'000'000;
 
  /// Veto (this delay has to be subtracted to the gate width ). Value is in ns 
  unsigned int vetoDelay = 0;

  /// Cryostat configuration
  std::array<CryoConfig, icarus::trigger::kNTriggerLocation> cryoConfig; 

  /// Majority trigger type (consider triggers from one cryostats, either cryostats, or both cryostats)
  std::string majorityTriggerType;

  /// Force the run to be fully a MinBias, if runType=="MinBias". If runType=="Majority" does a majority run with some prescaled minbias triggers depending on the gate selection in use
  std::string runType;

  /// TPCTriggerDelay: distance between the Global trigger time and the output for the TPC. NB: It is in units of 400 ns 
  unsigned int tpcTriggerDelay = 0;

  /// Gate Configuration 
  std::array<GateConfig, icarus::trigger::kNTriggerSource> gateConfig;

  // --- END ---- Data members -------------------------------------------------

  // --- BEGIN -- Derived quantities -------------------------------------------

  /**
   * @brief returns the effective gate width corrected for the veto delay in us 
   * @param source is the value of the sbn::bits::triggerSource enum type corresponding to the type of gate 
   *
   */
  float getGateWidth( std::size_t source ) const {

    // We really want the vetoDelay to be shorter than the gateWidth
    assert(!gateConfig[source].hasGate || (gateConfig[source].gateWidth >= vetoDelay));

    return gateConfig[source].hasGate ? 
      static_cast<float>( gateConfig[source].gateWidth - vetoDelay )/1000.  : 0.;
 
  }


  /**
   * @brief returns the width of the drift gate used for out-of-time light activity in us       
   * @param source is the value of the sbn::bits::triggerSource enum type corresponding to the type of gate 
   *
   */
  float getDriftGateWidth( std::size_t source ) const {

    return gateConfig[source].hasDriftGate ? 
      static_cast<float>( gateConfig[source].driftGateWidth )/1000. : 0U; 

  }

  /**
   * @brief returns the prescale value used to open the offbeam gates with respect to the total number of 
   * beam gates seen       
   * @param source is the value of the sbn::bits::triggerSource enum type corresponding to the type of gate 
   *
   */
  unsigned int getOffBeamRate( std::size_t source ) const {
     
    return gateConfig[source].hasGate ? gateConfig[source].offBeamGateRate : 0U;
   
  }

  /**
   * @brief returns the prescale value used to collect MinBias triggers with respect to the total number of 
   * gates seen of a particula type   
   * @param source is the value of the sbn::bits::triggerSource enum type corresponding to the type of gate 
   *
   */
  unsigned int getMinBiasPrescale( std::size_t source ) const {

    return gateConfig[source].hasGate ? gateConfig[source].prescaleMinBias : 0U;

  }


  // --- END ---- Derived quantities -------------------------------------------

  #if __cplusplus < 202004L
    //@{
    /// Comparison: all fields need to have the same values.
    bool operator== (TriggerConfiguration const& other) const noexcept;
    bool operator!= (TriggerConfiguration const& other) const noexcept
      { return ! this->operator== (other); }
    //@}
  #else
  # error "With C++20 support, enable the default comparison operators"
    // probably the compiler will be generating these anyway, so don't bother
    // bool operator== (ICARUSTriggerConfiguration const& other) const = default;
    // bool operator!= (ICARUSTriggerConfiguration const& other) const = default;
  #endif

  // -- BEGIN -- Dump facility -------------------------------------------------

  /// Maximum supported verbosity level supported by `dump()`.
  static constexpr unsigned int MaxDumpVerbosity = 2U;
  
  /// Default verbosity level for `dump()`.
  static constexpr unsigned int DefaultDumpVerbosity = MaxDumpVerbosity;
  
  
  /**
   * @brief Dumps the content of the configuration into `out` stream.
   * @param out stream to dump the information into
   * @param indent indentation string
   * @param firstIndent special indentation string for the first line
   * @param verbosity (default: `DefaultDumpVerbosity`) level of verbosity
   * 
   * The indentation string is prepended to each new line of the dump.
   * The first line indentation string is prepended before the first line of
   * the dump. The dump ends on a new empty line.
   * 
   * The amount of information printed depends on the `verbosity` level:
   * 
   * * `0`: Boardreader configuration
   * * `1`: FPGA configuration
   * * `2`: SPEXI configuration
   * 
   */

  void dump(std::ostream& out,
    std::string const& indent, std::string const& firstIndent,
    unsigned int verbosity = MaxDumpVerbosity
  ) const;
  
  /**
   * @brief Dumps the content of the configuration into `out` stream.
   * @param out stream to dump the information into
   * @param indent indentation level
   * @see `dump(std::ostream&, std::string const&, std::string const&, unsigned int) const`
   * 
   * Version of `dump()` with same first indentation level as the rest, and
   * default verbosity.
   */
   void dump(std::ostream& out, std::string const& indent = "") const
     { dump(out, indent, indent); }
  
  /**
   * @brief Dumps the content of the gate configuration into `out` stream.
   * @param out stream to dump the information into
   * @param gateConfig the gate to be dumped
   * @param indent (default: none) indentation string
   * @see `dump(std::ostream&, std::string const&, std::string const&, unsigned int) const`
   * 
   * Version of `dump()` with the specified `verbosity` level and same first
   * indentation level as the rest.
   */
  void dumpGateConfig(std::ostream& out, 
    icarus::TriggerConfiguration::GateConfig const& gateConfig, 
    std::string const& indent
  ) const;

  void dump(std::ostream& out,
    unsigned int verbosity,
    std::string const& indent = ""
    ) const
  { dump(out, indent, indent, verbosity); }
  
  // -- END ---- Dump facility -------------------------------------------------
  
  
}; // sbn::ICARUSTriggerConfiguration

//------------------------------------------------------------------------------

// Equality operators are incompete: fix for C++20 

// Equality operator for icarus::TriggerConfiguration
inline bool icarus::TriggerConfiguration::operator==
  (icarus::TriggerConfiguration const& other) const noexcept
{

 if ( useWrTime                             != other.useWrTime                             )   return false;
 if ( wrTimeOffset                          != other.wrTimeOffset                          )   return false;
 if ( vetoDelay                             != other.vetoDelay                             )   return false;
 if ( majorityTriggerType                   != other.majorityTriggerType                   )   return false;
 if ( runType                               != other.runType                               )   return false;
 if ( tpcTriggerDelay                       != other.tpcTriggerDelay                       )   return false; 

  
 return true;

} 


// Equality operator for icarus::TriggerConfiguration::CryoConfiguration
inline bool icarus::TriggerConfiguration::CryoConfig::operator==
  (icarus::TriggerConfiguration::CryoConfig const & other) const noexcept
{
    if ( majLevelInTime != other.majLevelInTime ) return false;
    if ( majLevelDrift  != other.majLevelDrift  ) return false;
    if ( slidingWindow  != other.slidingWindow  ) return false;

    return true;

}


// Equality operator for icarus::TriggerConfiguration::GateConfiguration
inline bool icarus::TriggerConfiguration::GateConfig::operator==
  (icarus::TriggerConfiguration::GateConfig const & other) const noexcept
  {
    if( hasGate                 != other.hasGate                 ) return false;
    if( hasDriftGate            != other.hasDriftGate            ) return false;
    if( hasMinBiasGate          != other.hasMinBiasGate          ) return false;
    if( gateWidth               != other.gateWidth               ) return false;
    if( driftGateWidth          != other.driftGateWidth          ) return false;
    if( prescaleMinBias         != other.prescaleMinBias         ) return false;
    if( offBeamGateRate         != other.offBeamGateRate         ) return false;
    if( earlyWarningOffset      != other.earlyWarningOffset      ) return false;
    if( earlyEarlyWarningOffset != other.earlyEarlyWarningOffset ) return false;
    if( period                  != other.period                  ) return false;

    return true;

  }


//------------------------------------------------------------------------------
inline std::ostream& icarus::operator<<(std::ostream& out, icarus::TriggerConfiguration const& config)
  { config.dump(out); return out; }


//------------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DATAPRODUCT_TRIGGERCONFIGURATION_H
