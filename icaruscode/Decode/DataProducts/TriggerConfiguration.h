
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
#include <vector>

#include "icaruscode/Decode/BeamBits.h" // sbn::triggerSource

namespace icarus {

	struct TriggerConfiguration;

  /// Prints the configuration into a stream with default verbosity.
  std::ostream& operator<<(std::ostream& out, icarus::TriggerConfiguration const& config);

}


struct icarus::TriggerConfiguration {

  struct CryoConfig {

    // Majority Level for in-time activity
    unsigned int majLevelInTime = 0U;

    // Majority Level for out-of-time activity
    unsigned int majLevelDrift = 0U;

    // Window selection "Fixed" or "Overlapping"
    std::string slidingWindow = "Fixed";

  };

  struct GateConfig {

    // Return gate activation
    bool hasGate = 0U;

    // Return drift gate activation status 
    bool hasDriftGate = 0U; 

    // Return MinBias gate activation status 
    bool hasMinBiasGate = 0U;

    // Return MinBias drift gate activation status 
    bool hasMinBiasDriftGate = 0U;

    // Duration of the coincidence gate for in-time activity
    unsigned int gateWidth = 0U;

    // Duration of the coincidence gate for the out-of-time activity
    unsigned int driftGateWidth = 0U;

    // Prescale for the MinBias triggers 
    unsigned long prescaleMinBias=0U;

    // Rate of gates opened outside the extraction
    unsigned long offBeamGateRate = 1U;

    // Early warning offset for the BNB (NuMI) GatedBES ($MIBS74)
    unsigned long earlyWarningOffset = 0U; 

    // Early Early warning offset for the BNB (NuMI) $1D ($AE)
    unsigned long earlyEarlyWarningOffset = 0U; 

    // Period of two consecutive pulses from the internal pulse generator (valid for calibration gate)
    unsigned int period = 0U;

  };

  // --- BEGIN -- Data members -------------------------------------------------
 
  // NOTE when adding data members, remember to add an element to the comparison

  // Use the WR time reference
  bool useWrTime = false;

  // Add an offset between the npt and tai time as used in the wr reference (normally it is 1 or 2 leap seconds)
  unsigned int wrTimeOffset = 1'000'000'000;
 
  // Veto (delay on the leading edge of the beam gate)
  unsigned int vetoDelay = 0;

  // Cryostat configuration
  std::array<CryoConfig, icarus::kNTriggerLocation> cryoConfig; 

 	// Majority trigger type ( consider trigger from one cryostats, either cryostats, or both cryostats )
 	std::string majorityTriggerType;

 	// Run type ( MinBias: ignoring the light in-time or Majority: which applies a logic based on combination of distriminated light signals )
 	std::string runType;

 	// TPCTriggerDelay: distance between the Global trigger time and the output for the TPC. NB: It is in units of 400 ns 
 	unsigned int tpcTriggerDelay = 0;

 	// GateSelection: available gates to produce triggers: see registers 0x00000 in SBNDOCDB: 
 	std::string gateSelection ; 

  // Gate Configuration 
  std::array<GateConfig, icarus::kNTriggerSource> gateConfig;

  // --- END ---- Data members -------------------------------------------------

  // --- BEGIN -- Derived quantities -------------------------------------------

  unsigned int getBeamGateWidth( std::size_t source ){

    // Get the duration of the Gate corrected by the Veto time
    using namespace std::string_literals;
    if( gateConfig[source].hasGate ){
      return gateConfig[source].gateWidth - vetoDelay;
    } else {
      return 0U;
    }

  }


  unsigned int getDriftGateWidth( std::size_t source ){

    // Get the duration of the Drift Gate 
    if( gateConfig[source].hasDriftGate ){
      return gateConfig[source].driftGateWidth;
    } else {
      return 0U;
    }

  }


  unsigned int getOffBeamRate( std::size_t source ){

    // Get the duration of the Drift Gate 
    if( gateConfig[source].hasGate ){
      return gateConfig[source].offBeamGateRate;
    } else {
      return 0U;
    }

  }


  unsigned int getMinBiasPrescale( std::size_t source ){

    // Get the duration of the Drift Gate 
    if( gateConfig[source].hasGate ){
      return gateConfig[source].prescaleMinBias;
    } else {
      return 0U;
    }

  }

  // --- END ---- Derived quantities -------------------------------------------

 	#if __cplusplus < 202004L
  	//@{
  		/// Comparison: all fields need to have the same values.
  		bool operator== (TriggerConfiguration const& other) const;
  		bool operator!= (TriggerConfiguration const& other) const
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
  		* @brief Dumps the content of the configuration into `out` stream.
   		* @param out stream to dump the information into
   		* @param indent (default: none) indentation string
   		* @see `dump(std::ostream&, std::string const&, std::string const&, unsigned int) const`
   		* 
   		* Version of `dump()` with the specified `verbosity` level and same first
   		* indentation level as the rest.
   */

  	void dump(std::ostream& out,
    	unsigned int verbosity,
    	std::string const& indent = ""
    	) const
    { dump(out, indent, indent, verbosity); }
  
  	// -- END ---- Dump facility -------------------------------------------------
  
 
}; // sbn::ICARUSTriggerConfiguration

//------------------------------------------------------------------------------
inline bool icarus::TriggerConfiguration::operator==
  (icarus::TriggerConfiguration const& other) const
{

 if ( useWrTime              != other.useWrTime              )   return false;
 if ( wrTimeOffset           != other.wrTimeOffset           )   return false;
 if ( vetoDelay              != other.vetoDelay              )   return false;
 //if ( cryoConfig             != other.cryoConfig             )   return false;
 //if ( gateConfig             != other.gateConfig             )   return false;
 if ( majorityTriggerType    != other.majorityTriggerType    )   return false;
 if ( runType                != other.runType                )   return false;
 if ( tpcTriggerDelay        != other.tpcTriggerDelay        )   return false; 

  
 return true;

} 


//------------------------------------------------------------------------------
inline std::ostream& icarus::operator<<(std::ostream& out, icarus::TriggerConfiguration const& config)
  { config.dump(out); return out; }


//------------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DATAPRODUCT_TRIGGERCONFIGURATION_H