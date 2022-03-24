
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
#include <vector>

namespace icarus {

	struct TriggerConfiguration;

  /// Prints the configuration into a stream with default verbosity.
  std::ostream& operator<<(std::ostream& out, icarus::TriggerConfiguration const& config);

}


struct icarus::TriggerConfiguration {

	// --- BEGIN -- Data members -------------------------------------------------
  
  // NOTE when adding data members, remember to add an element to the comparison

  // Use the WR time reference
  bool UseWrTime = true;

  // Add an offset between the npt and tai time as used in the wr reference (normally it is 1 or 2 leap seconds)
  unsigned int WrTimeOffset = 1e9;
 
 	// Veto (delay on the leading edge of the beam gate)
 	unsigned int VetoDelay = std::numeric_limits<unsigned int>::max();

 	// Majority level in-time cryo 1 (EAST)
 	unsigned int MajLevelBeamCryoEAST = std::numeric_limits<unsigned int>::max();

 	// Majority level out-of-time cryo 1 (EAST)
 	unsigned int MajLevelEnableCryoEAST = std::numeric_limits<unsigned int>::max();

 	// Sliding window option cryo 1 (EAST)
 	std::string SlidingWindowCryoEAST;

 	// Majority level in-time cryo 2 (WEST)
 	unsigned int MajLevelBeamCryoWEST = std::numeric_limits<unsigned int>::max();

 	// Majority level out-of-time cryo 2 (WEST)
 	unsigned int MajLevelEnableCryoWEST = std::numeric_limits<unsigned int>::max();

 	// Sliding window option cryo 2 (WEST)
 	std::string SlidingWindowCryoWEST;

 	// Majority trigger type ( consider trigger from one cryostats, either cryostats, or both cryostats )
 	std::string MajorityTriggerType;

 	// Run type ( MinBias: ignoring the light in-time or Majority: which applies a logic based on combination of distriminated light signals )
 	std::string RunType;

 	// TPCTriggerDelay: distance between the Global trigger time and the output for the TPC. NB: It is in units of 400 ns 
 	unsigned int TPCTriggerDelay = std::numeric_limits<unsigned int>::max();

 	// GateSelection: available gates to produce triggers: see registers 0x00000 in SBNDOCDB: 
 	std::string GateSelection ; 

 	// Duration of the conincidence gate synchronous with the BNB beam
 	unsigned int BNBBeamWidth = std::numeric_limits<unsigned int>::max();

 	// Duration of the drift window gate ( for out-of-time light activity ) synchronous with the BNB beam
 	unsigned int BNBEnableWidth = std::numeric_limits<unsigned int>::max();

 	// Duration of the conincidence gate synchronous with the NuMI beam
 	unsigned int NuMIBeamWidth = std::numeric_limits<unsigned int>::max();

 	// Duration of the drift window gate ( for out-of-time light activity ) synchronous with the NuMI beam
 	unsigned int NuMIEnableWidth = std::numeric_limits<unsigned int>::max();

 	// Ratio of spills to be collected ignoring the light coincidence ( MinBias ) compared to the total number of spills sent
 	std::string PreScaleBNBNuMI; 

 	// Duration of the conincidence gate synchronous gate opened 33 ms after a BNB extraction
 	unsigned int OffBeamBNBBeamWidth = std::numeric_limits<unsigned int>::max();

	// Duration of the drift gate ( for out-of-time light activity ) synchronous gate opened 33 ms after a BNB extraction
 	unsigned int OffBeamBNBEnableWidth = std::numeric_limits<unsigned int>::max();

 	// Duration of the conincidence gate synchronous gate opened 33 ms after a NuMI extraction
 	unsigned int OffBeamNuMIBeamWidth = std::numeric_limits<unsigned int>::max();

	// Duration of the drift gate ( for out-of-time light activity ) synchronous gate opened 33 ms after a NuMI extraction
 	unsigned int OffBeamNuMIEnableWidth = std::numeric_limits<unsigned int>::max();

 	// Rate of gates opened outside the extraction
 	std::string OffBeamGateRate;

 	// Ratio of gates to be collected ignoring the light coincidence ( MinBias ) compared to the total number of offbeam gate produced
 	std::string PreScaleOffBeam;

 	// Duration of the gate opened using a periodic signal (random trigger)
 	unsigned int ZeroBiasWidth = std::numeric_limits<unsigned int>::max();

 	// Duration of the drift gate (for out-of-time activity) opened using a periodic signal (random trigger)
 	unsigned int ZeroBiasEnableWidth = std::numeric_limits<unsigned int>::max();

 	// Frequency of the random trigger in ns
 	unsigned int ZeroBiasFreq = std::numeric_limits<unsigned int>::max();

 	// Ratio of gates to be collected ignoring the light coincidence ( MinBias ) compared to the total number of gates produced
 	std::string PrescaleZeroBias;

 	// Additional offset to be added to the Gated-BES delay set in the WR network
 	unsigned int BNBBESOffset = std::numeric_limits<unsigned int>::max();

 	// Additional offset to be added to the $1D delay set in the WR network
 	unsigned int BNB1DOffset = std::numeric_limits<unsigned int>::max();

 	// Additional offset to be added to the MIBS$74 delay set in the WR network
 	unsigned int NuMIMIBSOffset = std::numeric_limits<unsigned int>::max();

 	// Additional offset to be added to the $AD delay set in the WR network
 	unsigned int NuMIADOffset = std::numeric_limits<unsigned int>::max();
  
  	// --- END ---- Data members -------------------------------------------------


  	// --- BEGIN -- Derived quantities -------------------------------------------
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

 if ( UseWrTime              != other.UseWrTime            )   return false;
 if ( WrTimeOffset           != other.WrTimeOffset         )   return false;
 if ( VetoDelay              != other.VetoDelay            )   return false;
 if ( MajLevelBeamCryoEAST   !=other.MajLevelBeamCryoEAST )    return false;
 if ( MajLevelEnableCryoEAST !=other.MajLevelEnableCryoEAST )  return false;
 if ( SlidingWindowCryoEAST  !=other.SlidingWindowCryoEAST )   return false;
 if ( MajLevelBeamCryoWEST   !=other.MajLevelBeamCryoWEST )    return false;
 if ( MajLevelEnableCryoWEST !=other.MajLevelEnableCryoWEST )  return false;
 if ( SlidingWindowCryoWEST  !=other.SlidingWindowCryoWEST )   return false;
 if ( MajorityTriggerType    !=other.MajorityTriggerType )     return false;
 if ( RunType                !=other.RunType )                 return false;
 if ( TPCTriggerDelay        != other.TPCTriggerDelay )        return false; 
 if ( GateSelection 	     != other.GateSelection )          return false; 
 if ( BNBBeamWidth 			 != other.BNBBeamWidth )           return false; 
 if ( BNBEnableWidth 		 != other.BNBEnableWidth )         return false; 
 if ( NuMIBeamWidth 		 != other.NuMIBeamWidth )          return false; 
 if ( NuMIEnableWidth 		 != other.NuMIEnableWidth )        return false; 
 if ( PreScaleBNBNuMI 		 != other.PreScaleBNBNuMI )        return false; 
 if ( OffBeamBNBBeamWidth 	 != other.OffBeamBNBBeamWidth )    return false; 
 if ( OffBeamBNBEnableWidth  != other.OffBeamBNBEnableWidth )  return false; 
 if ( OffBeamNuMIBeamWidth 	 != other.OffBeamNuMIBeamWidth )   return false; 
 if ( OffBeamNuMIEnableWidth != other.OffBeamNuMIEnableWidth ) return false; 
 if ( OffBeamGateRate 		 != other.OffBeamGateRate )        return false; 
 if ( PreScaleOffBeam 		 != other.PreScaleOffBeam )        return false; 
 if ( ZeroBiasWidth 		 != other.ZeroBiasWidth )          return false; 
 if ( ZeroBiasEnableWidth 	 != other.ZeroBiasEnableWidth )    return false; 
 if ( ZeroBiasFreq 			 != other.ZeroBiasFreq )           return false; 
 if ( PrescaleZeroBias 		 != other.PrescaleZeroBias )       return false; 
 if ( BNBBESOffset 			 != other.BNBBESOffset )           return false; 
 if ( BNB1DOffset 			 != other.BNB1DOffset )            return false; 
 if ( NuMIMIBSOffset 		 != other.NuMIMIBSOffset )         return false; 
 if ( NuMIADOffset 			 != other.NuMIADOffset )           return false; 
  
 return true;

} 


//------------------------------------------------------------------------------
inline std::ostream& icarus::operator<<(std::ostream& out, icarus::TriggerConfiguration const& config)
  { config.dump(out); return out; }


//------------------------------------------------------------------------------


#endif // ICARUSCODE_DECODE_DATAPRODUCT_TRIGGERCONFIGURATION_H