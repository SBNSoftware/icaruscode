/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.h
 * @brief  Calibration utilities for adder board output and simulation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 5, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCALIBRATIONDATABASE_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCALIBRATIONDATABASE_H

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelID.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/Indenter.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microsecond

// framework libraries
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr
#include <ostream>
#include <sstream>
#include <string>
#include <utility> // std::move()
#include <valarray>


// -----------------------------------------------------------------------------
namespace icarus::trigger { class AdderCalibrationDatabase; }
/**
 * @brief Adder calibration database interface.
 * 
 * A calibration object returns a calibration factor for a specified
 * simulated adder waveform on a specified channel.
 * This class represents the common interface of all these calibration objects.
 * 
 * The interface consists of two stages.
 * In the first, the time period is chosen: the actual calibration object
 * is then returned, which can perform calibration on data relative to the
 * chosen run number.
 * This object is derived from the `RunCalibration` interface.
 * 
 */
class icarus::trigger::AdderCalibrationDatabase
  : protected icarus::ns::util::mfLoggingClass
{
  
    public:
  
  using microseconds = util::quantities::intervals::microseconds; // alias
  
  using RunNumber_t = int; ///< Type used to represent a run number.
  
  // --- BEGIN ---  Exceptions  ------------------------------------------------
  /// @name Exceptions
  /// @{
  
  /// Base exception for all calibration errors.
  struct Exception: cet::exception {
    Exception(std::string msg = "")
      : cet::exception{ "AdderCalibrationDatabase", std::move(msg) } {}
  };
  
  /// Requested run is not supported.
  struct RunNotSupportedError: Exception {
    
    RunNumber_t run; ///< Run the reported error is about.
    
    RunNotSupportedError(RunNumber_t run, std::string const& msg = "")
      : run{ run }
      {
        if (!msg.empty()) *this << msg;
        else *this << "This adder calibration database does not support run " << run;
      }
    
  }; // RunNotSupportedError
  
  /// @}
  // ---- END ----  Exceptions  ------------------------------------------------
  
  /**
   * @brief Interface for the calibration of adder waveforms in a specific run.
   * 
   * All queries to this object assume the data to belong to the run the object
   * was constructed for.
   * 
   * The interface is mainly the function `calibrationFactor()` to return the
   * calibration factor from the whole waveform information.
   * 
   * Note that at this point the interface only supports querying for a certain
   * channel and waveform, and it does not act on information as the waveform
   * timestamp, time relative to the beam gate, presence of clipped input and
   * such.
   * 
   */
  class RunCalibration: public icarus::ns::util::mfLoggingClass {
    
      public:
    using RunNumber_t = AdderCalibrationDatabase::RunNumber_t; // import type
    
    /// Type representing voltage, used internally for computation.
    using Voltage_t = float;
    
    /// Representation of samples in time domain for internal computation.
    using WaveformSamples_t = std::valarray<Voltage_t>;
    
    
    /// Constructor: records the requested run for future documentation.
    RunCalibration
      (RunNumber_t run, std::string logCategory = "AdderCalibrationDatabase")
      : icarus::ns::util::mfLoggingClass{ std::move(logCategory) }, fRun{ run }
      {}
    
    /// Destructor. Does nothing, but virtually so.
    virtual ~RunCalibration() = default;
    
    
    /// Returns the run number that was requested on construction (record only).
    RunNumber_t run() const { return fRun; }
    
    // --- BEGIN --- Calibration constants -------------------------------------
    /// @name Calibration constants
    /// @{
    
    /**
     * @brief Returns a single calibration factor for the specified waveform.
     * @param channel the number of adder channel the waveform is from
     * @param waveform the sequence of samples of the waveform [mV]
     * @return a calibration factor
     * @throw UnknownChannelError if `channel` is unknown or unsupported
     * 
     * A "linear" calibration factor is returned that is supposed to be
     * multiplied to all samples of the waveform.
     * 
     * The channel numbers are standardized via the `AdderChannelID` object.
     * 
     * The waveform samples are **measured in millivolt**.
     * They are specified as a `std::valarray`. It is assumed that they are the
     * result of computation and therefore are not directly from a
     * `raw::OpDetWaveform`.
     */
    double calibrationFactor
      (AdderChannelID channel, WaveformSamples_t const& waveform) const
      { return doCalibrationFactor(channel, waveform); }
    
    
    /**
     * @brief Returns a time offset to add to the waveform timestamp.
     * @param channel the number of adder channel the waveform is from
     * @param waveform the sequence of samples of the waveform [mV]
     * @return an offset as time interval
     * @throw UnknownChannelError if `channel` is unknown or unsupported
     * 
     * The returned time offset represents how later than the input PMT
     * waveforms the resulting adder waveform should be.
     * It includes delays that might be physically measured (signal propagation
     * though cables, electronics delays), but it should be considered as an
     * effective catch-all number.
     * 
     * The waveform may or may not be used to compute the delay.
     * 
     * The waveform samples are **measured in millivolt**.
     * They are specified as a `std::valarray`. It is assumed that they are the
     * result of computation and therefore are not directly from a
     * `raw::OpDetWaveform`.
     */
    microseconds timeOffset
      (AdderChannelID channel, WaveformSamples_t const& waveform) const
      { return doTimeOffset(channel, waveform); }
    
    /// @}
    // ---- END ---- Calibration constants -------------------------------------
    
    // --- BEGIN ---  Exceptions  ----------------------------------------------
    /// @name Exceptions
    /// @{
    
    /// Base exception for all calibration errors.
    struct UnknownChannelError: Exception {
      
      AdderChannelID channel; ///< The channel number the complain is about.
      
      UnknownChannelError(AdderChannelID channel, std::string const& msg = "");
      
    }; // UnknownChannelError
    
    
    /// @}
    // ---- END ----  Exceptions  ----------------------------------------------
    
      private:
    
    RunNumber_t fRun; ///< The run number that was requested at construction.
    
    // --- BEGIN --- Virtual implementation functions --------------------------
    /// @name Virtual implementation functions
    /// @{
    
    /**
     * @brief Implementation of the `calibrationFactor()` method.
     * @param channel the number of adder channel the waveform is from
     * @param waveform the sequence of samples of the waveform
     * @return a calibration factor
     * 
     * See `calibrationFactor()` for behavior details.
     */
    virtual double doCalibrationFactor
      (AdderChannelID channel, WaveformSamples_t const& waveform) const = 0;
    
    /**
     * @brief Implementation of the `timeOffset()` method.
     * @param channel the number of adder channel the waveform is from
     * @param waveform the sequence of samples of the waveform
     * @return a calibration factor
     * 
     * See `timeOffset()` for behavior details.
     * 
     * The default implementation returns no offset at all.
     */
    virtual microseconds doTimeOffset
      (AdderChannelID channel, WaveformSamples_t const& waveform) const
      { return microseconds { 0 }; }
    
    /// @}
    // ---- END ---- Virtual implementation functions --------------------------
    
  }; // RunCalibration
  
  
  /// Destructor. Does nothing, but virtually so.
  virtual ~AdderCalibrationDatabase() = default;
  
  
  /**
   * @brief Returns an independent calibration object for the specified run.
   * @param run the number of run to calibrate adder data of
   * @return a run-specific calibration object.
   * @see `RunCalibration`
   * @throw RunNotSupportedError if `run` is unknown, invalid or not supported
   * 
   * The returned object is an implementation of the `RunCalibration` interface
   * that holds all the information necessary to calibrate data from the
   * requested run.
   * 
   * The returned object is created anew and yielded to the caller for its
   * exclusive use. Repeated calls with the same run number will presumably
   * return new copies of the same object.
   */
  std::unique_ptr<RunCalibration> calibrationForRun(unsigned int run) const
    { return doCalibrationForRun(run); }
  
  
  // --- BEGIN --- Configuration dump ----------------------------------------
  /// @name Configuration dumping
  /// @{
  /**
   * @brief Dumps the configuration to a stream.
   * @tparam Stream type of stream to insert to
   * @param out the stream to dump the information into
   * @param indent string used to indent any new line
   * @param firstIndent string used to indent the starting line
   * 
   * The dump is in human-readable form and is intended for logging.
   * The output ends with a new line (and no indentation inserted).
   * 
   * The output can be multi-line, and it is indented using the `indent`
   * string, with the exception of the first line which uses `firstIndent`.
   * 
   * The `Stream` type must support insertion of stream buffers
   * (`operator<< (Stream&, std::streambuf*)`).
   */
  template <typename Stream>
  void dumpConfig
    (Stream& out, std::string indent, std::string firstIndent) const;
  
  /**
   * @brief Dumps the configuration to a stream.
   * @tparam Stream type of stream to insert to
   * @param out the stream to dump the information into
   * @param indent string used to indent all lines (including the first)
   * @see dumpConfig(Stream&, std::string, std::string)
   */
  template <typename Stream>
  void dumpConfig(Stream& out, std::string const& indent = "") const
    { dumpConfig(out, indent, indent); }
  

  /// @}
  // ---- END ---- Configuration dump ----------------------------------------
    protected:
  
  /// Constructor: initializes logging.
  AdderCalibrationDatabase(std::string logCategory = "AdderCalibrationDatabase")
    : icarus::ns::util::mfLoggingClass{ std::move(logCategory) }
    {}
  
  /// Returns the calibration factor for the specified waveform and channel.
  virtual std::unique_ptr<RunCalibration> doCalibrationForRun
    (unsigned int run) const = 0;
  
  /**
   * @brief Implementation of `dumpConfig`.
   * @param out stream to insert output to
   * @param nextLine an `Indenter` object to help with new lines
   * 
   * Implementation guidelines:
   *  * Use `out << nextLine` at the beginning of each output line
   *    (first included)
   *  * When passing the indenter around, pass it by reference.
   *  * Do not end the last line (it will be done by `dumpConfig()`).
   * 
   * See `SallenKeyFilter::doDumpConfig()` for an implementation example.
   */
  virtual void doDumpConfig
    (std::ostream& out, details::Indenter& nextline) const;
  
}; // icarus::trigger::AdderCalibrationDatabase


// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
template <typename Stream>
void icarus::trigger::AdderCalibrationDatabase::dumpConfig
  (Stream& out, std::string indent, std::string firstIndent) const
{
  details::Indenter indenter{ std::move(indent), std::move(firstIndent) };
  std::stringstream buffer; // needs to be I/O-enabled
  doDumpConfig(buffer, indenter);
  out << std::move(buffer).rdbuf() << std::endl;
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCALIBRATIONDATABASE_H
