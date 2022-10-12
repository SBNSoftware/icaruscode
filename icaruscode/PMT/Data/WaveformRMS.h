/**
 * @file   icaruscode/PMT/Data/WaveformRMS.h
 * @brief  A baseline RMS for a waveform.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2022
 * @see    icaruscode/PMT/Data/WaveformRMS.cxx
 */

#ifndef ICARUSCODE_PMT_DATA_WAVEFORMRMS_H
#define ICARUSCODE_PMT_DATA_WAVEFORMRMS_H


// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <cmath> // std::round()


//------------------------------------------------------------------------------
namespace icarus {
  
  struct WaveformRMS;
  
  /// Prints the value of the RMS into a stream.
  std::ostream& operator<<
    (std::ostream& out, icarus::WaveformRMS const& baseline);
  
} // namespace icarus

/**
 * @brief Class containing a waveform baseline RMS value.
 * @see `icarus::WaveformBaseline`
 * 
 * The baseline RMS is stored as a floating point value.
 * 
 * This class is a data product wrapper for a simple value, with some usability
 * candies attached.
 * 
 * 
 * Example of simple usage
 * ------------------------
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * 
 * icarus::WaveformRMS const RMS { 3.2f };
 * 
 * std::cout << "RMS: " << RMS << " ADC" << std::endl;
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will print `RMS: 3.2 ADC`.
 */
struct icarus::WaveformRMS {
  
  using BaselineRMS_t = float; ///< Type of baseline RMS value.
  
  BaselineRMS_t fRMS {}; ///< The current value of the baseline RMS.
  
  
  // --- BEGIN -- Constructors -------------------------------------------------
  
  /// Constructor: default baseline RMS (`0`).
  WaveformRMS() = default;
  
  /// Constructor: sets the baseline RMS.
  WaveformRMS(BaselineRMS_t RMS): fRMS(RMS) {}
  
  // --- END -- Constructors ---------------------------------------------------
  
  
  // --- BEGIN -- Access to the baseline RMS -----------------------------------
  /**
   * @name Access to the baseline RMS
   * 
   * In addition to the direct method (`RMS()`) a few candies are offered:
   * a function-like operator for baseline RMS access.
   * 
   */
  /// @{
  /// Returns the current baseline RMS value.
  BaselineRMS_t RMS() const { return fRMS; }
  
  /// Returns the current baseline RMS value.
  BaselineRMS_t operator() () const { return RMS(); }
  
  // --- END -- Access to the baseline RMS -------------------------------------
  
}; // icarus::trigger::WaveformRMS


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_DATA_WAVEFORMRMS_H
