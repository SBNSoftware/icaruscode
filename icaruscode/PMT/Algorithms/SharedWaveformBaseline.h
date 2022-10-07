/**
 * @file   icaruscode/PMT/Algorithms/SharedWaveformBaseline.h
 * @brief  Extracts and writes PMT waveform baselines.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 5, 2022
 * @see    icaruscode/PMT/Algorithms/SharedWaveformBaseline.cxx
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_SHAREDWAVEFORMBASELINE_H
#define ICARUSCODE_PMT_ALGORITHMS_SHAREDWAVEFORMBASELINE_H


// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <vector>
#include <utility> // std::move()
#include <string>
#include <ostream>
#include <cstdint> // std::size_t


// -----------------------------------------------------------------------------
namespace opdet { class SharedWaveformBaseline; }
/**
 * @class opdet::SharedWaveformBaseline
 * @brief Extracts a common baseline from waveforms.
 * 
 * This algorithm processes a group of waveforms at a time, and returns a common
 * baseline for them.
 * The baseline is learned by looking at a fixed size of the beginning of each
 * of the waveforms, as follows:
 * 
 * 1. the RMS of the first portion of each baseline is computed;
 * 2. the median of the samples on the same samples is also computed;
 * 3. an acceptance range is constructed, using the median of the sample as the
 *    center and `nRMS` times the median of the RMS as maximum distance from
 *    that center in either direction;
 * 4. as a second pass, if in the first portion of a waveform there are at least
 *    `nExcessSamples` samples in a row that are outside of the acceptance
 *    range, that waveform is excluded;
 * 5. all the samples in the first portion of the remaining waveforms are
 *    averaged to obtain the final estimation of the baseline; this last step
 *    should increase the resolution of the baseline beyond the median that was
 *    obtained at step 2.
 * 
 * The parameters are specified at algorithm construction time and are contained
 * in the `Params_t` object.
 * 
 */
class opdet::SharedWaveformBaseline {
    public:
  
  /// Algorithm configuration parameters.
  struct Params_t {
    
    std::size_t nSample; ///< Number of samples to use from each waveform.
    
    double nRMS; ///< Number of RMS from the baseline to discard a waveform.
    
    /// Number of samples out of range to discard a waveform.
    unsigned int nExcessSamples;
    
    /// Dumps this configuration into the output stream `out`.
    template <typename Stream>
    void dump(
      Stream& out,
      std::string const& indent, std::string const& firstIndent
      ) const;
    template <typename Stream>
    void dump(Stream& out, std::string const& indent = "") const
      { dump(out, indent, indent); }
    
  }; // Params_t
  
  
  /// Type for algorithm result.
  struct BaselineInfo_t {
    
    /// Magic value used to denote the lack of a (`double`) data item.
    static constexpr double NoInfo = std::numeric_limits<double>::max();
    
    /// Value of the baseline [ADC#]
    double baseline = NoInfo;
    
    /// The RMS found during the extraction.
    double RMS = NoInfo;
    
    /// Number of waveforms used for the extraction.
    unsigned int nWaveforms = 0U;
    
    /// Number of samples used for the extraction.
    unsigned int nSamples = 0U;
    
  }; // BaselineInfo_t
  
  
  SharedWaveformBaseline(Params_t params, std::string logCategory):
      fParams{ std::move(params) }
    , fLogCategory{ std::move(logCategory) }
    {}
  
  /// Returns a common baseline from all the specified waveforms.
  BaselineInfo_t operator()
    (std::vector<raw::OpDetWaveform const*> const& waveforms) const;
  
  /// Returns the set of configuration parameters of this algorithm.
  Params_t const& parameters() const { return fParams; }
  
    private:
  Params_t fParams; ///< Algorithm parameters.
  
  std::string fLogCategory; ///< Name of stream category for console messages.
  
}; // opdet::SharedWaveformBaseline


//------------------------------------------------------------------------------
namespace opdet {
  
  inline std::ostream& operator<<
    (std::ostream& out, SharedWaveformBaseline::Params_t const& params)
    { params.dump(out); return out; }

} // namespace opdet


//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void opdet::SharedWaveformBaseline::Params_t::dump(
  Stream& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out << firstIndent << "samples from each waveforms: " << nSample
    << "\n" << indent << "pedestal range: +/- " << nRMS << " x RMS"
    << "\n" << indent << "use only waveforms with less than "
      << nExcessSamples << " samples out of pedestal range"
    ;
} // opdet::SharedWaveformBaseline::Params_t::dump()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_SHAREDWAVEFORMBASELINE_H
