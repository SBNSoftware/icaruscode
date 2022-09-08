/**
 * @file   icaruscode/PMT/OpReco/Algorithms/PedAlgoFixed.h
 * @brief  Pedestal "algorithm" reading the pedestals from somewhere else.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 7, 2022
 *
 * @ingroup PulseReco
 */

#ifndef ICARUSCODE_PMT_OPRECO_ALGORITHMS_PEDALGOFIXED_H
#define ICARUSCODE_PMT_OPRECO_ALGORITHMS_PEDALGOFIXED_H

/// LArSoft libraries
#include "larana/OpticalDetector/OpHitFinder/PMTPedestalBase.h"
#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

// framework libraries
#ifndef ICARUS_NO_CONFIGURATION_VALIDATION
# include "fhiclcpp/types/Table.h"
# include "fhiclcpp/types/OptionalAtom.h"
# include "fhiclcpp/types/Atom.h"
#endif // ICARUS_NO_CONFIGURATION_VALIDATION

/// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair
#include <limits> // std::numeric_limits<>


// -----------------------------------------------------------------------------
namespace pmtana { class PedAlgoFixed; }
#ifdef ICARUS_NO_CONFIGURATION_VALIDATION
namespace fhicl { class ParameterSet; }
#else // ICARUS_NO_CONFIGURATION_VALIDATION

#endif // ICARUS_NO_CONFIGURATION_VALIDATION
/**
 * @class pmtana::PedAlgoFixed
 * @brief Pedestal "algorithm" reading the pedestals from somewhere else.
 * @addtogroup PulseReco
 * 
 * This is a big hack.
 * 
 * The optical reconstruction mini-framework passes to the pedestal algorithm
 * only a pointer to the waveform samples, with no context at all.
 * In this algorithm we need to figure out that context. We take the pointer of
 * every single possible input, because we need to know the index of the
 * waveform, which is required to match the index in the pedestal and RMS lists.
 * We can't even just take the address of the first waveform, because we don't
 * know the _actual_ type of the vector the waveform belongs too (in LArSoft,
 * that is `raw::OpDetWaveform`, which has different size and therefore pointer
 * math than the vector of samples it derives from).
 * 
 * 
 * Multithreading note
 * --------------------
 * 
 * This class is not multithread-ready. To process multiple events, the input
 * sets should be placed into a container with controlled access (one writes,
 * many read).
 * 
 */
class pmtana::PedAlgoFixed: public pmtana::PMTPedestalBase {
  
    public:

  /// Input information.
  struct InputSet_t {
    /// An unique identifier for this input set.
    unsigned int inputID = std::numeric_limits<unsigned int>::max();
    /// All the waveforms that are going to process.
    std::vector<pmtana::Waveform_t const*> waveforms;
    std::vector<float> pedestals; ///< In ADC counts by waveform.
    std::vector<float> RMSs; ///< In ADC counts by waveform.
    /// Returns the number of elements in the input.
    std::size_t size() const noexcept { return waveforms.size(); }
  }; // InputSet_t
  
#ifndef ICARUS_NO_CONFIGURATION_VALIDATION
  struct Config {
    fhicl::Atom<std::string> Name {
      fhicl::Name{ "Name" },
      fhicl::Comment
        { "Name of this algorithm [mandatory in the ophit miniframework]" },
      "Fixed" // default
      };
    fhicl::Atom<std::string> WaveformTag {
      fhicl::Name{ "WaveformTag" },
      fhicl::Comment{ "Tag of the data product with input waveforms." }
      };
    fhicl::Atom<std::string> PedestalTag {
      fhicl::Name{ "PedestalTag" },
      fhicl::Comment{ "Tag of the data product with input pedestal levels." }
      };
    fhicl::OptionalAtom<std::string> RMSTag {
      fhicl::Name{ "RMSTag" },
      fhicl::Comment{
        "Tag of the data product with input pedestal RMS [default: like "
        + PedestalTag.name() + "]." 
        }
      };
  }; // Config
  
  using Parameters = fhicl::Table<Config>;
  
  PedAlgoFixed(Parameters const& params, std::string const& name = "PedFixed");
  
#else // !ICARUS_NO_CONFIGURATION_VALIDATION
  /// Constructor from FHiCL configuration and optional class name
  PedAlgoFixed
    (fhicl::ParameterSet const& pset, std::string const& name = "PedFixed");
#endif // ICARUS_NO_CONFIGURATION_VALIDATION

  /// Returns the name of the configured waveform source.
  std::string const& waveformSourceName() const { return fWaveformTag; }

  /// Returns the name of the configured pedestal source.
  std::string const& pedestalSourceName() const { return fPedestalTag; }
  
  /// Returns the name of the configured RMS source.
  std::string const& pedestalRMSName() const { return fRMSTag; }
  
  // --- BEGIN -- Class-specific interface: setup ------------------------------
  
  /// Records the current pedestals and RMS per channel.
  void setParameters(InputSet_t inputSet);
  
  /// Removes the pedestal and RMS records.
  void clearParameters();
  
  // --- END ---- Class-specific interface: setup ------------------------------
  
    protected:
  
  // --- BEGIN -- Configuration ------------------------------------------------
  std::string fWaveformTag; ///< Name of the data source for waveforms.
  std::string fPedestalTag; ///< Name of the data source for pedestals.
  std::string fRMSTag; ///< Name of the data source for RMS.
  // --- END ---- Configuration ------------------------------------------------
  
  InputSet_t fInput; ///< The set of input waveforms currently registered.

  /**
   * @brief Computes the pedestal of the specified `waveform`.
   * @param[in] waveform the waveform to be processed
   * @param[out] mean_v the mean value of the pedestal, tick by tick
   * @param[out] sigma_v the RMS of the pedestal, tick by tick
   * @return whether the pedestal was successfully found
   * 
   * The algorithm is quite simple, since it picks the baseline that was
   * provided in the configuration.
   */
  virtual bool ComputePedestal(
    pmtana::Waveform_t const& waveform,
    pmtana::PedestalMean_t& mean_v,
    pmtana::PedestalSigma_t& sigma_v
    ) override;
  
  /// Returns the input set and index of the specified `waveform`.
  std::pair<InputSet_t const*, std::size_t> findWaveform
    (pmtana::Waveform_t const& waveform) const;

}; // pmtana::PedAlgoFixed

// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_OPRECO_ALGORITHMS_PEDALGOFIXED_H
