/**
 * @file   icaruscode/PMT/OpReco/Algorithms/PedAlgoFixed.cxx
 * @brief  Pedestal "algorithm" reading the pedestals from somewhere else.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 7, 2022
 *
 * @ingroup PulseReco
 */

// ICARUS libraries
#include "icaruscode/PMT/OpReco/Algorithms/PedAlgoFixed.h"

// framework libraries
#ifndef ICARUS_NO_CONFIGURATION_VALIDATION
# include "fhiclcpp/ParameterSet.h"
#endif // ICARUS_NO_CONFIGURATION_VALIDATION

// C/C++ standard libraries
#include <algorithm> // std::fill(), std::find()
#include <iterator> // std::distance()
#include <stdexcept> // std::runtime_error
#include <limits> // std::numeric_limits<>
#include <cstdint> // std::ptrdiff_t


// -----------------------------------------------------------------------------
#ifndef ICARUS_NO_CONFIGURATION_VALIDATION
pmtana::PedAlgoFixed::PedAlgoFixed
  (Parameters const& params, std::string const& name /* = "PedFixed" */)
  : pmtana::PMTPedestalBase(name)
  , fWaveformTag{ params().WaveformTag() }
  , fPedestalTag{ params().PedestalTag() }
  , fRMSTag{ params().RMSTag().value_or(fPedestalTag) }
  {}

#else // ICARUS_NO_CONFIGURATION_VALIDATION
pmtana::PedAlgoFixed::PedAlgoFixed(
  fhicl::ParameterSet const& pset,
  std::string const& name /* = "PedFixed" */
)
  : pmtana::PMTPedestalBase(name)
  , fWaveformTag{ pset.get<std::string>("WaveformTag") }
  , fPedestalTag{ pset.get<std::string>("PedestalTag") }
  , fRMSTag{ pset.get<std::string>("RMSTag", fPedestalTag) }
  {}

#endif // ICARUS_NO_CONFIGURATION_VALIDATION

// -----------------------------------------------------------------------------
void pmtana::PedAlgoFixed::setParameters(InputSet_t inputSet) {
  if (
    (inputSet.pedestals.size() != inputSet.waveforms.size())
    || (inputSet.RMSs.size() != inputSet.waveforms.size())
  ) {
    throw std::runtime_error{
      "pmtana::PedAlgoFixed::setParameters(): inconsistent sizes of input: "
      + std::to_string(inputSet.waveforms.size()) + " waveforms, "
      + std::to_string(inputSet.pedestals.size()) + " pedestals, "
      + std::to_string(inputSet.RMSs.size()) + " RMS.\n"
      };
  }
  
  fInput = std::move(inputSet);
  
} // pmtana::PedAlgoFixed::setParameters()


// -----------------------------------------------------------------------------
void pmtana::PedAlgoFixed::clearParameters() {
  fInput = {};
} // pmtana::PedAlgoFixed::clearParameters()


// -----------------------------------------------------------------------------
bool pmtana::PedAlgoFixed::ComputePedestal(
  const pmtana::Waveform_t& wf,
  pmtana::PedestalMean_t& mean_v,
  pmtana::PedestalSigma_t& sigma_v
) {
  
  auto const [ inputSet, waveformIndex ] = findWaveform(wf);
  if (!inputSet) {
    // this means that this work is hopeless
    throw std::logic_error{
      "pmtana::PedAlgoFixed::ComputePedestal(): input waveform (address="
      + std::to_string(reinterpret_cast<std::intptr_t>(&wf))
      + ") not registered in any input (ID: "
      + std::to_string(fInput.inputID) + ")"
      };
  }
  
  std::size_t const nSamples = inputSet->waveforms[waveformIndex]->size();
  
  mean_v.resize(nSamples);
  std::fill(mean_v.begin(), mean_v.end(), inputSet->pedestals[waveformIndex]);
  
  sigma_v.resize(nSamples);
  std::fill(sigma_v.begin(), sigma_v.end(), inputSet->RMSs[waveformIndex]);
  
  return true;

} // pmtana::PedAlgoFixed::ComputePedestal()


// -----------------------------------------------------------------------------
auto pmtana::PedAlgoFixed::findWaveform
  (pmtana::Waveform_t const& waveform) const
  -> std::pair<InputSet_t const*, std::size_t>
{
  InputSet_t const* inputSet = &fInput;
  auto const wbegin = inputSet->waveforms.begin();
  auto const wend = inputSet->waveforms.end();
  if (auto it = std::find(wbegin, wend, &waveform); it == wend)
    return { nullptr, std::numeric_limits<std::size_t>::max() };
  else
    return { inputSet, std::distance(wbegin, it) };
  
} // pmtana::PedAlgoFixed::findWaveform()


// -----------------------------------------------------------------------------
