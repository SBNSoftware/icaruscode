/**
 * @file   icaruscode/PMT/SinglePhotonPulseFunctionTool.h
 * @brief  Tool to create a pulse function for PMT single photon response.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h`
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_PMT_SINGLEPHOTONPULSEFUNCTIONTOOL_H
#define ICARUSCODE_PMT_SINGLEPHOTONPULSEFUNCTIONTOOL_H


// ICARUS libraries
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond

// C/C++ standard libraries
#include <memory> // std::unique_ptr()
#include <utility> // std::move()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::opdet { struct SinglePhotonPulseFunctionTool; }
/**
 * @brief Creates a `PhotoelectronPulseFunction` pulse shape.
 * 
 * Implementations of this tool create a function representing the single photon
 * response of a PMT, as an object derived from
 * `icarus::opdet::PhotoelectronPulseFunction<nanosecond>`.
 * 
 * The envisioned usage of this tool is from a temporary, one-time-only call:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * using util::quantities::nanosecond;
 * std::unique_ptr<icarus::opdet::PhotoelectronPulseFunction<nanosecond>> pulse
 *   = art::make_tool<icarus::opdet::SinglePhotonPulseFunctionTool>(config)
 *     ->getPulseFunction()
 *   ;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * where `pset` is the FHiCL parameter set configuration for the tool
 * (the type of `pulse` can be safely declared `auto`).
 * The tool object itself is useless after the first use, and should be
 * disposed of.
 * 
 * The configuration of the tool is implementation-dependent.
 * 
 */
struct icarus::opdet::SinglePhotonPulseFunctionTool {
  
  /// Convenience definition for time stored in nanoseconds.
  using nanoseconds = util::quantities::nanosecond;
  
  /// Type of function returned.
  using PulseFunction_t
    = icarus::opdet::PhotoelectronPulseFunction<nanoseconds>;
  
  
  /**
   * @brief Returns an instance of the pulse function.
   * 
   * This function is guaranteed to return a valid pulse function only once.
   * The return value on following calls is undefined, but it is considered
   * an error to return a null pointer.
   * 
   * See the description of the class for an usage example.
   */
  std::unique_ptr<PulseFunction_t> getPulseFunction();
  
  
    private:
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /**
   * @brief Returns an instance of the pulse function.
   * 
   * Only the first call is guaranteed to return a non-null pointer.
   * 
   */
  virtual std::unique_ptr<PulseFunction_t> doGetPulseFunction() = 0;

  
  // --- END -- Virtual interface ----------------------------------------------
  
  
}; // icarus::opdet::SinglePhotonPulseFunctionTool


//------------------------------------------------------------------------------
//--- icarus::opdet::SinglePhotonPulseFunctionTool implementation
//------------------------------------------------------------------------------
inline auto icarus::opdet::SinglePhotonPulseFunctionTool::getPulseFunction()
  -> std::unique_ptr<PulseFunction_t>
{
  auto ptr = doGetPulseFunction();
  assert(ptr);
  return std::move(ptr);
} // icarus::opdet::SinglePhotonPulseFunctionTool::getPulseFunction() &&


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_SINGLEPHOTONPULSEFUNCTIONTOOL_H
