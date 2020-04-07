/**
 * @file   icaruscode/PMT/AsymGaussPulseFunctionTool_tool.cc
 * @brief  Toolization of `icarus::opdet::AsymGaussPulseFunction<nanosecond>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/AsymGaussPulseFunction.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::SinglePhotonPulseFunctionTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/SinglePhotonPulseFunctionTool.h"
#include "icaruscode/PMT/Algorithms/AsymGaussPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electromagnetism.h" // picocoulomb
#include "lardataalg/Utilities/quantities_fhicl.h" // nanoseconds from FHiCL

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::opdet { struct AsymGaussPulseFunctionTool; }
/**
 * @brief Creates a `AsymGaussPulseFunction` pulse shape.
 * @see `icarus::opdet::SinglePhotonPulseFunctionTool`
 * 
 * This tool creates a `icarus::opdet::AsymGaussPulseFunction<nanosecond>`
 * function to describe a R5912 PMT pulse.
 * 
 * See `icarus::opdet::AsymGaussPulseFunction` for the details of the function.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description AsymGaussPulseFunctionTool` (or read `Config`
 * data structure) for a short explanation of the meaning of the parameters.
 * 
 * In addition, note that the actual amplitude in ADC counts of the pulse is
 * composed as the product of the amplitude in charge (`MeanAmplitude`) and
 * the charge-to-ADC conversion factor (`ADC`).
 * 
 */
struct icarus::opdet::AsymGaussPulseFunctionTool
  : public icarus::opdet::SinglePhotonPulseFunctionTool
{
  
  /// Configuration parameters.
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<nanoseconds> TransitTime {
      Name("TransitTime"),
      Comment("peak time from the beginning of the waveform [ns]")
      // mandatory
      };
    fhicl::Atom<util::quantities::picocoulomb> MeanAmplitude {
      Name("MeanAmplitude"),
      Comment("signal amplitude at peak [pC]")
      // mandatory
      };
    fhicl::Atom<nanoseconds> RaiseTime {
      Name("RaiseTime"),
      Comment("rise time (10% to 90%, sigma * ~1.687) [ns]")
      // mandatory
      };
    fhicl::Atom<nanoseconds> FallTime {
      Name("FallTime"),
      Comment("fall time (90% to 10%, sigma * ~1.687) [ns]")
      // mandatory
      };
    fhicl::Atom<float> ADC {
      Name("ADC"),
      Comment("Charge to ADC conversion factor [ADC counts/pC]")
      // mandatory
      };
    
  }; // struct Config

  
  /// Tool parameter configuration.
  using Parameters = art::ToolConfigTable<Config>;
  
  /// Constructor: sets the configuration.
  AsymGaussPulseFunctionTool(Parameters const& config)
    : fPulseFunction(makePulseFunction(config())) {}
  
  
    private:
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /// Returns the function that was created at construction time.
  virtual std::unique_ptr<PulseFunction_t> doGetPulseFunction() override
    { return std::move(fPulseFunction); }
  
  // --- END -- Virtual interface ----------------------------------------------
  
  /// Function stored while waiting to be delivered.
  std::unique_ptr<PulseFunction_t> fPulseFunction;
  
  
  /// Creates and returns a pulse function with the specified configuration.
  static std::unique_ptr<PulseFunction_t> makePulseFunction
    (Config const& config);
  
  
}; // icarus::opdet::AsymGaussPulseFunctionTool


//------------------------------------------------------------------------------
//--- icarus::opdet::AsymGaussPulseFunctionTool implementation
//------------------------------------------------------------------------------
auto icarus::opdet::AsymGaussPulseFunctionTool::makePulseFunction
  (Config const& config) -> std::unique_ptr<PulseFunction_t>
{
  
  using MyFunction_t = icarus::opdet::AsymGaussPulseFunction<nanoseconds>;
  using ADCcount = MyFunction_t::ADCcount;
  
  auto raiseTimeToRMS = [](auto raiseTime)
    {
      return raiseTime / (
        std::sqrt(2.0)
        * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9)))
        );
    };
  
  return std::make_unique<MyFunction_t>(
    // amplitude is a charge, so we have to twist the arm of the constructor to
    // accept it as ADC count (`value()` makes `meanAmplitude` lose its unit)
    ADCcount(config.ADC() * config.MeanAmplitude().value()), // amplitude
    config.TransitTime(),               // peakTime
    raiseTimeToRMS(config.RaiseTime()), // sigmaLeft
    raiseTimeToRMS(config.FallTime())   // sigmaRight
    );
  
} // icarus::opdet::AsymGaussPulseFunctionTool::makePulseFunction()


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::AsymGaussPulseFunctionTool)


//------------------------------------------------------------------------------

