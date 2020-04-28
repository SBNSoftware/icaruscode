/**
 * @file   icaruscode/PMT/AsymExpPulseFunctionTool_tool.cc
 * @brief  Toolization of `icarus::opdet::AsymExpPulseFunction<nanosecond>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/AsymExpPulseFunction.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::SinglePhotonPulseFunctionTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/SinglePhotonPulseFunctionTool.h"
#include "icaruscode/PMT/Algorithms/AsymExpPulseFunction.h"

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
namespace icarus::opdet { struct AsymExpPulseFunctionTool; }
/**
 * @brief Creates a `AsymExpPulseFunction` pulse shape.
 * @see `icarus::opdet::SinglePhotonPulseFunctionTool`
 * 
 * This tool creates a `icarus::opdet::AsymExpPulseFunction<nanosecond>`
 * function to describe a R5912 PMT pulse.
 * 
 * See `icarus::opdet::AsymExpPulseFunction` for the details of the function.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description AsymExpPulseFunctionTool` (or read `Config`
 * data structure) for a short explanation of the meaning of the parameters.
 * 
 * In addition, note that the actual amplitude in ADC counts of the pulse is
 * composed as the product of the amplitude in charge (`PeakCharge`) and
 * the charge-to-ADC conversion factor (`ADC`).
 * The latter can be considered as the product of the circuitry impedance
 * (transforming the charge into a voltage) and the digitization conversion
 * (full digitizer range divided by the largest output value).
 * 
 */
struct icarus::opdet::AsymExpPulseFunctionTool
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
    // TODO add milliampere to `electromagnetism.h`
    fhicl::Atom<float> PeakCurrent {
      Name("PeakCurrent"),
      Comment("measured charge current at peak [mA]")
      // mandatory
      };
    fhicl::Atom<nanoseconds> RaiseTimeConstant {
      Name("RaiseTimeConstant"),
      Comment("raise time constant (exponential raise) [ns]")
      // mandatory
      };
    fhicl::Atom<nanoseconds> FallTimeConstant {
      Name("FallTimeConstant"),
      Comment("fall time constant (exponential fall) [ns]")
      // mandatory
      };
    fhicl::Atom<float> ADC {
      Name("ADC"),
      Comment("Current-to-ADC conversion factor [ADC counts/mA]")
      // mandatory
      };
    
  }; // struct Config

  
  /// Tool parameter configuration.
  using Parameters = art::ToolConfigTable<Config>;
  
  /// Constructor: sets the configuration.
  AsymExpPulseFunctionTool(Parameters const& config)
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
  
  
}; // icarus::opdet::AsymExpPulseFunctionTool


//------------------------------------------------------------------------------
//--- icarus::opdet::AsymExpPulseFunctionTool implementation
//------------------------------------------------------------------------------
auto icarus::opdet::AsymExpPulseFunctionTool::makePulseFunction
  (Config const& config) -> std::unique_ptr<PulseFunction_t>
{
  
  using MyFunction_t = icarus::opdet::AsymExpPulseFunction<nanoseconds>;
  using ADCcount = MyFunction_t::ADCcount;
  
  return std::make_unique<MyFunction_t>(
    // amplitude is a charge, so we have to twist the arm of the constructor to
    // accept it as ADC count (`value()` makes `PeakCurrent()` lose its unit)
    // FIXME milliampere not available yet in `electromagnetism.h`
//     ADCcount(config.ADC() * config.PeakCurrent().value()), // amplitude
    ADCcount(config.ADC() * config.PeakCurrent()),         // amplitude
    config.TransitTime(),                                  // peakTime
    config.RaiseTimeConstant(),                            // raiseTau
    config.FallTimeConstant()                              // fallTau
    );
  
} // icarus::opdet::AsymExpPulseFunctionTool::makePulseFunction()


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::AsymExpPulseFunctionTool)


//------------------------------------------------------------------------------

