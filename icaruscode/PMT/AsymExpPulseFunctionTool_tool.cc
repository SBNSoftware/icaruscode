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
#include "canvas/Utilities/Exception.h"

// framework libraries
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr()
#include <cmath> // std::exp(), std::log()
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
 * Finally, note that it is also possible to specify the amplitude of the
 * function in terms of PMT gain (`Gain` parameter) instead of current
 * (`PeakAmplitude`), in which case the peak current will be calculated to
 * have a correct total charge.
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
    fhicl::OptionalAtom<float> PeakCurrent {
      Name("PeakCurrent"),
      Comment("measured charge current at peak [mA]")
      };
    fhicl::OptionalAtom<float> Gain {
      Name("Gain"),
      Comment("PMT amplification gain factor")
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
  
  /// Computes the peak current value out of the shape and gain parameters [mA].
  static float peakCurrentFromGain
    (float gain, nanoseconds raiseTau, nanoseconds fallTau);
  
}; // icarus::opdet::AsymExpPulseFunctionTool


//------------------------------------------------------------------------------
//--- icarus::opdet::AsymExpPulseFunctionTool implementation
//------------------------------------------------------------------------------
auto icarus::opdet::AsymExpPulseFunctionTool::makePulseFunction
  (Config const& config) -> std::unique_ptr<PulseFunction_t>
{
  
  using MyFunction_t = icarus::opdet::AsymExpPulseFunction<nanoseconds>;
  using ADCcount = MyFunction_t::ADCcount;
  
  // FIXME milliampere not available yet in `electromagnetism.h`
  double peakCurrent;
  if (config.Gain().has_value()) { // peak current from gain
    if (config.PeakCurrent().has_value()) {
      throw art::Exception(art::errors::Configuration)
        << "AsymExpPulseFunctionTool: Only one configuration parameter out of '"
        << config.PeakCurrent.name() << "' (set as " << (*config.PeakCurrent())
        << ") and '" << config.Gain.name() << "' (set as " << (*config.Gain())
        << ") can be specified at once!\n";
    }
    peakCurrent = peakCurrentFromGain(
      *config.Gain(),             // gain
      config.RaiseTimeConstant(), // raiseTau
      config.FallTimeConstant()   // fallTau
      );
  }
  else { // peak current directly specified
    if (!config.PeakCurrent().has_value()) {
      throw art::Exception(art::errors::Configuration)
        << "AsymExpPulseFunctionTool: either configuration parameter '"
        << config.PeakCurrent.name() << "' or '" << config.Gain.name()
        << "' must be specified!\n";
    }
    peakCurrent = *config.PeakCurrent();
  }
  
  return std::make_unique<MyFunction_t>(
    // amplitude is a charge, so we have to twist the arm of the constructor to
    // accept it as ADC count (`value()` makes `peakCurrent` lose its unit)
//     ADCcount(config.ADC() * config.peakCurrent.value()),   // amplitude
    ADCcount(config.ADC() * peakCurrent),                  // amplitude
    config.TransitTime(),                                  // peakTime
    config.RaiseTimeConstant(),                            // raiseTau
    config.FallTimeConstant()                              // fallTau
    );
  
} // icarus::opdet::AsymExpPulseFunctionTool::makePulseFunction()


//------------------------------------------------------------------------------
float icarus::opdet::AsymExpPulseFunctionTool::peakCurrentFromGain
  (float gain, nanoseconds raiseTau, nanoseconds fallTau)
{
  static constexpr double electronCharge = -1.602176634e-7; // picocoulomb
  
  nanoseconds const dTau = fallTau - raiseTau; // expected positive
  float const tauRatioLog = std::log(raiseTau / fallTau); // expected negative
  
  return gain * electronCharge / dTau.value() * (
      std::exp(raiseTau / dTau * tauRatioLog)
    - std::exp(fallTau  / dTau * tauRatioLog)
    );
  
} // icarus::opdet::AsymExpPulseFunctionTool::peakCurrentFromGain()


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::AsymExpPulseFunctionTool)


//------------------------------------------------------------------------------

