/**
 * @file   icaruscode/PMT/CustomPulseFunctionTool_tool.cc
 * @brief  Toolization of `icarus::opdet::AsymGaussPulseFunction<nanosecond>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/CustomPulseFunction.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::SinglePhotonPulseFunctionTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/SinglePhotonPulseFunctionTool.h"
#include "icaruscode/PMT/Algorithms/CustomPulseFunction.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electromagnetism.h" // picocoulomb
#include "lardataalg/Utilities/quantities_fhicl.h" // nanoseconds from FHiCL

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::opdet { struct CustomPulseFunctionTool; }
/**
 * @brief Creates a `CustomPulseFunction` pulse shape.
 * @see `icarus::opdet::SinglePhotonPulseFunctionTool`
 * 
 * This tool creates a `icarus::opdet::CustomPulseFunction<nanosecond>`
 * function to describe an arbitrary pulse shape.
 * 
 * See `icarus::opdet::CustomPulseFunction` for the details of the function.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description CustomPulseFunctionTool` (or read `Config`
 * data structure) for a short explanation of the meaning of the parameters.
 * 
 * * **ShapeFormula** (string, mandatory): an expression for the pulse shape;
 *   for example: `[A] * exp(-0.5*((x - [mu])/(sqrt2*[sigma]))**2)` is an
 *   extended way to describe a Gaussian pulse; the syntax is mostly C++ with a
 *   few ROOT extensions (see ROOT 6 `TFormula` documentation); `x` variable
 *   represents the time in nanoseconds, with `x = 0` the time of emission of
 *   the photoelectron;
 * * **PeakTime** (string, mandatory): evaluates to the time at which the pulse
 *   has its peak, in nanoseconds; the format can be an expression in the same
 *   fashion as `ShapeFormula` and it can use any of the parameters in
 *   `ShapeFormula`, but no extra parameters and with no variable, for example
 *   `[mu]` for the Gaussian pulse peak; or it can be a time, for example
 *   `55.1 ns`;
 * * **Parameters** (string): a table with one entry for each parameter of the
 *   pulse shape, and their numerical value (no unit is used); for example,
 *   `mu: 55.1 sigma: 2.4 A: -10.1` sets the parameters of the Gaussian shape
 *   in the previous example. This parameter must contain one entry for each
 *   parameter in `ShapeFormula` and no extra values; if `ShapeFormula` has
 *   no parameters, the `Parameters` table can be omitted.
 * 
 * @note Because of the limitations of FHiCL language in `Parameters`
 *       specification, the names of the parameters need to be simple (e.g.
 *       `[mu]` rather than `[#mu]`).
 */
struct icarus::opdet::CustomPulseFunctionTool
  : public icarus::opdet::SinglePhotonPulseFunctionTool
{
  
  /// Configuration parameters.
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> ShapeFormula {
      Name("ShapeFormula"),
      Comment("formuls for the pulse shape, in ROOT TFormula format")
      // mandatory
      };
    fhicl::Atom<std::string> PeakTime {
      Name("PeakTime"),
      Comment("time at which the pulse peaks [ns]")
      // mandatory
      };
    fhicl::OptionalDelegatedParameter Parameters {
      Name("Parameters"),
      Comment("collection of parameter names and their numerical values")
      };
    
  }; // struct Config

  
  /// Tool parameter configuration.
  using Parameters = art::ToolConfigTable<Config>;
  
  /// Constructor: sets the configuration.
  CustomPulseFunctionTool(Parameters const& config)
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
  
  
}; // icarus::opdet::CustomPulseFunctionTool


//------------------------------------------------------------------------------
//--- icarus::opdet::CustomPulseFunctionTool implementation
//------------------------------------------------------------------------------
auto icarus::opdet::CustomPulseFunctionTool::makePulseFunction
  (Config const& config) -> std::unique_ptr<PulseFunction_t>
{
  
  using MyFunction_t = icarus::opdet::CustomPulseFunction<nanoseconds>;
  
  std::string const& expression = config.ShapeFormula();
  
  // collect the parameters
  fhicl::ParameterSet configuredParameters; // will stay empty if not present
  config.Parameters.get_if_present(configuredParameters);
  MyFunction_t::PulseParameters_t parameters;
  for (auto const& parName: configuredParameters.get_names())
    parameters.emplace_back(parName, configuredParameters.get<double>(parName));
  
  std::string const& peakTimeStr = config.PeakTime();
  
  try {
    // try to convert `peakTimeStr` into a nanosecond quantity
    return std::make_unique<MyFunction_t>(
      expression,
      parameters,
      util::quantities::makeQuantity<nanoseconds>(peakTimeStr, true)
      );
  }
  catch (util::quantities::ValueError const&) {}
  catch (util::quantities::ExtraCharactersError const&) {}
  
  mf::LogDebug("CustomPulseFunctionTool")
    << "Parameter PeakTime ('" << peakTimeStr
    << "') does not seem to be a nanosecond quantity; will try as formula."
    ;
  
  return std::make_unique<MyFunction_t>(expression, parameters, peakTimeStr);
  
} // icarus::opdet::CustomPulseFunctionTool::makePulseFunction()


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::CustomPulseFunctionTool)


//------------------------------------------------------------------------------

