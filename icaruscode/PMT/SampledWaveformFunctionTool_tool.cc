/**
 * @file   icaruscode/PMT/SampledWaveformFunctionTool_tool.cc
 * @brief  Toolization of `icarus::opdet::SampledWaveformFunction<nanosecond>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/SampledWaveformFunction.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::SinglePhotonPulseFunctionTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/SinglePhotonPulseFunctionTool.h"
#include "icaruscode/PMT/Algorithms/SampledWaveformFunction.h"
#include "icaruscode/PMT/Algorithms/KeyValueParser.h"
#include "icaruscode/Decode/DecoderTools/details/KeyValuesData.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electromagnetism.h" // picocoulomb
#include "lardataalg/Utilities/quantities_fhicl.h" // nanoseconds from FHiCL
#include "lardataalg/Utilities/quantities.h" // util::quantities::makeQuantity()
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <fstream>
#include <memory> // std::unique_ptr()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::opdet { struct SampledWaveformFunctionTool; }
/**
 * @brief Creates a `SampledWaveformFunction` pulse shape.
 * @see `icarus::opdet::SinglePhotonPulseFunctionTool`
 * 
 * This tool creates a `icarus::opdet::SampledWaveformFunction<nanosecond>`
 * function to describe a R5912 PMT pulse attached to the ICARUS detector.
 * 
 * See `icarus::opdet::SampledWaveformFunction` for the details of the function.
 * 
 * 
 * Waveform specification file format
 * -----------------------------------
 * 
 * The response must be described in a plain text file following the syntax from
 * `icarus::details::KeyValueParser`.
 * The following fields are supported:
 * 
 * * `"FileFormat"` (integer, implied: `1`) represents the version of the data
 *   format, i.e. the list of the supported fields and their meaning.
 *   Currently it's only a formal parameter, which is ignored.
 * * `"Name"` (string): short name identifying this response
 * * `"Description"` (string): the description of this response; may be long and
 *   spanning multiple lines (always adhering `icarus::details::KeyValueParser`
 *   syntax)
 * * `"Date"` (string): the date of this response; free-form.
 * * `"Version"` (positive integer): the version of this response; it may
 *   describe updates for the same response `Name`.
 * * `"Tick"` (time quantity string, mandatory): the duration of one tick in the
 *   response sampling, in time quantity format (e.g. `"2 ns"` or "0.4 us").
 * * `"Gain"` (real number): the PMT gain associated to this response.
 *   If provided, it will allow rescaling to different gains.
 * * `"Samples"` (sequence of samples in mV, mandatory): values of the samples
 *   in the response, one per tick. No reference time is needed (the algorithms
 *   will look for the peak sample to be used as reference time).
 * * `"NSamples"` (positive integer): the number of samples in the `Samples`
 *   array; used to validate the input.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description SampledWaveformFunctionTool` (or read `Config`
 * data structure) for a short explanation of the meaning of the parameters.
 * 
 * * `TransitTime` (time, mandatory): time from the arrival of a photoelectron
 *     to the surface of the photodetector to the peak of the signal that the
 *     photodetector produces. It must include the unit (e.g. `"51.5 ns"`).
 * * `WaveformData` (path, mandatory): path to the file with the complete
 *     information about the single photoelectron response. The file is searched
 *     for in `FW_SEARCH_PATH` path
 * * `Gain` (real, optional): if specified, the input must provide the nominal
 *     gain of the response, and that response will be rescaled from its gain
 *     value to the one specified with this parameter. If not specified,
 *     the response is used as is.
 *
 */
struct icarus::opdet::SampledWaveformFunctionTool
  : public icarus::opdet::SinglePhotonPulseFunctionTool
{
  
  /// Configuration parameters.
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> WaveformData {
      Name("WaveformData"),
      Comment("Path to the data file with the SPR information")
      // mandatory
      };
    fhicl::Atom<nanoseconds> TransitTime {
      Name("TransitTime"),
      Comment("peak time from the beginning of the waveform [ns]")
      // mandatory
      };
    fhicl::OptionalAtom<float> Gain {
      Name("Gain"),
      Comment("PMT amplification gain factor")
      };
    
  }; // struct Config

  
  /// Tool parameter configuration.
  using Parameters = art::ToolConfigTable<Config>;
  
  /// Constructor: sets the configuration.
  SampledWaveformFunctionTool(Parameters const& config)
    : fPulseFunction(makePulseFunction(config())) {}
  
  
    private:
  
  /// The actual function type we offer.
  using MyFunction_t = icarus::opdet::SampledWaveformFunction<nanoseconds>;
  
  
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
  
  /// Parses the specified file and returns the information on the SPR waveform.
  static MyFunction_t::WaveformSpecs_t extractWaveformSpecification
    (std::string const& path);
  
}; // icarus::opdet::SampledWaveformFunctionTool


//------------------------------------------------------------------------------
//--- icarus::opdet::SampledWaveformFunctionTool implementation
//------------------------------------------------------------------------------

auto icarus::opdet::SampledWaveformFunctionTool::extractWaveformSpecification
  (std::string const& path) -> MyFunction_t::WaveformSpecs_t
{
  //
  // text file parsing
  //
  std::ifstream srcFile { path };
  if (!srcFile.is_open()) {
    // quite strange, actually, since the file was found by `cet::search_path`
    throw art::Exception{ art::errors::FileReadError }
      << "Can't open single photoelectron response file '" << path << "'\n";
  }
  
  icarus::details::KeyValueParser const parser;
  icarus::KeyValuesData const data { parser(srcFile) };
  srcFile.close();
  
  //
  // interpretation
  //
  
  auto makeException = [path]()
    {
      return cet::exception{ "SampledWaveformFunctionTool" }
        << "in '" << path << "': ";
    };
  
  MyFunction_t::WaveformSpecs_t specs;
  
  if (auto const* item = data.findItem("Name")) {
    if (item->nValues() != 1) {
      throw makeException() << "'Name' must have exactly 1 entry, not "
        << item->nValues() << "!\n";
    }
    specs.name = item->value();
  }
  
  if (auto const* item = data.findItem("Description")) {
    if (item->nValues() != 1) {
      throw makeException()
        << "'Description' must have exactly 1 entry (possibly quoted), not "
        << item->nValues() << "!\n";
    }
    specs.description = item->value();
  }
  
  if (auto const* item = data.findItem("Date")) {
    if (item->nValues() != 1) {
      throw makeException()
        << "'Date' must have exactly 1 entry (possibly quoted), not "
        << item->nValues() << "!\n";
    }
    specs.date = item->value();
  }
  
  if (auto const* item = data.findItem("Version")) {
    if (item->nValues() != 1) {
      throw makeException()
        << "'Version' must have exactly 1 entry, not " << item->nValues()
        << "!\n";
    }
    try {
      specs.version = item->getNumber<unsigned int>(0);
    }
    catch (icarus::KeyValuesData::Error const& e) {
      throw makeException() << "value in 'Version' ('" << item->value()
        << "') can't be interpreted as version number (unsigned int):\n"
        << e.what() << "\n";
    }
  }
  
  if (auto const* item = data.findItem("Tick")) {
    if (item->nValues() != 1) {
      throw makeException()
        << "'Date' must have exactly 1 entry (possibly quoted), not "
        << item->nValues() << "!\n";
    }
    try {
      specs.sampleDuration
        = util::quantities::makeQuantity<nanoseconds>(item->value());
    }
    catch(std::runtime_error const& e) {
      throw makeException()
        << "Failed to parse 'Tick' ('" << item->value() << "') as a time:\n"
        << e.what() << "\n";
    }
  }
  else throw makeException() << "'Tick' entry is mandatory.\n";
  
  if (auto const* item = data.findItem("Gain")) {
    if (item->nValues() != 1) {
      throw makeException()
        << "'Gain' must have exactly 1 entry, not " << item->nValues() << "!\n";
    }
    try {
      specs.gain = item->getNumber<float>(0);
    }
    catch (icarus::KeyValuesData::Error const& e) {
      throw makeException() << "value in 'Gain' ('" << item->value()
        << "') can't be interpreted as a gain factor:\n"
        << e.what() << "\n";
    }
  }
  
  
  if (auto const* item = data.findItem("Samples")) {
    if (item->nValues() < 2) {
      throw makeException()
        << "'Samples' has only " << item->nValues() << " values!!\n";
    }
    try {
      specs.samples = item->getVector<float>();
    }
    catch(icarus::KeyValuesData::Error const& e) {
      throw makeException()
        << "Error converting values in 'Samples' into voltage:\n"
        << e.what() << "\n";
    }
  }
  else throw makeException() << "'Samples' entry is mandatory.\n";
  
  
  if (auto const* item = data.findItem("NSamples")) {
    if (item->nValues() != 1) {
      throw makeException()
        << "'NSamples' must have exactly 1 entry, not " << item->nValues()
        << "!\n";
    }
    unsigned int nSamples;
    try {
      nSamples = item->getNumber<unsigned int>(0);
    }
    catch (icarus::KeyValuesData::Error const& e) {
      throw makeException() << "value in 'NSamples' ('" << item->value()
        << "') can't be interpreted as a sample number (unsigned int):\n"
        << e.what() << "\n";
    }
    if (specs.samples.size() != nSamples) {
      throw makeException() << "'Samples' has " << specs.samples.size()
        << " values, but according to 'NSamples' there should be " << nSamples
        << "!\n";
    }
  }
  
  
  return specs;
} // icarus::opdet::SampledWaveformFunctionTool::extractWaveformSpecification()

//------------------------------------------------------------------------------
auto icarus::opdet::SampledWaveformFunctionTool::makePulseFunction
  (Config const& config) -> std::unique_ptr<PulseFunction_t>
{
  
  //
  // find and process the waveform information file
  //
  cet::search_path searchPath{ "FW_SEARCH_PATH" };
  std::string waveformSpecsPath;
  try {
    waveformSpecsPath = searchPath.find_file(config.WaveformData());
  }
  catch (cet::exception& e) {
    throw art::Exception{ art::errors::Configuration, "", e }
      << "Error looking for the waveform data file '" << config.WaveformData()
      << "' (configured via: '" << config.WaveformData.name() << "')\n";
  }
  
  MyFunction_t::WaveformSpecs_t waveformSpecs
    = extractWaveformSpecification(waveformSpecsPath);
    
  //
  // fix the gain request
  //
  if ((config.Gain().value_or(0.0) != 0.0) && (waveformSpecs.gain == 0.0)) {
    throw art::Exception(art::errors::Configuration)
      << "The single photoelectron response '" << waveformSpecs.name
      << "' at '" << waveformSpecsPath
      << "' does not specify a base gain, so it can't be rescaled to "
      << config.Gain().value()
      << " ('" << config.Gain.name() << "' parameter)\n";
  }
  float const reqGain  = (config.Gain().value_or(0.0) != 0.0)
    ? config.Gain().value(): waveformSpecs.gain;
  
  //
  // create the algorithm
  //
  return std::make_unique<MyFunction_t>(
      std::move(waveformSpecs)      // waveformSpecs
    , config.TransitTime()          // peakTime
    , reqGain                       // gain
    );
  
} // icarus::opdet::SampledWaveformFunctionTool::makePulseFunction()


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::SampledWaveformFunctionTool)


//------------------------------------------------------------------------------
