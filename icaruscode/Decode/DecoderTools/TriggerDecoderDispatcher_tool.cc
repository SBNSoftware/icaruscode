/**
 * @file   icaruscode/Decode/DecodeTools/TriggerDecoderDispatched_tool.cc
 * @brief  Trigger decoder detecting which version of decoding to use.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu), Jacob Zettlemoyer
 * 
 * 
 */

// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"

// LArSoft and artDAQ libraries
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "artdaq-core/Data/Fragment.hh"

// framework libraries
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Run.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr
#include <regex>
#include <string>
#include <optional>
#include <vector>


// -----------------------------------------------------------------------------
namespace daq { class TriggerDecoderDispatched; }
/**
 * @brief Trigger decoding tool supporting multiple trigger versions.
 * 
 * This class delegates trigger decoding to another, single-purpose tool.
 * Which tool is chosen depends on the trigger configuration, and it may change
 * on each run (i.e. no access to any trigger data fragment is allowed to make
 * that choice).
 * 
 * A fragment generator is a part of the DAQ code that serves a hardware
 * component and generates data fragments out of it.
 * The trigger configuration contains the name of the generator used to create
 * the trigger fragment. That name is used to determine the configuration of the
 * trigger decoding tool ("worker" tool) that will be actually used. If there is
 * no configuration for that generator, an exception is typically thrown.
 * 
 * 
 * Configuration
 * --------------
 * 
 * * `TrigConfigLabel` (input tag, mandatory): the data product with the trigger
 *   configuration
 * * `Decoders` (list of generator configurations): each entry in this list is
 *   a generator configuration table structured as follows:
 *     * `Generator` (regular expression, mandatory): name of the generator this
 *       configuration deals with. It is a regular expression matched with
 *       `std::regex_match()`, but it can be specified as a simple name.
 *     * `FragmentsLabel` (input tag, mandatory): name of the data product with
 *       the trigger data fragment.
 *     * `ToolConfig` (configuration table): the configuration for the trigger
 *       decoder tool that can decode the data fragments written by the
 *       `Generator`.
 * * `LogCategory` (string, default: `"TriggerDecoderDispatched"`): name of
 *   the message facility stream where messages of this tool are sent to.
 * 
 * 
 * Input
 * ------
 * 
 * In addition to the input required by the delegated trigger decoder tools,
 * this tool also uses:
 * 
 * * `icarus::TriggerConfiguration` (`TrigConfigLabel`): configuration of
 *   trigger hardware
 *   (extracted e.g. by `icarus::TriggerConfigurationExtraction` module).
 * * `std::vector<artdaq::Fragment>` (`FragmentsLabel` of the worker tool(s)
 *   that is actually used).
 * 
 * 
 * Output
 * -------
 * 
 * This tool requires that each worker tool produces all (and only) the
 * following data products:
 * 
 * * `std::vector<raw::ExternalTrigger>` (empty instance name)
 * * `std::vector<raw::ExternalTrigger>` (instance name `"previous"`)
 * * `std::vector<raw::Trigger>`
 * * `std::vector<sim::BeamGateInfo>`
 * * `sbn::ExtraTriggerInfo`
 * 
 * The specifics of these data products should be described by the worker tools.
 * 
 * 
 * The `IDecoder` protocol
 * ------------------------
 * 
 * This object does not fully implement the `daq::IDecoder` protocol, nor it
 * supports all its functionalities.
 * 
 * The `preferredInput()` call will return the input data product tag for the
 * raw trigger data fragment. This is directly learnt from the tool
 * configuration (`FragmentsLabel`).
 * 
 * The `configure()` hook does not work and will throw an exception if invoked.
 * This tool only accepts configuration at construction time.
 * 
 * The `consumes()` functionality is not very robust. It must be invoked at
 * construction time, when it's not yet known which trigger generator(s) will
 * be encountered, and therefore which tools will be used. 
 * As a fallback, because the tool configurations in this module must include
 * the input tag of the data fragment they need, we can declare that all of them
 * _might_ be consumed, which is what this tool does. Most likely, all but one
 * of those declarations will not be confirmed.
 * 
 * Similarly, `produces()` is troublesome. This is not an optional protocol,
 * but at least the produce of the tools is well defined by LArSoft/ICARUS
 * conventions. So we stick to it and require that all worker tools we use will
 * honour that convention (documented in the `Outputs` section above).
 * 
 * 
 * ### Support for different triggers in the same job
 * 
 * The design of this tool should allow for the use of different trigger tools
 * in the same job. A limitation in the implementation in `icaruscode`
 * `v09_63_00` comes from the fact that the trigger decoding tools have all
 * the same fully qualified name, which likely makes it troublesome to load
 * more than one of them at the same time. This can be easily addressed by
 * renaming the tool C++ class names; since the tools are never called by name
 * in external code, this should have no effect in the big picture.
 * 
 * 
 * 
 */
class daq::TriggerDecoderDispatched: public daq::IDecoder {
  
    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    struct DecoderConfig {
      
      fhicl::Atom<std::string> Generator{
        Name{ "Generator" },
        Comment{ "the name of the trigger fragment generator (regex)" }
        };
      
      fhicl::Atom<art::InputTag> FragmentsLabel{
        Name{ "FragmentsLabel" },
        Comment{ "the trigger fragment input tag" }
        };
      
      fhicl::DelegatedParameter ToolConfig{
        Name{ "ToolConfig" },
        Comment{ "configuration of the tool for this generator" }
        };
      
    }; // DecoderConfig
    
    
    fhicl::Atom<art::InputTag> TrigConfigLabel{
      Name{ "TrigConfigLabel" },
      Comment{ "data product with the trigger hardware configuration" }
      };
    
    fhicl::Sequence<fhicl::Table<DecoderConfig>> Decoders{
      Name{ "Decoders" },
      Comment{ "Configurations of all supported decoders" }
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "name of stream for messages of this tool" },
      "TriggerDecoderDispatched"
      };
    
  }; // Config
  
  using Parameters = art::ToolConfigTable<Config>;
  
  explicit TriggerDecoderDispatched(Parameters const& params);
  
  /// Optionally declares which data products to read.
  virtual void consumes(art::ConsumesCollector& collector) override;
  
  /// Registers all the data products we promise to produce.
  virtual void produces(art::ProducesCollector& collector) override;
  
  virtual void configure(const fhicl::ParameterSet&) override;
  
  /// Delegates data product initialization to the current trigger decoder.
  virtual void initializeDataProducts() override;
  
  /// Constructs and prepares the tool according to the trigger in the run.
  virtual void setupRun(art::Run const& run) override;
  
  virtual std::optional<art::InputTag> preferredInput() const override;
  
  /// Delegates the fragment processing to the current trigger decoder.
  virtual void process_fragment(const artdaq::Fragment &fragment) override;
  
  /// Delegates data product output to the current trigger decoder.
  virtual void outputDataProducts(art::Event &event) override;
  
  
    private:
  
  /// Configuration of a single generator and tool to decode it.
  struct GeneratorParams_t {
    std::string generatorName; ///< Name of the generator (pattern text).
    std::regex generatorPattern; ///< Regular expression for generator name.
    art::InputTag fragmentTag; ///< Tag for the trigger input fragment.
    fhicl::ParameterSet toolConfig; ///< Configuration for the tool.
  };
  
  /// Standard name for data products referring to the current event.
  static std::string const CurrentTriggerInstanceName;
  
  /// Standard name for data products referring to the previous event.
  static std::string const PreviousTriggerInstanceName;
  
  
  // --- BEGIN -- Configuration ------------------------------------------------
  
  /// Data product with trigger hardware configuration.
  art::InputTag const fTriggerConfigTag;
  
  /// Configuration for each and all supported generators.
  std::vector<GeneratorParams_t> const fGenParams;
  
  std::string const fLogCategory; ///< Message facility stream for this tool.
  
  // --- END ---- Configuration ------------------------------------------------
  
  
  // --- BEGIN -- Current decoder ----------------------------------------------
  
  /// Pointer to the current generator and tool parameters.
  GeneratorParams_t const* fCurrentGenParam = nullptr;
  
  std::string fCurrentGenerator; ///< Name of the current generator.
  
  std::unique_ptr<IDecoder> fDecoder; ///< The decoder for the current run.
  
  // --- END ---- Current decoder ----------------------------------------------
  
  /// Sets up the object for the decoder for the specified `generator`.
  /// Throws exceptions if the setup fails.
  void createDecoder(std::string const& generator);
  
  /// Converts `DecoderConfig` into `GeneratorParams_t`.
  static std::vector<GeneratorParams_t> makeGeneratorParams
    (std::vector<Config::DecoderConfig> const& genConfigs);
  
}; // class daq::TriggerDecoderDispatched


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
std::string const daq::TriggerDecoderDispatched::CurrentTriggerInstanceName {};
std::string const daq::TriggerDecoderDispatched::PreviousTriggerInstanceName
  { "previous" };

// -----------------------------------------------------------------------------
daq::TriggerDecoderDispatched::TriggerDecoderDispatched
  (Parameters const& params)
    // configuration
  : fTriggerConfigTag{ params().TrigConfigLabel() }
  , fGenParams{ makeGeneratorParams(params().Decoders()) }
  , fLogCategory{ params().LogCategory() }
    // caches
{
  
  //
  // configuration checks
  //
  if (fGenParams.empty()) {
    throw art::Exception{ art::errors::Configuration }
      << "No tool configured for the trigger decoder ('"
      << params().Decoders.name() << "').\n";
  }
  
  //
  // configuration dump on screen
  //
  {
    mf::LogInfo log{ fLogCategory };
    log << "Trigger decoder dispatcher supports " << fGenParams.size()
      << " trigger generators:";
    for (GeneratorParams_t const& genParam: fGenParams) {
      log << "\n - '" << genParam.generatorName << "': ";
      if (genParam.toolConfig.is_key_to_atom("tool_type")) {
        log << " tool '" << genParam.toolConfig.get<std::string>("tool_type")
          << "' reading '" << genParam.fragmentTag.encode() << "'";
      }
      else {
        log << "not a tool?!";
      }
    } // for
    log << "\nTrigger generator is read from '" << fTriggerConfigTag.encode()
      << "'";
  }
  
} // daq::TriggerDecoderDispatched::TriggerDecoderDispatched()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::consumes(art::ConsumesCollector& collector)
{
  
  collector.consumes<icarus::TriggerConfiguration>(fTriggerConfigTag);
  
  // we don't ask the tools what they will/would consume, but rather have
  // a simple guess
  for (GeneratorParams_t const& genParam: fGenParams)
    collector.mayConsume<artdaq::Fragments>(genParam.fragmentTag);
  
} // daq::TriggerDecoderDispatched::consumes()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::produces(art::ProducesCollector& collector)
{
  collector.produces<std::vector<raw::ExternalTrigger>>
    (CurrentTriggerInstanceName);
  collector.produces<std::vector<raw::ExternalTrigger>>
    (PreviousTriggerInstanceName);
  collector.produces<std::vector<raw::Trigger>>
    (CurrentTriggerInstanceName);
  collector.produces<std::vector<sim::BeamGateInfo>>
    (CurrentTriggerInstanceName);
  collector.produces<sbn::ExtraTriggerInfo>
    (CurrentTriggerInstanceName);
} // daq::TriggerDecoderDispatched::produces()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::configure(const fhicl::ParameterSet&) {
  
  throw art::Exception{ art::errors::LogicError }
    << "TriggerDecoderDispatched tool does not support the `configure()` call:"
    " all configuration happened in the constructor already.\n";
  // ... delegated tools will be configured at construction time in setupRun()
  
} // daq::TriggerDecoderDispatched::configure()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::setupRun(art::Run const& run) {
  
  // determine the generator
  auto const& trigConfig
    = run.getProduct<icarus::TriggerConfiguration>(fTriggerConfigTag);
  if (fCurrentGenerator != trigConfig.generator)
    createDecoder(trigConfig.generator);

  // set up the tool for the run (delegated setup)
  fDecoder->setupRun(run);
  
} // daq::TriggerDecoderDispatched::setupRun()


// -----------------------------------------------------------------------------
std::optional<art::InputTag> daq::TriggerDecoderDispatched::preferredInput()
  const
{
  return fCurrentGenParam
    ? std::optional{ fCurrentGenParam->fragmentTag }: std::nullopt;
} // daq::TriggerDecoderDispatched::preferredInput()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::initializeDataProducts() {
  
  if (!fDecoder) {
    throw art::Exception{ art::errors::LogicError }
      << "No decoder tool configured when initializing the data product!\n";
  }
  
  fDecoder->initializeDataProducts();
  
} // daq::TriggerDecoderDispatched::initializeDataProducts()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::process_fragment
  (const artdaq::Fragment &fragment)
{
  
  if (!fDecoder) {
    throw art::Exception{ art::errors::LogicError }
      << "No decoder tool configured when processing the fragment!\n";
  }
  
  fDecoder->process_fragment(fragment);
  
} // daq::TriggerDecoderDispatched::process_fragment()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::outputDataProducts(art::Event &event) {
  
  if (!fDecoder) {
    throw art::Exception{ art::errors::LogicError }
      << "No decoder tool configured when outputting data products!\n";
  }
  
  fDecoder->outputDataProducts(event);
  
} // daq::TriggerDecoderDispatched::outputDataProducts()


// -----------------------------------------------------------------------------
void daq::TriggerDecoderDispatched::createDecoder(std::string const& generator)
{
  
  // determine the tool
  GeneratorParams_t const* newGenParam = nullptr;
  for (GeneratorParams_t const& genParam: fGenParams) {
    if (!std::regex_match(generator, genParam.generatorPattern)) continue;
    newGenParam = &genParam;
    break;
  } // for
  
  if (!newGenParam) {
    art::Exception e{ art::errors::NotFound };
    e << "Generator '" << generator << "' not supported (supported: ";
    auto itGenParam = fGenParams.cbegin();
    auto gend = fGenParams.cend();
    e << "'" << itGenParam->generatorName << "'";
    while (++itGenParam != gend) e << ", '" << itGenParam->generatorName << "'";
    e << ").\n";
    throw e;
  }
  
  // construct a new tool (only if not the same as the current one)
  fCurrentGenerator = generator;
  fCurrentGenParam = newGenParam;
  fDecoder = art::make_tool<IDecoder>(fCurrentGenParam->toolConfig);
  
  mf::LogDebug{ fLogCategory }
    << "Generator '" << fCurrentGenerator << "' matches '"
    << fCurrentGenParam->generatorName << "' tool configuration.";
  
} // daq::TriggerDecoderDispatched::createDecoder()


// -----------------------------------------------------------------------------
auto daq::TriggerDecoderDispatched::makeGeneratorParams
  (std::vector<Config::DecoderConfig> const& genConfigs)
  -> std::vector<GeneratorParams_t>
{
  std::vector<GeneratorParams_t> params;
  params.reserve(genConfigs.size());
  for (Config::DecoderConfig const& config: genConfigs) {
    params.push_back({ // C++20: use named elements
        config.Generator()                            // generatorName
      , std::regex{ config.Generator() }              // generatorPattern
      , config.FragmentsLabel()                       // fragmentTag
      , config.ToolConfig.get<fhicl::ParameterSet>()  // toolConfig
      });
  } // for
  return params;
} // daq::TriggerDecoderDispatched::makeGeneratorParams()


// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(daq::TriggerDecoderDispatched)


// -----------------------------------------------------------------------------
