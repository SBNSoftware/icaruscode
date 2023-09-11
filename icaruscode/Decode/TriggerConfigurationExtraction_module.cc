/**
 * @file   TriggerConfigurationExtraction_module.cc
 * @brief  Writes ICARUS configuration from FHiCL into data product.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 * @date   March 23, 2022
 */

// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/TriggerConfigurationExtractor.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"

// framework libraries
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus { class TriggerConfigurationExtraction; }
/**
 * @brief Writes trigger configuration from FHiCL into a data product.
 * 
 * This module reads the configuration related to the Trigger from the FHiCL
 * configuration of the input runs and puts it into each run as a data product.
 * 
 * Input
 * ------
 * 
 * This module requires any input with _art_ run objects in it.
 * The format expected for that configuration is defined in
 * `icarus::TriggerConfigurationExtractor`, which is the utility used for the actual
 * extraction.
 * 
 * 
 * Output
 * -------
 * 
 * A data product of type `icarus::TriggerConfiguration` is placed in each run.
 * Note that the module itself does not enforce any coherence in the
 * configuration.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * The following configuration parameters are supported:
 * 
 * * **TriggerFragmentType** (string, default: `ICARUSTrigger.*`):
 *     name of the type of trigger fragment, used to identify the configuration
 *     of the trigger.
 * * **Verbose** (flag, default: `false`): if set to `true`, it will print in
 *     full the configuration of the trigger the first time it is read and each time
 *     a different one is found.
 * * **LogCategory** (string, default: `ICARUSConfigurationExtraction`):
 *     category tag used for messages to message facility.
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * This module does not support multithreading, and _art_ does not provide
 * multithreading for its functionality anyway: the only action is performed
 * at run and input file level, and the only concurrency in _art_ is currently
 * (_art_ 3.7) at event level.
 * 
 */
class icarus::TriggerConfigurationExtraction: public art::EDProducer {
  
  /// Current Trigger configuration (may be still unassigned).
  std::optional<icarus::TriggerConfiguration> fTriggerConfig;
  
  /// Whether trigger configuration inconsistency is fatal.
  bool fRequireConsistency = true;
  
  std::string fTriggerFragmentType; ///< Name of the trigger fragment type.
  
  bool fVerbose = false; ///< Whether to print the configuration we read.
  
  std::string fLogCategory; ///< Category tag for messages.
  
    public:
  
  /// Configuration of the module.
  struct Config {
    
    fhicl::Atom<std::string> TriggerFragmentType {
      fhicl::Name("TriggerFragmentType"),
      fhicl::Comment("the name of the trigger generator in DAQ"),
      "ICARUSTrigger.*" // default
      };
    
    fhicl::Atom<bool> Verbose {
      fhicl::Name("Verbose"),
      fhicl::Comment("print in full each new Trigger configuration read"),
      false // default
      };
    
    fhicl::Atom<std::string> LogCategory {
      fhicl::Name("LogCategory"),
      fhicl::Comment("category tag used for messages to message facility"),
      "TriggerConfigurationExtraction" // default
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  
  /// Constructor: just reads the configuration.
  TriggerConfigurationExtraction(Parameters const& config);
  
  
  /// Action on new run: configuration is written.
  virtual void beginRun(art::Run& run) override;
  
  /// Mandatory method, unused.
  virtual void produce(art::Event&) override {}
  
  
    private:
  
  /// Throws an exception if the `newConfig` is not compatible with the current.
  bool checkConsistency
    (icarus::TriggerConfiguration const& newConfig, std::string const& srcName) const;
  
  
}; // icarus::TriggerConfigurationExtraction


// -----------------------------------------------------------------------------
icarus::TriggerConfigurationExtraction::TriggerConfigurationExtraction
  (Parameters const& config)
  : art::EDProducer(config)
  , fTriggerFragmentType(config().TriggerFragmentType())
  , fVerbose(config().Verbose())
  , fLogCategory(config().LogCategory())
{
  
  // no consummation here (except for FHiCL configuration)
  
  produces<icarus::TriggerConfiguration, art::InRun>();
  
} // icarus::TriggerConfigurationExtraction::ICARUSTriggerConfigurationExtraction()


// -----------------------------------------------------------------------------
void icarus::TriggerConfigurationExtraction::beginRun(art::Run& run) {
  
  icarus::TriggerConfiguration config = extractTriggerReadoutConfiguration
    (run, icarus::TriggerConfigurationExtractor{ fTriggerFragmentType });
  
  checkConsistency(config, "run " + std::to_string(run.run()));
  if (!fTriggerConfig.has_value() && fVerbose)
    mf::LogInfo(fLogCategory) << "Trigger readout:" << config;
  
  fTriggerConfig = std::move(config);
  
  // put a copy of the current configuration
  run.put(std::make_unique<icarus::TriggerConfiguration>(fTriggerConfig.value()), art::fullRun());
  
} // icarus::riggerConfigurationExtraction::beginRun()


// -----------------------------------------------------------------------------
bool icarus::TriggerConfigurationExtraction::checkConsistency
  (icarus::TriggerConfiguration const& config, std::string const& srcName) const
{
  if (!fTriggerConfig.has_value() || (*fTriggerConfig == config)) return true;
  
  // see debug information for more details
  if (fRequireConsistency) {
    throw cet::exception("icarus::TriggerConfigurationExtraction")
      << "Configuration from " << srcName
      << " is incompatible with the previously found configuration.\n"
      ;
  } // if consistency is required
  
  mf::LogWarning(fLogCategory)
    << "Configuration from " << srcName
    << " is incompatible with the previously found configuration.\n"
    ;
  if (fVerbose) mf::LogVerbatim(fLogCategory) << "Trigger readout:" << config;
  return false;
  
  
} // icarus::TriggerConfigurationExtraction::checkConsistency()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::TriggerConfigurationExtraction)


// -----------------------------------------------------------------------------
