/**
 * @file   SaveConfigurationIntoTFile_module.cc
 * @brief  Writes the configuration of this job into `TFileService` file.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 3, 2018
 *
 *
 *
 *
 *
 */

// framework libraries
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSetRegistry.h"

// ROOT libraries
#include "TNamed.h"

// C/C++ standard libraries
#include <set>
#include <string>


//------------------------------------------------------------------------------
/**
 * @brief Writes the _art_ configuration into the `TFileService` file.
 *
 * A `TStringObj` object is written into the directory assigned by `TFileService` 
 * to this analyzer for every selected process configuration. The `TObjString`
 * will have the same ROOT name as the process itself. In addition, the
 * following are saved:
 * * `current_process` (`TNamed`): the name of this process
 *
 * The selection of processes is controlled by the module configuration.
 * The current process is always selected.
 *
 * This module takes a simplistic approach, and prints only the first occurrence
 * of each process name. When saving configuration of processes from the input
 * files, and they contain different configurations for the same process name,
 * only one of them (the first one encountered) will be saved.
 *
 * @todo Currently the author does not know how to extract the configuration of the current process, which is therefore a dummy
 *
 *
 * Configuration parameters
 * =========================
 *
 * * **includePreviousProcesses** (boolean, default: `false`): saves also the
 *     configuration of all the processes the input file has seen; this has
 *     limitations when there are multiple input files (see above).
 *
 */
class SaveConfigurationIntoTFile: public art::EDAnalyzer {
  
    public:

  struct Config {

     using Name = fhicl::Name;
     using Comment = fhicl::Comment;
     
     fhicl::Atom<bool> includePreviousProcesses{
       Name("includePreviousProcesses"),
       Comment("Also save configuration of all previous processes in the input file"),
       false
       };
     
  }; // struct Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  /// Standard _art_ analyzer moduel constructor.
  SaveConfigurationIntoTFile(Parameters const& config);
  
  /// Writes the configuration information of the current process.
  virtual void beginJob() override;
  
  /// Writes the configuration information of the input processes.
  virtual void analyze(art::Event const& event) override;
  
  
    private:
  
  /// Whether to save the configuration of the input file.
  bool fIncludePreviousProcesses = false;
  
  // we might go by configuration IDs to store all possible configurations
  std::set<std::string> fProcessedProcesses; ///< Process names already saved.
  
  /// Writes the specified configuration into a ROOT directory.
  void saveProcessConfiguration
    (std::string const& procName, fhicl::ParameterSet const& config);
  
  /// Returns whether we have already saved configuration for `procName`.
  bool isProcessConfigurationSaved(std::string const& procName) const
    { return fProcessedProcesses.count(procName) > 0; }
  
  
  /// Extracts the configuration of the specified process from `event`.
  static fhicl::ParameterSet processConfiguration
    (art::Event const& event, std::string const& procName);

  /// Saves a string via `TFileService` (as `TNamed`).
  static void storeString(std::string const& key, std::string const& value);

  /// Returns the configuration of the current _art_ process.
  static fhicl::ParameterSet currentProcessConfiguration
    (std::string const& procName);
  
  /// Returns whether the specified one is a process configuration.
  static bool isProcessConfiguration(fhicl::ParameterSet const& pset);
  
  /// Returns whether `pset` is the configuration of `procName` process.
  static bool isProcessConfiguration
    (fhicl::ParameterSet const& pset, std::string const& procName);
  
}; // class SaveConfigurationIntoTFile


DEFINE_ART_MODULE(SaveConfigurationIntoTFile)


//------------------------------------------------------------------------------
namespace {
  
  bool has_table(fhicl::ParameterSet const& pset, std::string const& name)
    { return pset.has_key(name) && pset.is_key_to_table(name); }
  
  bool has_sequence(fhicl::ParameterSet const& pset, std::string const& name)
    { return pset.has_key(name) && pset.is_key_to_sequence(name); }
  
  bool has_atom(fhicl::ParameterSet const& pset, std::string const& name)
    { return pset.has_key(name) && pset.is_key_to_atom(name); }
  
  
} // local namespace

//------------------------------------------------------------------------------
SaveConfigurationIntoTFile::SaveConfigurationIntoTFile(Parameters const& config)
  : art::EDAnalyzer(config)
  , fIncludePreviousProcesses(config().includePreviousProcesses())
  {}


//------------------------------------------------------------------------------
void SaveConfigurationIntoTFile::beginJob() {
  
  std::string const procName = processName();
  //
  // current process
  //
  fhicl::ParameterSet pset = currentProcessConfiguration(procName);
  
  MF_LOG_DEBUG("SaveConfigurationIntoTFile")
    << "This process: '" << procName << "':\n"
    << std::string(80, '=') << '\n'
    << pset.to_indented_string()
    << '\n' << std::string(80, '=');
  saveProcessConfiguration(procName, pset);
  storeString("current_process", procName);
  
  mf::LogInfo("SaveConfigurationIntoTFile")
    << "Configuration of current process ('" << procName
    << "') stored via TFileService.";
  
} // SaveConfigurationIntoTFile::beginJob()
  

//------------------------------------------------------------------------------
void SaveConfigurationIntoTFile::analyze(art::Event const& event) {
  
  //
  // previous processes for input file
  //
  if (fIncludePreviousProcesses) { 
    for (auto const& procConfigInfo: event.processHistory()) {
      std::string const& procName = procConfigInfo.processName();
      if (isProcessConfigurationSaved(procName)) continue;
      
      saveProcessConfiguration
        (procName, processConfiguration(event, procName));
      mf::LogInfo("SaveConfigurationIntoTFile")
        << "Configuration of process '" << procName
        << "' from input file stored via TFileService.";
      
    } // for process config
  } // if include history
  
} // SaveConfigurationIntoTFile::analyze()


//------------------------------------------------------------------------------
void SaveConfigurationIntoTFile::saveProcessConfiguration(
  std::string const& procName, fhicl::ParameterSet const& config
) {
 
  if (fProcessedProcesses.count(procName) > 0) return; // already saved

  storeString(procName, config.to_indented_string());
  
  fProcessedProcesses.insert(procName);
  
} // SaveConfigurationIntoTFile::saveProcessConfiguration()


//------------------------------------------------------------------------------
fhicl::ParameterSet SaveConfigurationIntoTFile::processConfiguration
  (art::Event const& event, std::string const& procName)
{
  fhicl::ParameterSet pset;
  event.getProcessParameterSet(procName, pset);
  return pset;
} // SaveConfigurationIntoTFile::processConfiguration()


//------------------------------------------------------------------------------
void SaveConfigurationIntoTFile::storeString
  (std::string const& key, std::string const& value)
{
  /*
   * The following `TFileService::makeAndRegister()` call does two things:
   * 1) constructs a `TNamed` object with the _default constructor_
   * 2) assigns it the specified name and title
   * 3) adds the new object to the ROOT directory we are given.
   * 
   * The first two operations might be combined by calling the proper
   * constructor, but the generic `makeAndRegister()` will always explicitly
   * perform action (2) anyway, while `make()` will not perform action (3)
   * at all, which appears to be required for successfully writing `TNamed`.
   * 
   */
  auto const& fs = *(art::ServiceHandle<art::TFileService>());
  fs.makeAndRegister<TNamed>(key.c_str(), value.c_str());
  
} // SaveConfigurationIntoTFile::storeString()


//------------------------------------------------------------------------------
fhicl::ParameterSet SaveConfigurationIntoTFile::currentProcessConfiguration
  (std::string const& procName)
{

  /*
   * The strategy is to parse the global FHiCL registry (yes, there is such a
   * thing), looking for the table which contains 'process_name'.
   * This heuristic is not fool proof, since there might be other tables with
   * the same atom. So:
   * TODO find a way to ask art which is the root parameter set
   * The registry contains all root-level tables, including the main
   * configuration and other tables that are used to resolve references
   * (internally, references are not resolved to save space).
   * It turns out that when a fhicl::ParameterSet is created, all references
   * are resolved, so by the time we get the parameter sets we don't have to
   * worry to manually resolve the references. Which would not be fun.
   *
   */

  static std::string const rootKey { "process_name" };
  
  fhicl::ParameterSet const* rootConfig = nullptr;
  for (auto const& idAndSet: fhicl::ParameterSetRegistry::get()) {
    auto const& pset = idAndSet.second;
    if (isProcessConfiguration(pset, procName)) {
      if (rootConfig) {
        throw art::Exception(art::errors::LogicError)
          << "Found two candidate process parameter sets: with " << rootKey
          << " '" << rootConfig->get<std::string>(rootKey) << "' and '"
          << pset.get<std::string>(rootKey) << "'!\n";
      }
      rootConfig = &pset;
    } // if found
  }
  if (!rootConfig) {
    throw art::Exception(art::errors::LogicError)
      << "No parameter set with 'process_name' atom found!\n";
  }
  
  return *rootConfig;
} // SaveConfigurationIntoTFile::currentProcessConfiguration()


//------------------------------------------------------------------------------
bool SaveConfigurationIntoTFile::isProcessConfiguration
  (fhicl::ParameterSet const& pset)
{
  /*
   * We are forced to parse the parameter set heuristically.
   * If a configuration is something like:
   *     
   *     full_configuration: { ... }
   *
   *     @table::full_configuration
   *
   * the heuristic will match both the actual root configuration *and*
   * `full_configuration` parameter set (hint: put `full_configuration` within a
   * prolog). This may cause trouble downstream.
   * Barred that case, the more detailed the logic here is, the less likely
   * a false positive is to happen.
   *
   */
  if (!::has_table   (pset, "services"    )) return false;
  if (!::has_table   (pset, "physics"     )) return false;
  if (!::has_table   (pset, "source"      )) return false;
  if (!::has_atom    (pset, "process_name")) return false;
  if (!::has_sequence(pset, "physics.trigger_paths")) return false;
  if (!::has_sequence(pset, "physics.end_paths")) return false;
  
  return true;
} // SaveConfigurationIntoTFile::isProcessConfiguration()


//------------------------------------------------------------------------------
bool SaveConfigurationIntoTFile::isProcessConfiguration
  (fhicl::ParameterSet const& pset, std::string const& procName)
{
  if (!isProcessConfiguration(pset)) return false;
  return pset.get<std::string>("process_name") == procName;
} // SaveConfigurationIntoTFile::isProcessConfiguration(std::string)


//------------------------------------------------------------------------------

