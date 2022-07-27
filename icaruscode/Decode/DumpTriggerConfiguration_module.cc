/**
 * @file    DumpTriggerConfiguration_module.cc
 * @brief   Dumps to console the content of `icarus::TriggerConfiguration` data
 *          product.
 * @authors Andrea Scarpelli (ascarpell@bnl.gov),
 *          Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    June 11, 2022
 */

// SBN libraries
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <set>
#include <string>
#include <sstream> // std::ostringstream


//------------------------------------------------------------------------------
namespace sbn { class DumpTriggerConfiguration; }

/**
 * @brief Dumps on console the content of `icarus::TriggerConfiguration` data
 *        product.
 * 
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `icarus::TriggerConfiguration`: configuration of trigger in `art::Run` data
 *     product. See e.g. `icarus::TriggerConfigurationExtraction` module.
 * 
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description DumpTriggerConfiguration`.
 * 
 * * `TriggerConfigurationTag` (data product input tag): the tag identifying the
 *     data product to dump; this data product must be in `art::Run`.
 * * `Verbosity` (integral, default: maximum): verbosity level used in the
 *     dump; see `icarus::TriggerConfiguration::dump()` for details.
 * * `SkipDuplicateRuns` (flag, default: `true`): multiple files can contain
 *     information from the same run; with this flag set, only the first time
 *     a run is encountered its trigger configuration is dumped; otherwise, each
 *     time a run is opened by _art_, its configuration is printed.
 * * `OutputCategory` (string, default: `DumpTriggerConfiguration`): name of the
 *     message facility output stream to dump the information into
 * 
 */
class sbn::DumpTriggerConfiguration: public art::EDAnalyzer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> TriggerConfigurationTag {
      Name("TriggerConfigurationTag"),
      Comment("tag of trigger configuration data product (from art::Run)")
      };

    fhicl::Atom<unsigned int> Verbosity {
      Name("Verbosity"),
      Comment("verbosity level [default: maximum]"),
      icarus::TriggerConfiguration::MaxDumpVerbosity // default
      };

    fhicl::Atom<bool> SkipDuplicateRuns {
      Name("SkipDuplicateRuns"),
      Comment("print only one trigger configuration from each run"),
      true // default
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the category used for the output"),
      "DumpTriggerConfiguration"
      };

  }; // struct Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit DumpTriggerConfiguration(Parameters const& config);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Dumps the data product.
  virtual void beginRun(art::Run const& run) override;
  
  /// Does nothing, but it is mandatory.
  virtual void analyze(art::Event const& event) override {}
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Input trigger configuration tag.
  art::InputTag const fTriggerConfigurationTag;
  
  unsigned int const fVerbosity; ///< Verbosity level used for dumping.
  
  bool const fSkipDuplicateRuns; ///< Print only once from each run.

  /// Category used for message facility stream.
  std::string const fOutputCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  std::set<art::RunID> fRuns; ///< Set of runs already encountered.
  
}; // sbn::DumpTriggerConfiguration


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- sbn::DumpTriggerConfiguration
//------------------------------------------------------------------------------
sbn::DumpTriggerConfiguration::DumpTriggerConfiguration
  (Parameters const& config)
  : art::EDAnalyzer(config)
  // configuration
  , fTriggerConfigurationTag(config().TriggerConfigurationTag())
  , fVerbosity              (config().Verbosity())
  , fSkipDuplicateRuns      (config().SkipDuplicateRuns())
  , fOutputCategory         (config().OutputCategory())
{
  
  consumes<icarus::TriggerConfiguration, art::InRun>(fTriggerConfigurationTag);
  
} // sbn::DumpTriggerConfiguration::DumpTriggerConfiguration()


//------------------------------------------------------------------------------
void sbn::DumpTriggerConfiguration::beginRun(art::Run const& run) {
  
  if (fSkipDuplicateRuns) {
    art::RunID const& id = run.id();
    if (fRuns.count(id)) {
      mf::LogTrace(fOutputCategory) << id << " has already been encountered.";
      return;
    }
    fRuns.insert(id);
  } // if skip duplicates
  
  auto const& config
    = run.getProduct<icarus::TriggerConfiguration>(fTriggerConfigurationTag);
  
  std::ostringstream sstr;
  config.dump(sstr, "  ", "", fVerbosity);
  
  mf::LogVerbatim(fOutputCategory) << run.id() << ": " << sstr.str();
  

} // sbn::DumpTriggerConfiguration::beginRun()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::DumpTriggerConfiguration)


//------------------------------------------------------------------------------
