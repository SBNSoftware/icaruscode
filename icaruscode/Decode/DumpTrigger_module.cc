/**
 * @file   DumpTrigger_module.cc
 * @brief  Dumps on console the content of trigger data products.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 22, 2022
 */

// SBN libraries
#include "icaruscode/Decode/DataProducts/ExtraTriggerInfo.h"
#include "icaruscode/Decode/BeamBits.h" // sbn::triggerSource
#include "icaruscode/Utilities/StreamIndenter.h" // util::addIndent()

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <iomanip> // std::setfill(), std::setw()
#include <vector>
#include <optional>
#include <utility> // std::pair
#include <string>


//------------------------------------------------------------------------------
namespace sbn { class DumpTrigger; }

/**
 * @brief Dumps on console the content of trigger data products.
 * 
 * Supported trigger types are:
 * * `std::vector<raw::Trigger>`
 * * `std::vector<raw::ExternalTrigger>`
 * * `std::vector<sim::BeamGateInfo>`
 * * `std::vector<sbn::ExtraTriggerInfo>`
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::Trigger>`: simple trigger information
 * * `std::vector<sim::BeamGateInfo>`: beam gate information
 * * `std::vector<raw::ExternalTrigger>`: additional trigger information
 * * `sbn::ExtraTriggerInfo`: detailed trigger information
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * All data product tags are optional, but at least one needs to be specified.
 * If a tag is not specified, the same tag as `TriggerTag` will be attempted
 * (if specified).
 * A terse description of the parameters is printed by running
 * `lar --print-description DumpTrigger`.
 * 
 * * `TriggerTag` (data product input tag, optional): the tag identifying the
 *     data product of the simple trigger information to dump.
 * * `BeamGateTag` (data product input tag, optional): the tag identifying the
 *     data product of the beam gate; if explicitly empty, this type of data
 *     product will not be dumped; if omitted, an attempt to dump a data
 *     product with the same tag as `TriggerTag` will be performed, and in case
 *     of failure no message will be printed.
 * * `ExternalTriggerTag` (data product input tag, optional): the tag
 *     identifying the data product of the additional standard trigger
 *     information; this is a standard LArSoft data product that ICARUS abuses
 *     to store information in a portable way. If explicitly empty, this type of
 *     data product will not be dumped; if omitted, an attempt to dump a data
 *     product with the same tag as `TriggerTag` will be performed, and in case
 *     of failure no message will be printed.
 * * `ExtraTriggerTag` (data product input tag, optional): the tag identifying
 *     the data product of the detailed trigger information; this is
 *     SBN-specific. If explicitly empty, this type of data product will not be
 *     dumped; if omitted, an attempt to dump a data product with the same tag
 *     as `TriggerTag` will be performed, and in case of failure no message will
 *     be printed.
 * * `Verbosity` (integral, default: maximum): verbosity level used in the
 *     dump; see `sbn::PMTconfiguration::dump()` for details.
 * * `OutputCategory` (string, default: `DumpTrigger`): name of the
 *     message facility output stream to dump the information into.
 * 
 */
class sbn::DumpTrigger: public art::EDAnalyzer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::OptionalAtom<art::InputTag> TriggerTag {
      Name{ "TriggerTag" },
      Comment
        { "tag of the simple trigger data product (raw::Trigger collection)" }
      };

    fhicl::OptionalAtom<art::InputTag> BeamGateTag {
      Name{ "BeamGateTag" },
      Comment{ "tag of the beam gate information data product" }
      };

    fhicl::OptionalAtom<art::InputTag> ExternalTriggerTag {
      Name{ "ExternalTriggerTag" },
      Comment{
        "tag of additional trigger information data product"
        " (raw::ExternalTrigger collection)" 
        }
      };

    fhicl::OptionalAtom<art::InputTag> ExtraTriggerTag {
      Name{ "ExtraTriggerTag" },
      Comment{ "tag of detailed SBN-specific trigger information data product" }
      };

    fhicl::Atom<unsigned int> Verbosity {
      Name{ "Verbosity" },
      Comment{ "verbosity level [default: maximum]" },
      std::numeric_limits<unsigned int>::max() // default
      };

    fhicl::Atom<std::string> OutputCategory {
      Name{ "OutputCategory" },
      Comment{ "name of the category used for the output" },
      "DumpTrigger"
      };

  }; // struct Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit DumpTrigger(Parameters const& config);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Does the dumping.
  virtual void analyze(art::Event const& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  using electronics_time = detinfo::timescales::electronics_time; // alias
  using simulation_time = detinfo::timescales::simulation_time; // alias
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  std::optional<art::InputTag> const fTriggerTag; ///< Input trigger tag.
  
  std::optional<art::InputTag> const fBeamGateTag; ///< Input beam gate tag.
  
  /// Input additional trigger tag.
  std::optional<art::InputTag> const fExternalTag;
  
  std::optional<art::InputTag> const fExtraTag; /// Input extra trigger tag.
  
  unsigned int const fVerbosity; ///< Verbosity level used for dumping.
  
  /// Category used for message facility stream.
  std::string const fOutputCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  /// Dumps a simple LArSoft trigger data product.
  void dumpTrigger(std::vector<raw::Trigger> const& triggers) const;
  
  /// Dumps a LArSoft external trigger data product.
  void dumpTrigger(std::vector<raw::ExternalTrigger> const& triggers) const;
  
  /// Dumps a LArSoft beam gate information data product.
  void dumpBeamGate(std::vector<sim::BeamGateInfo> const& gates) const;
  
  /// Dumps a SBN trigger information data product.
  void dumpTrigger(sbn::ExtraTriggerInfo const& trigger) const;
  
  
  /// Returns the tag to try, and whether it is mandatory to find it.
  std::pair<art::InputTag, bool> inputTag
    (std::optional<art::InputTag> const& tag) const;
  
  
  /// Returns `1` if the data product of type `T` specified by `param` must
  /// be available, `0` otherwise. It also declares it will consume it.
  template <typename T>
  unsigned int countConsume(std::optional<art::InputTag> const& param);
  
}; // sbn::DumpTrigger


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  struct dumpTimestamp {
    raw::TriggerTimeStamp_t ts;
    dumpTimestamp(raw::TriggerTimeStamp_t ts): ts{ ts } {}
  }; // dumpTimestamp
  
  std::ostream& operator<< (std::ostream& out, dumpTimestamp tsw) {
    return out << (tsw.ts / 1'000'000'000) << "." << std::setfill('0')
      << std::setw(9) << (tsw.ts % 1'000'000'000) << std::setfill(' ');
  } // operator<< (dumpTimestamp)
  
  
  
} // local namespace

//------------------------------------------------------------------------------
//--- sbn::DumpTrigger
//------------------------------------------------------------------------------
sbn::DumpTrigger::DumpTrigger
  (Parameters const& config)
  : art::EDAnalyzer(config)
  // configuration
  , fTriggerTag    { config().TriggerTag() }
  , fBeamGateTag   { config().BeamGateTag() }
  , fExternalTag   { config().ExternalTriggerTag() }
  , fExtraTag      { config().ExtraTriggerTag() }
  , fVerbosity     { config().Verbosity() }
  , fOutputCategory{ config().OutputCategory() }
{
  
  unsigned int nRequiredTriggers = 0U; // count of required data products
  
  /*
   * data product usage declaration
   */
  nRequiredTriggers += countConsume<std::vector<raw::Trigger>>(fTriggerTag);
  nRequiredTriggers
    += countConsume<std::vector<sim::BeamGateInfo>>(fBeamGateTag);
  nRequiredTriggers
    += countConsume<std::vector<raw::ExternalTrigger>>(fExternalTag);
  nRequiredTriggers += countConsume<sbn::ExtraTriggerInfo>(fExtraTag);
  
  /*
   * parameter check
   */
  if (nRequiredTriggers == 0U) {
    throw art::Exception{ art::errors::Configuration }
      << "At least one trigger data product needs to be specified.\n";
  }
  
  /*
   * print configuration
   */
  mf::LogInfo log{ fOutputCategory };
  log << "Configuration:\n * dump:";
  
  if (fTriggerTag) log << "\n - '" << fTriggerTag->encode() << "' (raw::Trigger)";
  
  if (auto const [ tag, req ] = inputTag(fExternalTag); !tag.empty()) {
    log << "\n - '" << tag.encode() << "' (raw::ExternalTrigger)";
    if (!req) log << " (if available)";
  }
  
  if (auto const [ tag, req ] = inputTag(fExtraTag); !tag.empty()) {
    log << "\n - '" << tag.encode() << "' (sbn::ExtraTriggerInfo)";
    if (!req) log << " (if available)";
  }
  
  if (auto const [ tag, req ] = inputTag(fBeamGateTag); !tag.empty()) {
    log << "\n - '" << tag.encode() << "' (sim::BeamGateInfo)";
    if (!req) log << " (if available)";
  }
  
} // sbn::DumpTrigger::DumpTrigger()


//------------------------------------------------------------------------------
void sbn::DumpTrigger::analyze(art::Event const& event) {
  
  mf::LogVerbatim{ fOutputCategory }
    << "Trigger information in " << event.id() << ":";
  
  if (auto const [ tag, req ] = inputTag(fTriggerTag); !tag.empty()) {
    
    auto const& handle = event.getHandle<std::vector<raw::Trigger>>(tag);
    if (handle) {
      mf::LogVerbatim{ fOutputCategory }
        << "* raw::Trigger ('" << tag.encode() << "'):";
      dumpTrigger(*handle);
    }
    else if (req) throw *(handle.whyFailed());
  } // simple trigger
  
  
  if (auto const [ tag, req ] = inputTag(fExternalTag); !tag.empty()) {
    
    auto const& handle
      = event.getHandle<std::vector<raw::ExternalTrigger>>(tag);
    if (handle) {
      mf::LogVerbatim{ fOutputCategory }
        << "* raw::ExternalTrigger ('" << tag.encode() << "'):";
      dumpTrigger(*handle);
    }
    else if (req) throw *(handle.whyFailed());
  } // external trigger info
  
  
  if (auto const [ tag, req ] = inputTag(fExtraTag); !tag.empty()) {
    
    auto const& handle = event.getHandle<sbn::ExtraTriggerInfo>(tag);
    if (handle) {
      mf::LogVerbatim{ fOutputCategory }
        << "* sbn::ExtraTriggerInfo ('" << tag.encode() << "'):";
      dumpTrigger(*handle);
    }
    else if (req) throw *(handle.whyFailed());
  } // extra trigger info
  
  
  if (auto const [ tag, req ] = inputTag(fBeamGateTag); !tag.empty()) {
    
    auto const& handle = event.getHandle<std::vector<sim::BeamGateInfo>>(tag);
    if (handle) {
      mf::LogVerbatim{ fOutputCategory }
        << "* sim::BeamGateInfo ('" << tag.encode() << "'):";
      dumpBeamGate(*handle);
    }
    else if (req) throw *(handle.whyFailed());
  } // beam gate
  
  
} // sbn::DumpTrigger::analyze()


//------------------------------------------------------------------------------
void sbn::DumpTrigger::dumpTrigger
  (std::vector<raw::Trigger> const& triggers) const
{
  static std::string const indent(4U, ' ');
  
  if (triggers.empty()) {
    mf::LogVerbatim{ fOutputCategory } << indent << "no triggers.";
    return;
  }
  if (triggers.size() > 1) {
    mf::LogVerbatim{ fOutputCategory }
      << indent << "[" << triggers.size() << " triggers]";
  }
  
  for (auto const& [ iTrigger, trigger]: util::enumerate(triggers)) {
    mf::LogVerbatim log { fOutputCategory };
    log << indent;
    if (triggers.size() > 1) log << "[" << iTrigger << "] ";
    log << "trigger #" << trigger.TriggerNumber() << " at "
      << electronics_time{ trigger.TriggerTime() }
      << ", beam gate at " << electronics_time{ trigger.BeamGateTime() }
      << ", bits:";
    if (sbn::bits::triggerSourceMask const bitMask{ trigger.TriggerBits() }) {
      log << " {";
      for (std::string const& name: names(bitMask)) log << " " << name;
      log << " }";
    }
    else log << " none";
  
  } // for
  
} // sbn::DumpTrigger::dumpTrigger()


//------------------------------------------------------------------------------
void sbn::DumpTrigger::dumpTrigger
  (std::vector<raw::ExternalTrigger> const& triggers) const
{
  
  static std::string const indent(4U, ' ');
  
  if (triggers.empty()) {
    mf::LogVerbatim{ fOutputCategory } << indent << "no triggers.";
    return;
  }
  if (triggers.size() > 1) {
    mf::LogVerbatim{ fOutputCategory }
      << indent << "[" << triggers.size() << " triggers]";
  }
  
  for (auto const& [ iTrigger, trigger]: util::enumerate(triggers)) {
    mf::LogVerbatim log { fOutputCategory };
    log << indent;
    if (triggers.size() > 1) log << "[" << iTrigger << "] ";
    log << "trigger ID=" << trigger.GetTrigID()
      << " at " << dumpTimestamp(trigger.GetTrigTime());
  
  } // for
  
} // sbn::DumpTrigger::dumpTrigger(ExternalTrigger)


//------------------------------------------------------------------------------
void sbn::DumpTrigger::dumpBeamGate
  (std::vector<sim::BeamGateInfo> const& gates) const
{
  
  static std::string const indent(4U, ' ');
  
  if (gates.empty()) {
    mf::LogVerbatim{ fOutputCategory } << indent << "no beam gates.";
    return;
  }
  if (gates.size() > 1) {
    mf::LogVerbatim{ fOutputCategory }
      << indent << "[" << gates.size() << " beam gates]";
  }
  
  for (auto const& [ iGate, gate]: util::enumerate(gates)) {
    mf::LogVerbatim log { fOutputCategory };
    log << indent;
    simulation_time const start { gate.Start() };
    simulation_time const end { gate.Start() + gate.Width() };
    if (gates.size() > 1) log << "[" << iGate << "] ";
    log << "beam gate [ " << start << " -- " << end << " ] (duration: "
      << (end - start) << ") of type ";
    switch (gate.BeamType()) {
      case sim::kBNB:     log << "BNB"; break;
      case sim::kNuMI:    log << "NuMI"; break;
      case sim::kUnknown: log << "unknown"; break;
      default:
        log << "unsupported [code=" << static_cast<int>(gate.BeamType()) << "]";
    } // switch
  
  } // for
  
} // sbn::DumpTrigger::dumpBeamGate()


//------------------------------------------------------------------------------
void sbn::DumpTrigger::dumpTrigger(sbn::ExtraTriggerInfo const& trigger) const {
  
  static std::string const indent(4U, ' ');
  
  mf::LogVerbatim{ fOutputCategory } << util::addIndent(indent) << trigger;
  
} // sbn::DumpTrigger::dumpTrigger(ExtraTriggerInfo)


//------------------------------------------------------------------------------
std::pair<art::InputTag, bool> sbn::DumpTrigger::inputTag
  (std::optional<art::InputTag> const& tag) const
{
  return
    { tag.value_or(fTriggerTag.value_or(art::InputTag{})), tag.has_value() };
} // sbn::DumpTrigger::inputTag()


//------------------------------------------------------------------------------
template <typename T>
unsigned int sbn::DumpTrigger::countConsume
  (std::optional<art::InputTag> const& param)
{
  auto const [ tag, required ] = inputTag(param);
  if (tag.empty()) return 0;
  if (required) consumes<T>(tag);
  else          mayConsume<T>(tag);
  return required? 1: 0;
} // sbn::DumpTrigger::countConsume()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::DumpTrigger)


//------------------------------------------------------------------------------
