/**
 * @file   WriteBeamGateInfo_module.cc
 * @brief  Writes a fixed collection of sim::BeamGateInfo.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 18, 2021
 */

// LArSoft libraries
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <memory> // std::make_unique()


//------------------------------------------------------------------------------
namespace icarus::trigger { class WriteBeamGateInfo; }

/**
 * @brief Writes a set collection of beam gates into each event.
 * 
 * This module simply injects a list of `sim::BeamGateInfo` from the
 * configuration into each event.
 * 
 * It may be used as input to modules which require to operate on beam gates.
 * 
 * 
 * Output data products
 * =====================
 *
 * * `std::vector<sim::BeamGateInfo>`: a collection of as many
 *    `sim::BeamGateInfo` as the `BeamGates` configuration parameter entries,
 *    with the content from them.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse online description of the parameters is printed by running
 * `lar --print-description WriteBeamGateInfo`.
 * 
 * * `BeamGates` (list of configuration tables, mandatory): list of beam gate
 *     specifications; one output element is produced for each entry in this
 *     list. Each entry is a configuration table in the form:
 *     * `Duration` (time, _mandatory_): the duration of the beam
 *         gate; _the time requires the unit to be explicitly specified_: use
 *         `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *         `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *         `trigger_icarus.fcl`);
 *     * `Start` (time, default: `0_us`): how long after the
 *         @ref DetectorClocksBeamGateOpening "nominal beam gate opening time"
 *         the actual beam gate opens at;
 *     * `Type` (string, default: "unknown"): type of the gate; see online
 *         description for the configuration keys representing the values of
 *         `sim::BeamType_t`.
 * * `LogCategory` (string, default: `WriteBeamGateInfo`): name of the output
 *     stream category for console messages (managed by MessageFacility
 *     library).
 *
 * 
 */
class icarus::trigger::WriteBeamGateInfo: public art::SharedProducer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    using microseconds = util::quantities::intervals::microseconds;
    
    struct GateConfig {
      
      enum class BeamType_t { // we need to translate enum into a strong type
          kUnknown = sim::kUnknown
        , kBNB     = sim::kBNB
        , kNuMI    = sim::kNuMI
      };
      
      /// Selector for `Type` parameter.
      static util::MultipleChoiceSelection<BeamType_t> const BeamTypeSelector;
      
      fhicl::Atom<microseconds> Duration {
        Name("Duration"),
        Comment("length of gate time interval")
        };

      fhicl::Atom<microseconds> Start {
        Name("Start"),
        Comment("open the beam gate this long after the nominal beam gate time"),
        microseconds{ 0.0 }
        };

      fhicl::Atom<std::string> Type {
        Name("Type"),
        Comment("gate type: " + BeamTypeSelector.optionListString()),
        BeamTypeSelector.get(BeamType_t::kUnknown).name()
        };
      
      
      sim::BeamType_t getBeamType() const
        {
          try {
            return static_cast<sim::BeamType_t>
              (BeamTypeSelector.parse(Type()).value());
          }
          catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e)
          {
            throw art::Exception(art::errors::Configuration)
              << "Invalid value for '" << Type.name()
              << "' parameter: '" << e.label() << "'; valid options: "
              << BeamTypeSelector.optionListString() << ".\n";
          }
        } // getBeamType()
    
    };
    
    fhicl::Sequence<fhicl::TableAs<sim::BeamGateInfo, GateConfig>> BeamGates {
      Name("BeamGates"),
      Comment("list of gates to write")
      };


    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "WriteBeamGateInfo" // default
      };
    
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit WriteBeamGateInfo
    (Parameters const& config, art::ProcessingFrame const&);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Fills the plots. Also extracts the information to fill them with.
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  std::vector<sim::BeamGateInfo> const fBeamGates; ///< The gates to write.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
}; // class icarus::trigger::WriteBeamGateInfo


namespace icarus::trigger {
  
  // configuration conversion automatically picked by fhicl::TableAs<>
  sim::BeamGateInfo convert(WriteBeamGateInfo::Config::GateConfig const& config)
  {
    // nanoseconds is the standard unit for simulation and for sim::BeamGateInfo
    using nanoseconds = util::quantities::intervals::nanoseconds;
    return sim::BeamGateInfo{
        config.Start().convertInto<nanoseconds>().value()     // start
      , config.Duration().convertInto<nanoseconds>().value()  // width
      , config.getBeamType()                                  // type
      };
    
  } // convert(GateConfig)
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  util::MultipleChoiceSelection<WriteBeamGateInfo::Config::GateConfig::BeamType_t> const
  WriteBeamGateInfo::Config::GateConfig::BeamTypeSelector
    {
        { BeamType_t::kUnknown, "unknown" }
      , { BeamType_t::kBNB,     "BNB" }
      , { BeamType_t::kNuMI,    "NuMI" }
    };
  
} // namespace icarus::trigger

//------------------------------------------------------------------------------
icarus::trigger::WriteBeamGateInfo::WriteBeamGateInfo
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedProducer(config)
  // configuration
  , fBeamGates(config().BeamGates())
  , fLogCategory(config().LogCategory())
{
  
  //
  // output data declaration
  //
  produces<std::vector<sim::BeamGateInfo>>();
  
  async<art::InEvent>();
  
  //
  // configuration report (short)
  //
  {
    mf::LogInfo log { fLogCategory };
    log << "Writing " << fBeamGates.size() << " gates in each event:";
    for (sim::BeamGateInfo const& gate: fBeamGates) {
      auto const& beamType = Config::GateConfig::BeamTypeSelector.get
        (Config::GateConfig::BeamType_t(gate.BeamType()));
      log << "\n *  [ " << gate.Start() << " -- " << gate.Start() + gate.Width()
        << " ] ns (duration: " << gate.Width() << " ns), type: "
        << beamType.name() << " (" << gate.BeamType() << ")"
        ;
    } // for
  }
  
} // icarus::trigger::WriteBeamGateInfo::WriteBeamGateInfo()


//------------------------------------------------------------------------------
void icarus::trigger::WriteBeamGateInfo::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  // put a copy
  event.put(std::make_unique<std::vector<sim::BeamGateInfo>>(fBeamGates));
  
} // icarus::trigger::WriteBeamGateInfo::produce()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::WriteBeamGateInfo)


//------------------------------------------------------------------------------
