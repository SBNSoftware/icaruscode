/**
 *  @file   icaruscode/Decode/DecoderTools/PMTDecoder_tool.cc
 *
 *  @brief  This tool provides "standard" 3D hits built (by this tool) from 2D hits
 * 
 *  @author Andrea Scarpelli (ascarpell@bnl.gov)
 *
 */

// Framework Includes
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
// #include "cetlib/cpu_timer.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/intervals_fhicl.h" // for nanoseconds in FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/ExternalTrigger.h"

#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

#include "icaruscode/Decode/DecoderTools/details/PMTDecoderUtils.h"
#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()
#include "icarusalg/Utilities/BinaryDumpUtils.h" // icarus::ns::util::bin()
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration

// ROOT libraries
#include "TTree.h"

// std includes
#include <ostream>
#include <algorithm> // std::lower_bound()
#include <string>
#include <vector>
#include <tuple>
#include <optional>
#include <memory>
#include <cmath> // std::floor()
#include <cassert>


//------------------------------------------------------------------------------------------------------------------------------------------
using namespace util::quantities::time_literals;

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq { class PMTDecoder; }
/**
 * @brief Turns PMT readout fragments from DAQ into LArSoft data products.
 * 
 * The tool can read fragments from CAEN V1730 readout boards delivered by
 * artDAQ.
 * 
 * This decoder must support both a off-line mode (for storage and downstream
 * processing) and a on-line mode (for monitoring).
 * In particular, the on-line workflow is such that it may not be possible to
 * access the FHiCL configuration of the job and therefore the PMT configuration
 * data (see `icarus::PMTconfigurationExtraction` module).
 * 
 * 
 * Configuration
 * --------------
 * 
 * The set of supported parameters can be seen on command line by running
 * `lar --print-description PMTDecoder`.
 * 
 * Description of the configuration parameters:
 * * `DiagnosticOutput` (flag, default: `false`): enables additional console
 *     output, including dumping of the fragments (that is huge output).
 * * `PMTconfigTag` (data product tag, optional): if specified, the pre-trigger
 *     buffer duration is read from there; although optional, it is strongly
 *     recommended that this information be provided, since it is essential for
 *     the correct timing of the PMT waveforms (see
 *     @ref icarus_PMTDecoder_timestamps "the discussion on time stamps below").
 * * `BoardSetup` (list of board setup information): each entry specifies some
 *     information about a specific readout board; the boards are identified by
 *     their name; if a board is found in input that has no setup information,
 *     some time corrections are not applied (see
 *     @ref icarus_PMTDecoder_timestamps "the discussion on time stamps below").
 *     Each entry is in the form of a table:
 *     * `Name` (string, mandatory): the name of the board
 *        (e.g. `"icaruspmtwwtop01"`); this is use to match the setup
 *        information to a fragment ID in the PMT configuration.
 *     * `FragmentID` (integral, optional): if specified, allows the corrections
 *       using setup information to be applied even when no PMT configuration is
 *       provided (if neither PMT configuration nor setup information including
 *       `FragmentID` is available, no time correction is applied).
 *     * `TriggerDelay` (nanoseconds, default: 0 ns): measured delay from the
 *       primitive trigger time to the execution of the PMT trigger; specify
 *       the unit! (e.g. `"43 ns"`).
 * * `RequireKnownBoards` (flag, default: `true`): if set, the readout boards
 *     in input must each have a setup configuration (`BoardSetup`) *and* must
 *     be present in the PMT DAQ configuration, or an exception is thrown;
 *     if not set, readout boards in input will be processed even when their
 *     setup or DAQ configuration is not known, at the cost of fewer corrections
 *     on the timestamps (which should then be considered unreliable).
 * * `RequireBoardConfig` (flag, default: `true`): if set, the readout boards
 *     which have a setup (`BoardSetup`) are required to be included in the DAQ
 *     configuration of the input file, or an exception is thrown; if not set,
 *     missing readout boards are unnoticed.
 * * `TriggerTag` (data product tag, mandatory): tag for the information
 *     (currently required to be a collection of `raw::ExternalTrigger`,
 *     in the future it should become `raw::Trigger`);
 * * `DataTrees` (list of strings, default: none): list of data trees to be
 *     produced; if none (default), then `TFileService` is not required.
 * * `LogCategory` (string, default: `PMTDecoder`): name of the message facility
 *     category where the output is sent.
 * 
 * 
 * Requirements
 * -------------
 * 
 * Services required include:
 * 
 * * `IICARUSChannelMap` for the association of fragments to LArSoft channel ID;
 * * `DetectorClocksService` for the correct decoding of the time stamps
 *   (always required, even when dumbed-down timestamp decoding is requested);
 * * `TFileService` only if the production of trees or plots is requested.
 * 
 * 
 * Waveform time stamp
 * --------------------
 * 
 * @anchor icarus_PMTDecoder_timestamps
 * 
 * All waveforms on the same readout board share the same timestamp.
 * 
 * The time stamp of the waveform is defined as the time when the first sample
 * of the waveform started (that is, if the sample represent the value of the
 * signal in an interval of 2 ns, the time stamp is pointing at the beginning
 * of those 2 ns). Whether we can honour that definition, though, is a different
 * matter.
 * The representation of the time stamp is in the
 * @ref DetectorClocksElectronicsTime "electronics time scale".
 * 
 * There are two "types" of waveforms: the ones acquired at global trigger time,
 * and the ones acquired because of a "trigger primitive" which did not upgrade
 * to global (likely because not in coincidence with the beam gate).
 * In both cases, it is the same type of signal, a trigger primitive from
 * the NI7820 FPGA, which initializes the acquisition of the waveform.
 * Every delay between when that signal is emitted and when the PMT trigger is
 * executed shifts the time stamp of the waveform backward.
 * 
 * We assign the the time stamp of the waveforms matching the global trigger
 * as follow:
 * * the base time is the global trigger time; this effectively defines the
 *   electronics time scale, so its representation is a fixed number that is
 *   configured in LArSoft and can be accessed with
 *   `DetectorClocksData::TriggerTime()`;
 * * the delay of the propagation from the trigger board to the readout board
 *   is subtracted to the timestamp; this value must be independently measured
 *   and provided to this decoder via tool configuration as setup information
 *   (`TriggerDelay`); if not present in the setup, this delay is not considred;
 * * upon receiving the trigger, the readout board will keep some of the samples
 *   already digitized, in what we call pre-trigger buffer; the size of this
 *   buffer is a fixed number of samples which is specified in DAQ as a fraction
 *   of the complete buffer that is _post-trigger_; this amount, converted in
 *   time, is subtracted to the trigger time to point back to the beginning of
 *   the waveform instead that to the trigger primitive time. The necessary
 *   information is read from the PMT configuration (`PMTconfigTag`); if no
 *   configuration is available, this offset is not subtracted; note that this
 *   is a major shift (typically, a few microseconds) that should be always
 *   included.
 * 
 * We do not assign the the time stamp of the waveforms not matching the global
 * trigger because we have no clue how to do that. Ok, that is a to-do!
 * 
 * Each V1730 event record includes a trigger time tag (TTT), which is the value
 * of an internal counter of the board at the time the board received a trigger.
 * This can be used to relate the various waveforms (and the various fragments)
 * in the _art_ event.
 * 
 * 
 * Data trees
 * -----------
 * 
 * The tool supports the following ROOT trees production on demand:
 * 
 * * `PMTfragments`: data pertaining a single fragment; each entry is about a
 *   single fragment, and it includes full event ID, event time stamp (from
 *   _art_, i.e. the one assigned to the event by artDAQ), the ID of the
 *   fragment the entry is describing, and then the rest of the data, including:
 *     * `TTT` (64-bit integer): (Extended) Trigger Time Tag value, in readout
 *       board ticks (each worth 8 ns) from the last board reset;
 *       currently the value includes the 31-bit counter and in addition the
 *       overflow bit as MSB; the overflow bit is set by the readout board
 *       the first time the counter passes its limit (2^31) and wraps, and never
 *       cleared until the next board reset.
 *       While the tree has room for the 48-bit version (ETTT), the rest of the
 *       decoder does not yet support it.
 *     * `trigger` (64-bit signed integer): global trigger time (from the
 *       trigger decoder), in nanoseconds from The Epoch.
 *     * `triggerSec`, `triggerNS` (32-bit integer each): same time as `trigger`
 *       branch, split into second and nanosecond components.
 *     * `fragTime` (64-bit signed integer), `fragTimeSec` (32-bit signed
 *       integer): the timestamp of the PMT fragment, assigned by the board
 *       reader constructing the fragment.
 * 
 * 
 * 
 * Technical notes
 * ----------------
 * 
 * In order to correctly reconstruct the time stamp, this tool needs several
 * pieces of information.
 * These include the size of the pre-trigger buffer, which is set by the readout
 * board configuration, and the delay between the global trigger and the time
 * that trigger is received and acted upon in the readout board, which needs to
 * be measured.
 * The first category of information, from readout board configuration, are read
 * from the input file (`sbn::PMTconfiguration`), while the second category 
 * needs to be specified in the tool FHiCL configuration.
 * 
 * PMT configuration is optional, in the sense that it can be omitted; in that
 * case, some standard values will be used for it.
 * For a board to be served, an entry of that board must be present in the
 * tool configuration (`BoardSetup`). It is an error for a fragment in input not
 * to have an entry for the corresponding board setup.
 * 
 * The tool code extract the needed information and matches it into a
 * sort-of-database keyed by fragment ID, so that it can be quickly applied
 * when decoding a fragment. The matching is performed by board name.
 * 
 * 
 * Glossary
 * ---------
 * 
 * * **setup**, **[PMT] configuration**: this is jargon specific to this tool.
 *     Information about a readout board can come from two sources: the "setup"
 *     is information included in the `BoardSetup` configuration list of this
 *     tool; the "PMT configuration" is information included in the DAQ
 *     configuration that is delivered via `PMTconfigTag`.
 * * **TTT**: trigger time tag, from the V1730 event record (31 bits); may be:
 * * **ETTT**: extended trigger time tag, from the V1730 event record (48 bits).
 * * **trigger delay**: time point when a V1730 board processes a (PMT) trigger
 *     signal (and increments the TTT register) with respect to the time of the
 *     time stamp of the (SPEXi) global trigger that acquired the event.
 * 
 */
class daq::PMTDecoder: public IDecoder
{
  
  /// Enumerate the supported data trees.
  enum class DataTrees: std::size_t {
    Fragments, ///< Information about fragments
    N          ///< Counter.
  };
  using TreeNameList_t
    = std::array<std::string, static_cast<std::size_t>(DataTrees::N)>;
  static TreeNameList_t const TreeNames;
  
  /// Returns a string with all supported tree names.
  static std::string listTreeNames(std::string const& sep = "\n");
  
public:
    
    using nanoseconds = util::quantities::intervals::nanoseconds; ///< Alias.
    
    /// Data structure for trigger time.
    struct SplitTimestamp_t {
        static_assert(sizeof(int) >= 4U);
        struct Split_t {
            int seconds = std::numeric_limits<int>::min(); ///< The ongoing second.
            /// Nanoseconds from the start of current second.
            unsigned int nanoseconds = std::numeric_limits<unsigned int>::max();
        }; // Split_t
        
        /// Trigger time in nanoseconds from The Epoch.
        long long int time = std::numeric_limits<long long int>::min();
        /// Trigger time in nanoseconds from The Epoch (in components).
        Split_t split;
        
        constexpr SplitTimestamp_t() = default;
        constexpr SplitTimestamp_t(int sec, unsigned int ns);
        constexpr SplitTimestamp_t(long long int triggerTime);
    }; // SplitTimestamp_t
    
    
    // --- BEGIN -- FHiCL configuration ----------------------------------------
    
    /// Configuration of the V1730 readout board setup.
    struct BoardSetupConfig {
      
      fhicl::Atom<std::string> Name {
        fhicl::Name("Name"),
        fhicl::Comment("board name, as specified in the DAQ configuration")
        };
      
      fhicl::OptionalAtom<unsigned int> FragmentID {
        fhicl::Name("FragmentID"),
        fhicl::Comment("ID of the fragments associated with the board")
        };
      
      fhicl::Atom<nanoseconds> TriggerDelay {
        fhicl::Name("TriggerDelay"),
        fhicl::Comment
          ("from delay from the trigger timestamp to the PMT trigger [ns]"),
        0_ns
        };
      
    }; // struct BoardSetupConfig
    
    
    /// Main tool configuration.
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<bool> DiagnosticOutput {
        Name("DiagnosticOutput"),
        Comment("enable additional console output"),
        false // default
        };
      
      fhicl::Atom<bool> RequireKnownBoards {
        Name("RequireKnownBoards"),
        Comment("all readout boards in input must be known (setup+PMT configuration)"),
        true
        };
      
      fhicl::Atom<bool> RequireBoardConfig {
        Name("RequireBoardConfig"),
        Comment("all readout boards in setup must have a matching PMT configuration"),
        true
        };
      
      fhicl::OptionalAtom<art::InputTag> PMTconfigTag {
        Name("PMTconfigTag"),
        Comment("input tag for the PMT readout board configuration information")
        };
      
      fhicl::Sequence<fhicl::TableAs<details::BoardSetup_t, BoardSetupConfig>> BoardSetup {
        Name("BoardSetup"),
        Comment("list of the setup settings for all relevant V1730 boards")
        };
      
      fhicl::Atom<art::InputTag> TriggerTag {
        Name("TriggerTag"),
        Comment("input tag for the global trigger object (raw::ExternalTrigger)")
        };
      
      fhicl::Sequence<std::string> DataTrees {
        fhicl::Name("DataTrees"),
        fhicl::Comment("produces the specified ROOT trees (" + listTreeNames(",") + ")"),
        std::vector<std::string>{} // default
        };
      
    fhicl::Atom<std::string> LogCategory {
        Name("LogCategory"),
        Comment("name of the category for message stream"),
        "PMTDecoder" // default
        };
      
    }; // Config
    
    using Parameters = art::ToolConfigTable<Config>;
    
    // --- END ---- FHiCL configuration ----------------------------------------
  
  
    /**
     *  @brief  Constructor
     *
     *  @param  params configuration parameter set
     */
    explicit PMTDecoder(Parameters const& params);

    /// I hereby declare I will consume trigger and PMT configuration products.
    virtual void consumes(art::ConsumesCollector& consumerColl) override;
    
    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::ProducesCollector&) override;

    /// Reconfiguration is not supported: all configuration at construction time.
    virtual void configure(const fhicl::ParameterSet&) override;

    /// Reads the PMT configuration from the run.
    virtual void setupRun(art::Run const& run);

    /// Will read trigger information one day if needed.
    virtual void setupEvent(art::Event const& event);
    
    /**
     *  @brief Initialize any data products the tool will output
     *
     */
    virtual void initializeDataProducts() override;

    
    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param fragment            The artdaq fragment to process
     *  @param rawDigitColllection The output RawDigits
     */
    virtual void process_fragment(const artdaq::Fragment &fragment) override;

    /**
     *  @brief Output the data products to the event store
     * 
     *  @param event The event store objects
     */
    virtual void outputDataProducts(art::Event& event) override;

private:

    /// Information used in decoding from a board.
    struct NeededBoardInfo_t {
      std::string const name;
      nanoseconds preTriggerTime;
      nanoseconds PMTtriggerDelay;
    };
    
    // --- BEGIN -- Configuration parameters -----------------------------------
    bool const        fDiagnosticOutput; ///< If true will spew endless messages to output.
    
    bool const        fRequireKnownBoards; ///< Whether info on all input boards is required.
    bool const        fRequireBoardConfig; ///< Whether setup info on all boards is required.
    
    std::optional<art::InputTag> const fPMTconfigTag; ///< Input tag of the PMT configuration.
    
    art::InputTag const fTriggerTag; ///< Input tag of the global trigger.
    
    std::vector<details::BoardSetup_t> const fBoardSetup; ///< All board setup settings.
    
    std::string const fLogCategory; ///< Message facility category.
    
    // --- END ---- Configuration parameters -----------------------------------

    
    // --- BEGIN -- Services ---------------------------------------------------
    
    geo::GeometryCore const&           fGeometry; ///< Geometry service provider.

    /// Interface to LArSoft configuration for detector timing.
    detinfo::DetectorTimings const     fDetTimings;
    
    icarusDB::IICARUSChannelMap const& fChannelMap; ///< Fragment/channel mapping database.

    // --- END ---- Services ---------------------------------------------------

    
    // --- BEGIN -- Cached values ----------------------------------------------
    
    /// Duration of the optical detector readout sampling tick (i.e. 2 ns; shh!).
    nanoseconds const fOpticalTick;
    
    /// Trigger time as reported by `DetectorClocks` service.
    detinfo::timescales::electronics_time const fNominalTriggerTime;
    
    // --- END ---- Cached values ----------------------------------------------

    
    using OpDetWaveformCollection    = std::vector<raw::OpDetWaveform>;
    using OpDetWaveformCollectionPtr = std::unique_ptr<OpDetWaveformCollection>;

    OpDetWaveformCollectionPtr         fOpDetWaveformCollection;  ///< The output data collection pointer
    
    /// Find the information on a readout boards by fragment ID.
    std::optional<details::BoardInfoLookup> fBoardInfoLookup;
    
    
    /// Trigger timestamp for the current event.
    SplitTimestamp_t fTriggerTimestamp;
    
    
    /// Returns whether PMT configuration information is expected to be available.
    bool hasPMTconfiguration() const { return fPMTconfigTag.has_value(); }
    
    /// Updates the PMT configuration cache. How? Dunno. Placeholder.
    bool UpdatePMTConfiguration(sbn::PMTconfiguration const* PMTconfig);
    
    
    /**
     * @brief Returns a lookup object with board setup and configuration info.
     * @param PMTconfig the PMT configuration, if available
     * @return an object working like lookup table for all fragment information
     * 
     * This method merges the setup information from the tool configuration with
     * the PMT configuration specified in the argument, and returns an object
     * that can look up all the information as a single record, with the
     * fragment ID as key. In addition, a few intermediate quantities ("facts",
     * see `BoardFacts_t`) are computed and stored in this object.
     * 
     * If a fragment ID is missing, it means that no PMT configuration was
     * provided and that the setup information did not include a fragment ID.
     * If some information (configuration or setup) is missing, the "facts"
     * depending on the missing information will have default values.
     */
    details::BoardInfoLookup matchBoardConfigurationAndSetup
        (sbn::PMTconfiguration const* PMTconfig) const;
    
    /// Puts together all the needed information for a board.
    NeededBoardInfo_t fetchNeededBoardInfo
      (details::BoardInfoLookup::BoardInfo_t const* boardInfo, unsigned int fragmentID) const;


    /**
     * @brief Returns the count of set bits for each set bit.
     * @tparam NBits the number of bits to test
     * @param value the bit mask to be analyzed
     * @return a pair: `first` an array with the `NBits` values,
     *         `second` the number of bits set
     * 
     * Much better to go with an example.
     * Let's process `setBitIndices<16U>(0x6701)`: `value` is `0110 0111 0000 0001`,
     * and we request the `16` least significant bits (`NBits`).
     * There are six bit that are set, in position 0, 8, 9, 10, 13 and 14.
     * The returned value includes the number of set bits (6) as the `second`
     * element of the pair, an as `first` an array of 16 indices, one for each bit
     * position.
     * Its values are `0` at `[0]`, `1` at `[8]`, `2` at `[9]`, `3` at `[10]`,
     * `4` at `[13]` and `5` at `[14]`. That is, each value represents the number
     * of that bit in the list of set bits. All the other indices, associated to
     * the bits which are not set, are assigned a value that is equal or larger than
     * the number of bits set (i.e. `6` or larger).
     */
    template <std::size_t NBits, typename T>
    constexpr std::pair<std::array<std::size_t, NBits>, std::size_t>
      setBitIndices(T value) noexcept;

    
    // --- BEGIN -- Trees and their data ---------------------------------------
    
    /// Data structure for basic event information in simple ROOT trees.
    struct TreeData_EventID_t {
        unsigned int run;    ///< Run number.
        unsigned int subrun; ///< Subrun number.
        unsigned int event;  ///< Event number.
        SplitTimestamp_t timestamp; ///< Event timestamp (seconds from the epoch).
    }; // TreeData_EventID_t
    
    /// Structure collecting all data for a fragment ROOT tree.
    struct TreeFragment_t {
        struct Data_t: public TreeData_EventID_t {
            
            int fragmentID = 0; ///< ID of the fragment of this entry.
            unsigned long int TriggerTimeTag = 0; ///< Trigger time tag from the fragment.
            
            SplitTimestamp_t trigger; ///< Global trigger time.
            
            SplitTimestamp_t fragTime; ///< PMT fragment time stamp.
            
        }; // Data_t
        
        Data_t data;
        TTree* tree = nullptr;
    }; // TreeFragment_t
    
    
    std::unique_ptr<TreeData_EventID_t> fEventInfo; ///< Event ID for trees.
    std::unique_ptr<TreeFragment_t> fTreeFragment; ///< Tree with fragment information.
    
    // --- END ---- Trees and their data ---------------------------------------
    
    // --- BEGIN -- Tree-related methods ---------------------------------------
    
    /// Declares the use of event information.
    void usesEventInfo();
    
    /// Initializes all requested data trees.
    void initTrees(std::vector<std::string> const& treeNames);
    
    /// Initializes the event ID part of a tree.
    void initEventIDtree(TTree& tree, TreeData_EventID_t& data);
    
    /// Initializes the fragment data tree (`fTreeFragment`).
    void initFragmentsTree();

    /// Fills the base information of a tree data entry from an _art_ event.
    void fillTreeEventID
      (art::Event const& event, TreeData_EventID_t& treeData) const;

    /// Assigns the cached event information to the specified tree data.
    void assignEventInfo(TreeData_EventID_t& treeData) const;
  
    /// Static initialization.
    static TreeNameList_t initTreeNames();
    
    // --- END ---- Tree-related methods ---------------------------------------
    
}; // class daq::PMTDecoder


namespace daq {
  
    /// Special function `fhicl::TableAs` uses to convert BoardSetupConfig.
    details::BoardSetup_t convert(PMTDecoder::BoardSetupConfig const& config) {
        
        using ::util::fhicl::getOptionalValue;
        return {
            config.Name()                                     // name
          , getOptionalValue(config.FragmentID) 
              .value_or(details::BoardSetup_t::NoFragmentID)  // fragmentID
          , config.TriggerDelay()                             // triggerDelay
          };
    } // convert(BoardSetupConfig)
  
} // namespace daq


//-------------------------------------------------------------------------------------------                -----------------------------------------------
constexpr daq::PMTDecoder::SplitTimestamp_t::SplitTimestamp_t(int sec, unsigned int ns)
  : time { static_cast<long long int>(sec) * 1'000'000'000LL + ns }
  , split { sec, ns }
  {}

constexpr daq::PMTDecoder::SplitTimestamp_t::SplitTimestamp_t(long long int triggerTime)
  : time { triggerTime }
  , split {
      static_cast<int>(time / 1'000'000'000), // seconds
      static_cast<unsigned int>(time % 1'000'000'000) // nanoseconds
    }
  {}


namespace daq {
    std::ostream& operator<< (std::ostream& out, PMTDecoder::SplitTimestamp_t const& time) {
        out << time.split.seconds << '.'
          << std::setfill('0') << std::setw(9) << time.split.nanoseconds;
        return out;
    } // operator<< (std::ostream&, PMTDecoder::SplitTimestamp_t)
} // namespace daq


//------------------------------------------------------------------------------------------------------------------------------------------
daq::PMTDecoder::TreeNameList_t const daq::PMTDecoder::TreeNames
  = daq::PMTDecoder::initTreeNames();

auto daq::PMTDecoder::initTreeNames() -> TreeNameList_t {
  TreeNameList_t names;
  names[static_cast<std::size_t>(DataTrees::Fragments)] = "PMTfragments";
  return names;
} // daq::PMTDecoder::initTreeNames()

std::string daq::PMTDecoder::listTreeNames(std::string const& sep /* = " " */) {
  std::string l;
  for (std::string const& name: TreeNames) {
    if (!l.empty()) l += sep;
    l += '\'';
    l += name;
    l += '\'';
  } // for
  return l;
} // daq::PMTDecoder::listTreeNames()


//------------------------------------------------------------------------------------------------------------------------------------------
daq::PMTDecoder::PMTDecoder(Parameters const& params)
  : fDiagnosticOutput{ params().DiagnosticOutput() }
  , fRequireKnownBoards{ params().RequireKnownBoards() }
  , fRequireBoardConfig{ params().RequireBoardConfig() }
  , fPMTconfigTag{ ::util::fhicl::getOptionalValue(params().PMTconfigTag) }
  , fTriggerTag{ params().TriggerTag() }
  , fBoardSetup{ params().BoardSetup() }
  , fLogCategory{ params().LogCategory() }
  , fGeometry{ *(lar::providerFrom<geo::Geometry const>()) }
  , fDetTimings
    { art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob() }
  , fChannelMap{ *(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}) }
  , fOpticalTick{ fDetTimings.OpticalClockPeriod() }
  , fNominalTriggerTime{ fDetTimings.TriggerTime() }
{
  
  initTrees(params().DataTrees());
    
  mf::LogInfo log(fLogCategory);
  log << "Configuration:"
    << "\n * boards with setup: " << fBoardSetup.size();
  if (fPMTconfigTag)
    log << "\n * PMT configuration from '" << fPMTconfigTag->encode() << "'";
  else 
    log << "\n * PMT configuration not used (and some corrections will be skipped)";
  if (fRequireKnownBoards) {
    log << "\n * all readout boards in input must be known (from `"
      << params().BoardSetup.name() << "` or `"
      << params().PMTconfigTag.name() << "`)"
      ;
  }
  else {
    log << "\n * readout boards with no information (from neither `"
      << params().BoardSetup.name() << "` or `"
      << params().PMTconfigTag.name()
      << "`) are processed at the best we can (skipping corrections)"
      ;
  }
  if (fRequireBoardConfig) {
    log << "\n * all readout boards in `"
      << params().BoardSetup.name()
      << "` must appear in the PMT configuration from `"
      << params().PMTconfigTag.name() << "`"
      ;
  }
  else {
    log << "\n * all readout boards in `"
      << params().BoardSetup.name()
      << "` may lack a matching PMT configuration from `"
      << params().PMTconfigTag.name() << "`"
      ;
  }
  
} // daq::PMTDecoder::PMTDecoder()

//------------------------------------------------------------------------------------------------------------------------------------------


void daq::PMTDecoder::consumes(art::ConsumesCollector& consumerColl) {
    if (fPMTconfigTag)
        consumerColl.consumes<sbn::PMTconfiguration>(*fPMTconfigTag);
    consumerColl.consumes<std::vector<raw::ExternalTrigger>>(fTriggerTag);
}


void daq::PMTDecoder::produces(art::ProducesCollector& collector)
{
    collector.produces<OpDetWaveformCollection>();
}

//------------------------------------------------------------------------------------------------------------------------------------------
void daq::PMTDecoder::configure(fhicl::ParameterSet const&) {
  // Configuration all happens during construction
  throw cet::exception("PMTDecoder")
    << "This tool does not support reconfiguration.\n"; 
} // PMTDecoder::configure()


//------------------------------------------------------------------------------------------------------------------------------------------
void daq::PMTDecoder::setupRun(art::Run const& run)
{
    sbn::PMTconfiguration const* PMTconfig = fPMTconfigTag
      ? run.getPointerByLabel<sbn::PMTconfiguration>(*fPMTconfigTag): nullptr;
    
    UpdatePMTConfiguration(PMTconfig);
}


//------------------------------------------------------------------------------------------------------------------------------------------
void daq::PMTDecoder::setupEvent(art::Event const& event)
{
    //
    // global trigger
    //
    auto const& triggers
      = event.getByLabel<std::vector<raw::ExternalTrigger>>(fTriggerTag);
    if (triggers.size() != 1) {
      // if this is hit, the decoder needs some development to correctly deal
      // with input with no trigger, or more than one
      throw cet::exception("PMTDecoder")
        << "Found " << triggers.size() << " triggers from '"
        << fTriggerTag.encode() << "', can deal only with 1.\n";
    }
    fTriggerTimestamp = { triggers.front().GetTrigTime() };
    mf::LogTrace(fLogCategory)
      << "Trigger time ('" << fTriggerTag.encode() << "'): " << fTriggerTimestamp << " s";
    
    //
    // event ID
    //
      
    // if needed, fill the record with the basic information of the event
    if (fEventInfo) fillTreeEventID(event, *fEventInfo);
    
} // daq::PMTDecoder::setupEvent()


//------------------------------------------------------------------------------------------------------------------------------------------
void daq::PMTDecoder::initializeDataProducts()
{
    fOpDetWaveformCollection = OpDetWaveformCollectionPtr(new OpDetWaveformCollection);
}

void daq::PMTDecoder::process_fragment(const artdaq::Fragment &artdaqFragment)
{
    size_t const fragment_id = artdaqFragment.fragmentID();
    size_t const eff_fragment_id = artdaqFragment.fragmentID() & 0x0fff;

    // convert fragment to Nevis fragment
    sbndaq::CAENV1730Fragment         fragment(artdaqFragment);
    sbndaq::CAENV1730FragmentMetadata metafrag = *fragment.Metadata();
    sbndaq::CAENV1730Event            evt      = *fragment.Event();
    sbndaq::CAENV1730EventHeader      header   = evt.Header;

    size_t nChannelsPerBoard  = metafrag.nChannels; //fragment.nChannelsPerBoard();

    //std::cout << "Decoder:channels_per_board: " << nChannelsPerBoard << std::endl;

    uint32_t ev_size_quad_bytes         = header.eventSize;
    uint32_t evt_header_size_quad_bytes = sizeof(sbndaq::CAENV1730EventHeader)/sizeof(uint32_t);
    uint32_t data_size_double_bytes     = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
    uint32_t nSamplesPerChannel         = data_size_double_bytes/nChannelsPerBoard;
    uint16_t const enabledChannels      = header.ChannelMask();

    artdaq::Fragment::timestamp_t const fragmentTimestamp = artdaqFragment.timestamp();
    
    unsigned int const time_tag =  header.triggerTimeTag;

    size_t boardId = nChannelsPerBoard * eff_fragment_id;

    // chDataMap tells where in the buffer each digitizer channel is;
    // if nowhere, then the answer is a number no smaller than nEnabledChannels
    auto const [ chDataMap, nEnabledChannels ] = setBitIndices<16U>(enabledChannels);

    if (fDiagnosticOutput)
    {
        mf::LogVerbatim(fLogCategory)
          << "----> PMT Fragment ID: " << std::hex << fragment_id << std::dec
            << ", boardID: " << boardId
            << ", nChannelsPerBoard: " << nChannelsPerBoard
            << ", nSamplesPerChannel: " << nSamplesPerChannel
            << ", enabled (" << nEnabledChannels << "): " << icarus::ns::util::bin(enabledChannels)
          << "\n      size: " << ev_size_quad_bytes
            << ", data size: " << data_size_double_bytes
            << ", samples/channel: " << nSamplesPerChannel
            << ", trigger time tag: " << time_tag
            << ", time stamp: " << (fragmentTimestamp / 1'000'000'000UL)
              << "." << (fragmentTimestamp % 1'000'000'000UL) << " s"
          ;
    }

    const uint16_t* data_begin = reinterpret_cast<const uint16_t*>(artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));

    // Recover the information for this fragment
    if (fChannelMap.hasPMTDigitizerID(eff_fragment_id))
    {
        assert(fBoardInfoLookup);
    
        /*
         * The trigger time is always the nominal one, because that is the
         * reference time of the whole DAQ (PMT, TPC...).
         * We only need to know how sooner than the trigger the V1730 buffer
         * starts. Oh, and the delay from the global trigger time to when
         * the readout board receives and processes the trigger signal.
         */
        details::BoardInfoLookup::BoardInfo_t const* boardInfo = fBoardInfoLookup->findBoardInfo(fragment_id);
        if (!boardInfo) {
            if (fRequireKnownBoards) {
                cet::exception e("PMTDecoder");
                e << "Input fragment has ID " << fragment_id
                  << " which has no associated board information (`BoardSetup`";
                if (!hasPMTconfiguration()) e << " + `.FragmentID`";
                throw e << ").\n";
            }
        }
        else {
            assert(boardInfo->fragmentID == fragment_id);
            assert(boardInfo->setup);
        }
        
        NeededBoardInfo_t const info = fetchNeededBoardInfo(boardInfo, fragment_id);
        
        nanoseconds const preTriggerTime = info.preTriggerTime;
        nanoseconds const PMTtriggerDelay = info.PMTtriggerDelay;
        auto const timeStamp
          = fNominalTriggerTime - PMTtriggerDelay - preTriggerTime;
        mf::LogTrace(fLogCategory) << "V1730 board '" << info.name
          << "' has data starting at electronics time " << timeStamp
          << " = " << fNominalTriggerTime
          << " - " << PMTtriggerDelay << " - " << preTriggerTime
          ;

        const icarusDB::DigitizerChannelChannelIDPairVec& digitizerChannelVec
          = fChannelMap.getChannelIDPairVec(eff_fragment_id);
        
        std::optional<mf::LogVerbatim> diagOut;
        if (fDiagnosticOutput) {
          diagOut.emplace(fLogCategory);
          (*diagOut)
            << "      " << digitizerChannelVec.size() << " channels:";
        }
        // Allocate the vector outside the loop just since we'll reuse it over and over... 
        std::vector<uint16_t> wvfm(nSamplesPerChannel);

        // track what we do and what we want to
        uint16_t attemptedChannels = 0;
        
        for(auto const [ digitizerChannel, channelID ]: digitizerChannelVec)
        {
            if (diagOut)
              (*diagOut) << " " << digitizerChannel << " [=> " << channelID << "];";
            
            attemptedChannels |= (1 << digitizerChannel);
            
            std::size_t const channelPosInData = chDataMap[digitizerChannel];
            if (channelPosInData >= nEnabledChannels) {
                mf::LogTrace(fLogCategory)
                  << "Digitizer channel " << digitizerChannel
                  << " [=> " << channelID << "] skipped because not enabled.";
                continue;
            }
            
            std::size_t const ch_offset = channelPosInData * nSamplesPerChannel;

            std::copy_n(data_begin + ch_offset, nSamplesPerChannel, wvfm.begin());

            mf::LogTrace(fLogCategory)
              << "PMT channel " << channelID << " has " << wvfm.size()
              << " samples (read from entry #" << channelPosInData
              << " in fragment data) starting at electronics time " << timeStamp;
            fOpDetWaveformCollection->emplace_back(timeStamp.value(), channelID, wvfm);
        }
        if (diagOut) diagOut.reset(); // destroys and therefore prints out
        if (attemptedChannels != enabledChannels) {
            // this is mostly a warning; regularly, for example,
            // we effectively have 15 channels per board; but all 16 are enabled,
            // so one channel is not decoded at all
            mf::LogTrace log(fLogCategory);
            log << "Not all data read:";
            for (int const bit: util::counter(16U)) {
                uint16_t const mask = (1 << bit);
                bool const attempted = bool(attemptedChannels & mask);
                bool const enabled = bool(enabledChannels & mask);
                if (attempted == enabled) continue;
                if (!enabled) // and attempted
                    log << "\n  requested channel " << bit << " was not enabled";
                if (!attempted) // and enabled
                    log << "\n  data for enabled channel " << bit << " was ignored";
              } // for
        } // if request and availability did not match
        
    }
    else {
        mf::LogError(fLogCategory)
          << "*** PMT could not find channel information for fragment: "
            << fragment_id;
    }

    if (fDiagnosticOutput) {
        mf::LogVerbatim(fLogCategory)
          << "      - size of output collection: " << fOpDetWaveformCollection->size();
    }
    
    if (fTreeFragment) {
      fTreeFragment->data.fragmentID = fragment_id;
      fTreeFragment->data.TriggerTimeTag = time_tag;
      fTreeFragment->data.trigger = fTriggerTimestamp;
      fTreeFragment->data.fragTime = { static_cast<long long int>(fragmentTimestamp) };
      assignEventInfo(fTreeFragment->data);
      fTreeFragment->tree->Fill();
    } // if fTreeFragment
    
} // PMTDecoder::process_fragment()


void daq::PMTDecoder::outputDataProducts(art::Event& event)
{
    // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
    std::sort(fOpDetWaveformCollection->begin(),fOpDetWaveformCollection->end(),[](const auto& left,const auto&right){return left.ChannelNumber() < right.ChannelNumber();});

    // Now transfer ownership to the event store
    event.put(std::move(fOpDetWaveformCollection));

}

//------------------------------------------------------------------------------------------------------------------------------------------
bool daq::PMTDecoder::UpdatePMTConfiguration(sbn::PMTconfiguration const* PMTconfig) {
    
    fBoardInfoLookup.emplace(matchBoardConfigurationAndSetup(PMTconfig));
    
    mf::LogDebug(fLogCategory) << "Board information as cached:\n" << *fBoardInfoLookup;
    
    return true;
} // daq::PMTDecoder::UpdatePMTConfiguration()


//------------------------------------------------------------------------------------------------------------------------------------------
auto daq::PMTDecoder::matchBoardConfigurationAndSetup
  (sbn::PMTconfiguration const* PMTconfig) const -> details::BoardInfoLookup
{
    /*
     * We need to support the case where no PMT configuration is known
     * (that is the standard situation in the online monitor).
     * The "strategy" is that in such cases we give up the correct time stamp
     * decoding; if the setup information contains a fragment ID, it may be
     * possible to do a little better, that is to use the setup information
     * (this is not possible without knowing the fragment ID that each bit of
     * setup information pertains).
     * 
     * So the cases for a board are:
     * * setup information is not present: encountering such a board will cause
     *   an exception to be thrown (implemented elsewhere)
     * * PMT configuration and setup present: full configuration
     *     * exception thrown if setup fragment ID is present and inconsistent
     * * PMT configuration not present: a general warning is printed;
     *     * boards with setup fragment ID information: add setup information
     *       to the "database" for the board: it will be used for partial
     *       timestamp reconstruction
     *     * boards without setup fragment ID information: board will not be
     *       added into the database; no specific correction will be performed;
     *       a warning is printed for each board
     * 
     */
    
    // dictionary of board configurations (if any)
    std::vector<std::pair<std::string, sbn::V1730Configuration const*>> configByName;
    if (PMTconfig) {
        if (!PMTconfig->boards.empty())
            configByName.reserve(PMTconfig->boards.size());
        for (sbn::V1730Configuration const& boardConfig: PMTconfig->boards)
            configByName.emplace_back(boardConfig.boardName, &boardConfig);
        std::sort(configByName.begin(), configByName.end()); // sorted by board name
    } // if we have configuration
    
    
    auto findPMTconfig = [this,&configByName]
        (std::string const& name) -> sbn::V1730Configuration const*
      {
        if (!hasPMTconfiguration()) return nullptr;
        auto const* ppBoardConfig = details::binarySearch(configByName, name);
        if (!ppBoardConfig) {
            if (!fRequireBoardConfig) return nullptr;
            throw cet::exception("PMTDecoder")
              << "No DAQ configuration found for PMT readout board '"
              << name << "'\n"
              << "If this is expected, you may skip this check by setting "
              << "PMTDecoder tool configuration `RequireBoardConfig` to `false`.\n";
        }
        return ppBoardConfig->second;
      };
    
    // the filling is driven by boards configured in the tool
    // (which is how a setup entry is mandatory)
    details::BoardInfoLookup::Database_t boardInfoByFragment;
    
    for (details::BoardSetup_t const& boardSetup: fBoardSetup) {
        
        std::string const& boardName = boardSetup.name;
        
        sbn::V1730Configuration const* pBoardConfig = findPMTconfig(boardName);
        
        if (pBoardConfig) {
            // fragment ID from configuration and setup must match if both present
            if (boardSetup.hasFragmentID() && (boardSetup.fragmentID != pBoardConfig->fragmentID)) 
            {
                throw cet::exception("PMTDecoder")
                  << "Board '" << boardName << "' has fragment ID "
                  << std::hex << pBoardConfig->fragmentID << std::dec
                  << " but it is set up as "
                  << std::hex << boardSetup.fragmentID << std::dec
                  << "!\n";
            } // if fragment ID in setup
        }
        else {
            if (boardSetup.hasFragmentID()) {
                mf::LogPrint(fLogCategory)
                  << "Board '" << boardName
                  << "' has no configuration information:"
                    " some time stamp corrections will be skipped.";
                // to avoid this, make a PMT configuration available from input file
            }
            else {
                mf::LogPrint(fLogCategory)
                  << "Board '" << boardName
                  << "' can't be associated to a fragment ID: its time stamp corrections will be skipped.";
                // to avoid this, add a `BoardSetup.FragmentID` entry for it in the
                // configuration of this tool, or make a PMT configuration available
                continue; // no entry for this board at all
            }
        }
        
        unsigned int const fragmentID
          = pBoardConfig? pBoardConfig->fragmentID: boardSetup.fragmentID;
        assert(fragmentID != details::BoardSetup_t::NoFragmentID);
        
        nanoseconds const preTriggerTime
          = pBoardConfig
          ? (pBoardConfig->bufferLength * (1.0f - pBoardConfig->postTriggerFrac)) * fOpticalTick
          : nanoseconds{ 0.0 }
          ;
        
        details::BoardFacts_t boardFacts {
          preTriggerTime  // ditto
          };
        
        boardInfoByFragment.push_back({
          fragmentID,               // fragmentID
          &boardSetup,              // setup
          pBoardConfig,             // config
          std::move(boardFacts)     // facts
          });
    } // for
    
    return details::BoardInfoLookup{ std::move(boardInfoByFragment) };
    
} // daq::PMTDecoder::matchBoardConfigurationAndSetup()


//------------------------------------------------------------------------------------------------------------------------------------------
auto daq::PMTDecoder::fetchNeededBoardInfo
  (details::BoardInfoLookup::BoardInfo_t const* boardInfo, unsigned int fragmentID) const
  -> NeededBoardInfo_t
{
    
    return NeededBoardInfo_t{
      // name
        ((boardInfo && boardInfo->config)? boardInfo->config->boardName: ("<ID=" + std::to_string(fragmentID)))
      // preTriggerTime
      , (boardInfo? boardInfo->facts.preTriggerTime: nanoseconds{ 0.0 })
      // PMTtriggerDelay
      , ((boardInfo && boardInfo->setup)? boardInfo->setup->triggerDelay: nanoseconds{ 0.0 })
      };
      
} // daq::PMTDecoder::fetchNeededBoardInfo()


//------------------------------------------------------------------------------------------------------------------------------------------
template <std::size_t NBits, typename T>
constexpr std::pair<std::array<std::size_t, NBits>, std::size_t>
daq::PMTDecoder::setBitIndices(T value) noexcept {
    
    std::pair<std::array<std::size_t, NBits>, std::size_t> res;
    auto& [ indices, nSetBits ] = res;
    for (std::size_t& index: indices) {
        index = (value & 1)? nSetBits++: NBits;
        value >>= 1;
    } // for
    return res;
    
} // daq::PMTDecoder::setBitIndices()


//------------------------------------------------------------------------------------------------------------------------------------------
void daq::PMTDecoder::initTrees(std::vector<std::string> const& treeNames) {
    
    auto findTree = [](std::string const& name)
      {
        return static_cast<PMTDecoder::DataTrees>(
          std::distance(PMTDecoder::TreeNames.begin(),
          std::find(PMTDecoder::TreeNames.begin(), PMTDecoder::TreeNames.end(), name))
          );
      };
    
    for (std::string const& name: treeNames) {
        switch (findTree(name)) {
            case DataTrees::Fragments: initFragmentsTree(); break;
            case DataTrees::N:
            default:
                throw cet::exception("PMTdecoder")
                  << "initTrees(): no data tree supported with name '" << name << "'.\n";
        } // switch
    } // for names
    
} // daq::PMTDecoder::initTrees()


//------------------------------------------------------------------------------------------------------------------------------------------
void daq::PMTDecoder::initEventIDtree(TTree& tree, TreeData_EventID_t& data) {
    
    usesEventInfo(); // this tree includes event information
    
    tree.Branch("run", &data.run);
    tree.Branch("subrun", &data.subrun);
    tree.Branch("event", &data.event);
    tree.Branch("timestamp", &data.timestamp.time, "timestamp/L");
    
} // daq::PMTDecoder::initEventIDtree()


void daq::PMTDecoder::initFragmentsTree() {
    
    if (fTreeFragment) return;
    
    TTree* tree = art::ServiceHandle<art::TFileService>()
      ->make<TTree>("PMTfragments", "PMT fragment data");
    
    fTreeFragment = std::make_unique<TreeFragment_t>();
    fTreeFragment->tree = tree;
    auto& data = fTreeFragment->data;
    
    initEventIDtree(*tree, data);
    
    tree->Branch("fragmentID", &data.fragmentID);
    tree->Branch("fragTime", &data.fragTime.time, "fragTime/L"); // ROOT 6.24 can't detect 64-bit
    tree->Branch("fragTimeSec", &data.fragTime.split.seconds);
    tree->Branch("TTT", &data.TriggerTimeTag, "TTT/l"); // ROOT 6.24 can't detect 64-bit
    tree->Branch("trigger", &data.trigger.time, "trigger/L"); // ROOT 6.24 can't detect 64-bit
    tree->Branch("triggerSec", &data.trigger.split.seconds);
    tree->Branch("triggerNS", &data.trigger.split.nanoseconds);
    
} // daq::PMTDecoder::initFragmentsTree()


void daq::PMTDecoder::usesEventInfo() {
  
    // the allocation of fEventInfo is the flag for event information usage
    if (!fEventInfo) fEventInfo = std::make_unique<TreeData_EventID_t>();
  
} // daq::PMTDecoder::usesEventInfo()


void daq::PMTDecoder::assignEventInfo
  (TreeData_EventID_t& treeData) const
{
    
    assert(fEventInfo);
    treeData = *fEventInfo; // nice slicing
    
} // daq::PMTDecoder::assignEventInfo()


void daq::PMTDecoder::fillTreeEventID
  (art::Event const& event, TreeData_EventID_t& treeData) const
{
    art::EventID const& id = event.id();
    treeData.run    = id.run();
    treeData.subrun = id.subRun();
    treeData.event  = id.event();
    
    art::Timestamp const& timestamp = event.time();
    treeData.timestamp
      = { static_cast<int>(timestamp.timeHigh()), timestamp.timeLow() };
    
} // daq::PMTDecoder::fillTreeEventID()


//------------------------------------------------------------------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(daq::PMTDecoder)

