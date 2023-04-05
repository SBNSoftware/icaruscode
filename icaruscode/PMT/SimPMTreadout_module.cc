/**
 * @file   icaruscode/PMT/SimPMTreadout_module.cc
 * @brief  Produces PMT waveforms reflecting ICARUS readout.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 16, 2022
 * 
 */

// ICARUS libraries
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()
#include "icaruscode/PMT/PMTpedestalGeneratorTool.h"
#include "icaruscode/PMT/Algorithms/PedestalGeneratorAlg.h"
#include "icaruscode/PMT/Algorithms/PMTReadoutWindowMaker.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond
#include "lardataalg/Utilities/quantities/electronics.h" // counts_f
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/make_tool.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C/C++ standard library
#include <vector>
#include <algorithm> // std::min()
#include <atomic> // std::atomic_flag
#include <mutex> // std::mutex, std::lock_guard
#include <memory> // std::make_unique()
#include <tuple>
#include <utility> // std::move()
#include <cmath> // std::round()
#include <limits> // std::numeric_limits


namespace icarus::opdet { class SimPMTreadout; }
/**
 * @brief Simulates the digitization of ICARUS PMT response including readout.
 * 
 * This module operates on top of already digitised PMT waveforms, with the
 * purpose of simulating the effect of the selective readout of ICARUS PMT
 * system.
 * 
 * In the detector, PMT signals are first digitised and discriminated in real
 * time, then combined to determine trigger primitives. When a trigger
 * primitive is formed, it causes the PMT readout board to actually save the
 * PMT information.
 * 
 * In the simulation, this is done in steps:
 * 
 * 1. Digitization occurs first based on scintillation photons which converted
 *    at the PMT surface; this is the task of `icarus::opdet::SimPMTIcarus`
 *    module, which creates `raw::OpDetWaveform` with all scintillation signal
 *    (except that a low threshold is used), including noise.
 * 2. Trigger primitives are determined based on the waveforms from the
 *    previous step. This task can be performed by the proper configuration of
 *    `icarus::trigger::TriggerSimulationOnGates` module.
 * 3. Finally, the trigger primitives are used to determine the time intervals
 *    and the channels when and where the readout happens, and new waveforms
 *    are created combining the samples from the waveforms in the first step.
 * 
 * 
 * Readout details
 * ================
 * 
 * The current implementation operates mimicking the ICARUS PMT and triger
 * systems, although the delays are not simulated.
 * 
 * 1. Trigger primitives are expected to be generated with a time, precise at
 *    the PMT sample (2 nanoseconds), and a location, precise at the cryostat
 *    (conceivably it could be precise at the TPC, but currently the trigger
 *    compbines the two TPC in a cryostat).
 * 2. Optical detector waveforms used as input are expected to already contain
 *    all the signal.
 * 3. Readout parameters are learned from the module configuration, as well
 *    as the noise parameters.
 * 4. Artificial trigger primitives are added to reflect the ones sent by the
 *    hardware to ensure coverage of the beam gate time
 *    (see @ref SimPMTreadout_BeamGate "Readout during beam gate" below).
 * 5. For each trigger primitive, all PMT in the primitive cryostat are read
 *    out, with a time interval extending before the trigger primitive time
 *    by the amount of configured pre-trigger buffer, and with a fixed
 *    duration.
 * 6. Contiguous buffers are merged in a single waveform.
 * 
 * 
 * Waveform content creation
 * --------------------------
 * 
 * The algorithm proceeds with the following steps:
 * 
 * 1. Read all the trigger primitives from input.
 * 2. Generate additional, artificial trigger primitives as requested in the
 *    configuration
 *    (see @ref SimPMTreadout_BeamGate "Readout during beam gate" below)
 * 3. Surround each trigger primitive with a readout buffer interval, and
 *    determine the final readout intervals for all channels and waveforms.
 *    A single channel will have multiple intervals, and each interval will
 *    generate a waveform.
 * 4. For each interval, the source waveform data is copied into the new one.
 *    If there are times that are uncovered, the noise generator is tasked
 *    to full the gaps.
 * 
 * It is assumed that the noise will not contribute to further trigger
 * primitives. This is a solid assumption, considering that the noise is added
 * only on top of the baseline (no signal present), and that the noise is two
 * orders of magnitude smaller than the discrimination threshold the trigger
 * primitives utilize.
 * 
 * 
 * ### Readout during beam gate
 * 
 * @anchor SimPMTreadout_BeamGate
 * 
 * The trigger hardware is configured to always send a sequence of artificial
 * trigger primitives at around the time the beam gate is expected.
 * For example, in ICARUS Run1 three trigger primitives were sent to all the
 * PMT readout boards (i.e. triggering the readout of all 360 channels) at
 * the times of -4, +3 and +12 microseconds from the beam gate (which is also
 * the time the trigger hardware opens the gate to form a global trigger).
 * This, combined with a readout buffer of 10 microseconds, 3 before the
 * trigger primitive signal and 7 after it, would ensure the coverage by the
 * PMT of 26 microseconds, from -7 to +19 with respect to the beam gate, and
 * therefore a coverage of 7 microseconds before any beam interaction to more
 * than 9 microseconds after the latest of NuMI neutrino possible interaction.
 * 
 * In addition, a trigger primitive is still sent at the time of the global
 * trigger; with the aforementioned configuration, though, this does not add
 * any time coverage to the PMT.
 * 
 * 
 * ### Noise generation and absolute event time
 * 
 * A noise generator is necessary to extend the waveforms to intervals where no
 * original samples are available (as a reminder: the simulation fills a narrow
 * region of waveform around any minor activity, while the detector readout
 * acquires all PMT in a cryostat even when just a few of them show activity).
 * Noise generation is delegated to an algorithm belonging to
 * `icarus::opdet::NoiseGeneratorAlg` class. This class of generators may in
 * principle generate noise specific to a channel and to a certain data
 * acquisition time, and they require such information to be provided to them.
 * We have two options to learn about the absolute time of the event: from a
 * trigger object (`sbn::ExtraTriggerInfo` is produced by the trigger emulation)
 * from the time of the event generation.
 * In general, the status as of December 2022 in ICARUS is that no noise
 * generator actually uses that information, we don't have automatically
 * accessible records of different PMT noise conditions as function of time,
 * and the Monte Carlo production is not set up to reflect any specific data
 * taking condition. As a consequence, not only the choice here is largely
 * irrelevant, but it's also not possible to make one that is "future-proof".
 * The current implementation of this module utilizes the latter approach,
 * of using the event time, because it's the simplest and least demanding;
 * the former might be made optional too (for example, extracting the timestamp
 * from the `sbn::ExtraTriggerInfo` data product matching the first of the
 * `TriggerTags` tags).
 * 
 * 
 * Configuration
 * ==============
 * 
 * The configuration is also available with the usual
 * `lar --print-description SimPMTreadout`.
 * 
 * * `WaveformTags` (list of input tags, mandatory): list of data products
 *   with the waveform samples to be arranged in the new waveforms.
 * * `PrimitiveTags` (list of input tags, mandatory): list of trigger
 *   primitive data products to seed the PMT readout with.
 * * `TriggerTags` (list of input tags, default: empty): list of global
 *   triggers contributing additional trigger primitives to the whole
 *   detector.
 * * `FixedTriggerTimes` (list of times, optional): adds a global trigger (in
 *   the same fashion as the ones from `TriggerTags` collections) to the times
 *   specified in this list. Times here are relative to the beam gate opening
 *   time.
 * * `Readout` (table, mandatory): readout settings:
 *     * `WindowSize` (integer, mandatory): the size of a readout buffer in
 *       optical ticks.
 *     * `PreTrigFraction` (real, mandatory): the fraction of the readout buffer
 *       that is collected before the trigger primitive (e.g. a fraction `0.3`
 *       over a buffer of `5000` samples will have `1500` samples before the
 *       trigger primitive, and the remaining `3500` starting from the primitive
 *       time on).
 * * `Pedestal` (tool configuration table): configuration of the tool
 *   used to generate pedestal and noise in the intervals where no source
 *   samples are available. The tool is one derived from
 *   `icarus::opdet::PMTpedestalGeneratorTool`. The recommendation is to use the
 *   same generator and parameters as used in the original simulation.
 *   even if no noise addition is desired.
 * * `NoiseGeneratorSeed` (integer, optional): if specified, the value is used
 *   to seed the noise generation algorithm random engine; otherwise, the seed
 *   is assigned by `NuRandomService`.
 * * `CryostatFirstBit` (positive integer, mandatory): the trigger bit which is
 *   set on input triggers from cryostat 0; as many bits are expected as
 *   there are cryostats in the detector (thus, 2 for ICARUS). The value `0`
 *   represents the least significant bit.
 * * `LogCategory` (string, default: `SimPMTreadout`): message facility
 *   category to which the module messages are routed.
 * 
 * 
 * Input
 * ======
 * 
 * * `std::vector<raw::Trigger>` (`PrimitiveTags`): collection of trigger
 *   primitives. Each object must encode the location of the primitive in its
 *   trigger bits (see `CryostatFirstBit` configuration parameter).
 * * `std::vector<raw::Trigger>` (optional, `TriggerTags`): collections of
 *   global triggers; each global trigger will count as a trigger primitive
 *   for each cryostat. All triggers in all data products are added.
 * * `std::vector<raw::OpDetWaveform>` (`WaveformTags`): PMT waveforms already
 *   digitized, source of the data for the new ones.
 * 
 * 
 * Output
 * =======
 * 
 * * `std::vector<raw::OpDetWaveform>` with the rearranged readout samples.
 *   They are sorted by channel, then by time, like in the output of
 *   `icarus::opdet::SimPMTIcarus` module.
 * * `std::vector<sbn::OpDetWaveformMeta>` and associations with the optical
 *   waveform data product, with "summary" informations about each waveform.
 * 
 * 
 * Requirements
 * =============
 * 
 * This module currently requires LArSoft services:
 * * `Geometry` for discovery of optical detector channels and their cryostat;
 * * `DetectorClocksService` for timing conversions and settings;
 * * `NuRandomService` for the noise generation (even if a non-stochastic
 *   generator is chosen).
 * 
 * 
 */
class icarus::opdet::SimPMTreadout: public art::EDProducer {
    
  using electronics_time = detinfo::timescales::electronics_time;
  
  /// Type of trigger bits supported by `raw::Trigger`.
  using TriggerBits_t = decltype(std::declval<raw::Trigger>().TriggerBits());
  
  /// Number of bits supported by `raw::Trigger`.
  static constexpr std::size_t NTriggerBits = sizeof(TriggerBits_t) * 8U;
  
    public:
  
  using microseconds = util::quantities::intervals::microseconds; // alias
  
  /// Configuration parameters for PMT readout.
  struct ReadoutParams_t {
    
    unsigned int bufferSize; ///< Size of the full readout buffer [samples]
    unsigned int preSize; ///< Size of the buffer collected before trigger.
    
    unsigned int preTriggerSize() const { return preSize; }
    unsigned int postTriggerSize() const { return bufferSize - preSize; }
    
  }; // ReadoutParams_t
  
  
  /// Module configuration interface.
  struct Config {
    using Comment = fhicl::Comment;
    using Name = fhicl::Name;
    
    struct ReadoutConfig {
      
      fhicl::Atom<unsigned int> WindowSize{
        Name{ "WindowSize" },
        Comment{ "readout window size in samples" }
        };
      
      fhicl::Atom<float> PreTrigFraction{
        Name{ "PreTrigFraction" },
        Comment{ "fraction of readout coming before the trigger primitive" }
        };
      
    }; // ReadoutConfig
    
    
    fhicl::Sequence<art::InputTag> WaveformTags{
      Name{ "WaveformTags" },
      Comment{ "input PMT waveforms" }
      };
    
    fhicl::Sequence<art::InputTag> PrimitiveTags{
      Name{ "PrimitiveTags" },
      Comment{ "input trigger primitives" }
      };
    
    fhicl::Sequence<art::InputTag> TriggerTags{
      Name{ "TriggerTags" },
      Comment{ "additional global trigger(s) generating primitives" },
      std::vector<art::InputTag>{}
      };
    
    fhicl::Sequence<microseconds> FixedTriggerTimes{
      Name{ "FixedTriggerTimes" },
      Comment{ "additional global trigger(s) as time w.r.t. beam gate" },
      std::vector<microseconds>{}
      };
    
    fhicl::DelegatedParameter Pedestal {
      fhicl::Name{ "Pedestal" },
      fhicl::Comment{
        "parameters for the pedestal and electronics noise generation algorithm"
        " (PMTpedestalGeneratorTool tool)"
        }
      };
    
    fhicl::Atom<unsigned int> CryostatFirstBit{
      Name("CryostatFirstBit"),
      Comment(
        "number of the trigger bit representing the first cryostat (< "
        + std::to_string(NTriggerBits) + ")"),
      NTriggerBits
      };
    
    rndm::SeedAtom NoiseGeneratorSeed {
      Name("NoiseGeneratorSeed"),
      Comment("fix the seed for stochastic electronics noise generation")
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "category name for the message facility stream to be used" },
      "SimPMTreadout"
      };
    
    fhicl::TableAs<ReadoutParams_t, ReadoutConfig> Readout {
      Name{ "Readout" },
      Comment{ "Readout parameters" }
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;

  
  explicit SimPMTreadout(Parameters const& config);
  
  void produce(art::Event& event) override;
  
  
    private:
  
  using nanoseconds = util::quantities::intervals::nanoseconds; // alias
    
  /// Type of the pedestal and electronics noise generator algorithm.
  using PedestalGenerator_t
    = icarus::opdet::PMTpedestalGeneratorTool::Generator_t;
  
  /// Type of data product and its association(s).
  template <typename Data, typename... AssnsTo>
  using DataAndAssns
    = std::tuple<std::vector<Data>, art::Assns<Data, AssnsTo>...>;
  
  /// Class to manage the input waveforms.
  class WaveformManager {
    
    /// Pointers to the managed waveforms, sorted by channel.
    std::vector<std::vector<raw::OpDetWaveform const*>> fByChannel;
    
      public:
    
    WaveformManager(unsigned int nChannels)
      : fByChannel{ static_cast<std::size_t>(nChannels) }
      {}
    
    /// Returns (pointers to) the waveforms in the specified `channel`.
    std::vector<raw::OpDetWaveform const*> const& operator[]
      (raw::Channel_t channel) const
      { return fByChannel[static_cast<std::size_t>(channel)]; }
    
    /// Returns (pointers to) the waveforms in the specified `channel`.
    std::vector<raw::OpDetWaveform const*> const& at
      (raw::Channel_t channel) const
      { return fByChannel.at(static_cast<std::size_t>(channel)); }
    
    /// Adds to the manager references to all `waveforms` (not sorted).
    WaveformManager& add(std::vector<raw::OpDetWaveform> const& waveforms);
    
    /// Enforces order in increasing timestamp of all waveforms in each channel.
    void sortInTime();
    
  }; // class WaveformManager

  /// Tracks the required primitives per optical detector channel.
  class PrimitiveManager {
    
    /// Cryostat ID of each PMT channel.
    std::vector<geo::CryostatID> const fCryostatOf;
    
    /// Cached primitives per cryostat.
    mutable std::vector<std::vector<electronics_time>> fPrimitives;
    
    /// Changes happened since last cache refresh.
    mutable std::atomic_bool fModified{ false };
    
    mutable std::mutex fRefreshMutex; ///< Lock used while refreshing caches.
    
    void modify() { fModified = true; }
    
    /// Refreshes all caches. Note that this is a "const" method.
    void refreshCaches(bool force = false) const;
    
    /// Creates a map of all PMT channels into cryostat ID.
    static std::vector<geo::CryostatID> makeChannelMap
      (geo::GeometryCore const& geom);
    
      public:
    
    PrimitiveManager(geo::GeometryCore const& geom);
    
    /// Returns all the readout windows for the specified channel.
    std::vector<electronics_time> const& forChannel
      (raw::Channel_t channel) const;
    
    /// Adds a primitive at the specified `time` for all channels in `cryo`.
    void addPrimitiveForCryostat(electronics_time time, geo::CryostatID cryo);
    
    /// Adds a primitive at the specified `time` for all channels.
    void addPrimitiveToAllChannels(electronics_time time);
    
    /// Prepares for queries; not essential, but add deterministic flow.
    void prepare() { refreshCaches(); }
    
  }; // PrimitiveManager
  
  using Window_t = PMTReadoutWindowMaker<electronics_time>::Window_t;
  

  // --- BEGIN -- Configuration ------------------------------------------------
  /// Input tags for already digitized waveforms.
  std::vector<art::InputTag> const fWaveformTags;
  
  /// Input tags for trigger primitives.
  std::vector<art::InputTag> const fPrimitiveTags;
  
  /// Input tags for additional global triggers.
  std::vector<art::InputTag> const fTriggerTags;
  
  /// Time of the added, artificial global triggers.
  std::vector<microseconds> const fFixedTriggerTimes;
  
  ReadoutParams_t const fReadoutParams; ///< PMT readout parameters.
  
  TriggerBits_t const fCryostatZeroMask; ///< Bit mask for cryostat `0`.
  
  std::string const fLogCategory; ///< Message facility category name.
  
  // --- END ---- Configuration ------------------------------------------------
  
  
  // --- BEGIN -- Cached values ------------------------------------------------
  
  geo::GeometryCore const& fGeom; ///< Geometry service provider.
  detinfo::DetectorTimings const fDetTimings; ///< Timing conversion utility.
  
  unsigned int const fNOpChannels; ///< Number of optical detector channels.
  unsigned int const fNCryostats; ///< Number of cryostats in the detector.
  
  nanoseconds const fOpticalTick; ///< PMT digitization sampling time.
  
  /// Random engine for the electronics noise generation algorithm.
  CLHEP::HepRandomEngine& fNoiseGeneratorEngine;
  
  // --- END ---- Cached values ------------------------------------------------
  
  /// Pedestal and electronics noise generation algorithm.
  std::unique_ptr<PedestalGenerator_t> const fPedestalGenerator;
  
  
  /// Registers local primitives (affecting only part of the detector).
  void addLocalPrimitives(
    PrimitiveManager& manager,
    std::vector<raw::Trigger> const& primitives
    ) const;
  
  /// Registers global primitives (affecting the whole detector).
  void addGlobalPrimitives(
    PrimitiveManager& manager,
    std::vector<raw::Trigger> const& primitives
    ) const;
  
  /// Creates waveforms on a single channel, filling the specified readout
  /// windows, from the samples in `source` waveforms.
  std::vector<raw::OpDetWaveform> makeWaveforms(
    raw::Channel_t channel, std::uint64_t beamGateTimestamp,
    std::vector<raw::OpDetWaveform const*> const& source,
    std::vector<Window_t> const& readoutWindows
    ) const;
  
  /// Creates waveform metadata and associations for all `waveforms`.
  DataAndAssns<sbn::OpDetWaveformMeta, raw::OpDetWaveform> makeWaveformMetadata(
    std::vector<raw::OpDetWaveform> const& waveforms,
    detinfo::DetectorTimings const& detTimings,
    art::Event const& event
    ) const;
  
  /// Returns which cryostats the trigger belongs to (invalid if none/unknown).
  std::vector<geo::CryostatID> cryostatsOf(raw::Trigger const& trigger) const;

  
  
}; // class icarus::opdet::SimPMTreadout


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  template <typename T, typename=std::enable_if_t<std::is_rvalue_reference_v<T&&>>>
  std::unique_ptr<T> moveToUniquePtr(T&& data)
    { return std::make_unique<T>(std::move(data)); }
  
  
  // ---------------------------------------------------------------------------
  /// Moves the content of nested collection `src` flattening the outmost tier.
  template <typename Coll, typename=std::enable_if_t<std::is_rvalue_reference_v<Coll&&>>>
  std::vector<typename std::decay_t<Coll>::value_type::value_type> flatten
    (Coll&& src)
  {
    
    using SrcColl_t = typename std::decay_t<Coll>::value_type;
    using value_type = typename SrcColl_t::value_type;
    using DestColl_t = std::vector<value_type>;
    
    DestColl_t dest;
    for (SrcColl_t& coll: src) {
      for (value_type& value: coll) dest.push_back(std::move(value));
    } // for src
    
    return dest;
  } // flatten()
  
  
  // ---------------------------------------------------------------------------
  template <typename T>
  constexpr T bitMask(std::size_t i) {
    constexpr std::size_t nBits = sizeof(T) * 8;
    return (i >= nBits)? T{ 0 }: T{ 1 } << i;
  } // bitMask()
  
  
  // ---------------------------------------------------------------------------
  
} // local namespace


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  
  SimPMTreadout::ReadoutParams_t convert
    (icarus::opdet::SimPMTreadout::Config::ReadoutConfig const& config)
  {
    
    return {
        config.WindowSize()                                      // bufferSize
      , static_cast<unsigned int>                                // preSize
          (std::round(config.WindowSize() * config.PreTrigFraction()))
      };
    
  } // convert(SimPMTreadout::Config::ReadoutConfig)
  
} // namespace icarus::opdet


// -----------------------------------------------------------------------------
// --- icarus::opdet::SimPMTreadout implementation
// -----------------------------------------------------------------------------
// --- icarus::opdet::SimPMTreadout::WaveformManager
// -----------------------------------------------------------------------------
auto icarus::opdet::SimPMTreadout::WaveformManager::add
  (std::vector<raw::OpDetWaveform> const& waveforms) -> WaveformManager&
{
  
  for (raw::OpDetWaveform const& waveform: waveforms)
    fByChannel.at(waveform.ChannelNumber()).push_back(&waveform);
  
  sortInTime();
  
  return *this;
} // icarus::opdet::SimPMTreadout::WaveformManager::add()


// -----------------------------------------------------------------------------
void icarus::opdet::SimPMTreadout::WaveformManager::sortInTime() {
  
  using std::begin, std::end;
  
  // standard ordering for waveforms goes channel first, then timestamp
  for (std::vector<raw::OpDetWaveform const*>& waveforms: fByChannel) {
    std::sort
      (begin(waveforms), end(waveforms), [](auto a, auto b){ return *a < *b; });
  }
  
} // icarus::opdet::SimPMTreadout::WaveformManager::sortInTime()


// -----------------------------------------------------------------------------
// --- icarus::opdet::SimPMTreadout::PrimitiveManager
// -----------------------------------------------------------------------------
void icarus::opdet::SimPMTreadout::PrimitiveManager::refreshCaches
  (bool force /* = false */) const
{
  /*
   * All in all, this class is still not thread-safe, because it would need
   * explicit protection in the non-const methods.
   * I am cheap and am only guaranteeing thread-safety for const ones.
   */
  if (!force && !fModified) return;
  
  std::lock_guard guard{ fRefreshMutex };
  
  if (!force && !fModified) return; // again. just in case something happened
  
  // for correct multithreading, this would need a lock
  using std::begin, std::end;
  for (std::vector<electronics_time>& primitives: fPrimitives)
    std::sort(begin(primitives), end(primitives));
  
  fModified = false;
  
} // icarus::opdet::SimPMTreadout::PrimitiveManager::refreshCaches()


// -----------------------------------------------------------------------------
std::vector<geo::CryostatID>
icarus::opdet::SimPMTreadout::PrimitiveManager::makeChannelMap
  (geo::GeometryCore const& geom)
{
  std::vector<geo::CryostatID> cryos;
  for (auto channel: util::counter(geom.NOpChannels()))
    cryos.push_back(geom.OpDetGeoFromOpChannel(channel).ID());
  return cryos;
} // icarus::opdet::SimPMTreadout::PrimitiveManager::makeChannelMap()


// -----------------------------------------------------------------------------
icarus::opdet::SimPMTreadout::PrimitiveManager::PrimitiveManager
  (geo::GeometryCore const& geom)
  : fCryostatOf{ makeChannelMap(geom) }
  , fPrimitives{ geom.Ncryostats() }
  {}


// -----------------------------------------------------------------------------
auto icarus::opdet::SimPMTreadout::PrimitiveManager::forChannel
  (raw::Channel_t channel) const -> std::vector<electronics_time> const&
{
  refreshCaches();
  return fPrimitives[fCryostatOf.at(channel).Cryostat];
} // icarus::opdet::SimPMTreadout::PrimitiveManager::forChannel()


// -----------------------------------------------------------------------------
void icarus::opdet::SimPMTreadout::PrimitiveManager::addPrimitiveForCryostat
  (electronics_time time, geo::CryostatID cryo)
{
  if (!cryo) {
    // this is probably coming from some misunderstanding with the trigger bits
    // which made the cryostat extraction miserably fail
    throw art::Exception{ art::errors::LogicError }
      << "Added a trigger primitive to an invalid cryostat: " << cryo;
  }
  fPrimitives.at(cryo.Cryostat).push_back(time);
  modify();
} // icarus::opdet::SimPMTreadout::PrimitiveManager::addPrimitiveForCryostat()


// -----------------------------------------------------------------------------
void icarus::opdet::SimPMTreadout::PrimitiveManager::addPrimitiveToAllChannels
  (electronics_time time)
{
  for (std::vector<electronics_time>& primitives: fPrimitives)
    primitives.push_back(time);
  modify();
} // icarus::opdet::SimPMTreadout::PrimitiveManager::addPrimitiveToAllChannels()


// -----------------------------------------------------------------------------
// --- icarus::opdet::SimPMTreadout implementation
// -----------------------------------------------------------------------------
icarus::opdet::SimPMTreadout::SimPMTreadout(Parameters const& config)
  : art::EDProducer{ config }
  // configuration
  , fWaveformTags     { config().WaveformTags() }
  , fPrimitiveTags    { config().PrimitiveTags() }
  , fTriggerTags      { config().TriggerTags() }
  , fFixedTriggerTimes{ config().FixedTriggerTimes() }
  , fReadoutParams    { config().Readout() }
  , fCryostatZeroMask { bitMask<TriggerBits_t>(config().CryostatFirstBit()) }
  , fLogCategory      { config().LogCategory() }
  // caches
  , fGeom       { *(lar::providerFrom<geo::Geometry>()) }
  , fDetTimings { icarus::ns::util::makeDetTimings() }
  , fNOpChannels{ fGeom.NOpChannels() }
  , fNCryostats { fGeom.Ncryostats() }
  , fOpticalTick{ fDetTimings.OpticalClockPeriod() }
  , fNoiseGeneratorEngine{
      art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
        createEngine(0, "HepJamesRandom", "ElectronicsNoise"),
        "HepJamesRandom", "ElectronicsNoise", config().NoiseGeneratorSeed
        ).get()
      }
  , fPedestalGenerator{
      art::make_tool<icarus::opdet::PMTpedestalGeneratorTool>
        (config().Pedestal.get<fhicl::ParameterSet>())
        ->makeGenerator(fNoiseGeneratorEngine)
    }
{
  
  //
  // configuration checks
  //
  if (config().CryostatFirstBit() + fNCryostats > NTriggerBits) {
    throw art::Exception{ art::errors::Configuration }
      << config().CryostatFirstBit.name()
      << " must be small enough to accommodate all " << fNCryostats
      << " cryostats in the " << fGeom.DetectorName()
      << " detector, i.e. it must be no larger than "
      << (NTriggerBits - fNCryostats) << "\n";
  }
  
  //
  // input data product declaration
  //
  for (art::InputTag const& tag: fWaveformTags)
    consumes<std::vector<raw::OpDetWaveform>>(tag);
  for (art::InputTag const& tag: fPrimitiveTags)
    consumes<std::vector<raw::Trigger>>(tag);
  for (art::InputTag const& tag: fTriggerTags)
    consumes<std::vector<raw::Trigger>>(tag);
  
  //
  // output data product declaration
  //
  produces<std::vector<raw::OpDetWaveform>>();
  produces<std::vector<sbn::OpDetWaveformMeta>>();
  produces<art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>>();
  
  
  { // --- BEGIN -- configuration dump -----------------------------------------
    
    mf::LogInfo log{ fLogCategory };
    log << "Configuration parameters:"
      << "\n - waveforms:";
    if (fWaveformTags.empty()) log << " none?!";
    else {
      auto iTag = fWaveformTags.cbegin(), tend = fWaveformTags.cend();
      log << " " << iTag->encode();
      while (++iTag != tend) log << ", " << iTag->encode();
    }
    log << "\n - trigger primitives:";
    if (fPrimitiveTags.empty()) log << " none?!";
    else {
      auto iTag = fPrimitiveTags.cbegin(), tend = fPrimitiveTags.cend();
      log << " " << iTag->encode();
      while (++iTag != tend) log << ", " << iTag->encode();
    }
    if (!fTriggerTags.empty()) {
      log << "\n - additional global triggers:";
      auto iTag = fTriggerTags.cbegin(), tend = fTriggerTags.cend();
      log << " " << iTag->encode();
      while (++iTag != tend) log << ", " << iTag->encode();
    }
    
    unsigned int const firstBit = config().CryostatFirstBit();
    log
      << "\n - readout buffer: " << fReadoutParams.bufferSize << " samples ("
        << fReadoutParams.preTriggerSize() << " + "
        << fReadoutParams.postTriggerSize() << ")"
      << "\n - trigger bit for the first cryostat: #" << firstBit
      ;
    if (fNCryostats > 1) log << "-" << (firstBit + fNCryostats - 1);
    log << " (0x" << std::hex << fCryostatZeroMask;
    if (fNCryostats > 1)
      log << "-0x" << (fCryostatZeroMask << (fNCryostats-1));
    log << std::dec << ")";
    
    if (!fFixedTriggerTimes.empty()) {
      auto itTime = fFixedTriggerTimes.begin();
      auto const tend = fFixedTriggerTimes.end();
      log << "\n - additional global triggers at time (w.r.t. beam gate): "
        << *itTime;
      while (++itTime != tend) log << ", " << *itTime;
    }
    
    log << "\n - pedestal generation: " << *fPedestalGenerator;
    
  } // --- END ---- configuration dump -----------------------------------------
  
} // icarus::opdet::SimPMTreadout::SimPMTreadout()


// -----------------------------------------------------------------------------
void icarus::opdet::SimPMTreadout::produce(art::Event& event) {
  
  using detinfo::timescales::electronics_time;
  
  detinfo::DetectorTimings const detTimings
    = icarus::ns::util::makeDetTimings(event);
  
  //
  // determine the primitives and the channels they affect
  //
  PrimitiveManager primitives{ fGeom };
  
  // add the trigger primitives
  for (art::InputTag const& tag: fPrimitiveTags) {
    addLocalPrimitives
      (primitives, event.getProduct<std::vector<raw::Trigger>>(tag));
  }
  
  // add global triggers
  for (art::InputTag const& tag: fTriggerTags) {
    addGlobalPrimitives
      (primitives, event.getProduct<std::vector<raw::Trigger>>(tag));
  }
  
  // add artificial primitives
  electronics_time const beamGateTime = detTimings.BeamGateTime();
  for (microseconds const time: fFixedTriggerTimes)
    primitives.addPrimitiveToAllChannels(beamGateTime + time);
  
  //
  // read the input waveforms
  //
  WaveformManager simWaveforms{ fNOpChannels };
  for (art::InputTag const& tag: fWaveformTags)
    simWaveforms.add(event.getProduct<std::vector<raw::OpDetWaveform>>(tag));
  
  
  //
  // create the new waveforms
  //
  
  // transforms trigger primitives into readout windows
  // (we might do that once per cryostat rather than per channel...):
  PMTReadoutWindowMaker<electronics_time> const makeReadoutWindows{
    fReadoutParams.bufferSize * fOpticalTick,
    fReadoutParams.preTriggerSize() * fOpticalTick,
    fLogCategory + "_PMTReadoutWindowMaker"
    };
  
  std::uint64_t const beamGateTimestamp = event.time().value();
  std::vector<std::vector<raw::OpDetWaveform>> readoutWaveforms{ fNOpChannels };
  for (auto const channel: util::counter<raw::Channel_t>(fNOpChannels)) {
    
    std::vector<Window_t> const readoutWindows
      = makeReadoutWindows(primitives.forChannel(channel));
    std::vector<raw::OpDetWaveform const*> const waveforms
      = simWaveforms[channel];
    
    readoutWaveforms[channel]
      = makeWaveforms(channel, beamGateTimestamp, waveforms, readoutWindows);
    
  } // for channels
  
  
  //
  // sort the output and put it into the event
  //
  std::vector<raw::OpDetWaveform> allReadoutWaveforms
    = flatten(std::move(readoutWaveforms));
  
  auto [ allReadoutWaveformsMeta, waveToMeta ]
    = makeWaveformMetadata(allReadoutWaveforms, detTimings, event);
  
  event.put(moveToUniquePtr(std::move(allReadoutWaveforms)));
  event.put(moveToUniquePtr(std::move(allReadoutWaveformsMeta)));
  event.put(moveToUniquePtr(std::move(waveToMeta)));
  
} // icarus::opdet::SimPMTreadout::produce()


// ---------------------------------------------------------------------------
std::vector<geo::CryostatID> icarus::opdet::SimPMTreadout::cryostatsOf
  (raw::Trigger const& trigger) const
{
  TriggerBits_t const bits = trigger.TriggerBits();
  
  std::vector<geo::CryostatID> cryostats;
  TriggerBits_t cryoMask = fCryostatZeroMask;
  for (auto const c: util::counter<unsigned int>(fNCryostats)) {
    
    if (bits & cryoMask) cryostats.emplace_back(c);
    cryoMask <<= 1;
    
  } // for
  
  return cryostats;
} // geo::CryostatID icarus::opdet::SimPMTreadout::cryostatOf()


// ---------------------------------------------------------------------------
void icarus::opdet::SimPMTreadout::addLocalPrimitives
  (PrimitiveManager& manager, std::vector<raw::Trigger> const& primitives) const
{
  
  using microsecond = util::quantities::points::microsecond;
  
  for (raw::Trigger const& trigger: primitives) {
    
    electronics_time const time = microsecond{ trigger.TriggerTime() };
    
    mf::LogTrace log{ fLogCategory };
    log << "Adding trigger primitive at time " << time << " (bits: 0x"
      << std::hex << trigger.TriggerBits() << std::dec << ")";
    
    for (geo::CryostatID const& cid: cryostatsOf(trigger)) {
      log << " to " << cid;
      manager.addPrimitiveForCryostat(time, cid);
    }
    
  } // for
  
} // icarus::opdet::SimPMTreadout::addLocalPrimitives()


// ---------------------------------------------------------------------------
void icarus::opdet::SimPMTreadout::addGlobalPrimitives
  (PrimitiveManager& manager, std::vector<raw::Trigger> const& primitives) const
{
  
  using microsecond = util::quantities::points::microsecond;
  
  for (raw::Trigger const& trigger: primitives) {
    
    electronics_time const time = microsecond{ trigger.TriggerTime() };
    
    mf::LogTrace{ fLogCategory }
      << "Adding global trigger primitive at time " << time << " (bits: 0x"
      << std::hex << trigger.TriggerBits() << std::dec << ")";
    
    manager.addPrimitiveToAllChannels(time);
    
  } // for
  
} // icarus::opdet::SimPMTreadout::addGlobalPrimitives()


// ---------------------------------------------------------------------------
std::vector<raw::OpDetWaveform> icarus::opdet::SimPMTreadout::makeWaveforms(
  raw::Channel_t channel, std::uint64_t beamGateTimestamp,
  std::vector<raw::OpDetWaveform const*> const& source,
  std::vector<Window_t> const& readoutWindows
) const {
  
  // assumption: source waveforms are sorted by timestamp
  // assumption: readoutWindows sorted by timestamp and non-contiguous
  
  using microsecond = util::quantities::points::microsecond;
  
  auto const toTicks = [this](nanoseconds interval)
    { return static_cast<int>(std::round(interval / fOpticalTick)); };
  auto const toTimestamp = [this,beamGateTimestamp](electronics_time t)
    {
      return beamGateTimestamp + static_cast<std::int64_t>
        (std::round(nanoseconds(t - fDetTimings.BeamGateTime()).value()));
    };
  
  std::vector<raw::OpDetWaveform> waveforms;
  waveforms.reserve(readoutWindows.size());
  
  // --- BEGIN DEBUG -----------------------------------------------------------
  {
    mf::LogTrace log{ "SimPMTreadout" };
    log << "makeWaveforms(channel=" << channel << ", { " << source.size()
      << " waveforms }, { " << readoutWindows.size() << " windows })";
    for (auto const& [ i, w ]: util::enumerate(source)) {
      double const length = w->size() * fOpticalTick.value() / 1000.;
      log << "\n  [S" << i << "](" << ((void*) w)
        << ") ch=" << w->ChannelNumber() << " t=" << w->TimeStamp()
        << " -- " << (w->TimeStamp() + length) << " l=" << length
        << " s=" << w->size();
    } // for
    for (auto const& [ i, r ]: util::enumerate(readoutWindows)) {
      log << "\n  [R" << i << "]"
        << " " << r.start() << " -- " << r.stop() << " (" << r.width() << ")";
    } // for
  }
  // --- END   DEBUG -----------------------------------------------------------
  
  auto nextSource = source.begin();
  auto const sourceEnd = source.end();
  
  for (Window_t const& window: readoutWindows) {
    // --- BEGIN DEBUG ---------------------------------------------------------
    MF_LOG_TRACE("SimPMTreadout") << "Window: " << window.start() << " -- "
      << window.stop() << " (" << window.width() << ")";
    // --- END   DEBUG ---------------------------------------------------------
    
    if (window.empty()) continue;
    
    std::size_t const nSamples = toTicks(window.width());
    std::vector<raw::ADC_Count_t> samples(nSamples); // allocated, uninitialized
    auto itSample = samples.begin();
    auto const send [[maybe_unused]] = samples.end();
    
    electronics_time neededTime = window.start();
    
    // --- BEGIN DEBUG ---------------------------------------------------------
    auto sampleIndex [[maybe_unused]]
      = [&itSample,b=samples.begin()](){ return itSample - b; };
    MF_LOG_TRACE("SimPMTreadout")
      << "Expected waveform with " << samples.size() << " samples;"
      << " neededTime=" << neededTime << "  itSample=#" << sampleIndex();
    // --- END   DEBUG ---------------------------------------------------------
    while (neededTime < window.stop()) {
      
      // find the first waveform that might include the time we need
      while (nextSource != sourceEnd) {
        assert(*nextSource); // no null pointers please
        electronics_time const srcStart
          { microsecond{ (*nextSource)->TimeStamp() } };
        electronics_time const srcStop
          { srcStart + (*nextSource)->size() * fOpticalTick };
        
        // --- BEGIN DEBUG -----------------------------------------------------
        MF_LOG_TRACE("SimPMTreadout")
          << "Next source: " << ((void*) *nextSource) << " t=" << srcStart
          << " -- " << srcStop << " s=" << (*nextSource)->size();
        // --- END   DEBUG -----------------------------------------------------
        
        if (srcStop > neededTime) break;
        ++nextSource;
      } // while source available
      
      // start of the data of the next useful source waveform (may be never)
      electronics_time const srcStart{
        (nextSource == sourceEnd)
        ? std::numeric_limits<microsecond>::max()
        : microsecond{ (*nextSource)->TimeStamp() } 
        };
      
      // --- BEGIN DEBUG -------------------------------------------------------
      if (nextSource == sourceEnd) {
        MF_LOG_TRACE("SimPMTreadout")
          << "No source can help with neededTime=" << neededTime;
      }
      else {
        MF_LOG_TRACE("SimPMTreadout")
          << "Source might help with neededTime=" << neededTime
          << " (srcStart=" << srcStart << ")";
      }
      // --- END   DEBUG -------------------------------------------------------
      
      if (neededTime < srcStart) {
        // the data at needed time is not available: fill with noise
        electronics_time const noiseEnd = std::min(srcStart, window.stop());
        std::size_t const nNoiseSamples = toTicks(noiseEnd - neededTime);
        // --- BEGIN DEBUG -----------------------------------------------------
        MF_LOG_TRACE("SimPMTreadout") << "neededTime=" << neededTime
          << " < srcStart=" << srcStart << ": noise time";
        // --- END   DEBUG -----------------------------------------------------
        
        // while the pedestal generator is working in floating point,
        // our source and destination is in rounded ADC counts;
        // we need to convert (emitting directly in ADC counts is only good if
        // the pedestal is already integer, otherwise we stack yet another
        // rounding)
        std::vector<util::quantities::counts_f> extendedNoise(nNoiseSamples);
        fPedestalGenerator->fill
          (channel, toTimestamp(neededTime), extendedNoise);
        for (util::quantities::counts_f const noise: extendedNoise) {
          *(itSample++)
            = static_cast<raw::ADC_Count_t>(std::round(noise.value()));
        }
        assert(itSample <= send);
        
        neededTime = noiseEnd;
        // --- BEGIN DEBUG -----------------------------------------------------
        MF_LOG_TRACE("SimPMTreadout") << "Added " << nNoiseSamples
          << " samples of noise, now neededTime=" << neededTime
          << " itSample=#" << sampleIndex()
          << "  window: " << window.start() << " -- " << window.stop();
        // --- END   DEBUG -----------------------------------------------------
      } // if add noise
      
      if (neededTime < window.stop()) {
        // adding noise (if any) was not enough: there must be some data to copy
        // --- BEGIN DEBUG -----------------------------------------------------
        MF_LOG_TRACE("SimPMTreadout")
          << "neededTime=" << neededTime << "  window: "
          << window.start() << "--" << window.stop() << ": copy!";
        // --- END   DEBUG -----------------------------------------------------
        assert(nextSource != sourceEnd);
        electronics_time const srcStart
          { microsecond{ (*nextSource)->TimeStamp() } };
        electronics_time const srcStop
          { srcStart + (*nextSource)->size() * fOpticalTick };
        assert(neededTime < srcStop);
        assert(neededTime >= srcStart);
        
        electronics_time const dataEndTime = std::min(srcStop, window.stop());
        
        auto const dbegin
          = (*nextSource)->cbegin() + toTicks(neededTime - srcStart);
        auto const dend = dbegin + toTicks(dataEndTime - neededTime);
        
        // --- BEGIN DEBUG -----------------------------------------------------
        MF_LOG_TRACE("SimPMTreadout") << "Copy data from " << srcStart << " (#"
          << (dbegin - (*nextSource)->cbegin()) << ") to " << dataEndTime
          << " (#" << (dend - (*nextSource)->cbegin()) << ")"
          << " into itSample=#" << sampleIndex();
        // --- END   DEBUG -----------------------------------------------------
        
        itSample = std::copy(dbegin, dend, itSample);
        neededTime = dataEndTime;
        
        // if this source is over, point to next one for the next iteration
        if (dend >= (*nextSource)->cend()) ++nextSource;
        
        // --- BEGIN DEBUG -----------------------------------------------------
        MF_LOG_TRACE("SimPMTreadout") << "Copy is over, neededTime="
          << neededTime << " itSample=#" << sampleIndex();
        // --- END   DEBUG -----------------------------------------------------
        
      } // if copy more data
    
      // --- BEGIN DEBUG -------------------------------------------------------
      MF_LOG_TRACE("SimPMTreadout")
        << "Done with this source+window, neededTime=" << neededTime
        << "  window.stop()=" << window.stop();
      // --- END   DEBUG -------------------------------------------------------
    
    } // while needed time
    
    // --- BEGIN DEBUG ---------------------------------------------------------
    MF_LOG_TRACE("SimPMTreadout") << "Now neededTime=" << neededTime
      << " window.stop()=" << window.stop() << " => done";
    // --- END   DEBUG ---------------------------------------------------------
    
    assert(itSample == send);
    
    // add the new waveform (with a weird way dictated by a weird interface)
    waveforms.emplace_back(window.start().value(), channel);
    waveforms.back().Waveform() = std::move(samples);
    
    // --- BEGIN DEBUG ---------------------------------------------------------
    {
      raw::OpDetWaveform const& ow [[maybe_unused]] = waveforms.back();
      double const length [[maybe_unused]]
        = ow.size() * fOpticalTick.value() / 1000.;
      MF_LOG_TRACE("SimPMTreadout") 
        << "For window " << window.start() << " -- " << window.stop()
        << " added ch=" << ow.ChannelNumber() << " t=" << ow.TimeStamp()
        << " -- " << (ow.TimeStamp() + length) << " l=" << length
        << " s=" << ow.size() << std::endl;
    }
    // --- END   DEBUG ---------------------------------------------------------
    
  } // for
  
  return waveforms;
  
} // icarus::opdet::SimPMTreadout::makeWaveforms()


// ---------------------------------------------------------------------------
auto icarus::opdet::SimPMTreadout::makeWaveformMetadata(
  std::vector<raw::OpDetWaveform> const& waveforms,
  detinfo::DetectorTimings const& detTimings,
  art::Event const& event
) const -> DataAndAssns<sbn::OpDetWaveformMeta, raw::OpDetWaveform> {
  
  sbn::OpDetWaveformMetaMaker const metaMaker{ detTimings };
  
  art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr{ event };
  art::PtrMaker<sbn::OpDetWaveformMeta> const makeWaveformMetaPtr{ event };
  
  art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> waveToMeta;
  std::vector<sbn::OpDetWaveformMeta> meta;
  
  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
    meta.push_back(metaMaker(waveform));
    waveToMeta.addSingle
      (makeWaveformPtr(iWaveform), makeWaveformMetaPtr(iWaveform));
  }
  
  return { std::move(meta), std::move(waveToMeta) };
} // icarus::opdet::SimPMTreadout::makeWaveformMetadata()


// ---------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::opdet::SimPMTreadout)


// ---------------------------------------------------------------------------
