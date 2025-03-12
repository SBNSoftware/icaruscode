/**
 * @file   DiscriminatedAdderSignal_module.cc
 * @brief  Module to put together and discriminate adder signals.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 6, 2023
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h" // GateOps
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icarusalg/Utilities/TimeIntervalConfig.h" // icarus::...::TimeInterval
#include "icarusalg/Utilities/GroupByIndex.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/Utilities/intervals_fhicl.h" // for microseconds parameters
#include "lardataalg/Utilities/quantities/electromagnetism.h" // millivolt
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib/maybe_ref.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <iomanip>
#include <algorithm> // std::transform()
#include <numeric> // std::accumulate()
#include <iterator> // std::back_inserter(), std::iterator_traits
#include <utility> // std::move(), std::pair
#include <functional> // std::mem_fn()
#include <memory> // std::unique_ptr
#include <map>
#include <valarray>
#include <vector>
#include <string>
#include <optional>
#include <cassert>


//------------------------------------------------------------------------------
using namespace util::quantities::time_literals;

//------------------------------------------------------------------------------
namespace icarus::trigger { class DiscriminatedAdderSignal; }

/**
 * @brief Produces discriminated optical waveforms for PMT adders.
 * 
 * This module produces a "discriminated waveform" for each adder module and
 * for each configured threshold, via an algorithm of the class of
 * `icarus::trigger::TriggerGateBuilder` (`TriggerGateBuilder` configuration
 * parameter).
 * 
 * The input is PMT optical waveforms: this module composes them into added
 * analogue signals, then discriminates them.
 * No distortion is applied to the PMT waveforms before or after the addition.
 * 
 * Thresholds are ultimately chosen by the tool in charge of actually running
 * the discrimination algorithm. Out of the thresholds that this algorithm
 * produces, it is possible to choose only a subset of them
 * (`SelectThresholds`).
 * 
 * 
 * Output data products
 * =====================
 * 
 * * `std::vector<raw::OpDetWaveform>`, `std::vector<icarus::WaveformBaseline>`
 *   and their one-to-one associations (only if `SaveWaveforms` is set):
 *   waveforms from the adders; in this emulation, the waveforms are
 *   "reproduced" from the already digitized PMT waveforms, in contrast with
 *   the hardware which performs an analogue sum.
 *   The channel number is unique on each of them, but it has no set meaning
 *   beyond that.
 *   In addition, a baseline object is saved for each waveform (their baseline
 *   is `0` by construction since they originate from baseline-subtracted PMT
 *   waveforms), and an association
 *   `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>`between them
 *   (which is trivial since the two data products are
 *   @ref LArSoftProxyDefinitionParallelData "parallel").
 * * `std::vector<sbn::OpDetWaveformMeta>` (only if `SavePMTcoverage` is set):
 *   parallel to the `std::vector<raw::OpDetWaveform>` data product, each entry
 *   in the collection is matched with the corresponding one in the adder
 *   waveform collection; if `SaveWaveforms` is set, an association
 *   `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>` is also produced
 *   (which is also trivial since the two data products are
 *   @ref LArSoftProxyDefinitionParallelData "parallel").
 * * for each threshold _thr_ (in ADC counts), data products are created with
 *   an instance name equivalent to the threshold (e.g. for a threshold of 700
 *   ADC counts, the data product instance name will be `"700"`):
 *     * a data product of type
 *       `std::vector<icarus::trigger::OpticalTriggerGate::GateData_t>`,
 *       with one entry per adder module; each gate object has as index in the
 *       vector the number of the optical channel;
 *     * an _art_ association, one trigger gate data to many optical waveforms,
 *       of each trigger gate data (as above) with all the optical waveforms
 *       which contributed to it; the order of the association pairs is the same
 *       as the one of the trigger gate data, i.e. by channel, and some trigger
 *       gate data entries may be not associated to any waveform and therefore
 *       they can not be listed in the association; association pairs within
 *       the same optical channel are sorted by optical waveform timestamp;
 *       the type of the association is
 *       `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, raw::OpDetWaveform>`.
 *       If `SaveMetaAssns` is set, also an association
 *       art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, sbn::OpDetWaveformMeta>`
 *       (just like the previous one, but associating to the waveform metadata
 *       instead of the waveform itself).
 *     * if `SaveWaveforms` is set, an association
 *       `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, raw::OpDetWaveform>`
 *       to the adder waveforms;
 *       this will have the string `adder` appended to the instance name, to
 *       distinguish it from the associations to the PMT waveforms.
 *     * if `SavePMTcoverage` is set, an association
 *       `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, sbn::OpDetWaveformMeta>`
 *       to the adder waveform metadata; this will have the string `adders`
 *       appended to the instance name, for consistency with the associations
 *       of the same trigger gate objects to the PMT waveforms.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::OpDetWaveform>` (`WaveformTag`): a single waveform for
 *   each recorded optical detector activity; the activity belongs to a single
 *   channel, but there may be multiple waveforms on the same channel. The time
 *   stamp is expected to be from the
 *   @ref DetectorClocksElectronicsTime "electronics time scale" and therefore
 *   expressed in microseconds.
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>` (`BaselineTag`):
 *   the baseline associated to each of the input waveforms.
 * * `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>` (`WaveformTag`;
 *   only if `SaveMetaAssns` is set): associations between the input waveforms
 *   and their metadata (will be used to create associations to the
 *   discriminated waveforms).
 * 
 * 
 * Service requirements
 * ---------------------
 * 
 * The following services are _required_:
 * 
 * * `Geometry` (for window definition)
 * * an implementation of `detinfo::DetectorClocksService`
 * * _art_ message facility
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description DiscriminatedAdderSignal`.
 * 
 * * `WaveformTag` (input tag): the data product containing all optical
 *   detector waveforms.
 * * `BaselineTag` (input tag): the data product containing the association
 *   between each input waveform (`WaveformTag`) and its estimated baseline.
 * * `WaveformMetaTag` (input tag, optional): the data product containing the
 *   association between each input waveform (`WaveformTag`) and its metadata;
 *   if omitted, the tag from `WaveformTag` is used.
 * * `AmplitudeScale` (real number, default: `1.0`): how much to scale each
 *   input waveform.
 * * `TimeInterval` (table, optional): limits the time range of the adders to
 *   the specified time interval; all time points are _relative to the beam gate
 *   opening time_. No more than two of the following must be specified:
 *     * `Start` (time, default: ages ago): start from this time point.
 *     * `End` (time, default: ages in the future): end at this time point.
 *     * `Duration` (time, default: ages): duration of the time interval.
 * * `MissingChannels` (list of integers): list of the PMT channel numbers
 *   (as reported in the `WaveformTag` input collection) to be excluded.
 *   Note that in the hardware bad channels may be included in the adder input.
 * 
 * * `TriggerGateBuilder` (tool configuration): configuration of the _art_ tool
 *   used to discriminate the optional waveforms; the tool interface is
 *   `icarus::trigger::TriggerGateBuilder`.
 * * `SelectThresholds` (sequence of ADC thresholds): if specified, only the
 *   waveforms discriminated with the thresholds in this list will be saved
 *   by the module; if not specified, all the thresholds produced by the
 *   algorithm will be saved. @note If a requested threshold is not eventually
 *   produced by the algorithm, an exception will be thrown by the module.
 * 
 * * `ChannelOffset` (integer, default: `15000`): this offset is added to the
 *   channel number of the waveforms of the adders.
 * * `ChannelsPerAdder` (positive integer, default: `15`): PMT are added in
 *   groups of this many into each adder channel. The PMT are split in walls
 *   and then grouped longitudinally. This number must be a divisor of the
 *   number of PMT on the wall (e.g. in ICARUS there are 90 PMT, so this value
 *   can only be 1, 2, 3, 5, 6, 9, 10, 15, 18, 30, 45 or 90).
 * * `SaveWaveforms` (flag, default: `false`): if set, the waveforms from the
 *   adders are saved, otherwise only the discriminated waveforms are.
 * * `SavePMTcoverage` (flag, default: `true`): produces a collection of
 *   `sbn::OpDetWaveformMeta` representing each of the output waveforms; trigger
 *   tools can use this information in place for the more space-hungry waveforms
 *   for further processing.
 * * `SaveMetaAssns` (flag, default: `true`): if set, an association between the
 *   metadata of the input waveforms (from `WaveformTag`) and the discriminated
 *   waveforms is also saved.
 * * `LogCategory` (string, default: `"DiscriminatedAdderSignal"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 * 
 */
class icarus::trigger::DiscriminatedAdderSignal: public art::SharedProducer {
  
    public:
  
  using microseconds = util::quantities::intervals::microseconds;
  
  /// The type of data produced for a single channel.
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGate::GateData_t;
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> WaveformTag{
      Name{ "WaveformTag" },
      Comment{ "input tag of digitized optical waveform data product" }
      };
    
    fhicl::Atom<art::InputTag> BaselineTag{
      Name{ "BaselineTag" },
      Comment{ "input tag for associations between waveforms and baselines" }
      };
    
    fhicl::OptionalAtom<art::InputTag> WaveformMetaTag{
      Name{ "WaveformMetaTag" },
      Comment{
        "input tag for associations between waveforms and metadata"
        " [default: as WaveformTag]"
        }
      };
    
    fhicl::Atom<float> AmplitudeScale{
      Name{ "AmplitudeScale" },
      Comment{ "scale factor applied to each sample of all input waveforms" },
      1.0 // default value
      };
    
    fhicl::DelegatedParameter TriggerGateBuilder_{
      Name{ "TriggerGateBuilder" },
      Comment{
        "parameters for generating trigger gates from optical channel output"
        }
      };
    
    icarus::ns::fhicl::TimeIntervalOptionalTable<microseconds> TimeInterval{
      Name{ "TimeInterval" },
      Comment{ "limit the adder time to this interval, relative to beam gate" }
      };
    
    fhicl::OptionalSequence<raw::ADC_Count_t> SelectThresholds{
      Name{ "SelectThresholds" },
      Comment{ "thresholds to save (default: all those produced by algorithm)" }
      };
    
    fhicl::Sequence<raw::Channel_t> MissingChannels{
      Name{ "MissingChannels" },
      Comment("do not include the specified channels [default: all included]"),
      std::vector<raw::Channel_t>{}
      };
    
    fhicl::Atom<bool> SaveWaveforms{
      Name{ "SaveWaveforms" },
      Comment{ "write a raw::OpDetWaveform for each of the adder boards" },
      false
      };
    
    fhicl::Atom<bool> SavePMTcoverage{
      Name{ "SavePMTcoverage" },
      Comment
        { "write a sbn::OpDetWaveformMeta collection for adder waveforms" },
      true
      };
    
    fhicl::Atom<bool> SaveMetaAssns{
      Name{ "SaveMetaAssns" },
      Comment{
        "write associations between discriminated waveforms"
        " and their input metadata"
      },
      true
      };
    
    fhicl::Atom<int> ChannelOffset{
      Name{ "ChannelOffset" },
      Comment{ "offset added to the channel number of the adder waveforms" },
      15000
      };
    
    fhicl::Atom<unsigned int> ChannelsPerAdder{
      Name{ "ChannelsPerAdder" },
      Comment{ "number of PMT channels in each adder" },
      15
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "tag of the module output to console via message facility" },
      "DiscriminatedAdderSignal"
      };
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  explicit DiscriminatedAdderSignal
    (Parameters const& config, art::ProcessingFrame const& frame);
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Creates the data products.
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  using electronics_time = detinfo::timescales::electronics_time;
  using WaveformWithBaseline = icarus::trigger::WaveformWithBaseline;
  
  /// Time interval from an unspecified reference.
  using RelTimeInterval_t = icarus::ns::util::TimeInterval<microseconds>;
  
  /// Time interval in electronics time.
  using TimeInterval_t = icarus::ns::util::TimeInterval<electronics_time>;
  
  /// Map of waveform information records per channel.
  using WaveformsByChannel_t
    = icarus::ns::util::GroupByIndex<WaveformWithBaseline>;
  
  /// Conversion from ADC# to mV for an ADC sampling 2 V with 14 bits.
  static constexpr double mVperADC = 2'000. / (1<<14);
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fWaveformTag; ///< Input optical waveform tag.
  
  art::InputTag const fBaselineTag; ///< Input waveform baseline tag.
  
  art::InputTag const fWaveformMetaTag; ///< Input waveform metadata tag.
  
  float const fAmplitudeScale; ///< Scale factor for input waveform amplitude.
  
  RelTimeInterval_t const fTimeInterval; ///< Range of time to create adders on.
  
  /// Thresholds selected for saving, and their instance name.
  std::map<icarus::trigger::ADCCounts_t, std::string> fSelectedThresholds;
  
  std::vector<raw::Channel_t> const fMissingChannels; ///< Channels to skip.
  
  bool const fSaveWaveforms; ///< Whether to save `raw::OpDetWaveform`.
  
  bool const fSavePMTcoverage; ///< Whether to save `sbn::OpDetWaveformMeta`.
  
  bool const fSaveMetaAssns; ///< Whether to save input waveform associations.
  
  int const fAdderChannelOffset; ///< Channel number offset for adder waveforms.
  
  unsigned int const fChannelsPerAdder; ///< Channels per adder.
  
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  /// Algorithms to simulate trigger gates out of optical channel output.
  std::unique_ptr<icarus::trigger::TriggerGateBuilder> fTriggerGateBuilder;
  
  // --- END Algorithms --------------------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  geo::GeometryCore const& fGeom; ///< Cached geometry service provider.
  
  microseconds const fOpticalTick; ///< Duration of an optical tick.
  
  // --- END Service variables -------------------------------------------------
  
  
  /// Suffix used for the instance name of some data products.
  static std::string const AdderDataProductSuffix;
  
  /// Suffix used for the instance name of some data products.
  static std::string const PMTDataProductSuffix;
  
  
  /**
   * @brief Returns the samples of a waveform added from the specified ones.
   * @param timeInterval the interval to constrain the waveform into
   * @param waveforms the list of waveforms (plus baseline) to include
   * @param channel (default: max possible value) channel number to assign
   * @return the added waveform and the list of the contributing waveforms
   * 
   * The returned waveform starts at either `timeInterval.start` or the
   * earliest of the specified waveform start times. Likewise, the waveform ends
   * at either `timeInterval.end` or the latest of the waveform end times.
   * 
   * The input waveforms are assumed to be grouped by channel: each entry in
   * `waveforms` is assumed to be the list of pointers to information records
   * for all the relevant waveform on the same channel (empty lists are ignored,
   * and there is no way to know which channel they might have belonged to).
   * 
   * This function uses `computeTimeInterval()` to determine the actual time
   * interval that the added waveform needs to cover, and it performs the
   * addition with `addWaveformsInInterval()`. See the documentation of those
   * methods for additional details.
   * The returned list includes all and only the waveforms which have at least
   * one sample in the requested time interval.
   */
  std::pair<raw::OpDetWaveform, std::vector<WaveformWithBaseline const*>>
  makeAdderWaveform(
    TimeInterval_t const& timeInterval,
    std::vector<std::vector<WaveformWithBaseline const*>> const& waveforms,
    raw::Channel_t channel = std::numeric_limits<raw::Channel_t>::max()
    ) const;
  
  /**
   * @brief Removes the channels configured as missing from a `channelList`.
   * @param channels the list of channels (see the description)
   * @return `channels` list except for the "missing" ones
   */
  void removeChannels(std::vector<raw::Channel_t>& channels) const;
  
  /// Returns whether `channel` is in the configured missing channel list.
  bool isMissingChannel(raw::Channel_t channel) const;
  
  /**
   * @brief Returns the time interval where the `waveforms` should be added.
   * @param waveforms the input waveforms, grouped by channel
   * @param timeInterval the interval outside we are not interested in adding
   * @return a time interval where addition should happen
   * 
   * The `waveforms` are specified as lists of pointers to waveform information
   * records, each list pertaining a single channel.
   * The information used from the record is the timestamp and the number of
   * samples, both learnt from the optical waveform object.
   * 
   * The returned interval is guaranteed to be fully covered with data from all
   * the channels in the input, and it never extends beyond `timeInterval`.
   * It is possible that there is no sample available at all, in which case the
   * returned interval is empty (starting and ending at `timeInterval.start`).
   */
  TimeInterval_t computeTimeInterval(
    std::vector<std::vector<WaveformWithBaseline const*>> const& waveforms,
    TimeInterval_t const& timeInterval
    ) const;
  
  
  /**
   * @brief Returns the samples of a waveform added from the specified ones.
   * @param timeInterval the interval to fill
   * @param waveforms the list of waveforms (plus baseline) to include
   * @param channel (default: max possible value) channel number to assign
   * @return the added waveform and the list of the contributing waveforms
   * 
   * The returned waveform exactly covers the `timeInterval`.
   * The input waveforms are assumed to be grouped by channel: each entry in
   * `waveforms` is assumed to be the list of pointers to information records
   * for all the relevant waveform on the same channel (empty lists are ignored,
   * and there is no way to know which channel they might have belonged to).
   * The returned list includes all and only the waveforms which have at least
   * one sample in the requested time interval.
   * 
   * It is assumed that for each channel _one_ of the input `waveforms` covers
   * the whole `timeInterval`.
   */
  std::pair<std::valarray<float> , std::vector<WaveformWithBaseline const*>>
  addWaveformsInInterval(
    TimeInterval_t const& timeInterval,
    std::vector<std::vector<WaveformWithBaseline const*>> const& waveforms,
    raw::Channel_t channel = std::numeric_limits<raw::Channel_t>::max()
    ) const;

  /// Packs the specified `samples` at `startTime` into a `raw::OpDetWaveform`.
  template <typename Samples>
  raw::OpDetWaveform packWaveform
    (raw::Channel_t channel, electronics_time startTime, Samples const& samples)
    const;

  /// Returns the interval covered by the `waveform`, in electronics time.
  TimeInterval_t waveformTimeCoverage(raw::OpDetWaveform const& waveform) const;
  
  /**
   * @brief Returns the PMT channel composition of each adder window.
   * @param geom LArSoft geometry service provider
   * @return a list of adder windows, each with it list of contributing channels
   * 
   * Windows are defined as groups of 15 channels.
   * The `MissingChannels` are not included in any window.
   */
  icarus::trigger::TriggerWindowDefs_t createAdderWindows
    (geo::GeometryCore const& geom) const;
  
  /// Returns a map channel -> waveforms on that channel.
  /// If `waveformInfo` changes, the returned map becomes invalid.
  static WaveformsByChannel_t groupByChannel
    (std::vector<WaveformWithBaseline> const& waveformInfo);
  
  /// Returns the number of optical ticks from `start` to `stop`.
  template <typename StartTime, typename StopTime>
  std::ptrdiff_t tickDistance(StartTime start, StopTime time) const;
  
  /// Returns a crude baseline estimation at the first `interval` and at the
  /// last one. The baseline is just the arithmetic average.
  /// This method takes the `begin` and `end` iterator to the waveform data.
  template <typename BIter, typename EIter>
  std::pair<double, double> computeSimpleBaselines
    (BIter begin, EIter end, microseconds interval) const;
  
  /// Returns a crude baseline estimation at the first `interval` and at the
  /// last one. The baseline is just the arithmetic average.
  /// This method takes the whole `waveform` data collection.
  template <typename Waveform>
  std::pair<double, double> computeSimpleBaselines
    (Waveform const& added, microseconds interval) const;
  
  /// Converts an ADC count into voltage of a 14 bit ADC, 2 V peak-to-peak.
  template <typename T>
  static constexpr util::quantities::millivolt to_mV
    (util::quantities::counts_as<T> ADC)
    { return (double(ADC.value()) / (1<<14)) * util::quantities::volt{ 2.0 }; }
  
}; // icarus::trigger::DiscriminatedAdderSignal



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Appends a copy of all the elements of `src` into `dest`. Returns `dest`.
  template <typename DestColl, typename SrcColl>
  DestColl& append(DestColl& dest, SrcColl const& src);
  
  /// Moves all the elements of `src` at the end of `dest`. Returns `dest`.
  template <typename DestColl, typename SrcColl>
  DestColl& append(DestColl& dest, SrcColl&& src);
  
  /// Returns a `std::vector` with the results of `op` on each of `src`.
  // C++20: this is probably a single range operation
  template <typename SrcColl, typename Func>
  auto transformColl(SrcColl const& src, Func&& op)
    {
      using Result_t = std::decay_t
        <decltype(op(std::declval<typename SrcColl::value_type>()))>;
      using DestColl = std::vector<Result_t>;
      DestColl dest;
      std::transform(begin(src), end(src), back_inserter(dest), op);
      return dest;
    } // transformColl()
  
  
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& obj)
    { return std::make_unique<T>(std::move(obj)); }
  
  /// Returns a copy of the collection `coll`, sorted by `std::sort()`.
  template <typename Coll>
  Coll sorted(Coll coll)
    {
      Coll sortedColl{ coll };
      std::sort(begin(sortedColl), end(sortedColl));
      return sortedColl;
    }
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::trigger::DiscriminatedAdderSignal
//------------------------------------------------------------------------------
std::string const
icarus::trigger::DiscriminatedAdderSignal::AdderDataProductSuffix{ "adder" }; 
std::string const
icarus::trigger::DiscriminatedAdderSignal::PMTDataProductSuffix{ "" };


//------------------------------------------------------------------------------
icarus::trigger::DiscriminatedAdderSignal::DiscriminatedAdderSignal
  (Parameters const& config, art::ProcessingFrame const& frame)
  : art::SharedProducer{ config }
  // configuration
  , fWaveformTag       { config().WaveformTag()      }
  , fBaselineTag       { config().BaselineTag()      }
  , fWaveformMetaTag   { config().WaveformMetaTag().value_or(fWaveformTag) }
  , fAmplitudeScale    { config().AmplitudeScale()   }
  , fTimeInterval
    {
      icarus::ns::fhicl::makeTimeInterval(config().TimeInterval())
        .value_or(RelTimeInterval_t{})
    }
  , fMissingChannels   { sorted(config().MissingChannels()) }
  , fSaveWaveforms     { config().SaveWaveforms()    }
  , fSavePMTcoverage   { config().SavePMTcoverage()  }
  , fSaveMetaAssns     { config().SaveMetaAssns()    }
  , fAdderChannelOffset{ config().ChannelOffset()    }
  , fChannelsPerAdder  { config().ChannelsPerAdder() }
  , fLogCategory       { config().LogCategory()      }
  , fTriggerGateBuilder{
      art::make_tool<icarus::trigger::TriggerGateBuilder>
        (config().TriggerGateBuilder_.get<fhicl::ParameterSet>())
    }
  // service cache
  , fGeom              { *lar::providerFrom<geo::Geometry>() }
  , fOpticalTick       {
    detinfo::makeDetectorClocksWithUnits
      (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob())
      .OpticalClockPeriod()
    }
{
  async<art::InEvent>();
  
  //
  // optional configuration parameters
  //
  if (fTimeInterval.empty()) {
    throw art::Exception{ art::errors::Configuration }
      << "The '" << config().TimeInterval.name()
      << "' parameters configures an empty time interval [ "
      << fTimeInterval.start << " -- " << fTimeInterval.stop << " ]\n";
  }
  
  //
  // configuration post-processing
  //
  std::vector<raw::ADC_Count_t> selectedThresholds;
  if (!config().SelectThresholds(selectedThresholds)) {
    std::vector<ADCCounts_t> const& allThresholds
      = fTriggerGateBuilder->channelThresholds();
    
    if (allThresholds.empty()) {
      throw art::Exception{ art::errors::Configuration }
        << "Trigger building algorithm reports no threshold!\n"
        << "Check the configuration of `TriggerGateBuilder` tool.\n";
    }
    selectedThresholds = transformColl
      (allThresholds, std::mem_fn(&icarus::trigger::ADCCounts_t::value));
  }
  else if (selectedThresholds.empty()) {
    throw art::Exception{ art::errors::Configuration }
      << "Configured to save no thresholds!\n"
      << "Add values to the `SelectThresholds` configuration parameter.\n";
  }
  
  for (auto const& threshold: selectedThresholds) {
    fSelectedThresholds[icarus::trigger::ADCCounts_t{threshold}]
      = util::to_string(threshold);
  }
  
  {
    mf::LogInfo log(fLogCategory);
    log << "Configuration:"
      << "\n * time interval: [ "
        << fTimeInterval.start << " - " << fTimeInterval.stop << " ]"
      << "\n * selected thresholds:"
      ;
    for (auto const& [ threshold, instance ]: fSelectedThresholds) {
      log << " " << threshold << " ADC = " << to_mV(threshold)
        << " (\"" << instance << "\")";
    }
  } // nameless block
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  consumes<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>
    (fBaselineTag);
  if (fSaveMetaAssns) {
    consumes<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
      (fWaveformMetaTag);
  }
  
  //
  // declaration of output
  //
  
  if (fSaveWaveforms) {
    produces<std::vector<raw::OpDetWaveform>>();
    produces<std::vector<icarus::WaveformBaseline>>();
    produces<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>();
  }
  if (fSavePMTcoverage) {
    produces<std::vector<sbn::OpDetWaveformMeta>>();
    if (fSaveWaveforms)
      produces<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>();
  }
  for (std::string const& instanceName: util::const_values(fSelectedThresholds))
  {
    produces<std::vector<TriggerGateData_t>>(instanceName);
    produces<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
      (instanceName + PMTDataProductSuffix);
    if (fSaveMetaAssns) {
      produces<art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta>>
        (instanceName + PMTDataProductSuffix);
    }
    if (fSaveWaveforms) {
      produces<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
        (instanceName + AdderDataProductSuffix);
    } // if saving waveforms
    if (fSavePMTcoverage) {
      produces<art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta>>
        (instanceName + AdderDataProductSuffix);
    } // if saving waveform metadata
  } // for thresholds
  
} // icarus::trigger::DiscriminatedAdderSignal::DiscriminatedAdderSignal()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatedAdderSignal::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // set up the algorithm to create the trigger gates
  //
  detinfo::DetectorTimings const detTimings {
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  fTriggerGateBuilder->resetup(detTimings);
  mf::LogDebug{ fLogCategory }
    << "Trigger at " << detTimings.TriggerTime() << ", "
    << (detTimings.TriggerTime() - detTimings.BeamGateTime())
    << " after beam gate opened at " << detTimings.BeamGateTime();
  
  //
  // prepare input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  auto const& waveforms = *waveformHandle;
  
  // retrieve the baseline information
  art::FindOne<icarus::WaveformBaseline> const waveformBaselines
    { waveformHandle, event, fBaselineTag };
  
  // attach to each waveform additional information: baseline
  std::vector<WaveformWithBaseline> waveformInfo;
  waveformInfo.reserve(waveforms.size());
  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
    cet::maybe_ref<icarus::WaveformBaseline const> const baseline
      { waveformBaselines.at(iWaveform) };
    if (!baseline) {
      throw art::Exception{ art::errors::NotFound }
        << "No baseline ('" << fBaselineTag.encode()
        << "') found for waveform #" << iWaveform << " (channel="
        << waveform.ChannelNumber() << " time=" << waveform.TimeStamp()
        << ").\n";
    }
    waveformInfo.emplace_back(&waveform, &(baseline.ref()));
  }
  
  // returns a art pointer to the waveform of `wi`
  auto PMTwaveformInfoPtr
    = [handle=waveformHandle, firstWF=&(waveformHandle->front())]
      (WaveformWithBaseline const* wi)
    {
      return art::Ptr<raw::OpDetWaveform>(handle, wi->waveformPtr() - firstWF);
    };
  
  mf::LogDebug{ fLogCategory }
    << "Addition of " << waveforms.size() << " PMT waveforms from '"
    << fWaveformTag.encode() << "' with baselines from '"
    << fBaselineTag.encode() << "'";
  
  //
  // define the time interval
  //
  TimeInterval_t const timeInterval{
    detTimings.BeamGateTime() + fTimeInterval.start,
    detTimings.BeamGateTime() + fTimeInterval.stop
    };
  
  //
  // define the adder windows in terms of contributing channels;
  // this implicitly defines the channel number
  //
  icarus::trigger::TriggerWindowDefs_t const adderChannels
    = createAdderWindows(fGeom);
  
  //
  // create adder waveforms
  //
  // We need to create it in a format that `TriggerGateBuilder::build()`
  // is happy with, that is a pointer to waveform plus pointer to baseline;
  // even if the baselines are all `0`, we need them in a stable data object
  
  WaveformsByChannel_t const waveformInfoByChannel
    = groupByChannel(waveformInfo);
  
  auto const isMissing = [&missing=fMissingChannels](raw::Channel_t channel)
    { return std::binary_search(missing.cbegin(), missing.cend(), channel); };
  std::vector<raw::OpDetWaveform> adderWaveforms;
  std::vector<icarus::WaveformBaseline> adderWaveformBaselines; // quite a dummy
  std::vector<WaveformWithBaseline> adderWaveformInfo;
  std::vector<std::vector<art::Ptr<raw::OpDetWaveform>>>
    adderWaveformContribPtrs;
  // preallocate all the space we need
  adderWaveforms.reserve(adderChannels.size());
  adderWaveformBaselines.reserve(adderChannels.size());
  for (auto const& [ window, windowChannels ]: util::enumerate(adderChannels)) {
    
    raw::Channel_t const adderChannel = fAdderChannelOffset + window;
    
    mf::LogTrace log{ fLogCategory };
    log << "Collecting waveforms for adder channel " << adderChannel << " from "
      << windowChannels.size() << " channels on window " << window << ":";
    std::vector<std::vector<WaveformWithBaseline const*>> inputWaveforms;
    for (raw::Channel_t const ch: windowChannels) {
      if (isMissing(ch)) {
        // these should have been already filtered out by `createAdderWindows()`
        log << " (!" << ch << ")";
        continue;
      }
      inputWaveforms.push_back(waveformInfoByChannel[ch]);
      log << " " << ch;
      if (auto nCh = waveformInfoByChannel[ch].size(); nCh != 1)
        log << " [x" << nCh << "]";
    }
    
    auto [ waveform, waveformContribs ] = makeAdderWaveform
      (timeInterval, inputWaveforms, adderChannel);
    
    adderWaveforms.push_back(std::move(waveform));
    adderWaveformBaselines.emplace_back(0);
    // the preallocation of the two vectors guarantees that the pointers
    // are not going to be invalidated by a vector resize
    adderWaveformInfo.emplace_back
      (&adderWaveforms.back(), &adderWaveformBaselines.back());
    
    std::vector<art::Ptr<raw::OpDetWaveform>> waveformContribPtrs
      = transformColl(waveformContribs, PMTwaveformInfoPtr);
    
    adderWaveformContribPtrs.push_back(std::move(waveformContribPtrs));
    
  } // for adder windows
  
  //
  // discrimination (all thresholds at once)
  //
  
  // this is a collection where each entry (of type `TriggerGates`) contains
  // the complete set of trigger gates for an event for a given threshold.
  std::vector<icarus::trigger::TriggerGateBuilder::TriggerGates>
    triggerGatesByThreshold = fTriggerGateBuilder->build(adderWaveformInfo);
  
  using TriggerGates_t = icarus::trigger::TriggerGateBuilder::TriggerGates;
  
  { // nameless block
    mf::LogTrace log(fLogCategory);
    log << "Trigger gates from " << triggerGatesByThreshold.size()
      << " thresholds:\n";
    for (TriggerGates_t const& triggerGates: triggerGatesByThreshold)
      triggerGates.dump(log);
    
  } // nameless block
  
  // fix the channels in the gates: they contain the adder waveform channel,
  // we want the contributing PMT channels
  // (also because that will be needed when simulating a trigger)
  std::vector<std::pair<ADCCounts_t, std::vector<TriggerGateData_t>>>
    adderGatesAndThresholds;
  adderGatesAndThresholds.reserve(triggerGatesByThreshold.size());
  for (TriggerGates_t& triggerGates: triggerGatesByThreshold) {
    std::vector<TriggerGateData_t> thresholdAdders;
    for(
      auto const& [ adderWaveform, contribs ]
      : util::zip(adderWaveforms, adderWaveformContribPtrs)
    ) {
      // must already exist:
      assert(triggerGates.getGateFor(adderWaveform.ChannelNumber()));
      
      // tracking is used to keep a record of contributing waveforms from the
      // same channel; this does not apply here because only the adder "channel"
      // is contributing, and only with a single waveform (if we wished to
      // emulate trigger primitives from adders, it would be a  different story)
      // ... so we drop tracking without looking back
      
      TriggerGateData_t gate; // no channel list yet  // ... so we drop it
      // this assigns only the gate data ("levels"), not affecting its channels
      gate = std::move(triggerGates.gateFor(adderWaveform).gate().gateLevels());
      // we want to set the list of channels in the trigger gate
      // to the one of the contributing PMT channels (not the adder "channel")
      for (art::Ptr<raw::OpDetWaveform> const& waveform: contribs)
        gate.addChannel(waveform->ChannelNumber());
      
      thresholdAdders.push_back(std::move(gate));
    } // for discriminated adder
    
    adderGatesAndThresholds.emplace_back
      (triggerGates.threshold(), std::move(thresholdAdders));
    
  } // for threshold
  
  
  //
  // prepare output
  //
  auto const makeAdderPtr = fSaveWaveforms
    ? std::optional<art::PtrMaker<raw::OpDetWaveform>>{ event }: std::nullopt;
  auto const makeAdderMetaPtr = fSavePMTcoverage
    ? std::optional<art::PtrMaker<sbn::OpDetWaveformMeta>>{ event }
    : std::nullopt;
  
  // more input: needed only to save associations to input waveform metadata
  std::optional<art::FindOneP<sbn::OpDetWaveformMeta>> waveformToMeta
    = fSaveMetaAssns
    ? std::make_optional<art::FindOneP<sbn::OpDetWaveformMeta>>
      (waveformHandle, event, fWaveformMetaTag)
    : std::nullopt
    ;
  for (auto& [ thr, discrGates ]: adderGatesAndThresholds) {
    
    // find the threshold and its label, if we care of it
    auto const iThr = fSelectedThresholds.find(thr);
    if (iThr == fSelectedThresholds.end()) continue; // not interested
    
    std::string const& instanceName = iThr->second;
    art::PtrMaker<TriggerGateData_t> const makeGatePtr(event, instanceName);
    
    art::Assns<TriggerGateData_t, raw::OpDetWaveform> discrToAdderContrib;
    art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta>
      discrToAdderContribMeta;
    art::Assns<TriggerGateData_t, raw::OpDetWaveform> discrToAdder;
    art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta> discrToAdderMeta;
    for(
      auto const& [ iAdder, contribs ]
      : util::enumerate(adderWaveformContribPtrs)
    ) {
      art::Ptr<TriggerGateData_t> const gatePtr{ makeGatePtr(iAdder) };
      
      // associations: discriminated adder waveforms <=> input waveforms
      discrToAdderContrib.addMany(gatePtr, contribs);
      // associations: discriminated adder waveforms <=> input waveform metadata
      if (fSaveMetaAssns) {
        assert(waveformToMeta);
        for (art::Ptr<raw::OpDetWaveform> const& contribPtr: contribs) {
          discrToAdderContribMeta.addSingle
            (gatePtr, waveformToMeta->at(contribPtr.key()));
        }
      }
      
      // associations: discriminated adder waveforms <=> adder waveforms
      if (fSaveWaveforms) {
        assert(makeAdderPtr);
        discrToAdder.addSingle(gatePtr, (*makeAdderPtr)(iAdder));
      }
      
      // associations: discriminated adder waveforms <=> adder waveform metadata
      if (fSavePMTcoverage) {
        assert(makeAdderMetaPtr);
        discrToAdderMeta.addSingle(gatePtr, (*makeAdderMetaPtr)(iAdder));
      }
      
    } // for adder waveforms
    
    event.put(
      moveToUniquePtr(discrToAdderContrib), instanceName + PMTDataProductSuffix
      );
    
    if (fSaveMetaAssns) {
      event.put(
        moveToUniquePtr(discrToAdderContribMeta),
        instanceName + PMTDataProductSuffix
        );
    }
    
    if (fSaveWaveforms) {
      event.put
        (moveToUniquePtr(discrToAdder), instanceName + AdderDataProductSuffix);
    }
    
    if (fSavePMTcoverage) {
      event.put(
        moveToUniquePtr(discrToAdderMeta), instanceName + AdderDataProductSuffix
        );
    }

    // discriminated adder waveforms
    event.put(moveToUniquePtr(discrGates), instanceName);
    
  } // for all extracted thresholds
  
  if (fSavePMTcoverage) {
    assert(makeAdderMetaPtr);
    
    if (fSaveWaveforms) {
      assert(makeAdderPtr);
      
      art::PtrMaker<icarus::WaveformBaseline> const makeAdderBaselinePtr
        { event };
      
      // associations: adder waveforms <=> adder waveform metadata,
      //                               <=> adder waveform baselines
      art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> adderWaveformAssns;
      art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>
        adderWaveformBaselineAssns;
      for (std::size_t const i: util::counter(adderWaveforms.size())) {
        art::Ptr<raw::OpDetWaveform> const adderPtr = (*makeAdderPtr)(i);
        adderWaveformAssns.addSingle(adderPtr, (*makeAdderMetaPtr)(i));
        adderWaveformBaselineAssns.addSingle(adderPtr, makeAdderBaselinePtr(i));
      }
      event.put(moveToUniquePtr(adderWaveformBaselines));
      event.put(moveToUniquePtr(adderWaveformAssns));
      event.put(moveToUniquePtr(adderWaveformBaselineAssns));
      
    } // if also save adder waveforms
    
    // adder waveform metadata
    std::vector<sbn::OpDetWaveformMeta> adderWaveformMeta
     = transformColl(adderWaveforms, sbn::OpDetWaveformMetaMaker{ detTimings });
    assert(adderWaveformMeta.size() == adderWaveforms.size());
    event.put(moveToUniquePtr(adderWaveformMeta));
    
  } // if save metadata
  
  if (fSaveWaveforms) {
    // adder waveforms
    event.put(moveToUniquePtr(adderWaveforms));
  }
  
} // icarus::trigger::DiscriminatedAdderSignal::produce()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatedAdderSignal::removeChannels
  (std::vector<raw::Channel_t>& channels) const
{
  auto isMissing
    = [this](raw::Channel_t channel){ return isMissingChannel(channel); };
  auto itEnd = std::remove_if(channels.begin(), channels.end(), isMissing);
  channels.erase(itEnd, channels.end());
} // icarus::trigger::DiscriminatedAdderSignal::removeChannels()


//------------------------------------------------------------------------------
bool icarus::trigger::DiscriminatedAdderSignal::isMissingChannel
  (raw::Channel_t channel) const
{
  return std::binary_search
    (fMissingChannels.cbegin(), fMissingChannels.cend(), channel);
}


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::makeAdderWaveform(
  TimeInterval_t const& timeInterval,
  std::vector<std::vector<WaveformWithBaseline const*>> const& waveforms,
  raw::Channel_t channel /* = std::numeric_limits<raw::Channel_t>::max() */
) const
  -> std::pair<raw::OpDetWaveform, std::vector<WaveformWithBaseline const*>>
{
  
  // decide the time range: largest sub-interval of timeInterval where all the
  // input waveforms have information.
  TimeInterval_t const adderTimeInterval
    = computeTimeInterval(waveforms, timeInterval);
  auto [ added, contribs ]
    = addWaveformsInInterval(adderTimeInterval, waveforms, channel);
  
  raw::OpDetWaveform waveform
    = packWaveform(channel, timeInterval.start, added);
  
  {
    auto const [ startBaseline, endBaseline ]
      = computeSimpleBaselines(added, 500_ns);
    float const peak
      = startBaseline - *std::min_element(begin(added), end(added));
    mf::LogTrace log{ fLogCategory };
    log << "Adder waveform CH=" << waveform.ChannelNumber()
      << " created from " << contribs.size() << " waveforms: "
      << waveform.size() << " samples starting at timestamp="
      << waveform.TimeStamp();
    if (!waveform.empty()) {
      auto const [ itMinADC, itMaxADC ]
        = std::minmax_element(waveform.cbegin(), waveform.cend());
      electronics_time const minTime = timeInterval.start
        + fOpticalTick * std::distance(waveform.cbegin(), itMinADC);
      log << ", range [ " << *itMaxADC << " to " << *itMinADC
        << " ] ADC (peak at " << minTime << ")";
    }
    log << "\n  adder baseline at start " << startBaseline
      << " ADC, at end " << endBaseline << " ADC, peak: "
      << peak << " ADC = " << to_mV(util::quantities::counts_f{ peak });
  }
  
  return { std::move(waveform), std::move(contribs) };
  
} // icarus::trigger::DiscriminatedAdderSignal::makeAdderWaveform()


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::computeTimeInterval(
  std::vector<std::vector<WaveformWithBaseline const*>> const& waveforms,
  TimeInterval_t const& timeInterval
) const -> TimeInterval_t {
  /*
   * If we are willing to make up data (e.g. assigning samples of value `0` to
   * the ticks that are not covered by any waveform in a channel), we can then
   * exactly will the required `timeInterval`, with the caveat that we need
   * still to supply the limits that the `timeInterval` does not specify
   * (e.g. if `timeInterval.end` is "infinite", we may replace it with the last
   * waveform time end).
   * 
   * If we are not willing to make up data, then we need to trim the allowed
   * interval so that all channels have data in it. This is a more complicate
   * choice in that it requires an "and" of all the coverage in each channel,
   * in conjunction with the required `timeInterval`, and this and may also end
   * up being not contiguous, in which case a choice must be made since only a
   * single contiguous interval is accepted (with the premise above that we
   * don't want to fill the gaps with made up data).
   */
  using TimeCoverage_t = icarus::trigger::TriggerGateData
    <electronics_time, detinfo::timescales::time_interval>;
  
  //
  // 1. determine all the time sub-intervals within `timeInterval` covered by
  //    all channels
  //
  TimeCoverage_t coverage;
  coverage.setOpeningAt(TimeCoverage_t::MinTick, 1U); // start always open

  mf::LogTrace{ fLogCategory } << "Coverage within " << timeInterval
    << " (only channels narrowing it are listed):";
  for (std::vector<WaveformWithBaseline const*> const& chWaveforms: waveforms) {
    
    TimeCoverage_t chCoverage;
    
    for (WaveformWithBaseline const* wi: chWaveforms) {
      TimeInterval_t wInterval = waveformTimeCoverage(wi->waveform());
      assert(!wInterval.empty());
      
      // trim to the allowed` timeInterval`; we don't assume `wi` time ordering
      wInterval.intersect(timeInterval);
      if (wInterval.empty()) continue;
      
      chCoverage.openBetween(wInterval.start, wInterval.stop);
    } // for waveforms in channel
    
    if (chWaveforms.empty()) {
      mf::LogTrace{ fLogCategory } << " - no waveform in channel";
    }
    else if (chCoverage != coverage) {
      mf::LogTrace{ fLogCategory } << " - channel "
        << chWaveforms.front()->waveform().ChannelNumber()
        << " with " << chWaveforms.size() << " waveforms => " << chCoverage;
    }
    coverage = icarus::trigger::GateOps::Mul(coverage, chCoverage);
    
  } // for channels
  mf::LogTrace{ fLogCategory } << "Total coverage from " << waveforms.size()
    << " channels: " << coverage;
  
  //
  // 2. if there are multiple non-contiguous sub-intervals, pick the widest one
  //
  auto nextWindow
    = [](TimeCoverage_t const& coverage, TimeInterval_t const& prev)
    {
      TimeInterval_t next;
      next.start = coverage.findOpen(1, prev.stop);
      if (next.start == TimeCoverage_t::MaxTick) return TimeInterval_t{};
      next.stop = coverage.findClose(1, next.start);
      return next;
    };
  
  TimeInterval_t interval{ timeInterval.start };
  auto duration = interval.duration();
  TimeInterval_t runningWindow{ TimeCoverage_t::MinTick };
  while (true) {
    runningWindow = nextWindow(coverage, runningWindow);
    if (runningWindow.empty()) break; // no more openings left
    assert(runningWindow.stop != TimeCoverage_t::MaxTick);
    
    auto const windowDuration = runningWindow.duration();
    if (duration > windowDuration) continue; // smaller window, discarded
    duration = windowDuration;
    interval = runningWindow;
  } // while
  
  return interval;
} // icarus::trigger::DiscriminatedAdderSignal::computeTimeInterval()


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::addWaveformsInInterval(
  TimeInterval_t const& timeInterval,
  std::vector<std::vector<WaveformWithBaseline const*>> const& waveforms,
  raw::Channel_t channel /* = std::numeric_limits<raw::Channel_t>::max() */
) const
  -> std::pair<std::valarray<float> , std::vector<WaveformWithBaseline const*>>
{
  
  std::size_t const nSamples = timeInterval.empty()
    ? 0: tickDistance(timeInterval.start, timeInterval.stop);
  
  std::valarray<float> added(0.0f, nSamples);
  std::vector<WaveformWithBaseline const*> used;
  
  mf::LogTrace{ fLogCategory } << "Building adder waveform CH=" << channel
    << ", " << nSamples << " sample long, at interval " << timeInterval
    << " from " << waveforms.size() << " channels";
  
  for (std::vector<WaveformWithBaseline const*> const& chWaveforms: waveforms) {
    
    for (WaveformWithBaseline const* wi: chWaveforms) {
      
      TimeInterval_t const wfCoverage = waveformTimeCoverage(wi->waveform());
      
      if (!wfCoverage.contains(timeInterval.start)) continue;
      // in this approach, we assume that one of the selected waveforms covers
      // the whole interval (if it were more, there would be a gap, and this
      // algorithm should have carved an interval where no channel has gaps)
      assert(wfCoverage.contains(timeInterval.stop));
      
      // add (minus baseline)
      std::ptrdiff_t const startSample
        = tickDistance(wfCoverage.start, timeInterval.start);
      
      auto itSample = std::next(wi->waveform().cbegin(), startSample);
      for (float& sample: added) sample += *(itSample++);
      
      if (wi->hasBaseline()) added -= wi->baseline().baseline();
      
      { // --- BEGIN -- DEBUG --------------------------------------------------
        auto const [ startBaseline, endBaseline ]
          = computeSimpleBaselines(wi->waveform(), 500_ns);
        auto const itLowest
          = std::min_element(wi->waveform().begin(), wi->waveform().end());
        raw::ADC_Count_t const lowest = *itLowest;
        electronics_time const lowestTime = wfCoverage.start
          + fOpticalTick * std::distance(wi->waveform().begin(), itLowest);
        mf::LogTrace log{ fLogCategory };
        log << " - from ch=" << std::setw(3) << wi->waveform().ChannelNumber()
          << ", " << wfCoverage
          << ", sample=" << std::setw(5) << startSample << ", baseline: ";
        if (wi->hasBaseline()) {
          double const baseline = wi->baseline().baseline();
          log << std::setw(7) << std::setprecision(6) << std::left << baseline
            << " ADC, diff of local start="
            << std::showpos << std::setprecision(2) << std::setw(7) << std::left
              << (startBaseline - baseline) << " ADC, end="
            << std::showpos << std::setprecision(2) << std::setw(7) << std::left
              << (endBaseline - baseline) << " ADC, peak amplitude: "
            << std::defaultfloat << std::setprecision(6) << std::setw(8)
            << (baseline - lowest) << " ADC at time " << lowestTime;
        }
        else {
          log << "n/a, local start=" << std::setw(7) << std::setprecision(1)
            << startBaseline << " ADC, end="
            << std::setw(7) << std::setprecision(1) << endBaseline << " ADC"
            << " ADC, lowest point: " << lowest << " ADC at time "
            << lowestTime;
        }
      } // --- END ---- DEBUG --------------------------------------------------
      
      used.push_back(wi); // track
      break; // this waveform has full coverage as per above, no more needed
    } // for waveforms in channel
    
    if (chWaveforms.empty()) {
      mf::LogTrace{ fLogCategory } << " - empty input collection!";
    }
    else if (used.empty()
      || (
        chWaveforms.front()->waveform().ChannelNumber()
        != used.back()->waveform().ChannelNumber()
      )
    ) {
      mf::LogTrace{ fLogCategory } << " - none of the " << chWaveforms.size()
        << " waveforms of channel "
        << chWaveforms.front()->waveform().ChannelNumber()
        << " falls in the requested interval";
    }
    
  } // for channels
  
  added *= fAmplitudeScale;
  
  return { std::move(added), std::move(used) };
  
} // icarus::trigger::DiscriminatedAdderSignal::addWaveformsInInterval()


//------------------------------------------------------------------------------
template <typename Samples>
raw::OpDetWaveform icarus::trigger::DiscriminatedAdderSignal::packWaveform
  (raw::Channel_t channel, electronics_time startTime, Samples const& samples)
  const
{
  using std::size, std::begin, std::end;
  
  raw::TimeStamp_t const timeStamp
    = startTime.convertInto<util::quantities::points::microsecond>().value();
  raw::OpDetWaveform waveform{
      timeStamp     /* timestamp */
    , channel       /* channel */
    , size(samples)  /* length */
    };
  
  // float -> short
  std::transform
    (begin(samples), end(samples), back_inserter(waveform), std::roundf);
  
  return waveform;
} // icarus::trigger::DiscriminatedAdderSignal::packWaveform()


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::waveformTimeCoverage
  (raw::OpDetWaveform const& waveform) const -> TimeInterval_t
{
  electronics_time const start
    { util::quantities::microsecond{ waveform.TimeStamp() } };
  return TimeInterval_t{ start, start + waveform.size() * fOpticalTick };
}


//------------------------------------------------------------------------------
icarus::trigger::TriggerWindowDefs_t
icarus::trigger::DiscriminatedAdderSignal::createAdderWindows
  (geo::GeometryCore const& geom) const
{
  icarus::trigger::SlidingWindowDefinitionAlg makeWindows{ geom };
  
  icarus::trigger::TriggerWindowDefs_t windows
    = makeWindows(fChannelsPerAdder, fChannelsPerAdder);
  
  for (TriggerWindowChannels_t& window: windows)
    removeChannels(window);
  
  return windows;
} // icarus::trigger::DiscriminatedAdderSignal::createAdderWindows()


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::groupByChannel
  (std::vector<icarus::trigger::WaveformWithBaseline> const& waveformInfo)
  -> WaveformsByChannel_t
{
  auto waveformChannel = [](icarus::trigger::WaveformWithBaseline const& wi)
    { return static_cast<std::size_t>(wi.waveform().ChannelNumber()); };
  return icarus::ns::util::GroupByIndex{ waveformInfo, waveformChannel };
} // icarus::trigger::DiscriminatedAdderSignal::groupByChannel()


//------------------------------------------------------------------------------
template <typename StartTime, typename StopTime>
std::ptrdiff_t icarus::trigger::DiscriminatedAdderSignal::tickDistance
  (StartTime start, StopTime stop) const
{
  return static_cast<std::ptrdiff_t>(std::trunc((stop - start) / fOpticalTick));
}


//------------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::DiscriminatedAdderSignal::computeSimpleBaselines
  (BIter begin, EIter end, microseconds interval) const
  -> std::pair<double, double>
{
  assert(interval >= 0_us);
  std::ptrdiff_t const nSamples
    = std::min(tickDistance(0_us, interval), std::distance(begin, end));
  return {
    std::accumulate(begin, begin + nSamples, 0.0) / nSamples,
    std::accumulate(end - nSamples, end, 0.0) / nSamples
  };
} // icarus::trigger::DiscriminatedAdderSignal::computeSimpleBaselines(Iter)


//------------------------------------------------------------------------------
template <typename Waveform>
auto icarus::trigger::DiscriminatedAdderSignal::computeSimpleBaselines
  (Waveform const& waveform, microseconds interval) const
  -> std::pair<double, double>
{
  using std::begin, std::end;
  return computeSimpleBaselines(begin(waveform), end(waveform), interval);
} // icarus::trigger::DiscriminatedAdderSignal::computeSimpleBaselines(Waveform)


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::DiscriminatedAdderSignal)


//------------------------------------------------------------------------------
