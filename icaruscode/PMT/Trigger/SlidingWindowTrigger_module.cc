/**
 * @file   SlidingWindowTrigger_module.cc
 * @brief  Combines trigger channels into sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 6, 2020
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowCombinerAlg.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h"
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()

// LArSoft libraries
#include "lardata/Utilities/NestedIterator.h" // lar::double_fwd_const_iterator
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/CoreUtils/get_elements.h" // util::get_elements()
// #include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::::static_assert_on<>

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"


// C/C++ standard libraries
#include <map>
#include <vector>
#include <string>
#include <optional>
#include <memory> // std::unique_ptr
#include <utility> // std::move()
#include <cassert>
#if 0
#include <set>
#include <cstddef> // std::size_t

#endif // 0

//------------------------------------------------------------------------------
namespace icarus::trigger { class SlidingWindowTrigger; }

/**
 * @brief Combines trigger channels (V1730 LVDS gates) into sliding windows.
 * 
 * This module combines the input discriminated waveforms into a set of new
 * multi-level discriminated gates according to the geometric criterion of a
 * sliding window.
 *
 * The combination is a sum of all the discrminated waveforms within the window.
 * Different optical detector planes are treated independently and separately.
 * It is expected that the input is from V1730 LVDS gates
 * (see `icarus::trigger::LVDSgates` module) which combine (usually) two optical
 * detector channels each.
 *
 *
 * Sliding window configuration
 * -----------------------------
 *
 * The size of the window is specified in terms on number of optical detector
 * channels per window and of stride. The windows are allowed to overlap, but
 * they have all the same number of channels.
 * If the parameters yield an invalid configuration, the module will throw an
 * exception.
 *
 * The definition of a valid window configuration is quite restrictive.
 * Optical detector channels are grouped by the TPC width coordinate (_z_) into
 * "towers". In ICARUS, a tower is made of either two or three channels,
 * depending on the position in _z_. A window can't split a tower, but rather
 * has to include whole towers.
 * In addition, the input discriminated waveforms may include information from
 * multiple channels (that is the common case if the input is LVDS gates).
 * A window can't cover only part of the channels of any input, but rather must
 * extend wide enough so that all the channels in the included gates are
 * covered.
 * Finally, all the channels must be covered by at least one window.
 *
 * For the sliding window configuration to be valid, all the windows starting
 * with an offset multiple of the stride must be valid.
 *
 * For the standard pairing of LVDS gates in ICARUS, this translates into the
 * only allowed window sizes are by 10 and 15 channels or their multiples,
 * with a stride also in multiples of 10 or 15 respectively.
 * For example, the following configurations are all valid:
 * * 10 channels per window, with stride of 10 channels (9 windows per plane);
 * * 20 channels per window, with stride of 10 channels (8 windows per plane);
 * * 30 channels per window, with stride of 10 channels (7 windows per plane);
 * * 15 channels per window, with stride of 15 channels (6 windows per plane);
 * * 30 channels per window, with stride of 30 channels (3 windows per plane);
 * * 30 channels per window, with stride of 15 channels (5 windows per plane);
 * * 45 channels per window, with stride of 45 channels (2 windows per plane).
 *
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (labels out of
 *   `TriggerGatesTag` and `Thresholds`): full sets of discriminated waveforms,
 *   each waveform possibly covering multiple optical channels,
 *   and their associations to optical waveforms. One set per threshold.
 *
 * Requirements
 * -------------
 *
 * The algorithms require that in each set all optical channels are covered.
 *
 *
 * Output data products
 * =====================
 *
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (instance name:
 *   same as the input gates): sets of gates combined according to the window
 *   configuration; one set per input threshold.
 * * `art::Assns<icarus::trigger::OpticalTriggerGateData_t, raw::OpDetWaveform>`
 *   (instance name: same as the input gates): associations between each
 *   produced gate and the optical waveforms providing the original data.
 *
 * @note At the moment, associations between input and output gates is not
 *       produced. `art::Assns` would not support it, requiring distinct
 *       types for the associated objects.
 *
 *
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description SlidingWindowTrigger`.
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the discriminated waveforms; it must not include
 *     any instance name, as the instance names will be automatically added from
 *     `Thresholds` parameter.
 * * `Thresholds` (list of integers, mandatory): list of the discrimination
 *     thresholds to consider, in ADC counts. A data product containing a
 *     digital signal is read for each one of the thresholds, and the tag of the
 *     data product is expected to be the module label `TriggerGatesTag` with as
 *     instance name the value of the threshold (e.g. for a threshold of 6 ADC
 *     counts the data product tag might be `discrimopdaq:6`).
 * * `WindowSize` (integral, mandatory): the number of _optical detector
 *     channels_ that every window must include. Note that this is not
 *     equivalent to the number of input objects from the `TriggerGatesTag` data
 *     products: for example, in the standard ICARUS configuration, a window
 *     spanning 15 channels is generated out of only 8 LVDS input gates.
 * * `Stride` (integral, default: as `WindowSize`): the base offset, in number
 *     of optical detector channels, where to start a new window. Windows are
 *     started at all multiples of `Stride`. For example, windows of 30 channels
 *     with `Stride` 15 will start after 0, 15, 30, 45 and 60 channels on each
 *     optical detector plane.
 * * `LogCategory` (string): name of the output stream category for console
 *     messages (managed by MessageFacility library).
 *
 *
 * 
 */
class icarus::trigger::SlidingWindowTrigger: public art::EDProducer {
  
    public:
  
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGate::GateData_t;
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;


    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Sequence<raw::ADC_Count_t> Thresholds {
      Name("Thresholds"),
      Comment("thresholds to consider [ADC counts]")
      };

    fhicl::Atom<unsigned int> WindowSize {
      Name("WindowSize"),
      Comment("numer of optical channels to be included in each window")
      };
    fhicl::OptionalAtom<unsigned int> Stride {
      Name("Stride"),
      Comment(
        "number of optical channel used as offset for sliding window [as WindowSize]"
        )
      };

    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "SlidingWindowTrigger" // default
      };

  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit SlidingWindowTrigger(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  SlidingWindowTrigger(SlidingWindowTrigger const&) = delete;
  SlidingWindowTrigger(SlidingWindowTrigger&&) = delete;
  SlidingWindowTrigger& operator=(SlidingWindowTrigger const&) = delete;
  SlidingWindowTrigger& operator=(SlidingWindowTrigger&&) = delete;
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Fills the plots. Also extracts the information to fill them with.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  /// Definition of all windows.
  using WindowDefs_t = icarus::trigger::SlidingWindowCombinerAlg::Windows_t;
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<icarus::trigger::ADCCounts_t, art::InputTag> fADCthresholds;

  unsigned int const fWindowSize; ///< Sliding window size in number of channels.
  unsigned int const fWindowStride; ///< Sliding window base offset.

  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  /// Channel content of each window.
  WindowDefs_t const fWindowChannels;

  /// Combining algorithm.
  icarus::trigger::SlidingWindowCombinerAlg const fCombiner;


  /// Performs the combination for data with a specified threshold.
  void produceThreshold(
    art::Event& event,
    icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap,
    icarus::trigger::ADCCounts_t const threshold,
    art::InputTag const& dataTag
    ) const;

  /// Reads a set of input gates from the `event` and updates `waveformMap`.
  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> ReadTriggerGates(
    art::Event const& event,
    art::InputTag const& dataTag,
    icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap
    ) const;

  /// Defines the channels falling in each window.
  WindowDefs_t defineWindows() const;

  /// Adds the waveforms in the specified association to the waveform `map`.
  static void UpdateWaveformMap(
    icarus::trigger::OpDetWaveformDataProductMap_t& map,
    art::Assns<TriggerGateData_t, raw::OpDetWaveform> const& assns
    );

}; // class icarus::trigger::SlidingWindowTrigger



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
icarus::trigger::SlidingWindowTrigger::SlidingWindowTrigger
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fWindowSize(config().WindowSize())
  , fWindowStride
    (util::fhicl::getOptionalValue(config().Stride).value_or(fWindowSize))
  , fLogCategory(config().LogCategory())
  , fWindowChannels(defineWindows())
  , fCombiner(fWindowChannels)
{
  //
  // more complex parameter parsing
  //
  std::string const discrModuleLabel = config().TriggerGatesTag();
  for (raw::ADC_Count_t threshold: config().Thresholds()) {
    fADCthresholds[icarus::trigger::ADCCounts_t{threshold}]
      = art::InputTag{ discrModuleLabel, util::to_string(threshold) };
  }

  //
  // configuration report (short)
  //
  {
    mf::LogInfo log(fLogCategory);
    log << fWindowChannels.size() << " windows configured:";
    for (auto const& [ iWindow, channels ]: util::enumerate(fWindowChannels)) {
      log << "\n [#" << iWindow << "] " << channels.size() << " channels:";
      for (raw::Channel_t const channel: channels) log << " " << channel;
    } // for windows
    log << "\nConfigured " << fADCthresholds.size() << " thresholds:";
    for (auto const& [ threshold, dataTag ]: fADCthresholds)
      log << "\n * " << threshold << " ADC (from '" << dataTag.encode() << "')";
  } // local block


  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  //
  // input data declaration
  //
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    consumes<std::vector<OpticalTriggerGateData_t>>(inputDataTag);
    consumes<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag);
  } // for


  //
  // output data declaration
  //
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    produces<std::vector<OpticalTriggerGateData_t>>(inputDataTag.instance());
    produces<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag.instance());
  } // for

  
} // icarus::trigger::SlidingWindowTrigger::SlidingWindowTrigger()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTrigger::produce(art::Event& event) {
  
  icarus::trigger::OpDetWaveformDataProductMap_t waveformMap;

  for (auto const& [ threshold, dataTag ]: fADCthresholds) {

    produceThreshold(event, waveformMap, threshold, dataTag);
    
  } // for all thresholds
  
} // icarus::trigger::SlidingWindowTrigger::produce()


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowTrigger::defineWindows() const
  -> WindowDefs_t
{
  /*
   * 1. compute the vertical PMT towers in each separate optical detector plane
   * 2. fill the windows by counting channels (i.e. op. det.)
   */
  using icarus::trigger::PMTverticalSlicingAlg;

  //
  // 1. compute the vertical PMT towers in each separate optical detector plane
  //
  geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
  PMTverticalSlicingAlg slicerAlg(fLogCategory);
  PMTverticalSlicingAlg::Slices_t slices;
  for (geo::CryostatGeo const& cryo: geom.IterateCryostats())
    slicerAlg.appendCryoSlices(slices, cryo);

  //
  // 2. fill the windows by counting channels (i.e. op. det.)
  //
  WindowDefs_t windows;

  for (PMTverticalSlicingAlg::PMTtowerOnPlane_t const& planeSlices: slices) {

    auto itSlice = planeSlices.begin();
    auto const send = planeSlices.end();
    while (itSlice != send) {

      mf::LogTrace(fLogCategory) << "Assembling window #" << windows.size();

      WindowDefs_t::value_type window;
      window.reserve(fWindowSize);

      std::optional<decltype(itSlice)> nextStart;
      unsigned int nChannels = 0U;
      while (nChannels < fWindowSize) {

        // aside: check if this is the right place to start the next window
        if (nChannels == fWindowStride) {
          mf::LogTrace(fLogCategory)
            << "  (next window will start from this slice)";
          nextStart = itSlice;
        }
        else if ((nChannels > fWindowStride) && !nextStart) {
          throw cet::exception("SlidingWindowTrigger")
            << "Unable to start a new window " << fWindowStride
            << " channels after window #" << windows.size()
            << " (next slice starts " << nChannels << " channels after)\n";
        }

        mf::LogTrace(fLogCategory)
          << "  adding " << itSlice->size() << " channels to existing "
          << nChannels;
        for (geo::OpDetGeo const* opDet: *itSlice) {
          geo::OpDetID const& id = opDet->ID();
          raw::Channel_t const channel
            = geom.OpDetFromCryo(id.OpDet, id.Cryostat);
          mf::LogTrace(fLogCategory)
            << "   * " << id << " (channel " << channel << ")";
          window.push_back(channel);
        } // for channels in slice
        nChannels += (itSlice++)->size();
      } // while
      if (nChannels == fWindowStride) nextStart = itSlice;
      assert(nextStart);

      if (nChannels != fWindowSize) {
        throw cet::exception("SlidingWindowTrigger")
          << "Definition of one window yielded " << nChannels
          << " elements (window should be of size " << fWindowSize
          << " and with stride " << fWindowStride << ").\n";
      }

      windows.push_back(std::move(window));

      itSlice = nextStart.value();
    } // for all slices
  } // for all windows
  mf::LogTrace(fLogCategory)
    << "SlidingWindowTrigger defined " << windows.size() << " windows.";

  return windows;
} // icarus::trigger::SlidingWindowTrigger::defineWindows()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTrigger::produceThreshold(
  art::Event& event,
  icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap,
  icarus::trigger::ADCCounts_t const threshold,
  art::InputTag const& dataTag
) const {
  
  mf::LogDebug(fLogCategory)
    << "Processing threshold " << threshold
    << " from '" << dataTag.encode() << "'";

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience
  
  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> const& gates
    = ReadTriggerGates(event, dataTag, waveformMap);
  
  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> combinedGates
    = fCombiner.combine(gates);

  // transform the data; after this line, `gates` is not usable any more
  art::PtrMaker<OpticalTriggerGateData_t> const makeGatePtr
    (event, dataTag.instance());
  auto&& [ outputGates, outputAssns ]
    = icarus::trigger::transformIntoOpticalTriggerGate
    (combinedGates, makeGatePtr, waveformMap);

  mf::LogTrace(fLogCategory)
    << "Threshold " << threshold << " ('" << dataTag.encode() << "'): "
    << gates.size() << " input channels, "
    << outputGates.size() << " output channels (+"
    << outputAssns.size() << " associations to waveforms) into '"
    << moduleDescription().moduleLabel() << ":" << dataTag.instance() << "'"
    ;

  event.put(
    std::make_unique<std::vector<OpticalTriggerGateData_t>>
      (std::move(outputGates)),
    dataTag.instance()
    );
  event.put(
    std::make_unique<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (std::move(outputAssns)),
    dataTag.instance()
    );

} // icarus::trigger::SlidingWindowTrigger::produceThreshold()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTrigger::UpdateWaveformMap(
  icarus::trigger::OpDetWaveformDataProductMap_t& map,
  art::Assns<TriggerGateData_t, raw::OpDetWaveform> const& assns
) {
  
  for (art::Ptr<raw::OpDetWaveform> const& wave: util::get_elements<1U>(assns))
    map.emplace(wave.get(), wave);
  
} // icarus::trigger::SlidingWindowTrigger::UpdateWaveformMap()


//------------------------------------------------------------------------------
std::vector<icarus::trigger::MultiChannelOpticalTriggerGate>
icarus::trigger::SlidingWindowTrigger::ReadTriggerGates(
  art::Event const& event,
  art::InputTag const& dataTag,
  icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap
) const {
  using icarus::trigger::OpticalTriggerGateData_t;

  auto const& assns =
    *(event.getValidHandle<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>(dataTag));

  UpdateWaveformMap(waveformMap, assns);

  return icarus::trigger::FillTriggerGates<icarus::trigger::MultiChannelOpticalTriggerGate>
    (*(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>(dataTag)), assns);

} // icarus::trigger::SlidingWindowTrigger::ReadTriggerGates()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::SlidingWindowTrigger)


//------------------------------------------------------------------------------
