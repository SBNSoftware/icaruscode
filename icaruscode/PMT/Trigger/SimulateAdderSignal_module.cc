/**
 * @file   SimulateAdderSignal_module.cc
 * @brief  Executes the simulation of adder board output from PMT signals.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 1, 2025
 */

#undef NDEBUG // FIXME

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulationTypes.h" // adder::types
#include "icaruscode/PMT/Trigger/Algorithms/AmplitudeAdderCalibration.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.h"
#include "icaruscode/PMT/Trigger/Algorithms/SallenKeyFilter.h"
#include "icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelID.h"
#include "icaruscode/PMT/Trigger/Algorithms/WaveformWithBaseline.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icaruscode/IcarusObj/ChannelToChannelMap.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "sbnalg/Utilities/ToStr.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/SharedResource.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib/maybe_ref.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/TableFragment.h"
#include "fhiclcpp/types/OptionalTableAs.h"
#include "fhiclcpp/types/OptionalTable.h"

// ROOT libraries
#include "TDirectoryFile.h"

// // C/C++ standard libraries
#include <algorithm> // std::sort(), std::max()
#include <cassert>
#include <memory> // std::unique_ptr
#include <ostream>
#include <string>
#include <type_traits> // std::remove_reference_t
#include <utility> // std::move(), std::pair
#include <vector>


//------------------------------------------------------------------------------
namespace icarus::trigger { class SimulateAdderSignal; }
/**
 * @brief Simulates the waveform output of ICARUS adder boards.
 * 
 * This module produces a waveform for each adder module.
 * 
 * The input is PMT optical waveforms: this module composes them into added
 * analogue signals, reshapes and finally "calibrates" them.
 * 
 * 
 * Missing channels
 * -----------------
 * 
 * The module allows for the specification of a list of missing PMT channels.
 * A missing channel is in general a channel whose output should be ignored.
 * In practice, depending on the type of input and of channel issue, the channel
 * may be completely absent from the input data product, be present with only
 * noise, normal or high (typical of detector data), or be present with physical
 * activity (in simulation where these channels are not masked).
 * 
 * The correct treatment of these channels depends on the input and this module
 * is not spending time in figuring it out.
 * The current implementation is that waveforms from a "missing" channels, as
 * flagged by the configuration parameter `MissingChannels`, are used if present
 * but not required (as opposed to throwing an exception when an expected and
 * good channel of an adder is not found).
 * 
 * 
 * On calibration
 * ---------------
 * 
 * The calibration is bridging the gap between the first principle, approximate
 * simulation and the observed data.
 * 
 * While an interface exists for calibration algorithms, only one is currently
 * implemented and hard-plugged: `icarus::trigger::AmplitudeAdderCalibration`.
 * This calibration just rescales the whole waveform so that the peak amplitude
 * reflects the one in data.
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * Enabling the writing of plots into the ROOT file with `art::TFileService`
 * (`DebugPlots` configuration parameter) explicitly inhibits multithreading.
 * 
 * 
 * Output
 * =======
 * 
 * Data products
 * --------------
 * 
 * In each run:
 * * `icarus::ChannelToChannelMap<raw::Channel_t>` (if `WriteChannelMap` is
 *   set): the mapping between adder and PMT channels.
 * 
 * In each event:
 * * `std::vector<raw::OpDetWaveform>`:
 *   waveforms from the adders; in this emulation, the waveforms are
 *   "reproduced" from the already digitized PMT waveforms, in contrast with
 *   the hardware which performs an analogue sum.
 *   The channel number is unique on each of them, and it is following a
 *   convention documented in `icarus::trigger::AdderChannelMaps`.
 *   Waveforms are sorted by channel, and within each channel by starting time.
 * * `std::vector<icarus::WaveformBaseline>`: the adder baseline; it is `0` by
 *   construction.
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>`: the association
 *   between each adder waveform and its baseline object (which is trivial since
 *   the two data products are
 *   @ref LArSoftProxyDefinitionParallelData "parallel").
 * * `std::vector<sbn::OpDetWaveformMeta>`: parallel to the
 *   `std::vector<raw::OpDetWaveform>` data product, each entry in the
 *   collection is matched with the corresponding one in the adder waveform
 *   collection; an association
 *   `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>` is also produced
 *   (which is also trivial since the two data products are
 *   @ref LArSoftProxyDefinitionParallelData "parallel").
 * 
 * 
 * Plots
 * ------
 * 
 * If `DebugPlots` is enabled, intermediate waveforms are saved in the file
 * managed by `TFileService` as ROOT `TGraph` objects.
 * The plots are organised as usual under a directory named after the label of
 * this module, and then in a few directory levels, the top-most being the
 * event ID (run and event number). See the section "Plots" of
 * `icarus::trigger::AdderChannelSimulator` documentation for further details.
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
 * * an implementation of `detinfo::DetectorClocksService`
 * * `IICARUSChannelMap` (ICARUS hardware channel mapping service)
 * * _art_ message facility
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description SimulateAdderSignal`.
 * 
 * * `WaveformTag` (input tag): the data product containing all optical
 *   detector waveforms.
 * * `BaselineTag` (input tag): the data product containing the association
 *   between each input waveform (`WaveformTag`) and its estimated baseline.
 * * `MissingChannels` (list of integers): list of the PMT channel numbers
 *   (as reported in the `WaveformTag` input collection) flagged as "missing".
 *   Note that in the hardware bad channels may be included in the adder input.
 * * `SplitToAdders` (real, default: 5%): the fraction of PMT signal sent to the
 *   adder boards. It is assumed that the PMT signal in `WaveformTag` received
 *   the complement of this.
 * * `ReshapingCircuitParameters`: configuration of the reshaping filter.
 *   Currently the choice of filter is hard-coded to be a Sallen-Key filter; see
 *   `icarus::trigger::SallenKeyFilter` for the configuration details.
 *   If omitted entirely (default), no reshaping is performed, which is only
 *   useful for debugging the PMT summation code.
 * * `Calibration`: configuration of the calibration algorithm. Currently the
 *   choice of calibration algorithm is hard-coded: see
 *   `icarus::trigger::AmplitudeAdderCalibration` for the configuration details.
 *   If omitted entirely (default), no calibratio is applied, which is useful
 *   for extracting one.
 * * `DebugPlots` (flag, default: `false`): write ROOT plots of each
 *   intermediate waveform in the simulation (PMT sum, its reshaping, and the
 *   final, calibrated result); see the `Plots` section above.
 *   Do not activate this flag in production jobs: it can produce a very
 *   sizeable output.
 * * `WriteChannelMap` (flag, default: `false`): if set, the adder channel
 *   mapping data object will be written into the run data. This module always
 *   creates a mapping anew (i.e. it does not read it from the run even if it
 *   would be available). Another way to write the mapping is to execute the
 *   module `WriteAdderChannelMap`.
 * * `LogCategory` (string, default: `"SimulateAdderSignal"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 * 
 */
class icarus::trigger::SimulateAdderSignal
  : public art::SharedProducer, private icarus::ns::util::mfLoggingClass
{
  
    public:
  
  /// Type of channel map optionally serialized.
  using ChannelMap_t = icarus::ChannelToChannelMap<raw::Channel_t>;

  /// FHiCL configuration of the module.
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
    
    fhicl::TableFragment<AdderSignalSimulation::Config> simulation;
    
    fhicl::OptionalTableAs<SallenKeyFilter::Parameters, SallenKeyFilterConfig>
    ReshapingCircuitParameters {
      Name{ "ReshapingCircuitParameters" },
      Comment{ "parameters of the reshaping (Sallen-Key filter) circuit" }
      };
    
    // should become a fhicl::DelegatedParameter
    fhicl::OptionalTable<icarus::trigger::AmplitudeAdderCalibration::Config>
    Calibration{
      Name{ "Calibration" },
      Comment{ "configuration of the chosen adder calibration database" }
      };
    
    fhicl::Atom<bool> DebugPlots {
      Name{ "DebugPlots" },
      Comment{ "write ROOT plots of each intermediate waveform in the simulation" },
      false
      };
    
    fhicl::Atom<bool> WriteChannelMap {
      Name{ "WriteChannelMap" },
      Comment{ "write the adder channel map in a run-level data product" },
      false
      };
    
  }; // Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  
  SimulateAdderSignal(Parameters const& params, art::ProcessingFrame const&);
  
  /// Loads run-specific database information.
  virtual void beginRun(art::Run& run, art::ProcessingFrame const&) override;
  
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  
    private:
  
  // imported types
  using microseconds = util::quantities::intervals::microseconds;
  using electronics_time = detinfo::timescales::electronics_time;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fWaveformTag; ///< Input optical waveform tag.
  
  art::InputTag const fBaselineTag; ///< Input waveform baseline tag.
  
  bool const fDebugPlots; ///< Whether to write debugging plots.
  
  /// Whether to put the adder channel map into the run.
  bool const fWriteChannelMap;
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  icarusDB::IICARUSChannelMapProvider const& fPMTchannelMap;
  
  microseconds const fOpticalTick; ///< Duration of an optical tick. TODO is this still necessary?
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  /// Reshaping algorithm.
  std::unique_ptr<FilteringCircuit> const fReshapingFilter;
  
  /// Adder calibration algorithm.
  std::unique_ptr<AdderCalibrationDatabase> const fCalibrationDB;
  
  // this algorithm requires by-run setup so it's not declared `const`
  AdderSignalSimulation fAdderSimulation; ///< Adder simulation algorithm.
  
  // --- END Algorithms --------------------------------------------------------
  
  // --- BEGIN Caches ----------------------------------------------------------
  
  /// Calibration for the current run.
  std::unique_ptr<AdderCalibrationDatabase::RunCalibration> fRunCalibration;
  
  /// Adder channel mapping for the current run.
  AdderChannelMap fAdderChannels;
  
  // --- END Caches ------------------------------------------------------------
  
  /// Returns a set of ADC settings matching the configuration.
  adder::types::ADCsettings_t makeADCsettings() const;
  
  /**
   * @brief Combines waveforms and their baselines.
   * @param waveforms waveform data product (or must match data product order)
   * @param waveformBaselines waveform-to-baseline association
   * @return a collection of `WaveformWithBaseline` (same order as `waveforms`)
   * @throw art::Exception (code `art::errors::NotFound`) on missing baseline
   */
  std::vector<WaveformWithBaseline> fetchWaveformsAndBaselines(
    std::vector<raw::OpDetWaveform> const& waveforms,
    art::FindOne<icarus::WaveformBaseline> const& waveformBaselines
    ) const;
  
  
  /// Creates and returns a directory where to put all the plots for this event.
  TDirectoryFile* makeEventPlotDir(art::EventID const& ID) const;
  
  
  /// Sorts the `waveforms` in place, by channel and then by start time.
  static void sortByChannelAndTime(std::vector<raw::OpDetWaveform>& waveforms);


  /**
   * @brief Returns metadata and associations out of a list of `waveforms`.
   * @param event the event where the metadata is going to be stored
   * @param detTimings detector timing information
   * @param waveforms the waveforms to extract metadata out of
   * @return a pair: metadata and associations between metadata and waveforms
   * 
   * Metadata and associations are
   * @ref LArSoftProxyDefinitionParallelData "parallel" to the `waveforms`
   * collection.
   */
  std::pair<
    std::vector<sbn::OpDetWaveformMeta>,
    art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>
    >
  createMetadata(
    art::Event const& event, detinfo::DetectorTimings const& detTimings,
    std::vector<raw::OpDetWaveform> const& waveforms
    ) const;
  
  /**
   * @brief Creates a simple parallel association between new data products.
   * @tparam Left type of the data product on the left side of the association
   * @tparam Right type of the data product on the right side of the association
   * @param event the _art_ event where the data products are going to be stored
   * @param left the data product on the left side of the association
   * @param right the data product on the right side of the association
   * @return the _art_ association
   * 
   * @ref LArSoftProxyDefinitionParallelData "Parallel data products" are
   * collections of objects associated one-to-one in order (the first object of
   * one collection with the first object of the other collection, the second
   * with the second, etc.).
   * While this does not require the existence of an explicit association
   * object, this function creates it anyway.
   * 
   * There are strong assumptions here:
   *  * the two collections have the same size (asserted, not checked);
   *  * both collections are new data products not yet in the event;
   *  * both collections will become data products an empty instance name;
   *  * no metadata is needed.
   */
  template <typename Left, typename Right>
  static art::Assns<Left, Right> makeParallelAssociation(
    art::Event const& event,
    std::vector<Left> const& lefts, std::vector<Right> const& rights
    );
    
}; // icarus::trigger::SimulateAdderSignal


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Returns a copy of the collection `coll`, sorted by `std::sort()`.
  template <typename Coll>
  Coll sorted(Coll coll)
    {
      Coll sortedColl{ coll };
      std::sort(begin(sortedColl), end(sortedColl));
      return sortedColl;
    }
  
  /// Shoves the content of `data` into a `std::unique_ptr`.
  template <typename T>
  std::unique_ptr<std::remove_reference_t<T>> moveToUniquePtr(T&& data)
    { return std::make_unique<std::remove_reference_t<T>>(std::move(data)); }
  
  /// Creates an algorithm object (type `Algo`) only if `config` is present.
  template <typename Algo, typename Config>
  std::unique_ptr<Algo> configureAlgo(Config const& config)
    { return config? std::make_unique<Algo>(*config): nullptr; }

} // local namespace

namespace raw {
  
  /// Prints a short summary of the waveform.
  std::ostream& operator<< (std::ostream& out, raw::OpDetWaveform const& waveform)
  {
    raw::Channel_t const channel = waveform.ChannelNumber();
    out << "CH=" << icarus::trigger::AdderChannelID{ channel }
      << " T=" << waveform.TimeStamp()
      << " us (" << waveform.Waveform().size() << " samples)";
    return out;
  }
  
} // namespace raw


//------------------------------------------------------------------------------
//--- icarus::trigger::SimulateAdderSignal
//------------------------------------------------------------------------------
icarus::trigger::SimulateAdderSignal::SimulateAdderSignal
  (Parameters const& params, art::ProcessingFrame const&)
  // base classes
  : art::SharedProducer{ params }
  , icarus::ns::util::mfLoggingClass{ params().simulation().LogCategory() }
  // configuration parameters
  , fWaveformTag    { params().WaveformTag()     }
  , fBaselineTag    { params().BaselineTag()     }
  , fDebugPlots     { params().DebugPlots()      }
  , fWriteChannelMap{ params().WriteChannelMap() }
  // service cache
  , fPMTchannelMap
    { *art::ServiceHandle<icarusDB::IICARUSChannelMap>()->provider() }
  , fOpticalTick {
    detinfo::makeDetectorClocksWithUnits
      (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob())
      .OpticalClockPeriod()
    }
  // algorithms
  , fReshapingFilter{
      configureAlgo<SallenKeyFilter>(params().ReshapingCircuitParameters())
    }
  // TODO make this into a tool to freely choose which calibration to use
  , fCalibrationDB
    { configureAlgo<AmplitudeAdderCalibration>(params().Calibration()) }
  , fAdderSimulation { params().simulation(), makeADCsettings() }
{
  
  if (fDebugPlots)
    serialize(art::SharedResource<art::TFileService>);
  else
    async<art::InEvent>(); // FIXME is FFTW algorithm reentrant? [A: it is if we call it in a way we did not]
  
  consumes<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  consumes<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>
    (fBaselineTag);
  
  if (fWriteChannelMap) produces<ChannelMap_t, art::InRun>();
  produces<std::vector<raw::OpDetWaveform>>();
  produces<std::vector<sbn::OpDetWaveformMeta>>();
  produces<art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>>();
  produces<std::vector<icarus::WaveformBaseline>>();
  produces<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>();
  
  //
  // dump configuration
  //
  {
    auto log = mfLogInfo();
    log << "SimulateAdderSignal configuration:"
      << "\n * PMT waveforms from '" << fWaveformTag.encode() << "'"
      << "\n     * baselines from '" << fBaselineTag.encode() << "'"
      << "\n * " << fAdderSimulation.configurationDump("     ", "")
      ;
    
    if (fReshapingFilter) {
      log << "\n * filter configuration:"
        << "\n" << fReshapingFilter->printConfig("   ");
    }
    else log << "\n * no filter applied";
    if (fCalibrationDB) {
      log << "\n * reshaped waveform calibration ";
      fCalibrationDB->dumpConfig(log, "   ", "");
    }
    else log << "\n * no calibration applied";
    log
      << "\n * " << (fDebugPlots? "": "do not ") << "produce debugging plots";
    if (fWriteChannelMap)
      log << "\n * write the extracted adder channel map into the run data";
  }
  
} // icarus::trigger::SimulateAdderSignal::SimulateAdderSignal()


//------------------------------------------------------------------------------
auto icarus::trigger::SimulateAdderSignal::makeADCsettings() const
  -> adder::types::ADCsettings_t
{
  adder::types::ADCsettings_t settings; // with all default values;
  settings.samplingTime = fOpticalTick.quantity();
  return settings;
} // icarus::trigger::SimulateAdderSignal::makeADCsettings()


//------------------------------------------------------------------------------
void icarus::trigger::SimulateAdderSignal::beginRun
  (art::Run& run, art::ProcessingFrame const&)
{
  
  AdderChannelMapBuilder const channelMapBuilder{ logCategory() };
  fAdderChannels = channelMapBuilder.build(fPMTchannelMap, run.run());
  
  if (fCalibrationDB)
    fRunCalibration = fCalibrationDB->calibrationForRun(run.run());
  
  fAdderSimulation.setup
    (&fAdderChannels, fReshapingFilter.get(), fRunCalibration.get());
  
  if (fWriteChannelMap) {
    // build and write the mapping data product
    ChannelMap_t map;
    for (auto const& adderInfo: util::const_values(fAdderChannels)) {
      map.addChannel
        (adderInfo.channel.channel(), { begin(adderInfo.PMTs), end(adderInfo.PMTs) });
    } // for
    
    run.put(moveToUniquePtr(map), art::fullRun());
  } // if fWriteChannelMap
  
} // icarus::trigger::SimulateAdderSignal::beginRun()


//------------------------------------------------------------------------------
void icarus::trigger::SimulateAdderSignal::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // prepare input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  auto const& waveforms = *waveformHandle;
  
  // retrieve and integrate the baseline information
  art::FindOne<icarus::WaveformBaseline> const waveformBaselines
    { waveformHandle, event, fBaselineTag };
  std::vector<WaveformWithBaseline> const waveformInfo
    = fetchWaveformsAndBaselines(waveforms, waveformBaselines);
  
  mfLogDebug()
    << "Addition of " << waveforms.size() << " PMT waveforms from '"
    << fWaveformTag.encode() << "' with baselines from '"
    << fBaselineTag.encode() << "'";
  
  //
  // run simulation
  //
  detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings
    (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event))
    ;
  
  std::vector<raw::OpDetWaveform> adderWaveforms = fAdderSimulation.simulate
    (waveformInfo, detTimings, makeEventPlotDir(event.id()));
  
  sortByChannelAndTime(adderWaveforms);
  
  //
  // waveform baselines
  //
  
  // this is kind of embarrassing, since they are expected to be `0` by
  // construction and we abuse that assumption here;
  // `icarus::WaveformBaseline` default-constructs to 0
  std::vector<icarus::WaveformBaseline> adderWaveformBaselines
    (adderWaveforms.size());
  art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline> adderWaveformToBaseline
    = makeParallelAssociation(event, adderWaveforms, adderWaveformBaselines);
  
  //
  // get the timing information for this event
  //
  //
  // waveform metadata
  //
  {
    electronics_time const triggerTime = detTimings.TriggerTime();
    electronics_time const beamGateTime = detTimings.BeamGateTime();
    mfLogDebug()
      << "Event " << event.id() << " has beam gate starting at " << beamGateTime
      << " and trigger at " << triggerTime << ".";
  }
  
  // std::vector<sbn::OpDetWaveformMeta> AdderInfo
  // art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform> metaToWaveform
  auto [ adderInfo, metaToWaveform ]
    = createMetadata(event, detTimings, adderWaveforms);
  
  //
  // store the output into the event
  //
  event.put(moveToUniquePtr(adderWaveforms));
  event.put(moveToUniquePtr(adderInfo));
  event.put(moveToUniquePtr(metaToWaveform));
  event.put(moveToUniquePtr(adderWaveformBaselines));
  event.put(moveToUniquePtr(adderWaveformToBaseline));
  
} // icarus::trigger::SimulateAdderSignal::produce()


//------------------------------------------------------------------------------
auto icarus::trigger::SimulateAdderSignal::fetchWaveformsAndBaselines(
  std::vector<raw::OpDetWaveform> const& waveforms,
  art::FindOne<icarus::WaveformBaseline> const& waveformBaselines
) const
  -> std::vector<WaveformWithBaseline>
{

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
  } // for
  
  return waveformInfo;
} // icarus::trigger::SimulateAdderSignal::fetchWaveformsAndBaselines()


//------------------------------------------------------------------------------
TDirectoryFile* icarus::trigger::SimulateAdderSignal::makeEventPlotDir
  (art::EventID const& ID) const
{
  if (!fDebugPlots) return nullptr;
  
  using util::toStr;
  
  auto const paddedStr = [](int v, int padding, char padChar = '0'){
    std::string sv = std::to_string(v);
    return std::string(std::max(padding - (int) sv.length(), 0), padChar) + sv;
  };
  
  std::string const name
    = "R" + paddedStr(ID.run(), 5, '0') + "_E" + paddedStr(ID.event(), 7, '0');
  std::string const title
    = toStr << "Run " << ID.run() << " event " << ID.event();
  
  return art::ServiceHandle<art::TFileService>()->make<TDirectoryFile>
    (name.c_str(), title.c_str());

} // icarus::trigger::SimulateAdderSignal::makeEventPlotDir()


//------------------------------------------------------------------------------
namespace {
  template <typename Func>
  struct memberFuncLess {
    Func funcPtr;
    memberFuncLess(Func funcPtr): funcPtr{ funcPtr } {}
    template <typename T>
    bool operator() (T const& a, T const& b) const
      { return (a.*funcPtr)() < (b.*funcPtr)(); }
  }; // memberFuncLess
} // local namespace

void icarus::trigger::SimulateAdderSignal::sortByChannelAndTime
  (std::vector<raw::OpDetWaveform>& waveforms)
{
  std::stable_sort(begin(waveforms), end(waveforms),
    memberFuncLess{ &raw::OpDetWaveform::TimeStamp });
  std::stable_sort(begin(waveforms), end(waveforms),
    memberFuncLess{ &raw::OpDetWaveform::ChannelNumber });
}


//------------------------------------------------------------------------------
std::pair
  <std::vector<sbn::OpDetWaveformMeta>, art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>>
icarus::trigger::SimulateAdderSignal::createMetadata(
  art::Event const& event, detinfo::DetectorTimings const& detTimings,
  std::vector<raw::OpDetWaveform> const& waveforms
) const {
  
  sbn::OpDetWaveformMetaMaker const makeOpDetWaveformMeta { detTimings };
  
  std::vector<sbn::OpDetWaveformMeta> adderInfo;
  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
    assert(iWaveform == adderInfo.size());
    
    adderInfo.push_back(makeOpDetWaveformMeta(waveform));
    
    {
      sbn::OpDetWaveformMeta const& info = adderInfo.back();
      auto log = mfLogTrace();
      log << "Coverage for waveform #" << iWaveform
        << " on channel " << AdderChannelID{ info.channel } << ": "
        << info.startTime << " -- " << info.endTime;
      if (info.withTrigger()) log << "; includes trigger";
      if (info.withBeamGate()) log << "; includes beam gate start";
    }
    
  } // for
  
  art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform> infoToWaveform
    = makeParallelAssociation(event, adderInfo, waveforms);
  
  return { std::move(adderInfo), std::move(infoToWaveform) };
  
} // icarus::trigger::SimulateAdderSignal::createMetadata()


//------------------------------------------------------------------------------
template <typename Left, typename Right>
art::Assns<Left, Right>
icarus::trigger::SimulateAdderSignal::makeParallelAssociation(
  art::Event const& event,
  std::vector<Left> const& lefts, std::vector<Right> const& rights
) {
  
  std::vector<sbn::OpDetWaveformMeta> adderInfo;
  art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform> infoToWaveform;
  
  // assuming the products are new, and assuming they have empty instance name
  art::PtrMaker<Left> const makeLPtr { event };
  art::PtrMaker<Right> const makeRPtr { event };
  
  assert(size(lefts) == size(rights));
  art::Assns<Left, Right> assns;
  for (std::size_t const i: util::counter(size(lefts)))
    assns.addSingle(makeLPtr(i), makeRPtr(i));
  
  return assns;
  
} // icarus::trigger::SimulateAdderSignal::makeParallelAssociation()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::SimulateAdderSignal)


//------------------------------------------------------------------------------
