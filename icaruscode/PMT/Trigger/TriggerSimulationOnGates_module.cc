/**
 * @file   TriggerSimulationOnGates_module.cc
 * @brief  Plots of efficiency for triggers based on independent window gates.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 27, 2021
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.h"
#include "icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.h" // WindowTopologyManager
#include "icaruscode/PMT/Trigger/Algorithms/WindowPatternConfig.h"
#include "icaruscode/PMT/Trigger/Algorithms/WindowPattern.h"
#include "icaruscode/PMT/Trigger/Algorithms/ApplyBeamGate.h"
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h"
#include "icaruscode/PMT/Trigger/Algorithms/LVDSbitMaps.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // FillTriggerGates()
#include "icaruscode/PMT/Trigger/Utilities/DiscrThresholdParser.h"
#include "icaruscode/PMT/Algorithms/ADCsettings.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()...
#include "icaruscode/Utilities/CacheCounter.h" // util::CacheGuard
#include "icarusalg/Utilities/BeamBitChoiceSelectors.h" // sbn::bits::triggerLogic
#include "icarusalg/Utilities/CommonChoiceSelectors.h" // util::TimeScale...
#include "icarusalg/Utilities/PlotSandbox.h"
#include "icarusalg/Utilities/ROOTutils.h" // util::ROOT
#include "icarusalg/Utilities/BinningSpecs.h"
#include "icarusalg/Utilities/FixedBins.h"
#include "icarusalg/Utilities/PassCounter.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "icarusalg/Utilities/ChangeMonitor.h" // ThreadSafeChangeMonitor
#include "icarusalg/Utilities/rounding.h" // icarus::ns::util::roundup()
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnobj/Common/Trigger/BeamBits.h" // sbn::triggerSource, mask()...

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_tick...
#include "lardataalg/Utilities/quantities/electromagnetism.h" // millivolt
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds, ...
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/get_elements.h" // util::get_elements()
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger
#include "lardataobj/RawData/OpDetWaveform.h" // raw::ADC_Count_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::CryostatID

// framework libraries
#if 0
#  include "art_root_io/TFileService.h"
#endif // 0
#  include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TEfficiency.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

// C/C++ standard libraries
#include <ostream>
#include <algorithm> // std::fill(), std::any_of(), std::find()
#include <map>
#include <set>
#include <vector>
#include <iterator> // std::make_move_iterator()
#include <memory> // std::make_unique()
#include <string>
#include <atomic>
#include <optional>
#include <functional> // std::mem_fn()
#include <utility> // std::pair<>, std::move()
#include <cmath> // std::ceil()
#include <cstddef> // std::size_t
#include <cassert>


//------------------------------------------------------------------------------
using namespace util::quantities::time_literals;


//------------------------------------------------------------------------------
namespace icarus::trigger { class TriggerSimulationOnGates; }
/**
 * @brief Simulates a trigger at specified gate times.
 * 
 * This module produces `raw::Trigger` objects each representing the outcome of
 * some trigger logic applied to discriminated optical detector input
 * ("trigger primitives").
 * The logic is applied to each event at multiple times, according to a list
 * of time intervals read from each event. In addition, the module finds all the
 * possible triggers _within_ each gate (honouring a configurable dead time).
 * 
 * The main purpose of this module is to simulate the trigger logic at times of
 * special interest, typically the times some track is believed to have crossed
 * the detector.
 * 
 * This module applies the configured window pattern to the input: the pattern
 * consists of a requirement on the main window and optional additional
 * requirements on the neighbouring windows. This module rebases the configured
 * pattern on each of the available windows, evaluates the requirement of the
 * pattern in that configuration, and decides whether those requirements are
 * met. The general trigger is considered passed if _any_ of the rebased
 * patterns satisfies the requirement at any time, and no special treatment is
 * performed in case multiple windows fulfil them, except that the trigger time
 * is driven by the earliest of the satisfied patterns.
 * 
 * A single trigger pattern is configured for each instance of the module,
 * while multiple input sets (e.g. with different discrimination thresholds)
 * can be processed on the same pattern by the same module instance.
 * Conversely, testing a different pattern requires the instantiation of a new
 * module.
 * 
 * The input of this module is one discretized waveform _per window_. In the
 * typical ICARUS sliding window configuration, each window is a 6 metre area
 * on one side of the TPC (16 LVDS pairs, 30 PMT). In the typical adder trigger
 * configuration, each window is one adder signal (corresponding to 15 PMT).
 * The input discretized waveforms must preserve some information of their
 * physical locations, that is they need to remember the PMT channels they
 * ultimately originate from.
 * 
 * 
 * Configuration
 * ==============
 * 
 * Run `lar --print-description TriggerSimulationOnGates` for a configuration
 * example and terse explanations (it will also include the full list of
 * available options for multiple choice parameters like `BeamGateReference`).
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module instance which
 *     produced the trigger primitives to be used as input; it must not include
 *     any instance name, as the instance names will be automatically added from
 *     `Thresholds` parameter.
 *     The typical trigger primitives used as input are LVDS discriminated
 *     output combined into trigger windows (e.g. from
 *     `icarus::trigger::SlidingWindowTrigger` module).
 * * threshold specifications: `Thresholds` (list of names) and `VThresholds`
 *     (list of voltages): list of the discrimination thresholds to include.
 *     A data product containing a digital signal is read for each of the
 *     thresholds in `Thresholds`, and the tag of the data product is
 *     expected to be the instance name in this configuration parameter for the
 *     module label set in `TriggerGatesTag` (e.g. for a threshold of
 *     `"400"`, supposedly 400 ADC counts, and with `TriggerGatesTag` set to
 *     `"TrigSlidingWindows"`, the data product tag would be
 *     `TrigSlidingWindows:60`). The thresholds in `VThresholds` are specified
 *     as voltage quantities (e.g. `"60 mV"`), and they are converted to ADC
 *     counts (integral number) according to the CAEN V1730B board settings
 *     (`icarus::ADCsettings`). The resulting number (e.g. `491` for `"60 mV"`)
 *     is used as if it had been specified in `Thresholds` instead, that is
 *     directly as tag for input data products.
 *     At least one threshold must be specified, either with `Thresholds` or
 *     with `VThresholds`. Duplicate entries will cause a configuration error.
 * * `LVDSgatesTag` (string, default: empty): if non-empty, input tags will be
 *     created in the same way as for `TriggerGatesTag`, and the resulting
 *     data products will be used to fill the
 *     `sbn::ExtraTriggerInfo::CryostatInfo::LVDSstatus` information
 *     (see @ref TriggerSimulationOnGates_Output "Output data products").
 * * `KeepThresholdName` (flag, optional): by default, output data products have
 *     each an instance name according to their threshold (from the `Threshold`
 *     parameter), unless there is only one threshold specified. If this
 *     parameter is specified as `true`, the output data product always
 *     includes the threshold instance name, even when there is only one
 *     threshold specified. If this parameter is specified as `false`, if there
 *     is only one threshold the default behaviour (of not adding an instance
 *     name) is confirmed; otherwise, it is a configuration error to have this
 *     parameter set to `false`.
 * * `Pattern` (configuration table, mandatory): describes the sliding window
 *     pattern; the configuration format for a pattern is described under
 *     `icarus::trigger::ns::fhicl::WindowPatternConfig`.
 * * `BeamGates` (input tag, _mandatory_): the data product with the beam gates
 *     to run the simulation in;
 * * `BeamBits` (bitmask as 32-bit integral number, optional): bits to be set in
 *     the produced `raw::Trigger` objects. If omitted, bits are created from
 *     a `sbn::bits::triggerSource` mask reflecting the content of the input
 *     beam gate bits, which basically can distinguish only between BNB, NuMI
 *     and others, "unknown" (see also `daq::TriggerDecoder` tool).
 * * `BeamGateReference` (text, default: `BeamGate`): time reference of the
 *     input gates in the data product defined by `BeamGates`. By LArSoft
 *     standard, that is simulation time in Monte Carlo, which is mirrored in
 *     data by the hardware beam gate time (although the latter allows a margin
 *     around the beam spill opening a bit before neutrinos are expected).
 * * `TriggerOnTransition` (flag, default: `False`): if set, a trigger candidate
 *     is issued only when there is a transition from non-satisfying the trigger
 *     requirements to satisfying them. This implicitly excludes all candidates
 *     coming earlier than the beam gate opening. The default setting of `false`
 *     reflects ICARUS trigger behaviour until Run2 (included), when only the
 *     current status of the trigger was considered. The setting of `true`
 *     reflects the trigger from run `11649` on (and from Run3 on), which
 *     effectively removes the peak of triggers at the time of the opening of
 *     the beam gate.
 * * `DeadTime` (time, default: forever): when looking for sequences of
 *     triggers, ignore this much time after every trigger found. This applies
 *     only within each requested gate: a trigger at the very end of a gate
 *     will not suppress one at the beginning of the next one even if that is
 *     within the set dead time. If set to a very large value (as the default),
 *     only one trigger will be found per input gate.
 * * `TriggerDelay` (time, default: `0 ns`): fixed time to add to the time of
 *     all the triggers.
 * * `EmitEmpty` (flag, default: `true`): if set, each gate gets at least a
 *     trigger object, and if there is no trigger during a gate its trigger
 *     object will be marked by having no bit set. If unset, when there is no
 *     trigger in the gate, no trigger object will be produced, but the counts
 *     will still proceed.
 * * `ExtraInfo` (flag, default: `false`): also produces a data product
 *     `sbn::ExtraTriggerInfo` with reduced information from the _first_ of the
 *     triggers from the _first_ of the gates. If the first gate did not trigger
 *     the object will have fields marked invalid.
 * * `LVDSstatusDelay` (nanoseconds, default: `0`): when extraction of LVDS bits
 *     for `LVDSstatus` is requested (see
 *     @ref TriggerSimulationOnGates_Output "Output data products" section),
 *     save in the bits the state of the gates from `LVDSgatesTag` at the time
 *     the trigger conditions are met plus the delay specified in this
 *     parameter.
 * * `RetriggeringBit` (positive integer, default: `17`): the bit to set for all
 *     the triggers found after the first one within each gate. The value `0`
 *     represents the least significant bit; the default value is `17`, bit mask
 *     `0x20000`. Using a value larger than the size of the trigger bit field
 *     will disable this mark.
 * * `TriggerLogicBit` (bit name from `sbn::bits::triggerLogic`, default:
 *     disabled): the name of the bit to set in
 *     `sbn::ExtraTriggerInfo::CryostatInfo::triggerLogicBits` when a trigger is
 *     found. It describes the type of input to the trigger logic and the type
 *     of its combination. If left empty, no bit is assigned.
 * * `TriggerClock` (time, default: 25 ns): the period of the clock applying
 *     the trigger logic. This value is used as unit for some of logic (e.g.,
 *     `sbn::ExtraTriggerInfo::CryostatInfo::beamToTrigger` and the meaning of
 *     `OverlapTicks`). This parameter replaces `TriggerTimeResolution`.
 * * `OverlapTicks` (integer, default: 0): if the difference between the tick
 *     of two triggers from different sources is this many trigger clock ticks
 *     or less, the triggers are merged (the first trigger determines almost
 *     everything, but `sbn::ExtraTriggerInfo::CryostatInfo` will have
 *     information from both and the source bit mask will be set for both.
 * * `CryostatFirstBit` (positive integer, default: disabled): the bit to set
 *     for triggers from windows in cryostat 0; as many bits will be used as
 *     there are cryostats in the detector (thus, 2 for ICARUS). The value `0`
 *     represents the least significant bit. Using a value larger than the size
 *     of the trigger bit field (which is the default) will disable this mark.
 * * `LogCategory` (string, default `TriggerSimulationOnGates`): name of
 *     category used to stream messages from this module into message facility.
 * 
 * An example job configuration is provided as
 * `simulate_sliding_window_trigger_icarus.fcl`.
 * 
 * 
 * Dependencies
 * =============
 * 
 * The following services are required:
 *  * The usual message facility.
 *  * `geo::Geometry` and `geo::WireReadout` to build trigger window topology.
 *  * `detinfo::DetectorClocksService` for timing translations.
 *  * `icarusDB::IICARUSChannelMap` if output of LVDS status is requested.
 *  * `art::TFileService` in principle for plots, but currently there is none
 *     (the service is required nonetheless).
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `TriggerGatesTag` + `Thresholds`: input discretized waveform ("gate")
 *    collections. Each gate, of type `icraus::trigger::OpticalTriggerGateData_t`,
 *    needs to remember or be associated to physical PMT.
 * * `BeamGates` (`std::vector<sim::BeamGateInfo>`): the beam gate intervals
 *     to run the simulation on; one trigger result is produced and saved for
 *     each of the gates in this data product. By default the gates are
 *     interpreted following LArSoft convention for the simulation, with the
 *     times in nanoseconds and in
 *     @ref DetectorClocksSimulationTime "simulation time reference".
 *     The parameter `BeamGateReference` can change that interpretation.
 * * `LVDSgatesTag` + `Thresholds`: LVDS input gate collections (if LVDS status
 *     output is requested: see
 *     @ref TriggerSimulationOnGates_Output "Output data products" section).
 * 
 * 
 * Output data products
 * =====================
 * 
 * @anchor TriggerSimulationOnGates_Output
 *
 * For each ADC threshold there will be an instance of the following:
 * 
 * * `std::vector<raw::Trigger>`:
 *     list of triggers fired according to the configured trigger definition;
 *     there is one collection (and data product) per ADC threshold, and the
 *     data product has the same instance name as the input data one, unless
 *     there is only one threshold (see `TriggerGatesTag`, `Thresholds` and
 *     `KeepThresholdName` configuration parameters);
 *     at least one trigger object is produced for each of the beam gates found
 *     in the input data product specified by the `BeamGates` parameter,
 *     unless `EmitEmpty` is set, in which case beam gates with no trigger fired
 *     will not contribute any object.
 *     Each trigger object has the time stamp matching the time when the trigger
 *     criteria were satisfied, plus the fixed delay in the `TriggerDelay`
 *     configuration parameter. All triggers feature the bits specified in
 *     `BeamBits` configuration parameter, with the following exceptions:
 *     if there was no trigger found, the bits will all be cleared; and if there
 *     was more than one trigger found in the same gate, all triggers except the
 *     first one (which all share the same ID) will have the `RetriggeringBit`
 *     set.
 * * `sbn::ExtraTriggerInfo` (if `ExtraInfo` configuration parameter is set):
 *     always present, and reflecting the triggers on the first beam gate: if a
 *     trigger fired in that beam gate, a _reduced_ version of
 *     `sbn::ExtraTriggerInfo` is provided; currently it is guaranteed to have:
 *     * `triggerTimestamp`: based on the _art_ event timestamp, or invalid
 *       timestamp if the trigger did not fire.
 *     * `beamGateTimestamp`: set accordingly to the relative time of the
 *       trigger vs. simulated beam gate opening, if the trigger happened.
 *       Otherwise, it will be set at the event time.
 *     * `sourceType`: interpreted according to the trigger bits of the input
 *       beam gate object, it may be `sbn::triggerSource::BNB`,
 *       `sbn::bits::triggerSource::NuMI` or `sbn::bits::triggerSource::Unknown`
 *       (note that in any case it's not left to its default value
 *       `sbn::bits::triggerSource::NBits`).
 *     * `triggerType`: majority type is always set.
 *     * `triggerLocationBits`: set to either cryostat east or west depending on
 *        where trigger was (bit mask for `sbn::bits::triggerLocation::CryoEast`
 *        or `sbn::bits::triggerLocation::CryoWest`).
 *     * `triggerID` and `triggerCount` match the trigger number in
 *       `raw::Trigger::TriggerNumber()`, except in case no trigger fired,
 *       in which case `triggerID` is `sbn::ExtraTriggerInfo::NoID` and
 *       `triggerCount` is `0`, both in their default values.
 *     * `gateID` and `gateCount` match the event number.
 *     * `cryostats`: information per cryostat:
 *         * `beamToTrigger`: time from beam gate opening to the time trigger
 *           conditions are met. This interval does not include the
 *           `TriggerDelay`. The resolution of this time only reflects the
 *           digitized input resolution, without any quantization (e.g., from
 *           hardware clocks).
 *         * `triggerLogicBits`: value from `TriggerLogicBit` configuration
 *           parameter (`0` if no trigger happened).
 *         * `triggerCount`: number of triggers found in the gate. This is the
 *           actual number of triggers found, and there is no attempt to
 *           emulate the equivalent count in the hardware.
 *         * `LVDSstatus` (if `LVDSgatesTag` is set): state of each channel pair
 *           at trigger time. A data product tag derived from `LVDSgatesTag`
 *           configuration parameter is used to read the state of all PMT pairs.
 *           The state of each pair is evaluated at the tick of the emulated
 *           trigger time (without the `TriggerDelay`) plus a specific freezing
 *           delay (from `LVDSstatusDelay` configuration parameter).
 *           This state is assigned to a "logic LVDS" bit according to the LVDS
 *           bit mapping read from `IICARUSChannelMap` service
 *           (see `icarus::trigger::LVDSbitMaps`).
 *   
 *     If the first gate did not emit a trigger, the object will be left
 *     default-constructed, noticeably with an invalid trigger timestamp
 *     (`sbn::ExtraTriggerInfo::NoTimestamp`). Note that there may be exceptions
 *     to this rule (for example, `triggerSource`, `beamGateTimestamp` and the
 *     gate counts described above). Also note that
 *     `sbn::ExtraTriggerInfo::isValid()`, bound to the source of the trigger,
 *     will return `true` even when there was no trigger firing.
 * 
 * 
 * ### Trigger bits in `raw::Trigger`
 * 
 * The trigger bits may be overridden by the `BeamBits` configuration parameter.
 * If they are not, the bits from an input "beam" gate will be used verbatim for
 * all the triggers within that gate.
 * 
 * On top of this base set of bits, this module may add a few additional bits.
 * 
 * The module can be configured so that multiple trigger objects are created for
 * a single input gate. In that case, the triggers after the first one may be
 * marked as "retriggering" by setting a bit whose position is determined by
 * `RetriggeringBit` configuration parameter.
 * 
 * In addition, the module can mark a trigger according to which cryostat it
 * originated from. The position of the bit for cryostat `0` is set by the
 * `CryostatFirstBit`, and the next cryostat(s) will follow in their order with
 * more significant bits.
 * 
 * For the options where a bit position can be specified, specifying a position
 * beyond the number of available bits effectively disables that option.
 * 
 * 
 * ### Trigger timing
 * 
 * This module replicates the logic used by the hardware to provide a trigger
 * response, but it does not attempt to emulate the delays that occur in the
 * hardware. The trigger time is by default the instant when in the input
 * waveforms the conditions are met, assuming that their synchronization already
 * reflects the one on the hardware input.
 * The only allowance to the hardware timing is the option to add a _fixed_
 * delay (`TriggerDelay` configuration parameter) to all trigger times produced
 * by this module. The `beamToTrigger` delay, on the other end, is deliberately
 * excluded from that delay, and no specific delay is offered for it.
 * 
 * 
 * 
 * 
 * Trigger logic algorithm
 * ========================
 * 
 * @anchor TriggerSimulationOnGates_Algorithm
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::TriggerSimulationOnGates` and its assumptions.
 * Nevertheless, more up-to-date information can be found in
 * `SlidingWindowTrigger` module (for the combination of the LVDS signals into
 * window-wide gates) and in `icarus::trigger::SlidingWindowPatternAlg`,
 * which applies the configured pattern logic to the input.
 * 
 * The module receives as input a multi-level trigger gate for each of the
 * windows to be considered.
 * On the first input (i.e. the first event), that input is parsed to learn
 * the windows and their relative position from the input trigger gates.
 * This topology will be used to apply the configured patterns. On the following
 * events, their input is checked to confirm the compatibility of the
 * composition of its windows with the one from that first event (both aspects
 * are handled by an `icarus::trigger::WindowTopologyManager` object).
 * 
 * All multi-level gates are set in coincidence with the beam gate by
 * multiplying the multi-level and the beam gates. Because of this, trigger
 * gates are suppressed everywhere except than during the beam gate (see below).
 * The reference time for the beam gates is the time configured in
 * `DetectorClocks` service provider
 * (`detinfo::DetectorClocks::BeamGateTime()`).
 * 
 * The algorithm handles independently multiple trigger patterns.
 * On each input, each configured pattern is applied based on the window
 * topology. Each pattern describes a minimum level of the trigger
 * gate in the window, that usually means the number of LVDS signals in
 * coincidence at any given time ("majority"). A pattern may have requirements
 * on the neighbouring windows in addition to the main one. The pattern is
 * satisfied if all involved windows pass their specific requirements at the
 * same time (coincidence between windows).
 * Each pattern is applied in turn to each of the windows (which is the "main"
 * window). The neighborhood described in the pattern is applied with respect to
 * that main window. The trigger fires if one or more of the windows satisfy
 * the pattern, and the trigger time is the one of the earliest satisfied
 * pattern (more precisely, the earliest tick when the coincidence required
 * by that pattern is satisfied).* 
 * All windows in the detector are considered independently, but the supported
 * patterns may only include components in the same cryostat. Therefore,
 * triggers are effectively on a single cryostat.
 * An object of class `icarus::trigger::SlidingWindowPatternAlg` applies this
 * logic: see its documentation for the most up-to-date details.
 * 
 * Eventually, for each event there are as many different trigger responses as
 * how many different patterns are configured (`Patterns` configuration
 * parameter), _times_ how many ADC thresholds are provided in input,
 * configured in `Thresholds`.
 * 
 * 
 * Beam gates
 * -----------
 * 
 * A single instance of this module can perform the simulation on several beam
 * gates. The values of these beam gates are picked from the data product
 * specified in `BeamGates`, event by event. The specified beam gate times are
 * on beam gate time scale, i.e. their reference time `0` is the time of the
 * beam gate as known by `detinfo::DetectorClocks::BeamGateTime()`.
 * In case the same beam gate is desired for all events, such data product can
 * be produced by `icarus::trigger::WriteBeamGateInfo` module.
 * The trigger data product collection produced by this module has the same
 * number of entries as the beam gates in the data product, and they match
 * one-to-one.
 * 
 * 
 * 
 * Technical aspects of the module
 * --------------------------------
 * 
 * @anchor TriggerSimulationOnGates_Tech
 * 
 * This module does not build the trigger gates of the sliding windows, but
 * rather it takes them as input (see e.g. `SlidingWindowTrigger` module).
 * Window topology (size of the windows and their relations) is stored in
 * `icarus::trigger::WindowChannelMap` objects, and its construction is
 * delegated to `icarus::trigger::WindowTopologyAlg` (under the hood of the
 * `WindowTopologyManager` class) which learns it from the actual trigger gate
 * input rather than on explicit configuration. Pattern definitions and
 * configuration are defined in `icarus::trigger::WindowPattern` and
 * `icarus::trigger::ns::fhicl::WindowPatternConfig` respectively. Trigger
 * simulation is delegated to `icarus::trigger::SlidingWindowPatternAlg`.
 * 
 * 
 * @todo Plots need to be thought and implemented.
 * 
 */
class icarus::trigger::TriggerSimulationOnGates
  : public art::EDProducer
  , private lar::UncopiableAndUnmovableClass
{
  
  /// Type of trigger bits supported by `raw::Trigger`.
  using TriggerBits_t = decltype(std::declval<raw::Trigger>().TriggerBits());
  
  /// Number of bits supported by `raw::Trigger`.
  static constexpr std::size_t NTriggerBits = sizeof(TriggerBits_t) * 8U;
  
    public:
  
  using microseconds = util::quantities::intervals::microseconds;
  using nanoseconds = util::quantities::intervals::nanoseconds;
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Sequence<std::string> Thresholds {
      Name("Thresholds"),
      Comment("tags of the thresholds to consider"),
      std::vector<std::string>{}
      };

    fhicl::Sequence<util::quantities::millivolt> VThresholds {
      Name("VThresholds"),
      Comment("thresholds to convert to ADC tags (exclusive with `Thresholds`)"),
      std::vector<util::quantities::millivolt>{}
      };

    fhicl::Atom<std::string> LVDSgatesTag {
      Name("LVDSgatesTag"),
      Comment(
        "data product to read LVDS state from (same rules as `TriggerGatesTag`)"
        ),
      "" // default
      };

    fhicl::OptionalAtom<bool> KeepThresholdName {
      Name("KeepThresholdName"),
      Comment
        ("add threshold to output product tag even with only one threshold")
      };

    icarus::trigger::ns::fhicl::WindowPatternTable Pattern {
      Name("Pattern"),
      Comment("trigger requirements as a trigger window pattern")
      };
 
    fhicl::Atom<art::InputTag> BeamGates {
      Name("BeamGates"),
      Comment("data product with all beam gates to run simulation into")
      };
    
    fhicl::OptionalAtom<std::uint32_t> BeamBits {
      Name("BeamBits"),
      Comment("bits to be set in the trigger object as beam identifier")
      };

    fhicl::Atom<util::TimeScale> BeamGateReference {
      Name{ "BeamGateReference" },
      Comment{ "time scale the beam gates refer to" },
      util::TimeScale::BeamGate
      };
    
    fhicl::Atom<bool> EmitEmpty {
      Name("EmitEmpty"),
      Comment("produce a trigger object even for gates with no firing trigger"),
      true
      };

    fhicl::Atom<bool> ExtraInfo {
      Name("ExtraInfo"),
      Comment("produce a sbn::ExtraTriggerInfo object out of the first gate"),
      false
      };

    fhicl::Atom<bool> TriggerOnTransition {
      Name("TriggerOnTransition"),
      Comment("only emit triggers on the first tick the requirements are met"),
      false
      };

    fhicl::Atom<nanoseconds> DeadTime {
      Name("DeadTime"),
      Comment("veto time after each trigger found"),
      std::numeric_limits<nanoseconds>::max()
      };
    
    fhicl::Atom<nanoseconds> TriggerDelay {
      Name("TriggerDelay"),
      Comment("trigger response delay, added to the trigger times"),
      0_ns // default
      };
    
    fhicl::Atom<nanoseconds> LVDSstatusDelay {
      Name("LVDSstatusDelay"),
      Comment("delay to freeze LVDS state bits, added to the trigger times"),
      0_ns // default
      };
    
    fhicl::Atom<unsigned int> RetriggeringBit {
      Name("RetriggeringBit"),
      Comment(
        "number of the bit to set for triggers from retriggering (< "
        + std::to_string(NTriggerBits) + ")"),
      17U
      };
    
    fhicl::Atom<sbn::bits::triggerLogic> TriggerLogicBit {
      Name("TriggerLogicBit"),
      Comment("trigger logic bit to set for these triggers"),
      sbn::bits::triggerLogic::NBits
      };
    
    fhicl::Atom<unsigned int> CryostatFirstBit {
      Name("CryostatFirstBit"),
      Comment(
        "number of the bit to set for triggers from the first cryostat (< "
        + std::to_string(NTriggerBits) + ")"),
      NTriggerBits
      };
    
    fhicl::OptionalAtom<nanoseconds> TriggerTimeResolution {
      Name("TriggerTimeResolution"),
      Comment("deprecated: use `TriggerClock`")
      };
    
    fhicl::OptionalAtom<nanoseconds> TriggerClock {
      Name("TriggerClock"),
      Comment("tick duration of the trigger logic hardware clock [25 ns]")
      };
    
    fhicl::Atom<unsigned int> OverlapTicks {
      Name("OverlapTicks"),
      Comment
        ("triggers from different cryostat closer than this ticks are merged"),
      0
      };
    
    fhicl::Atom<double> EventTimeBinning {
      Name("EventTimeBinning"),
      Comment("binning for the trigger time plot [second]"),
      300 // 5 minutes
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "TriggerSimulationOnGates" // default
      };
    
  }; // struct Config

  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit TriggerSimulationOnGates(Parameters const& config);

  // --- END Constructors ------------------------------------------------------


  // --- BEGIN Framework hooks -------------------------------------------------

  /// Initializes the plots.
  virtual void beginJob() override;
  
  /// Updates time-dependent caches.
  virtual void beginRun(art::Run&) override;
    
  /// Runs the simulation and saves the results into the _art_ event.
  virtual void produce(art::Event& event) override;
  
  /// Prints end-of-job summaries.
  virtual void endJob() override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // for convenience:
  using OpticalTriggerGateData_t = icarus::trigger::OpticalTriggerGateData_t;
  using optical_tick = detinfo::timescales::optical_tick;
  
  /// Information about the input for a single threshold.
  struct InputInfo_t {
    art::InputTag triggerGatesTag; ///< Gates for trigger response evaluation.
    art::InputTag LVDSgatesTag; ///< Original gates for LVDS state saving.
  };
  
  using TriggerInfo_t = details::TriggerInfo_t; ///< Type alias.
  
  /// Data structure to communicate internally a trigger response.
  using WindowTriggerInfo_t
    = icarus::trigger::SlidingWindowPatternAlg::AllTriggerInfo_t;
  
  /// Event-level information.
  struct EventAux_t {
    std::uint64_t time; ///< Event timestamp [ns]
    unsigned int event; ///< Event number.
  };

  /// Content for future histograms, binned.
  using BinnedContent_t = icarus::ns::util::FixedBins<double>;
  
  /// All information needed to generate plots for a specific threshold.
  struct ThresholdPlotInfo_t {
    BinnedContent_t eventTimes;
    BinnedContent_t HWtrigTimeVsBeam;
    BinnedContent_t triggerTimesVsHWtrig;
    BinnedContent_t triggerTimesVsBeam;
  };
  
  /// Type of list of gates to simulate trigger into.
  using BeamGates_t = std::vector<sim::BeamGateInfo>;
  
  
  /// Type for all PMT pair bits, per cryostat and per PMT wall.
  using LVDSbitArrays_t = std::array<
    std::array<std::uint64_t, sbn::ExtraTriggerInfo::MaxWalls>,
    sbn::ExtraTriggerInfo::MaxCryostats
    >;
  
  
  /// Utility to carry beam bits along with the beam gates.
  struct ApplyBeamGateClassWithBits: icarus::trigger::ApplyBeamGateClass {
    sbn::triggerSourceMask source; ///< Where trigger will look to come from.
  };
  
  
  using PlotSandbox_t = icarus::ns::util::PlotSandbox<art::TFileDirectory>;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Name of ADC thresholds to read, and the input connected to their data.
  std::map<std::string, InputInfo_t> fInputInfo;
  
  /// Configured sliding window requirement pattern.
  WindowPattern const fPattern;
  
  art::InputTag const fBeamGateTag; ///< Data product of beam gates to simulate.
  
  util::TimeScale const fBeamGateReference; ///< Reference time of beam gates.
  
  /// Bits for the beam gate being simulated.
  std::optional<std::uint32_t> fBeamBits;
  
  /// Whether to produce trigger objects for gates without triggers.
  bool const fEmitEmpty;
  
  bool const fExtraInfo; ///< Whether to produce a `sbn::ExtraTriggerInfo`.
  
  /// Whether to trigger only the moment requirements are met.
  bool const fTriggerOnTransition;
  
  nanoseconds const fDeadTime; ///< Veto time after a trigger in a gate.
  
  nanoseconds const fTriggerDelay; ///< Time to be added to the trigger time.
  
  nanoseconds const fLVDSstatusDelay; ///< Delay to freeze LVDS state bits.
  
  /// Bit mask set for triggers after the first one in a gate.
  TriggerBits_t const fRetriggeringMask;
  
  /// Bit describing the logic of this trigger (for `sbn::ExtraTriggerInfo`).
  sbn::bits::triggerLogicMask const fTriggerLogicMask;
  
  TriggerBits_t const fCryostatZeroMask; ///< Bit mask for cryostat `0`.
  
  nanoseconds const fTriggerClock; ///< Trigger resolution in time.
  
  unsigned int const fOverlapTicks; ///< Merge triggers closer than this.
  
  double const fEventTimeBinning; ///< Trigger time plot binning [s]
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------

  /// PMT channel mapping service provider.
  icarusDB::IICARUSChannelMapProvider const* fChannelMap = nullptr;
  
  ///< Tracks the cache of `IICARUSChannelMapProvider`.
  util::CacheGuard fChannelMapCacheGuard;
  
  // --- END Service variables -------------------------------------------------

  
  // --- BEGIN Internal variables ----------------------------------------------
  
  bool const fSaveLVDSbits; ///< Whether we need to save LVDS state.
  
  /// Output data product instance names (same order as `fInputInfo`).
  std::vector<std::string> fOutputInstances;
  
  /// Mapping of each sliding window with location and topological information.
  // mutable = not thread-safe
  mutable icarus::trigger::WindowTopologyManager fWindowMapMan;
  
  /// Pattern algorithm.
  std::optional<icarus::trigger::SlidingWindowPatternAlg> fPatternAlg;
  
  /// Cached LVDS bit mappings (`optional` because needs reinitialization).
  std::optional<icarus::trigger::LVDSbitMaps> fLVDSmaps;
  
#if 0
  /// All plots in one practical sandbox.
  PlotSandbox_t fPlots;
#endif // 0
  
  /// Proto-histogram information in a convenient packet; event-wide.
  ThresholdPlotInfo_t fEventPlotInfo;
  
  /// Proto-histogram information in a not-so-practical array; per threshold.
  std::vector<ThresholdPlotInfo_t> fThresholdPlots;

  ///< Count of fired triggers, per threshold.
  std::vector<std::atomic<unsigned int>> fTriggerCount;
  std::atomic<unsigned int> fTotalGates { 0U }; ///< Count of opened gates.
  
  
  // --- END Internal variables ------------------------------------------------
  

  // --- BEGIN --- Plot infrastructure -----------------------------------------
  
  /// @brief Initializes the full set of plots (all ADC thresholds).
  void initializePlots();
  
  /// Creates summary plots from proto-histogram data.
  void finalizePlots();
  
  /// Creates in `plots` sandbox all plots for threshold `threshold` from data
  /// in `plotInfo`.
  void makeThresholdPlots(
    std::string const& threshold,
    PlotSandbox_t& plots,
    ThresholdPlotInfo_t const& plotInfo
    );
  
  /// Creates in the main sandbox all event-wide plots.
  void makeEventPlots();
  
  /// Fills event-wide plots.
  void plotEvent(
    art::Event const& event, detinfo::DetectorTimings const& detTimings,
    std::vector<sim::BeamGateInfo> const& gates
    );
  
  /// Fills the plots for threshold index `iThr` with trigger information.
  void plotTriggerResponse(
    std::size_t iThr, std::string const& thrTag,
    std::vector<WindowTriggerInfo_t> const& triggerInfos,
    detinfo::DetectorTimings const& detTimings
    );
  
  // --- END ----- Plot infrastructure -----------------------------------------
  
  /// Creates the window map manager.
  icarus::trigger::WindowTopologyManager makeWindowMapManager() const;
  
  
  /**
   * @brief Performs the simulation for the specified ADC threshold.
   * @param event _art_ event to read data from and put results into
   * @param detTimings detector clocks service provider proxy
   * @param beamGates list of all beam gates to evaluate
   * @param iThr index of the threshold in the configuration
   * @param thr value of the threshold (ADC counts)
   * @param firstTriggerNumber the next unassigned trigger number
   * @return the trigger response information
   * 
   * For the given threshold, the simulation of the configured trigger is
   * performed.
   * The input data is read from the event (the source tag is from the module
   * configuration), simulation is performed, auxiliary plots are drawn and
   * a `raw::Trigger` collection is stored into the event.
   * 
   * The return value contains one trigger list per gate. Depending on the
   * configuration (and especially `fDeadTime`), one gate may be any number of
   * triggers, including none.
   * 
   * The simulation itself is performed by the `fPatternAlg` algorithm.
   */
  std::vector<std::vector<WindowTriggerInfo_t>> produceForThreshold(
    art::Event& event,
    detinfo::DetectorTimings const& detTimings,
    std::vector<sim::BeamGateInfo> const& beamGates,
    std::size_t const iThr, std::string const& thrTag,
    unsigned int firstTriggerNumber
    );
  
  /**
   * @brief Returns the first trigger in the specified gate.
   * @param start time start for the gate, in optical ticks
   * @param stop time end for the gate, in optical ticks
   * @param gates the trigger gates from all the windows needed in the pattern
   * @param detTimings detector timings helper
   * @return information about the first time trigger requirements were met
   */
  WindowTriggerInfo_t findFirstTrigger(
    optical_tick start, optical_tick stop,
    icarus::trigger::SlidingWindowPatternAlg::TriggerGates_t const& gates,
    detinfo::DetectorTimings const& detTimings
    ) const;
  
  /**
   * @brief Converts the trigger information into trigger objects.
   * @param detTimings detector clocks service provider proxy
   * @param beamGates list of all beam gates to evaluate
   * @param eventInfo event-wide information to be transferred to trigger data
   * @param triggerNumber the "unique" number to assign to this trigger
   * @param info the information about the fired triggers
   * @return objects with all the trigger information encoded
   * 
   * The triggers described by `info` are encoded into `raw::Trigger` objects.
   * If no trigger is present, one object may be made up depending of the
   * value configured in `fEmitEmpty`.
   * A `sbn::ExtraTriggerInfo` object is also returned representing the first
   * trigger in the beam gate (see the
   * @ref TriggerSimulationOnGates_Output "Output section of module documentation"
   * for details).
   */
  std::pair<std::vector<raw::Trigger>, sbn::ExtraTriggerInfo>
  triggerInfoToTriggerData(
    detinfo::DetectorTimings const& detTimings,
    sim::BeamGateInfo const& beamGate,
    EventAux_t const& eventInfo,
    unsigned int triggerNumber, std::vector<WindowTriggerInfo_t> const& info,
    std::vector<OpticalTriggerGateData_t> const* PMTpairGates
    ) const;
  
  /// Converts trigger bits from `beamInfo` into a `sbn::triggerSourceMask`.
  static sbn::triggerSourceMask makeTriggerBits
    (sim::BeamGateInfo const& beamInfo);
  
  /// Converts a `sim::BeamType_t` value into a `sbn::bits::triggerSource` one.
  static sbn::triggerSource beamTypeToTriggerSource(sim::BeamType_t beamType);
  
  /// Converts a cryostat ID into a bitmask pointing to that cryostat.
  static sbn::bits::triggerLocationMask cryoIDtoTriggerLocation
    (geo::CryostatID const& cid);
  
  /// Builds the PMT pair words from the status of the `gates` at `triggerTick`.
  LVDSbitArrays_t extractLVDSstatus(
    optical_tick triggerTick,
    std::vector<OpticalTriggerGateData_t> const& PMTgates
    ) const;
  
  /// Returns whether the `gate` channel list completely matches `channels`.
  static bool matchChannelList(
    OpticalTriggerGateData_t const& gate,
    std::vector<raw::Channel_t> const& channels
    );
  
  /// Returns a map from channel number to pointer to the gate covering it.
  static std::vector<OpticalTriggerGateData_t const*> makeGateMap
    (std::vector<OpticalTriggerGateData_t> const& gates);
  
  
  /// Prints the summary of fired triggers on screen.
  void printSummary() const;
  
  
  /**
   * @brief Converts a time into electronics time scale.
   * @param time the time to be converted
   * @param detTimings detector timings helper
   * @return `time` in electronics time scale
   * 
   * The time value is assumed to be in the reference as `fBeamGateReference`
   * and converted into the electronics time scale.
   */
  detinfo::timescales::electronics_time toElectronicsTime(
    util::quantities::nanosecond time,
    detinfo::DetectorTimings const& detTimings
    ) const;
  
  /**
   * @brief Converts a time into beam gate time.
   * @param time the time to be converted
   * @param detTimings detector timings helper
   * @return delay from the start of beam gate opening time to `time`
   * 
   * The time value is assumed to be in the reference as `fBeamGateReference`
   * and converted into the beam gate time scale.
   */
  nanoseconds toBeamGateTime(
    util::quantities::nanosecond time,
    detinfo::DetectorTimings const& detTimings
    ) const;
  
  
  /// Creates and returns a 1D histogram filled with `binnedContent`.
  TH1* makeHistogramFromBinnedContent(
    PlotSandbox_t& plots,
    std::string const& name, std::string const& title,
    BinnedContent_t const& binnedContent
    ) const;
  
  
  //@{
  /// Returns the time of the event in seconds from The Epoch.
  static double eventTimestampInSeconds(art::Timestamp const& time);
  static double eventTimestampInSeconds(art::Event const& event);
  //@}
  
  /// Fills an `EventAux_t` from the information found in the argument.
  static EventAux_t extractEventInfo(art::Event const& event);
  
  /// Converts a standard _art_ timestamp into an UTC time [ns]
  static std::uint64_t TimestampToUTC(art::Timestamp const& ts);
  
  /// Returns the ID of the cryostat the specified window belongs to.
  static geo::CryostatID WindowCryostat
    (icarus::trigger::WindowChannelMap::WindowInfo_t const& winfo)
    { return winfo.composition.cryoid; }

}; // icarus::trigger::TriggerSimulationOnGates



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Moves all the elements of `src` to the end of `dest`.
  /// The status of `src` after the call is undefined but destroyable.
  template <typename DestColl, typename SrcColl>
  DestColl& append(DestColl& dest, SrcColl&& src) {
    using std::begin, std::end;
    dest.insert(end(dest),
      std::make_move_iterator(begin(src)), std::make_move_iterator(end(src))
      );
    return dest;
  } // append()
  
  // ---------------------------------------------------------------------------
  
  /// Adds the numbers avoiding overflow. Only really safe if `T` == `U`.
  template <typename T, typename U, typename... Others>
  T cappedSum(T a, U&& b, Others&&... others) {
    
    static constexpr T max = std::numeric_limits<T>::max();
    
    auto const max_d = max - a;
    a = (b > max_d)? max: a + b;
    
    if constexpr(sizeof...(Others) > 0) {
      return cappedSum(a, std::forward<Others>(others)...);
    }
    else {
      return a;
    }
  } // cappedSum()
  
  
  // ---------------------------------------------------------------------------
  template <typename T>
  constexpr T bitMask(std::size_t i) {
    constexpr std::size_t nBits = sizeof(T) * 8;
    return (i >= nBits)? T{ 0 }: T{ 1 } << i;
  } // bitMask()
  
  
  // ---------------------------------------------------------------------------
  template <typename Coll>
  std::string join(Coll const& values, std::string const& sep = ", ") {
    std::string res;
    auto it = begin(values), vend = end(values);
    if (it == vend) return {};
    res = *it;
    while (++it != vend) {
      res += sep;
      res += *it;
    }
    return res;
  } // join()
  
  
  // ---------------------------------------------------------------------------
  
} // local namespace


//------------------------------------------------------------------------------
//---  icarus::trigger::TriggerSimulationOnGates
//------------------------------------------------------------------------------
icarus::trigger::TriggerSimulationOnGates::TriggerSimulationOnGates
  (Parameters const& config)
  : art::EDProducer       (config)
  // configuration
  , fPattern              (config().Pattern())
  , fBeamGateTag          (config().BeamGates())
  , fBeamGateReference    (config().BeamGateReference())
  , fBeamBits             (config().BeamBits())
  , fEmitEmpty            (config().EmitEmpty())
  , fExtraInfo            (config().ExtraInfo())
  , fTriggerOnTransition  (config().TriggerOnTransition())
  , fDeadTime             (config().DeadTime())
  , fTriggerDelay         (config().TriggerDelay())
  , fLVDSstatusDelay      (config().LVDSstatusDelay())
  , fRetriggeringMask     (bitMask<TriggerBits_t>(config().RetriggeringBit()))
  , fTriggerLogicMask     (mask(config().TriggerLogicBit()))
  , fCryostatZeroMask     (bitMask<TriggerBits_t>(config().CryostatFirstBit()))
  , fTriggerClock         (config().TriggerClock().value_or
                            (config().TriggerTimeResolution().value_or(25_ns)))
  , fOverlapTicks         (config().OverlapTicks())
  , fEventTimeBinning     (config().EventTimeBinning())
  , fLogCategory          (config().LogCategory())
  // services
  , fChannelMapCacheGuard{ "PMT" } // track the PMT cache only
  // internal and cached
  , fSaveLVDSbits         (fExtraInfo && !config().LVDSgatesTag().empty())
  , fWindowMapMan{ makeWindowMapManager() }
#if 0
  , fPlots{
    *art::ServiceHandle<art::TFileService>(),
    "", "requirement: " + fPattern.description()
    }
#endif //0
  , fEventPlotInfo{
        BinnedContent_t{ fEventTimeBinning }  // eventTimes
      , BinnedContent_t{                      // HWtrigTimeVsBeam
          fTriggerClock.convertInto
            <detinfo::timescales::trigger_time::interval_t>().value()
        }
      , BinnedContent_t{                      // triggerTimesVsHWtrig
          fTriggerClock.convertInto
            <detinfo::timescales::trigger_time::interval_t>().value()
        }
      , BinnedContent_t{                      // triggerTimesVsBeam
          fTriggerClock.convertInto
            <detinfo::timescales::trigger_time::interval_t>().value()
        }
    }
{
  //
  // simple parameter validation
  //
  if (config().TriggerTimeResolution()) {
    if (config().TriggerClock()
      && config().TriggerTimeResolution().value() != fTriggerClock
    ) {
      throw art::Exception{ art::errors::Configuration }
        << "Both `TriggerClock` and `TriggerTimeResolution` specified:"
           " remove the latter (which is deprecated).";
    }
    else {
      mf::LogWarning{ fLogCategory }
        << "TriggerSimulationOnGates: "
           " parameter 'TriggerTimeResolution' is deprecated:"
           " use `TriggerClock` instead."
        ;
    }
  } // if TriggerTimeResolution specified
  
  //
  // more complex parameter parsing
  //
  std::string const& discrModuleLabel = config().TriggerGatesTag();
  
  constexpr icarus::ADCsettings<double> ADCsettings;
  std::set<std::string> thresholds
    { cbegin(config().Thresholds()), cend(config().Thresholds()) };
  for (util::quantities::millivolt const& Vthr: config().VThresholds()) {
    std::string const thr = std::to_string(ADCsettings.to_ADC(Vthr));
    if (thresholds.count(thr) > 0) {
      throw art::Exception{ art::errors::Configuration }
        << "Threshold '" << Vthr << "' in " << config().VThresholds.name()
        << " was already specified as '" << thr << "' in "
        << config().Thresholds.name() << "!";
    }
    thresholds.insert(thr);
  } // for
  if (thresholds.empty()) {
    throw art::Exception{ art::errors::Configuration }
        << "No threshold tag configured in " << config().Thresholds.name()
        << " nor " << config().VThresholds.name() << "!";
  }
  
  std::string const& LVDSmoduleLabel = config().LVDSgatesTag();
  for (std::string const& threshold: thresholds) {
    InputInfo_t inputInfo;
    inputInfo.triggerGatesTag = art::InputTag{ discrModuleLabel, threshold };
    if (!LVDSmoduleLabel.empty())
      inputInfo.LVDSgatesTag = art::InputTag{ LVDSmoduleLabel, threshold };
    fInputInfo[threshold] = std::move(inputInfo);
  }
  
  // initialization of a vector of atomic is not as trivial as it sounds...
  fTriggerCount = std::vector<std::atomic<unsigned int>>(fInputInfo.size());
  std::fill(fTriggerCount.begin(), fTriggerCount.end(), 0U);
  
  if ((fTriggerClock <= 0_ns) && (fDeadTime <= 0_ns)) {
    throw art::Exception{ art::errors::Configuration }
      << "Either trigger resolution or dead time must be larger than 0.\n";
  }

  //
  // input data declaration
  //

  // trigger primitives
  for (InputInfo_t const& inputInfo: util::const_values(fInputInfo)) {
    icarus::trigger::TriggerGateReader<>{ inputInfo.triggerGatesTag }
      .declareConsumes(consumesCollector());
    if (fSaveLVDSbits) {
      icarus::trigger::TriggerGateReader<>{ inputInfo.LVDSgatesTag }
        .declareConsumes(consumesCollector());
    }
  } // for
  
  //
  // output data declaration
  //
  // keepThresholdName is true if we write instance name in output data products
  bool const keepThresholdName
    = config().KeepThresholdName().value_or(config().Thresholds().size() > 1);
  if (!keepThresholdName && (config().Thresholds().size() > 1)) {
    throw art::Exception(art::errors::Configuration)
      << config().KeepThresholdName.name()
      << " can be set to `true` only when a single threshold is specified ("
      << config().Thresholds.name() << " has " << config().Thresholds().size()
      << ")";
  }
  
  for (InputInfo_t const& inputInfo: util::const_values(fInputInfo)) {
    std::string const outputInstance
      = keepThresholdName? inputInfo.triggerGatesTag.instance(): "";
    produces<std::vector<raw::Trigger>>(outputInstance);
    if (fExtraInfo) produces<sbn::ExtraTriggerInfo>(outputInstance);
    fOutputInstances.push_back(outputInstance);
  }
  
  //
  // service access
  //
  if (fSaveLVDSbits) {
    try {
      fChannelMap
        = art::ServiceHandle<icarusDB::IICARUSChannelMap>()->provider();
    }
    catch (art::Exception const& e) { // print a more verbose error message
      if (e.categoryCode() == art::errors::ServiceNotFound) {
        mf::LogError(fLogCategory)
          << "PMT channel mapping service is required to save LVDS states.";
      }
      throw;
    }
    fChannelMapCacheGuard.setCache(*fChannelMap);
  }
  assert(!fSaveLVDSbits || fChannelMap);
  
  //
  // configuration dump
  //
  {
    mf::LogInfo log(fLogCategory);
    log << "\nConfigured " << fInputInfo.size() << " thresholds (ADC):";
    for (auto const& [ thresholdTag, inputInfo ]: fInputInfo) {
      log << "\n * " << thresholdTag << " (from '"
        << inputInfo.triggerGatesTag.encode() << "')";
    }
    log << "\nOther parameters:"
      << "\n * trigger time clock period: " << fTriggerClock
      << "\n * trigger response delay: " << fTriggerDelay
      << "\n * input beam gate: '" << fBeamGateTag.encode()
        << "', reference time: " << util::StandardSelectorFor<util::TimeScale>{}
          .get(fBeamGateReference).name()
      << "\n * emit triggers"
        << (fTriggerOnTransition
          ? " at the instant requirements become satisfied"
          : " any time the requirements are met"
          )
      ;
    if (fDeadTime == std::numeric_limits<nanoseconds>::max())
      log << "\n * only one trigger per beam gate (infinite dead time)";
    else {
      log << "\n * veto time after a trigger in one gate: " << fDeadTime
        << "\n * retriggering bit:";
      if (config().RetriggeringBit() >= NTriggerBits) {
        log << " none (#" << config().RetriggeringBit()
          << ", 0x" << std::hex << fRetriggeringMask << std::dec << ")";
      }
      else {
        log << " #" << config().RetriggeringBit() << " (0x"
          << std::hex << fRetriggeringMask << std::dec << ")";
      }
    } // if no dead time ... else
    log << "\n * bits for triggering cryostat:";
    if (config().CryostatFirstBit() >= NTriggerBits) 
      log << " disabled (#" << config().CryostatFirstBit() << ")";
    else {
      geo::GeometryCore const& geom = *lar::providerFrom<geo::Geometry>();
      unsigned int const firstBit = config().CryostatFirstBit();
      unsigned int const NCryostats = geom.Ncryostats();
      log << " #" << firstBit;
      if (NCryostats > 1) log << "-" << (firstBit + NCryostats - 1);
      log << " (0x" << std::hex << fCryostatZeroMask;
      if (NCryostats > 1)
        log << "-0x" << (fCryostatZeroMask << (NCryostats-1));
      log << std::dec << ")";
    }
    log << "\n * cryostat trigger overlap interval: <" << fOverlapTicks
      << " ticks = " << (fOverlapTicks*fTriggerClock);
    if (fExtraInfo) {
      log
        << "\n * will produce a sbn::ExtraTriggerInfo from the first gate"
        << "\n * trigger logic bit: " << join(names(fTriggerLogicMask))
        ;
    }
    if (fSaveLVDSbits) {
      assert(fExtraInfo);
      log << "\n     * will save frozen LVDS states from '" << LVDSmoduleLabel
        << "'";
      if (fLVDSstatusDelay != 0_ns)
        log << " after a delay of " << fLVDSstatusDelay;
    }
    
  } // local block
  
  
} // icarus::trigger::TriggerSimulationOnGates::TriggerSimulationOnGates()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::beginJob() {
  
  initializePlots();
  
} // icarus::trigger::TriggerSimulationOnGates::beginJob()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::beginRun(art::Run&) {
  
  // refresh the LVDS bit map
  if (fChannelMap && fChannelMapCacheGuard.update()) {
    // old versions of the database do not have all needed information;
    // we won't provide LVDS bits for those
    fLVDSmaps.emplace(*fChannelMap); // on failure, !fLVDSmaps.has_value()
    if (!fLVDSmaps->hasMap(icarus::trigger::LVDSbitMaps::Map::PMTpairs)) {
      mf::LogWarning(fLogCategory)
        << "PMT pair mapping could not be extracted."
        << " LVDS bits will not be simulated.";
    }
  }
  
} // icarus::trigger::TriggerSimulationOnGates::beginRun()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::produce(art::Event& event)
{
  
  //
  // prepare all the gates to run the simulation on
  //
  detinfo::DetectorTimings const detTimings {
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  
  std::vector<sim::BeamGateInfo> const& beamGates
    = event.getProduct<std::vector<sim::BeamGateInfo>>(fBeamGateTag);

  
  { // BEGIN local block
    mf::LogDebug log(fLogCategory);
    log << "Trigger simulation for " << beamGates.size() << " gates";
    if (!beamGates.empty()) {
      log << " ('" << fBeamGateTag.encode() << "'):";
      for (auto const& [iGate, gate]: util::enumerate(beamGates)) {
        log << "\n [" << iGate << "]  from " << gate.Start() << " to "
          << (gate.Start() + gate.Width()) << " ns (" << gate.Width() << "ns)";
      } // for
    } // if
  } // END local block
  
  
  //
  // run the simulation on each threshold in turn
  //
  mf::LogDebug log(fLogCategory); // this will print at the end of produce()
  log << "Event " << event.id() << ":";
  
  // FIXME these two operations should be atomic
  unsigned int const firstTriggerNumber
    = fTotalGates.fetch_add(beamGates.size());
  
  for (auto const& [ iThr, thrTag ]
    : util::enumerate(util::get_elements<0U>(fInputInfo))
  ) {
    
    std::vector<std::vector<WindowTriggerInfo_t>> const triggers
      = produceForThreshold
        (event, detTimings, beamGates, iThr, thrTag, firstTriggerNumber)
      ;
    
    log << "\n * threshold " << thrTag << ": ";
    icarus::ns::util::PassCounter gateResults;
    for (std::vector<WindowTriggerInfo_t> const& triggerInfo: triggers) {
      gateResults.add(std::any_of(triggerInfo.cbegin(), triggerInfo.cend(),
        std::mem_fn(&WindowTriggerInfo_t::info)));
    }
    log << gateResults.passed() << "/" << gateResults.total()
      << " gates triggered"; // ... at least once
    
  } // for
  
  //
  // event-level plots
  //
  plotEvent(event, detTimings, beamGates);
  
} // icarus::trigger::TriggerSimulationOnGates::produce()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::endJob() {
  
  finalizePlots();
  
  printSummary();
  
} // icarus::trigger::TriggerSimulationOnGates::endJob()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::initializePlots() {
  
#if 0
  //
  // overview plots with different settings
  //
  
  std::vector<std::string> thresholdLabels;
  thresholdLabels.reserve(size(fInputInfo));
  for (std::string thr: util::get_elements<0U>(fInputInfo))
    thresholdLabels.push_back(std::move(thr));
  
  auto const beamGate = makeMyBeamGate();
  mf::LogInfo(fLogCategory)
    << "Beam gate for plots: " << beamGate.asSimulationTime()
    << " (simulation time), " << beamGate.tickRange()
    << " (optical ticks)"
    ;

  //
  // Triggering efficiency vs. ADC threshold.
  //
  auto* NTriggers = fPlots.make<TH1F>(
    "NTriggers",
    "Number of triggering events"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";events",
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );
  util::ROOT::applyAxisLabels(NTriggers->GetXaxis(), thresholdLabels);
  
  auto* Eff = fPlots.make<TEfficiency>(
    "Eff",
    "Triggering pass fraction"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";trigger pass fraction",
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );
  // people are said to have earned hell for things like this;
  // but TEfficiency really does not expose the interface to assign labels to
  // its axes, which supposedly could be done had we chosen to create it by
  // histograms instead of directly as recommended.
  util::ROOT::applyAxisLabels
    (const_cast<TH1*>(Eff->GetTotalHistogram())->GetXaxis(), thresholdLabels);
  
  detinfo::timescales::optical_time_ticks const triggerResolutionTicks{
    icarus::ns::util::makeDetTimings().toOpticalTicks(fTriggerClock)
    };
  
  auto const& beamGateTicks = beamGate.tickRange();
  auto* TrigTime = fPlots.make<TH2F>(
    "TriggerTick",
    "Trigger time tick"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";events",
    static_cast<int>(std::ceil(beamGate.lengthTicks()/triggerResolutionTicks)),
    beamGateTicks.start().value(),
    icarus::ns::util::roundup
     (beamGateTicks.start() + beamGate.lengthTicks(), triggerResolutionTicks)
     .value(),
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );
  util::ROOT::applyAxisLabels(TrigTime->GetYaxis(), thresholdLabels);
  
  
  // we allow some fixed margin in the plot, just in case:
  constexpr microseconds beamPlotPadding { 4_us };
  
  using detinfo::timescales::trigger_time;
  
  // hardware trigger may happen at any place within the beam gate,
  // and in this plot range I want to include the full beam gate;
  // so I take a beam gate worth of time before the trigger time,
  // and as much after it; since this plot is relative to the hardware trigger,
  // the hardware trigger time itself is 0
  icarus::ns::util::BinningSpecs const HWtrigBinning = alignBinningTo(
    icarus::ns::util::BinningSpecs{
      (- beamGate.length() - beamPlotPadding).value(),
      (beamGate.length() + beamPlotPadding).value(),
      fTriggerClock.convertInto<trigger_time::interval_t>().value()
      },
      0.0
    );
  fPlots.make<TH2F>(
    "TriggerTimeVsHWTrig",
    "Time of the trigger"
      ";trigger time (relative to hardware trigger)  [ #mus ]"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";triggers",
    HWtrigBinning.nBins(), HWtrigBinning.lower(), HWtrigBinning.upper(),
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );

  icarus::ns::util::BinningSpecs const beamGateBinning = alignBinningTo(
    icarus::ns::util::BinningSpecs{
      (-beamPlotPadding).value(),
      (beamGate.length() + beamPlotPadding).value(),
      fTriggerClock.convertInto<trigger_time::interval_t>().value()
      },
      0.0
    );
  fPlots.make<TH2F>(
    "TriggerTimeVsBeamGate",
    "Time of the trigger"
      ";trigger time (relative to beam gate opening)  [ #mus ]"
      ";PMT discrimination threshold  [ ADC counts ]"
      ";triggers",
    beamGateBinning.nBins(), beamGateBinning.lower(), beamGateBinning.upper(),
    thresholdLabels.size(), 0.0, double(thresholdLabels.size())
    );

  // 
  // per-threshold plots; should this initialization be set into its own method?
  // 
  for (auto const& [ thr, info ]
    : util::zip(util::get_elements<0U>(fInputInfo), fThresholdPlots))
  {
    PlotSandbox_t& plots
      = fPlots.addSubSandbox("Thr" + thr, "Threshold: " + thr);
    
    plots.make<TGraph>(
      "TriggerTimeVsHWTrigVsBeam",
      "Time of the trigger: emulated vs. hardware"
        ";hardware trigger time (relative to beam gate opening)  [ #mus ]"
        ";emulated trigger time (relative to beam gate opening)  [ #mus ]"
      );
    
  } // for thresholds
  
  fThresholdPlots.resize(
    size(fInputInfo),
    {
      BinnedContent_t{ fEventTimeBinning },         // eventTimes
      BinnedContent_t{ HWtrigBinning.binWidth() },  // HWtrigTimeVsBeam
      BinnedContent_t{ HWtrigBinning.binWidth() },  // triggerTimesVsHWtrig
      BinnedContent_t{ beamGateBinning.binWidth() } // triggerTimesVsBeam
    }
    );
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::initializePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::finalizePlots() {
  
#if 0
  
  for (auto const& [ thr, info ]
    : util::zip(util::get_elements<0U>(fInputInfo), fThresholdPlots))
  {
    PlotSandbox_t& plots = fPlots.demandSandbox("Thr" + thr);
    makeThresholdPlots(thr, plots, info);
    if (plots.empty()) fPlots.deleteSubSandbox(plots.name());
  } // for thresholds
  
  makeEventPlots();
  
#endif // 0

} // icarus::trigger::TriggerSimulationOnGates::finalizePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::makeThresholdPlots(
  std::string const& threshold,
  PlotSandbox_t& plots,
  ThresholdPlotInfo_t const& plotInfo
) {
  
#if 0
  
  BinnedContent_t const* content;
  
  content = &(plotInfo.eventTimes);
  makeHistogramFromBinnedContent(plots,
    "TriggerTime",
    "Time of the triggered events"
      ";time"
      ";triggered events  [ / " + std::to_string(content->binWidth())
        + "\" ]",
    *content
    );
  
  content = &(plotInfo.HWtrigTimeVsBeam);
  makeHistogramFromBinnedContent(plots,
    "HWTrigVsBeam",
    "Time of the hardware trigger"
      ";trigger time (relative to beam gate)  [ #mus ]"
      ";events  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
  content = &(plotInfo.triggerTimesVsHWtrig);
  makeHistogramFromBinnedContent(plots,
    "TriggerTimeVsHWTrig",
    "Time of the trigger"
      ";trigger time (relative to hardware trigger)  [ #mus ]"
      ";triggers  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
  content = &(plotInfo.triggerTimesVsBeam);
  makeHistogramFromBinnedContent(plots,
    "TriggerTimeVsBeamGate",
    "Time of the trigger"
      ";trigger time (relative to beam gate opening)  [ #mus ]"
      ";triggers  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::makeThresholdPlots()



//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::makeEventPlots() {
  
#if 0
  
  BinnedContent_t const* content;
  
  content = &(fEventPlotInfo.eventTimes);
  makeHistogramFromBinnedContent(fPlots,
    "EventTime",
    "Time of the events"
      ";time"
      ";events  [ / " + std::to_string(content->binWidth())
        + "\" ]",
    *content
    );
  
  content = &(fEventPlotInfo.HWtrigTimeVsBeam);
  makeHistogramFromBinnedContent(fPlots,
    "HWTrigVsBeam",
    "Time of the hardware trigger"
      ";trigger time (relative to beam gate)  [ #mus ]"
      ";events  [ / " + std::to_string(content->binWidth())
        + " #mus ]",
    *content
    );
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::makeEventPlots()


//------------------------------------------------------------------------------
icarus::trigger::WindowTopologyManager
icarus::trigger::TriggerSimulationOnGates::makeWindowMapManager() const {
  return {
    *lar::providerFrom<geo::Geometry>(),
    art::ServiceHandle<geo::WireReadout const>()->Get(),
    fLogCategory + "_WindowMapManager"
    };
} // icarus::trigger::TriggerSimulationOnGates::makeWindowMapManager()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::produceForThreshold(
  art::Event& event,
  detinfo::DetectorTimings const& detTimings,
  std::vector<sim::BeamGateInfo> const& beamGates,
  std::size_t const iThr, std::string const& thrTag,
  unsigned int firstTriggerNumber
) -> std::vector<std::vector<WindowTriggerInfo_t>> {
  
//   auto& plotInfo = fThresholdPlots[iThr];
  
  //
  // get the input
  //
  art::InputTag const& dataTag = fInputInfo.at(thrTag).triggerGatesTag;
  auto const& gates = icarus::trigger::ReadTriggerGates(event, dataTag);
  
  // attempt to read LVDSgatesTag only if needed
  auto const* PMTpairGates = fSaveLVDSbits
    ? &(event.getProduct<std::vector<OpticalTriggerGateData_t>>
      (fInputInfo.at(thrTag).LVDSgatesTag))
    : nullptr
    ;
  
  // extract or verify the topology of the trigger windows
  if (fWindowMapMan(gates))
    fPatternAlg.emplace(*fWindowMapMan, fPattern, fLogCategory);
  assert(fPatternAlg);
  
  //
  // simulate the trigger response on all beam gates
  //
  std::vector<std::vector<WindowTriggerInfo_t>> allTriggerInfo; // one per gate
  auto triggers = std::make_unique<std::vector<raw::Trigger>>();
  std::unique_ptr<sbn::ExtraTriggerInfo> firstExtraInfo;
  unsigned int triggerNumber = firstTriggerNumber;
  
  EventAux_t const eventInfo = extractEventInfo(event);
  
  // after each trigger, skip
  // at least one trigger resolution unit, and at least the dead time
  optical_time_ticks const skipTicks
    = (fDeadTime == std::numeric_limits<nanoseconds>::max())
    ? std::numeric_limits<optical_time_ticks>::max()
    : detTimings.toOpticalTicks(std::max(fDeadTime, fTriggerClock));
  
  auto const nextTrigger
    = [this, &gates, &detTimings](optical_tick start, optical_tick stop)
    { return findFirstTrigger(start, stop, gates, detTimings); };
  
  for (sim::BeamGateInfo const& beamGate: beamGates) {
    std::vector<WindowTriggerInfo_t> triggerInfos;
    
    // relative to the beam gate time (also simulation time for MC);
    optical_tick start = detTimings.toOpticalTick(
      toElectronicsTime
        (util::quantities::nanosecond{ beamGate.Start() }, detTimings)
      );
    optical_tick const stop = start + detTimings.toOpticalTicks
      (util::quantities::nanosecond{ beamGate.Width() });
    
    // we want to start a tick earlier, to see if a trigger is already there:
    if (fTriggerOnTransition) start -= optical_time_ticks{ 1 };
    
    while (start < stop) {
      
      WindowTriggerInfo_t triggerInfo = nextTrigger(start, stop);
      if (!triggerInfo) break;
      
      // if detecting transitions, `start` will be kept at the closing of the
      // current trigger;
      // otherwise, after a dead time from the previous trigger
      if (fTriggerOnTransition && (triggerInfo.info.atTick() == start)) {
        // trigger at the very start means no transition occurred,
        // need to look for another one
        // shouldn't happen after the first trigger:
        assert(triggerInfos.empty());
        start = triggerInfo.info.endTick();
        triggerInfo = nextTrigger(start, stop);
        if (!triggerInfo) break;
      }
      
      // FIXME what do we do with statistics and plots?
//       plotInfo.eventTimes.add(eventTimestampInSeconds(event));
      
      mf::LogTrace(fLogCategory) << "Found a trigger at optical tick "
        << triggerInfo.info.atTick() << " (" << start
        << ") from window #" << triggerInfo.extra.windowIndex;
      
      // restart after some dead time (if may be infinite, so needs cap);
      // if triggering on transitions, skip at least to the end of this trigger
      start = fTriggerOnTransition
        ? std::max(triggerInfo.info.endTick(), cappedSum(start, skipTicks))
        : cappedSum(start, skipTicks);
      
      triggerInfos.push_back(std::move(triggerInfo));
      
    } // while
    
    if (!triggerInfos.empty()) ++fTriggerCount[iThr]; // keep the unique count
    
    //
    // fill the plots
    //
    plotTriggerResponse(iThr, thrTag, triggerInfos, detTimings);

    //
    // create and store the data product
    //
    auto [ gateTriggers, extraInfo ] = triggerInfoToTriggerData(
      detTimings, beamGate, eventInfo, triggerNumber++,
      triggerInfos, PMTpairGates
      );
    
    append(*triggers, std::move(gateTriggers));
    
    if (!firstExtraInfo && fExtraInfo) {
      firstExtraInfo
        = std::make_unique<sbn::ExtraTriggerInfo>(std::move(extraInfo));
    }
    
    allTriggerInfo.push_back(std::move(triggerInfos));
    
  } // for beam gates

  if (!firstExtraInfo && fExtraInfo) {
    // even with no gates at all our contract prescribes we emit extra info
    firstExtraInfo = std::make_unique<sbn::ExtraTriggerInfo>();
  }
  
  event.put(std::move(triggers), fOutputInstances[iThr]);
  if (fExtraInfo) event.put(std::move(firstExtraInfo), fOutputInstances[iThr]);
  
  return allTriggerInfo;
  
} // icarus::trigger::TriggerSimulationOnGates::produceForThreshold()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::plotEvent(
  art::Event const& event, detinfo::DetectorTimings const& detTimings,
  std::vector<sim::BeamGateInfo> const& gates
) {
  
#if 0
  
  detinfo::timescales::trigger_time const beamGateTime
    { detTimings.toTriggerTime(detTimings.BeamGateTime()) };
  
  fEventPlotInfo.eventTimes.add(eventTimestampInSeconds(event));
  fEventPlotInfo.HWtrigTimeVsBeam.add(-beamGateTime.value());
  
  // `gates` is currently unused; it may be used e.g. to show how many gates
  // were tested in each event
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::plotEvent()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::plotTriggerResponse(
  std::size_t iThr, std::string const& thrTag,
  std::vector<WindowTriggerInfo_t> const& triggerInfos,
  detinfo::DetectorTimings const& detTimings
) {
  
#if 0
  
  bool const fired = triggerInfo.info.fired();
  
  fPlots.demand<TEfficiency>("Eff").Fill(fired, iThr);
  
  if (fired) {
    using namespace detinfo::timescales;
    
    // time of the beam gate in hardware trigger time scale
    trigger_time const beamGateTime
      { detTimings.toTriggerTime(detTimings.BeamGateTime()) };
    
    optical_tick const thisTriggerTick { triggerInfo.info.atTick() };
    trigger_time const thisTriggerTimeVsHWtrig
      { detTimings.toTriggerTime(thisTriggerTick) };
    time_interval const thisTriggerTimeVsBeamGate
      { thisTriggerTimeVsHWtrig - beamGateTime };
    
    mf::LogTrace(fLogCategory)
      << "Trigger " << fPattern.tag() << " at tick " << thisTriggerTick
      << " (" << thisTriggerTimeVsHWtrig << " vs. HW trigger, "
      << thisTriggerTimeVsBeamGate << " vs. beam gate)"
      ;
    
    fPlots.demand<TH1>("NTriggers").Fill(iThr);
    fPlots.demand<TH2>("TriggerTick").Fill(thisTriggerTick.value(), iThr);
    fPlots.demand<TH2>("TriggerTimeVsHWTrig").Fill
      (thisTriggerTimeVsHWtrig.value(), iThr);
    fPlots.demand<TH2>("TriggerTimeVsBeamGate").Fill
      (thisTriggerTimeVsBeamGate.value(), iThr);
    
    PlotSandbox_t& plots{ fPlots.demandSandbox("Thr" + thrTag) };
//     plots.demand<TGraph>("TriggerTimeVsHWTrigVsBeam").AddPoint( // ROOT 6.24?
    TGraph& graph = plots.demand<TGraph>("TriggerTimeVsHWTrigVsBeam");
    graph.SetPoint(graph.GetN(),
      -beamGateTime.value(), thisTriggerTimeVsBeamGate.value()
      );
    
    ThresholdPlotInfo_t& plotInfo { fThresholdPlots[iThr] };
    plotInfo.HWtrigTimeVsBeam.add(-beamGateTime.value());
    plotInfo.triggerTimesVsHWtrig.add(thisTriggerTimeVsHWtrig.value());
    plotInfo.triggerTimesVsBeam.add(thisTriggerTimeVsBeamGate.value());
    
  }
  
#endif // 0
  
} // icarus::trigger::TriggerSimulationOnGates::plotTriggerResponse()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerSimulationOnGates::printSummary() const {
  
  //
  // summary from our internal counters
  //
  mf::LogInfo log(fLogCategory);
  log
    << "Summary of triggers for " << fTriggerCount.size()
    << " thresholds (ADC) with pattern: " << fPattern.description()
    ;
  for (auto const& [ count, thr ]
    : util::zip(fTriggerCount, util::get_elements<0U>(fInputInfo)))
  {
    log << "\n  threshold " << thr
      << ": " << count;
    if (fTotalGates > 0U) {
      log << "/" << fTotalGates
        << " (" << (double(count) / fTotalGates * 100.0) << "%)";
    }
    else log << " gates triggered";
  } // for
  
} // icarus::trigger::TriggerSimulationOnGates::printSummary()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::toElectronicsTime(
  util::quantities::nanosecond time, detinfo::DetectorTimings const& detTimings
) const
  -> detinfo::timescales::electronics_time
{
  // currently (LArSoft v09_77_00) `detinfo::DetectorTimings` does not support
  // beam gate timescale conversion, so we need to do it "by hand" from...
  // electronics time, as usual
  
  switch (fBeamGateReference) {
    case util::TimeScale::Electronics:
      return detinfo::timescales::electronics_time{ time };
    case util::TimeScale::BeamGate:
      return detTimings.BeamGateTime() + nanoseconds{ time };
    case util::TimeScale::Trigger:
      return detTimings.toElectronicsTime
        (detinfo::timescales::trigger_time{ time });
      break;
    case util::TimeScale::Simulation:
      return detTimings.toElectronicsTime
        (detinfo::timescales::simulation_time{ time });
    default:
      throw art::Exception{ art::errors::Configuration }
        << "Conversion of times from reference '"
        << util::StandardSelectorFor<util::TimeScale>{}
          .get(fBeamGateReference).name()
        << "' not supported.\n";
  } // switch
  
} // icarus::trigger::TriggerSimulationOnGates::toElectronicsTime()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::toBeamGateTime(
  util::quantities::nanosecond time, detinfo::DetectorTimings const& detTimings
) const
  -> nanoseconds
{
  // currently (LArSoft v09_77_00) `detinfo::DetectorTimings` does not support
  // beam gate timescale conversion, so we need to do it "by hand" from...
  // electronics time, as usual
  
  return fBeamGateReference == util::TimeScale::BeamGate
    ? time: toElectronicsTime(time, detTimings) - detTimings.BeamGateTime();
  
} // icarus::trigger::TriggerSimulationOnGates::toBeamGateTime()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::extractEventInfo
  (art::Event const& event) -> EventAux_t
{
  return {
      TimestampToUTC(event.time()) // time
    , event.event()                // event
    };
} // icarus::trigger::TriggerSimulationOnGates::extractEventInfo()


//------------------------------------------------------------------------------
std::uint64_t icarus::trigger::TriggerSimulationOnGates::TimestampToUTC
  (art::Timestamp const& ts)
{
  return static_cast<std::uint64_t>(ts.timeHigh())
    + static_cast<std::uint64_t>(ts.timeLow()) * 1'000'000'000ULL;
} // icarus::trigger::TriggerSimulationOnGates::TimestampToUTC()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::findFirstTrigger(
    optical_tick start, optical_tick stop,
  icarus::trigger::SlidingWindowPatternAlg::TriggerGates_t const& gates,
  detinfo::DetectorTimings const& detTimings
) const -> WindowTriggerInfo_t {
  
  detinfo::DetectorClocksData const& detClocks = detTimings.clockData();
  
  icarus::trigger::ApplyBeamGateClass const applyBeamGate
   = makeApplyBeamGate(start, stop - start, detClocks, fLogCategory);
  
  mf::LogTrace(fLogCategory) << "Applying gate: " << applyBeamGate;
  
  return fPatternAlg->simulateResponse(applyBeamGate.applyToAll(gates));
  
} // icarus::trigger::TriggerSimulationOnGates::findFirstTrigger()


//------------------------------------------------------------------------------
std::pair<std::vector<raw::Trigger>, sbn::ExtraTriggerInfo>
icarus::trigger::TriggerSimulationOnGates::triggerInfoToTriggerData(
  detinfo::DetectorTimings const& detTimings,
  sim::BeamGateInfo const& beamGate,
  EventAux_t const& eventInfo,
  unsigned int triggerNumber, std::vector<WindowTriggerInfo_t> const& info,
  std::vector<OpticalTriggerGateData_t> const* PMTpairGates
) const {
  
  static constexpr auto NoTime = std::numeric_limits<unsigned int>::max();
  struct PreviousTriggerInfo_t {
    unsigned int beamToTriggerTicks = NoTime;
    geo::CryostatID cryo;
    constexpr operator bool() const { return beamToTriggerTicks != NoTime; }
    constexpr bool operator!() const { return !bool(*this); }
  };
  
  detinfo::timescales::electronics_time const beamTime
    = detTimings.BeamGateTime() + nanoseconds{ beamGate.Start() };
  TriggerBits_t const beamBits
    = fBeamBits.value_or(makeTriggerBits(beamGate.BeamType()));
  
  std::vector<raw::Trigger> triggers;
  
  sbn::ExtraTriggerInfo extraInfo;
  assert(extraInfo.triggerID == sbn::ExtraTriggerInfo::NoID);
  assert(extraInfo.triggerLocationBits == 0);
  assert(extraInfo.triggerTimestamp == sbn::ExtraTriggerInfo::NoTimestamp);
  // these are fixed by the gate and set this way whether trigger fired or not:
  extraInfo.beamGateTimestamp = eventInfo.time;
  extraInfo.triggerType       = sbn::bits::triggerType::Majority;
  extraInfo.sourceType        = beamTypeToTriggerSource(beamGate.BeamType());
  extraInfo.gateID            = eventInfo.event;
  extraInfo.gateCount         = eventInfo.event;
  
  PreviousTriggerInfo_t prevTrigger;
  for (WindowTriggerInfo_t const& trInfo: info) {
    
    if (trInfo.info.fired()) { // trigger fired
      
      detinfo::timescales::electronics_time const triggerTime
        = detTimings.toElectronicsTime(trInfo.info.atTick());
      
      // include the delay from timestamping (and all other fixed delays)
      detinfo::timescales::electronics_time const timestampedTriggerTime
        = triggerTime + fTriggerDelay;
      
      // find the location and set the bits accordingly
      geo::CryostatID const triggeringCryo
        = WindowCryostat(fWindowMapMan->info(trInfo.extra.windowIndex));
      TriggerBits_t const sourceMask = triggeringCryo.isValid
        ? (fCryostatZeroMask << triggeringCryo.Cryostat): 0U;
      
      // round to the next tick (because that's what the hardware must do):
      nanoseconds const beamToTrigger = triggerTime - beamTime;
      auto const beamToTriggerTicks
        = static_cast<unsigned int>(std::ceil(beamToTrigger / fTriggerClock));
      
      bool isFirstTrigger = triggers.empty(); // might still be merging tough
      
      // if the previous trigger is from a different cryostat and close in time,
      if (prevTrigger && (prevTrigger.cryo != triggeringCryo)
        && (beamToTriggerTicks - prevTrigger.beamToTriggerTicks <= fOverlapTicks)
      ) { // ... merge
        assert(!triggers.empty());
        isFirstTrigger = (triggers.size() <= 1);
        
        prevTrigger = PreviousTriggerInfo_t{}; // overlap only one trigger
        
        raw::Trigger prevTrigger = triggers.back();
        triggers.back() = raw::Trigger{
          prevTrigger.TriggerNumber(),             // counter
          prevTrigger.TriggerTime(),               // trigger time
          prevTrigger.BeamGateTime(),              // beam gate
          prevTrigger.TriggerBits() | sourceMask,  // bits
          };
        
      } // if merging trigger
      else { // new trigger
        prevTrigger = { beamToTriggerTicks, triggeringCryo };
        
        TriggerBits_t bits = beamBits | sourceMask;
        if (!isFirstTrigger) bits |= fRetriggeringMask;
        
        triggers.emplace_back(
          triggerNumber,                   // counter
          double(timestampedTriggerTime),  // trigger time
          double(beamTime),                // beam gate
          bits                             // bits
          );
        
        if (fExtraInfo && isFirstTrigger) {
          nanoseconds const relTrigTime = timestampedTriggerTime - beamTime;
          
          extraInfo.triggerTimestamp  = eventInfo.time;
          extraInfo.beamGateTimestamp = extraInfo.triggerTimestamp
            - static_cast<std::int64_t>(std::round(relTrigTime.value()));
          
          extraInfo.triggerID    = triggerNumber;
          extraInfo.triggerCount = triggerNumber;
          
          if (fSaveLVDSbits) {
            assert(PMTpairGates);
            optical_tick const LVDSfreezeTick
              = detTimings.toOpticalTick(triggerTime + fLVDSstatusDelay);
            LVDSbitArrays_t const LVDSbits
              = extractLVDSstatus(LVDSfreezeTick, *PMTpairGates);
            for (std::size_t const cryo
              : util::counter(sbn::ExtraTriggerInfo::MaxCryostats)
            ) {
              extraInfo.cryostats[cryo].LVDSstatus = LVDSbits[cryo];
            } // for cryostat
          } // if saving LVDS bits
        }
      } // if new trigger
      
      // more extra information (only if this is the first trigger)
      if (fExtraInfo && isFirstTrigger) {
        extraInfo.triggerLocationBits |= cryoIDtoTriggerLocation(triggeringCryo);
        if (triggeringCryo.isValid) {
          sbn::ExtraTriggerInfo::CryostatInfo& cryoInfo = extraInfo.cryostats[
            triggeringCryo.Cryostat == 0
            ? sbn::ExtraTriggerInfo::EastCryostat
            : sbn::ExtraTriggerInfo::WestCryostat
            ];
          if (cryoInfo.triggerCount++ == 0) { // first trigger in the cryostat
            cryoInfo.triggerLogicBits = fTriggerLogicMask;
            cryoInfo.beamToTrigger
              = static_cast<unsigned int>(beamToTrigger.value());
          }
        } // cryostat info
      } // if extra info
      
    }
    else { // trigger did not fire
      
      // nothing to do for raw::Trigger: we just don't produce any
      
    } // if fired... else
    
  } // for
  
  if (triggers.empty() && fEmitEmpty) {
    
    triggers.emplace_back(
      triggerNumber,                          // counter
      std::numeric_limits<double>::lowest(),  // trigger time
      double(beamTime),                   // beam gate in electronics time scale
      0                                       // bits
      );
    
  } // if
  
  return { std::move(triggers), std::move(extraInfo) };
  
} // icarus::trigger::TriggerSimulationOnGates::triggerInfoToTriggerData()


//------------------------------------------------------------------------------
sbn::triggerSourceMask
icarus::trigger::TriggerSimulationOnGates::makeTriggerBits
  (sim::BeamGateInfo const& beamInfo)
  { return mask(beamTypeToTriggerSource(beamInfo.BeamType())); }


//------------------------------------------------------------------------------
sbn::triggerSource
icarus::trigger::TriggerSimulationOnGates::beamTypeToTriggerSource
  (sim::BeamType_t beamType)
{
  
  switch (beamType) {
    case sim::kBNB: return sbn::bits::triggerSource::BNB;
    case sim::kNuMI: return sbn::bits::triggerSource::NuMI;
    case sim::kUnknown:
    default:
      return sbn::bits::triggerSource::Unknown;
  } // switch
  
} // icarus::trigger::TriggerSimulationOnGates::beamTypeToTriggerSource()


//------------------------------------------------------------------------------
sbn::bits::triggerLocationMask
icarus::trigger::TriggerSimulationOnGates::cryoIDtoTriggerLocation
  (geo::CryostatID const& cid)
{
  
  if      (cid == geo::CryostatID{ 0 })
    return mask(sbn::bits::triggerLocation::CryoEast);
  else if (cid == geo::CryostatID{ 1 })
    return mask(sbn::bits::triggerLocation::CryoWest);
  else
    return mask(sbn::bits::triggerLocation::NBits);
  
} // icarus::trigger::TriggerSimulationOnGates::cryoIDtoTriggerLocation()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::extractLVDSstatus(
  optical_tick triggerTick,
  std::vector<OpticalTriggerGateData_t> const& PMTgates
) const -> LVDSbitArrays_t {
  /*
   * For each "logic" bit of each cryostat and PMT wall,
   * we ask the mapping which channels that bit is associated to,
   * and find the gate with those gates.
   * The state of that gate at `triggerTick` becomes the value of the bit.
   */
  
  LVDSbitArrays_t LVDSbits;
  for (auto& cryoBits: LVDSbits) cryoBits.fill(0);
  
  assert(fLVDSmaps);
  if (!fLVDSmaps->hasMap(icarus::trigger::LVDSbitMaps::Map::PMTpairs))
    return LVDSbits;
  
  // map channel -> gate containing it
  std::vector<OpticalTriggerGateData_t const*> const gateMap
    = makeGateMap(PMTgates);
  
  auto const mask
    = [](auto bitNo){ return 1ULL << static_cast<std::uint64_t>(bitNo); };
  
  using icarus::trigger::PMTpairBitID;
  for (std::size_t const cryostat:
    { sbn::ExtraTriggerInfo::EastCryostat, sbn::ExtraTriggerInfo::WestCryostat }
  ) {
    
    for (std::size_t const PMTwall:
      { sbn::ExtraTriggerInfo::EastPMTwall, sbn::ExtraTriggerInfo::WestPMTwall }
    ) {
      
      std::uint64_t bits = 0;
      
      for (auto const bit: util::counter<PMTpairBitID::StatusBit_t>(64)) {
        
        PMTpairBitID const bitID{ cryostat, PMTwall, bit };
        std::vector<raw::Channel_t> const& channels
          = fLVDSmaps->bitSource(bitID).channels;
        
        if (channels.empty()) continue; // bit not mapped to anything
        
        // find the relevant gate
        raw::Channel_t const matchingChannel = channels.front();
        OpticalTriggerGateData_t const* LVDSgate = gateMap.at(matchingChannel);
        
        if (!LVDSgate) {
          mf::LogWarning log{ fLogCategory };
          log << "LVDS status " << bitID << " (" << channels.size()
            << " channels: ";
          for (raw::Channel_t const channel: channels) log << " " << channel;
          log << ") is not covered by any LVDS gate.";
          continue;
        }
        
        // check the full matching
        if (!matchChannelList(*LVDSgate, channels)) {
          art::Exception e{ art::errors::Unknown };
          e << "LVDS status " << bitID << " (" << channels.size() <<
            " channels: ";
          for (raw::Channel_t const channel: channels) e << " " << channel;
          e << ") matched to a LVDS gate from a different list of channels (";
          for (raw::Channel_t const channel: LVDSgate->channels())
            e << " " << channel;
          e << "); matching channel: " << matchingChannel << ".\n";
          throw e;
        } // if not matching channels
        
        // extract the value
        bool const bitValue = LVDSgate->isOpen(triggerTick.value());
        if (bitValue) bits |= mask(bit);
        
      } // bit
      
      LVDSbits[cryostat][PMTwall] = bits;
      
    } // PMT wall
  } // cryostat
  
  return LVDSbits;
} // icarus::trigger::TriggerSimulationOnGates::extractLVDSstatus()


//------------------------------------------------------------------------------
bool icarus::trigger::TriggerSimulationOnGates::matchChannelList(
  OpticalTriggerGateData_t const& gate,
  std::vector<raw::Channel_t> const& channels
) {
  
  // brute force
  auto const channelIsInList = [&chList=channels](raw::Channel_t ch)
    { return std::find(chList.begin(), chList.end(), ch) != chList.end(); };
  for (raw::Channel_t channel: gate.channels())
    if (!channelIsInList(channel)) return false;
  return true;
  
} // icarus::trigger::TriggerSimulationOnGates::matchChannelList()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerSimulationOnGates::makeGateMap
  (std::vector<OpticalTriggerGateData_t> const& gates)
  -> std::vector<OpticalTriggerGateData_t const*>
{
  // map channel -> gate containing it
  std::vector<OpticalTriggerGateData_t const*> gateMap
    (icarus::trigger::LVDSbitMaps::NChannels, nullptr);
  for (OpticalTriggerGateData_t const& gate: gates) {
    for (raw::Channel_t const channel: gate.channels()) {
      gateMap.at(channel) = &gate;
    }
  }
  return gateMap;
} // icarus::trigger::TriggerSimulationOnGates::makeGateMap()


//------------------------------------------------------------------------------
TH1*
icarus::trigger::TriggerSimulationOnGates::makeHistogramFromBinnedContent(
  PlotSandbox_t& plots,
  std::string const& name, std::string const& title,
  BinnedContent_t const& binnedContent
) const {
  
  if (binnedContent.empty()) return nullptr;
  
  TH1* hist = plots.make<TH1F>(
    name, title,
    binnedContent.nBins(), binnedContent.min(), binnedContent.max()
    );
  
  // directly transfer the content bin by bin
  unsigned int total = 0U;
  for (auto [ iBin, count ]: util::enumerate(binnedContent)) {
    hist->SetBinContent(iBin + 1, count);
    total += count;
  }
  hist->SetEntries(static_cast<double>(total));
  return hist;
} // icarus::trigger::TriggerSimulationOnGates::makeHistogramFromBinnedContent


//------------------------------------------------------------------------------
double icarus::trigger::TriggerSimulationOnGates::eventTimestampInSeconds
  (art::Timestamp const& time)
{
  // high value: seconds from the Epoch (Jan 1, 1970 UTC?);
  // low value: nanoseconds after that the start of that second
  return static_cast<double>(time.timeHigh())
    + static_cast<double>(time.timeHigh()) * 1e-9;
} // icarus::trigger::TriggerSimulationOnGates::eventTimestampInSeconds()


//------------------------------------------------------------------------------
double icarus::trigger::TriggerSimulationOnGates::eventTimestampInSeconds
  (art::Event const& event)
  { return eventTimestampInSeconds(event.time()); }


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::TriggerSimulationOnGates)


//------------------------------------------------------------------------------
