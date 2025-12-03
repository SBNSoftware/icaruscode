/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.cxx
 * @brief  Provides a adder simulation manager.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.h
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.h"

// ICARUS/SBN libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelSimulator.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/Indenter.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateDataIteration.h" // intervalsOver()
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "icaruscode/PMT/Trigger/Algorithms/AmplitudeAdderCalibration.h"
#include "sbnalg/Utilities/ToStr.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/TriggerGateData.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/spacetime.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()

// // ROOT libraries
#include "TDirectory.h"

// C/C++ standard libraries
#include <algorithm> // std::sort(), std::binary_search()
#include <iterator> // std::distance()
#include <sstream>
#include <stdexcept> // std::runtime_error()
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace {
  
  /// Returns a copy of the collection `coll`, sorted by `std::sort()`.
  template <typename Coll>
  Coll sorted(Coll coll) {
    Coll sortedColl{ coll };
    std::sort(begin(sortedColl), end(sortedColl));
    return sortedColl;
  }
  
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


// -----------------------------------------------------------------------------
icarus::trigger::AdderSignalSimulation::AdderSignalSimulation(
  Config const& config,
  adder::types::ADCsettings_t ADCsettings /* = {} */,
  icarus::trigger::AdderChannelMap const* adderChannels /* = nullptr */,
  FilteringCircuit const* reshapingFilter /* = nullptr */,
  AdderCalibrationDatabase::RunCalibration const* calibration /* = nullptr */
  )
  // base classes
  : icarus::ns::util::mfLoggingClass{ config.LogCategory() }
  // configuration parameters
  , fMissingChannels { sorted(config.MissingChannels()) }
  , fSplitToAdders   { config.SplitToAdders() }
  , fAdderBaseline   { config.AdderBaseline() }
  , fADCsettings     { std::move(ADCsettings) }
{
  
  setup(adderChannels, reshapingFilter, calibration);
  
} // icarus::trigger::AdderSignalSimulation::AdderSignalSimulation()


// -----------------------------------------------------------------------------
void icarus::trigger::AdderSignalSimulation::setup(
  icarus::trigger::AdderChannelMap const* adderChannels,
  FilteringCircuit const* reshapingFilter,
  AdderCalibrationDatabase::RunCalibration const* calibration
) {
  fAdderChannels = adderChannels;
  fReshapingFilter = reshapingFilter;
  fCalibration = calibration;
}


// -----------------------------------------------------------------------------
std::string icarus::trigger::AdderSignalSimulation::configurationDump
  (std::string const& indent, std::string const& firstIndent) const
{
  // the settings from `setup()` call are not dumped
  
  std::ostringstream log;
  details::Indenter nextLine{ indent, firstIndent };
  
  log << nextLine << "AdderSignalSimulation algorithm configuration:"
    << nextLine << " * PMT signal fraction: " << (1-fSplitToAdders) << " to digitizer/"
      << fSplitToAdders << " to adders"
    ;
  if (fMissingChannels.empty())
    log << nextLine << " * no missing channels configured";
  else {
    log << nextLine << " * configured " << fMissingChannels.size() << " missing channels:";
    for (raw::Channel_t const channel: fMissingChannels)
      log << " " << channel;
  }

  return log.str();
} // icarus::trigger::AdderSignalSimulation::configurationDump()


// -----------------------------------------------------------------------------
std::vector<raw::OpDetWaveform>
icarus::trigger::AdderSignalSimulation::simulate(
  std::vector<WaveformWithBaseline> const& waveformInfo,
  TDirectory* plotDir /* = nullptr */
) const {
  
  using namespace util::quantities::time_literals;
  
  if (!fAdderChannels) {
    throw std::runtime_error
      { "AdderSignalSimulation::simulate(): adder channel map not specified!" };
  }
  
  WaveformsByChannel_t const waveformInfoByChannel
    = groupByChannel(waveformInfo);
  
  std::vector<raw::OpDetWaveform> adderWaveforms;
  
  auto makeChannelPlotDir = [plotDir](AdderChannelID channel) -> TDirectory*
    {
      if (!plotDir) return nullptr;
      using util::toStr;
      std::string const dirName
        = toStr << plotDir->GetName() << "_CH" << channel;
      std::string const dirTitle
        = toStr << plotDir->GetTitle() << " channel " << channel;
      return plotDir->mkdir(dirName.c_str(), dirTitle.c_str());
    };
  
  for (auto const& [ adderChannel, PMTs ]: util::const_values(*fAdderChannels)) {
    
    TDirectory* const chPlotDir = makeChannelPlotDir(adderChannel);
    
    // TODO needs more configuration and more centralization
    AdderChannelSimulator const adderSimulator{
      adderChannel, fADCsettings, fSplitToAdders,
      fReshapingFilter, fCalibration, chPlotDir, logCategory() };
    
    mfLogTrace() << "Processing adder channel " << adderChannel;
    
    // std::vector<std::vector<WaveformWithBaseline const*>> inputWaveforms
    // std::set<raw::Channel_t> requiredChannels
    auto const& [ inputWaveforms, requiredChannels ]
      = collectAdderWindowWaveforms(waveformInfoByChannel, PMTs, adderChannel);
    
    
    // identify all the time intervals where we can operate
    std::vector<std::pair<TimeInterval_t, std::vector<WaveformWithBaseline const*>>>
    timeIntervals
      = findOverlappingTimeIntervals(inputWaveforms, requiredChannels);
    
    for (auto const& [ timeInterval, inputWaveforms ]: timeIntervals) {
      
      auto log = mfLogTrace();
      log << "Processed adder channel " << adderChannel
        << " time interval " << timeInterval;
      // AdderChannelSimulator::AdderSignal_t adderWaveformInfo
      // std::vector<WaveformWithBaseline const*> contributions
      auto const& [ adderWaveformInfo, contributions ]
        = adderSimulator.simulate(timeInterval, inputWaveforms);
      
      unsigned int const nClippedInputs = adderWaveformInfo.nClipped;
      
      raw::OpDetWaveform adderWaveform = adderWaveformInfo.makeOpDetWaveform
        (adderSimulator.ADCsettings(), fAdderBaseline);
      
      // AmplitudeAdderCalibration is not necessarily used for calibration;
      // here its computePeakAmplitude() is invoked for convenience
      log << " => adder channel " << adderWaveform
        << " from " << contributions.size() << "/" << inputWaveforms.size()
        << " PMT waveforms; peak amplitude: "
        << AmplitudeAdderCalibration::computePeakAmplitude(adderWaveform.Waveform());
      if (nClippedInputs > 0)
        log << " (" << nClippedInputs << " input waveforms could be clipped)";
      
      adderWaveforms.push_back(std::move(adderWaveform));
      
    } // for time intervals
    
  } // for windows
  
  return adderWaveforms;
} // icarus::trigger::AdderSignalSimulation::simulate()


//------------------------------------------------------------------------------
auto icarus::trigger::AdderSignalSimulation::groupByChannel
  (std::vector<icarus::trigger::WaveformWithBaseline> const& waveformInfo)
  -> WaveformsByChannel_t
{
  auto waveformChannel = [](icarus::trigger::WaveformWithBaseline const& wi)
    { return static_cast<std::size_t>(wi.waveform().ChannelNumber()); };
  return icarus::ns::util::GroupByIndex{ waveformInfo, waveformChannel };
} // icarus::trigger::AdderSignalSimulation::groupByChannel()


//------------------------------------------------------------------------------
auto icarus::trigger::AdderSignalSimulation::collectAdderWindowWaveforms(
  WaveformsByChannel_t const& waveformInfoByChannel,
  icarus::trigger::AdderChannelInfo_t::PMTchannelList_t const& windowChannels,
  AdderChannelID adderChannel /* = NoChannel */
) const 
  -> std::pair<std::vector<std::vector<WaveformWithBaseline const*>>, std::set<raw::Channel_t>>
{
  
  auto log = mfLogTrace();
  log << "Collecting waveforms for adder channel " << adderChannel << " from "
    << windowChannels.size() << " channels:";
  
  std::set<raw::Channel_t> requiredChannels;
  std::vector<std::vector<WaveformWithBaseline const*>> inputWaveforms;
  for (raw::Channel_t const ch: windowChannels) {
    
    // in data and MC with PMT readout simulation, all channels are present;
    // in the other cases, all photoelectron signals are present
    // (actually, only above a certain low threshold defined in `ThresholdADC`
    //  configuration parameter of `icarus::opdet::PMTsimulationAlg`)
    // TODO treatment of "missing channels" needs to be decided and checked;
    //      in the current code, a missing channel is not required;
    //      if this is changed, update class documentation too
    bool const bMissing = isMissingChannel(ch);
    if (bMissing) log << " (!" << ch << ")";
    else {
      log << " " << ch;
      requiredChannels.insert(ch);
    }
    if (auto nCh = waveformInfoByChannel[ch].size(); nCh != 1)
      log << " [x" << nCh << "]";
    
    inputWaveforms.push_back(waveformInfoByChannel[ch]);
    
  } // for
  
  return { std::move(inputWaveforms), std::move(requiredChannels) };
} // icarus::trigger::AdderSignalSimulation::collectAdderWindowWaveforms()


//------------------------------------------------------------------------------
auto icarus::trigger::AdderSignalSimulation::findOverlappingTimeIntervals(
  std::vector<std::vector<WaveformWithBaseline const*>> const& waveformsPerChannel,
  std::set<raw::Channel_t> const& requiredChannels
) const -> std::vector<std::pair<TimeInterval_t, std::vector<WaveformWithBaseline const*>>>
{
  
  /* The promise:
   *
   * For each returned item, the time interval is included together with all the
   * input waveforms that have some overlap with it.
   * By construction, in each time interval there will be only one waveform from
   * each of the required channels.
   * The behaviour in corner cases where there are non-required channels that
   * overlap a single interval is not well defined.
   * 
   * It may happen that the same waveform contributes to different intervals.
   */
  
  //
  // overlap the coverage of all mandatory channels
  // (and collect the non-mandatory ones)
  //
  // type to track overlaps in time (it was designed to have time like ticks):
  using Levels_t = icarus::trigger::TriggerGateData<
    electronics_time, util::quantities::concepts::interval_of<electronics_time>
    >;
    
  Levels_t overlapTimeline;
  std::vector<std::vector<WaveformWithBaseline const*> const*> optionalChannels;
  
  for (std::vector<WaveformWithBaseline const*> const& channelWaveforms
    : waveformsPerChannel)
  {
    if (channelWaveforms.empty()) continue;
    
    auto log = mfLogTrace();
    
    raw::Channel_t const channel
      = channelWaveforms.front()->waveform().ChannelNumber(); // any will do
    if (!requiredChannels.count(channel)) {
      log << "Channel " << channel << " is not required.";
      optionalChannels.push_back(&channelWaveforms);
      continue;
    }
    
    log << "Channel " << channel
      << " (" << channelWaveforms.size() << " waveforms):";
    
    for (WaveformWithBaseline const* wnb: channelWaveforms) {
      assert(wnb);
      assert(wnb->waveformPtr());
      
      raw::OpDetWaveform const& waveform = wnb->waveform();
      
      TimeInterval_t const interval = waveformInterval(waveform);
      log << " " << interval << " (" << interval.duration() << ")";
      assert(!interval.empty());
      
      overlapTimeline.openBetween(interval.start, interval.stop);
      
    } // for waveform in channel
    
  } // for all channels
  
  mfLogTrace() << "Timeline: " << overlapTimeline;
  
  // shortcut if there is no overlap at all (there should be some though)
  if (overlapTimeline.findOpen(requiredChannels.size()) == overlapTimeline.MaxTick) {
    mfLogTrace() << "No overlap between all " << requiredChannels.size()
      << " required channels.";
    return {};
  }
  
  //
  // for each interval, collect all the contributing waveforms
  //
  std::vector<std::pair<TimeInterval_t, std::vector<WaveformWithBaseline const*>>>
    intervals;
  
  // track which waveform should be considered next for each channel
  std::vector<std::vector<WaveformWithBaseline const*>::const_iterator> nextWaveform;
  for (std::vector<WaveformWithBaseline const*> const& channelWaveforms
    : waveformsPerChannel)
  {
    nextWaveform.push_back(channelWaveforms.begin());
  }
  
  for (auto [ start, stop ]
    : intervalsOver(overlapTimeline, requiredChannels.size()))
  {
    TimeInterval_t const interval(start, stop);
    
    auto log = mfLogTrace();
    log << "Interval " << interval << ":";
    
    std::vector<WaveformWithBaseline const*> contrib;
    
    for (auto const& [ channelIndex, channelWaveforms ]
      : util::enumerate(waveformsPerChannel))
    {
      std::vector<WaveformWithBaseline const*>::const_iterator& itNextWaveform
        = nextWaveform[channelIndex];
      
      std::vector<WaveformWithBaseline const*>::const_iterator const wend
        = channelWaveforms.end();
      
      // collect all waveforms overlapping `interval`;
      // assumes within each channel waveform list, waveforms are time-sorted
      log << "\n CH=" << (*itNextWaveform)->waveform().ChannelNumber() << ":";
      while (itNextWaveform != wend) {
        TimeInterval_t const wInterval = waveformInterval(**itNextWaveform);
        log << " -> vs. [#" << std::distance(channelWaveforms.begin(), itNextWaveform)
          << "] " << wInterval;
        if (wInterval.after(interval)) break; // we passed the interval
        
        if (wInterval.contains(interval)) {
          contrib.push_back(*itNextWaveform);
          log << " -> " << wInterval << ";";
          break;
        }
        
        ++itNextWaveform;
      } // while waveforms in channel
      
    } // for each channel
    
    log << " total: " << contrib.size() << " waveforms";
    assert(contrib.size() >= requiredChannels.size());
    assert(contrib.size() <= waveformsPerChannel.size());
    intervals.emplace_back(interval, std::move(contrib));
    
  } // while
  
  return intervals;
  
} // icarus::trigger::AdderSignalSimulation::findOverlappingTimeIntervals()


//------------------------------------------------------------------------------
bool icarus::trigger::AdderSignalSimulation::isMissingChannel
  (raw::Channel_t channel) const
{
  return std::binary_search
    (fMissingChannels.cbegin(), fMissingChannels.cend(), channel);
}


//------------------------------------------------------------------------------
auto icarus::trigger::AdderSignalSimulation::waveformInterval
  (raw::OpDetWaveform const& waveform) const -> TimeInterval_t
{
  electronics_time const start{ waveform.TimeStamp() };
  return { start, start + opticalTick() * waveform.size() };
} // icarus::trigger::AdderSignalSimulation::waveformInterval(WaveformWithBaseline)


auto icarus::trigger::AdderSignalSimulation::waveformInterval
  (WaveformWithBaseline const& waveform) const -> TimeInterval_t
  { return waveformInterval(waveform.waveform()); }


// -----------------------------------------------------------------------------
