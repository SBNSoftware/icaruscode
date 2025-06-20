////////////////////////////////////////////////////////////////////////
/**
 * @file icaruscode/Timing/PMTBeamSignalsExtractor_module.cc
 * @brief `icarus::timing::PMTBeamSignalsExtractor` producer module.
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date  Sun Feb 11 11:37:14 2024
 **/
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/IcarusObj/PMTWaveformTimeCorrection.h"
#include "icaruscode/IcarusObj/PMTBeamSignal.h"
#include "icaruscode/Timing/Tools/PulseStartExtractor.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/quantities/spacetime.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"

#include "TTree.h"
#include "TMath.h"

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <utility>
#include <cstddef>
#include <array>

namespace icarus::timing
{
  class PMTBeamSignalsExtractor;
}

/**
 * @brief Extracts RWM and EW times from special waveforms.
 *
 * This module extracts the RWM and EW from the waveforms recorded in the
 * spare PMT channels and associates these times with all the channels
 * in the same PMT readout crate.
 *
 * Input parameters
 * ------------------------
 *
 * This modules requires the following products:
 *
 * * `TriggerLabel` (input tag): tag for the `sbn::ExtraTriggerInfo`
 *    data product. This is required to check if it's a beam event or not.
 * * `RWMlabel` (input tag): tag for the `std::vector<raw::OpDetWaveform>`
 *    product containing the RWM special waveforms.
 * * `EWlabel` (input tag): tag for the `std::vector<raw::OpDetWaveform>`
 *    product containign the EW special waveforms.
 * *  `TriggerCorrectionLabel` (input tag): tag for the
 *    std::vector<icarus::timing::PMTWaveformTimeCorrection> product that
 *    stores the channel-by-channel waveform timing corrections
 * *  `BoardSetup` (fhicl parameter set): description of the current
 *    V1730 setup from `CAEN_V1730_setup_icarus.fcl` mapping the special signals.
 *    It is meant to be the same configuration as used by the PMT decoding
 *    (see `BoardSetup` configuration parameter in `icarus::DaqDecoderICARUSPMT`).
 * *  `TimeExtractionMethod` (string): pulse start time extraction method. 
 *     @see icarus::timing::PulseStartExtractor for details and options.
 * *  `ADCThreshold` (int): detection threshold to avoid cross-talk
 *    noise if one signal is missing from its waveform.
 * *  `ApplyConsensusFiltering` (bool): if true, filter bad or missing times using 
 *    the consensus among other crates in the same cryostat.
 * *  `DebugTree` (bool): flag to produce plain ROOT tree for debugging.
 * *  `SaveWaveforms` (bool): flag to save full waveforms in the debug tree.
 *
 * Signal timing extraction
 * -------------------------
 *
 * The RWM/EW signals are expected to be a sharp square wave in negative
 * polarity. The time of the correction is based on the front side of that wave.
 * The extraction is perfomed by `icarus::timing::PulseStartExtractor` by
 * specifying one of the available extraction methods (constant-fraction discrimation,
 * logitstic function fit) and an ADC threshold for the signal identification.
 * @see `icarus::timing::PulseStartExtractor` for details.
 *
 * The extracted times are then corrected using the standard waveform time corrections
 * based on the digitized trigger signals (`icarus::timing::PMTWaveformTimeCorrection`).
 * After this step, times should match between crates of the same module apart from
 * small differences in the fit results. 
 
 * Unfortunately, this is not always the case: a crate might be completely missing a signal 
 * or seeing it distorted. While these issues are eventually corrected in the hardware, 
 * the redundancy due to the other crates is exploited to recover a value for the affected crate.
 * All the valid (non-missing) signals are used to estimate the mean, median and standard
 * deviation of the extracted times in each cryostat. Large outliers from the median
 * (more than 2 &sigma away) are then replaced with the median itself.  
 *
 * Output products
 * -------------------------
 *
 * This modules produces two `std::vector<icarus::timing::PMTBeamSignal>`
 * with 360 elements each, representing the relevant RWM or EW time for
 * the corresponding PMT channel. When the discrimination algorithm decides
 * the "peak" on a waveform to be noise, the corresponding entries are placed
 * into the vector with an invalid value (`isValid()` returns `false`,
 * but the identification data fields are correctly set).
 *
 * If the event is offbeam, these vectors are produced empty.
 *
 * 
 * Debugging tree
 * ---------------
 *
 * The debugging tree, enabled by `DebugTree`, contains both Early Warning (EW) signals,
 * named after the instance name of `EWlabel`, and Resistive Wall Monitor (RWM) signals,
 * named after the instance name of `RWMlabel`. Each entry in the tree represents what a
 * PMT crate sees in an event (thus, typically for ICARUS there are 8 entries per event).
 * The content of an entry includes:
 *   * `run`, `event`: identifier of the run/event.
 *   * `timestamp`: timestamp of the event (in UTC), truncated at the second.
 *   * `n_crates`: number of crates with special signals in the event.
 *   * `crate`: string identifier for the PMT crate (eg, EE-BOT).
 *   * `rmw_channel`: special RWM channel number assigned to this signal.
 *   * `rwm_wfstart`: (uncorrected) RWM waveform time start.
 *   * `rwm_sample`: the sample number in the RWM waveform at which the signal was found.
 *   * `rwm_utime_abs`: uncorrected RWM signal time [&micro;s].
 *   * `rwm_itime_abs`: intermediate (corrected) RWM signal time [&micro;s].
 *   * `rwm_time_abs`:  corrected and adjusted RWM signal time [&micro;s].
 *   * `rwm_time`: corrected and adjusted RWM signal time (relative to the trigger time) [&micro;s].
 *   * `rwm_nsize` (only if `SaveWaveforms` is set): number of samples in RWM waveform.
 *   * `rwm_wf` (only if `SaveWaveforms` is set): full RWM waveform.
 *   * `ew_channel`: special EW channel number assigned to this signal.
 *   * `ew_wfstart`: (uncorrected) EW waveform time start.
 *   * `ew_sample`: the sample number in the EW waveform at which the signal was found.
 *   * `ew_utime_abs`: uncorrected EW signal time [&micro;s].
 *   * `ew_itime_abs`: intermediate (corrected) EW signal time [&micro;s].
 *   * `ew_time_abs`:  corrected and adjusted EW signal time [&micro;s].
 *   * `ew_time`: corrected and adjusted EW signal time (relative to the trigger time) [&micro;s].
 *   * `ew_nsize` (only if `SaveWaveforms` is set): number of samples in EW waveform.
 *   * `ew_wf` (only if `SaveWaveforms` is set): full EW waveform.
 */

class icarus::timing::PMTBeamSignalsExtractor : public art::EDProducer
{
public:
  explicit PMTBeamSignalsExtractor(fhicl::ParameterSet const &pset);

  // process waveforms
  void extractBeamSignalTime(art::Event &e, art::InputTag const &label);

  // filter times based on consensus among crates
  void filterBeamSignalsByCryostat(std::string const &l, std::vector<double> const &times, std::string const &c);
  double stdDev(std::vector<double> data, double mean);

  // associate times to PMT channels
  void associateBeamSignalsToChannels(art::InputTag const &label);

  // trigger-hardware timing correction
  double getTriggerCorrection(unsigned int channel) const;

  // quick mapping conversions
  std::string getDigitizerLabel(unsigned int channel) const;
  std::string getCrate(unsigned int channel) const;

  // unpack the V1730 special channels settings in a useful way
  static std::map<int, std::string> extractBoardBySpecialChannel(std::vector<fhicl::ParameterSet> const &setup);

  // fill debug trees
  void fillDebugTree();

  // Plugins should not be copied or assigned.
  PMTBeamSignalsExtractor(PMTBeamSignalsExtractor const &) = delete;
  PMTBeamSignalsExtractor(PMTBeamSignalsExtractor &&) = delete;
  PMTBeamSignalsExtractor &operator=(PMTBeamSignalsExtractor const &) = delete;
  PMTBeamSignalsExtractor &operator=(PMTBeamSignalsExtractor &&) = delete;

  void beginJob() override;
  void beginRun(art::Run &run) override;
  void produce(art::Event &e) override;

private:
  /// Channel mappping
  icarusDB::IICARUSChannelMap const &fChannelMap;
  /// Save plain ROOT TTree for debugging
  bool const fDebugTree;
  /// Save raw waveforms in debug TTrees
  bool const fSaveWaveforms;
  /// Trigger instance label
  art::InputTag const fTriggerLabel;
  /// RWM waveform instance label
  art::InputTag const fRWMlabel;
  /// EW waveform instance label
  art::InputTag const fEWlabel;
  /// Trigger-hardware correction instance
  art::InputTag const fTriggerCorrectionLabel;
  /// Manager for pulse start extraction
  icarus::timing::PulseStartExtractor const fPulseStartExtractor;
  /// Special channel to board association in a map
  std::map<int, std::string> const fBoardBySpecialChannel;
  /// Apply filtering for bad/missing signals
  bool const fApplyConsensusFiltering;

  /// PMT sample duration [&micro;s]
  static constexpr double fPMTsamplingTick = 0.002;
  /// Number of PMT channels
  static constexpr std::size_t fNPMTChannels = 360;

  std::vector<icarus::timing::PMTWaveformTimeCorrection> fCorrections;
  std::map<std::string, int> fBoardEffFragmentID;
  double ftriggerTime;

  // output TTrees
  TTree* fOutTree;
  // event info
  int m_run;
  int m_event;
  int m_timestamp;
  int m_n_crates;
  std::string m_crate;
  // rwm signal info
  unsigned int m_rwm_channel;
  double m_rwm_wfstart;
  double m_rwm_sample;
  double m_rwm_time;
  double m_rwm_time_abs;
  double m_rwm_utime_abs;
  double m_rwm_itime_abs;
  std::size_t m_rwm_nsize;
  std::vector<short> m_rwm_wf;
  // ew signal info
  unsigned int m_ew_channel;
  double m_ew_wfstart;
  double m_ew_sample;
  double m_ew_time;
  double m_ew_time_abs;
  double m_ew_utime_abs;
  double m_ew_itime_abs;
  std::size_t m_ew_nsize;
  std::vector<short> m_ew_wf;

  /// Struct for storing data
  struct SignalData {
    unsigned int channel;
    std::string diglabel;
    std::string crate;
    double wfstart;
    double sample;
    double utime_abs;
    double itime_abs;
    double time_abs;
    double time;
    std::vector<short> wf;
  };

  // prepare pointers for data products
  using BeamSignalCollection = std::vector<icarus::timing::PMTBeamSignal>;
  using BeamSignalCollectionPtr = std::unique_ptr<BeamSignalCollection>;

  std::map<std::string, std::map<std::string, SignalData>> fBeamSignals;
  std::map<std::string, BeamSignalCollectionPtr> fSignalCollection;
};

// -----------------------------------------------------------------------------

icarus::timing::PMTBeamSignalsExtractor::PMTBeamSignalsExtractor(fhicl::ParameterSet const &pset)
    : EDProducer{pset},
      fChannelMap(*(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{})),
      fDebugTree(pset.get<bool>("DebugTree")),
      fSaveWaveforms(pset.get<bool>("SaveWaveforms")),
      fTriggerLabel(pset.get<art::InputTag>("TriggerLabel")),
      fRWMlabel(pset.get<art::InputTag>("RWMlabel")),
      fEWlabel(pset.get<art::InputTag>("EWlabel")),
      fTriggerCorrectionLabel(pset.get<art::InputTag>("TriggerCorrectionLabel")),
      fPulseStartExtractor{ icarus::timing::stringToExtractionMethod.at(pset.get<std::string>("TimeExtractionMethod")), 
                            pset.get<double>("ADCThreshold")},
      fBoardBySpecialChannel(extractBoardBySpecialChannel(pset.get<std::vector<fhicl::ParameterSet>>("BoardSetup"))),
      fApplyConsensusFiltering(pset.get<bool>("ApplyConsensusFiltering"))
{
  // Call appropriate consumes<>() functions here.
  consumes<std::vector<raw::OpDetWaveform>>(fEWlabel);
  consumes<std::vector<raw::OpDetWaveform>>(fRWMlabel);

  // Call appropriate produces<>() functions here.
  produces<std::vector<PMTBeamSignal>>("RWM");
  produces<std::vector<PMTBeamSignal>>("EW");
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::beginJob()
{
  // prepare outupt TTree if requested
  if (!fDebugTree)
    return;

  art::ServiceHandle<art::TFileService> tfs;
  std::string name = fRWMlabel.instance() + "_" + fEWlabel.instance() + "_tree";
  std::string desc = fRWMlabel.instance() + " " + fEWlabel.instance() + " info";
  fOutTree = tfs->make<TTree>(name.c_str(), desc.c_str());
  fOutTree->Branch("run", &m_run);
  fOutTree->Branch("event", &m_event);
  fOutTree->Branch("timestamp", &m_timestamp);
  fOutTree->Branch("n_crates", &m_n_crates);
  fOutTree->Branch("crate", &m_crate);
  fOutTree->Branch("rwm_channel", &m_rwm_channel);
  fOutTree->Branch("rwm_wfstart", &m_rwm_wfstart);
  fOutTree->Branch("rwm_sample", &m_rwm_sample);
  fOutTree->Branch("rwm_utime_abs", &m_rwm_utime_abs);
  fOutTree->Branch("rwm_itime_abs", &m_rwm_itime_abs);
  fOutTree->Branch("rwm_time_abs", &m_rwm_time_abs);
  fOutTree->Branch("rwm_time", &m_rwm_time);
  fOutTree->Branch("ew_channel", &m_ew_channel);
  fOutTree->Branch("ew_wfstart", &m_ew_wfstart);
  fOutTree->Branch("ew_sample", &m_ew_sample);
  fOutTree->Branch("ew_utime_abs", &m_ew_utime_abs);
  fOutTree->Branch("ew_itime_abs", &m_ew_itime_abs);
  fOutTree->Branch("ew_time_abs", &m_ew_time_abs);
  fOutTree->Branch("ew_time", &m_ew_time);

  // add std::vector with full waveforms
  // this can quickly make the TTrees quite heavy
  if (fSaveWaveforms)
  {
    fOutTree->Branch("rwm_nsize", &m_rwm_nsize);
    fOutTree->Branch("rwm_wf", &m_rwm_wf);
    fOutTree->Branch("ew_nsize", &m_ew_nsize);
    fOutTree->Branch("ew_wf", &m_ew_wf);
  }
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::beginRun(art::Run &run)
{
  // pre-save the association between digitizer_label and effective fragment ID
  // needs to be done at the begin of each run in case mapping changed
  fBoardEffFragmentID.clear();

  for (unsigned int fragid = 0; fragid < fChannelMap.nPMTfragmentIDs(); fragid++)
  {
    auto const &pmtinfo = fChannelMap.getPMTchannelInfo(fragid)[0]; // pick first pmt on board
    fBoardEffFragmentID[pmtinfo.digitizerLabel] = fragid;
  }
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::produce(art::Event &e)
{

  // initialize the data products
  fBeamSignals.clear();
  fSignalCollection[fRWMlabel.instance()] = std::make_unique<BeamSignalCollection>();
  fSignalCollection[fEWlabel.instance()] = std::make_unique<BeamSignalCollection>();

  // extract meta event information
  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh();                                                         // precision to the second
  ftriggerTime = e.getProduct<std::vector<raw::Trigger>>(fTriggerLabel).at(0).TriggerTime(); // us

  // this module should run on beam-only events
  // check the current beam gate
  auto const &triggerInfo = e.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel);
  sbn::triggerSource const gateType = triggerInfo.sourceType;

  switch (gateType)
  {
  case sbn::triggerSource::BNB:
    break;
  case sbn::triggerSource::NuMI:
    break;
  default:
    mf::LogTrace("PMTBeamSignalsExtractor") << "Skipping offbeam gate '" << name(gateType) << "'";
    e.put(std::move(fSignalCollection[fRWMlabel.instance()]), "RWM");
    e.put(std::move(fSignalCollection[fEWlabel.instance()]), "EW");
    return;
  }

  // reside the collections to match the expected 360 PMT channels
  // this is to avoid dynamic resizing later on
  for (auto &collection : fSignalCollection)
    collection.second->resize(fNPMTChannels);

  // get the trigger-hardware corrections that are applied on all signal waveforms
  // if EW or RWM is to be compared to the PMT signals, it must be applied on them as well
  // it also takes care of board-to-board offsets (see SBN-doc-34631, slide 5)
  // Note: this is a vector of 360 elements, one correction for each signal channel
  fCorrections = e.getProduct<std::vector<icarus::timing::PMTWaveformTimeCorrection>>(fTriggerCorrectionLabel);
  std::size_t ntrig = fCorrections.size();

  if (ntrig < 1)
    mf::LogError("PMTBeamSignalsExtractor") << "Not found PMTWaveformTimeCorrections with label '"
                                            << fTriggerCorrectionLabel.instance() << "'";
  else if (ntrig < fNPMTChannels)
    mf::LogError("PMTBeamSignalsExtractor") << "Missing " << fNPMTChannels - ntrig << " PMTWaveformTimeCorrections with label '"
                                            << fTriggerCorrectionLabel.instance() << "'";

  // now the main course: getting EW and RWM waveforms
  // information is stored by PMT crate name
  extractBeamSignalTime(e, fRWMlabel);
  extractBeamSignalTime(e, fEWlabel);

  // fix the cable swap for part of Run 2 right here!!
  // see SBN-doc-34631 for details
  if (gateType == sbn::triggerSource::BNB && m_run > 9704 && m_run < 11443)
    std::swap(fBeamSignals[fRWMlabel.instance()], fBeamSignals[fEWlabel.instance()]);

  if (fDebugTree) fillDebugTree();

  // associating the proper RWM and EW time to each PMT channel
  // collections are vectors of 360 elements (one value for each channel)
  associateBeamSignalsToChannels(fRWMlabel);
  associateBeamSignalsToChannels(fEWlabel);

  // place data products in the stream
  e.put(std::move(fSignalCollection[fRWMlabel.instance()]), "RWM");
  e.put(std::move(fSignalCollection[fEWlabel.instance()]), "EW");
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::extractBeamSignalTime(art::Event &e, art::InputTag const &label)
{

  std::string const &l = label.instance();
  auto const &waveforms = e.getProduct<std::vector<raw::OpDetWaveform>>(label);
  int n_channels = waveforms.size();

  if (n_channels < 1)
    mf::LogError("PMTBeamSignalsExtractor") << "Not found raw::OpDetWaveform with label '" << label.encode() << "'";
  else if (n_channels < 8)
    mf::LogError("PMTBeamSignalsExtractor") << "Missing " << 8 - n_channels << " raw::OpDetWaveform with label '"
                                            << label.encode() << "'";

  m_rwm_wf.clear();
  m_ew_wf.clear();

  // prepare to separe out east/west signals
  std::vector<double> east;
  std::vector<double> west;

  // get the start sample of the signals, one instance per PMT crate
  // threshold is used to skip spikes from electric crosstalk see SBN-doc-34928, slides 4-5.
  // if no signal is found, set time to icarus::timing::NoTime
  for (auto const &wave : waveforms)
  {

    detinfo::timescales::electronics_time tstart = util::quantities::points::microsecond{wave.TimeStamp()};

    unsigned int channel = wave.ChannelNumber();
    double wfstart = tstart.value();
    std::string crate = getCrate(channel);

    // if nothing is found, first sample is returned (0)
    double sample = fPulseStartExtractor.extractStart(wave.Waveform());
    if (sample < 1){
      sample = icarus::timing::NoSample;
      mf::LogTrace("PMTBeamSignalsExtractor") << "No " << label.encode() << " signal found in channel " << channel
                                              << " (crate " << crate << ")";
    }

    double utime_abs = (sample != icarus::timing::NoSample) ? tstart.value() + fPMTsamplingTick * sample : icarus::timing::NoTime;
    double itime_abs = (sample != icarus::timing::NoSample) ? utime_abs + getTriggerCorrection(channel) : icarus::timing::NoTime;
    double itime = (sample != icarus::timing::NoSample) ? itime_abs - ftriggerTime : icarus::timing::NoTime;

    // store signal information, final time might be updated later
    SignalData beamTime{channel, getDigitizerLabel(channel), crate, wfstart, sample, utime_abs, itime_abs, itime_abs, itime,
                        fSaveWaveforms ? wave.Waveform() : std::vector<short>{}};
    fBeamSignals[l].emplace(crate, std::move(beamTime));

    // skip invalid times (signal completely missing)
    if( sample == icarus::timing::NoSample ) continue; 

    // save valid times in east/west separately to estimate consensus
    // this will be used to computed median and RMS
    if(crate[0] == 'E') east.push_back(itime_abs);
    else if( crate[0] == 'W') west.push_back(itime_abs);
    
  }

  // filter out and correct missing signals or large outliers
  // using the consensus from the majority of crates in the same cryostat
  if( fApplyConsensusFiltering ){
    filterBeamSignalsByCryostat(l, east, "E");
    filterBeamSignalsByCryostat(l, west, "W");
  }
}

// -----------------------------------------------------------------------------

double icarus::timing::PMTBeamSignalsExtractor::stdDev(std::vector<double> data, double mean)
{
  double sum = 0;
  for(auto d : data) sum += (d-mean)*(d-mean);
  return TMath::Sqrt(sum/data.size());
}

void icarus::timing::PMTBeamSignalsExtractor::filterBeamSignalsByCryostat(
  std::string const &l, std::vector<double> const &times, std::string const &c)
{
  
  // if there is only one valid time in this cryostat,
  // do no attempt any filtering/adjustments
  if (times.size() < 2) return;

  // compute stddev and median among valid signals (consensus)
  double median = TMath::Median(times.size(), times.data());
  double mean = TMath::Mean(times.size(), times.data());
  double stddev = stdDev(times, mean);

  for (auto const& [crate, signal] : fBeamSignals[l]){
    
    // skip the other cryostat
    if( crate[0] != c) continue; 
    double t = signal.itime_abs;
    
    // if the current signal is more than 2 standard deviations away, substitute it with the median
    // this affects both missing signals (whose value is currently icarus::timing::NoTime)
    // and distorted signal for which the time extraction is skewed
    if ( std::fabs(t-median) > 2*stddev )
    { 
      fBeamSignals[l][crate].time_abs = median; 
      fBeamSignals[l][crate].time = median - ftriggerTime;
    }

  } 

}

// -----------------------------------------------------------------------------

std::map<int, std::string> icarus::timing::PMTBeamSignalsExtractor::extractBoardBySpecialChannel(
    std::vector<fhicl::ParameterSet> const &setup)
{

  // map from special PMT channel to corresponding board
  std::map<int, std::string> boardBySpecialChannel;

  // unpack the V1730 special channels settings in a useful way
  for (fhicl::ParameterSet const &s : setup)
  {
    auto const &innerSet = s.get<std::vector<fhicl::ParameterSet>>("SpecialChannels");
    boardBySpecialChannel[innerSet[0].get<int>("Channel")] = s.get<std::string>("Name");
  }

  return boardBySpecialChannel;
}

// -----------------------------------------------------------------------------

std::string icarus::timing::PMTBeamSignalsExtractor::getDigitizerLabel(unsigned int channel) const
{

  // get the board name, convert to digitizer_label
  std::string board = fBoardBySpecialChannel.at(channel); // eg. icaruspmtewtop02

  std::string head = "icaruspmt";
  char letter = board.back() - '1' + 'A'; // converts 01,02,03 to A,B,C

  board.erase(board.find(head), head.size());                           // eg. ewtop02
  std::transform(board.begin(), board.end(), board.begin(), ::toupper); // eg. EWTOP02
  board.insert(2, 1, '-');                                              // insert dash at position 2, e.g: EW-TOP02
  board.insert(6, 1, '-');                                              // insert dash at position 6, e.g: EW-TOP-02

  return board.substr(0, board.size() - 2) + letter; // e.g: EW-TOP-B
}

// -----------------------------------------------------------------------------

std::string icarus::timing::PMTBeamSignalsExtractor::getCrate(unsigned int channel) const
{

  std::string digitizer_label = getDigitizerLabel(channel);
  return digitizer_label.substr(0, digitizer_label.size() - 2);
}

// -----------------------------------------------------------------------------

double icarus::timing::PMTBeamSignalsExtractor::getTriggerCorrection(unsigned int channel) const
{

  std::string digitizer_label = getDigitizerLabel(channel);

  // trigger-hardware corrections are shared by all channels on the same board
  // we can pick the first channel on the desired board
  int fragID = fBoardEffFragmentID.at(digitizer_label);
  auto pmtinfo = fChannelMap.getPMTchannelInfo(fragID)[0]; // pick first ch on board
  int pmtch = pmtinfo.channelID;

  // trigger-hardware correction are in order
  // index of vector is pmtch
  return fCorrections.at(pmtch).startTime;
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::associateBeamSignalsToChannels(art::InputTag const &label)
{

  std::string const &l = label.instance();

  // loop through the signals which are one per PMT crate
  // for each crate, find the corresponding digitizers
  for (auto const &signal : fBeamSignals[l])
  {

    // convert SignalData into icarus::timing::PMTBeamSignal
    icarus::timing::PMTBeamSignal beamTime{ signal.second.channel, signal.second.diglabel, signal.first,
                                            signal.second.sample, signal.second.time_abs, signal.second.time};

    // build the PMT digitizer labels that live in this crate
    // then convert it into fragment id
    std::array const letters{"-A", "-B", "-C"};
    for (auto letter : letters)
    {

      std::string digitizer_label = signal.first + letter;
      int fragID = fBoardEffFragmentID[digitizer_label];

      // use the fragment id to access the PMT channels
      // mapping works via fragment id (it's very annoying..)
      // loop through the PMT channels and set their RWM/EW time
      for (auto pmtinfo : fChannelMap.getPMTchannelInfo(fragID))
      {

        std::size_t channel = pmtinfo.channelID;
        // collections are resized to fNPMTChannels
        // always room in the collection (channel -> vector index)
        fSignalCollection[l]->at(channel) = beamTime;

      } // for each channel
    } // for each board
  } // for each crate
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::fillDebugTree()
{
  
  auto RWMSignals = fBeamSignals[fRWMlabel.instance()];
  auto EWSignals = fBeamSignals[fEWlabel.instance()];

  std::unordered_set<std::string> crates;
  for (const auto& [crate, _] : RWMSignals ) crates.insert(crate);
  for (const auto& [crate, _] : EWSignals ) crates.insert(crate);

  m_n_crates = crates.size();
  
  for (const auto& crate : crates)
  {
    auto itRWM = RWMSignals.find(crate);
    auto itEW = EWSignals.find(crate);
    m_crate = crate;

    if (itRWM != RWMSignals.end()) 
    {
      m_rwm_channel   = itRWM->second.channel;
      m_rwm_wfstart   = itRWM->second.wfstart;
      m_rwm_sample    = itRWM->second.sample;
      m_rwm_utime_abs = itRWM->second.utime_abs; //uncorrected
      m_rwm_itime_abs = itRWM->second.itime_abs; //pre-filtering
      m_rwm_time_abs  = itRWM->second.time_abs;  //post-filtering
      m_rwm_time      = itRWM->second.time;
      if (fSaveWaveforms)
      {
        m_rwm_nsize = itRWM->second.wf.size();
        m_rwm_wf    = itRWM->second.wf;
      }
    } else {
      m_rwm_channel = 0;
      m_rwm_wfstart   = icarus::timing::NoTime;
      m_rwm_sample    = icarus::timing::NoTime;
      m_rwm_utime_abs = icarus::timing::NoTime; 
      m_rwm_itime_abs = icarus::timing::NoTime; 
      m_rwm_time_abs  = icarus::timing::NoTime;
      m_rwm_time      = icarus::timing::NoTime;
    }

    if (itEW != EWSignals.end()) 
    {
      m_ew_channel   = itEW->second.channel;
      m_ew_wfstart   = itEW->second.wfstart;
      m_ew_sample    = itEW->second.sample;
      m_ew_utime_abs = itEW->second.utime_abs; //uncorrected
      m_ew_itime_abs = itEW->second.itime_abs; //pre-filtering
      m_ew_time_abs  = itEW->second.time_abs;  //post-filtering
      m_ew_time      = itEW->second.time;
      if (fSaveWaveforms)
      {
        m_ew_nsize = itEW->second.wf.size();
        m_ew_wf    = itEW->second.wf;
      }
    } else {
      m_ew_channel = 0;
      m_ew_wfstart   = icarus::timing::NoTime;
      m_ew_sample    = icarus::timing::NoTime;
      m_ew_utime_abs = icarus::timing::NoTime; 
      m_ew_itime_abs = icarus::timing::NoTime; 
      m_ew_time_abs  = icarus::timing::NoTime;
      m_ew_time      = icarus::timing::NoTime;
    }
    
    fOutTree->Fill();
  }
}

DEFINE_ART_MODULE(icarus::timing::PMTBeamSignalsExtractor)

