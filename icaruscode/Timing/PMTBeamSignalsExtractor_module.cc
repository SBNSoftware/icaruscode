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
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"

#include "TTree.h"

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <utility>
#include <cstddef>

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

 * *  `ADCThreshold` (int): detection threshold to avoid cross-talk
 *    noise if one signal is missing from its waveform.
 * *  `DebugTrees` (bool): flag to produce plain ROOT trees for debugging.
 * *  `SaveWaveforms` (bool): flag to save full waveforms in the debug trees.
 *
 * Signal timing extraction
 * -------------------------
 *
 * The algorithm is quite unsophisticated and follows what is done for
 * the digitized trigger signal in `PMTWaveformTimeCorrectionExtractor`.
 * The reference signal is expected to be a sharp square wave in negative
 * polarity. The time is based on the front side of that wave:
 *
 *  * the absolute minimum of the waveform is found
 *  * an interval starting 20 ticks before that minimum is considered
 *  * the baseline level is defined as the value at the start of that interval
 *  * if the baseline-minimum difference is below a threshold, it is assumed to be noise
 *    and no time is returned
 *  * if not, the start time is set to the exact tick with an amplitude exceeding 20%
 *    of the minimum of the signal from the baseline
 *
 *
 * Output products
 * -------------------------
 *
 * This modules produces two `std::vector<icarus::timing::PMTBeamSignal>`
 * with 360 elements each, representing the relevent RWM or EW time for
 * the corresponding PMT channel. When the discrimination algorithm decides
 * the "peak" on a waveform to be noise, the corresponding entries are placed
 * into the vector with an invalid value (`isValid()` returns `false`,
 * but the identification data fields are correctly set).
 *
 * If the event is offbeam, these vectors are produced empty.
 *
 */

class icarus::timing::PMTBeamSignalsExtractor : public art::EDProducer
{
public:
  explicit PMTBeamSignalsExtractor(fhicl::ParameterSet const &pset);

  // process waveforms
  void extractBeamSignalTime(art::Event &e, art::InputTag label);
  template <typename T>
  static T Median(std::vector<T> data);
  template <typename T>
  static std::size_t getMaxBin(std::vector<T> const &vv, std::size_t startElement, std::size_t endElement);
  template <typename T>
  static std::size_t getMinBin(std::vector<T> const &vv, std::size_t startElement, std::size_t endElement);
  template <typename T>
  static std::size_t getStartSample(std::vector<T> const &vv, T thres);

  // unpack the V1730 special channels settings in a useful way
  static std::map<int, std::string> extractBoardBySpecialChannel(std::vector<fhicl::ParameterSet> const &setup);

  // associate times to PMT channels
  void associateBeamSignalsToChannels(art::InputTag label);

  // trigger-hardware timing correction
  double getTriggerCorrection(int channel);

  // quick mapping conversions
  std::string getDigitizerLabel(int channel);
  std::string getCrate(int channel);

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
  /// Save plain ROOT TTrees for debugging
  bool const fDebugTrees;
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
  /// Threshold for pulse selection
  short int const fADCThreshold;
  /// Special channel to board association in a map
  std::map<int, std::string> const fBoardBySpecialChannel;

  std::vector<icarus::timing::PMTWaveformTimeCorrection> fCorrections;
  std::map<std::string, int> fBoardEffFragmentID;
  double ftriggerTime;

  // output TTrees
  std::map<std::string, TTree *> fOutTree;
  // event info
  int m_run;
  int m_event;
  int m_timestamp;
  // special signal info
  int m_n_channels;
  unsigned int m_channel;
  double m_wfstart;
  std::size_t m_sample;
  double m_time;
  double m_time_abs;
  double m_utime_abs;
  std::vector<short> m_wf;

  // prepare pointers for data products
  using BeamSignalCollection = std::vector<icarus::timing::PMTBeamSignal>;
  using BeamSignalCollectionPtr = std::unique_ptr<BeamSignalCollection>;

  std::map<std::string, std::map<std::string, icarus::timing::PMTBeamSignal>> fBeamSignals;
  std::map<std::string, BeamSignalCollectionPtr> fSignalCollection;
};

// -----------------------------------------------------------------------------

icarus::timing::PMTBeamSignalsExtractor::PMTBeamSignalsExtractor(fhicl::ParameterSet const &pset)
    : EDProducer{pset},
      fChannelMap(*(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{})),
      fDebugTrees(pset.get<bool>("DebugTrees")),
      fSaveWaveforms(pset.get<bool>("SaveWaveforms")),
      fTriggerLabel(pset.get<art::InputTag>("TriggerLabel")),
      fRWMlabel(pset.get<art::InputTag>("RWMlabel")),
      fEWlabel(pset.get<art::InputTag>("EWlabel")),
      fTriggerCorrectionLabel(pset.get<art::InputTag>("TriggerCorrectionLabel")),
      fADCThreshold(pset.get<short int>("ADCThreshold")),
      fBoardBySpecialChannel(extractBoardBySpecialChannel(pset.get<std::vector<fhicl::ParameterSet>>("BoardSetup")))
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
  // prepare outupt TTrees if requested
  if (!fDebugTrees)
    return;

  art::ServiceHandle<art::TFileService> tfs;
  std::vector<std::string> labels = {fRWMlabel.instance(), fEWlabel.instance()};
  for (auto l : labels)
  {
    std::string name = l + "tree";
    std::string desc = l + " info";
    TTree *tree = tfs->make<TTree>(name.c_str(), desc.c_str());
    tree->Branch("run", &m_run);
    tree->Branch("event", &m_event);
    tree->Branch("timestamp", &m_timestamp);
    tree->Branch("n_channels", &m_n_channels);
    tree->Branch("channel", &m_channel);
    tree->Branch("wfstart", &m_wfstart);
    tree->Branch("sample", &m_sample);
    tree->Branch("utime_abs", &m_utime_abs);
    tree->Branch("time_abs", &m_time_abs);
    tree->Branch("time", &m_time);

    // add std::vector with full waveforms
    // this can quickly make the TTrees quite heavy
    if (fSaveWaveforms)
    {
      tree->Branch("wf", &m_wf);
    }
    fOutTree[l] = tree;
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

  // get the trigger-hardware corrections that are applied on all signal waveforms
  // if EW or RWM is to be compared to the PMT signals, it must be applied on them as well
  // it also takes care of board-to-board offsets (see SBN-doc-34631, slide 5)
  // Note: this is a vector of 360 elements, one correction for each signal channel
  fCorrections = e.getProduct<std::vector<icarus::timing::PMTWaveformTimeCorrection>>(fTriggerCorrectionLabel);
  int ntrig = fCorrections.size();

  if (ntrig < 1)
    mf::LogError("PMTBeamSignalsExtractor") << "Not found PMTWaveformTimeCorrections with label '"
                                            << fTriggerCorrectionLabel.instance() << "'";
  else if (ntrig < 360)
    mf::LogError("PMTBeamSignalsExtractor") << "Missing " << 360 - ntrig << " PMTWaveformTimeCorrections with label '"
                                            << fTriggerCorrectionLabel.instance() << "'";

  // now the main course: getting EW and RWM waveforms
  // information is stored by PMT crate name
  extractBeamSignalTime(e, fRWMlabel);
  extractBeamSignalTime(e, fEWlabel);

  // associating the proper RWM and EW time to each PMT channel
  // collections are vectors of 360 elements (one value for each channel)
  associateBeamSignalsToChannels(fRWMlabel);
  associateBeamSignalsToChannels(fEWlabel);

  // place data products in the stream
  // fix the cable swap for part of Run 2 right here!!
  // see SBN-doc-34631 for details
  if (gateType == sbn::triggerSource::BNB && m_run > 9704 && m_run < 11443)
  {

    e.put(std::move(fSignalCollection[fRWMlabel.instance()]), "EW");
    e.put(std::move(fSignalCollection[fEWlabel.instance()]), "RWM");
  }
  else
  { // STANDARD BEHAVIOR

    e.put(std::move(fSignalCollection[fRWMlabel.instance()]), "RWM");
    e.put(std::move(fSignalCollection[fEWlabel.instance()]), "EW");
  }
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::extractBeamSignalTime(art::Event &e, art::InputTag label)
{

  std::string l = label.instance();
  auto const &waveforms = e.getProduct<std::vector<raw::OpDetWaveform>>(label);
  m_n_channels = waveforms.size();

  if (m_n_channels < 1)
    mf::LogError("PMTBeamSignalsExtractor") << "Not found raw::OpDetWaveform with label '" << l << "'";
  else if (m_n_channels < 8)
    mf::LogError("PMTBeamSignalsExtractor") << "Missing " << 8 - m_n_channels << " raw::OpDetWaveform with label '" << l << "'";

  m_wf.clear();

  // get the start sample of the signals, one instance per PMT crate
  // use threshold to skip spikes from electric crosstalk see SBN-doc-34928, slides 4-5.
  // if no signal is found, set both time to "NoTime" so that it triggers the isValid() call of PMTBeamSignal

  for (auto const &wave : waveforms)
  {

    detinfo::timescales::electronics_time tstart = util::quantities::points::microsecond{wave.TimeStamp()};

    // if nothing is found, first sample is returned (0)
    m_sample = getStartSample(wave.Waveform(), fADCThreshold);

    m_channel = wave.ChannelNumber();
    m_wfstart = tstart.value();
    m_utime_abs = (m_sample != icarus::timing::NoSample) ? tstart.value() + 0.002 * m_sample : icarus::timing::NoTime;
    m_time_abs = (m_sample != icarus::timing::NoSample) ? m_utime_abs + getTriggerCorrection(m_channel) : icarus::timing::NoTime;
    m_time = (m_sample != icarus::timing::NoSample) ? m_time_abs - ftriggerTime : icarus::timing::NoTime;

    std::string crate = getCrate(m_channel);
    icarus::timing::PMTBeamSignal beamTime{m_channel, getDigitizerLabel(m_channel), crate, m_sample, m_time_abs, m_time};
    fBeamSignals[l].insert(std::make_pair(crate, beamTime));

    if (fSaveWaveforms)
      m_wf = wave.Waveform();
    if (fDebugTrees)
      fOutTree[l]->Fill();
  }
}

// ---------------------------------------------------------------------------

template <typename T>
T icarus::timing::PMTBeamSignalsExtractor::Median(std::vector<T> data)
{

  std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
  return data[data.size() / 2];
}

// -----------------------------------------------------------------------------

template <typename T>
std::size_t icarus::timing::PMTBeamSignalsExtractor::getMinBin(
    std::vector<T> const &vv, std::size_t startElement, std::size_t endElement)
{

  auto minel =
      std::min_element(vv.begin() + startElement, vv.begin() + endElement);
  std::size_t minsample = std::distance(vv.begin() + startElement, minel);

  return minsample;
}

// -----------------------------------------------------------------------------

template <typename T>
std::size_t icarus::timing::PMTBeamSignalsExtractor::getMaxBin(
    std::vector<T> const &vv, std::size_t startElement, std::size_t endElement)
{

  auto maxel =
      std::max_element(vv.begin() + startElement, vv.begin() + endElement);

  std::size_t maxsample = std::distance(vv.begin() + startElement, maxel);

  return maxsample;
}

// -----------------------------------------------------------------------------

template <typename T>
std::size_t icarus::timing::PMTBeamSignalsExtractor::getStartSample(std::vector<T> const &vv, T thres)
{

  // We are thinking in inverted polarity
  std::size_t minbin = getMinBin(vv, 0, vv.size());

  // Search only a cropped region of the waveform backward from the min
  std::size_t maxbin = minbin - 20;

  // Now we crawl betweem maxbin and minbin and we stop when:
  // bin value  > (maxbin value - bin value )*0.2
  std::size_t startbin = maxbin;
  auto delta = vv[maxbin] - vv[minbin];

  if (delta < thres)                 // just noise
    return icarus::timing::NoSample; // return no sample

  for (std::size_t bin = maxbin; bin < minbin; bin++)
  {
    auto val = vv[maxbin] - vv[bin];
    if (val >= 0.2 * delta)
    { // 20%
      startbin = bin - 1;
      break;
    }
  }

  if (startbin < maxbin)
  {
    startbin = maxbin;
  }

  return startbin;
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

std::string icarus::timing::PMTBeamSignalsExtractor::getDigitizerLabel(int channel)
{

  // get the board name, convert to digitizer_label
  std::string board = fBoardBySpecialChannel.at(channel);

  std::string head = "icaruspmt";
  std::string dash = "-";
  std::string letter = (board.substr(board.size() - 2, board.size()) == "02") ? "B" : "C";

  board.erase(board.find(head), head.size());
  std::transform(board.begin(), board.end(), board.begin(), ::toupper);
  board.insert(2, dash);
  board.insert(6, dash);

  return board.substr(0, board.size() - 2) + letter;
}

// -----------------------------------------------------------------------------

std::string icarus::timing::PMTBeamSignalsExtractor::getCrate(int channel)
{

  std::string digitizer_label = getDigitizerLabel(channel);
  return digitizer_label.substr(0, digitizer_label.size() - 2);
}

// -----------------------------------------------------------------------------

double icarus::timing::PMTBeamSignalsExtractor::getTriggerCorrection(int channel)
{

  std::string digitizer_label = getDigitizerLabel(channel);

  // trigger-hardware corrections are shared by all channels on the same board
  // we can pick the first channel on the desired board
  int fragID = fBoardEffFragmentID[digitizer_label];
  auto pmtinfo = fChannelMap.getPMTchannelInfo(fragID)[0]; // pick first ch on board
  int pmtch = pmtinfo.channelID;

  // trigger-hardware correction are in order
  // index of vector is pmtch
  return fCorrections.at(pmtch).startTime;
}

// -----------------------------------------------------------------------------

void icarus::timing::PMTBeamSignalsExtractor::associateBeamSignalsToChannels(art::InputTag label)
{

  std::string l = label.instance();

  // loop through the signals which are one per PMT crate
  // for each crate, find the corresponding digitizers
  for (auto signal : fBeamSignals[l])
  {

    // build the PMT digitizer labels that live in this crate
    // then convert it into fragment id
    std::vector<std::string> letters = {"-A", "-B", "-C"};
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
        // make sure there is enough room in the collection (channel -> vector index)
        if (channel >= fSignalCollection[l]->size())
          fSignalCollection[l]->resize(channel + 1);

        fSignalCollection[l]->at(channel) = signal.second;

      } // for each channel
    } // for each board
  } // for each crate
}

DEFINE_ART_MODULE(icarus::timing::PMTBeamSignalsExtractor)
