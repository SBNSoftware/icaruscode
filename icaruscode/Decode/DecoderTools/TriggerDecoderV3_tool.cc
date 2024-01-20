////////////////////////////////////////////////
//   
//    File: TriggerDecoderV3_tool.cc
//       
//    Description: Extracting ICARUS trigger fragment information from new fragment necessary after trigger information  
//
//    Author: Jacob Zettlemoyer, FNAL
//
///////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // util::quantities::nanosecond
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerV3Fragment.hh"

#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "sbnobj/Common/Trigger/BeamBits.h"
#include "icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h" // dumpFragment()
#include "icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h"
#include "icaruscode/Decode/DecoderTools/details/KeyValuesData.h"
#include "icarusalg/Utilities/BinaryDumpUtils.h" // hexdump() DEBUG

#include <cstdlib>
#include <iostream>
#include <iomanip> // std::setw(), std::setfill()
#include <string_view>
#include <memory>
#include <array>


using namespace std::string_literals;

namespace daq 
{
  
  /**
   * @brief Tool decoding the trigger information from DAQ.
   * 
   * Produces:
   * * `std::vector<raw::ExternalTrigger>` containing the absolute trigger time
   *     stamp from the White Rabbit system, and a trigger count;
   *     it always includes a single entry (zero _might_ be supported).
   * * `std::vector<raw::Trigger>` containing:
   *     * `TriggerTime()`: the relative time of the trigger as reported in the
   *         White Rabbit timestamp, measured in the
   *         @ref DetectorClocksElectronicsTime "electronics time scale" (for
   *         ICARUS it will always be
   *         `detinfo::DetectorClocksData::TriggerTime()`).
   *     * `BeamGateTime()`: relative time of the announced arrival of the beam,
   *         also in @ref DetectorClocksElectronicsTime "electronics time scale".
   *     * `TriggerCount()`: the trigger count from the beginning of the run.
   *     * `TriggerBits()`: includes the beam(s) with an open gate when the
   *         trigger happened (currently only one beam gate per trigger);
   *         definitions are in `sbn::beamType` namespace.
   * 
   *     It always includes a single entry (zero _might_ be supported).
   * * `std::vector<sim::BeamGateInfo>` containing information on the "main"
   *     beam gate associated to each trigger (a way to say that if by any
   *     chance there are more than one bits set for the trigger, this gate
   *     will pick only one of them):
   *     * `Start()`: relative time of the announced arrival of the beam, in
   *         @ref DetectorClocksSimulationTime "simulation time scale".
   *     * `Width()`: duration of the gate, in nanoseconds; read from trigger
   *         configuration if specified (`TrigConfigLabel`), set to `0`
   *         otherwise.
   *     * `BeamType()`: the type of the beam gate being described (BNB, NuMI).
   * * `sbn::ExtraTriggerInfo`: the most complete collection of information,
   *     duplicating also some from the other data products.
   *     Note that differently from the usual, this is a _single object_, not
   *     a collection; also, this data product has no instance name.
   *     The information already available includes:
   *     * `sourceType`: the type of beam or trigger source, a value from
   *         `sbn::triggerSource` (equivalent to `raw::Trigger::TriggerBits()`,
   *         but in the form of an enumerator rather than a bit mask).
   *         Also called gate type.
   *     * `triggerType`: the type of trigger logic applied, a value from
   *         `sbn::triggerType`, e.g. minimum bias or majority.
   *     * `triggerTimestamp`: same as `raw::ExternalTrigger::GetTrigTime()`
   *         (nanoseconds from the Epoch, Coordinated Universal Time).
   *     * `beamGateTimestamp`: absolute time of the beam gate opening as
   *         reported by the trigger hardware, directly comparable with
   *         `triggerTimestamp` (same scale and unit).
   *     * `enableGateTimestamp`: absolute time of the enable gate opening as
   *         reported by the trigger hardware, directly comparable with
   *         `triggerTimestamp` (same scale and unit). This is the gate that
   *         enables the generation of trigger primitives and the off-spill
   *         readout of PMT. Its duration can be found in the trigger
   *         configuration data product (`icarus::TriggerConfiguration`).
   *     * `triggerID`: same as `raw::ExternalTrigger::GetTrigID()`. Should
   *         match the event number.
   *     * `gateID`: the count of gates since the beginning of the run, as
   *         reported by the trigger hardware.
   *     * `gateCountFromPreviousTrigger`: number of gates since the last
   *         trigger: specifically, if this is e.g. a majority trigger is from 
   *         off-beam BNB gate, this is the number of off-beam BNB gates
   *         considered for majority triggers (i.e. excluding minimum bias) from
   *         the last off-beam BNB majority trigger (minimum number is `1`,
   *         since that gate is included). The same logic applies for minimum
   *         bias triggers, which excludes from the count all minimum bias
   *         triggers (and therefore for minimum bias this count is always
   *         expected to be `1`).
   *     * `previousTriggerTimestamp`: absolute timestamp of the previous
   *         trigger from this source. For example, if this is a majority
   *         trigger from an off-beam BNB gate, this represents the previous
   *         off-beam BNB majority trigger.
   *     * `gateCount`: total number of gates of this type (triggered or not)
   *         from the beginning of the run (minimum `1` since this one is
   *         included). For example, if this trigger is from an off-beam BNB
   *         gate, this is the number of off-beam BNB gates from the beginning
   *         of the run. Note that this count, differently from the
   *         source-specific _trigger_ counts below, mixes the types of triggers
   *         (majority and minimum bias).
   *     * `triggerCount`: total number of triggers from gates of this type
   *         from the beginning of the run (minimum `1` since this one is
   *         included). For example, if this trigger is from an off-beam BNB
   *         majority trigger, this is the number of off-beam BNB majority
   *         triggers from the beginning of the run.
   *     * `anyTriggerCountFromPreviousTrigger`: number of triggers that
   *         occurred since the last trigger of this type (the one with
   *         timestamp `previousTriggerTimestamp`).
   *     * `anyGateCountFromAnyPreviousTrigger`: _not available_. It would hold
   *         how many gates have passed since the last trigger (reported by
   *         `anyPreviousTriggerTimestamp`), with minimum value `1`.
   *     * `anyPreviousTriggerTimestamp`: absolute timestamp of the previous
   *         trigger (from any source). For example, if this trigger is from an
   *         off-beam BNB gate, and the previous was from a NuMI gate, this
   *         represents that previous (NuMI) trigger.
   *     * `anyPreviousTriggerSourceType`: the type of gate of the previous
   *         trigger (the one reported by `anyPreviousTriggerTimestamp`).
   *     * `WRtimeToTriggerTime` (nanoseconds): correction added to the
   *         GPS/White Rabbit time to obtain the trigger timestamp (normally
   *         it's the conversion from TAI to NTP timestamps).
   *     * `triggerLocationBits`: whether the trigger came from the east or west
   *         cryostat (currently the triggers combine the two opposite PMT
   *         walls, so there is no TPC-level granularity).
   *     * `cryostats`: information per cryostat (east first, then west).
   *         Note however that only the cryostats that are mentioned in
   *         `triggerLocationBits` have the information filled.
   *        * `triggerCount`: number of triggers in this cryostat fired during
   *            the trigger window; only the first one becomes the global
   *            trigger, but we still keep the count of how many happen.
   *            Its value is `0` when trigger happened from elsewhere.
   *        * `LVDSstatus`: information per PMT wall (i.e. TPC; east first, then
   *            west) of the LVDS signals of the discriminated PMT pairs at the
   *            time of the global trigger. All bits are `0` when trigger
   *            happened elsewhere. Otherwise, the encoding is currently
   *            implemented in terms of hardware connectors as follows:
   *            * east wall:  `00<C3P2><C3P1><C3P0>00<C2P2><C2P1><C2P0>`
   *            * west wall:  `00<C1P2><C1P1><C1P0>00<C0P2><C0P1><C0P0>`
   *            
   *            For the expected matching with PMT, see the documentation of
   *            `sbn::ExtraTriggerInfo::CryostatInfo::LVDSstatus`.
   *     
   *     Information may be missing. If a count is not available, its value is
   *     set to `0` (which is an invalid value because their valid range starts
   *     with `1` since they include the current event), and if a timestamp is
   *     not available it is set to `sbn::ExtraTriggerInfo::NoTimestamp`; these
   *     two conditions can be checked with static methods
   *     `sbn::ExtraTriggerInfo::isValidTimestamp()`
   *     and `sbn::ExtraTriggerInfo::isValidCount()` respectively.
   * 
   * Besides the main data product (empty instance name) an additional
   * `std::vector<raw::ExternalTrigger>` data product with instance name
   * `"previous"` is also produced, which relays the same kind of information
   * but for the _previous_ trigger. This information also comes from the
   * trigger DAQ. If no previous trigger is available, this collection will be
   * empty.
   * 
   * 
   * Timestamps and corrections
   * ---------------------------
   * 
   * The reference trigger time is driven by the trigger fragment time, which
   * is expected to have been derived from the actual trigger time from the
   * White Rabbit system properly corrected to UTC by the board reader.
   * 
   * All absolute timestamps are corrected to be on that same scale.
   * The absolute timestamps related to the White Rabbit time are added an
   * offset to achieve this correction; this offset is stored in the data
   * product (`sbn::ExtraTriggerInfo::WRtimeToTriggerTime`).
   * 
   * 
   * Configuration
   * --------------
   * 
   * * `TrigConfigLabel` (input tag, mandatory): tag of the trigger
   *     configuration data product (see `icarus::TriggerConfigurationExtractor`
   *     module) to be used. Specifying its tag is mandatory, but if it is
   *     explicitly specified empty, the decoder will try to work around its
   *     absence.
   * * `DiagnosticOutput` (flag, default: `false`): prints on console trigger
   *     data diagnostics (including a full dump of the parsed content).
   * * `Debug` (flag, default: `false`): prints on console decoding debug
   *     information, including a dump of the trigger data fragment.
   * 
   */
  class TriggerDecoderV3 : public IDecoder
  {
    using microseconds = util::quantities::microsecond;
    using nanoseconds = util::quantities::nanosecond;
  public:
    explicit TriggerDecoderV3(fhicl::ParameterSet const &pset);
    
    virtual void consumes(art::ConsumesCollector& collector) override;
    virtual void produces(art::ProducesCollector&) override;
    virtual void configure(const fhicl::ParameterSet&) override;
    virtual void initializeDataProducts() override;
    virtual void setupRun(art::Run const& run) override;
    virtual void process_fragment(const artdaq::Fragment &fragment) override;
    virtual void outputDataProducts(art::Event &event) override;
   
  private: 
    using TriggerCollection = std::vector<raw::ExternalTrigger>;
    using TriggerPtr = std::unique_ptr<TriggerCollection>;
    using RelativeTriggerCollection = std::vector<raw::Trigger>;
    using BeamGateInfoCollection = std::vector<sim::BeamGateInfo>;
    using BeamGateInfoPtr = std::unique_ptr<BeamGateInfoCollection>;
    using ExtraInfoPtr = std::unique_ptr<sbn::ExtraTriggerInfo>;
    TriggerPtr fTrigger;
    TriggerPtr fPrevTrigger;
    std::unique_ptr<RelativeTriggerCollection> fRelTrigger;
    ExtraInfoPtr fTriggerExtra;
    BeamGateInfoPtr fBeamGateInfo;
    art::InputTag fTriggerConfigTag; ///< Data product with hardware trigger configuration.
    bool fDiagnosticOutput; ///< Produces large number of diagnostic messages, use with caution!
    bool fDebug; ///< Use this for debugging this tool
    
    long fLastEvent = 0;
    
    detinfo::DetectorTimings const fDetTimings; ///< Detector clocks and timings.
    
    /// Cached pointer to the trigger configuration of the current run, if any.
    icarus::TriggerConfiguration const* fTriggerConfiguration = nullptr;
    
    
    /// Creates a `ICARUSTriggerInfo` from a generic fragment.
    icarus::ICARUSTriggerV3Fragment makeTriggerFragment
      (artdaq::Fragment const& fragment) const;
    
    /// Parses the trigger data packet with the "standard" parser.
    icarus::ICARUSTriggerInfo parseTriggerString(std::string_view data) const;
    
    /// Parses the trigger data packet with a CSV parser.
    icarus::KeyValuesData parseTriggerStringAsCSV
      (std::string const& data) const;
    
    /// Name of the data product instance for the current trigger.
    static std::string const CurrentTriggerInstanceName;
    
    /// Name of the data product instance for the previous trigger.
    static std::string const PreviousTriggerInstanceName;
    
    static constexpr double UnknownBeamTime = std::numeric_limits<double>::max();
    
    /// Codes of gate types from the trigger hardware.
    struct TriggerGateTypes {
      static constexpr int BNB { 1 };
      static constexpr int NuMI { 2 };
      static constexpr int OffbeamBNB { 3 };
      static constexpr int OffbeamNuMI { 4 };
      static constexpr int Calib { 5 };
    }; 
    
    static std::string_view firstLine
      (std::string const& s, std::string const& endl = "\0\n\r"s);
    
    /// Combines second and nanosecond counts into a 64-bit timestamp.
    static std::uint64_t makeTimestamp(unsigned int s, unsigned int ns)
      { return s * 1000000000ULL + ns; }
    /// Returns the difference `a - b`.
    static long long int timestampDiff(std::uint64_t a, std::uint64_t b)
      { return static_cast<long long int>(a) - static_cast<long long int>(b); }

    /// Fills one `cryostat` information of `sbn::ExtraTriggerInfo`.
    static sbn::ExtraTriggerInfo::CryostatInfo unpackPrimitiveBits(
      std::size_t cryostat, bool firstEvent, unsigned long int counts,
      std::uint64_t connectors01, std::uint64_t connectors23
      );

    /// Encodes the `connectorWord` LVDS bits from the specified `cryostat`
    /// and `connector` into the format required by `sbn::ExtraTriggerInfo`.
    static std::uint64_t encodeLVDSbits
      (short int cryostat, short int connector, std::uint64_t connectorWord);
    
    static std::uint16_t encodeSectorBits
      (short int cryostat, short int connector, std::uint64_t connectorWord);
    
    /// Returns the `nBits` bits of `value` from `startBit` on.
    template <unsigned int startBit, unsigned int nBits, typename T>
    static constexpr T bits(T value);

    /// Returns the beam type corresponding to the specified trigger `source`.
    static sim::BeamType_t simGateType(sbn::triggerSource source);
    
  };


  std::string const TriggerDecoderV3::CurrentTriggerInstanceName {};
  std::string const TriggerDecoderV3::PreviousTriggerInstanceName { "previous" };
  

  TriggerDecoderV3::TriggerDecoderV3(fhicl::ParameterSet const &pset)
    : fDetTimings
      { art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob() }
  {
    this->configure(pset);
  }

  
  void TriggerDecoderV3::consumes(art::ConsumesCollector& collector) {
    collector.consumes<icarus::TriggerConfiguration, art::InRun>
      (fTriggerConfigTag);
  }
  
  void TriggerDecoderV3::produces(art::ProducesCollector& collector) 
  {
    collector.produces<TriggerCollection>(CurrentTriggerInstanceName);
    collector.produces<TriggerCollection>(PreviousTriggerInstanceName);
    collector.produces<RelativeTriggerCollection>(CurrentTriggerInstanceName);
    collector.produces<BeamGateInfoCollection>(CurrentTriggerInstanceName);
    collector.produces<sbn::ExtraTriggerInfo>(CurrentTriggerInstanceName);
  }
    

  void TriggerDecoderV3::configure(fhicl::ParameterSet const &pset) 
  {
    fTriggerConfigTag = pset.get<std::string>("TrigConfigLabel");
    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);
    fDebug = pset.get<bool>("Debug", false);
    if (pset.has_key("TimeOffset")) {
      if (auto offset = pset.get<long long int>("TimeOffset"); offset == 0) {
        mf::LogWarning("TriggerDecoder")
          << "Configuration parameter 'TimeOffset' has been dropped.\n";
      }
      else { // unforgiving: dropping non-zero offset changes output
        throw art::Exception(art::errors::Configuration)
          << "Adding offset (configuration parameter 'TimeOffset', here set to "
          << offset << " seconds) is not supported any more.\n";
      }
    }
    return;
  }
  
  void TriggerDecoderV3::initializeDataProducts()
  {
    //use until different object chosen 
    //fTrigger = new raw::Trigger();
    fTrigger = std::make_unique<TriggerCollection>();
    fPrevTrigger = std::make_unique<TriggerCollection>();
    fRelTrigger = std::make_unique<RelativeTriggerCollection>();
    fBeamGateInfo = BeamGateInfoPtr(new BeamGateInfoCollection);
    fTriggerExtra = std::make_unique<sbn::ExtraTriggerInfo>();
    return;
  }
  
  
  icarus::ICARUSTriggerV3Fragment TriggerDecoderV3::makeTriggerFragment
    (artdaq::Fragment const& fragment) const
  {
    try {
      return icarus::ICARUSTriggerV3Fragment { fragment };
    }
    catch(std::exception const& e) {
      mf::LogSystem("TriggerDecoder")
        << "Error while creating trigger fragment from:\n"
          << sbndaq::dumpFragment(fragment)
        << "\nError message: " << e.what();
      throw;
    }
    catch(...) {
      mf::LogSystem("TriggerDecoder")
        << "Unidentified exception while creating trigger fragment from:"
          << sbndaq::dumpFragment(fragment);
      throw;
    }
  } // TriggerDecoderV3::makeTriggerFragment()

  
  icarus::ICARUSTriggerInfo TriggerDecoderV3::parseTriggerString
    (std::string_view data) const
  {
    try {
      return icarus::parse_ICARUSTriggerV3String(data.data());
    }
    catch(std::exception const& e) {
      mf::LogSystem("TriggerDecoder")
        << "Error while running standard parser on " << data.length()
        << "-char long trigger string:\n==>|" << data
        << "|<==\nError message: " << e.what();
      throw;
    }
    catch(...) {
      mf::LogSystem("TriggerDecoder")
        << "Unidentified exception while running standard parser on "
        << data.length() << "-char long trigger string:\n==>|" << data << "|.";
      throw;
    }
  } // TriggerDecoderV3::parseTriggerString()


  icarus::KeyValuesData TriggerDecoderV3::parseTriggerStringAsCSV
    (std::string const& data) const
  {
    icarus::details::KeyedCSVparser parser;
    parser.addPatterns({
        { "Cryo. (EAST|WEST) Connector . and .", 1U }
        , { "Trigger Type", 1U }
      });
    std::string_view const dataLine = firstLine(data);
    try {
      return parser(dataLine);
    }
    catch(icarus::details::KeyedCSVparser::Error const& e) {
      mf::LogError("TriggerDecoder")
        << "Error parsing " << dataLine.length()
        << "-char long trigger string:\n==>|" << dataLine
        << "|<==\nError message: " << e.what() << std::endl;
      throw;
    }
  } // TriggerDecoderV3::parseTriggerStringAsCSV()
  

  void TriggerDecoderV3::setupRun(art::Run const& run) {
    
    fTriggerConfiguration = fTriggerConfigTag.empty()
      ? nullptr
      : &(run.getProduct<icarus::TriggerConfiguration>(fTriggerConfigTag))
      ;
    
  } // TriggerDecoderV3::setupRun()
  
  
  void TriggerDecoderV3::process_fragment(const artdaq::Fragment &fragment)
  {
    // artdaq_ts is reworked by the trigger board reader to match the corrected
    // trigger time; to avoid multiple (potentially inconsistent) corrections,
    // the decoder trusts it and references all the times with respect to it.
    uint64_t const artdaq_ts = fragment.timestamp();
    icarus::ICARUSTriggerV3Fragment frag { makeTriggerFragment(fragment) };
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());

    icarus::ICARUSTriggerInfo datastream_info = parseTriggerString(buffer);
    uint64_t const raw_wr_ts // this is raw, unadultered, uncorrected
      = makeTimestamp(frag.getWRSeconds(), frag.getWRNanoSeconds());
    
    // correction (explicitly converted to signed)
    int64_t const WRtimeToTriggerTime
      = static_cast<int64_t>(artdaq_ts) - raw_wr_ts;
    auto const correctWRtime = [WRtimeToTriggerTime](uint64_t time)
      { return time + WRtimeToTriggerTime; };
    assert(correctWRtime(raw_wr_ts) == artdaq_ts);
    
    //
    // we parse again the trigger string for information that was not saved
    // by the board reader in the trigger fragment nor in `datastream_info`
    //
    auto const parsedData = parseTriggerStringAsCSV(data); 
    
    unsigned int beamgate_count { std::numeric_limits<unsigned int>::max() };
    std::uint64_t beamgate_ts { artdaq_ts }; // we cheat
    /* [20210717, petrillo@slac.stanford.edu] `(pBeamGateInfo->nValues() == 3)`:
     * this is an attempt to "support" a few Run0 runs (6017 to roughly 6043)
     * which have the beam gate information truncated; this workaround should
     * be removed when there is enough ICARUS data that these runs become
     * uninteresting.
     */
    if (auto pBeamGateInfo = parsedData.findItem("Beam_TS")) {
      /*
       * The Veto Business:
       * 
       * to better cover the pre-spill time, a trigger primitive is issued
       * in advance of the actual beam gate, but during the additional time
       * the global trigger is vetoed.
       * So the physical beam gate signal starts in advance of the beam gate,
       * and teh actual beam gate is expected after the veto time has passed.
       * The hardware tagging the beam gate time does not know of this, and it
       * marks the beam gate to start at the opening of the physical signal
       * (at the beginning of the vetoed time).
       * This is all a trick for PMT readout, nothing should depend on it,
       * so we correct that as soon as we can and then forget it.
       * The amount of veto time is written in the trigger configuration;
       * if that configuration is not available, this code assumes there was
       * no veto.
       */
      int64_t const triggerVetoDurationNS
        = fTriggerConfiguration? fTriggerConfiguration->vetoDelay: 0LL;
      
      
      // if gate information is found, it must be complete
      beamgate_count = pBeamGateInfo->getNumber<unsigned int>(0U);
      
      uint64_t const raw_bg_ts = makeTimestamp( // raw and uncorrected too...
        pBeamGateInfo->getNumber<unsigned int>(1U),
        pBeamGateInfo->getNumber<unsigned int>(2U)
        )
        + triggerVetoDurationNS // ... but remove the veto time
        ;
      
      // assuming the raw times from the fragment are on the same time scale
      // (same offset corrections)
      beamgate_ts += raw_bg_ts - raw_wr_ts;
      
    } // if has gate information
    std::uint64_t enablegate_ts { artdaq_ts };
    if (auto pEnableGateInfo = parsedData.findItem("Enable_TS"))
    {
      // if gate information is found, it must be complete
      //enablegate_count = pEnableGateInfo->getNumber<unsigned int>(0U);

      uint64_t const raw_en_ts = makeTimestamp( // raw and uncorrected too
                                               pEnableGateInfo->getNumber<unsigned int>(1U),
                                               pEnableGateInfo->getNumber<unsigned int>(2U)
                                                );

      // assuming the raw times from the fragment are on the same time scale 
      // (same offset corrections)

      enablegate_ts += raw_en_ts - raw_wr_ts;
    } // if has gate information

    // --- END ---- TEMPORARY --------------------------------------------------
    
    if(fDiagnosticOutput || fDebug)
    {
      std::cout << "Full Timestamp = " << artdaq_ts
        << "\nBeam gate " << beamgate_count << " at "
        << (beamgate_ts/1'000'000'000) << "." << std::setfill('0')
        << std::setw(9) << (beamgate_ts%1'000'000'000) << std::setfill(' ')
        << " s (" << timestampDiff(beamgate_ts, artdaq_ts)
        << " ns relative to trigger)"
        << "\nParsed data (from " << data.size() << " characters): "
        << parsedData << std::endl;
      
      if (fDebug) { // this grows tiresome quickly when processing many events
        std::cout << "Trigger packet content:\n" << data
          << "\nFull trigger fragment dump:"
          << sbndaq::dumpFragment(fragment) << std::endl;
      }
    }
    //
    // extra trigger info
    //
    unsigned int const triggerID = datastream_info.wr_event_no;
    sbn::triggerSource beamGateBit;
    switch (datastream_info.gate_type) {
      case TriggerGateTypes::BNB:{
        beamGateBit = sbn::triggerSource::BNB;
      if(datastream_info.trigger_type == 0)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesBNBMaj();
        fTriggerExtra->previousTriggerTimestamp = frag.getLastTimestampBNBMaj();
        fTriggerExtra->gateCount = datastream_info.gate_id_BNB;
        fTriggerExtra->triggerCount = frag.getTotalTriggerBNBMaj();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerBNBMaj();
      }
      else if(datastream_info.trigger_type == 1)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesBNBMinbias();
          fTriggerExtra->previousTriggerTimestamp = frag.getLastTimestampBNBMinbias();
          fTriggerExtra->gateCount = datastream_info.gate_id_BNB;
          fTriggerExtra->triggerCount = frag.getTotalTriggerBNBMinbias();
          fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerBNBMinbias();
      }
        break;
      }
      case TriggerGateTypes::NuMI:{
        beamGateBit = sbn::triggerSource::NuMI;
      if(datastream_info.trigger_type == 0)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesNuMIMaj();
        fTriggerExtra->previousTriggerTimestamp = frag.getLastTimestampNuMIMaj();
        fTriggerExtra->gateCount = datastream_info.gate_id_NuMI;
        fTriggerExtra->triggerCount = frag.getTotalTriggerNuMIMaj();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerNuMIMaj();
      }
      else if(datastream_info.trigger_type == 1)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesNuMIMinbias();
        fTriggerExtra->previousTriggerTimestamp = frag.getLastTimestampNuMIMinbias();
        fTriggerExtra->gateCount = datastream_info.gate_id_NuMI;
        fTriggerExtra->triggerCount = frag.getTotalTriggerNuMIMinbias();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerNuMIMinbias();
      }
        break;
      }
      case TriggerGateTypes::OffbeamBNB:{
        beamGateBit = sbn::triggerSource::OffbeamBNB;
      if(datastream_info.trigger_type == 0)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesBNBOffMaj();
        fTriggerExtra->previousTriggerTimestamp= frag.getLastTimestampBNBOffMaj();
        fTriggerExtra->gateCount = datastream_info.gate_id_BNBOff;
        fTriggerExtra->triggerCount = frag.getTotalTriggerBNBOffMaj();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerBNBOffMaj();
      }
      else if(datastream_info.trigger_type == 1)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesBNBOffMinbias();
        fTriggerExtra->previousTriggerTimestamp= frag.getLastTimestampBNBOffMinbias();
        fTriggerExtra->gateCount = datastream_info.gate_id_BNBOff;
        fTriggerExtra->triggerCount = frag.getTotalTriggerBNBOffMinbias();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerBNBOffMinbias();
      }
        break;
      }
      case TriggerGateTypes::OffbeamNuMI:{
        beamGateBit = sbn::triggerSource::OffbeamNuMI;
      if(datastream_info.trigger_type == 0)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesNuMIOffMaj();
        fTriggerExtra->previousTriggerTimestamp= frag.getLastTimestampNuMIOffMaj();
        fTriggerExtra->gateCount = datastream_info.gate_id_NuMIOff;
        fTriggerExtra->triggerCount = frag.getTotalTriggerNuMIOffMaj();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerNuMIOffMaj();
      }
      if(datastream_info.trigger_type == 1)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesNuMIOffMinbias();
        fTriggerExtra->previousTriggerTimestamp= frag.getLastTimestampNuMIOffMinbias();
        fTriggerExtra->gateCount = datastream_info.gate_id_NuMIOff;
        fTriggerExtra->triggerCount = frag.getTotalTriggerNuMIOffMinbias();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerNuMIOffMinbias();
      }
        break;
      }
      case TriggerGateTypes::Calib:{
        beamGateBit = sbn::triggerSource::Calib;
      if(datastream_info.trigger_type == 0)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesCalibMaj();
        fTriggerExtra->previousTriggerTimestamp = frag.getLastTimestampCalibMaj();
        //fTriggerExtra->gateCount = datastream_info.gate_id_calib;
        fTriggerExtra->triggerCount = frag.getTotalTriggerCalibMaj();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerCalibMaj();
      }
      if(datastream_info.trigger_type == 1)
      {
        fTriggerExtra->gateCountFromPreviousTrigger = frag.getDeltaGatesCalibMinbias();
        fTriggerExtra->previousTriggerTimestamp = frag.getLastTimestampCalibMinbias();
        //fTriggerExtra->gateCount = datastream_info.gate_id_calib;
        fTriggerExtra->triggerCount = frag.getTotalTriggerCalibMinbias();
        fTriggerExtra->anyTriggerCountFromPreviousTrigger = triggerID - frag.getLastTriggerCalibMinbias();
      }
      break;
      }
      default:                            beamGateBit = sbn::triggerSource::Unknown;
    } // switch gate type
    
    fTriggerExtra->sourceType = beamGateBit;
    fTriggerExtra->triggerType = static_cast<sbn::triggerType>(datastream_info.trigger_type);
    fTriggerExtra->triggerTimestamp = artdaq_ts;
    fTriggerExtra->beamGateTimestamp = beamgate_ts;
    fTriggerExtra->enableGateTimestamp = enablegate_ts;
    fTriggerExtra->triggerID = triggerID; //all triggers (event ID)
    fTriggerExtra->gateID = datastream_info.gate_id; //all gate types (gate ID)
    fTriggerExtra->anyGateCountFromAnyPreviousTrigger = frag.getDeltaGates();
    fTriggerExtra->anyPreviousTriggerTimestamp = frag.getLastTimestamp();
    sbn::triggerSource previousTriggerSourceBit;
    if(frag.getLastTriggerType() == 1)
      previousTriggerSourceBit = sbn::triggerSource::BNB;
    else if(frag.getLastTriggerType() == 2)
      previousTriggerSourceBit = sbn::triggerSource::NuMI;
    else if(frag.getLastTriggerType() == 3)
      previousTriggerSourceBit = sbn::triggerSource::OffbeamBNB;
    else if(frag.getLastTriggerType() == 4)
      previousTriggerSourceBit = sbn::triggerSource::OffbeamNuMI;
    else if(frag.getLastTriggerType() == 5)
      previousTriggerSourceBit = sbn::triggerSource::Calib;
    else
      previousTriggerSourceBit = sbn::triggerSource::Unknown;
    fTriggerExtra->anyPreviousTriggerSourceType = previousTriggerSourceBit;

    fTriggerExtra->WRtimeToTriggerTime = WRtimeToTriggerTime;
    sbn::bits::triggerLocationMask locationMask;
    // trigger location: 0x01=EAST; 0x02=WEST; 0x07=ALL
    int const triggerLocation = parsedData.getItem("Trigger Source").getNumber<int>(0);
    if(triggerLocation == 1)
      locationMask = mask(sbn::triggerLocation::CryoEast);
    else if(triggerLocation == 2)
      locationMask = mask(sbn::triggerLocation::CryoWest);
    else if(triggerLocation == 7)
      locationMask = mask(sbn::triggerLocation::CryoEast, sbn::triggerLocation::CryoWest);
    fTriggerExtra->triggerLocationBits = locationMask;
    
    //
    // fill sbn::ExtraTriggerInfo::cryostats
    //
    auto setCryoInfo = [&extra=*fTriggerExtra,isFirstEvent=(triggerID <= 1)]
      (std::size_t cryo, icarus::KeyValuesData const& data)
      {
        std::string const SIDE = (cryo == sbn::ExtraTriggerInfo::EastCryostat)
          ? "Cryo1 EAST": "Cryo2 WEST";
        extra.cryostats[cryo] = unpackPrimitiveBits(
          cryo, isFirstEvent,
          data.getItem(SIDE + " counts").getNumber<unsigned long int>(0),
          data.getItem(SIDE + " Connector 0 and 1").getNumber<std::uint64_t>(0, 16),
          data.getItem(SIDE + " Connector 2 and 3").getNumber<std::uint64_t>(0, 16)
          );
      };
    
    if (triggerLocation & 1)
      setCryoInfo(sbn::ExtraTriggerInfo::EastCryostat, parsedData);
    if (triggerLocation & 2)
      setCryoInfo(sbn::ExtraTriggerInfo::WestCryostat, parsedData);
    
    //
    // absolute time trigger (raw::ExternalTrigger)
    //
    fTrigger->emplace_back
      (fTriggerExtra->triggerID, artdaq_ts);//fTriggerExtra->triggerTimestamp);
    
    //
    // previous absolute time trigger (raw::ExternalTrigger)
    //
    uint64_t lastTrigger = 0;
    if(fTriggerExtra->triggerID == 1)
    {
      fLastEvent = 0;
    }
    else 
    {
      fLastEvent = fTriggerExtra->triggerID - 1;
      lastTrigger = fTriggerExtra->anyPreviousTriggerTimestamp;
      fPrevTrigger->emplace_back(fLastEvent, lastTrigger);
    }
    
    //
    // beam gate
    //
    
    // beam gate width (read in microseconds, but we use it in nanoseconds);
    // we need to protect from some defective configurations where the (unused)
    // beam gate duration is smaller than the veto time
    icarus::TriggerConfiguration::GateConfig const* gateConfig
      = fTriggerConfiguration
      ? &(fTriggerConfiguration->gateConfig.at(value(beamGateBit))): nullptr
      ;
    nanoseconds const gateWidth = microseconds{
      (gateConfig && (gateConfig->gateWidth > fTriggerConfiguration->vetoDelay))
      ? fTriggerConfiguration->getGateWidth(value(beamGateBit))
      : 0.0
      };
    
    // beam gate, defined in simulation time, which should match beam gate time
    // ... trivial:
    fBeamGateInfo->emplace_back(0, gateWidth.value(), simGateType(beamGateBit));
    
    //
    // relative time trigger (raw::Trigger)
    //
    nanoseconds const gateStartFromTrigger{
      static_cast<double>(timestampDiff
        (fTriggerExtra->beamGateTimestamp, fTriggerExtra->triggerTimestamp)
        )
      }; // narrowing!!
    auto const elecGateStart = fDetTimings.TriggerTime() + gateStartFromTrigger;
    fRelTrigger->emplace_back(
      fTriggerExtra->triggerID,                               // counter
      fDetTimings.TriggerTime().value(),                      // trigger_time
      elecGateStart.value(),                                  // beamgate_time
      mask(beamGateBit)                                       // bits
      );
    
    //Once we have full trigger data object, set up and place information into there
    return;
  }

  void TriggerDecoderV3::outputDataProducts(art::Event &event)
  {
    //Place trigger data object into raw data store 
    event.put(std::move(fTrigger), CurrentTriggerInstanceName);
    event.put(std::move(fRelTrigger), CurrentTriggerInstanceName);
    event.put(std::move(fPrevTrigger), PreviousTriggerInstanceName);
    event.put(std::move(fBeamGateInfo), CurrentTriggerInstanceName);
    event.put(std::move(fTriggerExtra));
    return;
  }

  std::string_view TriggerDecoderV3::firstLine
    (std::string const& s, std::string const& endl /* = "\0\n\r" */)
  {
    return { s.data(), std::min(s.find_first_of(endl), s.size()) };
  }
  
  
  sbn::ExtraTriggerInfo::CryostatInfo TriggerDecoderV3::unpackPrimitiveBits(
    std::size_t cryostat, bool firstEvent, unsigned long int counts,
    std::uint64_t connectors01, std::uint64_t connectors23
  ) {
    sbn::ExtraTriggerInfo::CryostatInfo cryoInfo;
    
    // there is (or was?) a bug on the first event in the run,
    // which would make this triggerCount wrong
    cryoInfo.triggerCount = firstEvent? 0UL: counts,
    
    cryoInfo.LVDSstatus = {
      encodeLVDSbits(cryostat, 2 /* any of the connectors */, connectors23),
      encodeLVDSbits(cryostat, 0 /* any of the connectors */, connectors01)
      };
    
    cryoInfo.SectorStatus = {
      encodeSectorBits(cryostat, 2 /* any of the connectors */, connectors23),
      encodeSectorBits(cryostat, 0 /* any of the connectors */, connectors01)
      };
    
    return cryoInfo;
  } // TriggerDecoderV3::unpackPrimitiveBits()
  
  
  std::uint64_t TriggerDecoderV3::encodeLVDSbits
    (short int cryostat, short int connector, std::uint64_t connectorWord)
  {
    /*
     * Encoding of the LVDS channels from the trigger:
     * * east wall:  `00<C0P2><C0P1><C0P0>00<C1P2><C1P1><C1P0>`
     * * west wall:  `00<C2P2><C2P1><C2P0>00<C3P2><C3P1><C3P0>`
     * The prescription from `sbn::ExtraTriggerInfo` translates into:
     * * east wall:  `00<C3P2><C3P1><C3P0>00<C2P2><C2P1><C2P0>`
     * * west wall:  `00<C1P2><C1P1><C1P0>00<C0P2><C0P1><C0P0>`
     * Therefore, the two 32-bit half-words need to be swapped
     * This holds for both cryostats, and both walls.
     * (note that the `00` bits may actually contain information from adders)
     */

    std::uint64_t lsw = connectorWord & 0x00FFFFFFULL;
    std::uint64_t msw = (connectorWord >> 32ULL) & 0x00FFFFFFULL;
    assert((connectorWord & 0x00FF'FFFF'00FF'FFFF) == ((msw << 32ULL) | lsw));
    std::swap(lsw, msw);
    return (msw << 32ULL) | lsw;
  } // TriggerDecoderV3::encodeLVDSbits()
  
  
  template <unsigned int startBit, unsigned int nBits, typename T>
  constexpr T TriggerDecoderV3::bits(T value) {
    constexpr T nBitsMask = (nBits == sizeof(T)*8)? ~T{0}: T((1 << nBits) - 1);
    return (value >> T(startBit)) & nBitsMask;
  }
  
  
  std::uint16_t TriggerDecoderV3::encodeSectorBits
    (short int cryostat, short int connector, std::uint64_t connectorWord)
  {
    /*
     * Encoding of the LVDS channels from the trigger (cf. `encodeLVDSbits()`):
     *  * connector 0: east wall south, LSB the south-most one
     *  * connector 1: east wall north, LSB the south-most one
     *  * connector 2: west wall south, LSB the south-most one
     *  * connector 3: west wall north, LSB the south-most one
     * Target:
     *  * [0]: 00000000 00nnnsss (east wall)
     *  * [1]: 00000000 00nnnsss (west wall)
     */
    
    // connector (32 bit): 00000SSS LVDSLVDS LVDSLVDS LVDSLVDS
    constexpr std::size_t BitsPerConnector = 32;
    constexpr std::size_t FirstSectorBit = 27;
    constexpr std::size_t SectorBitsPerConnector = 3;
    
    // TODO check the order of the bits: [ ] the blocks of 3 and [ ] within each block
    return static_cast<std::uint16_t>(
        (bits<BitsPerConnector+FirstSectorBit, SectorBitsPerConnector>(connectorWord) << SectorBitsPerConnector)
      | (bits<FirstSectorBit,                  SectorBitsPerConnector>(connectorWord)                          )
      );
  } // TriggerDecoderV3::encodeSectorBits()
  
  
  sim::BeamType_t TriggerDecoderV3::simGateType(sbn::triggerSource source)
  {
    switch (source) {
      case sbn::triggerSource::BNB:
      case sbn::triggerSource::OffbeamBNB:
        return sim::kBNB;
      case sbn::triggerSource::NuMI:
      case sbn::triggerSource::OffbeamNuMI:
        return sim::kNuMI;
      case sbn::triggerSource::Calib:
        return sim::kUnknown;
      default:
        mf::LogWarning("TriggerDecoder") << "Unsupported trigger source " << name(source);
        return sim::kUnknown;
    } // switch source
  } // TriggerDecoderV3::simGateType()
  
  
  DEFINE_ART_CLASS_TOOL(TriggerDecoderV3)

}
