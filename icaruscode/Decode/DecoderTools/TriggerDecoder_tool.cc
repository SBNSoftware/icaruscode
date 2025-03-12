////////////////////////////////////////////////
//   
//    File: TriggerDecoder_tool.cc
//       
//    Description: Starting point for extracting ICARUS trigger fragment information into LArSoft object TBD 
//
//    Author: Jacob Zettlemoyer, FNAL
//
///////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
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

#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerUDPFragment.hh"

#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/Decode/DecoderTools/IDecoder.h"
#include "sbnobj/Common/Trigger/BeamBits.h" 
#include "icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h" // dumpFragment()
#include "icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h"
#include "icarusalg/Utilities/BinaryDumpUtils.h" // hexdump() DEBUG

#include <cstdlib>
#include <iostream>
#include <iomanip> // std::setw(), std::setfill()
#include <string_view>
#include <memory>


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
   *     * `BeamGateTime()`: relative time of the announced arrival of the beam
   *         (currently not available) also in
   *         @ref DetectorClocksElectronicsTime "electronics time scale".
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
   *     * `Start()`: relative time of the announced arrival of the beam
   *         (currently not available), in
   *         @ref DetectorClocksSimulationTime "simulation time scale".
   *     * `Width()`: duration of the gate, in nanoseconds; currently set to a
   *         nominal value.
   *     * `BeamType()`: the type of the beam gate being described (BNB, NuMI).
   * * `sbn::ExtraTriggerInfo`: the most complete collection of information,
   *     duplicating also some from the other data products. Some of the
   *     information is not available yet: if a count is not available, its
   *     value is set to `0` (which is an invalid value because their valid
   *     range starts with `1` since they include the current event), and if a
   *     timestamp is not available it is set to
   *     `sbn::ExtraTriggerInfo::NoTimestamp`; these two conditions can be
   *     checked with static methods 
   *     `sbn::ExtraTriggerInfo::isValidTimestamp()` and 
   *     `sbn::ExtraTriggerInfo::isValidCount()` respectively.
   *     Note that differently from the usual, this is a _single object_, not
   *     a collection; also, this data product has no instance name.
   *     The information already available includes:
   *     * `sourceType`: the type of beam or trigger source, a value from
   *         `sbn::triggerSource` (equivalent to `raw::Trigger::TriggerBits()`,
   *         but in the form of an enumerator rather than a bit mask).
   *     * `triggerTimestamp`: same as `raw::ExternalTrigger::GetTrigTime()`
   *         (nanoseconds from the Epoch, Coordinated Universal Time).
   *     * `beamGateTimestamp`: absolute time of the beam gate opening as
   *         reported by the trigger hardware, directly comparable with
   *         `triggerTimestamp` (same scale and unit).
   *     * `triggerID`: same as `raw::ExternalTrigger::GetTrigID()`. Should
   *         match the event number.
   *     * `gateID`: the count of gates since the beginning of the run, as
   *         reported by the trigger hardware.
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
   */
  class TriggerDecoder : public IDecoder
  {
    using nanoseconds = util::quantities::nanosecond;
  public:
    explicit TriggerDecoder(fhicl::ParameterSet const &pset);
    
    virtual void produces(art::ProducesCollector&) override;
    virtual void configure(const fhicl::ParameterSet&) override;
    virtual void initializeDataProducts() override;
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
    bool fDiagnosticOutput; //< Produces large number of diagnostic messages, use with caution!
    bool fDebug; //< Use this for debugging this tool
    //Add in trigger data member information once it is selected, current LArSoft object likely not enough as is
    
    // uint64_t fLastTimeStamp = 0;
   
    long fLastEvent = 0;
    
    detinfo::DetectorTimings const fDetTimings; ///< Detector clocks and timings.
    
    /// Creates a `ICARUSTriggerInfo` from a generic fragment.
    icarus::ICARUSTriggerUDPFragment makeTriggerFragment
      (artdaq::Fragment const& fragment) const;
    
    /// Parses the trigger data packet with the "standard" parser.
    icarus::ICARUSTriggerInfo parseTriggerString(std::string_view data) const;
    
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
    }; 
    
    static constexpr nanoseconds BNBgateDuration { 1600. };
    static constexpr nanoseconds NuMIgateDuration { 9500. };
    
    static std::string_view firstLine
      (std::string const& s, std::string const& endl = "\0\n\r"s);
    
    /// Combines second and nanosecond counts into a 64-bit timestamp.
    static std::uint64_t makeTimestamp(unsigned int s, unsigned int ns)
      { return s * 1000000000ULL + ns; }
    /// Returns the difference `a - b`.
    static long long int timestampDiff(std::uint64_t a, std::uint64_t b)
      { return static_cast<long long int>(a) - static_cast<long long int>(b); }
    
  };


  std::string const TriggerDecoder::CurrentTriggerInstanceName {};
  std::string const TriggerDecoder::PreviousTriggerInstanceName { "previous" };
  

  TriggerDecoder::TriggerDecoder(fhicl::ParameterSet const &pset)
    : fDetTimings
      { art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob() }
  {
    this->configure(pset);
  }

  
  void TriggerDecoder::produces(art::ProducesCollector& collector) 
  {
    collector.produces<TriggerCollection>(CurrentTriggerInstanceName);
    collector.produces<TriggerCollection>(PreviousTriggerInstanceName);
    collector.produces<RelativeTriggerCollection>(CurrentTriggerInstanceName);
    collector.produces<BeamGateInfoCollection>(CurrentTriggerInstanceName);
    collector.produces<sbn::ExtraTriggerInfo>(CurrentTriggerInstanceName);
  }
    

  void TriggerDecoder::configure(fhicl::ParameterSet const &pset) 
  {
    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);
    fDebug = pset.get<bool>("Debug", false);
    return;
  }
  
  void TriggerDecoder::initializeDataProducts()
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
  
  
  icarus::ICARUSTriggerUDPFragment TriggerDecoder::makeTriggerFragment
    (artdaq::Fragment const& fragment) const
  {
    try {
      return icarus::ICARUSTriggerUDPFragment { fragment };
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
  } // TriggerDecoder::parseTriggerString()

  
  icarus::ICARUSTriggerInfo TriggerDecoder::parseTriggerString
    (std::string_view data) const
  {
    try {
      return icarus::parse_ICARUSTriggerString(data.data());
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
  } // TriggerDecoder::parseTriggerString()

  

  void TriggerDecoder::process_fragment(const artdaq::Fragment &fragment)
  {
    // artdaq_ts is reworked by the trigger board reader to match the corrected
    // trigger time; to avoid multiple (potentially inconsistent) corrections,
    // the decoder trusts it and references all the times with respect to it.
    uint64_t const artdaq_ts = fragment.timestamp();
    icarus::ICARUSTriggerUDPFragment frag { makeTriggerFragment(fragment) };
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
    
    // --- END ---- TEMPORARY --------------------------------------------------
    int gate_type = datastream_info.gate_type;
    long delta_gates_bnb [[maybe_unused]] = frag.getDeltaGatesBNB();
    long delta_gates_numi [[maybe_unused]] = frag.getDeltaGatesOther(); //this is actually NuMI due to abrupt changes in trigger board logic
    long delta_gates_other [[maybe_unused]] = frag.getDeltaGatesNuMI();
    uint64_t lastTrigger = 0;
    
    // --- BEGIN -- TEMPORARY --------------------------------------------------
    // remove this part when the beam gate timestamp is available via fragment
    // or via the parser
    auto const parsedData = icarus::details::KeyedCSVparser{}(firstLine(data));
    unsigned int beamgate_count { std::numeric_limits<unsigned int>::max() };
    std::uint64_t beamgate_ts { artdaq_ts }; // we cheat
    /* [20210717, petrillo@slac.stanford.edu] `(pBeamGateInfo->nValues() == 3)`:
     * this is an attempt to "support" a few Run0 runs (6017 to roughly 6043)
     * which have the beam gate information truncated; this workaround should
     * be removed when there is enough ICARUS data that these runs become
     * uninteresting.
     */
    if (auto pBeamGateInfo = parsedData.findItem("Beam_TS");
      pBeamGateInfo && (pBeamGateInfo->nValues() == 3)
    ) {
      // if gate information is found, it must be complete
      beamgate_count = pBeamGateInfo->getNumber<unsigned int>(0U);
      
      uint64_t const raw_bg_ts = makeTimestamp( // raw and uncorrected too
        pBeamGateInfo->getNumber<unsigned int>(1U),
        pBeamGateInfo->getNumber<unsigned int>(2U)
        );
      
      // assuming the raw times from the fragment are on the same time scale
      // (same offset corrections)
      beamgate_ts += raw_bg_ts - raw_wr_ts;
      
    } // if has gate information
    // --- END ---- TEMPORARY --------------------------------------------------
    
    if(fDiagnosticOutput || fDebug)
    {
      std::cout << "Full Timestamp = " << artdaq_ts
        << "\nBeam gate " << beamgate_count << " at "
        << (beamgate_ts/1'000'000'000) << "." << std::setfill('0')
        << std::setw(9) << (beamgate_ts%1'000'000'000) << std::setfill(' ')
        << " s (" << timestampDiff(beamgate_ts, artdaq_ts)
        << " ns relative to trigger)" << std::endl;
      
      // note that this parsing is independent from the one used above
      std::string_view const dataLine = firstLine(data);
      try {
        auto const parsedData = icarus::details::KeyedCSVparser{}(dataLine);
        std::cout << "Parsed data (from " << dataLine.size() << " characters): "
          << parsedData << std::endl;
      }
      catch(icarus::details::KeyedCSVparser::Error const& e) {
        mf::LogError("TriggerDecoder")
          << "Error parsing " << dataLine.length()
          << "-char long trigger string:\n==>|" << dataLine
          << "|<==\nError message: " << e.what() << std::endl;
        throw;
      }
      
      if (fDebug) { // this grows tiresome quickly when processing many events
        std::cout << "Trigger packet content:\n" << dataLine
          << "\nFull trigger fragment dump:"
          << sbndaq::dumpFragment(fragment) << std::endl;
      }
    }
    
    //
    // extra trigger info
    //
    sbn::triggerSource beamGateBit;
    switch (gate_type) {
      case TriggerGateTypes::BNB:         beamGateBit = sbn::triggerSource::BNB;         break;
      case TriggerGateTypes::NuMI:        beamGateBit = sbn::triggerSource::NuMI;        break;
      case TriggerGateTypes::OffbeamBNB:  beamGateBit = sbn::triggerSource::OffbeamBNB;  break;
      case TriggerGateTypes::OffbeamNuMI: beamGateBit = sbn::triggerSource::OffbeamNuMI; break;
      default:                            beamGateBit = sbn::triggerSource::Unknown;
    } // switch gate_type
    
    fTriggerExtra->sourceType = beamGateBit;
    fTriggerExtra->triggerTimestamp = artdaq_ts;
    fTriggerExtra->beamGateTimestamp = beamgate_ts;
    fTriggerExtra->triggerID = datastream_info.wr_event_no;
    fTriggerExtra->gateID = datastream_info.gate_id;
    /* TODO (may need to add WRtimeToTriggerTime to some timestamps):
    fTriggerExtra->triggerCount
    fTriggerExtra->gateCount
    fTriggerExtra->gateCountFromPreviousTrigger
    fTriggerExtra->anyTriggerCountFromPreviousTrigger
    fTriggerExtra->anyGateCountFromPreviousTrigger
    fTriggerExtra->anyPreviousTriggerSourceType
    fTriggerExtra->anyGateCountFromAnyPreviousTrigger
    fTriggerExtra->previousTriggerTimestamp
    fTriggerExtra->anyPreviousTriggerTimestamp
    */
    fTriggerExtra->WRtimeToTriggerTime = WRtimeToTriggerTime;
    
    //
    // absolute time trigger (raw::ExternalTrigger)
    //
    fTrigger->emplace_back
      (fTriggerExtra->triggerID, fTriggerExtra->triggerTimestamp);
    
    //
    // previous absolute time trigger (raw::ExternalTrigger)
    //
    if(fTriggerExtra->triggerID == 1)
    {
      fLastEvent = 0;
    }
    else 
    {
      fLastEvent = fTriggerExtra->triggerID - 1;
      lastTrigger = frag.getLastTimestampBNB();
      fPrevTrigger->emplace_back(fLastEvent, lastTrigger);
    }
    
    //
    // beam gate
    //
    // defined in simulation time, which should match beam gate time... trivial:
    nanoseconds const gateStartFromTrigger{
      static_cast<double>(timestampDiff
        (fTriggerExtra->beamGateTimestamp, fTriggerExtra->triggerTimestamp)
        )
      }; // narrowing!!
    switch (gate_type) {
      case TriggerGateTypes::BNB:
        fBeamGateInfo->emplace_back(0.0, BNBgateDuration.value(), sim::kBNB);
        break;
      case TriggerGateTypes::NuMI:
        fBeamGateInfo->emplace_back(0.0, NuMIgateDuration.value(), sim::kNuMI);
        break;
      case TriggerGateTypes::OffbeamBNB:
        fBeamGateInfo->emplace_back(0.0, BNBgateDuration.value(), sim::kBNB);
        break;
      case TriggerGateTypes::OffbeamNuMI:
        fBeamGateInfo->emplace_back(0.0, NuMIgateDuration.value(), sim::kNuMI);
        break;
      default:
        mf::LogWarning("TriggerDecoder") << "Unsupported gate type #" << gate_type;
    } // switch gate_type
    
    //
    // relative time trigger (raw::Trigger)
    //
    auto const elecGateStart = fDetTimings.TriggerTime() + gateStartFromTrigger;
    fRelTrigger->emplace_back(
      static_cast<unsigned int>(datastream_info.wr_event_no), // counter
      fDetTimings.TriggerTime().value(),                      // trigger_time
      elecGateStart.value(),                                  // beamgate_time
      mask(beamGateBit)                                       // bits
      );
    
    //Once we have full trigger data object, set up and place information into there
    return;
  }

  void TriggerDecoder::outputDataProducts(art::Event &event)
  {
    //Place trigger data object into raw data store 
    event.put(std::move(fTrigger), CurrentTriggerInstanceName);
    event.put(std::move(fRelTrigger), CurrentTriggerInstanceName);
    event.put(std::move(fPrevTrigger), PreviousTriggerInstanceName);
    event.put(std::move(fBeamGateInfo), CurrentTriggerInstanceName);
    event.put(std::move(fTriggerExtra));
    return;
  }

  std::string_view TriggerDecoder::firstLine
    (std::string const& s, std::string const& endl /* = "\0\n\r" */)
  {
    return { s.data(), std::min(s.find_first_of(endl), s.size()) };
  }
  
  
  DEFINE_ART_CLASS_TOOL(TriggerDecoder)

}





  
