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
#include "lardataobj/RawData/ExternalTrigger.h" //JCZ: TBD, placeholder for now to represent the idea
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger
#include "lardataobj/Simulation/BeamGateInfo.h" //JCZ:, another placeholder I am unsure if this and above will be combined at some point into a dedicated object 

#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerUDPFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"
// #include "sbnobj/Common/Trigger/BeamBits.h" // maybe future location of:
#include "icaruscode/Decode/BeamBits.h" // sbn::triggerSource
#include "icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h" // dumpFragment()
#include "icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h"
#include "icarusalg/Utilities/BinaryDumpUtils.h" // hexdump() DEBUG

#include <cstdlib>
#include <iostream>
#include <memory>

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
   * 
   * Besides the main data product (empty instance name) an additional
   * `std::vector<raw::ExternalTrigger>` data product with instance name
   * `"previous"` is also produced, which relays the same kind of information
   * but for the _previous_ trigger. This information also comes from the
   * trigger DAQ. If no previous trigger is available, this collection will be
   * empty.
   * 
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
    TriggerPtr fTrigger;
    TriggerPtr fPrevTrigger;
    std::unique_ptr<RelativeTriggerCollection> fRelTrigger;
    BeamGateInfoPtr fBeamGateInfo; 
    bool fDiagnosticOutput; //< Produces large number of diagnostic messages, use with caution!
    bool fDebug; //< Use this for debugging this tool
    int fOffset; //< Use this to determine additional correction needed for TAI->UTC conversion from White Rabbit timestamps. Needs to be changed if White Rabbit firmware is changed and the number of leap seconds changes! 
    //Add in trigger data member information once it is selected, current LArSoft object likely not enough as is
    
    // uint64_t fLastTimeStamp = 0;
   
    long fLastEvent = 0;
    
    detinfo::DetectorTimings const fDetTimings; ///< Detector clocks and timings.
    
    /// Name of the data product instance for the current trigger.
    static std::string const CurrentTriggerInstanceName;
    
    /// Name of the data product instance for the previous trigger.
    static std::string const PreviousTriggerInstanceName;
    
    static constexpr double UnknownBeamTime = std::numeric_limits<double>::max();
    
    /// Codes of gate types from the trigger hardware.
    struct TriggerGateTypes {
      static constexpr int BNB { 1 };
      static constexpr int NuMI { 2 };
    }; 
    
    static constexpr nanoseconds BNBgateDuration { 1600. };
    static constexpr nanoseconds NuMIgateDuration { 9500. };
    
    static std::string_view firstLine(std::string const& s, const char* endl = "\0\n\r");
    
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
  }
    

  void TriggerDecoder::configure(fhicl::ParameterSet const &pset) 
  {
    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);
    fDebug = pset.get<bool>("Debug", false);
    fOffset = pset.get<int>("TAIOffset", 2'000'000'000);
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
    return;
  }

  void TriggerDecoder::process_fragment(const artdaq::Fragment &fragment)
  {
    uint64_t artdaq_ts = fragment.timestamp();
    icarus::ICARUSTriggerUDPFragment frag(fragment);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerString(buffer);
    uint64_t wr_ts = datastream_info.getNanoseconds_since_UTC_epoch() + fOffset;
    int gate_type = datastream_info.gate_type;
    long delta_gates_bnb = frag.getDeltaGatesBNB();
    long delta_gates_numi = frag.getDeltaGatesOther(); //this is actually NuMI due to abrupt changes in trigger board logic
    long delta_gates_other = frag.getDeltaGatesNuMI();
    uint64_t lastTrigger = 0;
    if(fDiagnosticOutput || fDebug)
    {
      std::cout << "Full Timestamp = " << wr_ts << std::endl;
      auto const cross_check = static_cast<signed long long int>(wr_ts)
        - static_cast<signed long long int>(artdaq_ts);
      if(abs(cross_check) > 1000)
        std::cout << "Loaded artdaq TS and fragment data TS are > 1 ms different! They are " << cross_check << " nanoseconds different!" << std::endl;
      std::cout << delta_gates_bnb << " " << delta_gates_numi << " " << delta_gates_other << std::endl; // nonsensical print statement to avoid error that I don't use these...until we have an object to store them in...    
      
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
        std::cout << "Full trigger fragment dump:"
          << sbndaq::dumpFragment(fragment) << std::endl;
      }
    }
    
    //
    // absolute time trigger (raw::ExternalTrigger)
    //
    fTrigger->emplace_back(datastream_info.wr_event_no, wr_ts);
    
    //
    // previous absolute time trigger (raw::ExternalTrigger)
    //
    if(datastream_info.wr_event_no == 1)
    {
      fLastEvent = 0;
    }
    else 
    {
      fLastEvent = datastream_info.wr_event_no - 1;
      lastTrigger = frag.getLastTimestampBNB();
      fPrevTrigger->emplace_back(fLastEvent, lastTrigger);
    }
    
    //
    // beam gate
    //
    // TODO UnknownBeamTime: beam gate time is not available yet, here assuming it's same as trigger
    nanoseconds gateStartFromTrigger { 0 }; // beam gate - trigger: hope it's negative...
    auto const gateStart
      = fDetTimings.toSimulationTime(fDetTimings.TriggerTime() + gateStartFromTrigger);
    switch (gate_type) {
      case TriggerGateTypes::BNB:
        fBeamGateInfo->emplace_back
          (gateStart.value(), BNBgateDuration.value(), sim::kBNB);
        break;
      case TriggerGateTypes::NuMI:
        fBeamGateInfo->emplace_back
          (gateStart.value(), NuMIgateDuration.value(), sim::kNuMI);
        break;
      default:
        mf::LogWarning("TriggerDecoder") << "Unsupported gate type #" << gate_type;
    } // switch gate_type
    
    //
    // relative time trigger (raw::Trigger)
    //
    sbn::triggerSource beamGateBit;
    switch (gate_type) {
      case TriggerGateTypes::BNB:  beamGateBit = sbn::triggerSource::BNB;  break;
      case TriggerGateTypes::NuMI: beamGateBit = sbn::triggerSource::NuMI; break;
      default:                     beamGateBit = sbn::triggerSource::Unknown;
    } // switch gate_type
    
    fRelTrigger->emplace_back(
      static_cast<unsigned int>(datastream_info.wr_event_no), // counter
      fDetTimings.TriggerTime().value(),                      // trigger_time
      UnknownBeamTime,                                        // beamgate_time TODO
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
    return;
  }

  std::string_view TriggerDecoder::firstLine
    (std::string const& s, const char* endl /* = "\0\n\r" */)
  {
    return { s.data(), std::min(s.find_first_of(endl), s.size()) };
  }
  
  
  DEFINE_ART_CLASS_TOOL(TriggerDecoder)

}





  
