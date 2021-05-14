////////////////////////////////////////////////
//   
//    File: TriggerDecoder_tool.cc
//       
//    Description: Starting point for extracting ICARUS trigger fragment information into LArSoft object TBD 
//
//    Author: Jacob Zettlemoyer, FNAL
//
///////////////////////////////////////////////

#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// #include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
// #include "cetlib/cpu_timer.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/ExternalTrigger.h" //JCZ: TBD, placeholder for now to represent the idea

#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerUDPFragment.hh"
#include "icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h"
#include "icaruscode/Decode/DecoderTools/TriggerPayloadParser.h"

#include "icaruscode/Decode/DecoderTools/IDecoder.h"

#include <algorithm> // std::replace()
#include <utility> // std::move()
#include <vector>
#include <optional>
#include <stdexcept> // std::logic_error
#include <ctime> // std::ctime()
#include <cstdlib> // std::size_t
#include <memory>

// FIXME in some of the runs, the board reader appears to return nonsense
//       (or the interpretation this side of the online border is wrong);
//       use the local parser to store trigger information?
#define TRIGGERDECODER_TriggerFromDirectParsing 1

// -----------------------------------------------------------------------------
namespace daq 
{
  
  class TriggerDecoder: public IDecoder
  {
  public:
    
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<bool> DiagnosticOutput {
        Name("DiagnosticOutput"),
        Comment("prints additional diagnostic output"),
        false
        };
      
      fhicl::Atom<bool> Debug {
        Name("Debug"),
        Comment("enables debug mode"),
        false
        };
      
    }; // Config
    
    using Parameters = art::ToolConfigTable<Config>;
    
    explicit TriggerDecoder(Parameters const &pset);
    
    virtual void produces(art::ProducesCollector&) override;
    virtual void configure(fhicl::ParameterSet const&) override;
    virtual void initializeDataProducts() override;
    virtual void process_fragment(const artdaq::Fragment &fragment) override;
    virtual void outputDataProducts(art::Event &event) override;
   
  private: 
    using TriggerCollection = std::vector<raw::ExternalTrigger>;
    using TriggerPtr = std::unique_ptr<TriggerCollection>;
    
    TriggerPtr fTrigger;
    
    /// Produces large number of diagnostic messages, use with caution!
    bool const fDiagnosticOutput;
    
    /// Use this for debugging this tool
    bool const fDebug;
    
    /// Information of the previous processed event.
    std::optional<TriggerPayloadParser::TriggerData_t> fLastEventData;
    
    //Add in trigger data member information once it is selected, current LArSoft object likely not enough as is
  }; // class TriggerDecoder

  
  TriggerDecoder::TriggerDecoder(Parameters const &pset)
    : fDiagnosticOutput{ pset().DiagnosticOutput() }
    , fDebug{ pset().Debug() }
  {
  }

  void TriggerDecoder::produces(art::ProducesCollector& collector) 
  {
    collector.produces<TriggerCollection>();
  }

  void TriggerDecoder::configure [[noreturn]] (fhicl::ParameterSet const&) {
    throw
      std::logic_error{ "TriggerDecoder does not support reconfiguration." };
  }


  void TriggerDecoder::initializeDataProducts()
  {
    //use until different object chosen 
    //fTrigger = new raw::Trigger();
    fTrigger = std::make_unique<TriggerPtr::element_type>();
  }

  void TriggerDecoder::process_fragment(const artdaq::Fragment &fragment)
  {
    size_t fragmentID = fragment.fragmentID();
    icarus::ICARUSTriggerUDPFragment frag(fragment);
    //artdaq::Fragment::timestamp_t ts = frag.timestamp();

    int TriggerNumber = 0;
    std::uint64_t GlobalTriggerTimestamp = 0;

#if TRIGGERDECODER_TriggerFromDirectParsing
    std::string payloadAsText
      = reinterpret_cast<char const*>(fragment.dataBeginBytes());
    std::replace(payloadAsText.begin(), payloadAsText.end(), '\x0D', '\n');
    TriggerPayloadParser parser;
    auto const& triggerData = parser(payloadAsText);
    TriggerNumber = triggerData.WR_TS1->eventNo;
    GlobalTriggerTimestamp
      = static_cast<std::uint64_t>(triggerData.WR_TS1->timeStampHigh) * 1'000'000'000UL
      + triggerData.WR_TS1->timeStampLow;
#endif // TRIGGERDECODER_TriggerFromDirectParsing

    if(fDiagnosticOutput) {
#if !TRIGGERDECODER_TriggerFromDirectParsing
      std::string payloadAsText
        = reinterpret_cast<char const*>(fragment.dataBeginBytes());
      std::replace(payloadAsText.begin(), payloadAsText.end(), '\x0D', '\n');
#endif // !TRIGGERDECODER_TriggerFromDirectParsing
      icarus::ICARUSTriggerUDPFragmentMetadata meta = *frag.Metadata();
      mf::LogVerbatim log { "TriggerDecoder" };
      log
          << meta
        << "\nFragment dump:\n"
        << sbndaq::dumpFragment(fragment)
        << "Trigger data as text:\n"
        << " --- BEGIN " << std::string(60, '-') << "\n"
        << payloadAsText
        << "\n --- END --" << std::string(60, '-')
        ;
#if !TRIGGERDECODER_TriggerFromDirectParsing
      TriggerPayloadParser parser;
      auto const& triggerData = parser(payloadAsText);
#endif // !TRIGGERDECODER_TriggerFromDirectParsing
      
      std::time_t dateTime = triggerData.WR_TS1->timeStampHigh;
      std::string dateTimeStr = std::ctime(&dateTime);
      if (dateTimeStr.back() == '\n') dateTimeStr.pop_back();
      mf::LogVerbatim("TriggerDecoder")
        << "Parsed data:"
        << "\nLocal timestamp: "
        << "\n  - event number: " << triggerData.Local_TS1->eventNo
        << "\n  - time:         " << triggerData.Local_TS1->timeStampHigh
          << " s + " << triggerData.Local_TS1->timeStampLow << " ns"
        << "\nWhite Rabbit timestamp: "
        << "\n  - event number: " << triggerData.WR_TS1->eventNo
        << "\n  - time:         " << triggerData.WR_TS1->timeStampHigh
          << " s + " << triggerData.WR_TS1->timeStampLow << " ns ("
          << dateTimeStr << ")"
        ;
      
      if (fLastEventData) { // no previous information
        
        int const dWRevents
          = triggerData.WR_TS1->eventNo - fLastEventData->WR_TS1->eventNo;
        auto const dWRtime
          = triggerData.WR_TS1->nanoseconds() - fLastEventData->WR_TS1->nanoseconds();
        int const dLocalEvents
          = triggerData.Local_TS1->eventNo - fLastEventData->Local_TS1->eventNo;
        auto const dLocalTime
          = triggerData.Local_TS1->nanoseconds() - fLastEventData->Local_TS1->nanoseconds();
        
        mf::LogVerbatim("TriggerDecoder") << std::showpos
          << "Compared with the previous event:"
          << "\n* local timestamp: "
          << "\n  - event number change: " << dLocalEvents
          << "\n  - time change:         " << dLocalTime << " ns"
          << "\n* White Rabbit timestamp: "
          << "\n  - event number change: " << dWRevents
          << "\n  - time change:         " << dWRtime << " ns"
          ;
      }
      fLastEventData = triggerData;
      
    } // if diagnostics
    int local_trig = frag.getName();
    int local_event_no = frag.getEventNo();
    int local_seconds = frag.getSeconds();
    long local_nsec = frag.getNanoSeconds();
    int wr_trig = frag.getWRName();
    int wr_event_no = frag.getWREventNo();
    long wr_seconds = frag.getWRSeconds();
    long wr_nsec = frag.getWRNanoSeconds();
    uint64_t trigger_ts = wr_seconds*1e9 + wr_nsec;
    
#if !TRIGGERDECODER_TriggerFromDirectParsing
    TriggerNumber = wr_event_no;
    GlobalTriggerTimestamp = trigger_ts;
#endif // !TRIGGERDECODER_TriggerFromDirectParsing
    
    fTrigger->emplace_back(TriggerNumber, GlobalTriggerTimestamp);
    
    if(fDiagnosticOutput || fDebug)
    {
      mf::LogInfo log { "TriggerDecoder" };
      
      log 
        <<   "Fragment ID: " << fragmentID
          << " Local Trigger ID: " << local_trig
          << " Event Number: " << local_event_no
          << " Local Seconds: " << local_seconds
          << " Local Nanoseconds: " << local_nsec
        << "\nWhite Rabbit Trigger ID: " << wr_trig
          << " Event Number: " << wr_event_no 
          << " White Rabbit POSIX Seconds: " << std::fixed << wr_seconds
          << " White Rabbit Timestamp Nanoseconds: " << wr_nsec
        << "\nFull Timestamp = " << std::fixed << trigger_ts
        ;
    }
    
    //Once we have full trigger data object, set up and place information into there
    
  } // TriggerDecoder::process_fragment()

  void TriggerDecoder::outputDataProducts(art::Event &event)
  {
    //Place trigger data object into raw data store 
    mf::LogInfo log { "TriggerDecoder" };
    log << "Decoded " << fTrigger->size() << " triggers:";
    for (raw::ExternalTrigger const& trigger: *fTrigger) {
      log << "\n - ID=" << trigger.GetTrigID() << " timestamp=" << trigger.GetTrigTime();
    }
    
    event.put(std::move(fTrigger));
  }

  DEFINE_ART_CLASS_TOOL(TriggerDecoder)

} // namespace daq


