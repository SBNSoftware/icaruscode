/**
 * @file   OpDetWaveformMetaMaker_module.cc
 * @brief  Writes a collection of sbn::OpDetWaveformMeta from PMT waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 22, 2021
 */


// ICARUS libraries
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <memory> // std::make_unique()
#include <utility> // std::move()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::trigger { class OpDetWaveformMetaMaker; }

/**
 * @brief Extracts and saves the time coverage of optical detector waveforms.
 * 
 * This module writes a list of `sbn::OpDetWaveformMeta` objects matching the
 * information of each optical detector waveform.
 * 
 * It may be used as input to modules which require all the information of a
 * PMT waveform except the actual content of the waveform. For such uses,
 * the large waveforms may be dropped and this summary information be kept
 * instead.
 * 
 * 
 * Input data products
 * ====================
 * 
 * This module acts on a _selection_ of tracks, which implies that the input is
 * a set of _pointers_ to tracks rather than to an actual track collection.
 * For each track, an associated time is required.
 * 
 * * `std::vector<raw::OpDetWaveform>` (tag from `Waveforms`):
 *   all optical detector waveforms to extract the information from
 * 
 * 
 * Output data products
 * =====================
 *
 * * `std::vector<sbn::OpDetWaveformMeta>`: a collection parallel to the input
 *   one (from data product configured by `Waveforms`) with the summary
 *   information on each of them; also an explicit courtesy association
 *   `art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>` for uses where
 *   order is not preserved.
 *   The times in these objects are on the same scale as the ones in the source
 *   data product, which is expected to be
 *   @ref DetectorClocksElectronicsStartTime "electronics time scale [us]".
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse online description of the parameters is printed by running
 * `lar --print-description OpDetWaveformMetaMaker`.
 * 
 * * `Waveforms` (input tag, mandatory): the list of optical detector waveforms
 *   to be processed.
 * * `LogCategory` (string, default: `OpDetWaveformMetaMaker`): name of the
 *     output stream category for console messages (managed by MessageFacility
 *     library).
 * 
 * 
 */
class icarus::trigger::OpDetWaveformMetaMaker: public art::SharedProducer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
    fhicl::Atom<art::InputTag> Waveforms {
      Name("Waveforms"),
      Comment("tag of input optical detector waveforms")
      // mandatory
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "OpDetWaveformMetaMaker" // default
      };
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit OpDetWaveformMetaMaker
    (Parameters const& config, art::ProcessingFrame const&);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Fills the plots. Also extracts the information to fill them with.
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fWaveformTag; ///< Input waveforms.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
}; // class icarus::trigger::OpDetWaveformMetaMaker


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& data)
    { return std::make_unique<T>(std::move(data)); }

} // local namespace


//------------------------------------------------------------------------------
icarus::trigger::OpDetWaveformMetaMaker::OpDetWaveformMetaMaker
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedProducer(config)
  // configuration
  , fWaveformTag(config().Waveforms())
  , fLogCategory(config().LogCategory())
{
  
  async<art::InEvent>();
  
  //
  // output data declaration
  //
  produces<std::vector<sbn::OpDetWaveformMeta>>();
  produces<art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>>();
  
  //
  // configuration report (short)
  //
  
  mf::LogInfo{ fLogCategory }
    << "Configuration:"
    << "\n - input waveforms: '" << fWaveformTag.encode() << '\''
    ;
  
} // icarus::trigger::OpDetWaveformMetaMaker::OpDetWaveformMetaMaker()


//------------------------------------------------------------------------------
void icarus::trigger::OpDetWaveformMetaMaker::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  using detinfo::timescales::electronics_time;
  
  //
  // get the timing information for this event
  //
  detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings
    (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event))
    ;
  
  //
  // fetch input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fWaveformTag);

  {
    electronics_time const triggerTime = detTimings.TriggerTime();
    electronics_time const beamGateTime = detTimings.BeamGateTime();
    mf::LogDebug(fLogCategory)
      << "Event " << event.id() << " has beam gate starting at " << beamGateTime
      << " and trigger at " << triggerTime << "."
      << "\nNow extracting information from " << waveformHandle->size()
        << " waveforms."
      ;
  }

  //
  // create the content
  //
  sbn::OpDetWaveformMetaMaker makeOpDetWaveformMeta{ detTimings };
  
  std::vector<sbn::OpDetWaveformMeta> PMTinfo;
  art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform> infoToWaveform;
  
  art::PtrMaker<sbn::OpDetWaveformMeta> const makeInfoPtr { event };
  art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr
    { event, waveformHandle.id() };
  
  for (auto const& [ iWaveform, waveform ]: util::enumerate(*waveformHandle)) {
    assert(iWaveform == PMTinfo.size());
    
    PMTinfo.push_back(makeOpDetWaveformMeta(waveform));
    
    {
      sbn::OpDetWaveformMeta const& info = PMTinfo.back();
      mf::LogTrace log{ fLogCategory };
      log << "Coverage for waveform #" << iWaveform
        << " on channel " << info.channel << ": "
        << info.startTime << " -- " << info.endTime;
      if (info.withTrigger()) log << "; includes trigger";
      if (info.withBeamGate()) log << "; includes beam gate start";
    }
    
    infoToWaveform.addSingle
      (makeInfoPtr(iWaveform), makeWaveformPtr(iWaveform));
    
  } // for
  
  //
  // store output
  //
  event.put(moveToUniquePtr(PMTinfo));
  event.put(moveToUniquePtr(infoToWaveform));
  
} // icarus::trigger::OpDetWaveformMetaMaker::produce()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::OpDetWaveformMetaMaker)


//------------------------------------------------------------------------------
