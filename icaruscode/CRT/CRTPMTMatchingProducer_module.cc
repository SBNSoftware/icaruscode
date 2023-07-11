/**
 * @file   icaruscode/CRT/CRTPMTMatchingProducer_module.cc
 * @brief  Producer to extract and store time matches between CRT hits and PMT.
 * @author Francesco Poppi ( poppi@bo.infn.it )  and Anna Heggestuen (aheggest@colostate.edu)
 * @date   April 26 2023
 */

// icaruscode includes
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h"
//#include "icaruscode/IcarusObj/CRTPMTMatching.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
// Framework includes
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <memory>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <tuple>
#include <utility>

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"


using std::vector;

//namespace icarus::crt {
namespace sbn::crt {
  /**
   * @brief Extracts and saves matches between CRT hits and PMT flashes.
   * 
   * Matching between reconstructed CRT hits and optical flashes is driven by
   * the time difference between them.
   * 
   * Only flashes with a minimum amount (`nOpHitToTrigger`) of contributing PMT
   * above a certain amplitude threshold (`PMTADCThresh`) are considered for
   * matching.
   * Matches are considered successful if they are not farther than a time
   * cutoff (`TimeOfFlightInterval`). Either way, all candidate flashes have
   * matching information saved, even just to say that the matching failed.
   * A single CRT hit may in principle be matched to multiple flashes.
   * In practice this is unlikely since it is rare to have multiple flashes
   * close enough in time to match the same CRT hit time.
   * 
   * This module also works for simulated events, provided that trigger
   * simulation results are available, a proper CRT time offset is specified
   * in the configuration and beam gate limits are also configured.
   * 
   * 
   * ### Timing
   * 
   * The time of a flash used to match it to CRT hits is not taken from the data
   * product, but rather computed anew, as the earliest peak time
   * (`rceob::OpHit::PeakTime()`) among the hits associated to the flash.
   * 
   * The time of a CRT hit from a _data_ sample is taken directly from the
   * timestamp 1 (`sbn::crt::CRTHit::ts1_ns`). For a simulated sample, instead,
   * the absolute timestamp 0 is used (`sbn::crt::CRTHit::ts0()`), compared with
   * a time offset specified in the configuration (`GlobalT0Offset`).
   * 
   * 
   * Configuration parameters
   * -------------------------
   * 
   * All parameters for which no default value is reported are mandatory.
   * 
   * * `OpFlashModuleLabels` (sequence of input tags):
   *   list of reconstructed PMT flash data products to be matched to CRT hits.
   * * `CrtHitModuleLabel` (input tag, default: `crthit`): data product with
   *   the reconstructed CRT hits to be matched to the flashes above.
   * * `TriggerLabel` (input tag, default: `daqTrigger`): data product
   *   containing hardware trigger information. It may be from simulation.
   *   To omit the trigger information, an empty tag can be specified;
   *   if the specified trigger data product is not available, a warning will
   *   be printed on console, but the processing will continue.
   * * `TimeOfFlightInterval` (real number, nanoseconds): the maximum time
   *   between a CRT hit and a flash in order to consider them matched.
   * * `PMTADCThresh` (integral number, PMT ADC counts): optical hits with
   *   amplitude lower than this do not contribute to the contributing hit count
   *   of a flash.
   * * `nOpHitToTrigger` (integral number): minimum number of contributing
   *   reconstructed optical hits to a flash; flashes contributed by fewer hits
   *   will not be matched.
   * * `GlobalT0Offset` (real number, nanoseconds): offset used in simulation,
   *   equivalent to the CRT hit timestamp (TS0) at the time of the reference
   *   trigger.
   * * `MatchBottomCRT` (flag): if unset, the CRT hits from the bottom regions
   *   will be ignored.
   * * `BNBBeamGateMin`, `BNBBeamGateMax`, `NuMIBeamGateMin`, `NuMIBeamGateMax`
   *   (real numbers, nanoseconds): the start and stop of the beam gates,
   *   i.e. the time while the trigger is active for an event. This interval is
   *   used to determine the value `flashInGate` of the matching.
   *   The intervals are specified independently for BNB and NuMI gates
   *   (on-beam and off-beam gates share the same gate interval).
   *   See below for the definition of the time scale.
   * * `BNBinBeamMin`, `BNBinBeamMax`, `NuMIinBeamMin`, `NuMIinBeamMax`
   *   (real numbers, nanoseconds): the start and stop of the beam spill time,
   *   i.e. the time at which neutrinos from the accelerator pass through the
   *   detector. This interval is used to determine the value `flashInBeam` of
   *   the matching. The intervals are specified independently for BNB and NuMI
   *   gates (on-beam and off-beam gates share the same gate interval).
   *   The times are relative to the beam gate opening timestamp; see below for
   *   an explanation of the scale of these times.
   * 
   * 
   * ### Time intervals for the determination of flashes inside beam gates.
   * 
   * The gate interval parameters are on a time scale with as reference the
   * beam gate opening time from the hardware.
   * Ideally, this means that the `BNBBeamGateMin` and `NuMIBeamGateMin`
   * parameters would evaluate to `0`. In practice the beam gate time (and also
   * the trigger time) as stamped by the hardware may not match the actual beam
   * gate opening, and these values need to be aligned with some calibration.
   * 
   * Therefore, effectively `BNBBeamGateMin` represents how much before
   * (negative sign) or after (positive sign) a flash happening exactly at the
   * time of the beam gate opening would be reconstructed with respect to the
   * beam gate timestamp.
   * The best way to determine these values is empirically from a distribution
   * of the time of reconstructed flashes in a sample of majority-triggered
   * events, which should appear like a continuum with a step due to the opening
   * of the beam gate (which is biassed because of the presence of the trigger)
   * and a further step due to the presence of the neutrinos.
   * The first plateau defines the beam gate interval (`XxxxBeamGateXxx`), the
   * second one defined the spill interval (`XxxxinBeamXxx`).
   * 
   * 
   * Input
   * ------
   * 
   * * `std::vector<recob::OpFlash>` (all tags from `OpFlashModuleLabels`):
   *   the collections of flashes to match. All the flashes are treated the
   *   same, independently of which of the collections they come from.
   * * `std::vector<sbn::crt::CRTHit>` (tag from `CrtHitModuleLabel`): the CRT
   *   hits to match to the flashes.
   * * `sbn::ExtraTriggerInfo` (tag from `TriggerLabel`): if available,
   *   information whether the flash is at beam time and in the beam gate time
   *   are saved. It is otherwise optional.
   * 
   * 
   * Output
   * -------
   * 
   * * `std::vector<sbn::crt::CRTPMTMatching>`: an entry for each matched
   *   flash. Entries are sorted by flast time, from the earliest to the latest,
   *   irregardless of which of the input flash collections they come from.
   *   Each entry contains information of all the matched CRT hits and the type
   *   of the matching:
   *    * `flashID`: an ID made up out of the flash content
   *         (`sbn::crt::makeFlashID()`).
   *    * `flashTime`: from `recob::OpFlash::Time()`.
   *    * `flashGateTime`: time of the flash from the beam gate opening.
   *    * `firstOpHitPeakTime`: first `recob::OpHit::PeakTime()` in the flash.
   *    * `firstOpHitStartTime`: from `recob::OpHit::StartTime()` in the flash;
   *         only filled if present in the hit (`HasStartTime()`).
   *    * `flashInGate`: whether the flash is in the beam gate interval as
   *       configured via `BNBinBeamMin`/`BNBinBeamMax` or the corresponding
   *       settings for NuMI beam.
   *    * `flashInBeam`: whether the flash is in the beam gate interval as
   *       configured via `BNBBeamGateMin`/`BNBBeamGateMax` or the corresponding
   *       settings for NuMI beam.
   *    * `flashPE`: `recob::OpFlash::TotalPE()` saved in the flash.
   *    * `flashPosition`: centroid of the flash, recomputed as the average of
   *         the position of the contributing PMTs, weighted by the amplitude of
   *         their hits.
   *    * `flashYWidth`: from `recob::OpFlash::YWidth()`.
   *    * `flashZWidth`: from `recob::OpFlash::ZWidth()`.
   *    * `flashClassification`: the topology of this track, according to the
   *      categories in `sbn::crt::MatchType`; i.e. if it is a particle that
   *      seems to be entering from the top, or exiting from within, etc.
   *    * `matchedCRTHits` (`sbn::crt::MatchedCRT`): for each CRT hit matched
   *      to the flash, information on its location and time relative to the
   *      flash.
   *        *`pos`: the location of the hit in space.
   *        * `PMTTimeDiff` [&micro;s]: the time of flight computed using
   *          the flash time as described above (not `recob::OpFlash::Time()`);
   *          it is also the time on which the matching decision was taken,
   *          and the time that settles the relative time of hit and flash.
   *        * `time` [&micro;s]: the time of the matched CRT hit.
   *        * `region`: the number of CRT region where the matched hit is.
   *        * `sys`: which subdetector the hit is in; `0` for top CRT, `1`
   *          for side CRT.
   *    * `nTopCRTHitsBefore`: not saved yet (left to default value).
   *    * `nTopCRTHitsAfter`: not saved yet (left to default value).
   *    * `nSideCRTHitsBefore`: not saved yet (left to default value).
   *    * `nSideCRTHitsAfter`: not saved yet (left to default value).
   * * `art::Assns<sbn::crt::CRTPMTMatching, recob::OpFlash>`:
   *   associations linking the matched flash and the match information;
   *   this is a one-to-one association; all `sbn::crt::CRTPMTMatching`
   *   objects in the association come from the collection data product
   *   documented above, and all have exactly one associated flash. Flashes
   *   that were not matched are not present in this association.
   * * `art::Assns<sbn::crt::CRTPMTMatching, sbn::crt::CRTHit>`:
   *   associations linking the matched CRT hits and the match information;
   *   this is a one-to-many association; all `sbn::crt::CRTPMTMatching`
   *   objects in the association come from the collection data product
   *   documented above, and they may have one or more CRT hits associated to
   *   them. CRT hits that were not matched are not present in this association.
   *   Also failed matches, with no CRT hit, are not present.
   * * `art::Assns<recob::OpFlash, sbn::crt::CRTHit, sbn::crt::CRTPMTMatchingInfo>`:
   *   direct association between a flash and its matched CRT hit (or hits);
   *   this is a one-to-many association. This information is redundant with
   *   the other two associations produced by this module. Metadata includes:
   *     * `timeOfFlight`: time of flight between hit and flash [&micro;s]
   *     * `direction`, derived from the time of flight, of the flash compared
   *       to the hit.
   *     * `distance` between the reconstructed CRT hit position and the
   *       centroid of the reconstructed flash [cm]
   * 
   * 
   * Services
   * ---------
   * 
   * * `Geometry` for the determination of the position of the PMT
   * 
   * 
   * Multithreading
   * ---------------
   * 
   * The module currently supports _art_ multithreading, i.e. the processing of
   * multiple events in the same run (and subrun) at the same time.
   * It does not directly use any multithreading within its code.
   * 
   */
  class CRTPMTMatchingProducer : public art::SharedProducer {
  public:
 
    using CRTHit = sbn::crt::CRTHit;
    using CRTPMTMatching = sbn::crt::CRTPMTMatching;
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Sequence<art::InputTag> OpFlashModuleLabels{
        Name{ "OpFlashModuleLabels" },
        Comment{ "reconstructed PMT flash data products to be matched to CRT hits" }
        };
      
      fhicl::Atom<art::InputTag> CrtHitModuleLabel{
        Name{ "CrtHitModuleLabel" },
        Comment{ "reconstructed CRT hits to be matched to the flashes" },
        "crthit"
        };
      
      fhicl::Atom<art::InputTag> TriggerLabel{
        Name{ "TriggerLabel" },
        Comment{ "data product containing hardware trigger information" },
        "daqTrigger"
        };
      
      fhicl::Atom<double> TimeOfFlightInterval{
        Name{ "TimeOfFlightInterval" },
        Comment{ "maximum time between a CRT hit and a flash"
          " in order to consider them matched [ns]" }
        };
      
      fhicl::Atom<int> PMTADCThresh{
        Name{ "PMTADCThresh" },
        Comment{ "optical hits with amplitude lower than this do not contribute"
          " to the contributing hit count of a flash [ADC]" }
        };
      
      fhicl::Atom<int> nOpHitToTrigger{
        Name{ "nOpHitToTrigger" },
        Comment{ "minimum number of contributing reconstructed optical hits to a flash" }
        };
      
      fhicl::Atom<double> GlobalT0Offset{
        Name{ "GlobalT0Offset" },
        Comment{ "offset used in simulation, equivalent to the CRT hit timestamp (TS0)"
          " at the time of the reference trigger [ns]" }
        };
      fhicl::Atom<bool> MatchBottomCRT{
        Name{ "MatchBottomCRT" },
        Comment{ "bool fcl parameter to perform CRT-PMT Matching on bottom CRT data" }
        };
      
      fhicl::Atom<double> BNBBeamGateMin{
        Name{ "BNBBeamGateMin" },
        Comment{ "start of the BNB beam gate w.r.t. beam gate timestamp [ns]" }
        };
      
      fhicl::Atom<double> BNBBeamGateMax{
        Name{ "BNBBeamGateMax" },
        Comment{ "end of the BNB beam gate w.r.t. beam gate timestamp [ns]" }
        };
      
      fhicl::Atom<double> BNBinBeamMin{
        Name{ "BNBinBeamMin" },
        Comment{ "start of the BNB spill w.r.t. beam gate timestamp [ns]" }
        };
      
      fhicl::Atom<double> BNBinBeamMax{
        Name{ "BNBinBeamMax" },
        Comment{ "end of the BNB spill w.r.t. beam gate timestamp [ns]" }
        };
      
      fhicl::Atom<double> NuMIBeamGateMin{
        Name{ "NuMIBeamGateMin" },
        Comment{ "start of the NuMI beam gate w.r.t. beam gate timestamp [ns]" }
        };
      
      fhicl::Atom<double> NuMIBeamGateMax{
        Name{ "NuMIBeamGateMax" },
        Comment{ "end of the NuMI beam gate w.r.t. beam gate timestamp [ns]" }
        };
      
      fhicl::Atom<double> NuMIinBeamMin{
        Name{ "NuMIinBeamMin" },
        Comment{ "start of the NuMI spill w.r.t. beam gate timestamp [ns]" }
        };
      
      fhicl::Atom<double> NuMIinBeamMax{
        Name{ "NuMIinBeamMax" },
        Comment{ "end of the NuMI spill w.r.t. beam gate timestamp [ns]" }
        };
      
    };
    
    using Parameters = art::SharedProducer::Table<Config>;
    
    explicit CRTPMTMatchingProducer
      (Parameters const & p, art::ProcessingFrame const&);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Required functions.
    void produce(art::Event & e, art::ProcessingFrame const&) override;

  private:

    // Params from fcl file.......

    std::vector<art::InputTag> const fFlashLabels;   ///< All optical flash input tags.
    art::InputTag const fCrtHitModuleLabel;          ///< Name of CRT producer.
    art::InputTag const fTriggerLabel;               ///< Name of trigger producer.
    
    double const fTimeOfFlightInterval;              ///< CRT-PMT time difference interval to find the match.
    int const fPMTADCThresh;                         ///< ADC amplitude for a PMT to be considered above threshold.
    int const fnOpHitToTrigger;                      ///< Number of OpHit above threshold to mimic the triggering PMT.
    double const fGlobalT0Offset;                    ///< 1.6 ms delay to shift CRT Hit T0, the CRT Timing variable we use in MC.
    bool const fMatchBottomCRT;                      ///< Bool fcl value to perform matching on bottom CRT Hits if true, currently set to false.
    double const fBNBBeamGateMin;
    double const fBNBBeamGateMax;
    double const fBNBinBeamMin;
    double const fBNBinBeamMax;
    double const fNuMIBeamGateMin;
    double const fNuMIBeamGateMax;
    double const fNuMIinBeamMin;
    double const fNuMIinBeamMax;

    geo::GeometryCore const* const fGeometryService;  ///< Pointer to Geometry provider.
    
    /**
     * @brief Returns a `MatchedCRT` out of the provided information.
     * @param hit the matched CRT hit
     * @param tflash the time of the matched flash [us]
     * @param isRealData whether the data is real (as opposed to Monte Carlo)
     * @return a complete `sbn::crt::MatchedCRT` record
     */
    sbn::crt::MatchedCRT makeMatchedCRT
      (sbn::crt::CRTHit const& hit, double tflash, bool isRealData) const;
    
  }; // class CRTPMTMatchingProducer

  CRTPMTMatchingProducer::CRTPMTMatchingProducer
    (Parameters const & p, art::ProcessingFrame const&)
  : SharedProducer{p},
    fFlashLabels{p().OpFlashModuleLabels()},
    fCrtHitModuleLabel{p().CrtHitModuleLabel()},
    fTriggerLabel{p().TriggerLabel()},
    fTimeOfFlightInterval{p().TimeOfFlightInterval()},
    fPMTADCThresh{p().PMTADCThresh()},
    fnOpHitToTrigger{p().nOpHitToTrigger()},
    fGlobalT0Offset{p().GlobalT0Offset()},
    fMatchBottomCRT{p().MatchBottomCRT()},
    fBNBBeamGateMin{p().BNBBeamGateMin()},
    fBNBBeamGateMax{p().BNBBeamGateMax()},
    fBNBinBeamMin{p().BNBinBeamMin()},
    fBNBinBeamMax{p().BNBinBeamMax()},
    fNuMIBeamGateMin{p().NuMIBeamGateMin()},
    fNuMIBeamGateMax{p().NuMIBeamGateMax()},
    fNuMIinBeamMin{p().NuMIinBeamMin()},
    fNuMIinBeamMax{p().NuMIinBeamMax()},
    fGeometryService{lar::providerFrom<geo::Geometry>()}
  {
    async<art::InEvent>();
    
    produces< std::vector<CRTPMTMatching> >();
    
    produces< art::Assns<sbn::crt::CRTPMTMatching, recob::OpFlash> >();
    produces< art::Assns<sbn::crt::CRTPMTMatching, sbn::crt::CRTHit> >();
    produces< art::Assns<recob::OpFlash, sbn::crt::CRTHit, CRTPMTMatchingInfo> >();
    
  } // CRTPMTMatchingProducer()


  void CRTPMTMatchingProducer::produce(art::Event & e, art::ProcessingFrame const&)
  {
    mf::LogDebug("CRTPMTMatchingProducer") << "beginning CRTPMTProducer";
    // Trigger data product variables
    // add trigger info
    
    sbn::ExtraTriggerInfo const* trigInfo = nullptr;
    if (!fTriggerLabel.empty()) {
      if (auto const trigger_handle = e.getHandle<sbn::ExtraTriggerInfo>(fTriggerLabel))
      {
        trigInfo = trigger_handle.product();
      } else {
        mf::LogDebug("CRTPMTMatchingProducer") << "No sbn::ExtraTriggerInfo associated to label: " << fTriggerLabel.encode() << "\n";
      }
    }
    else {
      mf::LogDebug("CRTPMTMatchingProducer") << "No TriggerLabel in : " << fTriggerLabel.encode() << "\n";
    }
    sbn::triggerSource const gateType
      = trigInfo? trigInfo->sourceType: sbn::triggerSource::Unknown;
    int64_t const triggerGateDiff
      = trigInfo? trigInfo->triggerFromBeamGate(): 0;
    
    double BeamGateMin = 0.0, BeamGateMax = 0.0;
    double inBeamMin = 0.0, inBeamMax = 0.0;
    switch (gateType) {
      case sbn::triggerSource::BNB:
      case sbn::triggerSource::OffbeamBNB:
        BeamGateMin = fBNBBeamGateMin;
        BeamGateMax = fBNBBeamGateMax;
        inBeamMin = fBNBinBeamMin;
        inBeamMax = fBNBinBeamMax;
        break;
      case sbn::triggerSource::NuMI:
      case sbn::triggerSource::OffbeamNuMI:
        BeamGateMin = fNuMIBeamGateMin;
        BeamGateMax = fNuMIBeamGateMax;
        inBeamMin = fNuMIinBeamMin;
        inBeamMax = fNuMIinBeamMax;
        break;
      case sbn::triggerSource::Unknown:
        mf::LogWarning("CRTPMTMatchingProducer")
          << "Unknown beam gate type: in-gate flags won't be available.";
        break;
      default:
        mf::LogWarning("CRTPMTMatchingProducer")
          << "Unsupported beam gate type: " << name(gateType)
          << "; in-gate flags won't be available.";
        break;
    } // switch
    
    auto CRTPMTMatchesColl = std::make_unique<std::vector<CRTPMTMatching>>();

    art::PtrMaker<CRTPMTMatching> const makeInfoPtr(e);
    auto FlashAssociation = std::make_unique<art::Assns<CRTPMTMatching, recob::OpFlash>>();
    auto CRTAssociation = std::make_unique<art::Assns<CRTPMTMatching, sbn::crt::CRTHit>>();
    auto FlashCRTAssociation = std::make_unique<art::Assns<recob::OpFlash, sbn::crt::CRTHit, CRTPMTMatchingInfo>>();
    
    // add CRTHits
    std::vector<art::Ptr<CRTHit>> crtHitList;
    if (auto crtHitListHandle = e.getHandle<std::vector<CRTHit>>(fCrtHitModuleLabel))
      art::fill_ptr_vector(crtHitList, crtHitListHandle);
    {
      mf::LogTrace log("CRTPMTMatchingProducer");
      log << fCrtHitModuleLabel.encode() << " has " << crtHitList.size() << " CRT hits";
      for (art::Ptr<sbn::crt::CRTHit> const& hitPtr: crtHitList) {
        log << "\n [#" << hitPtr.key() << "] at T0=" << hitPtr->ts0() << " ns, T1="
          << hitPtr->ts1()/1000.0 << " us, (" << hitPtr->x_pos << ", "
          << hitPtr->y_pos << ", " << hitPtr->z_pos << ") cm, " << hitPtr->peshit
          << " p.e., on " << hitPtr->tagger;
      } // for
    }
    bool const isRealData = e.isRealData();
    int n_entering_matches = 0;
    int n_exiting_matches = 0;
    mf::LogTrace("CRTPMTMatchingProducer") << "is this real data? " << std::boolalpha << isRealData;
    
    //
    // prepare the flashes (sorted by time)
    //
    std::vector
      <std::tuple<art::Ptr<recob::OpFlash>, std::vector<recob::OpHit const*>>>
      flashesAndHits;
    for (art::InputTag const& flashLabel : fFlashLabels) {
      auto const flashHandle =
        e.getHandle<std::vector<recob::OpFlash>>(flashLabel);
      
      mf::LogTrace("CRTPMTMatchingProducer")
        << "Matching " << flashHandle->size() << " flashes from " << flashLabel.encode();
      
      art::FindMany<recob::OpHit> const findManyHits(flashHandle, e, flashLabel);
      for (std::size_t const iFlash: util::counter(flashHandle->size())) {
        flashesAndHits.emplace_back(
          art::Ptr<recob::OpFlash>{ flashHandle, iFlash },
          findManyHits.at(iFlash)
          );
      } // for flashes
    } // for flash labels
    
    // sort the flashes by Time(); their associated hits follow
    std::sort(flashesAndHits.begin(), flashesAndHits.end(),
      [](auto const& flashAndHitsA, auto const& flashAndHitsB)
        {
          art::Ptr<recob::OpFlash> const& flashPtrA = std::get<0U>(flashAndHitsA);
          art::Ptr<recob::OpFlash> const& flashPtrB = std::get<0U>(flashAndHitsB);
          if (flashPtrA == flashPtrB) return false; // or we would emit warnings
          if (!flashPtrA || !flashPtrB) return flashPtrA < flashPtrB;
          if (flashPtrA->Time() != flashPtrB->Time())
            return flashPtrA->Time() < flashPtrB->Time();
          // this should be very, very rare
          mf::LogPrint("CRTPMTMatchingProducer")
            << "Flashes " << flashPtrA << " and " << flashPtrB
            << " both are at time " << flashPtrA->Time() << " us";
          if (flashPtrA->TotalPE() != flashPtrB->TotalPE())
            return flashPtrA->TotalPE() > flashPtrB->TotalPE();
          mf::LogPrint("CRTPMTMatchingProducer")
            << "  ... and they both have " << flashPtrA->TotalPE() << "?!?";
          return flashPtrA < flashPtrB;
        }
      );
    
    // 
    // process all the flashes
    // 
    for (auto const& [ flashPtr, hits ]: flashesAndHits) {

      double const tflash = flashPtr->Time();
      double const totPe  = flashPtr->TotalPE();
      double const flashYWidth = flashPtr->YWidth();
      double const flashZWidth = flashPtr->ZWidth();

      int nPMTsTriggering = 0;
      double firstOpHitPeakTime = 999999;
      double firstOpHitStartTime = 999999;

      geo::vect::MiddlePointAccumulator flashCentroid;
      for (auto const& hit : hits) {
        if (hit->Amplitude() > fPMTADCThresh) nPMTsTriggering++;
        if (firstOpHitPeakTime > hit->PeakTime()) firstOpHitPeakTime = hit->PeakTime();
        if(hit->HasStartTime()){
          if (firstOpHitStartTime > hit->StartTime()) firstOpHitStartTime = hit->StartTime(); 
          //mf::LogTrace("CRTPMTMatchingProducer") << "OpHit has a starttime " << hit->StartTime() << ", saving " <<  firstOpHitStartTime << " to firstOpHitStartTime.\n";
        }

        geo::Point_t const pos =
          fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel())
          .GetCenter();
        double const amp = hit->Amplitude();
        flashCentroid.add(pos, amp);
      }
      geo::Point_t const flash_pos = flashCentroid.middlePoint();
      {
        mf::LogTrace log("CRTPMTMatchingProducer");
        log << "Now matching flash #" << flashPtr.key()
          << " at " << flashPtr->Time() << " us and (";
        if (flashPtr->hasXCenter()) log << flashPtr->XCenter();
        else log << "n/a";
        log << ", " << flashPtr->YCenter() << ", " << flashPtr->ZCenter()
          << ") cm [" << hits.size() << " op.hits] -> first time: "
          << firstOpHitPeakTime << " us, centroid: " << flash_pos << " cm";
      }
      
      if (nPMTsTriggering < fnOpHitToTrigger) {
        mf::LogTrace("CRTPMTMatchingProducer")
          << "  => skipped (only " << nPMTsTriggering << " < " << fnOpHitToTrigger
          << " hits above threshold)";
        continue;
      }
      
      double const thisRelGateTime = triggerGateDiff + tflash * 1e3; // ns
      bool const thisInTime_gate
        = thisRelGateTime > BeamGateMin && thisRelGateTime < BeamGateMax;
      bool const thisInTime_beam
        = thisRelGateTime > inBeamMin && thisRelGateTime < inBeamMax;
      
      icarus::crt::CRTMatches const crtMatches = icarus::crt::CRTHitmatched(
        firstOpHitPeakTime, flash_pos, crtHitList, fTimeOfFlightInterval, isRealData, fGlobalT0Offset, fMatchBottomCRT);
      
      std::vector<MatchedCRT> thisFlashCRTmatches;
        std::vector<art::Ptr<sbn::crt::CRTHit>> CRTPtrs; // same order as thisFlashCRTmatches
      MatchType const eventType = crtMatches.flashType;
      if (!crtMatches.entering.empty()) {
        mf::LogTrace("CRTPMTMatchingProducer")
          << "Entering matches (" << crtMatches.entering.size() << "):";
        for (auto const& entering : crtMatches.entering) {
          n_entering_matches++;
          CRTPtrs.push_back(entering.CRTHit);
          thisFlashCRTmatches.push_back(makeMatchedCRT(*entering.CRTHit, tflash, isRealData));
        }
      }
      if (!crtMatches.exiting.empty()) {
        mf::LogTrace("CRTPMTMatchingProducer")
          << "Exiting matches (" << crtMatches.exiting.size() << "):";
        for (auto const& exiting : crtMatches.exiting) {
          n_exiting_matches++;
          CRTPtrs.push_back(exiting.CRTHit);
          thisFlashCRTmatches.push_back(makeMatchedCRT(*exiting.CRTHit, tflash, isRealData));
        }
      }
      if (!thisFlashCRTmatches.empty() ) {
        mf::LogTrace("CRTPMTMatchingProducer") << "pushing back flash with "
                                                << thisFlashCRTmatches.size() << " CRT Matches. --> Match classification = " << std::to_string(static_cast<int>(eventType))<< "\n" ;
        mf::LogTrace("CRTPMTMatchingProducer") << "\tfirstOpHitPeakTime " << firstOpHitPeakTime << ", firstOpHitStartTime = " << firstOpHitStartTime << " (us)\n\ttotPe = " << totPe << "\n\tflashYWidth = " << flashYWidth << ", flashZWidth = " << flashZWidth << " (cm)\n";

      }
      
      int const flashID = icarus::crt::makeFlashID(flashPtr.get());
      mf::LogTrace("CRTPMTMatchingProducer") << "saving Flash with FlashID = " << flashID << "\n";
      FlashType thisFlashType = { /* .flashID = */ flashID, // C++20: restore initializers
                                  /* .flashPos = */ flash_pos,
                                  /* .flashTime = */ tflash,
                                  /* .flashGateTime = */ thisRelGateTime / 1000.0, // -> us
                                  /* .firstOpHitPeakTime = */ firstOpHitPeakTime,
                                  /* .firstOpHitStartTime = */ firstOpHitStartTime,
                                  /* .inBeam = */ thisInTime_beam,
                                  /* .inGate = */ thisInTime_gate,
                                  /* .flashPE = */ totPe, 
                                  /* /.flashYWidth = */ flashYWidth,
                                  /* /.flashZWidth = */ flashZWidth,
                                  /* .classification = */ eventType,
                                  /* .CRTmatches = */ std::move(thisFlashCRTmatches)};
      
      // add to the data products
      
      // pointer to the matching info object we are going to add to the collection:
      art::Ptr<CRTPMTMatching> const infoPtr
        = makeInfoPtr(CRTPMTMatchesColl->size());
      
      FlashAssociation->addSingle(infoPtr, flashPtr);
      
      for (icarus::crt::CRTPMT const& CRTmatch: crtMatches.entering) {
        CRTAssociation->addSingle(infoPtr, CRTmatch.CRTHit);
        FlashCRTAssociation->addSingle(
          flashPtr, CRTmatch.CRTHit,
          CRTPMTMatchingInfo{
            /* .direction =    */ CRTPMTMatchingInfo::Dir::entering, // C++20: restore initializers
            /* .timeOfFlight = */ CRTmatch.tof,
            /* .distance =     */ CRTmatch.distance
          }
          );
      }
      for (icarus::crt::CRTPMT const& CRTmatch: crtMatches.exiting) {
        CRTAssociation->addSingle(infoPtr, CRTmatch.CRTHit);
        FlashCRTAssociation->addSingle(
          flashPtr, CRTmatch.CRTHit,
          CRTPMTMatchingInfo{
            /* .direction =    */ CRTPMTMatchingInfo::Dir::exiting,  // C++20: restore initializers
            /* .timeOfFlight = */ CRTmatch.tof,
            /* .distance =     */ CRTmatch.distance
          }
        );
      }
      
      CRTPMTMatching matchInfo = sbn::crt::FillCRTPMT(thisFlashType);
      CRTPMTMatchesColl->push_back(std::move(matchInfo));
      
    } // for flash data products
    mf::LogTrace("CRTPMTMatchingProducer") << "in total " << n_entering_matches << " entering CRT Hit matches, " << n_exiting_matches << " exiting CRT Hit matches\n";
    mf::LogTrace("CRTPMTMatchingProducer") <<"This Event has "<<CRTPMTMatchesColl->size()<<"  Flashes candidate for CRT matching."<<std::endl;
    e.put(std::move(CRTPMTMatchesColl));
    e.put(std::move(FlashAssociation));
    e.put(std::move(CRTAssociation));
    e.put(std::move(FlashCRTAssociation));

  } // CRTPMTMatchingProducer::produce()


  // ---------------------------------------------------------------------------
  sbn::crt::MatchedCRT CRTPMTMatchingProducer::makeMatchedCRT
    (sbn::crt::CRTHit const& hit, double tflash, bool isRealData) const
  {
    geo::Point_t const thisCRTpos {hit.x_pos, hit.y_pos, hit.z_pos};
    // fGlobalT0Offset = 1.6e6 ns = 1600000 ns.
    double const CRTtime  // ns -> us
      = icarus::crt::CRTHitTime(hit, fGlobalT0Offset, isRealData) / 1000.0;
    double const CRTTof_opflash = CRTtime - tflash; // us
    int const CRTRegion = hit.plane;
    // Very lazy way to determine if the Hit is a Top or a Side, bottom CRT is only set if fMatchBottomCRT true; 
    int const CRTSys = (fMatchBottomCRT && CRTRegion >= 49) ? 2 : ((CRTRegion >= 36)? 1: 0);

    mf::LogTrace("CRTPMTMatchingProducer")
      << "  Match: tof = crtTime - tflash  = " << CRTtime
      << " - "<< tflash << " = " << CRTTof_opflash << " (us) in region " << CRTRegion << "";
    return { /* .position =*/  thisCRTpos, // C++20: restore initializers
             /* .PMTTimeDiff = */ CRTTof_opflash,
             /* .time = */ CRTtime,
             /* .sys = */ CRTSys,
             /* .region = */ CRTRegion};
  }


  DEFINE_ART_MODULE(CRTPMTMatchingProducer)

} //end namespace
