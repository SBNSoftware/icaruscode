/**
 * @file   icaruscode/CRT/CRTPMTMatchingProducer_module.cc
 * @brief  Producer to extract and store time matches between CRT hits and PMT.
 * @author Francesco Poppi ( poppi@bo.infn.it )  and Anna Heggestuen 
 * @date   April 26 2023
 */

// icaruscode includes
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/CRT/CRTUtils/CRTPMTMatchingUtils.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <memory>
#include <iomanip>
#include <vector>
#include <utility>

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"


using std::vector;

namespace icarus::crt {

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
   * * `OpFlashModuleLabels` (sequence of input tags, mandatory):
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
   * * `std::vector<icarus::crt::CRTPMTMatching>`: an entry for each matched
   *   flash; the entry contains information of all the matched CRT hits and
   *   the type of the matching:
   *    * `flashID`: not saved yet (set to `0`).
   *    * `flashTime`: from `recob::OpFlash::Time()`.
   *    * `flashGateTime`: time of the flash from the beam gate opening.
   *    * `firstOpHitPeakTime`: not saved yet (set to `0`).
   *    * `firstOpHitStartTime`: not saved yet (set to `0`).
   *    * `flashInGate`: whether the flash is in the beam gate interval as
   *       configured via `BNBinBeamMin`/`BNBinBeamMax` or the corresponding
   *       settings for NuMI beam.
   *    * `flashInBeam`: whether the flash is in the beam gate interval as
   *       configured via `BNBBeamGateMin`/`BNBBeamGateMax` or the corresponding
   *       settings for NuMI beam.
   *    * `flashAmplitude_pe`: not saved yet (set to `0`).
   *    * `flashPosition`: centroid of the flash, recomputed as the average of
   *         the position of the contributing PMTs, weighted by the amplitude of
   *         their hits.
   *    * `flashYWidth`: not saved yet (set to `0`).
   *    * `flashZWidth`: not saved yet (set to `0`).
   *    * `flashClassification`: the topology of this track, according to the
   *      categories in `icarus::crt::MatchType`; i.e. if it is a particle that
   *      seems to be entering from the top, or exiting from within, etc.
   *    * `matchedCRTHits` (`icarus::crt::MatchedCRT`): for each CRT hit matched
   *      to the flash, information on its location and time relative to the
   *      flash.
   *        *`CRTHitPos`: the location of the hit in space.
   *        * `CRTPMTTimeDiff` [&micro;s]: the time of flight computed using
   *          the flash time as described above (not `recob::OpFlash::Time()`);
   *          it is also the time on which the matching decision was taken,
   *          and the time that settles the relative time of hit and flash.
   *        * `CRTTime` [&micro;s]: the time of the matched CRT hit.
   *        * `CRTRegion`: the number of CRT region where the matched hit is.
   *        * `CRTSys`: which subdetector the hit is in; `0` for top CRT, `1`
   *          for side CRT.
   *    * `topCRTBefore`: not saved yet (set to `0`).
   *    * `topCRTAfter`: not saved yet (set to `0`).
   *    * `sideCRTBefore`: not saved yet (set to `0`).
   *    * `sideCRTAfter`: not saved yet (set to `0`).
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
   * The module currently does not support multithreading, even if it is ready
   * to.
   * 
   * 
   */
  class CRTPMTMatchingProducer : public art::EDProducer {
  public:
 
    using CRTHit = sbn::crt::CRTHit;
    
    explicit CRTPMTMatchingProducer(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Required functions.
    void produce(art::Event & e) override;

  private:

    // Params from fcl file.......

    std::vector<art::InputTag> const fFlashLabels;   ///< All optical flash input tags.
    art::InputTag const fCrtHitModuleLabel;          ///< Name of CRT producer.
    art::InputTag const fTriggerLabel;               ///< Name of trigger producer.
    
    double const fTimeOfFlightInterval;              ///< CRT-PMT time difference interval to find the match.
    int const fPMTADCThresh;                         ///< ADC amplitude for a PMT to be considered above threshold.
    int const fnOpHitToTrigger;                      ///< Number of OpHit above threshold to mimic the triggering PMT.
    double const fGlobalT0Offset;                    ///< 1.6 ms delay to shift CRT Hit T0, the CRT Timing variable we use in MC.

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
     * @return a complete `icarus::crt::MatchedCRT` record
     */
    icarus::crt::MatchedCRT makeMatchedCRT
      (sbn::crt::CRTHit const& hit, double tflash, bool isRealData) const;
    

  }; // class CRTPMTMatchingProducer

  CRTPMTMatchingProducer::CRTPMTMatchingProducer(fhicl::ParameterSet const & p)
  : EDProducer{p},
    fFlashLabels(p.get<std::vector<art::InputTag>>("OpFlashModuleLabels")),
    fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit")),
    fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")),
    fTimeOfFlightInterval(p.get<double>("TimeOfFlightInterval")),
    fPMTADCThresh(p.get<int>("PMTADCThresh")),
    fnOpHitToTrigger(p.get<int>("nOpHitToTrigger")),
    fGlobalT0Offset(p.get<int>("GlobalT0Offset")),
    fBNBBeamGateMin(p.get<double>("BNBBeamGateMin")),
    fBNBBeamGateMax(p.get<double>("BNBBeamGateMax")),
    fBNBinBeamMin(p.get<double>("BNBinBeamMin")),
    fBNBinBeamMax(p.get<double>("BNBinBeamMax")),
    fNuMIBeamGateMin(p.get<double>("NuMIBeamGateMin")),
    fNuMIBeamGateMax(p.get<double>("NuMIBeamGateMax")),
    fNuMIinBeamMin(p.get<double>("NuMIinBeamMin")),
    fNuMIinBeamMax(p.get<double>("NuMIinBeamMax")),
    fGeometryService(lar::providerFrom<geo::Geometry>())
  {
    produces< std::vector<CRTPMTMatching> >();
  } // CRTPMTMatchingProducer()


  void CRTPMTMatchingProducer::produce(art::Event & e)
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

    //auto FlashAssociation = std::make_unique<art::Assns<CRTPMTMatching, recob::OpFlash>>();
    
    // add CRTHits
    std::vector<art::Ptr<CRTHit>> crtHitList;
    if (auto crtHitListHandle = e.getHandle<std::vector<CRTHit>>(fCrtHitModuleLabel))
      art::fill_ptr_vector(crtHitList, crtHitListHandle);
    bool const isRealData = e.isRealData();
    mf::LogTrace("CRTPMTMatchingProducer") << "is this real data? " << std::boolalpha << isRealData;
    // add optical flashes
    for (art::InputTag const& flashLabel : fFlashLabels) {
    auto const flashHandle =
      e.getHandle<std::vector<recob::OpFlash>>(flashLabel);
    art::FindMany<recob::OpHit> const findManyHits(flashHandle, e, flashLabel);

    std::vector<FlashType> thisEventFlashes;

    for (auto const& [iflash, flash] : util::enumerate(*flashHandle)) {
      double const tflash = flash.Time();
      vector<recob::OpHit const*> const& hits = findManyHits.at(iflash);
      int nPMTsTriggering = 0;
      double firstTime = 999999;
      geo::vect::MiddlePointAccumulator flashCentroid;
      for (auto const& hit : hits) {
        if (hit->Amplitude() > fPMTADCThresh) nPMTsTriggering++;
        if (firstTime > hit->PeakTime()) firstTime = hit->PeakTime();
        geo::Point_t const pos =
          fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel())
          .GetCenter();
        double const amp = hit->Amplitude();
        flashCentroid.add(pos, amp);
      }
      geo::Point_t const flash_pos = flashCentroid.middlePoint();
      if (nPMTsTriggering < fnOpHitToTrigger) continue;

      double const thisRelGateTime = triggerGateDiff + tflash * 1e3; // ns
      bool const thisInTime_gate
        = thisRelGateTime > BeamGateMin && thisRelGateTime < BeamGateMax;
      bool const thisInTime_beam
        = thisRelGateTime > inBeamMin && thisRelGateTime < inBeamMax;
      
      icarus::crt::CRTMatches const crtMatches = CRTHitmatched(
        firstTime, flash_pos, crtHitList, fTimeOfFlightInterval, isRealData, fGlobalT0Offset);
      auto const nCRTHits = crtMatches.entering.size() + crtMatches.exiting.size();
      
      std::vector<MatchedCRT> thisFlashCRTmatches;
      MatchType const eventType = crtMatches.flashType;
      if (nCRTHits > 0) {
          
          mf::LogTrace("CRTPMTMatchingProducer")
            << "Entering matches (" << crtMatches.entering.size() << "):";
          for (auto const& entering : crtMatches.entering) {
            thisFlashCRTmatches.push_back(makeMatchedCRT(*entering.CRTHit, tflash, isRealData));
          }
          
          mf::LogTrace("CRTPMTMatchingProducer")
            << "Exiting matches (" << crtMatches.exiting.size() << "):";
          for (auto const& exiting : crtMatches.exiting) {
            thisFlashCRTmatches.push_back(makeMatchedCRT(*exiting.CRTHit, tflash, isRealData));
          }
          
        }
        if (!thisFlashCRTmatches.empty() )
          mf::LogTrace("CRTPMTMatchingProducer") << "pushing back flash with " << thisFlashCRTmatches.size() << " CRT Matches.";
        FlashType thisFlashType = { /* .flashPos = */ flash_pos, // C++20: restore initializers
                                    /* .flashTime = */ tflash,
                                    /* .flashGateTime = */ thisRelGateTime / 1000.0, // -> us
                                    /* .inBeam = */ thisInTime_beam,
                                    /* .inGate = */ thisInTime_gate,
                                    /* .classification = */ eventType,
                                    /* .CRTmatches = */ std::move(thisFlashCRTmatches)};
        thisEventFlashes.push_back(std::move(thisFlashType));
      } // end of this flash
      mf::LogTrace("CRTPMTMatchingProducer") << "Event has " << thisEventFlashes.size() << " flashes in " << flashLabel.encode();
      for (auto const& theseFlashes : thisEventFlashes){
        CRTPMTMatching ProducedFlash = FillCRTPMT(theseFlashes);
        CRTPMTMatchesColl->push_back(std::move(ProducedFlash));
      }
    }
    //art::PtrMaker<sbn::crt::CRTHit> makeHitPtr(event);
    mf::LogTrace("CRTPMTMatchingProducer") <<"This Event has "<<CRTPMTMatchesColl->size()<<"  Flashes candidate for CRT matching."<<std::endl;
    e.put(std::move(CRTPMTMatchesColl));
    //e.put(std::move(FlashAssociation));

  } // CRTPMTMatchingProducer::produce()


  // ---------------------------------------------------------------------------
  icarus::crt::MatchedCRT CRTPMTMatchingProducer::makeMatchedCRT
    (sbn::crt::CRTHit const& hit, double tflash, bool isRealData) const
  {
    geo::Point_t const thisCRTpos {hit.x_pos, hit.y_pos, hit.z_pos};
    // fGlobalT0Offset = 1.6e6 ns = 1600000 ns.
    double const CRTtime  // ns -> us
      = icarus::crt::CRTHitTime(hit, fGlobalT0Offset, isRealData) / 1000.0;
    double const CRTTof_opflash = CRTtime - tflash; // us
    int const CRTRegion = hit.plane;
    // Very lazy way to determine if the Hit is a Top or a Side;
    // Will update it when bottom CRT will be available:
    int const CRTSys = (CRTRegion >= 36)? 1: 0;
    mf::LogTrace("CRTPMTMatchingProducer")
      << "  Match: tof = crtTime - tflash  = " << CRTtime
      << " - "<< tflash << " = " << CRTTof_opflash << " (us)";
    return { /* .CRTHitPos =*/  thisCRTpos, // C++20: restore initializers
             /* .CRTPMTTimeDiff = */ CRTTof_opflash,
             /* .CRTTime = */ CRTtime,
             /* .CRTSys = */ CRTSys,
             /* .CRTRegion = */ CRTRegion};
  }

  DEFINE_ART_MODULE(CRTPMTMatchingProducer)

} //end namespace
