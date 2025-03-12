////////////////////////////////////////////////////////////////////////
// Class:       TPCPMTBarycenterMatchProducer
// Plugin Type: producer (Unknown Unknown)
// File:        TPCPMTBarycenterMatchProducer_module.cc
//
// Generated at Sun Oct 22 14:43:16 2023 by John Smedley using cetskelgen
// from  version .
//
//  @file   icaruscode/PMT/OpReco/TPCPMTBarycenterMatchProducer_module.cc
//  @brief  Producer to match Pandora slices to their best match OpFlash by minimizing barycenter distance, as well as compare slices to the triggering OpFlash
//  @author Jack Smedley ( jsmedley@fnal.gov )
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "icarusalg/Utilities/TrackTimeInterval.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

//Data product includes
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/TriggerData.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

//ROOT includes
#include "TTree.h"
#include "TVector3.h"

#include <cmath> // std::hypot(), std::abs(), std::sqrt()
#include <iostream>
#include <memory>
#include <string>
#include <utility> // std::move()
#include <vector>

using microseconds = util::quantities::intervals::microseconds;
using electronics_time = detinfo::timescales::electronics_time;


/**
 * @brief Matches optical flashes and charge slices based on their location.
 * 
 * 
 * This algorithm associates slices of charge from the TPC (`recob::Slice`) to
 * reconstructed scintillation light flashes (`recob::OpFlash`).
 * 
 * For each one of the `InputTags` tags, one slice and one flash set are
 * identified and matching is done between them.
 * 
 * Each slice is associated independently to the closest among the flashes.
 * Distance is computed between the charge centroid of the slice and the one 
 * of the flash.
 * 
 * All flashes are candidate for matching with any given slice. If
 * `UseTimeRange` is specified, though, an allowed time range is determined
 * for the slice based on the drift time of its TPC hits (in the extreme case,
 * a slice covering from anode to cathode has an allowed time interval reduced
 * to a single instant). In this case, only flashes that are no farther from
 * this interval (with a set margin, `TimeRangeMargin`, accommodating for some
 * reconstruction biases) are kept in the candidate pool.
 * 
 * Matching of a slice immediately fails if no charge is associated to it, or,
 * in the case `UseTimeRange` is set, if there are no reconstructed flash
 * compatible with the time interval of the slice. Whenever the candidate pool
 * contains at least a flash, the matching is "successful", meaning that it
 * yields a result. It is the privilege of the analyzer to decide whether that
 * is a good match or not, based on the information in the matching data
 * product associated to the slice. If the matching did not succeed, matching
 * information will still be associated with the slice, but most or all of its
 * information will be placeholders.
 * 
 * It is possible that a flash ends being associated with multiple slices.
 * 
 * In addition, if a trigger time is available, the algorithm will attempt to
 * identify the triggering flash and, if found, will provide the distance of
 * each slice from that flash.
 * 
 * Content of the matching information:
 * 
 * * `ChargeT0` (microseconds): time associated to the slice via particle flow
 *     object; it is left in the original time reference. If no `anab::T0`
 *     object is associated, this variable is assigned the magic value `-9999`.
 * 
 * ### Definition of the centroids
 * 
 * The centroid of the slice is defined as the position of all the valid space
 * points associated with the slice, weighted by their charge.
 * Space points (`recob::SpacePoint`) are associated each to a TPC hit
 * (`recob::Hit`); the charge of the hit is used as a weight, with no further
 * calibration applied.
 * If `CollectionOnly` is specified, only space points associated to hits on
 * a collection plane are included.
 * 
 * The centroid of a flash is defined only in its projection on the PMT plane.
 * Its definition is delegated to the flash reconstruction algorithm
 * (`recob::OpFlash::YCenter()` and `recob::OpFlash::ZCenter()`).
 * In particular note that there is no attempt to enforce a location based on
 * the earliest light. In fact, the standard flash reconstruction algorithm
 * builds the center of the flash from _all_ the associated hits no matter
 * their time.
 * 
 * ### Determination of the trigger flash
 * 
 * If a trigger time is found from the `TriggerLabel` data product, the module
 * will attempt to identify the triggering flash.
 * The identification is based solely on the time of the trigger and of the
 * flash.
 * In ICARUS data, the time reference is the trigger time, therefore the time
 * of the global trigger is a fixed value by definition. Nevertheless, to
 * accommodate relative delays between the trigger and the light systems, an
 * expected time difference between the reconstrcuted flash (`TriggerDelay`)
 * is accounted for with a cushion ('TriggerTolerance'). These parameters are
 * measured in microseconds and controlled directly by the module caller.
 * In ICARUS simulation the mechanism is the same; the times are stored with
 * respect to the simulation time, but the module will take care of
 * the appropriate conversion, so the configuration parameter holds the same
 * meaning and scale as in data.
 * 
 * 
 * Input
 * ------
 * 
 * * `std::vector<recob::OpFlash>` (based on `OpFlashLabel`): collection of
 *     reconstructed flashes to be associated to the charge; also their
 *     associations to `recob::OpHit` (and those hits themselves).
 * * `std::vector<recob::Slice>` (based on `PandoraLabel`): collection of
 *     reconstructed TPC charge slices to be associated to the light; also their
 *     associations to `recob::Hit` and `recob::PFParticle`, and
 *     `recob::SpacePoint` objects associated to those hits
 *     (and also the associated objects as well).
 * * `std::vector<raw::Trigger>` (`TriggerLabel`): collection of the global
 *     triggers (only at most one is expected).
 * 
 * 
 * Output
 * -------
 * 
 * A single collection, merging all the input, is produced for each of the
 * following data products in the _art_/ROOT output:
 * 
 * * `std::vector<sbn::TPCPMTBarycenterMatch>`: collection of matching information;
 *     one matching information object is present for each slice in the input
 *     collections, in the same order as the input.
 * * `art::Assns<sbn::TPCPMTBarycenterMatch, recob::Slice>`: association of each
 *     matched slice with its matching information.
 * * `art::Assns<sbn::TPCPMTBarycenterMatch, recob::OpFlash>`: association of each
 *     matched light flash with its matching information.
 * 
 * Note that while there is currently no direct association between the slice
 * and the flash, the information contained in `sbn::TPCPMTBarycenterMatch` is very
 * detailed.
 * 
 * In addition, if `FillMatchTree` is set, a ROOT tree called `matchTree` will
 * be written via `TFileService`. The tree contains one entry per slice,
 * whether matched or not,
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * * `InputTags` (list of strings, mandatory): suffixes to be added to the
 *     labels of the reconstructed data products to be used in the matching.
 *     Two typical choices are to process independently parts of the detector,
 *     e.g. `[ "cryoE", "cryoW" ]` to process east cryostat first and west one
 *     next, or to have a single empty suffix (`[ "" ]`) to use as input tags
 *     exactly the labels as specified by the following parameters.
 * * `OpFlashLabel` (string, mandatory): base of the tag of the input flashes.
 * * `PandoraLabel` (string, mandatory): base of the tag of input slices, and
 *     their associations to hits, space points, particle flow objects etc.
 * * `TriggerLabel` (input tag, mandatory): information on the global trigger
 *     (hardware or synthetic).
 * * `CollectionOnly` (flag, default: `true`): if set, only hits from
 *     collection planes will contribute to the centroid of the TPC slice.
 * * `UseTimeRange` (flag, default: `true`): if set, a slice is assigned an
 *     allowed time interval based on the drift time of its hits and only
 *     flashes falling in that interval are considered for matching. Otherwise,
 *     all flashes are game for matching.
 * * `Verbose` (flag, default: `false`): enables verbose output directly to
 *     console standard output.
 * * `FillMatchTree` (flag, default: `false`): if set to `true`, a ROOT tree
 *     with detailed matching information called `"matchTree"` will be written
 *     via `TFileService`.
 * * `TriggerDelay` (real, mandatory, in microseconds): the expected time
 *     difference between the scintillation flash and the global trigger.
 * * `TriggerTolerance` (real, mandatory, in microseconds): the tolerance used
 *     to identify a light flash associated with the trigger.
 * * `TimeRangeMargin` (real, microseconds; default: `0`): when `UseTimeRange`
 *     is set, the allowed time interval for each slice is extended on both
 *     sides by this amount of time.
 * 
 */
class TPCPMTBarycenterMatchProducer : public art::EDProducer {
public:
  explicit TPCPMTBarycenterMatchProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCPMTBarycenterMatchProducer(TPCPMTBarycenterMatchProducer const&) = delete;
  TPCPMTBarycenterMatchProducer(TPCPMTBarycenterMatchProducer&&) = delete;
  TPCPMTBarycenterMatchProducer& operator=(TPCPMTBarycenterMatchProducer const&) = delete;
  TPCPMTBarycenterMatchProducer& operator=(TPCPMTBarycenterMatchProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  void InitializeSlice();                                                                     ///< Re-initialize all slice-level data members
  double CentroidOverlap(double center1, double center2, double width1, double width2) const; ///< Return overlap between charge and light centroids OR distance apart if no overlap
  double CalculateAsymmetry(art::Ptr<recob::OpFlash> flash, int cryo);                        ///< Return the east-west asymmetry of PEs in a given OpFlash
  void updateChargeVars(double sumCharge, TVector3 const& sumPos, TVector3 const& sumPosSqr, std::array<double, 2> const& triggerFlashCenter); ///< Update slice-level data members with charge and trigger match info
  void updateFlashVars(art::Ptr<recob::OpFlash> flash, double firstHit);                      ///< Update slice-level data members with best match info
  void updateMatchInfo(sbn::TPCPMTBarycenterMatch& matchInfo);                                      ///< Update match product with slice-level data members
 
  // Input parameters
  std::vector<std::string>  fInputTags;            ///< Suffix added onto fOpFlashLabel and fPandoraLabel, used by ICARUS for separate cryostat labels but could be empty
  std::string               fOpFlashLabel;         ///< Label for PMT reconstruction products
  std::string               fPandoraLabel;         ///< Label for Pandora output products
  std::string               fTriggerLabel;         ///< Label for trigger product
  bool                      fCollectionOnly;       ///< Only use TPC spacepoints from the collection plane
  bool                      fUseTimeRange;         ///< Reject impossible matches based on allowed time range of TPC hits relative to trigger 
  bool                      fVerbose;              ///< Print extra info
  bool                      fFillMatchTree;        ///< Fill an output TTree in the supplemental file
  double                    fTriggerDelay;         ///< Typical time of triggering flash, EYEBALLED (us)
  double                    fTriggerTolerance;     ///< Spread of triggering flash times, EYEBALLED (us)
  microseconds              fTimeRangeMargin;      ///< Symmetric acceptable margin for allowed time range of TPC hits (us)
  
  // Event-level data members
  int                       fRun;                  ///< Number of the run being processed
  int                       fEvent;                ///< Number of the event being processed
  int                       fCryo;                 ///< Cryostat this event occured in
  int                       fSliceNum;             ///< Number of slice in the event
  // Slice-level data members 
  double                    fChargeT0;             ///< Start time for cathode-crossing PFPs, not always available (us)
  double                    fChargeTotal;          ///< Total charge in slice (integrated ADC counts)
  double                    fChargeCenterXGlobal;  ///< Weighted mean X position of spacepoints (cm)
  double                    fChargeCenterXLocal;   ///< Weighted mean X position of spacepoints, measured with respect to the cathode (cm)
  double                    fChargeCenterY;        ///< Weighted mean Y position of spacepoints (cm)
  double                    fChargeCenterZ;        ///< Weighted mean Z position of spacepoints (cm)
  double                    fChargeWidthX;         ///< Weighted standard deviation of X position of spacepoints (cm)
  double                    fChargeWidthY;         ///< Weighted standard deviation of Y position of spacepoints (cm)
  double                    fChargeWidthZ;         ///< Weighted standard deviation of Z position of spacepoints (cm)
  double                    fFlashFirstHit;        ///< Earliest OpHit time in matched OpFlash (us)
  double                    fFlashTime;            ///< Matched OpFlash time (us)
  double                    fFlashPEs;             ///< Brightness of matched flash (photoelectrons)
  double                    fFlashAsymmetry;       ///< East-West asymmetry of PEs in matched flash
  double                    fFlashCenterY;         ///< Weighted mean Y postion of hit PMTs (cm)
  double                    fFlashCenterZ;         ///< Weighted mean Z postion of hit PMTs (cm)
  double                    fFlashWidthY;          ///< Weighted standard deviation of Y postion of hit PMTs (cm)
  double                    fFlashWidthZ;          ///< Weighted standard deviation of Z postion of hit PMTs (cm)
  double                    fDeltaT;               ///< | Matched flash time - charge T0 | when available (us)
  double                    fDeltaY;               ///< | Matched flash Y center - charge Y center | (cm)
  double                    fDeltaZ;               ///< | Matched flash Z center - charge Z center | (cm)
  double                    fRadius;               ///< Hypotenuse of DeltaY and DeltaZ *parameter minimized by matching* (cm)
  double                    fOverlapY;             ///< Spatial overlap of flash and charge centroids in Y [>0] OR distance apart if no overlap [<0] (cm)
  double                    fOverlapZ;             ///< Spatial overlap of flash and charge centroids in Z [>0] OR distance apart if no overlap [<0] (cm)
  double                    fDeltaY_Trigger;       ///< | Triggering flash Y center - charge Y center | (cm)
  double                    fDeltaZ_Trigger;       ///< | Triggering flash Z center - charge Z center | (cm)
  double                    fRadius_Trigger;       ///< Hypotenuse of DeltaY_Trigger and DeltaZ_Trigger (cm)
  
  TTree*                    fMatchTree;            ///< Tree to store all match information
  
  // Detector geometry and properties
  geo::GeometryCore const&                  fGeom;
  detinfo::DetectorClocksService const&     fDetClocks;
  detinfo::DetectorPropertiesService const& fDetProp; 
  lar::util::TrackTimeIntervalMaker const   fTimeIntervalMaker;
};


TPCPMTBarycenterMatchProducer::TPCPMTBarycenterMatchProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  // More initializers here.
  fInputTags(p.get<std::vector<std::string>>("InputTags")),
  fOpFlashLabel(p.get<std::string>("OpFlashLabel")),
  fPandoraLabel(p.get<std::string>("PandoraLabel")),
  fTriggerLabel(p.get<std::string>("TriggerLabel")),
  fCollectionOnly(p.get<bool>("CollectionOnly", true)),
  fUseTimeRange(p.get<bool>("UseTimeRange", true)),
  fVerbose(p.get<bool>("Verbose", false)),
  fFillMatchTree(p.get<bool>("FillMatchTree", false)),
  fTriggerDelay(p.get<double>("TriggerDelay", 0.0)),
  fTriggerTolerance(p.get<double>("TriggerTolerance", 25.0)),
  fTimeRangeMargin(fUseTimeRange? microseconds{p.get<double>("TimeRangeMargin")}: microseconds{0.0}),
  fGeom(*lar::providerFrom<geo::Geometry>()),
  fDetClocks(*art::ServiceHandle<detinfo::DetectorClocksService>()),
  fDetProp(*art::ServiceHandle<detinfo::DetectorPropertiesService>()),
  fTimeIntervalMaker{ fGeom }
{
  // Call appropriate produces<>() functions here.

  produces< std::vector<sbn::TPCPMTBarycenterMatch> >();
  produces< art::Assns<sbn::TPCPMTBarycenterMatch, recob::Slice> >();
  produces< art::Assns<sbn::TPCPMTBarycenterMatch, recob::OpFlash> >();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<raw::Trigger>>(fTriggerLabel);
  for ( const std::string& inputTag : fInputTags ) {
    consumes<std::vector<recob::OpFlash>>(fOpFlashLabel + inputTag);
    consumes<std::vector<recob::Slice>>(fPandoraLabel + inputTag);
    
    // via art::FindMany:
    consumes<art::Assns<recob::OpFlash, recob::OpHit>>(fOpFlashLabel + inputTag);
    consumes<art::Assns<recob::Slice, recob::Hit>>(fPandoraLabel + inputTag);
    consumes<art::Assns<recob::Slice, recob::PFParticle>>(fPandoraLabel + inputTag);
    consumes<art::Assns<recob::Hit, recob::SpacePoint>>(fPandoraLabel + inputTag);
    consumes<art::Assns<recob::PFParticle, anab::T0>>(fPandoraLabel + inputTag);
  } // for

  if ( fFillMatchTree ) {
    art::ServiceHandle<art::TFileService> tfs;
    fMatchTree = tfs->make<TTree>("matchTree","TPC Slice - OpFlash Matching Analysis");

    //Event Info
    fMatchTree->Branch("run",                 &fRun,                 "run/I"                );
    fMatchTree->Branch("event",               &fEvent,               "event/I"              );
    fMatchTree->Branch("cryo",                &fCryo,                "cryo/I"               );
    fMatchTree->Branch("sliceNum",            &fSliceNum,            "sliceNum/I"           );

    //Charge Info
    fMatchTree->Branch("chargeT0",            &fChargeT0,            "chargeT0/d"           );
    fMatchTree->Branch("chargeTotal",         &fChargeTotal,         "chargeTotal/d"        );
    fMatchTree->Branch("chargeCenterXGlobal", &fChargeCenterXGlobal, "chargeCenterXGlobal/d");
    fMatchTree->Branch("chargeCenterXLocal",  &fChargeCenterXLocal,  "chargeCenterXLocal/d" );
    fMatchTree->Branch("chargeCenterY",       &fChargeCenterY,       "chargeCenterY/d"      );
    fMatchTree->Branch("chargeCenterZ",       &fChargeCenterZ,       "chargeCenterZ/d"      );
    fMatchTree->Branch("chargeWidthX",        &fChargeWidthX,        "chargeWidthX/d"       );
    fMatchTree->Branch("chargeWidthY",        &fChargeWidthY,        "chargeWidthY/d"       );
    fMatchTree->Branch("chargeWidthZ",        &fChargeWidthZ,        "chargeWidthZ/d"       );

    //Matched Flash Info
    fMatchTree->Branch("flashFirstHit",       &fFlashFirstHit,       "flashFirstHit/d"      );
    fMatchTree->Branch("flashTime",           &fFlashTime,           "flashTime/d"          );
    fMatchTree->Branch("flashPEs",            &fFlashPEs,            "flashPEs/d"           );
    fMatchTree->Branch("flashAsymmetry",      &fFlashAsymmetry,      "flashAsymmetry/d"     );
    fMatchTree->Branch("flashCenterY",        &fFlashCenterY,        "flashCenterY/d"       );
    fMatchTree->Branch("flashCenterZ",        &fFlashCenterZ,        "flashCenterZ/d"       );
    fMatchTree->Branch("flashWidthY",         &fFlashWidthY,         "flashWidthY/d"        );
    fMatchTree->Branch("flashWidthZ",         &fFlashWidthZ,         "flashWidthZ/d"        );

    //Match Quality Info
    fMatchTree->Branch("deltaT",              &fDeltaT,              "deltaT/d"             );
    fMatchTree->Branch("deltaY",              &fDeltaY,              "deltaY/d"             );
    fMatchTree->Branch("deltaZ",              &fDeltaZ,              "deltaZ/d"             );
    fMatchTree->Branch("radius",              &fRadius,              "radius/d"             );
    fMatchTree->Branch("overlapY",            &fOverlapY,            "overlapY/d"           );
    fMatchTree->Branch("overlapZ",            &fOverlapZ,            "overlapZ/d"           );
    fMatchTree->Branch("deltaZ_Trigger",      &fDeltaZ_Trigger,      "deltaZ_Trigger/d"     );
    fMatchTree->Branch("deltaY_Trigger",      &fDeltaY_Trigger,      "deltaY_Trigger/d"     );
    fMatchTree->Branch("radius_Trigger",      &fRadius_Trigger,      "dadius_Trigger/d"     );

  } //End MatchTree

}

void TPCPMTBarycenterMatchProducer::produce(art::Event& e)
{
  // Implementation of required member function here.
  fEvent  = e.id().event();
  fRun    = e.run();

  //Fetch trigger info and if MC check whether this event triggered
  art::Handle const triggerHandle
    = e.getHandle<std::vector<raw::Trigger>>(fTriggerLabel);
  double triggerTime = -9999.;
  if ( triggerHandle.isValid() && triggerHandle->size() >= 1 ) {
    raw::Trigger const& trigger = triggerHandle->at(0);
    triggerTime = trigger.TriggerTime();
    if ( fVerbose ) std::cout << "Valid trigger product found. Trigger time: " << trigger.TriggerTime() << ", Beam gate time: " << trigger.BeamGateTime() << ", Difference: " << trigger.TriggerTime() - trigger.BeamGateTime() << std::endl;
  }
  else if ( fVerbose ) std::cout << "No valid trigger product found for this event!"  << std::endl;

  //Infrastructure for checking allowed time range of a slice
  detinfo::DetectorTimings const detTimings{ fDetClocks.DataFor(e) };
  detinfo::DetectorPropertiesData const& detProp { fDetProp.DataFor(e, detTimings.clockData()) };
  lar::util::TrackTimeInterval const timeIntervals = fTimeIntervalMaker(detProp, detTimings);

  //Initialize new data products
  auto matchInfoVector = std::make_unique< std::vector<sbn::TPCPMTBarycenterMatch> >();
  art::PtrMaker< sbn::TPCPMTBarycenterMatch > const makeInfoPtr(e); 
  auto sliceAssns = std::make_unique< art::Assns<sbn::TPCPMTBarycenterMatch, recob::Slice> >();
  auto flashAssns = std::make_unique< art::Assns<sbn::TPCPMTBarycenterMatch, recob::OpFlash> >();

  //For InputTag...
  for ( const std::string& inputTag : fInputTags ) {
    //East-->0, West-->1
    fCryo = ( inputTag.find("W") != std::string::npos ) ? 1 : 0;


/* ~~~~~~~~~~~~~~~~~~~~ Flash Section
 *
 * Here we gather the OpFlashes found in this cryostat and their OpHits
 * We iterate through the flashes to identify a triggering flash
 */

    //Fetch the flashes and their associated hits, pointer vector needed for generating associations
    art::Handle const flashHandle
      = e.getHandle<std::vector<recob::OpFlash>>(fOpFlashLabel + inputTag);
    art::FindMany<recob::OpHit> fmOpHits(flashHandle, e, fOpFlashLabel + inputTag);
    int nFlashes = (*flashHandle).size();

    std::array<double, 2> triggerFlashCenter = {-9999., -9999.};
    //If triggered event..
    if ( triggerTime >= 0. ) {
      double minTimeDiff = 1.6;
      //For flash...
      for ( int i = 0; i < nFlashes; i++ ) {
        const recob::OpFlash &flash = (*flashHandle)[i];

        //Identify the triggering flash
        double timeDiff = abs( (triggerTime - flash.AbsTime()) - fTriggerDelay );
        if ( timeDiff < fTriggerTolerance && timeDiff < minTimeDiff ) {
          triggerFlashCenter = {flash.YCenter(), flash.ZCenter()};
          minTimeDiff = timeDiff;
        }

      } //End for flash
    } //End if triggered event

    if ( fVerbose ) std::cout << "Event: " << fEvent << ", Cryo: " << inputTag << ", nFlashes: " << nFlashes << ", Triggering flash center Y: " << triggerFlashCenter[0] << ", Triggering flash center Z: " << triggerFlashCenter[1]  << std::endl;


/* ~~~~~~~~~~~~~~~~~~~~ TPC Section
 * Here we start by gathering the Slices in the event
 * For each slice, the charge centroid is first calculated
 * Then we iterate through flashes to identify the best match flash
 * If a triggering flash was found in this cyrostat, the barycenter distance to the triggering flash is also stored
 */

    //Fetch slices, TPC hits, and PFPs; pointer vector needed for generating associations
    art::Handle const sliceHandle
      = e.getHandle<std::vector<recob::Slice>>(fPandoraLabel + inputTag);
    art::FindManyP<recob::Hit> fmTPCHits(sliceHandle, e, fPandoraLabel + inputTag);
    art::FindManyP<recob::PFParticle> fmPFPs(sliceHandle, e, fPandoraLabel + inputTag);

    unsigned nSlices = (*sliceHandle).size();

    //For slice...
    for ( unsigned j = 0; j < nSlices; j++ ) {
      fSliceNum = j;
      const art::Ptr<recob::Slice> slicePtr { sliceHandle, j };
      InitializeSlice();
      sbn::TPCPMTBarycenterMatch sliceMatchInfo;
      updateMatchInfo(sliceMatchInfo);

      const std::vector<art::Ptr<recob::Hit>> &tpcHitsVec = fmTPCHits.at(j);
      const std::vector<art::Ptr<recob::PFParticle>> &pfpsVec = fmPFPs.at(j);
      art::FindOne<recob::SpacePoint> f1SpacePoint(tpcHitsVec, e, fPandoraLabel + inputTag);

      int nHits = tpcHitsVec.size();
      int nPFPs = pfpsVec.size();

      //Establish possible time range for this slice
      lar::util::TrackTimeInterval::TimeRange const& timeRange = timeIntervals.timeRangeOfHits(tpcHitsVec);
      const bool rangeIsValid = timeRange.isValid();

      //Retrieve Pandora's T0 for this slice if available, same for every PFP in slice so we only need one
      if ( nPFPs != 0 ) {
        art::FindOne<anab::T0> f1T0( {pfpsVec.at(0)}, e, fPandoraLabel + inputTag);
        if ( f1T0.at(0).isValid() ) {
          fChargeT0 = f1T0.at(0).ref().Time() / 1e3;
        }
      }

      double thisCharge;
      double sumCharge = 0.;
      TVector3 sumPos {0.,0.,0.};
      TVector3 sumPosSqr {0.,0.,0.};

      //For hit...
      for ( int k = 0; k < nHits; k++ ) {
        const art::Ptr<recob::Hit> &tpcHit = tpcHitsVec.at(k);

        //Only use hits with associated SpacePoints, and optionally only collection plane hits
        if ( fCollectionOnly && tpcHit->SignalType() != geo::kCollection ) continue;
        if ( !f1SpacePoint.at(k).isValid() ) continue;

        const recob::SpacePoint point = f1SpacePoint.at(k).ref();
        thisCharge = tpcHit->Integral();
        TVector3 const thisPoint = point.XYZ();
        TVector3 const thisPointSqr {thisPoint.X()*thisPoint.X(), thisPoint.Y()*thisPoint.Y(), thisPoint.Z()*thisPoint.Z()};
        sumCharge += thisCharge;
        sumPos += thisPoint * thisCharge;
        sumPosSqr += thisPointSqr * thisCharge;
      } //End for hit

      //No charge found in slice...
      if ( sumCharge == 0. ) {
        if ( fFillMatchTree ) fMatchTree->Fill();
        art::Ptr<sbn::TPCPMTBarycenterMatch> const infoPtr = makeInfoPtr(matchInfoVector->size());
        sliceAssns->addSingle(infoPtr, slicePtr);
        matchInfoVector->push_back(std::move(sliceMatchInfo));
        if ( fVerbose ) std::cout << "No charge found in Event: " << fEvent << " Slice: " << j << "! Continuing..."  << std::endl;
        continue;
      }

      //Update charge variables
      updateChargeVars(sumCharge, sumPos, sumPosSqr, triggerFlashCenter);
      updateMatchInfo(sliceMatchInfo);

      int matchIndex = -5;
      double minDistance = 1e6;
      double thisFlashCenterY, thisFlashCenterZ, thisDistance;

      //For flash...
      for ( int m = 0; m < nFlashes; m++ ) {
        const recob::OpFlash &flash = (*flashHandle).at(m);

        //Skip over flashes that are very out of time with respect to the slice
        if ( fUseTimeRange && rangeIsValid ) {
          electronics_time eTime (flash.AbsTime());
          if ( !timeRange.contains(eTime, fTimeRangeMargin) ) continue;
        }

        //TODO: if ( flash has entering CRT match ) continue? Or at least just store that as a bool?

        //Find index of flash that minimizes barycenter distance in YZ place
        thisFlashCenterY = flash.YCenter();
        thisFlashCenterZ = flash.ZCenter();
        thisDistance = std::hypot( (thisFlashCenterY - fChargeCenterY), (thisFlashCenterZ - fChargeCenterZ) );
        if ( thisDistance < minDistance ) {
          minDistance = thisDistance;
          matchIndex = m;
        }
      } //End for flash

      //No valid match found...
      if ( matchIndex == -5 ) {
        if ( fFillMatchTree ) fMatchTree->Fill();
        art::Ptr<sbn::TPCPMTBarycenterMatch> const infoPtr = makeInfoPtr(matchInfoVector->size());
        sliceAssns->addSingle(infoPtr, slicePtr);
        matchInfoVector->push_back(std::move(sliceMatchInfo));
        if ( fVerbose ) std::cout << "No matching flash found for Event: " << fEvent << " Slice: " << j << "! Continuing..."  << std::endl;
        continue;
      }

      //Best match flash pointer
      unsigned unsignedMatchIndex = matchIndex;
      const art::Ptr<recob::OpFlash> flashPtr { flashHandle, unsignedMatchIndex };

      //Find time of first OpHit in matched flash
      const std::vector<recob::OpHit const*> &opHitsVec = fmOpHits.at(matchIndex);
      double minTime = 1e6;
      for (const recob::OpHit *opHit : opHitsVec ) { if ( opHit->PeakTime() < minTime ) minTime = opHit->PeakTime(); }

      //Update match info
      updateFlashVars(flashPtr, minTime);
      updateMatchInfo(sliceMatchInfo);
      art::Ptr<sbn::TPCPMTBarycenterMatch> const infoPtr = makeInfoPtr(matchInfoVector->size());
      sliceAssns->addSingle(infoPtr, slicePtr);
      flashAssns->addSingle(infoPtr, flashPtr);
      matchInfoVector->push_back(std::move(sliceMatchInfo));
      if ( fFillMatchTree ) fMatchTree->Fill();

    } //End for slice

  } //End for InputTag

  //Store new products at the end of the event
  e.put(std::move(matchInfoVector));
  e.put(std::move(sliceAssns));
  e.put(std::move(flashAssns));

} //End produce()

void TPCPMTBarycenterMatchProducer::InitializeSlice() {
  fChargeT0 = -9999.;
  fChargeTotal = -9999.;
  fChargeCenterXGlobal = -9999.;
  fChargeCenterXLocal = -9999.;
  fChargeCenterY = -9999.;
  fChargeCenterZ = -9999.;
  fChargeWidthX = -9999.;
  fChargeWidthY = -9999.;
  fChargeWidthZ = -9999.;
  fFlashFirstHit = -9999.;
  fFlashTime = -9999.;
  fFlashPEs = -9999.;
  fFlashCenterY = -9999.;
  fFlashCenterZ = -9999.;
  fFlashWidthY = -9999.;
  fFlashWidthZ = -9999.;
  fDeltaT = -9999.;
  fDeltaY = -9999.;
  fDeltaZ = -9999.;
  fRadius = -9999.;
  fOverlapY = -9999.;
  fOverlapZ = -9999.;
  fDeltaZ_Trigger = -9999.;
  fDeltaY_Trigger = -9999.;
  fRadius_Trigger = -9999.;
} //End InitializeSlice()


double TPCPMTBarycenterMatchProducer::CentroidOverlap(double center1, double center2, double width1, double width2) const {
  //Centroid 2 is contained within Centroid 1, so overlap is the whole Centroid 2
  if ( (center1 - width1 < center2 - width2) && (center1 + width1 > center2 + width2) ) return (2 * width2);

  //Centroid 1 is contained within Centroid 2, so overlap is the whole Centroid 1
  else if ( (center1 - width1 > center2 - width2) && (center1 + width1 < center2 + width2) ) return (2 * width1);

  double difference = center1 - center2;
  if ( center1 > center2 ) difference *= -1;
  return difference + width1 + width2;
} //End CentroidOverlap()


double TPCPMTBarycenterMatchProducer::CalculateAsymmetry(art::Ptr<recob::OpFlash> flash, int cryo) {
  double sumEast = 0.;
  double sumWest = 0.;

  //East Cryo flashes have a 180-element PE vector; 0-89 --> East wall, 90-179 --> West wall
  //West Cryo flashes have a 360-element PE vector; 0-179 --> ALL 0, 180-269 --> East wall, and 270-359 --> West wall
  int countingOffset = 0;
  if ( cryo == 1 ) countingOffset += 180;

  for ( int PMT = 0; PMT < 180; PMT++ ) {
    if ( PMT <= 89 ) sumEast += flash->PEs().at(PMT + countingOffset);
    else sumWest += flash->PEs().at(PMT + countingOffset);
  }

  return (sumWest - sumEast) / (sumWest + sumEast);
} //End CalculateAsymmetry()


//TODO: Get the cathode position and shift global X to local X in a less hacky way
//According to a geometrydump, the cathode X positions are +/-(210.14, 210.29), depending on the TPC. Here I just averaged those...
void TPCPMTBarycenterMatchProducer::updateChargeVars(double sumCharge, TVector3 const& sumPos, TVector3 const& sumPosSqr, std::array<double, 2> const& triggerFlashCenter) {
  fChargeCenterXGlobal = sumPos[0] / sumCharge;
  fChargeCenterXLocal = fChargeCenterXGlobal - 210.215 * (2*fCryo - 1);
  fChargeCenterY = sumPos[1] / sumCharge;
  fChargeCenterZ = sumPos[2] / sumCharge;
  fChargeWidthX = std::sqrt( sumPosSqr[0]/sumCharge - (sumPos[0]/sumCharge)*(sumPos[0]/sumCharge) );
  fChargeWidthY = std::sqrt( sumPosSqr[1]/sumCharge - (sumPos[1]/sumCharge)*(sumPos[1]/sumCharge) );
  fChargeWidthZ = std::sqrt( sumPosSqr[2]/sumCharge - (sumPos[2]/sumCharge)*(sumPos[2]/sumCharge) );
  if ( triggerFlashCenter[1] != -9999. ) {
    fDeltaY_Trigger = abs(triggerFlashCenter[0] - fChargeCenterY);
    fDeltaZ_Trigger = abs(triggerFlashCenter[1] - fChargeCenterZ);
    fRadius_Trigger = std::hypot(fDeltaY_Trigger, fDeltaZ_Trigger);
  }
} //End updateChargeVars()


void TPCPMTBarycenterMatchProducer::updateFlashVars(art::Ptr<recob::OpFlash> flash, double firstHit) {
  double matchedTime = flash->Time();
  double matchedYCenter = flash->YCenter();
  double matchedZCenter = flash->ZCenter();
  double matchedYWidth = flash->YWidth();
  double matchedZWidth = flash->ZWidth();

  fFlashFirstHit = firstHit;
  fFlashTime = matchedTime;
  fFlashPEs =  flash->TotalPE();
  fFlashAsymmetry = CalculateAsymmetry(flash, fCryo);
  fFlashCenterY = matchedYCenter;
  fFlashCenterZ = matchedZCenter;
  fFlashWidthY = matchedYWidth;
  fFlashWidthZ = matchedZWidth;
  if ( fChargeT0 != -9999 ) fDeltaT = abs(matchedTime - fChargeT0);
  fDeltaY = abs(matchedYCenter - fChargeCenterY);
  fDeltaZ = abs(matchedZCenter - fChargeCenterZ);
  fRadius = std::hypot(fDeltaY, fDeltaZ);
  fOverlapY = CentroidOverlap(matchedYCenter, fChargeCenterY, matchedYWidth, fChargeWidthY);
  fOverlapZ = CentroidOverlap(matchedZCenter, fChargeCenterZ, matchedZWidth, fChargeWidthZ);
} //End updateFlashVars()


void TPCPMTBarycenterMatchProducer::updateMatchInfo(sbn::TPCPMTBarycenterMatch& matchInfo) {
  matchInfo.chargeTotal = fChargeTotal;
  matchInfo.chargeCenterXLocal = fChargeCenterXLocal;
  matchInfo.chargeCenter = {fChargeCenterXGlobal, fChargeCenterY, fChargeCenterZ};
  matchInfo.chargeWidth = {fChargeWidthX, fChargeWidthY, fChargeWidthZ};
  matchInfo.flashFirstHit = fFlashFirstHit;
  matchInfo.flashTime = fFlashTime;
  matchInfo.flashPEs = fFlashPEs;
  matchInfo.flashCenter = {-9999., fFlashCenterY, fFlashCenterZ};
  matchInfo.flashWidth = {-9999., fFlashWidthY, fFlashWidthZ};
  matchInfo.deltaT = fDeltaT;
  matchInfo.deltaY = fDeltaY;
  matchInfo.deltaZ = fDeltaZ;
  matchInfo.radius = fRadius;
  matchInfo.overlapY = fOverlapY;
  matchInfo.overlapZ = fOverlapZ;
  matchInfo.deltaY_Trigger = fDeltaY_Trigger;
  matchInfo.deltaZ_Trigger = fDeltaZ_Trigger;
  matchInfo.radius_Trigger = fRadius_Trigger;
} //End updateMatchInfo()


DEFINE_ART_MODULE(TPCPMTBarycenterMatchProducer)
