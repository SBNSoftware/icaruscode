////////////////////////////////////////////////////////////////////////
// Class:       BarycenterMatchProducer
// Plugin Type: producer (Unknown Unknown)
// File:        BarycenterMatchProducer_module.cc
//
// Generated at Sun Oct 22 14:43:16 2023 by John Smedley using cetskelgen
// from  version .
//
//  @file   icaruscode/PMT/OpReco/BarycenterMatchProducer_module.cc
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
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "icarusalg/Utilities/TrackTimeInterval.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
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
#include "sbnobj/Common/Reco/BarycenterMatch.h"

//ROOT includes
#include "TTree.h"
#include "TVector3.h"

#include <memory>

using microseconds = util::quantities::intervals::microseconds;
using electronics_time = detinfo::timescales::electronics_time;

class BarycenterMatchProducer;


class BarycenterMatchProducer : public art::EDProducer {
public:
  explicit BarycenterMatchProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BarycenterMatchProducer(BarycenterMatchProducer const&) = delete;
  BarycenterMatchProducer(BarycenterMatchProducer&&) = delete;
  BarycenterMatchProducer& operator=(BarycenterMatchProducer const&) = delete;
  BarycenterMatchProducer& operator=(BarycenterMatchProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  void InitializeSlice();                                                                                  ///< Re-initialize all slice-level data members
  double CentroidOverlap(double center1, double center2, double width1, double width2);                    ///< Return overlap between charge and light centroids OR distance apart if no overlap
  double CalculateAsymmetry(art::Ptr<recob::OpFlash> flash, int cryo);                                     ///< Return the east-west asymmetry of PEs in a given OpFlash
  void updateChargeVars(double sumCharge, TVector3 sumPos, TVector3 sumPosSqr, double triggerFlashCenter); ///< Update slice-level data members with charge and trigger match info
  void updateFlashVars(art::Ptr<recob::OpFlash> flash, double firstHit);                     ///< Update slice-level data members with best match info
  void updateMatchInfo(sbn::BarycenterMatch& matchInfo);                                                   ///< Update match product with slice-level data members
 
  // Input parameters
  std::vector<std::string>  fInputTags;            ///< Suffix added onto fOpFlashLabel and fPandoraLabel, used by ICARUS for separate cryostat labels but could be empty
  std::string               fOpFlashLabel;         ///< Label for PMT reconstruction products
  std::string               fPandoraLabel;         ///< Label for Pandora output products
  std::string               fTriggerLabel;         ///< Label for trigger product
  bool                      fCollectionOnly;       ///< Only use TPC spacepoints from the collection plane
  bool                      fUseTimeRange;         ///< Reject impossible matches based on allowed time range of TPC hits relative to trigger 
  bool                      fVerbose;              ///< Print extra info
  bool                      fFillMatchTree;        ///< Fill an output TTree in the supplemental file
  double                    fNominalTrigTime;      ///< Typical time of triggering flash, EYEBALLED (us)
  double                    fTriggerTolerance;     ///< Spread of triggering flash times, EYEBALLED (us)
  double                    fTimeRangeMargin;      ///< Symmetric acceptable margin for allowed time range of TPC hits (us)
  
  // Event-level data members
  int                       fRun;                  ///< Number of the run being processed
  int                       fEvent;                ///< Number of the event being processed
  int                       fCryo;                 ///< Cryostat this event occured in
  int                       fSliceNum;             ///< Number of slice in the event
  // Slice-level data members 
  double                    fChargeT0;             ///< Start time for cathode-crossing PFPs, not always available (us)
  double                    fChargeTotal;          ///< Total charge in slice
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
  double                    fRadius;               ///< Hypotenuse of DeltaY and DeltaZ, PARAMETER MINIMIZED BY MATCHING (cm)
  double                    fOverlapY;             ///< Spacial overlap of flash and charge centroids in Y [>0] OR distance apart if no overlap [<0] (cm)
  double                    fOverlapZ;             ///< Spacial overlap of flash and charge centroids in Z [>0] OR distance apart if no overlap [<0] (cm)
  double                    fDeltaZ_Trigger;       ///< | Triggering flash Z center - charge Z center | (cm)
  
  TTree*                    fMatchTree;            ///< Tree to store all match information
  
  // Detector geometry and properties
  geo::GeometryCore const&                  fGeom;
  detinfo::DetectorClocksService const&     fDetClocks;
  detinfo::DetectorPropertiesService const& fDetProp; 
  lar::util::TrackTimeIntervalMaker const   fTimeIntervalMaker;
};


BarycenterMatchProducer::BarycenterMatchProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  // More initializers here.
  fInputTags(p.get<std::vector<std::string>>("InputTags")),
  fOpFlashLabel(p.get<std::string>("OpFlashLabel")),
  fPandoraLabel(p.get<std::string>("PandoraLabel")),
  fTriggerLabel(p.get<std::string>("TriggerLabel")),
  fCollectionOnly(p.get<bool>("CollectionOnly")),
  fUseTimeRange(p.get<bool>("UseTimeRange")),
  fVerbose(p.get<bool>("Verbose")),
  fFillMatchTree(p.get<bool>("FillMatchTree")),
  fNominalTrigTime(p.get<double>("NominalTrigTime")),
  fTriggerTolerance(p.get<double>("TriggerTolerance")),
  fTimeRangeMargin(p.get<double>("TimeRangeMargin")),
  fGeom(*lar::providerFrom<geo::Geometry>()),
  fDetClocks(*art::ServiceHandle<detinfo::DetectorClocksService>()),
  fDetProp(*art::ServiceHandle<detinfo::DetectorPropertiesService>()),
  fTimeIntervalMaker{ fGeom }
{
  // Call appropriate produces<>() functions here.

  produces< std::vector<sbn::BarycenterMatch> >();
  produces< art::Assns<sbn::BarycenterMatch, recob::Slice> >();
  produces< art::Assns<sbn::BarycenterMatch, recob::OpFlash> >();

  // Call appropriate consumes<>() for any products to be retrieved by this module.

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

  } //End MatchTree

}

void BarycenterMatchProducer::produce(art::Event& e)
{
  // Implementation of required member function here.
  fEvent  = e.id().event();
  fRun    = e.run();
  const bool isData = e.isRealData();

  //Fetch trigger info and if MC check whether this event triggered
  art::Handle const triggerHandle
    = e.getHandle<std::vector<raw::Trigger>>(fTriggerLabel);
  double triggerWithinGate = 0.;
  bool triggeredEvt = false;
  if ( triggerHandle.isValid() && triggerHandle->size() >= 1 ) {
    raw::Trigger const& trigger = triggerHandle->at(0);
    if ( trigger.TriggerTime() >= 0 ) {
      triggerWithinGate = trigger.TriggerTime() - trigger.BeamGateTime();
      triggeredEvt = true;
    }
    if ( fVerbose ) std::cout << "Valid trigger product found. Trigger time: " << trigger.TriggerTime() << ", Beam gate time: " << trigger.BeamGateTime() << ", Difference: " << triggerWithinGate  << std::endl;
  }
  else if ( fVerbose ) std::cout << "No valid trigger product found for this event!"  << std::endl;

  //Infrastructure for checking allowed time range of a slice
  microseconds margin(fTimeRangeMargin);
  detinfo::DetectorTimings const detTimings{ fDetClocks.DataFor(e) };
  detinfo::DetectorPropertiesData const& detProp { fDetProp.DataFor(e, detTimings.clockData()) };
  lar::util::TrackTimeInterval const timeIntervals = fTimeIntervalMaker(detProp, detTimings);

  //Initialize new data products
  auto matchInfoVector = std::make_unique< std::vector<sbn::BarycenterMatch> >();
  art::PtrMaker< sbn::BarycenterMatch > const makeInfoPtr(e); 
  auto sliceAssns = std::make_unique< art::Assns<sbn::BarycenterMatch, recob::Slice> >();
  auto flashAssns = std::make_unique< art::Assns<sbn::BarycenterMatch, recob::OpFlash> >();

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

    double triggerFlashCenter = -9999.;
    double flashTime_Trigger;
    //For flash...
    for ( int i = 0; i < nFlashes; i++ ) {
      const recob::OpFlash &flash = (*flashHandle).at(i);

      //Is this a triggering flash?
      flashTime_Trigger = flash.Time();
      if ( !isData ) flashTime_Trigger += -1.*triggerWithinGate;
      if ( triggeredEvt && abs(flashTime_Trigger - fNominalTrigTime) < fTriggerTolerance ) triggerFlashCenter = flash.ZCenter();

    } //End for flash

    if ( fVerbose ) std::cout << "Event: " << fEvent << ", Cryo: " << inputTag << ", nFlashes: " << nFlashes << ", Triggering flash center: " << triggerFlashCenter  << std::endl;


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
    std::vector<art::Ptr<recob::Slice>> sliceVector;
    art::fill_ptr_vector(sliceVector, sliceHandle);

    int nSlices = (*sliceHandle).size();

    //For slice...
    for ( int j = 0; j < nSlices; j++ ) {
      fSliceNum = j;
      const art::Ptr<recob::Slice>& slicePtr = sliceVector.at(j);
      InitializeSlice();
      sbn::BarycenterMatch sliceMatchInfo;
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
      TVector3 thisPoint, thisPointSqr;
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
        thisPoint = point.XYZ();
        thisPointSqr = {thisPoint[0]*thisPoint[0], thisPoint[1]*thisPoint[1], thisPoint[2]*thisPoint[2]};
        sumCharge += thisCharge;
        sumPos += thisPoint * thisCharge;
        sumPosSqr += thisPointSqr * thisCharge;
      } //End for hit

      //No charge found in slice...
      if ( sumCharge == 0. ) {
        if ( fFillMatchTree ) fMatchTree->Fill();
        art::Ptr<sbn::BarycenterMatch> const infoPtr = makeInfoPtr(matchInfoVector->size());
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
          if ( !timeRange.contains(eTime, margin) ) continue;
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
        art::Ptr<sbn::BarycenterMatch> const infoPtr = makeInfoPtr(matchInfoVector->size());
        sliceAssns->addSingle(infoPtr, slicePtr);
        matchInfoVector->push_back(std::move(sliceMatchInfo));
        if ( fVerbose ) std::cout << "No matching flash found for Event: " << fEvent << " Slice: " << j << "! Continuing..."  << std::endl;
        continue;
      }

      //Best match flash pointer
      const art::Ptr<recob::OpFlash> flashPtr { flashHandle, matchIndex };

      //Find time of first OpHit in matched flash
      const std::vector<recob::OpHit const*> &opHitsVec = fmOpHits.at(matchIndex);
      double minTime = 1e6;
      for (const recob::OpHit *opHit : opHitsVec ) { if ( opHit->PeakTime() < minTime ) minTime = opHit->PeakTime(); }

      //Update match info
      updateFlashVars(flashPtr, minTime);
      updateMatchInfo(sliceMatchInfo);
      art::Ptr<sbn::BarycenterMatch> const infoPtr = makeInfoPtr(matchInfoVector->size());
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

void BarycenterMatchProducer::InitializeSlice() {
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
} //End InitializeSlice()


double BarycenterMatchProducer::CentroidOverlap(double center1, double center2, double width1, double width2) {
  //Centroid 2 is contained within Centroid 1, so overlap is the whole Centroid 2
  if ( (center1 - width1 < center2 - width2) && (center1 + width1 > center2 + width2) ) return (2 * width2);

  //Centroid 1 is contained within Centroid 2, so overlap is the whole Centroid 1
  else if ( (center1 - width1 > center2 - width2) && (center1 + width1 < center2 + width2) ) return (2 * width1);

  double difference = center1 - center2;
  if ( center1 > center2 ) difference *= -1;
  return difference + width1 + width2;
} //End CentroidOverlap()


double BarycenterMatchProducer::CalculateAsymmetry(art::Ptr<recob::OpFlash> flash, int cryo) {
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
void BarycenterMatchProducer::updateChargeVars(double sumCharge, TVector3 sumPos, TVector3 sumPosSqr, double triggerFlashCenter) {
  fChargeCenterXGlobal = sumPos[0] / sumCharge;
  fChargeCenterXLocal = fChargeCenterXGlobal - 210.215 * (2*fCryo - 1);
  fChargeCenterY = sumPos[1] / sumCharge;
  fChargeCenterZ = sumPos[2] / sumCharge;
  fChargeWidthX = std::sqrt( sumPosSqr[0]/sumCharge - (sumPos[0]/sumCharge)*(sumPos[0]/sumCharge) );
  fChargeWidthY = std::sqrt( sumPosSqr[1]/sumCharge - (sumPos[1]/sumCharge)*(sumPos[1]/sumCharge) );
  fChargeWidthZ = std::sqrt( sumPosSqr[2]/sumCharge - (sumPos[2]/sumCharge)*(sumPos[2]/sumCharge) );
  if ( triggerFlashCenter != -9999 ) fDeltaZ_Trigger = abs(triggerFlashCenter - fChargeCenterZ);
} //End updateChargeVars()


void BarycenterMatchProducer::updateFlashVars(art::Ptr<recob::OpFlash> flash, double firstHit) {
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


void BarycenterMatchProducer::updateMatchInfo(sbn::BarycenterMatch& matchInfo) {
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
  matchInfo.deltaZ_Trigger = fDeltaZ_Trigger;
} //End updateMatchInfo()


DEFINE_ART_MODULE(BarycenterMatchProducer)