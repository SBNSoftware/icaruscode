////////////////////////////////////////////////////////////////////////
// Class:       SimpleScanFilter
// Plugin Type: filter (art v3_06_03)
// File:        SimpleScanFilter_module.cc
//
// Generated at Tue Feb 23 11:08:33 2021 by Bruce Howard using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art (FindMany{P}, TFileService, etc.)
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Includes for PFParticle, Vertex, Track, Geometry
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// ROOT
#include "TH1D.h"
#include "TH2D.h"

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

class SimpleScanFilter;


class SimpleScanFilter : public art::EDFilter {
public:
  explicit SimpleScanFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleScanFilter(SimpleScanFilter const&) = delete;
  SimpleScanFilter(SimpleScanFilter&&) = delete;
  SimpleScanFilter& operator=(SimpleScanFilter const&) = delete;
  SimpleScanFilter& operator=(SimpleScanFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.

  // fcl params
  art::InputTag fPFParticleModuleLabel; // Label for PFParticles
  art::InputTag fTrackModuleLabel;      // Label for tracks
  bool          fCompleteLoop;          // Sets whether to end the loop early or make debug plots
  bool          fContainTracks;         // Sets whether to check for track containment
  double        fTolUsDs;               // Upstream/Downstream tolerance for vertex
  double        fTolAnode;              // Anode tolerance for vertex
  double        fTolTopBot;             // Top/Bottom tolerance for vertex
  double        fTolCathode;            // Cathode tolerance for vertex
  double        fTolContainTrack;       // Tolerance for track containment on all sides

  // histograms
  TH1D *pfpCounter;
  TH1D *pfpNuCounter;
  TH1D *vtxDetX;
  TH1D *vtxDetY;
  TH1D *vtxDetZ;
  TH1D *trkMaxLength;
  TH1D *trkMaxLengthX;
  TH1D *trkLengthAll;
  TH1D *trkLengthAllX;
  TH1D *trkLengthAllY;
  TH1D *trkMaxDirZ;
  TH2D *trkMaxLengthVsTrkMaxDirZ;

  // Active volume map - will be filled out in constructor
  std::map< geo::TPCID, short int > mapDriftDir;
  std::map< geo::TPCID, double >    mapPosAnode;
  std::map< geo::TPCID, double >    mapPosCathode;
  std::map< geo::TPCID, double >    mapPosTopY;
  std::map< geo::TPCID, double >    mapPosBotY;
  std::map< geo::TPCID, double >    mapPosUsZ;
  std::map< geo::TPCID, double >    mapPosDsZ;

  // TPC min and max combined for containment volume
  double tpcContXMin;
  double tpcContXMax;
  double tpcContYMin;
  double tpcContYMax;
  double tpcContZMin;
  double tpcContZMax;

  // Services
  geo::GeometryCore const* fGeom;
};


SimpleScanFilter::SimpleScanFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fPFParticleModuleLabel( p.get< art::InputTag >("PFParticleModuleLabel") ),
    fTrackModuleLabel( p.get< art::InputTag >("TrackModuleLabel") ),
    fCompleteLoop( p.get< bool >("CompleteLoop") ),
    fContainTracks( p.get< bool >("ContainTracks") ),
    fTolUsDs( p.get< double >("TolUsDs") ),
    fTolAnode( p.get< double >("TolAnode") ),
    fTolTopBot( p.get< double >("TolTopBot") ),
    fTolCathode( p.get< double >("TolCathode") ),
    fTolContainTrack( p.get< double >("TolContainTrack") )
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  art::ServiceHandle<art::TFileService> tfs;
  tfs->mkdir("SimpleScanFilter");
  pfpCounter = tfs->make<TH1D>("SimpleScan_pfpCounter",";PFParticles;Events",1000,0,1000);
  pfpNuCounter = tfs->make<TH1D>("SimpleScan_pfpNuCounter",";PFP Neutrinos;Events",1000,0,1000);
  vtxDetX = tfs->make<TH1D>("SimpleScan_vtxDetX",";Detector X [cm];Events",100,-500,500);
  vtxDetY = tfs->make<TH1D>("SimpleScan_vtxDetY",";Detector Y [cm];Events",40,-200,200);
  vtxDetZ = tfs->make<TH1D>("SimpleScan_vtxDetZ",";Detector Z [cm];Events",100,-1000,1000);
  trkMaxLength = tfs->make<TH1D>("SimpleScan_trkMaxLength",";Longest track length, vertex contained [cm];Events",503,-6,1000);
  trkMaxLengthX = tfs->make<TH1D>("SimpleScan_trkMaxLengthX",";Longest track X length, vertex contained [cm];Events",503,-6,1000);
  trkLengthAll = tfs->make<TH1D>("SimpleScan_trkLengthAll",";All track lengths, vertex contained [cm];Events",503,-6,1000);
  trkLengthAllX = tfs->make<TH1D>("SimpleScan_trkLengthAllX",";All track X lengths, vertex contained [cm];Events",503,-6,1000);
  trkLengthAllY = tfs->make<TH1D>("SimpleScan_trkLengthAllY",";All track Y lengths, vertex contained [cm];Events",503,-6,1000);
  trkMaxDirZ = tfs->make<TH1D>("SimpleScan_trkMaxDirZ",";Z component of direction - longest accepted track;Events",43,-1.15,1);
  trkMaxLengthVsTrkMaxDirZ = tfs->make<TH2D>("SimpleScan_trkMaxLengthVsTrkMaxDirZ",";Z component of direction - longest accepted track;Track Length",9,-1.25,1,9,-50,400);

  // Set up geometry
  //  art::ServiceHandle<geo::Geometry> geom;
  fGeom = lar::providerFrom<geo::Geometry>();

  for( geo::TPCID const& thisTPC : fGeom->IterateTPCIDs() ){
    auto thisTPCPtr = fGeom->TPCPtr(thisTPC);
    short int thisDriftDir = thisTPCPtr->DetectDriftDirection();

    mapDriftDir[ thisTPC ]   = thisDriftDir;
    mapPosAnode[ thisTPC ]   = (thisDriftDir==1 ? thisTPCPtr->ActiveBoundingBox().MaxX() : thisTPCPtr->ActiveBoundingBox().MinX()); // if drift +x, then anode = x max
    mapPosCathode[ thisTPC ] = (thisDriftDir==-1 ? thisTPCPtr->ActiveBoundingBox().MaxX() : thisTPCPtr->ActiveBoundingBox().MinX()); // opposite of anode
    mapPosTopY[ thisTPC ]    = thisTPCPtr->ActiveBoundingBox().MaxY(); // top y of AV
    mapPosBotY[ thisTPC ]    = thisTPCPtr->ActiveBoundingBox().MinY(); // bottom y of AV
    mapPosUsZ[ thisTPC ]     = thisTPCPtr->ActiveBoundingBox().MinZ(); // upstream z
    mapPosDsZ[ thisTPC ]     = thisTPCPtr->ActiveBoundingBox().MaxZ(); // downstream z
  }

  // this time use the map to get the min and max x,y,z of the TPCs
  // TODO: probably a better way to do this...
  int counterTPC = 0;
  for( geo::TPCID const& thisTPC : fGeom->IterateTPCIDs() ){
    if(counterTPC==0){
      tpcContXMin = std::min( mapPosAnode[thisTPC], mapPosCathode[thisTPC] );
      tpcContXMax = std::max( mapPosAnode[thisTPC], mapPosCathode[thisTPC] );
      tpcContZMax = mapPosDsZ[thisTPC];
      tpcContZMin = mapPosUsZ[thisTPC];
      tpcContYMax = mapPosTopY[thisTPC];
      tpcContYMin = mapPosBotY[thisTPC];
    }
    else{
      if( mapPosAnode[thisTPC] < tpcContXMin )        tpcContXMin = mapPosAnode[thisTPC];
      else if( mapPosAnode[thisTPC] > tpcContXMax )   tpcContXMax = mapPosAnode[thisTPC];

      if( mapPosCathode[thisTPC] < tpcContXMin )      tpcContXMin = mapPosCathode[thisTPC];
      else if( mapPosCathode[thisTPC] > tpcContXMax ) tpcContXMax = mapPosCathode[thisTPC];

      if( mapPosTopY[thisTPC] < tpcContYMin )         tpcContYMin = mapPosTopY[thisTPC];
      else if( mapPosTopY[thisTPC] > tpcContYMax )    tpcContYMax = mapPosTopY[thisTPC];

      if( mapPosBotY[thisTPC] < tpcContYMin )         tpcContYMin = mapPosBotY[thisTPC];
      else if( mapPosBotY[thisTPC] > tpcContYMax )    tpcContYMax = mapPosBotY[thisTPC];

      if( mapPosDsZ[thisTPC] < tpcContZMin )          tpcContZMin = mapPosDsZ[thisTPC];
      else if( mapPosDsZ[thisTPC] > tpcContZMax )     tpcContZMax = mapPosDsZ[thisTPC];

      if( mapPosUsZ[thisTPC] < tpcContZMin )          tpcContZMin = mapPosUsZ[thisTPC];
      else if( mapPosUsZ[thisTPC] > tpcContZMax )     tpcContZMax = mapPosUsZ[thisTPC];
    }
    counterTPC+=1;
  }
}

bool SimpleScanFilter::filter(art::Event& e)
{
  ///////////////////////////////
  // Get the PFParticles
  ///////////////////////////////
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  //std::vector< art::Ptr<recob::PFParticle> > pfparticles;
  lar_pandora::PFParticleVector pfparticles;
  if( e.getByLabel(fPFParticleModuleLabel,pfparticleHandle) ){
    art::fill_ptr_vector(pfparticles,pfparticleHandle);
  }
  else{
    mf::LogWarning("SimpleScanFilter") << "Event failed to find recob::PFParticle. Are you sure you gave the right label and are running this at the right time?";
    return false;
  }

  if( fCompleteLoop ) pfpCounter->Fill( pfparticles.size() ); 
  if( pfparticles.size()==0 ){
    // If no PFParticles then skip, but give a message
    mf::LogWarning("SimpleScanFilter") << "No PFParticles in the considered event.";
    return false;
  }

  // PFParticle Map
  lar_pandora::PFParticleMap pfpMap;
  lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticles,pfpMap);

  ///////////////////////////////
  // Get the vertex association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::Vertex> fmvtx(pfparticleHandle, e, fPFParticleModuleLabel);
  if( !fmvtx.isValid() ){
    // Check against validity of fmvtx (using isValid, as used in an analyzer module from SBND)
    mf::LogError("SimpleScanFilter") << "Error in validity of fmvtx. Retruning false.";
    return false;
  }
  if( fmvtx.size()==0 ){
    // Check if any vertices, if not, skip event
    if( pfparticles.size()==0 )
      mf::LogWarning("SimpleScanFilter") << "Event had no vertices, but no PFParticles. Events with no PFParticles were meant to be filtered earlier... Returning false.";
    else mf::LogError("SimpleScanFilter") << "Event had no vertices, but had PFParticles. Returning false.";
    return false;
  }

  ///////////////////////////////
  // Get the track association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::Track> fmtrk(pfparticleHandle, e, fTrackModuleLabel);
  if( !fmtrk.isValid() ){
    // Check against validity of fmtrk (using isValid, as used in an analyzer module from SBND)
    mf::LogError("SimpleScanFilter") << "Error in validity of fmtrk. Returning false.";
    return false;
  }

  ///////////////////////////////
  // Get geometry
  ///////////////////////////////
  fGeom = lar::providerFrom<geo::Geometry>();

  ///////////////////////////////
  // Loop through the PFPs and do the checks
  ///////////////////////////////
  // Note: If efficiency matters, then probably you want to return true as soon as you find a
  // a track of >50cm meeting the vertex criteria. But, to make some debug plots, then set
  //   CompleteLoop: true
  // in the fcl settings and it will find the maximum track length for tracks with vertices
  // meeting the criteria, and will fill additional stats in the vertex plots.

  double maxTrkLengthVal = -5.;
  double maxTrkLengthXVal = -5.;
  double maxTrkDirZVal = -1.14;
  int pfpNus = 0;

  for( auto const& iPfp : pfparticles ){
    // FOR CHECKING JUST TRACKS
    //    if(abs(iPfp->PdgCode())!=13 ) continue;

    // CHECKING PFP NUs
    if( abs(iPfp->PdgCode()!=12 && iPfp->PdgCode()!=14) ) continue;
    pfpNus+=1;

    // debug print out useful to look at pfparticle hierarchies
    mf::LogDebug("SimpleScanFilter") << "PFParticle Nu with ID: " << iPfp.key();

    // CHECK VERTEX
    // -------------------
    // access the associated vertex
    std::vector< art::Ptr<recob::Vertex> > pfpVtx = fmvtx.at(iPfp.key());
    if( pfpVtx.size()!=1 ){
      mf::LogWarning("SimpleScanFilter") << "Zero/multiple vertices (" << pfpVtx.size() << ") found for given PFParticle! Skipping this PFParticle but take note...";
      continue;
    }
    double vtxX = pfpVtx.at(0)->position().X();
    double vtxY = pfpVtx.at(0)->position().Y();
    double vtxZ = pfpVtx.at(0)->position().Z();
    // debug plots
    if( fCompleteLoop ){
      vtxDetX->Fill( vtxX );
      vtxDetY->Fill( vtxY );
      vtxDetZ->Fill( vtxZ );
    }

    // debug print out useful to look at pfparticle hierarchies
    mf::LogDebug("SimpleScanFilter") << "PFParticle Nu Vtx: " << vtxX << ", " << vtxY << ", " << vtxZ;

    // check against the positions in the map and the tolerances
    auto thisVtxTPC = fGeom->PositionToTPCID( pfpVtx.at(0)->position() );
    if( !thisVtxTPC.isValid ){
      // TODO: is this doing anything...
      mf::LogWarning("SimpleScanFilter") << "fGeom->PositionToTPCID not valid. Vertex: (" << vtxX << "," << vtxY << "," << vtxZ <<  "). Skipping this PFParticle";
      continue; 
    }

    // Us/Ds fixed here I think... need to propagate to the analyzer module too...
    if( !( mapPosDsZ[thisVtxTPC]>0.5 ? (mapPosDsZ[thisVtxTPC]-vtxZ)>fTolUsDs : true ) ||
        !( mapPosUsZ[thisVtxTPC]<-0.5 ? (vtxZ-mapPosUsZ[thisVtxTPC])>fTolUsDs : true ) ||
        !( mapDriftDir[thisVtxTPC]==1 ? ( (mapPosAnode[thisVtxTPC]-vtxX)>fTolAnode && (vtxX-mapPosCathode[thisVtxTPC])>fTolCathode ) :
	   ( (mapPosCathode[thisVtxTPC]-vtxX)>fTolCathode && (vtxX-mapPosAnode[thisVtxTPC])>fTolAnode ) ) ||
        !( (mapPosTopY[thisVtxTPC]-vtxY)>fTolTopBot ) ||
        !( (vtxY-mapPosBotY[thisVtxTPC])>fTolTopBot ) )
      continue;

    // FOR PFP NUs, loop through the daughters and require at least one track of >50cm (could also try 2...)
    for( auto const& iDaughter : iPfp->Daughters() ){
      // TODO: Probably want to add a check and a throw here if it can't find the daughter (not just an error or warning, 
      //since that would not be good I think). See how it was implemented in the bool lar_pandora::LArPandoraHelper::IsFinalState function in doxygem
      art::Ptr< recob::PFParticle > iPfpDaughter = pfpMap.find(iDaughter)->second;

      // debug print out useful to look at pfparticle hierarchies
      mf::LogDebug("SimpleScanFilter") << "PFParticle Child with ID: " << iDaughter << "(assumption check: " << iPfpDaughter.key() << ", " << iPfpDaughter->Self() << ")" << "\n"
				       << "                     PDG: " << iPfpDaughter->PdgCode();

      if( abs(iPfpDaughter->PdgCode())!=13 ) continue; // skip non-track PFPs for now

      //////////////////////////////////////////
      // NOTE: When looping through tracks only and not the PFP Nus, then the following 
      //       lines should be back with the main loop through PFParticles
      //// START
      //////////////////////////////////////////

      // CHECK TRACK
      // -------------------
      // access the associated track, if PFParticle's PDG Code == +/- 13
      std::vector< art::Ptr<recob::Track> > pfpTrk = fmtrk.at(iPfpDaughter.key()); // for putting back in main pfp loop, iPfpDaughter.key() -> iPfp.key()
      if( pfpTrk.size()!=1 ){
	mf::LogWarning("SimpleScanFilter") << "Zero/multiple tracks (" << pfpTrk.size() << ") found for given PFParticle (after checking it was a track). Skipping this PFParticle but take note...";
	continue;
      }

      // debug print out useful to look at pfparticle hierarchies
      mf::LogDebug("SimpleScanFilter") << "PFParticle Child with ID: " << iDaughter << "(assumption check: " << iPfpDaughter.key() << ")" << "\n"
				       << "            Track Length: " << pfpTrk.at(0)->Length();

      // check track containment -- currently hard coded 2cm outside boundaries
      if( fContainTracks ){
	double trkStartX = pfpTrk.at(0)->Start().X();
	double trkEndX = pfpTrk.at(0)->End().X();
	double trkStartY = pfpTrk.at(0)->Start().Y();
	double trkEndY = pfpTrk.at(0)->End().Y();
	double trkStartZ = pfpTrk.at(0)->Start().Z();
	double trkEndZ = pfpTrk.at(0)->End().Z();

	// check start
	if( !( (tpcContXMax-trkStartX > fTolContainTrack) && (trkStartX-tpcContXMin > fTolContainTrack) && 
	       (tpcContYMax-trkStartY > fTolContainTrack) && (trkStartY-tpcContYMin > fTolContainTrack) &&
	       (tpcContZMax-trkStartZ > fTolContainTrack) && (trkStartZ-tpcContZMin > fTolContainTrack) ) )
	  continue;
	// check end
	if( !( (tpcContXMax-trkEndX > fTolContainTrack) && (trkEndX-tpcContXMin > fTolContainTrack) &&
	       (tpcContYMax-trkEndY > fTolContainTrack) && (trkEndY-tpcContYMin > fTolContainTrack) &&
	       (tpcContZMax-trkEndZ > fTolContainTrack) && (trkEndZ-tpcContZMin > fTolContainTrack) ) )
	  continue;
      }

      double trkLength = pfpTrk.at(0)->Length();

      if( !fCompleteLoop && trkLength>50. ) return true; // if not making debug plots and we found our criteria, then go!
      if( fCompleteLoop ){
	trkLengthAll->Fill(trkLength);
	trkLengthAllX->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X()) );
	trkLengthAllY->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().Y()) );
	if( trkLength > maxTrkLengthVal ){
	  maxTrkLengthVal  = trkLength;
	  maxTrkLengthXVal = trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X());
	  maxTrkDirZVal    = pfpTrk.at(0)->StartDirection().Unit().Z();
	}
      }
      //// END
    } // end loop through daughters of this PFPNeutrino
  } // end loop through PFPs

  // if doing fCompleteLoop then you should reach this point. Make one final decision about if the track is passing
  if( fCompleteLoop ){
    trkMaxLength->Fill( maxTrkLengthVal );
    trkMaxLengthX->Fill( maxTrkLengthXVal );
    trkMaxDirZ->Fill( maxTrkDirZVal );
    trkMaxLengthVsTrkMaxDirZ->Fill( maxTrkDirZVal, maxTrkLengthVal );
    pfpNuCounter->Fill( pfpNus );
    //    if( maxTrkLengthVal>50. ) return true;
  }

  ///////////////////////////////
  // If you get here, return false
  ///////////////////////////////
  return false;
}

DEFINE_ART_MODULE(SimpleScanFilter)
