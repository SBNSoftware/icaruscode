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

// ROOT
#include "TH1D.h"

#include <memory>

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
  double        fTolUsDs;               // Upstream/Downstream tolerance for vertex
  double        fTolAnode;              // Anode tolerance for vertex
  double        fTolTopBot;             // Top/Bottom tolerance for vertex
  double        fTolCathode;            // Cathode tolerance for vertex

  // histograms
  TH1D *pfpCounter;
  TH1D *vtxDetX;
  TH1D *vtxDetY;
  TH1D *vtxDetZ;
  TH1D *trkMaxLength;
  TH1D *trkMaxLengthX;
  TH1D *trkLengthAll;
  TH1D *trkLengthAllX;
  TH1D *trkLengthAllY;
};


SimpleScanFilter::SimpleScanFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fPFParticleModuleLabel( p.get< art::InputTag >("PFParticleModuleLabel") ),
    fTrackModuleLabel( p.get< art::InputTag >("TrackModuleLabel") ),
    fCompleteLoop( p.get< bool >("CompleteLoop") ),
    fTolUsDs( p.get< double >("TolUsDs") ),
    fTolAnode( p.get< double >("TolAnode") ),
    fTolTopBot( p.get< double >("TolTopBot") ),
    fTolCathode( p.get< double >("TolCathode") )
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  art::ServiceHandle<art::TFileService> tfs;
  tfs->mkdir("SimpleScanFilter");
  pfpCounter = tfs->make<TH1D>("SimpleScan-pfpCounter",";PFParticles;Events",300,0,300);
  vtxDetX = tfs->make<TH1D>("SimpleScan_vtxDetX",";Detector X [cm];Events",100,-500,500);
  vtxDetY = tfs->make<TH1D>("SimpleScan_vtxDetY",";Detector Y [cm];Events",40,-200,200);
  vtxDetZ = tfs->make<TH1D>("SimpleScan_vtxDetZ",";Detector Z [cm];Events",100,-1000,1000);
  trkMaxLength = tfs->make<TH1D>("SimpleScan_trkMaxLength",";Longest track length, vertex contained [cm];Events",503,-6,1000);
  trkMaxLengthX = tfs->make<TH1D>("SimpleScan_trkMaxLengthX",";Longest track X length, vertex contained [cm];Events",503,-6,1000);
  trkLengthAll = tfs->make<TH1D>("SimpleScan_trkLengthAll",";All track lengths, vertex contained [cm];Events",503,-6,1000);
  trkLengthAllX = tfs->make<TH1D>("SimpleScan_trkLengthAllX",";All track X lengths, vertex contained [cm];Events",503,-6,1000);
  trkLengthAllY = tfs->make<TH1D>("SimpleScan_trkLengthAllY",";All track Y lengths, vertex contained [cm];Events",503,-6,1000);
}

bool SimpleScanFilter::filter(art::Event& e)
{
  ///////////////////////////////
  // Get the PFParticles
  ///////////////////////////////
  art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > pfparticles;
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
  art::ServiceHandle<geo::Geometry> geom;

  ///////////////////////////////
  // Loop through the PFPs and do the checks
  ///////////////////////////////
  // Note: If efficiency matters, then probably you want to return true as soon as you find a
  // a track of >50cm meeting the vertex criteria. But, to make some debug plots, then set
  //   CompleteLoop: true 
  // in the fcl settings and it will find the maximum track length for tracks with vertices 
  // meeting the criteria, and will fill additional stats in the vertex plots.

  //  std::cout << "Size of pfparticles " << fPFParticleModuleLabel << ": " << pfparticles.size() << std::endl;
  float maxTrkLength = -5.;
  float maxTrkLengthX = -5.;
  for( auto const& iPfp : pfparticles ){
    ////    std::cout << iPfp->PdgCode() << std::endl;
    if(abs(iPfp->PdgCode())!=13 ) continue;

    // CHECK VERTEX
    // -------------------
    // access the associated vertex
    std::vector< art::Ptr<recob::Vertex> > pfpVtx = fmvtx.at(iPfp.key());
    if( pfpVtx.size()!=1 ){
      mf::LogWarning("SimpleScanFilter") << "Zero/multiple vertices found for given PFParticle! Skipping this PFParticle but take note...";
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

    // get positions to check against
    auto vtxTPC = geom->PositionToTPCptr( pfpVtx.at(0)->position() ); // tpc containing the vertex
    if( vtxTPC == nullptr ){
      mf::LogDebug("SimpleScanFilter") << "geom->PositionToTPCptr not valid. Vertex: (" << vtxX << "," << vtxY << "," << vtxZ <<  "). Skipping this PFParticle";
      continue;
    }
    double posAnode = (vtxTPC->DetectDriftDirection()==1 ? vtxTPC->ActiveBoundingBox().MaxX() : vtxTPC->ActiveBoundingBox().MinX()); // if drift +x, then anode = x max
    double posCathode = (vtxTPC->DetectDriftDirection()==-1 ? vtxTPC->ActiveBoundingBox().MaxX() : vtxTPC->ActiveBoundingBox().MinX()); // opposite of anode
    double posTopY = vtxTPC->ActiveBoundingBox().MaxY(); // top y of AV
    double posBotY = vtxTPC->ActiveBoundingBox().MinY(); // bottom y of AV
    double posUsZ = vtxTPC->ActiveBoundingBox().MinZ(); // upstream z
    double posDsZ = vtxTPC->ActiveBoundingBox().MaxZ(); // downstream z

    ////    std::cout << "TPC" << vtxTPC->TPCInfo() << std::endl;
    ////    std::cout << "Anode: " << posAnode << ", Cathode: " << posCathode << std::endl;
    ////    std::cout << "Top Y: " << posTopY << ", Bot Y: " << posBotY << std::endl;
    ////    std::cout << "Upstream Z: " << posUsZ << ", Downstream Z: " << posDsZ << std::endl;
    ////    std::cout << "" << std::endl;

    // check that the vertex is within allowed tolerances
    if( !( posUsZ>0.5 ? (posUsZ-vtxZ)>fTolUsDs : true ) || 
	!( posDsZ<-0.5 ? (vtxZ-posDsZ)>fTolUsDs : true ) || 
	!( vtxTPC->DetectDriftDirection()==1 ? ( (posAnode-vtxX)>fTolAnode && (vtxX-posCathode)>fTolCathode ) : 
	                                      ( (posCathode-vtxX)>fTolCathode && (vtxX-posAnode)>fTolAnode ) ) ||
	!( (posTopY-vtxY)>fTolTopBot ) ||
	!( (vtxY-posBotY)>fTolTopBot ) )
      continue;

    // CHECK TRACK
    // -------------------
    // access the associated track, if PFParticle's PDG Code == +/- 13
    std::vector< art::Ptr<recob::Track> > pfpTrk = fmtrk.at(iPfp.key());
    ////    std::cout << pfpTrk.size() << std::endl;
    if( pfpTrk.size()!=1 ){
      mf::LogWarning("SimpleScanFilter") << "Zero/multiple tracks found for given PFParticle (after checking it was a track). Skipping this PFParticle but take note...";
      continue;
    }
    double trkLength = pfpTrk.at(0)->Length();

    std::cout << "length = " << trkLength
              << " ... Dir X = " << pfpTrk.at(0)->StartDirection().Unit().X()
              << " ... Dir Y = " << pfpTrk.at(0)->StartDirection().Unit().Y() << std::endl;

    if( !fCompleteLoop && trkLength>50. ) return true; // if not making debug plots and we found our criteria, then go!
    if( fCompleteLoop ){
      trkLengthAll->Fill(trkLength);
      trkLengthAllX->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X()) );
      trkLengthAllY->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().Y()) );
      if( trkLength > maxTrkLength ){
	maxTrkLength = trkLength;
	maxTrkLengthX = trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X());
      }
    }
  }

  // if doing fCompleteLoop then you should reach this point. Make one final decision about if the track is passing
  if( fCompleteLoop ){
    trkMaxLength->Fill( maxTrkLength );
    trkMaxLengthX->Fill( maxTrkLengthX );
    //    if( maxTrkLength>50. ) return true;
  }

  ///////////////////////////////
  // If you get here, return false
  ///////////////////////////////
  return false;
}

DEFINE_ART_MODULE(SimpleScanFilter)
