////////////////////////////////////////////////////////////////////////
// Class:       SimpleTrackAndPFParticleAna
// Plugin Type: analyzer (art v3_06_03)
// File:        SimpleTrackAndPFParticleAna_module.cc
//
// Generated at Mon Mar  1 12:55:06 2021 by Bruce Howard using cetskelgen
// from cetlib version v3_11_01.
//
// Starting by making most of what is in a filter tool 
// under development (SimpleScanFilter) into analyzer
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional Includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"

#include "TH1D.h"

class SimpleTrackAndPFParticleAna;


class SimpleTrackAndPFParticleAna : public art::EDAnalyzer {
public:
  explicit SimpleTrackAndPFParticleAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleTrackAndPFParticleAna(SimpleTrackAndPFParticleAna const&) = delete;
  SimpleTrackAndPFParticleAna(SimpleTrackAndPFParticleAna&&) = delete;
  SimpleTrackAndPFParticleAna& operator=(SimpleTrackAndPFParticleAna const&) = delete;
  SimpleTrackAndPFParticleAna& operator=(SimpleTrackAndPFParticleAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

  // fcl params
  art::InputTag fPFParticleModuleLabel; // Label for PFParticles
  art::InputTag fTrackModuleLabel;      // Label for tracks
  bool          fUseFV;                 // Use a fiducial volume defined by the Tolerances
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


SimpleTrackAndPFParticleAna::SimpleTrackAndPFParticleAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFParticleModuleLabel( p.get< art::InputTag >("PFParticleModuleLabel") ),
  fTrackModuleLabel( p.get< art::InputTag >("TrackModuleLabel") ),
  fUseFV( p.get< bool >("UseFV") ),
  fTolUsDs( p.get< double >("TolUsDs") ),
  fTolAnode( p.get< double >("TolAnode") ),
  fTolTopBot( p.get< double >("TolTopBot") ),
  fTolCathode( p.get< double >("TolCathode") )
{
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

void SimpleTrackAndPFParticleAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.

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
    return;
  }
  pfpCounter->Fill( pfparticles.size() ); 

  if( pfparticles.size()==0 ){
    // If no PFParticles then skip, but give a message
    mf::LogWarning("SimpleScanFilter") << "No PFParticles in the considered event.";
    return;
  }

  ///////////////////////////////
  // Get the vertex association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::Vertex> fmvtx(pfparticleHandle, e, fPFParticleModuleLabel);
  if( !fmvtx.isValid() ){
    // Check against validity of fmvtx (using isValid, as used in an analyzer module from SBND)
    mf::LogError("SimpleScanFilter") << "Error in validity of fmvtx. Retruning false.";
    return;
  }
  if( fmvtx.size()==0 ){
    // Check if any vertices, if not, skip event
    if( pfparticles.size()==0 )
      mf::LogWarning("SimpleScanFilter") << "Event had no vertices, but no PFParticles. Events with no PFParticles were meant to be filtered earlier... Returning false.";
    else mf::LogError("SimpleScanFilter") << "Event had no vertices, but had PFParticles. Returning false.";
    return;
  }

  ///////////////////////////////
  // Get the track association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::Track> fmtrk(pfparticleHandle, e, fTrackModuleLabel);
  if( !fmtrk.isValid() ){
    // Check against validity of fmtrk (using isValid, as used in an analyzer module from SBND)
    mf::LogError("SimpleScanFilter") << "Error in validity of fmtrk. Returning false.";
    return;
  }

  ///////////////////////////////
  // Get geometry
  /////////////////////////////// 
  art::ServiceHandle<geo::Geometry> geom;

  ///////////////////////////////
  // Loop through PFParticles and do analysis
  ///////////////////////////////
  float maxTrkLength = -5.;
  float maxTrkLengthX = -5.;
  for( auto const& iPfp : pfparticles ){
    // For now look at track-like PFParticles
    if(abs(iPfp->PdgCode())!=13 ) continue;

    // VERTEX
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

    vtxDetX->Fill( vtxX );
    vtxDetY->Fill( vtxY );
    vtxDetZ->Fill( vtxZ );

    // IF CHECKING FIDUCIAL VOLUME
    // get positions to check against
    if( fUseFV ){
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

      // check that the vertex is within allowed tolerances
      if( !( posUsZ>0.5 ? (posUsZ-vtxZ)>fTolUsDs : true ) || 
	  !( posDsZ<-0.5 ? (vtxZ-posDsZ)>fTolUsDs : true ) || 
	  !( vtxTPC->DetectDriftDirection()==1 ? ( (posAnode-vtxX)>fTolAnode && (vtxX-posCathode)>fTolCathode ) : 
	     ( (posCathode-vtxX)>fTolCathode && (vtxX-posAnode)>fTolAnode ) ) ||
	  !( (posTopY-vtxY)>fTolTopBot ) ||
	  !( (vtxY-posBotY)>fTolTopBot ) )
	continue;
    }

    // TRACK
    // -------------------
    // access the associated track
    std::vector< art::Ptr<recob::Track> > pfpTrk = fmtrk.at(iPfp.key());
    if( pfpTrk.size()!=1 ){
      mf::LogWarning("SimpleScanFilter") << "Zero/multiple tracks found for given PFParticle (after checking it was a track). Skipping this PFParticle but take note...";
      continue;
    }
    double trkLength = pfpTrk.at(0)->Length();

    std::cout << "length = " << trkLength 
	      << " ... Dir X = " << pfpTrk.at(0)->StartDirection().Unit().X() 
	      << " ... Dir Y = " << pfpTrk.at(0)->StartDirection().Unit().Y() << std::endl;

    trkLengthAll->Fill( trkLength );
    trkLengthAllX->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X()) );
    trkLengthAllY->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().Y()) );
    if( trkLength > maxTrkLength ){
      maxTrkLength = trkLength;
      maxTrkLengthX = trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X());
    }
  }

  trkMaxLength->Fill( maxTrkLength );
  trkMaxLengthX->Fill( maxTrkLengthX );

  return;
}

DEFINE_ART_MODULE(SimpleTrackAndPFParticleAna)
