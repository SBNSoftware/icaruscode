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

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"

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
  bool          fPrintAnodeToCathode;   // Print out anode-cathode crossing candidates
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
  TH2D *vtxDetZX;
  TH2D *vtxDetXY;

  TH1D *trkMaxLength;
  TH1D *trkMaxLengthX; // projected length using direction and length
  TH1D *trkLengthAll;
  TH1D *trkLengthAllX; // projected length using direction and length
  TH1D *trkLengthAllY; // projected length using direction and length

  TH1D *trkStartX;
  TH1D *trkEndX;
  TH2D *trkStartEndX;
  TH1D *trkDistX;      // actual End_x - Start_x
  TH1D *trkDistY;      // actual End_y - Start_y
  TH1D *trkDistZ;      // actual End_z - Start_z
};


SimpleTrackAndPFParticleAna::SimpleTrackAndPFParticleAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFParticleModuleLabel( p.get< art::InputTag >("PFParticleModuleLabel") ),
  fTrackModuleLabel( p.get< art::InputTag >("TrackModuleLabel") ),
  fPrintAnodeToCathode( p.get< bool >("PrintAnodeToCathode") ),
  fUseFV( p.get< bool >("UseFV") ),
  fTolUsDs( p.get< double >("TolUsDs") ),
  fTolAnode( p.get< double >("TolAnode") ),
  fTolTopBot( p.get< double >("TolTopBot") ),
  fTolCathode( p.get< double >("TolCathode") )
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  art::ServiceHandle<art::TFileService> tfs;
  tfs->mkdir("SimpleTrackAndPFParticleAna");
  pfpCounter = tfs->make<TH1D>("SimpleAnalysis_pfpCounter",";PFParticles;Events",1000,0,1000); //300,0,300);
  vtxDetX = tfs->make<TH1D>("SimpleAnalysis_vtxDetX",";Detector X [cm];Events",100,-500,500);
  vtxDetY = tfs->make<TH1D>("SimpleAnalysis_vtxDetY",";Detector Y [cm];Events",40,-200,200);
  vtxDetZ = tfs->make<TH1D>("SimpleAnalysis_vtxDetZ",";Detector Z [cm];Events",100,-1000,1000);
  vtxDetZX = tfs->make<TH2D>("SimpleAnalysis_vtxDetZX",";Detector Z [cm];Detector X [cm]",100,-1000,1000,100,-500,500);
  vtxDetXY = tfs->make<TH2D>("SimpleAnalysis_vtxDetXY",";Detector X [cm];Detector Y [cm]",100,-500,500,40,-200,200);

  trkMaxLength = tfs->make<TH1D>("SimpleAnalysis_trkMaxLength",";Longest track length [cm];Events",503,-6,1000);
  trkMaxLengthX = tfs->make<TH1D>("SimpleAnalysis_trkMaxLengthX",";Longest track X (projected) length [cm];Events",503,-6,1000);
  trkLengthAll = tfs->make<TH1D>("SimpleAnalysis_trkLengthAll",";All considered track lengths [cm];Events",503,-6,1000);
  trkLengthAllX = tfs->make<TH1D>("SimpleAnalysis_trkLengthAllX",";All considered track X (projected) lengths [cm];Events",503,-6,1000);
  trkLengthAllY = tfs->make<TH1D>("SimpleAnalysis_trkLengthAllY",";All considered track Y (projected) lengths [cm];Events",503,-6,1000);

  trkStartX = tfs->make<TH1D>("SimpleAnalysis_trkStartX",";Track Start X [cm];Events",100,-500,500);
  trkEndX = tfs->make<TH1D>("SimpleAnalysis_trkEndX",";Track End X [cm];Events",100,-500,500);
  trkStartEndX = tfs->make<TH2D>("SimpleAnalysis_trkStartEndX",";Track Start X [cm];Track End X",100,-500,500,100,-500,500);
  trkDistX = tfs->make<TH1D>("SimpleAnalysis_trkDistX",";Track Dist X [cm];Events",252,-8,1000);
  trkDistY = tfs->make<TH1D>("SimpleAnalysis_trkDistY",";Track Dist Y [cm];Events",252,-8,1000);
  trkDistZ = tfs->make<TH1D>("SimpleAnalysis_trkDistZ",";Track Dist Z [cm];Events",502,-8,2000);
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
    mf::LogWarning("SimpleTrackAndPFParticleAna") << "Event failed to find recob::PFParticle. Are you sure you gave the right label and are running this at the right time?";
    return;
  }
  pfpCounter->Fill( pfparticles.size() );

  std::cout << "pfparticles.size() = " << pfparticles.size() << std::endl;

  if( pfparticles.size()==0 ){
    // If no PFParticles then skip, but give a message
    mf::LogWarning("SimpleTrackAndPFParticleAna") << "No PFParticles in the considered event.";
    return;
  }

  ///////////////////////////////
  // Get the vertex association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::Vertex> fmvtx(pfparticleHandle, e, fPFParticleModuleLabel);
  if( !fmvtx.isValid() ){
    // Check against validity of fmvtx (using isValid, as used in an analyzer module from SBND)
    mf::LogError("SimpleTrackAndPFParticleAna") << "Error in validity of fmvtx. Retruning false.";
    return;
  }
  if( fmvtx.size()==0 ){
    // Check if any vertices, if not, skip event
    if( pfparticles.size()==0 )
      mf::LogWarning("SimpleTrackAndPFParticleAna") << "Event had no vertices, but no PFParticles. Events with no PFParticles were meant to be filtered earlier... Returning false.";
    else mf::LogError("SimpleTrackAndPFParticleAna") << "Event had no vertices, but had PFParticles. Returning false.";
    return;
  }

  ///////////////////////////////
  // Get the track association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::Track> fmtrk(pfparticleHandle, e, fTrackModuleLabel);
  if( !fmtrk.isValid() ){
    // Check against validity of fmtrk (using isValid, as used in an analyzer module from SBND)
    mf::LogError("SimpleTrackAndPFParticleAna") << "Error in validity of fmtrk. Returning false.";
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
      mf::LogWarning("SimpleTrackAndPFParticleAna") << "Zero/multiple vertices (" << pfpVtx.size() << ") found for given PFParticle! Skipping this PFParticle but take note...";
      continue;
    }
    double vtxX = pfpVtx.at(0)->position().X();
    double vtxY = pfpVtx.at(0)->position().Y();
    double vtxZ = pfpVtx.at(0)->position().Z();

    vtxDetX->Fill( vtxX );
    vtxDetY->Fill( vtxY );
    vtxDetZ->Fill( vtxZ );
    vtxDetZX->Fill( vtxZ, vtxX );
    vtxDetXY->Fill( vtxX, vtxY );

    // IF CHECKING FIDUCIAL VOLUME
    // get positions to check against
    if( fUseFV ){
      auto vtxTPC = geom->PositionToTPCptr( pfpVtx.at(0)->position() ); // tpc containing the vertex
      if( vtxTPC == nullptr ){
	mf::LogDebug("SimpleTrackAndPFParticleAna") << "geom->PositionToTPCptr not valid. Vertex: (" << vtxX << "," << vtxY << "," << vtxZ <<  "). Skipping this PFParticle";
	continue;
      }
      double posAnode = (vtxTPC->DetectDriftDirection()==1 ? vtxTPC->ActiveBoundingBox().MaxX() : vtxTPC->ActiveBoundingBox().MinX()); // if drift +x, then anode = x max
      double posCathode = (vtxTPC->DetectDriftDirection()==-1 ? vtxTPC->ActiveBoundingBox().MaxX() : vtxTPC->ActiveBoundingBox().MinX()); // opposite of anode
      double posTopY = vtxTPC->ActiveBoundingBox().MaxY(); // top y of AV
      double posBotY = vtxTPC->ActiveBoundingBox().MinY(); // bottom y of AV
      double posUsZ = vtxTPC->ActiveBoundingBox().MinZ(); // upstream z
      double posDsZ = vtxTPC->ActiveBoundingBox().MaxZ(); // downstream z

      // check that the vertex is within allowed tolerances
      if( !( posDsZ>0.5 ? (posDsZ-vtxZ)>fTolUsDs : true ) ||
	  !( posUsZ<-0.5 ? (vtxZ-posUsZ)>fTolUsDs : true ) ||
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
      mf::LogWarning("SimpleTrackAndPFParticleAna") << "Zero/multiple (" << pfpTrk.size() << ") tracks found for PFParticle (after checking it was a track). Skipping this PFParticle but take note...";
      continue;
    }
    double trkLength = pfpTrk.at(0)->Length();
    double trkX0 = pfpTrk.at(0)->Start().X();
    double trkXf = pfpTrk.at(0)->End().X();
    double trkY0 = pfpTrk.at(0)->Start().Y();
    double trkYf = pfpTrk.at(0)->End().Y();
    double trkZ0 = pfpTrk.at(0)->Start().Z();
    double trkZf = pfpTrk.at(0)->End().Z();

    //    std::cout << "length = " << trkLength
    //	      << " ... Dir X = " << pfpTrk.at(0)->StartDirection().Unit().X()
    //	      << " ... Dir Y = " << pfpTrk.at(0)->StartDirection().Unit().Y() << std::endl;

    trkLengthAll->Fill( trkLength );
    trkLengthAllX->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X()) );
    trkLengthAllY->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().Y()) );

    trkStartX->Fill( trkX0 );
    trkEndX->Fill( trkXf );
    trkStartEndX->Fill( trkX0, trkXf );
    trkDistX->Fill( fabs(trkXf-trkX0) );
    trkDistY->Fill( fabs(trkYf-trkY0) );
    trkDistZ->Fill( fabs(trkZf-trkZ0) );

    if( trkLength > maxTrkLength ){
      maxTrkLength = trkLength;
      maxTrkLengthX = trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X());
    }

    if( fPrintAnodeToCathode && fabs(trkXf-trkX0)>=142.5 ){
      mf::LogInfo("SimpleTrackAndPFParticleAna") << "Potential anode-to-cathode crossing track with length (" << fabs(trkXf-trkX0)
						 << ") and start point (" << trkX0 << ", " << trkY0 << ", " << trkZ0 << ")";
    }
  }

  trkMaxLength->Fill( maxTrkLength );
  trkMaxLengthX->Fill( maxTrkLengthX );

  return;
}

DEFINE_ART_MODULE(SimpleTrackAndPFParticleAna)
