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
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcore/Geometry/Geometry.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TVector3.h"

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
  double        fTrkMinLenCut;          // Cut on minimum track length for plots (default in the analysis fcl is 0.)

  // histograms
  TH1D *pfpCounter;
  TH1D *vtxDetX;
  TH1D *vtxDetY;
  TH1D *vtxDetZ;
  TH2D *vtxDetZX;
  TH2D *vtxDetXY;

  TH1D *vtxDetXAll;
  TH1D *vtxDetYAll;
  TH1D *vtxDetZAll;

  TH1D *trkMaxLength;
  //TH1D *trkMaxLengthX; // projected length using direction and length
  TH1D *trkLengthAll;
  //TH1D *trkLengthAllX; // projected length using direction and length
  //TH1D *trkLengthAllY; // projected length using direction and length

  TH1D *pfpNumSps;
  TH1D *trkNumHits;
  TH1D *trkTheta;
  TH1D *trkPhi;
  TH2D *trkThetaPhi;
  TH1D *trkThetaXZ;
  TH1D *trkThetaYZ;
  TH2D *trkThetaXZThetaYZ;
  TH2D *trkLengthNumHits;

  TH1D *trkStartX;
  TH1D *trkStartY;
  TH1D *trkStartZ;
  TH1D *trkEndX;
  TH1D *trkEndY;
  TH1D *trkEndZ;
  TH2D *trkStartEndX;
  TH2D *trkStartEndY;
  TH2D *trkStartEndZ;
  TH2D *trkStartXLength;
  TH2D *trkStartYLength;
  TH2D *trkStartZLength;
  TH1D *trkDistX;      // actual End_x - Start_x
  TH1D *trkDistY;      // actual End_y - Start_y
  TH1D *trkDistZ;      // actual End_z - Start_z

  TH1D *trkThetaToWire0;
  TH1D *trkThetaToWire1;
  TH1D *trkThetaToWire2;
  TH1D *trkNumHits0;
  TH1D *trkNumHits1;
  TH1D *trkNumHits2;
  TH1D *trkNumWires0;
  TH1D *trkNumWires1;
  TH1D *trkNumWires2;
  TH1D *trkNumHitsPerWire0;
  TH1D *trkNumHitsPerWire1;
  TH1D *trkNumHitsPerWire2;
  TH2D *trkNumWiresNumHits0;
  TH2D *trkNumWiresNumHits1;
  TH2D *trkNumWiresNumHits2;

  // No pass requirements
  TH1D *trkStartXAll;
  TH1D *trkStartYAll;
  TH1D *trkStartZAll;

  // Pass Track Length
  TH1D *trkStartXPassTrkLen;
  TH1D *trkStartYPassTrkLen;
  TH1D *trkStartZPassTrkLen;

  // Fail the AV cut
  TH1D *trkThetaXZFailAV;
  TH1D *trkThetaYZFailAV;

  // Plots for anode to cathode crossing candidates
  TH1D *trkLengthAllA2CC;
  TH1D *trkDistXA2CC;
  TH1D *trkDistYA2CC;
  TH1D *trkDistZA2CC;
  TH1D *trkNumHitsA2CC;
  TH1D *trkThetaA2CC;
  TH1D *trkPhiA2CC;
  TH2D *trkThetaPhiA2CC;
  TH1D *trkThetaXZA2CC;
  TH1D *trkThetaYZA2CC;
  TH2D *trkThetaXZThetaYZA2CC;
  TH1D *trkStartXA2CC;
  TH1D *trkStartYA2CC;
  TH1D *trkStartZA2CC;
  TH1D *trkThetaToWireA2CC0;
  TH1D *trkThetaToWireA2CC1;
  TH1D *trkThetaToWireA2CC2;
  TH1D *trkNumHitsA2CC0;
  TH1D *trkNumHitsA2CC1;
  TH1D *trkNumHitsA2CC2;
  TH1D *trkNumWiresA2CC0;
  TH1D *trkNumWiresA2CC1;
  TH1D *trkNumWiresA2CC2;
  TH1D *trkNumHitsPerWireA2CC0;
  TH1D *trkNumHitsPerWireA2CC1;
  TH1D *trkNumHitsPerWireA2CC2;
  TH2D *trkNumWiresNumHitsA2CC0;
  TH2D *trkNumWiresNumHitsA2CC1;
  TH2D *trkNumWiresNumHitsA2CC2;

  // Plots for anode to anode crossing candidates
  TH1D *trkLengthAllA2AC;
  TH1D *trkDistXA2AC;
  TH1D *trkDistYA2AC;
  TH1D *trkDistZA2AC;
  TH1D *trkNumHitsA2AC;
  TH1D *trkThetaA2AC;
  TH1D *trkPhiA2AC;
  TH2D *trkThetaPhiA2AC;
  TH1D *trkThetaXZA2AC;
  TH1D *trkThetaYZA2AC;
  TH2D *trkThetaXZThetaYZA2AC;
  TH1D *trkStartXA2AC;
  TH1D *trkStartYA2AC;
  TH1D *trkStartZA2AC;
  TH1D *trkThetaToWireA2AC0;
  TH1D *trkThetaToWireA2AC1;
  TH1D *trkThetaToWireA2AC2;
  TH1D *trkNumHitsA2AC0;
  TH1D *trkNumHitsA2AC1;
  TH1D *trkNumHitsA2AC2;
  TH1D *trkNumWiresA2AC0;
  TH1D *trkNumWiresA2AC1;
  TH1D *trkNumWiresA2AC2;
  TH1D *trkNumHitsPerWireA2AC0;
  TH1D *trkNumHitsPerWireA2AC1;
  TH1D *trkNumHitsPerWireA2AC2;
  TH2D *trkNumWiresNumHitsA2AC0;
  TH2D *trkNumWiresNumHitsA2AC1;
  TH2D *trkNumWiresNumHitsA2AC2;
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
  fTolCathode( p.get< double >("TolCathode") ),
  fTrkMinLenCut( p.get< double >("TrkMinLenCut") )
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

  vtxDetXAll = tfs->make<TH1D>("SimpleAnalysis_vtxDetXAll",";Detector X [cm];Events",100,-500,500);
  vtxDetYAll = tfs->make<TH1D>("SimpleAnalysis_vtxDetYAll",";Detector Y [cm];Events",40,-200,200);
  vtxDetZAll = tfs->make<TH1D>("SimpleAnalysis_vtxDetZAll",";Detector Z [cm];Events",100,-1000,1000);

  trkMaxLength = tfs->make<TH1D>("SimpleAnalysis_trkMaxLength",";Longest track length [cm];Events",503,-6,1000);
  //trkMaxLengthX = tfs->make<TH1D>("SimpleAnalysis_trkMaxLengthX",";Longest track X (projected) length [cm];Events",503,-6,1000);
  trkLengthAll = tfs->make<TH1D>("SimpleAnalysis_trkLengthAll",";All considered track lengths [cm];Events",503,-6,1000);
  //trkLengthAllX = tfs->make<TH1D>("SimpleAnalysis_trkLengthAllX",";All considered track X (projected) lengths [cm];Events",503,-6,1000);
  //trkLengthAllY = tfs->make<TH1D>("SimpleAnalysis_trkLengthAllY",";All considered track Y (projected) lengths [cm];Events",503,-6,1000);

  pfpNumSps = tfs->make<TH1D>("SimpleAnalysis_pfpNumSps",";Number SpacePoints associated to PFParticle;Events",600,0,6000);
  trkNumHits = tfs->make<TH1D>("SimpleAnalysis_trkNumHits",";Number Hits associated to Track;Events",600,0,6000);
  trkTheta = tfs->make<TH1D>("SimpleAnalysis_trkTheta",";Theta [rad];Events",100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkPhi = tfs->make<TH1D>("SimpleAnalysis_trkPhi",";Phi [rad];Events",100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkThetaPhi = tfs->make<TH2D>("SimpleAnalysis_trkThetaPhi",";Theta [rad];Phi [rad]",50,-2.*TMath::Pi(),2.*TMath::Pi(),50,-2.*TMath::Pi(),2.*TMath::Pi());
  trkThetaXZ = tfs->make<TH1D>("SimpleAnalysis_trkThetaXZ",";ThetaXZ [deg];Events",100,-180,180);
  trkThetaYZ = tfs->make<TH1D>("SimpleAnalysis_trkThetaYZ",";ThetaYZ [deg];Events",100,-180,180);
  trkThetaXZThetaYZ = tfs->make<TH2D>("SimpleAnalysis_trkThetaXZThetaYZ",";ThetaXZ [deg];ThetaYZ [deg]",50,-180,180,50,-180,180);
  trkLengthNumHits = tfs->make<TH2D>("SimpleAnalysis_trkLengthNumHits",";Track Length [cm];Number Hits associated to Track",100,0,500,100,0,2000);

  trkStartX = tfs->make<TH1D>("SimpleAnalysis_trkStartX",";Track Start X [cm];Events",100,-500,500);
  trkStartY = tfs->make<TH1D>("SimpleAnalysis_trkStartY",";Track Start Y [cm];Events",40,-200,200);
  trkStartZ = tfs->make<TH1D>("SimpleAnalysis_trkStartZ",";Track Start Z [cm];Events",100,-1000,1000);
  trkEndX = tfs->make<TH1D>("SimpleAnalysis_trkEndX",";Track End X [cm];Events",100,-500,500);
  trkEndY = tfs->make<TH1D>("SimpleAnalysis_trkEndY",";Track End Y [cm];Events",40,-200,200);
  trkEndZ = tfs->make<TH1D>("SimpleAnalysis_trkEndZ",";Track End Z [cm];Events",100,-1000,1000);
  trkStartEndX = tfs->make<TH2D>("SimpleAnalysis_trkStartEndX",";Track Start X [cm];Track End X",100,-500,500,100,-500,500);
  trkStartEndY = tfs->make<TH2D>("SimpleAnalysis_trkStartEndY",";Track Start Y [cm];Track End Y",40,-200,200,40,-200,200);
  trkStartEndZ = tfs->make<TH2D>("SimpleAnalysis_trkStartEndZ",";Track Start Z [cm];Track End Z",100,-1000,1000,100,-1000,1000);
  trkStartXLength = tfs->make<TH2D>("SimpleAnalysis_trkStartXLength",";Track Start X [cm];Track Length",100,-500,500,100,-500,500);
  trkStartYLength = tfs->make<TH2D>("SimpleAnalysis_trkStartYLength",";Track Start Y [cm];Track Length",40,-200,200,40,-200,200);
  trkStartZLength = tfs->make<TH2D>("SimpleAnalysis_trkStartZLength",";Track Start Z [cm];Track Length",100,-1000,1000,100,-1000,1000);
  trkDistX = tfs->make<TH1D>("SimpleAnalysis_trkDistX",";Track Dist X [cm];Events",252,-8,1000);
  trkDistY = tfs->make<TH1D>("SimpleAnalysis_trkDistY",";Track Dist Y [cm];Events",252,-8,1000);
  trkDistZ = tfs->make<TH1D>("SimpleAnalysis_trkDistZ",";Track Dist Z [cm];Events",502,-8,2000);

  trkThetaToWire0 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWire0",";Theta To Wire Plane 0 [deg];Events",37,0,185);
  trkThetaToWire1 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWire1",";Theta To Wire Plane 1 [deg];Events",37,0,185);
  trkThetaToWire2 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWire2",";Theta To Wire Plane 2 [deg];Events",37,0,185);
  trkNumHits0 = tfs->make<TH1D>("SimpleAnalysis_trkNumHits0",";Number Plane 0 Hits associated to Track;Events",150,0,3000);
  trkNumHits1 = tfs->make<TH1D>("SimpleAnalysis_trkNumHits1",";Number Plane 1 Hits associated to Track;Events",150,0,3000);
  trkNumHits2 = tfs->make<TH1D>("SimpleAnalysis_trkNumHits2",";Number Plane 2 Hits associated to Track;Events",150,0,3000);
  trkNumWires0 = tfs->make<TH1D>("SimpleAnalysis_trkNumWires0",";Number Expected Plane 0 Wires for Track;Events",150,0,3000);
  trkNumWires1 = tfs->make<TH1D>("SimpleAnalysis_trkNumWires1",";Number Expected Plane 1 Wires for Track;Events",150,0,3000);
  trkNumWires2 = tfs->make<TH1D>("SimpleAnalysis_trkNumWires2",";Number Expected Plane 2 Wires for Track;Events",150,0,3000);
  trkNumHitsPerWire0 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWire0",";Number of Hits / Number of Wires [Plane0];Events",300,0,15);
  trkNumHitsPerWire1 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWire1",";Number of Hits / Number of Wires [Plane1];Events",300,0,15);
  trkNumHitsPerWire2 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWire2",";Number of Hits / Number of Wires [Plane2];Events",300,0,15);
  trkNumWiresNumHits0 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHits0",";Number of Wires [Plane 0];Number of Hits [Plane0]",150,0,3000,150,0,3000);
  trkNumWiresNumHits1 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHits1",";Number of Wires [Plane 1];Number of Hits [Plane1]",150,0,3000,150,0,3000);
  trkNumWiresNumHits2 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHits2",";Number of Wires [Plane 2];Number of Hits [Plane2]",150,0,3000,150,0,3000);

  // No pass requirements
  trkStartXAll = tfs->make<TH1D>("SimpleAnalysis_trkStartXAll",";Track Start X [cm];Events",100,-500,500);
  trkStartYAll = tfs->make<TH1D>("SimpleAnalysis_trkStartYAll",";Track Start Y [cm];Events",40,-200,200);
  trkStartZAll = tfs->make<TH1D>("SimpleAnalysis_trkStartZAll",";Track Start Z [cm];Events",100,-1000,1000);

  // Pass Track Length
  trkStartXPassTrkLen = tfs->make<TH1D>("SimpleAnalysis_trkStartXPassTrkLen",";Track Start X [cm];Events",100,-500,500);
  trkStartYPassTrkLen = tfs->make<TH1D>("SimpleAnalysis_trkStartYPassTrkLen",";Track Start Y [cm];Events",40,-200,200);
  trkStartZPassTrkLen = tfs->make<TH1D>("SimpleAnalysis_trkStartZPasstrkLen",";Track Start Z [cm];Events",100,-1000,1000);

  // Fail the AV cut
  trkThetaXZFailAV = tfs->make<TH1D>("SimpleAnalysis_trkThetaXZFailAV",";ThetaXZ [deg];Events",100,-180,180);
  trkThetaYZFailAV = tfs->make<TH1D>("SimpleAnalysis_trkThetaYZFailAV",";ThetaYZ [deg];Events",100,-180,180);

  // Plots for Anode To Cathode Candidates (A2C Candidates)
  trkLengthAllA2CC = tfs->make<TH1D>("SimpleAnalysis_trkLengthAllA2CC",";Considered track lengths [cm];Events",503,-6,1000);
  trkDistXA2CC = tfs->make<TH1D>("SimpleAnalysis_trkDistXA2CC",";Track Dist X [cm];Events",252,-8,1000);
  trkDistYA2CC = tfs->make<TH1D>("SimpleAnalysis_trkDistYA2CC",";Track Dist Y [cm];Events",252,-8,1000);
  trkDistZA2CC = tfs->make<TH1D>("SimpleAnalysis_trkDistZA2CC",";Track Dist Z [cm];Events",502,-8,2000);
  trkNumHitsA2CC = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2CC",";Number Hits associated to Track;Events",600,0,6000);
  trkThetaA2CC = tfs->make<TH1D>("SimpleAnalysis_trkThetaA2CC",";Theta [rad];Events",100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkPhiA2CC = tfs->make<TH1D>("SimpleAnalysis_trkPhiA2CC",";Phi [rad];Events",100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkThetaPhiA2CC = tfs->make<TH2D>("SimpleAnalysis_trkThetaPhiA2CC",";Theta [rad];Phi [rad]",100,-2.*TMath::Pi(),2.*TMath::Pi(),100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkThetaXZA2CC = tfs->make<TH1D>("SimpleAnalysis_trkThetaXZA2CC",";ThetaXZ [deg];Events",100,-180,180);
  trkThetaYZA2CC = tfs->make<TH1D>("SimpleAnalysis_trkThetaYZA2CC",";ThetaYZ [deg];Events",100,-180,180);
  trkThetaXZThetaYZA2CC = tfs->make<TH2D>("SimpleAnalysis_trkThetaXZThetaYZA2CC",";ThetaXZ [deg];ThetaYZ [deg]",100,-180,180,100,-180,180);
  trkStartXA2CC = tfs->make<TH1D>("SimpleAnalysis_trkStartXA2CC",";Track Start X [cm];Events",100,-500,500);
  trkStartYA2CC = tfs->make<TH1D>("SimpleAnalysis_trkStartYA2CC",";Track Start Y [cm];Events",40,-200,200);
  trkStartZA2CC = tfs->make<TH1D>("SimpleAnalysis_trkStartZA2CC",";Track Start Z [cm];Events",100,-1000,1000);
  trkThetaToWireA2CC0 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWireA2CC0",";Theta To Wire Plane 0 [deg];Events",37,0,185);
  trkThetaToWireA2CC1 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWireA2CC1",";Theta To Wire Plane 1 [deg];Events",37,0,185);
  trkThetaToWireA2CC2 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWireA2CC2",";Theta To Wire Plane 2 [deg];Events",37,0,185);
  trkNumHitsA2CC0 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2CC0",";Number Plane 0 Hits associated to Track;Events",150,0,3000);
  trkNumHitsA2CC1 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2CC1",";Number Plane 1 Hits associated to Track;Events",150,0,3000);
  trkNumHitsA2CC2 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2CC2",";Number Plane 2 Hits associated to Track;Events",150,0,3000);
  trkNumWiresA2CC0 = tfs->make<TH1D>("SimpleAnalysis_trkNumWiresA2CC0",";Number Expected Plane 0 Wires for Track;Events",150,0,3000);
  trkNumWiresA2CC1 = tfs->make<TH1D>("SimpleAnalysis_trkNumWiresA2CC1",";Number Expected Plane 1 Wires for Track;Events",150,0,3000);
  trkNumWiresA2CC2 = tfs->make<TH1D>("SimpleAnalysis_trkNumWiresA2CC2",";Number Expected Plane 2 Wires for Track;Events",150,0,3000);
  trkNumHitsPerWireA2CC0 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWireA2CC0",";Number of Hits / Number of Wires [Plane0];Events",300,0,15);
  trkNumHitsPerWireA2CC1 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWireA2CC1",";Number of Hits / Number of Wires [Plane1];Events",300,0,15);
  trkNumHitsPerWireA2CC2 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWireA2CC2",";Number of Hits / Number of Wires [Plane2];Events",300,0,15);
  trkNumWiresNumHitsA2CC0 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHitsA2CC0",";Number of Wires [Plane 0];Number of Hits [Plane0]",150,0,3000,150,0,3000);
  trkNumWiresNumHitsA2CC1 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHitsA2CC1",";Number of Wires [Plane 1];Number of Hits [Plane1]",150,0,3000,150,0,3000);
  trkNumWiresNumHitsA2CC2 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHitsA2CC2",";Number of Wires [Plane 2];Number of Hits [Plane2]",150,0,3000,150,0,3000);

  // Plots for Anode To Anode Candidates (A2A Candidates)
  trkLengthAllA2AC = tfs->make<TH1D>("SimpleAnalysis_trkLengthAllA2AC",";Considered track lengths [cm];Events",503,-6,1000);
  trkDistXA2AC = tfs->make<TH1D>("SimpleAnalysis_trkDistXA2AC",";Track Dist X [cm];Events",252,-8,1000);
  trkDistYA2AC = tfs->make<TH1D>("SimpleAnalysis_trkDistYA2AC",";Track Dist Y [cm];Events",252,-8,1000);
  trkDistZA2AC = tfs->make<TH1D>("SimpleAnalysis_trkDistZA2AC",";Track Dist Z [cm];Events",502,-8,2000);
  trkNumHitsA2AC = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2AC",";Number Hits associated to Track;Events",600,0,6000);
  trkThetaA2AC = tfs->make<TH1D>("SimpleAnalysis_trkThetaA2AC",";Theta [rad];Events",100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkPhiA2AC = tfs->make<TH1D>("SimpleAnalysis_trkPhiA2AC",";Phi [rad];Events",100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkThetaPhiA2AC = tfs->make<TH2D>("SimpleAnalysis_trkThetaPhiA2AC",";Theta [rad];Phi [rad]",100,-2.*TMath::Pi(),2.*TMath::Pi(),100,-2.*TMath::Pi(),2.*TMath::Pi());
  trkThetaXZA2AC = tfs->make<TH1D>("SimpleAnalysis_trkThetaXZA2AC",";ThetaXZ [deg];Events",100,-180,180);
  trkThetaYZA2AC = tfs->make<TH1D>("SimpleAnalysis_trkThetaYZA2AC",";ThetaYZ [deg];Events",100,-180,180);
  trkThetaXZThetaYZA2AC = tfs->make<TH2D>("SimpleAnalysis_trkThetaXZThetaYZA2AC",";ThetaXZ [deg];ThetaYZ [deg]",100,-180,180,100,-180,180);
  trkStartXA2AC = tfs->make<TH1D>("SimpleAnalysis_trkStartXA2AC",";Track Start X [cm];Events",100,-500,500);
  trkStartYA2AC = tfs->make<TH1D>("SimpleAnalysis_trkStartYA2AC",";Track Start Y [cm];Events",40,-200,200);
  trkStartZA2AC = tfs->make<TH1D>("SimpleAnalysis_trkStartZA2AC",";Track Start Z [cm];Events",100,-1000,1000);
  trkThetaToWireA2AC0 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWireA2AC0",";Theta To Wire Plane 0 [deg];Events",37,0,185);
  trkThetaToWireA2AC1 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWireA2AC1",";Theta To Wire Plane 1 [deg];Events",37,0,185);
  trkThetaToWireA2AC2 = tfs->make<TH1D>("SimpleAnalysis_trkThetaToWireA2AC2",";Theta To Wire Plane 2 [deg];Events",37,0,185);
  trkNumHitsA2AC0 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2AC0",";Number Plane 0 Hits associated to Track;Events",150,0,3000);
  trkNumHitsA2AC1 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2AC1",";Number Plane 1 Hits associated to Track;Events",150,0,3000);
  trkNumHitsA2AC2 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsA2AC2",";Number Plane 2 Hits associated to Track;Events",150,0,3000);
  trkNumWiresA2AC0 = tfs->make<TH1D>("SimpleAnalysis_trkNumWiresA2AC0",";Number Expected Plane 0 Wires for Track;Events",150,0,3000);
  trkNumWiresA2AC1 = tfs->make<TH1D>("SimpleAnalysis_trkNumWiresA2AC1",";Number Expected Plane 1 Wires for Track;Events",150,0,3000);
  trkNumWiresA2AC2 = tfs->make<TH1D>("SimpleAnalysis_trkNumWiresA2AC2",";Number Expected Plane 2 Wires for Track;Events",150,0,3000);
  trkNumHitsPerWireA2AC0 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWireA2AC0",";Number of Hits / Number of Wires [Plane0];Events",300,0,15);
  trkNumHitsPerWireA2AC1 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWireA2AC1",";Number of Hits / Number of Wires [Plane1];Events",300,0,15);
  trkNumHitsPerWireA2AC2 = tfs->make<TH1D>("SimpleAnalysis_trkNumHitsPerWireA2AC2",";Number of Hits / Number of Wires [Plane2];Events",300,0,15);
  trkNumWiresNumHitsA2AC0 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHitsA2AC0",";Number of Wires [Plane 0];Number of Hits [Plane0]",150,0,3000,150,0,3000);
  trkNumWiresNumHitsA2AC1 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHitsA2AC1",";Number of Wires [Plane 1];Number of Hits [Plane1]",150,0,3000,150,0,3000);
  trkNumWiresNumHitsA2AC2 = tfs->make<TH2D>("SimpleAnalysis_trkNumWiresNumHitsA2AC2",";Number of Wires [Plane 2];Number of Hits [Plane2]",150,0,3000,150,0,3000);
}

void SimpleTrackAndPFParticleAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // TODO: should I actually be doing exceptions instead of returns in some of these places...

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

  //  std::cout << "pfparticles.size() = " << pfparticles.size() << std::endl;
  mf::LogInfo("SimpleTrackAndPFParticleAna") << "pfparticles.size() = " << pfparticles.size() << std::endl;

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
      mf::LogWarning("SimpleTrackAndPFParticleAna") << "Event had no vertices, but no PFParticles. Events with no PFParticles were meant to be filtered earlier... Returning.";
    else mf::LogError("SimpleTrackAndPFParticleAna") << "Event had no vertices, but had PFParticles. Returning.";
    return;
  }

  ///////////////////////////////
  // Get the track association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::Track> fmtrk(pfparticleHandle, e, fTrackModuleLabel);
  if( !fmtrk.isValid() ){
    // Check against validity of fmtrk (using isValid, as used in an analyzer module from SBND)
    mf::LogError("SimpleTrackAndPFParticleAna") << "Error in validity of fmtrk. Returning.";
    return;
  }

  ///////////////////////////////
  // Get the spacepoint association to the PFPs
  ///////////////////////////////
  art::FindManyP<recob::SpacePoint> fmsps(pfparticleHandle, e, fPFParticleModuleLabel);
  if( !fmsps.isValid() ){
    // Check against validity of fmsps (using isValid)
    mf::LogError("SimpleTrackAndPFParticle") << "Error in validity of fmsps. Returning";
    return;
  }

  ///////////////////////////////
  // Get the hit assocation to the tracks
  ///////////////////////////////
  // get the tracks
  art::Handle< std::vector<recob::Track> > trackHandle;
  std::vector< art::Ptr<recob::Track> > tracks;
  if( e.getByLabel(fTrackModuleLabel,trackHandle) ){
    art::fill_ptr_vector(tracks,trackHandle);
  }
  else{
    mf::LogWarning("SimpleTrackAndPFParticleAna") << "Event failed to find recob::Track. Are you sure you gave the right label and are running this at the right time?";
    return;
  }
  // now findmany
  art::FindManyP<recob::Hit> fmtrkhit(trackHandle,e,fTrackModuleLabel);
  if( !fmtrkhit.isValid() ){
    // Check against validity of fmtrkhit (using isValid)
    mf::LogError("SimpleTrackAndPFParticleAna") << "Error in validity of fmtrkhit. Returning.";
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
  //float maxTrkLengthX = -5.;
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

    vtxDetXAll->Fill( vtxX );
    vtxDetYAll->Fill( vtxY );
    vtxDetZAll->Fill( vtxZ );

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
    // Necessary to get some of this info now since will use it in plots before and after cuts...
    double trkLength = pfpTrk.at(0)->Length();
    double trkX0 = pfpTrk.at(0)->Start().X();
    double trkXf = pfpTrk.at(0)->End().X();
    double trkY0 = pfpTrk.at(0)->Start().Y();
    double trkYf = pfpTrk.at(0)->End().Y();
    double trkZ0 = pfpTrk.at(0)->Start().Z();
    double trkZf = pfpTrk.at(0)->End().Z();
    double trkAngleTheta = pfpTrk.at(0)->Theta();
    double trkAnglePhi = pfpTrk.at(0)->Phi();
    double trkAngleThetaXZ = (180./TMath::Pi())*TMath::ATan2(pfpTrk.at(0)->StartDirection().X(),pfpTrk.at(0)->StartDirection().Z());
    double trkAngleThetaYZ = (180./TMath::Pi())*TMath::ATan2(pfpTrk.at(0)->StartDirection().Y(),pfpTrk.at(0)->StartDirection().Z());

    trkStartXAll->Fill(trkX0);
    trkStartYAll->Fill(trkY0);
    trkStartZAll->Fill(trkZ0);

    // continue if track length < min length cut
    if( trkLength < fTrkMinLenCut ) continue;

    trkStartXPassTrkLen->Fill(trkX0);
    trkStartYPassTrkLen->Fill(trkY0);
    trkStartZPassTrkLen->Fill(trkZ0);

    // get the TPC from which to iterate planes from PositionToTPCID
    // -- if start is null try end
    // -- if both are null continue
    //
    // TODO: Is this okay to do?
    //
    // We wont use this till later when looping through the planes,
    // but if we're going to continue, let's figure out now
    auto trkTPCID = geom->PositionToTPCID( pfpTrk.at(0)->Start() ); // tpc containing the start point
    //std::cout << pfpTrk.at(0)->Start().X() << " " << pfpTrk.at(0)->Start().Y() << " " << pfpTrk.at(0)->Start().Z() << std::endl;
    if( !trkTPCID.isValid ){
      trkTPCID = geom->PositionToTPCID( pfpTrk.at(0)->End() ); // fallback to tpc containing the end point
      //std::cout << pfpTrk.at(0)->End().X() << " " << pfpTrk.at(0)->End().Y() << " " << pfpTrk.at(0)->End().Z() << std::endl;
    }
    if( !trkTPCID.isValid ) {
      trkThetaXZFailAV->Fill( trkAngleThetaXZ );
      trkThetaYZFailAV->Fill( trkAngleThetaYZ );

      mf::LogInfo("SimpleTrackAndPFParticleAna") << "Neither the track start nor end returned a valid TPCID. Skipping this PFParticle";
      continue;
    }

    // Alright, if we're not going to skip the PFParticle, then let's start calculating things and filling things
    vtxDetX->Fill( vtxX );
    vtxDetY->Fill( vtxY );
    vtxDetZ->Fill( vtxZ );
    vtxDetZX->Fill( vtxZ, vtxX );
    vtxDetXY->Fill( vtxX, vtxY );

    // for now just fill the number of spacepoints
    std::vector< art::Ptr<recob::SpacePoint> > pfpSps = fmsps.at(iPfp.key());
    pfpNumSps->Fill( pfpSps.size() );

    //    std::cout << "length = " << trkLength
    //	      << " ... Dir X = " << pfpTrk.at(0)->StartDirection().Unit().X()
    //	      << " ... Dir Y = " << pfpTrk.at(0)->StartDirection().Unit().Y() << std::endl;

    trkLengthAll->Fill( trkLength );
    //trkLengthAllX->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X()) );
    //trkLengthAllY->Fill( trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().Y()) );

    trkStartX->Fill( trkX0 );
    trkStartY->Fill( trkY0 );
    trkStartZ->Fill( trkZ0 );
    trkEndX->Fill( trkXf );
    trkEndY->Fill( trkYf );
    trkEndZ->Fill( trkZf );
    trkStartEndX->Fill( trkX0, trkXf );
    trkStartEndY->Fill( trkY0, trkYf );
    trkStartEndZ->Fill( trkZ0, trkZf );
    trkStartXLength->Fill( trkX0, trkLength );
    trkStartYLength->Fill( trkY0, trkLength );
    trkStartZLength->Fill( trkZ0, trkLength );
    trkDistX->Fill( fabs(trkXf-trkX0) );
    trkDistY->Fill( fabs(trkYf-trkY0) );
    trkDistZ->Fill( fabs(trkZf-trkZ0) );

    trkTheta->Fill( trkAngleTheta );
    trkPhi->Fill( trkAnglePhi );
    trkThetaPhi->Fill( trkAngleTheta, trkAnglePhi );
    trkThetaXZ->Fill( trkAngleThetaXZ );
    trkThetaYZ->Fill( trkAngleThetaYZ );
    trkThetaXZThetaYZ->Fill( trkAngleThetaXZ, trkAngleThetaYZ );

    if( trkLength > maxTrkLength ){
      maxTrkLength = trkLength;
      //maxTrkLengthX = trkLength * fabs(pfpTrk.at(0)->StartDirection().Unit().X());
    }

    // HITS
    // -------------------
    // access the associated hits
    std::vector< art::Ptr<recob::Hit> > trkHits = fmtrkhit.at( pfpTrk.at(0).key() );
    trkNumHits->Fill( trkHits.size() );
    trkLengthNumHits->Fill( trkLength, trkHits.size() );
    // pieces for calculation of hit efficiency
    unsigned long trkHitsByPlane[3]        = {0,0,0};                 // count per plane
    int           trkNumWiresByPlane[3]    = {0,0,0};                 // expected number of wires to cross
    double        trkThetaToWireByPlane[3] = {0.,0.,0.};              // angle of track to the wire direction
    double        trkHitsPerWireByPlane[3] = {0.,0.,0.};              // ratio of number of hits to number of wires
    TVector3      trkLenVec( trkXf-trkX0, trkYf-trkY0, trkZf-trkZ0 ); // vector for the track
    TVector3      trkLenVecYZ( 0., trkYf-trkY0, trkZf-trkZ0 );        // just the portion in the YZ plane to compare to wire direction
    // loop through hits, for now just sum up hits on plane
    for( auto const& iHit : trkHits ){
      trkHitsByPlane[ iHit->WireID().Plane ]+=1;
    }
    // following the example as in Doxygen, iterate over PlaneIDs - use the trkTPCID (found above)
    for( geo::PlaneID const& iPlaneID : geom->IteratePlaneIDs(trkTPCID) ){
      geo::PlaneGeo const& iPlane = geom->Plane(iPlaneID);

      // can use PlaneID to get plane index (thanks Tracy!)
      unsigned int iPlaneNum = iPlaneID.Plane;

      // get the expected number of wires and angle of track w.r.t. wire (assumes 3 mm between wires!)
      trkNumWiresByPlane[iPlaneNum] = int( fabs( trkLenVec.Dot(iPlane.GetIncreasingWireDirection())/0.3 ) );
      trkThetaToWireByPlane[iPlaneNum] = (180./TMath::Pi())*trkLenVecYZ.Angle(iPlane.GetWireDirection());

      if( trkNumWiresByPlane[iPlaneNum]>0 )
	trkHitsPerWireByPlane[iPlaneNum] = double(trkHitsByPlane[iPlaneNum])/double(trkNumWiresByPlane[iPlaneNum]);
    }

    // FOR THE CASE of "All" tracks, i.e. not anode-to-cathode candidates, put a minimum 
    // length requirement of 50cm and play around with the angle restrictions a bit...
    if( trkLength>50. ){
      // fill the track wires/hit info
      trkThetaToWire0->Fill( trkThetaToWireByPlane[0] );
      trkThetaToWire1->Fill( trkThetaToWireByPlane[1] );
      trkThetaToWire2->Fill( trkThetaToWireByPlane[2] );
      // fill if the track meets angular criteria
      if( (trkAngleThetaXZ>-30. && trkAngleThetaXZ<30.) || trkAngleThetaXZ<-150. || trkAngleThetaXZ>150. ){
	if( trkThetaToWireByPlane[0]>20. && trkThetaToWireByPlane[0]<160. && trkNumWiresByPlane[0]>0 ){
	  trkNumHits0->Fill( trkHitsByPlane[0] );
	  trkNumWires0->Fill( trkNumWiresByPlane[0] );
	  trkNumWiresNumHits0->Fill( trkNumWiresByPlane[0],trkHitsByPlane[0] );
	  trkNumHitsPerWire0->Fill( trkHitsPerWireByPlane[0] );
	}
	if( trkThetaToWireByPlane[1]>20. && trkThetaToWireByPlane[1]<160. && trkNumWiresByPlane[1]>0 ){
	  trkNumHits1->Fill( trkHitsByPlane[1] );
	  trkNumWires1->Fill( trkNumWiresByPlane[1] );
	  trkNumWiresNumHits1->Fill( trkNumWiresByPlane[1],trkHitsByPlane[1] );
	  trkNumHitsPerWire1->Fill( trkHitsPerWireByPlane[1] );
	}
	if( trkThetaToWireByPlane[2]>20. && trkThetaToWireByPlane[2]<160. && trkNumWiresByPlane[2]>0 ){
	  trkNumHits2->Fill( trkHitsByPlane[2] );
	  trkNumWires2->Fill( trkNumWiresByPlane[2] );
	  trkNumWiresNumHits2->Fill( trkNumWiresByPlane[2],trkHitsByPlane[2] );
	  trkNumHitsPerWire2->Fill( trkHitsPerWireByPlane[2] );
	}
      }
    }

    // now look at if the track is a possible anode-cathode or anode-anode crosser and fill up some variables
    if( fPrintAnodeToCathode && fabs(trkXf-trkX0)>=293. ){
      mf::LogInfo("SimpleTrackAndPFParticleAna") << "Potential anode-to-anode crossing track with length (" << fabs(trkXf-trkX0)
                                                 << ") and start point (" << trkX0 << ", " << trkY0 << ", " << trkZ0 << ")";
    }
    else if( fPrintAnodeToCathode && fabs(trkXf-trkX0)>=146.5 ){
      mf::LogInfo("SimpleTrackAndPFParticleAna") << "Potential anode-to-cathode crossing track with length (" << fabs(trkXf-trkX0)
                                                 << ") and start point (" << trkX0 << ", " << trkY0 << ", " << trkZ0 << ")";
    }

    // anode to cathode candidates
    if( fabs(trkXf-trkX0) > 146.5 ){
      //std::cout << "FILLING A2CC" << std::endl;
      trkLengthAllA2CC->Fill( trkLength );
      trkDistXA2CC->Fill( fabs(trkXf-trkX0) );
      trkDistYA2CC->Fill( fabs(trkYf-trkY0) );
      trkDistZA2CC->Fill( fabs(trkZf-trkZ0) );
      trkStartXA2CC->Fill( trkX0 );
      trkStartYA2CC->Fill( trkY0 );
      trkStartZA2CC->Fill( trkZ0 );
      trkThetaA2CC->Fill( trkAngleTheta );
      trkPhiA2CC->Fill( trkAnglePhi );
      trkThetaPhiA2CC->Fill( trkAngleTheta, trkAnglePhi );
      trkThetaXZA2CC->Fill( trkAngleThetaXZ );
      trkThetaYZA2CC->Fill( trkAngleThetaYZ );
      trkThetaXZThetaYZA2CC->Fill( trkAngleThetaXZ, trkAngleThetaYZ );
      trkNumHitsA2CC->Fill( trkHits.size() );

      // fill the track wires/hit info
      trkThetaToWireA2CC0->Fill( trkThetaToWireByPlane[0] );
      trkThetaToWireA2CC1->Fill( trkThetaToWireByPlane[1] );
      trkThetaToWireA2CC2->Fill( trkThetaToWireByPlane[2] );
      // fill if the track meets angular criteria
      if( (trkAngleThetaXZ>-20. && trkAngleThetaXZ<20.) || trkAngleThetaXZ<-160. || trkAngleThetaXZ>160. ){
	if( trkThetaToWireByPlane[0]>60. && trkThetaToWireByPlane[0]<120. && trkNumWiresByPlane[0]>0 ){
          trkNumHitsA2CC0->Fill( trkHitsByPlane[0] );
          trkNumWiresA2CC0->Fill( trkNumWiresByPlane[0] );
	  trkNumWiresNumHitsA2CC0->Fill( trkNumWiresByPlane[0],trkHitsByPlane[0] );
          trkNumHitsPerWireA2CC0->Fill( trkHitsPerWireByPlane[0] );
	}
	if( trkThetaToWireByPlane[1]>60. && trkThetaToWireByPlane[1]<120. && trkNumWiresByPlane[1]>0 ){
	  trkNumHitsA2CC1->Fill( trkHitsByPlane[1] );
	  trkNumWiresA2CC1->Fill( trkNumWiresByPlane[1] );
	  trkNumWiresNumHitsA2CC1->Fill( trkNumWiresByPlane[1],trkHitsByPlane[1] );
	  trkNumHitsPerWireA2CC1->Fill( trkHitsPerWireByPlane[1] );
	}
	if( trkThetaToWireByPlane[2]>60. && trkThetaToWireByPlane[2]<120. && trkNumWiresByPlane[2]>0 ){
	  trkNumHitsA2CC2->Fill( trkHitsByPlane[2] );
	  trkNumWiresA2CC2->Fill( trkNumWiresByPlane[2] );
	  trkNumWiresNumHitsA2CC2->Fill( trkNumWiresByPlane[2],trkHitsByPlane[2] );
	  trkNumHitsPerWireA2CC2->Fill( trkHitsPerWireByPlane[2] );
	}
      }
    }

    // anode to anode candidates
    if( fabs(trkXf-trkX0) > 293. ){
      //std::cout << "FILLING A2AC" << std::endl;
      trkLengthAllA2AC->Fill( trkLength );
      trkDistXA2AC->Fill( fabs(trkXf-trkX0) );
      trkDistYA2AC->Fill( fabs(trkYf-trkY0) );
      trkDistZA2AC->Fill( fabs(trkZf-trkZ0) );
      trkStartXA2AC->Fill( trkX0 );
      trkStartYA2AC->Fill( trkY0 );
      trkStartZA2AC->Fill( trkZ0 );
      trkThetaA2AC->Fill( trkAngleTheta );
      trkPhiA2AC->Fill( trkAnglePhi );
      trkThetaPhiA2AC->Fill( trkAngleTheta, trkAnglePhi );
      trkThetaXZA2AC->Fill( trkAngleThetaXZ );
      trkThetaYZA2AC->Fill( trkAngleThetaYZ );
      trkThetaXZThetaYZA2AC->Fill( trkAngleThetaXZ, trkAngleThetaYZ );
      trkNumHitsA2AC->Fill( trkHits.size() );

      // fill the track wires/hit info
      trkThetaToWireA2AC0->Fill( trkThetaToWireByPlane[0] );
      trkThetaToWireA2AC1->Fill( trkThetaToWireByPlane[1] );
      trkThetaToWireA2AC2->Fill( trkThetaToWireByPlane[2] );
      // fill if the track meets angular criteria
      if( (trkAngleThetaXZ>-20. && trkAngleThetaXZ<20.) || trkAngleThetaXZ<-160. || trkAngleThetaXZ>160. ){
        if( trkThetaToWireByPlane[0]>60. && trkThetaToWireByPlane[0]<120. && trkNumWiresByPlane[0]>0 ){
          trkNumHitsA2AC0->Fill( trkHitsByPlane[0] );
          trkNumWiresA2AC0->Fill( trkNumWiresByPlane[0] );
	  trkNumWiresNumHitsA2AC0->Fill( trkNumWiresByPlane[0],trkHitsByPlane[0] );
          trkNumHitsPerWireA2AC0->Fill( trkHitsPerWireByPlane[0] );
        }
        if( trkThetaToWireByPlane[1]>60. && trkThetaToWireByPlane[1]<120. && trkNumWiresByPlane[1]>0 ){
          trkNumHitsA2AC1->Fill( trkHitsByPlane[1] );
          trkNumWiresA2AC1->Fill( trkNumWiresByPlane[1] );
	  trkNumWiresNumHitsA2AC1->Fill( trkNumWiresByPlane[1],trkHitsByPlane[1] );
          trkNumHitsPerWireA2AC1->Fill( trkHitsPerWireByPlane[1] );
        }
        if( trkThetaToWireByPlane[2]>60. && trkThetaToWireByPlane[2]<120. && trkNumWiresByPlane[2]>0 ){
          trkNumHitsA2AC2->Fill( trkHitsByPlane[2] );
          trkNumWiresA2AC2->Fill( trkNumWiresByPlane[2] );
	  trkNumWiresNumHitsA2AC2->Fill( trkNumWiresByPlane[2],trkHitsByPlane[2] );
          trkNumHitsPerWireA2AC2->Fill( trkHitsPerWireByPlane[2] );
        }
      }
    }
  } // End looop through PFParticle

  trkMaxLength->Fill( maxTrkLength );
  //trkMaxLengthX->Fill( maxTrkLengthX );

  return;
}

DEFINE_ART_MODULE(SimpleTrackAndPFParticleAna)
