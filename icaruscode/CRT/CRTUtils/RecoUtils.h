#ifndef FDSELECTIONUTILS_H_SEEN
#define FDSELECTIONUTILS_H_SEEN


///////////////////////////////////////////////
// RecoUtils.h
//
// A few reco utilities like truth matching 
// D Brailsford (adapted from work by D Brailsford and M Wallbank), October 2017
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
//#include "lardataobj/RecoBase/Track.h"
//#include "lardataobj/RecoBase/Shower.h"
//#include "lardataobj/AnalysisBase/MVAPIDResult.h"
//#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larcore/Geometry/Geometry.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
namespace detinfo { class DetectorClocksData; }


// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"


namespace icarus::crt::dataTools{

using ModuleCenter = geo::Point_t;

struct AffineTrans
{
  double A11;
  double A12;
  double A21;
  double A22;
  double B1;
  double B2;
  double Accuracy;
  double N;
};

using FebIndex_t = int;

using TopCRTCorrectionMap = std::map<FebIndex_t, AffineTrans>;

struct TopCRTTransformations
{
  TopCRTCorrectionMap EE;
  TopCRTCorrectionMap EW;
  TopCRTCorrectionMap EastCC;
  TopCRTCorrectionMap WE;
  TopCRTCorrectionMap WW;
  TopCRTCorrectionMap WestCC;
  bool imported;
};

using TopCRTCentersMap = std::map<FebIndex_t, ModuleCenter>;

/// The transformed CRT Hits are in cm
using TransformedCRTHit = std::pair<double, double>;

/// This function loads the Top CRT modules centers.
TopCRTCentersMap LoadTopCRTCenters();

/// This function performs the affine transformation of the CRT hit points.
/// The AffineTransformation requires input variables in cm.
TransformedCRTHit AffineTransformation(double DX, double DZ, AffineTrans affine);

/// This functions loads the Affine Transformation TXT files.
TopCRTTransformations LoadTopCRTTransformations();

/// This function transforms a CRTHit with AffineTransformationFunctions.
geo::Point_t ApplyTransformation(const sbn::crt::CRTHit& crthit, const TopCRTCorrectionMap& TopCRTCorrections, const TopCRTCentersMap& TopCRTCenters);

}

namespace RecoUtils{

// Returns the geant4 ID which contributes the most to a single reco hit.  
// The matching method looks for true particle which deposits the most true energy in the reco hit.  
// If rollup_unsaved_ids is set to true, any unsaved daughter than 
// contributed energy to the hit has its energy included in its closest ancestor that was saved.

int TrueParticleID(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Hit> hit, bool rollup_unsaved_ids=1); 
  
// Returns the geant4 ID which contributes the most to the vector of hits.  
// The matching method looks for which true particle deposits the most true energy in the reco hits

int TrueParticleIDFromTotalTrueEnergy(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids=1); 

// Returns the geant4 ID which contributes the most to the vector of hits.  
// The matching method looks for which true particle contributes the most reconstructed charge to the hit selection 
// (the reco charge of each hit is correlated with each maximally contributing true particle and summed)
  
int TrueParticleIDFromTotalRecoCharge(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids=1); 

// Returns the geant4 ID which contributes the most to the vector of hits.  
// The matching method looks for which true particle maximally contributes to the most reco hits
  
int TrueParticleIDFromTotalRecoHits(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids=1); 

// Checks if a position is within any of the TPCs in the geometry (user can define some distance buffer from the TPC walls)
  
bool IsInsideTPC(TVector3 position, double distance_buffer); 

// Calculates the total length of a recob::track by summing up the distances between adjacent traj. points

double CalculateTrackLength(const art::Ptr<recob::Track> track); 
  
std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels, const geo::GeometryCore &geo);

std::map<int, std::vector<art::Ptr<recob::Hit>>> buildTrackIDtoHitsMap(const std::vector<art::Ptr<recob::Hit>> &allHits, const detinfo::DetectorClocksData &clockData, const cheat::BackTrackerService &backtracker);
}

#endif
