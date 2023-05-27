////////////////////////////////////////////////////////////////////////
// Class:       HARPS (Hit Activity Removal from Particles for Systematics)
// Plugin Type: producer (Unknown Unknown)
// File:        HARPS_module.cc
//
// Generated at Tue May 23 10:46:39 2023 by Bruce Howard using cetskelgen
// from  version .
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

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// For random implementation: see
// https://github.com/SBNSoftware/icaruscode/blob/7b71f0213a1a66913cfe1cd795f49a9215570145/icaruscode/PMT/OpReco/FakePhotoS_module.cc#L16
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>

class HARPS;


class HARPS : public art::EDProducer {
public:
  explicit HARPS(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HARPS(HARPS const&) = delete;
  HARPS(HARPS&&) = delete;
  HARPS& operator=(HARPS const&) = delete;
  HARPS& operator=(HARPS&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  bool fOnlyOnePerEvent;

  std::vector<art::InputTag> fHitModuleLabels;
  std::vector<art::InputTag> fPFParticleModuleLabels;
  std::vector<art::InputTag> fTrackModuleLabels;

  std::string fInterestFileName;

  // TODO: how to handle the hits without spacepoints...
  double fMaskBeginningDist; ///< How far from beginning of track to mask (in cm)
  double fMaskEndDist;       ///< How far from end of track to mask (in cm)

  // Allow for one to choose random groups of wires to skip on the different planes...
  // TODO: Attempt made to handle things crossing boundaries, but possibly could be improved.
  std::vector< unsigned int > fMaskNBreak;
  std::vector< int > fMaskNWires;

  bool fTPCHitsWireAssn;
  std::string fTPCHitCreatorInstanceName;

  std::map< std::string, std::vector<size_t> > fParticleListCryo0;
  std::map< std::string, std::vector<size_t> > fParticleListCryo1;

  CLHEP::HepRandomEngine &fFlatEngine;
  CLHEP::RandFlat *fFlatRand; ///< Random number generator as in PMT/OpReco/FakePhotoS_module.cc
};


HARPS::HARPS(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fOnlyOnePerEvent           ( p.get< bool >("OnlyOnePerEvent", true) ),
  fHitModuleLabels           ( p.get< std::vector<art::InputTag> >("HitModuleLabels") ),
  fPFParticleModuleLabels    ( p.get< std::vector<art::InputTag> >("PFParticleModuleLabels") ),
  fTrackModuleLabels         ( p.get< std::vector<art::InputTag> >("TrackModuleLabels") ),
  fInterestFileName          ( p.get<std::string>("InterestFileName") ),
  fMaskBeginningDist         ( p.get<double>("MaskBeginningDist", -1.) ),
  fMaskEndDist               ( p.get<double>("MaskEndDist", -1.) ),
  fMaskNBreak                ( p.get< std::vector<unsigned int> >("MaskNBreak", {0}) ),
  fMaskNWires                ( p.get< std::vector<int> >("MaskNWires", {-1}) ),
  fTPCHitsWireAssn           ( p.get< bool >("TPCHitsWireAssn", true) ),
  fTPCHitCreatorInstanceName ( p.get<std::string>("TPCHitCreatorInstanaceName","") ),
  fFlatEngine                ( art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "Gen", p, "Seed") )
{
  if ( fPFParticleModuleLabels.size()!=fTrackModuleLabels.size() || fPFParticleModuleLabels.size()!=fHitModuleLabels.size() )
    throw cet::exception("HARPS") << "Error... InputTag vectors need to be same size..." << std::endl; // borrow from elsewhere

  if ( fMaskBeginningDist>0. && fMaskEndDist>0. )
    throw cet::exception("HARPS") << "Error... Should not try to mask the beginning and ALSO some distance from the end..." << std::endl;

  recob::HitCollectionCreator::declare_products(producesCollector(), fTPCHitCreatorInstanceName, fTPCHitsWireAssn, false);

  // Load in the PFP list of interest and put contents into a map for fast lookup. BASED ON larevt/Filters/EventFilter_module.cc
  std::ifstream in;
  in.open(fInterestFileName.c_str());
  char line[1024];
  while (1) {
    in.getline(line, 1024);
    if (!in.good()) break;
    unsigned int nR, nS, nE, nC, nP;
    sscanf(line, "%u %u %u %u %u", &nR, &nS, &nE, &nC, &nP);
    std::string evtID = std::to_string(nR) + ":" + std::to_string(nS) + ":" + std::to_string(nE);
    if ( fOnlyOnePerEvent && (fParticleListCryo0.find(evtID) != fParticleListCryo0.end() ||
                              fParticleListCryo1.find(evtID) != fParticleListCryo1.end()) ) continue;
    if ( nC == 0 ) fParticleListCryo0[evtID].push_back(nP);
    else           fParticleListCryo1[evtID].push_back(nP);
  }
  in.close();

  fFlatRand = new CLHEP::RandFlat(fFlatEngine,0.,1.);
}

void HARPS::produce(art::Event& e)
{
  // Return empty hit container if no particles of interest.
  std::string evtID = std::to_string(e.run()) + ":" + std::to_string(e.subRun()) + ":" + std::to_string(e.event());
  if ( fParticleListCryo0.find(evtID) == fParticleListCryo0.end() &&
       fParticleListCryo1.find(evtID) == fParticleListCryo1.end() )
    return;

  // As in Gauss Hit Finder code (larreco/HitFinder/GausHitFinder_module.cc)
  recob::HitCollectionCreator hitCol(e, fTPCHitCreatorInstanceName, fTPCHitsWireAssn, false);

  // Load in the PFParticles, Tracks, and Hits & the necessary FindManyP's
  for ( unsigned int iCryo=0; iCryo<2; ++iCryo ) {
    if ( iCryo==0 && fParticleListCryo0.find(evtID) == fParticleListCryo0.end() ) continue;
    if ( iCryo==1 && fParticleListCryo1.find(evtID) == fParticleListCryo1.end() ) continue;

    // as in previous code
    art::Handle< std::vector<recob::PFParticle> > pfpHandle;
    std::vector< art::Ptr<recob::PFParticle> > pfps;
    if ( e.getByLabel(fPFParticleModuleLabels[iCryo], pfpHandle) ) {
      art::fill_ptr_vector(pfps, pfpHandle);
    }
    else {
      mf::LogError("HARPS") << "Error pulling in PFParticle product... Skipping cryostat.";
    }

    art::Handle< std::vector<recob::Track> > trkHandle;
    std::vector< art::Ptr<recob::Track> > trks;
    if ( e.getByLabel(fTrackModuleLabels[iCryo], trkHandle) ) {
      art::fill_ptr_vector(trks, trkHandle);
    }
    else {
      mf::LogError("HARPS") << "Error pulling in Track product... Skipping cryostat.";
    }

    art::Handle< std::vector<recob::Hit> > hitHandle;
    std::vector< art::Ptr<recob::Hit> > hits;
    if ( e.getByLabel(fHitModuleLabels[iCryo], hitHandle) ) {
      art::fill_ptr_vector(hits, hitHandle);
    }
    else {
      mf::LogError("HARPS") << "Error pulling in Hit product... Skipping cryostat.";
    }

    // For the Wire Assns in HitCreator
    art::FindManyP<recob::Wire> fmWire(hitHandle, e, fHitModuleLabels[iCryo]);
    if ( !fmWire.isValid() && fTPCHitsWireAssn ){
      throw cet::exception("HARPS") << "Module wants to use Hit-Wire Associations, but fmWire invalid." << std::endl;
    }

    art::FindManyP<recob::Track>      fmTrk(pfpHandle, e, fTrackModuleLabels[iCryo]);
    art::FindManyP<recob::Hit>        fmHit(trkHandle, e, fTrackModuleLabels[iCryo]);
    art::FindManyP<recob::SpacePoint> fmSps(hitHandle, e, fPFParticleModuleLabels[iCryo]);
    if ( !fmTrk.isValid() || !fmHit.isValid() || !fmSps.isValid() ){
      throw cet::exception("HARPS") << "Module wants to use invalid associations." << std::endl;
    }

    // Find the PFParticles we care about and put their hits into the collection
    for ( auto const& pfp : pfps ) {
      bool keep = false;
      for ( auto const& id : (iCryo==0 ? fParticleListCryo0[evtID] : fParticleListCryo1[evtID]) ) {
        if ( pfp->Self() == id ) {
          keep = true;
          break;
        }
      }
      if( !keep ) continue;

      const std::vector< art::Ptr<recob::Track> > trkFromPFP = fmTrk.at(pfp.key());
      if ( trkFromPFP.size() != 1 ) {
        mf::LogError("HARPS") << "fmTrk gave us more than one track for this PFParticle... Skipping.";
        continue;
      }

      double trkStartX = trkFromPFP.at(0)->Start().X();
      double trkStartY = trkFromPFP.at(0)->Start().Y();
      double trkStartZ = trkFromPFP.at(0)->Start().Z();
      double trkEndX   = trkFromPFP.at(0)->End().X();
      double trkEndY   = trkFromPFP.at(0)->End().Y();
      double trkEndZ   = trkFromPFP.at(0)->End().Z();

      const std::vector< art::Ptr<recob::Hit> > hitsFromTrk = fmHit.at(trkFromPFP.at(0).key());

      // Make a vector of the hits for each plane, and while doing it, compile a list of idx = Nw * 4 * 3 * Cryo + Nw * 3 * TPC + Nw * Plane + Wire, for now hardcode Nw = 5000 which should be long enough...
      std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > planeToHitsMap;
      std::map< unsigned int, std::set< unsigned int > > planeToIdxMap;

      // Loop as in my overlay module, copy hits passing the criteria into the maps for checking against other criteria...
      for( auto const& iHitPtr : hitsFromTrk ) {
        // Check if we want to remove this hit because we are masking X distance from the beginning:
        // TODO: as a start, use the Pandora spacepoints<->hit association... is this safe or should we do something in 2d distances?
        if ( fMaskBeginningDist > 0. ) {
          const std::vector< art::Ptr<recob::SpacePoint> > spsFromHit = fmSps.at(iHitPtr.key());
          if ( spsFromHit.size() != 1 ) {
            mf::LogError("HARPS") << "fmSps gave us something other than one spacepoint for this hit (" << spsFromHit.size() << ")... Skipping.";
            continue;
          }
          const auto thisXYZ = spsFromHit.at(0)->XYZ();
          if ( std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] ) < fMaskBeginningDist ) continue;
        }
        else if ( fMaskEndDist > 0. ) {
          const std::vector< art::Ptr<recob::SpacePoint> > spsFromHit = fmSps.at(iHitPtr.key());
          if ( spsFromHit.size() != 1 ) {
            mf::LogError("HARPS") << "fmSps gave us something other than one spacepoint for this hit (" << spsFromHit.size() << ")... Skipping.";
            continue;
          }
          const auto thisXYZ = spsFromHit.at(0)->XYZ();
          if ( std::hypot( trkEndX-thisXYZ[0], trkEndY-thisXYZ[1], trkEndZ-thisXYZ[2] ) > fMaskEndDist ) continue;
        }

        unsigned int plane = iHitPtr->WireID().Plane;
        unsigned int idx = 5000*4*3*iHitPtr->WireID().Cryostat + 5000*3*iHitPtr->WireID().TPC + 5000*iHitPtr->WireID().Plane + iHitPtr->WireID().Wire;
        planeToHitsMap[plane].push_back( iHitPtr );
        planeToIdxMap[plane].emplace( idx );
      }

      // Pick wires to remove!
      std::map< unsigned int, std::set< unsigned int > > planeToIdxRemoval;
      if ( fMaskNWires.size() >= 1 && fMaskNBreak.size() >= 1 && fMaskNWires.size() == fMaskNBreak.size() ) {
        for ( unsigned int idxPlane = 0; idxPlane < fMaskNWires.size(); ++idxPlane ) {
          if( fMaskNWires[idxPlane] < 0 || fMaskNBreak[idxPlane] == 0 ) continue;
          if ( planeToIdxMap.find(idxPlane) == planeToIdxMap.end() ) continue;
          std::vector<unsigned int> planeIdxVec(planeToIdxMap[idxPlane].begin(),planeToIdxMap[idxPlane].end()); // vector copy for easy access... thanks Stack Overflow
          unsigned int ithBreak = 0;
          while ( ithBreak < fMaskNBreak[idxPlane] ) {
            unsigned int locToBreak = planeToIdxMap[idxPlane].size()*fFlatRand->fire(0.,1.);
            for ( unsigned int idxSubBreak=0; idxSubBreak < (unsigned int)fMaskNWires[idxPlane]; ++idxSubBreak ) {
              if ( locToBreak+idxSubBreak >= planeToIdxMap[idxPlane].size() ) continue;
              planeToIdxRemoval[idxPlane].emplace( planeIdxVec.at(locToBreak+idxSubBreak) );
            } // loop number of wires to break in this group
            ithBreak+=1;
          } // loop through number of groups of breaks to make in this plane
        } // loop through number of planes on which to add breaks
      }

      // Now loop over the hits we have and remove them if its wire didn't pass the test...
      for ( auto const& [ plane, hitVec ] : planeToHitsMap ) {
        for ( auto const& hit : hitVec ) {
          unsigned int plane = hit->WireID().Plane;
          unsigned int idx = 5000*4*3*hit->WireID().Cryostat + 5000*3*hit->WireID().TPC + 5000*hit->WireID().Plane + hit->WireID().Wire;
          if ( planeToIdxRemoval.find(plane) != planeToIdxRemoval.end() && planeToIdxRemoval[plane].count( idx ) ) continue;

          // If we're not skipping this hit then write it to the file!
          recob::Hit theHit = *hit;
          if ( fTPCHitsWireAssn ) {
            std::vector< art::Ptr<recob::Wire> > hitWires = fmWire.at(hit.key());
            if ( hitWires.size() == 0 ) throw cet::exception("HARPS") << "Hit found with no associated wires...\n";
            else if ( hitWires.size() > 1 ) mf::LogWarning("HARPS") << "Hit with >1 recob::Wire associated...";
            hitCol.emplace_back(theHit,hitWires[0]);
          }
          else hitCol.emplace_back(theHit);
        } // loop hits
      } // loop planes with hits
    } // loop pfps
  } // loop cryos

  hitCol.put_into(e);
}

DEFINE_ART_MODULE(HARPS)
