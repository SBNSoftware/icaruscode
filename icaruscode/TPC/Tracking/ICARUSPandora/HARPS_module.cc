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

  double fMaskBeginningDist; ///< How far from beginning of track to mask (in cm)

  bool fTPCHitsWireAssn;
  std::string fTPCHitCreatorInstanceName;

  std::map< std::string, std::vector<size_t> > fParticleListCryo0;
  std::map< std::string, std::vector<size_t> > fParticleListCryo1;
};


HARPS::HARPS(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fOnlyOnePerEvent           ( p.get< bool >("OnlyOnePerEvent", true) ),
  fHitModuleLabels           ( p.get< std::vector<art::InputTag> >("HitModuleLabels") ),
  fPFParticleModuleLabels    ( p.get< std::vector<art::InputTag> >("PFParticleModuleLabels") ),
  fTrackModuleLabels         ( p.get< std::vector<art::InputTag> >("TrackModuleLabels") ),
  fInterestFileName          ( p.get<std::string>("InterestFileName") ),
  fMaskBeginningDist         ( p.get<double>("MaskBeginningDist", -1.) ),
  fTPCHitsWireAssn           ( p.get< bool >("TPCHitsWireAssn", true) ),
  fTPCHitCreatorInstanceName ( p.get<std::string>("TPCHitCreatorInstanaceName","") )
{
  if ( fPFParticleModuleLabels.size()!=fTrackModuleLabels.size() || fPFParticleModuleLabels.size()!=fHitModuleLabels.size() )
    throw cet::exception("HARPS") << "Error... InputTag vectors need to be same size..." << std::endl; // borrow from elsewhere

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

      const std::vector< art::Ptr<recob::Hit> > hitsFromTrk = fmHit.at(trkFromPFP.at(0).key());

      //art::FindManyP<recob::SpacePoint> fmSps(hitsFromTrk, e, fPFParticleModuleLabels[iCryo]);
      //if ( !fmSps.isValid() ){
      //  throw cet::exception("HARPS") << "Module wants to use invalid associations (fmSps)." << std::endl;
      //}

      // Copied from my overlay module
      for ( auto const& iHitPtr : hitsFromTrk ) {
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

        // TODO: add some sort of random ability to miss some wires.

        // If we're not skipping this hit then write it to the file!
        recob::Hit theHit = *iHitPtr;
        if ( fTPCHitsWireAssn ) {
          std::vector< art::Ptr<recob::Wire> > hitWires = fmWire.at(iHitPtr.key());
          if ( hitWires.size() == 0 ) throw cet::exception("HARPS") << "Hit found with no associated wires...\n";
          else if ( hitWires.size() > 1 ) mf::LogWarning("HARPS") << "Hit with >1 recob::Wire associated...";
          hitCol.emplace_back(theHit,hitWires[0]);
        }
        else hitCol.emplace_back(theHit);
      } // loop hits
    } // loop pfps
  } // loop cryos

  hitCol.put_into(e);
}

DEFINE_ART_MODULE(HARPS)
