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
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// For random implementation: see
// https://github.com/SBNSoftware/icaruscode/blob/7b71f0213a1a66913cfe1cd795f49a9215570145/icaruscode/PMT/OpReco/FakePhotoS_module.cc#L16
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include <numeric>

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

  double FindNearestTrajPointRR( const std::map<std::tuple<double,double,double>, double>& trjPtToRR, const double* in_xyz );

  typedef std::pair<unsigned int, float> WireTimePair_t;
  typedef std::pair<WireTimePair_t, WireTimePair_t> PairOfWireTimePair_t;

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

  bool fTagDaughters; ///< Are we tagging daughter particles?
  bool fKeepContext;  ///< Are we doing the inverse operation (ie keeping only the context)?
  bool fUseTrajPointForDist; ///< Are we making a map for track trajectory points to start/end distance to use
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
  fFlatEngine                ( art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "Gen", p, "Seed") ),
  fTagDaughters              ( p.get< bool >("TagDaughters", false) ),
  fKeepContext               ( p.get< bool >("KeepContext", false) ),
  fUseTrajPointForDist       ( p.get< bool >("UseTrajPointForDist", false) )
{
  if ( fPFParticleModuleLabels.size()!=fTrackModuleLabels.size() || fPFParticleModuleLabels.size()!=fHitModuleLabels.size() )
    throw cet::exception("HARPS") << "Error... InputTag vectors need to be same size..." << std::endl; // borrow from elsewhere

  if ( fMaskBeginningDist>0. && fMaskEndDist>0. )
    throw cet::exception("HARPS") << "Error... Should not try to mask the beginning and ALSO some distance from the end..." << std::endl;

  //////
  // TODO: there is something a little strange here with saving the daughters and keeping the end or start dist... FOR NOW, if
  //       we do fTagDaughters, do *NOT* allow fMaskBeginningDist or fMaskEndDist... and in this case the vertex we want to return
  //       is the parent particle's start/vertex position.
  if ( (fKeepContext || fTagDaughters) && (fMaskBeginningDist>0. || fMaskEndDist>0. || fUseTrajPointForDist) )
    throw cet::exception("HARPS") << "Error... something strange about trying to keep daughters and remove based on start/end positions..." << std::endl;

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
    std::map< std::string, std::vector<size_t> >& aliasParticleList = (nC == 0) ? fParticleListCryo0 : fParticleListCryo1;
    if ( fOnlyOnePerEvent && (aliasParticleList.find(evtID) != aliasParticleList.end()) )
      continue;
    aliasParticleList[evtID].push_back(nP);
  }
  in.close();

  fFlatRand = new CLHEP::RandFlat(fFlatEngine,0.,1.);

  produces<std::vector<recob::Vertex>>();
}

void HARPS::produce(art::Event& e)
{
  // As in Gauss Hit Finder code (larreco/HitFinder/GausHitFinder_module.cc)
  recob::HitCollectionCreator hitCol(e, fTPCHitCreatorInstanceName, fTPCHitsWireAssn, false);

  // Make the vector of vertices: if using fKeepContext then we will output an EMPTY vector of vertices
  std::unique_ptr< std::vector<recob::Vertex> > vtxVec( std::make_unique< std::vector<recob::Vertex> >() );

  // Return empty hit container if no particles of interest.
  std::string evtID = std::to_string(e.run()) + ":" + std::to_string(e.subRun()) + ":" + std::to_string(e.event());
  if ( fParticleListCryo0.find(evtID) == fParticleListCryo0.end() &&
       fParticleListCryo1.find(evtID) == fParticleListCryo1.end() ) {
    e.put(std::move(vtxVec));
    hitCol.put_into(e);
    return;
  }

  //std::cout << "!!!! --- Got here 1" << std::endl;

  // Load in the PFParticles, Tracks, and Hits & the necessary FindManyP's
  for ( unsigned int iCryo=0; iCryo<2; ++iCryo ) {
    std::map< std::string, std::vector<size_t> >& aliasParticleList = (iCryo == 0) ? fParticleListCryo0 : fParticleListCryo1;
    if ( aliasParticleList.find(evtID) == aliasParticleList.end() ) continue;

    //std::cout << "Currently looking at Cryostat " << iCryo << " which has " << aliasParticleList[evtID].size() << " particles." << std::endl;

    art::InputTag prtLabel = (fPFParticleModuleLabels.size() > 1 ? fPFParticleModuleLabels[iCryo] : fPFParticleModuleLabels[0]);
    art::InputTag trkLabel = (fTrackModuleLabels.size() > 1 ? fTrackModuleLabels[iCryo] : fTrackModuleLabels[0]);
    art::InputTag hitLabel = (fHitModuleLabels.size() > 1 ? fHitModuleLabels[iCryo] : fHitModuleLabels[0]);

    // as in previous code
    art::Handle< std::vector<recob::PFParticle> > pfpHandle;
    std::vector< art::Ptr<recob::PFParticle> > pfps;
    if ( e.getByLabel(prtLabel, pfpHandle) ) {
      art::fill_ptr_vector(pfps, pfpHandle);
    }
    else {
      mf::LogError("HARPS") << "Error pulling in PFParticle product... Skipping cryostat.";
    }

    // How many/check
    //std::cout << "!!!! --- Got here 2 --> there are " << pfps.size() << " pfps." << std::endl;

    // Get the find many for vertices, we only really want this if fTagDaughters is true at moment...
    art::FindManyP<recob::Vertex> fmVtx(pfpHandle, e, prtLabel);
    if ( !fmVtx.isValid() ){
      throw cet::exception("HARPS") << "Module wants to use invalid association to recob::Vertex." << std::endl;
    }

    // sort the PFParticles by id
    // annoying as recob::PFParticle has a < opperator, but we have art::Ptr<recob::PFParticle> 
    std::sort(pfps.begin(), pfps.end(), [](art::Ptr<recob::PFParticle> pfpA, art::Ptr<recob::PFParticle> pfpB){ return pfpA->Self() < pfpB->Self(); });

    // if we are tagging daughters find them here
    if (fTagDaughters)
    {
      // if we're going to tag daughters, save the parent particles' vertices first
      for ( auto const& pfpId : aliasParticleList[evtID] ) {
        const std::vector< art::Ptr<recob::Vertex> > vtxFromPFP = fmVtx.at(pfps[pfpId]->Self());
        if ( vtxFromPFP.size() != 1 ) continue;
        double vtxPos[3] = {vtxFromPFP.at(0)->position().X(), vtxFromPFP.at(0)->position().Y(), vtxFromPFP.at(0)->position().Z()};
        vtxVec->push_back( recob::Vertex(vtxPos) );
      }

      auto idItr = aliasParticleList[evtID].begin();
      while (idItr != aliasParticleList[evtID].end())
      {
        std::vector<size_t> dghts = pfps[*idItr]->Daughters();
        aliasParticleList[evtID].insert(aliasParticleList[evtID].end(), dghts.begin(), dghts.end());
        ++idItr;
      }
    }

    // if we are keeping context, get the complement of our indices
    if (fKeepContext)
    {
      // start with {0, 1, ..., pfps.size() - 1}
      std::vector<size_t> complement(pfps.size());
      std::iota(complement.begin(), complement.end(), 0);

      // make sure the initial pfp ids are sorted
      // we want descending order so when we erase from complement we don't effect the indexing
      std::sort(aliasParticleList[evtID].begin(), aliasParticleList[evtID].end(), [](size_t a, size_t b){ return a > b; });

      // erase the indices in aliasParticleList[evtID] from complement
      for (auto const& id : aliasParticleList[evtID])
        complement.erase(complement.begin() + id);
      
      // set the particle list to the complement
      aliasParticleList[evtID] = complement;
    }

    art::Handle< std::vector<recob::Track> > trkHandle;
    std::vector< art::Ptr<recob::Track> > trks;
    if ( e.getByLabel(trkLabel, trkHandle) ) {
      art::fill_ptr_vector(trks, trkHandle);
    }
    else {
      mf::LogError("HARPS") << "Error pulling in Track product... Skipping cryostat.";
    }

    art::Handle< std::vector<recob::Hit> > hitHandle;
    std::vector< art::Ptr<recob::Hit> > hits;
    if ( e.getByLabel(hitLabel, hitHandle) ) {
      art::fill_ptr_vector(hits, hitHandle);
    }
    else {
      mf::LogError("HARPS") << "Error pulling in Hit product... Skipping cryostat.";
    }

    // For the Wire Assns in HitCreator
    art::FindManyP<recob::Wire> fmWire(hitHandle, e, hitLabel);
    if ( !fmWire.isValid() && fTPCHitsWireAssn ){
      throw cet::exception("HARPS") << "Module wants to use Hit-Wire Associations, but fmWire invalid." << std::endl;
    }

    art::FindManyP<recob::Track>      fmTrk(pfpHandle, e, trkLabel);
    art::FindManyP<recob::Hit>        fmHit(trkHandle, e, trkLabel);
    art::FindManyP<recob::SpacePoint> fmSps(hitHandle, e, prtLabel);
    if ( !fmTrk.isValid() || !fmHit.isValid() || !fmSps.isValid() ){
      throw cet::exception("HARPS") << "Module wants to use invalid associations." << std::endl;
    }

    //std::cout << "!!!! --- Got here 3 --> About to loop over the pfpIds... We have "
    //          << aliasParticleList[evtID].size() << " to care about." << std::endl;

    // Find the PFParticles we care about and put their hits into the collection
    for (auto const& id : aliasParticleList[evtID]) {
      art::Ptr<recob::PFParticle> pfp = pfps[id];

      const std::vector< art::Ptr<recob::Track> > trkFromPFP = fmTrk.at(pfp.key());
      if ( trkFromPFP.size() != 1 ) {
        mf::LogError("HARPS") << "fmTrk gave us more (or less) than one track for this PFParticle... Skipping.";
        continue;
      }

      // If using Traj Points for distance calculation (via RR)
      std::vector<recob::tracking::Point_t> trkFilterTrajPt;
      std::vector<double> trkTrajPtLengths;
      std::map<std::tuple<double,double,double>, double> mapTrajPtEndDist;

      if ( fUseTrajPointForDist ) {
        // First let's make a vector that is iterable (I think we can assume they are sorted?)
        double prev_x(0.), prev_y(0.), prev_z(0.);
        for ( unsigned int idxTrajPt=0; idxTrajPt < trkFromPFP.at(0)->NumberTrajectoryPoints(); ++idxTrajPt ) {
          if ( !trkFromPFP.at(0)->HasValidPoint(idxTrajPt) ) continue;
          trkFilterTrajPt.push_back( trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position );
          if ( idxTrajPt == 0 ) {
            trkTrajPtLengths.push_back(0.);
            prev_x = trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position.x();
            prev_y = trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position.y();
            prev_z = trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position.z();
          }
          else {
            double this_x = trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position.x();
            double this_y = trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position.y();
            double this_z = trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position.z();
            trkTrajPtLengths.push_back( std::hypot( this_x-prev_x, this_y-prev_y, this_z-prev_z) );
            prev_x = this_x;
            prev_y = this_y;
            prev_z = this_z;
          }
        }
        // Now let's build the maps
        unsigned int idxTrajPt = 0;
        for ( std::vector<double>::iterator it=trkTrajPtLengths.begin(); it!=trkTrajPtLengths.end(); ++it ) {
          const double rr = std::accumulate( it, trkTrajPtLengths.end(), 0. );
          //const recob::tracking::Point_t trajPtPoint = trkFromPFP.at(0)->TrajectoryPoint(idxTrajPt).position;
          const std::tuple<double,double,double> trajPtPos = {trkFilterTrajPt[idxTrajPt].x(),trkFilterTrajPt[idxTrajPt].y(),trkFilterTrajPt[idxTrajPt].z()};
          if ( mapTrajPtEndDist.find(trajPtPos) == mapTrajPtEndDist.end() ) {
            mapTrajPtEndDist[trajPtPos] = rr;
          }
          else {
            std::cout << "Found same traj point already? SKIPPING BUT PRINTING OUT THAT THIS HAS HAPPENED..." << std::endl;
          }
          idxTrajPt+=1;
        }
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

      // map for the interpolation
      // first elem is PlaneID --> combination of cryo, tpc, plane
      // second elem is pair of <<lo wire, lo time>, <hi wire, hi time>>
      std::map< geo::PlaneID, PairOfWireTimePair_t > wireToTimeRangeMap;

      std::map< unsigned int, geo::PlaneID > planeNumToPlaneID_ForMaskBeginningDist; // planeID for begin point from which we'll search a box
      std::map< unsigned int, WireTimePair_t > planeNumToWireTime_ForMaskBeginningDist; // WireTimePair_t for begin point from which we'll search a box
      std::map< unsigned int, double > planeNumToPtDist_ForMaskBeginningDist; // min distance for begin point from which we'll search a box
      std::map< unsigned int, geo::PlaneID > planeNumToPlaneID_ForMaskEndDist; // planeID for end point from which we'll search a box
      std::map< unsigned int, WireTimePair_t > planeNumToWireTime_ForMaskEndDist; // WireTimePair_t for end point from which we'll search a box
      std::map< unsigned int, double > planeNumToPtDist_ForMaskEndDist; // min distance for end point from which we'll search a box

      // Vertex coordinates: only used if fMaskBeginningDist or fMaskEndDist are > 0.
      double minmaxDist = ( fMaskBeginningDist > 0. ? 9999. : ( fMaskEndDist > 0. ? -9999. :  0.) );
      double vtxX = -9999.;
      double vtxY = -9999.;
      double vtxZ = -9999.;

      // Loop as in my overlay module, copy hits passing the criteria into the maps for checking against other criteria...
      for( auto const& iHitPtr : hitsFromTrk ) {
        // Check if we want to remove this hit because we are masking X distance from the beginning:
        // TODO: as a start, use the Pandora spacepoints<->hit association... is this safe or should we do something in 2d distances?
        // UPDATE: loop through twice and save the max and min wires (and times?) the first time for things with a match and use these to
        //         save some of the 0 match cases the second time through via interpolation... only in the MaskBeginning or MaskEnd case...
        if ( fMaskBeginningDist > 0. ) {
          const std::vector< art::Ptr<recob::SpacePoint> > spsFromHit = fmSps.at(iHitPtr.key());
          if ( spsFromHit.size() != 1 ) {
            if ( spsFromHit.size() > 1 ) mf::LogError("HARPS") << "fmSps gave us something more than one spacepoint for this hit (" << spsFromHit.size() << ")... Skipping.";
            continue;
          }
          const auto thisXYZ = spsFromHit.at(0)->XYZ();
          if ( std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] ) < fMaskBeginningDist ) continue;

          // if this spacepoint is closest to the track start, save it as our 3D Vertex
          if ( std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] ) < minmaxDist ) {
            minmaxDist = std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] );
            vtxX = thisXYZ[0];
            vtxY = thisXYZ[1];
            vtxZ = thisXYZ[2];
          }

          // Check wireToTimeRangeMap
          const auto planeid = iHitPtr->WireID().asPlaneID();
          const unsigned int thisplane = iHitPtr->WireID().Plane;
          const unsigned int thiswire = iHitPtr->WireID().Wire;
          const float thistime = iHitPtr->PeakTime();
          if ( wireToTimeRangeMap.count(planeid) == 0 ) {
            wireToTimeRangeMap[planeid] = std::make_pair<WireTimePair_t, WireTimePair_t>( std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime()),
                                                                                          std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime()) );
          }
          else {
            if ( thiswire < wireToTimeRangeMap[planeid].first.first ) {
              wireToTimeRangeMap[planeid].first.first = thiswire;
              wireToTimeRangeMap[planeid].first.second = thistime;
            }
            else if ( thiswire == wireToTimeRangeMap[planeid].first.first && thistime < wireToTimeRangeMap[planeid].first.second ) {
              wireToTimeRangeMap[planeid].first.second = thistime;
            }
            else if ( thiswire == wireToTimeRangeMap[planeid].second.first && thistime > wireToTimeRangeMap[planeid].second.second ) {
              wireToTimeRangeMap[planeid].second.second = thistime;
            }
            if ( thiswire > wireToTimeRangeMap[planeid].second.first ) {
              wireToTimeRangeMap[planeid].second.first = thiswire;
              wireToTimeRangeMap[planeid].second.second = thistime;
            }
          } // filling up the hit min/max map

          // Also deal with where to put a box to search near the beginning of track for hits unassociated with spacepoints
          if ( planeNumToPtDist_ForMaskBeginningDist.count(thisplane) == 0 ) {
            planeNumToPtDist_ForMaskBeginningDist[thisplane] = std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] );
            planeNumToPlaneID_ForMaskBeginningDist[thisplane] = planeid;
            planeNumToWireTime_ForMaskBeginningDist[thisplane] = std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime());
          }
          else if ( std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] ) < planeNumToPtDist_ForMaskBeginningDist[thisplane] ) {
            planeNumToPtDist_ForMaskBeginningDist[thisplane] = std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] );
            planeNumToPlaneID_ForMaskBeginningDist[thisplane] = planeid;
            planeNumToWireTime_ForMaskBeginningDist[thisplane] = std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime());
          }
        } // fMaskBeginningDist
        else if ( fMaskEndDist > 0. ) {
          const std::vector< art::Ptr<recob::SpacePoint> > spsFromHit = fmSps.at(iHitPtr.key());
          if ( spsFromHit.size() != 1 ) {
            if ( spsFromHit.size() > 1 ) mf::LogError("HARPS") << "fmSps gave us something more than one spacepoint for this hit (" << spsFromHit.size() << ")... Skipping.";
            continue;
          }
          const auto thisXYZ = spsFromHit.at(0)->XYZ();

          if ( !fUseTrajPointForDist ){
            // Just uses stright-line distance
            if ( std::hypot( trkEndX-thisXYZ[0], trkEndY-thisXYZ[1], trkEndZ-thisXYZ[2] ) > fMaskEndDist ) continue;

            // if this spacepoint is furthest from the track end, save it as our 3D Vertex
            if ( std::hypot( trkEndX-thisXYZ[0], trkEndY-thisXYZ[1], trkEndZ-thisXYZ[2] ) > minmaxDist ) {
              minmaxDist = std::hypot( trkEndX-thisXYZ[0], trkEndY-thisXYZ[1], trkEndZ-thisXYZ[2] );
              vtxX = thisXYZ[0];
              vtxY = thisXYZ[1];
              vtxZ = thisXYZ[2];
            }
          }
          else {
            // Uses the track's residual range by finding the closest traj point to this spacepoint
            double nearestTrajPtRR = FindNearestTrajPointRR(mapTrajPtEndDist,thisXYZ);
            if ( nearestTrajPtRR > fMaskEndDist || nearestTrajPtRR < 0. ) continue;
            if ( nearestTrajPtRR > minmaxDist ) {
              minmaxDist = nearestTrajPtRR;
              vtxX = thisXYZ[0];
              vtxY = thisXYZ[1];
              vtxZ = thisXYZ[2];
            }
          }

          // Check wireToTimeRangeMap
          const auto planeid = iHitPtr->WireID().asPlaneID();
          const unsigned int thisplane = iHitPtr->WireID().Plane;
          const unsigned int thiswire = iHitPtr->WireID().Wire;
          const float thistime = iHitPtr->PeakTime();
          if ( wireToTimeRangeMap.count(planeid) == 0 ) {
            wireToTimeRangeMap[planeid] = std::make_pair<WireTimePair_t, WireTimePair_t>( std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime()),
                                                                                          std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime()) );
          }
          else {
            if ( thiswire < wireToTimeRangeMap[planeid].first.first ) {
              wireToTimeRangeMap[planeid].first.first = thiswire;
              wireToTimeRangeMap[planeid].first.second = thistime;
            }
            else if ( thiswire == wireToTimeRangeMap[planeid].first.first && thistime < wireToTimeRangeMap[planeid].first.second ) {
              wireToTimeRangeMap[planeid].first.second = thistime;
            }
            else if ( thiswire == wireToTimeRangeMap[planeid].second.first && thistime > wireToTimeRangeMap[planeid].second.second ) {
              wireToTimeRangeMap[planeid].second.second = thistime;
            }
            if ( thiswire > wireToTimeRangeMap[planeid].second.first ) {
              wireToTimeRangeMap[planeid].second.first = thiswire;
              wireToTimeRangeMap[planeid].second.second = thistime;
            }
          } // filling up the hit min/max map

          // Also deal with where to put a box to search near the beginning of track for hits unassociated with spacepoints
          if ( planeNumToPtDist_ForMaskEndDist.count(thisplane) == 0 ) {
            planeNumToPtDist_ForMaskEndDist[thisplane] = std::hypot( trkEndX-thisXYZ[0], trkEndY-thisXYZ[1], trkEndZ-thisXYZ[2] );
            planeNumToPlaneID_ForMaskEndDist[thisplane] = planeid;
            planeNumToWireTime_ForMaskEndDist[thisplane] = std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime());
          }
          else if ( std::hypot( trkEndX-thisXYZ[0], trkEndY-thisXYZ[1], trkEndZ-thisXYZ[2] ) < planeNumToPtDist_ForMaskEndDist[thisplane] ) {
            planeNumToPtDist_ForMaskEndDist[thisplane] = std::hypot( trkStartX-thisXYZ[0], trkStartY-thisXYZ[1], trkStartZ-thisXYZ[2] );
            planeNumToPlaneID_ForMaskEndDist[thisplane] = planeid;
            planeNumToWireTime_ForMaskEndDist[thisplane] = std::make_pair<unsigned int, float>(iHitPtr->WireID().Wire, iHitPtr->PeakTime());
          }
        } // fMaskEndDist

        unsigned int plane = iHitPtr->WireID().Plane;
        unsigned int idx = 5000*4*3*iHitPtr->WireID().Cryostat + 5000*3*iHitPtr->WireID().TPC + 5000*iHitPtr->WireID().Plane + iHitPtr->WireID().Wire;
        planeToHitsMap[plane].push_back( iHitPtr );
        planeToIdxMap[plane].emplace( idx );
      } // loop hits

      // If we are masking by distance we may have skipped hits due to not being associated with a spacepoint...
      // Let's try to add those back to the output by interpolating wires/times
      if ( fMaskBeginningDist > 0. || fMaskEndDist > 0. ) {
        for( auto const& iHitPtr : hitsFromTrk ) {
          const std::vector< art::Ptr<recob::SpacePoint> > spsFromHit = fmSps.at(iHitPtr.key());
          if ( spsFromHit.size() != 0 ) {
            continue; // only want hits without spacepoints here
          }
          // The info to compare to map
          const auto planeid = iHitPtr->WireID().asPlaneID();
          const unsigned int thisplane = iHitPtr->WireID().Plane;
          const unsigned int thiswire = iHitPtr->WireID().Wire;
          const float thistime = iHitPtr->PeakTime();
          if ( wireToTimeRangeMap.count(planeid) == 0 ) {
            continue; // nothing on this plane in the map
          }
          // Check if it's within the interpolation and add it to the save map if so
          if ( thiswire >= wireToTimeRangeMap[planeid].first.first && thistime >= wireToTimeRangeMap[planeid].first.second &&
               thiswire <= wireToTimeRangeMap[planeid].second.first && thistime <= wireToTimeRangeMap[planeid].second.second ) {
            unsigned int plane = iHitPtr->WireID().Plane;
            unsigned int idx = 5000*4*3*iHitPtr->WireID().Cryostat + 5000*3*iHitPtr->WireID().TPC + 5000*iHitPtr->WireID().Plane + iHitPtr->WireID().Wire;
            planeToHitsMap[plane].push_back( iHitPtr );
            planeToIdxMap[plane].emplace( idx );
          } // check on interpolation
          // Also check in a box near the beginning (if fMaskBeginningDist) or the end (if fMaskEndDist): let's start with ~3 cm
          // --> so a box of about 10 wires and 50 ticks
          else if ( ( fMaskBeginningDist > 0. &&
                      planeid == planeNumToPlaneID_ForMaskBeginningDist[thisplane] &&
                      abs(int(thiswire)-int(planeNumToWireTime_ForMaskBeginningDist[thisplane].first)) < 10 &&
                      fabs(thistime-planeNumToWireTime_ForMaskBeginningDist[thisplane].second) < 50. ) ||
                    ( fMaskEndDist > 0. &&
                      planeid == planeNumToPlaneID_ForMaskEndDist[thisplane] &&
                      abs(int(thiswire)-int(planeNumToWireTime_ForMaskEndDist[thisplane].first)) < 10 &&
                      fabs(thistime-planeNumToWireTime_ForMaskEndDist[thisplane].second) < 50. ) ){
            unsigned int plane = iHitPtr->WireID().Plane;
            unsigned int idx = 5000*4*3*iHitPtr->WireID().Cryostat + 5000*3*iHitPtr->WireID().TPC + 5000*iHitPtr->WireID().Plane + iHitPtr->WireID().Wire;
            planeToHitsMap[plane].push_back( iHitPtr );
            planeToIdxMap[plane].emplace( idx );
          } // check "extrapolation box"
        } // loop hits

        // PUT THE 3D VERTEX BACK INTO EVENT
        double vtxPos[3] = {vtxX, vtxY, vtxZ};
        vtxVec->push_back( recob::Vertex(vtxPos) );
      } // only run if masking by distance

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

  e.put(std::move(vtxVec));
  hitCol.put_into(e);
}

double HARPS::FindNearestTrajPointRR( const std::map<std::tuple<double,double,double>, double>& trjPtToRR, const double* in_xyz )
{
  std::cout << "The hit spacepoint at (" << in_xyz[0] << ", " << in_xyz[1] << ", " << in_xyz[2] << ") ";

  if ( trjPtToRR.empty() ) {
    std::cout << " ... has no match." << std::endl;
    return -5.;
  }

  double minDistToPt = std::numeric_limits<double>::max();
  double rrAtMinDist = -5.;
  //double xAtMinDist = -5.;
  //double yAtMinDist = -5.;
  //double zAtMinDist = -5.;

  for ( auto const& [pt, rr] : trjPtToRR ) {
    double ptx = std::get<0>(pt);
    double pty = std::get<1>(pt);
    double ptz = std::get<2>(pt);
    if ( std::hypot(in_xyz[0]-ptx, in_xyz[1]-pty, in_xyz[2]-ptz) < minDistToPt ) {
      minDistToPt = std::hypot(in_xyz[0]-ptx, in_xyz[1]-pty, in_xyz[2]-ptz);
      rrAtMinDist = rr;
      //xAtMinDist = ptx;
      //yAtMinDist = pty;
      //zAtMinDist = ptz;
    }
  }

  //std::cout << " ... has matched traj point at (" << xAtMinDist << ", " << yAtMinDist << ", " << zAtMinDist << "), with rr " << rrAtMinDist << std::endl;

  return rrAtMinDist;
}


DEFINE_ART_MODULE(HARPS)
