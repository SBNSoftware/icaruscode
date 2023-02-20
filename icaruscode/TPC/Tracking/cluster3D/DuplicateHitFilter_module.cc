////////////////////////////////////////////////////////////////////////
// Class:       DuplicateHitFilter
// Plugin Type: producer (Unknown Unknown)
// File:        DuplicateHitFilter_module.cc
//
// Generated at Fri Feb 17 00:02:19 2023 by Bruce Howard using cetskelgen
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
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"

#include <memory>

class DuplicateHitFilter;


class DuplicateHitFilter : public art::EDProducer {
public:
  explicit DuplicateHitFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DuplicateHitFilter(DuplicateHitFilter const&) = delete;
  DuplicateHitFilter(DuplicateHitFilter&&) = delete;
  DuplicateHitFilter& operator=(DuplicateHitFilter const&) = delete;
  DuplicateHitFilter& operator=(DuplicateHitFilter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::string   fHitCreatorInstanceName;
  bool          fHitWireAssns;
  art::InputTag fInputHitLabel;

  int           fDebugLevel;
};


DuplicateHitFilter::DuplicateHitFilter(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fHitCreatorInstanceName ( p.get< std::string >("HitCreatorInstanceName","") ),
    fHitWireAssns           ( p.get<bool>("HitWireAssns", false) ),
    fInputHitLabel          ( p.get< art::InputTag >("InputHitLabel", "cluster3D") ),
    fDebugLevel             ( p.get<int>("DebugLevel", 0) )
{
  recob::HitCollectionCreator::declare_products(producesCollector(), fHitCreatorInstanceName, fHitWireAssns, false);
}

void DuplicateHitFilter::produce(art::Event& e)
{
  // As in Gauss hit finder code (larreco/HitFinder/GausHitFinder_module.cc) and copying from the application in other modules/test modules
  recob::HitCollectionCreator hitCol( e, fHitCreatorInstanceName, fHitWireAssns, false);

  art::Handle< std::vector<recob::Hit> > hitsHandle;
  std::vector< art::Ptr<recob::Hit> > hits;
  if ( e.getByLabel(fInputHitLabel,hitsHandle) )
    art::fill_ptr_vector(hits,hitsHandle);

  art::FindManyP<recob::Wire> fmwire(hitsHandle, e, fInputHitLabel);
  if ( !fmwire.isValid() && fHitWireAssns ) {
    mf::LogError("DuplicateHitFilter") << "Error in validity of fmwire when specifying to save the association. Returning.";
    return;
  }

  // converting the float peaktime to int for second item in pair. Is this specific enough or should we use other hit info?
  std::map< std::pair<raw::ChannelID_t, int>, unsigned int > channelTimeHits;

  for ( auto const& hit : hits ) {
    recob::Hit theHit = *hit;

    if ( channelTimeHits.find( std::make_pair(theHit.Channel(), int(theHit.PeakTime())) ) == channelTimeHits.end() )
      channelTimeHits[ std::make_pair(theHit.Channel(), int(theHit.PeakTime())) ] = 1;
    else {
      if ( fDebugLevel > 1 ) mf::LogWarning("DuplicateHitFilter") << "Duplicate hit found, not considering it moving forward.";
      continue;
    }

    if ( fHitWireAssns ) {
      std::vector< art::Ptr<recob::Wire> > hitWires = fmwire.at( hit.key() );
      if ( hitWires.size() == 0 || hitWires.size() > 1 ) {
	mf::LogWarning("DuplicateHitFilter") << "No associated wires when they were expected, or > 1 associated wire. Returning.";
	return;
      }
      hitCol.emplace_back( theHit, hitWires[0] );
    }
    else hitCol.emplace_back( theHit );
  }

  if ( fDebugLevel > 0 ) std::cout << "INITIAL HITS = " << hits.size() << " ... KEPT = " << hitCol.size() << std::endl;

  hitCol.put_into(e);
}

DEFINE_ART_MODULE(DuplicateHitFilter)
