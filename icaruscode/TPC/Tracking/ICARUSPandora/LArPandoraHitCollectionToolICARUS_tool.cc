/**
 *  @file  larpandora/LArPandoraInterface/Tools/LArPandoraHitCollectionToolICARUS_tool.cc
 * 
 *  @brief Define class for hit collection tools and implements ICARUS tool, with a duplication filter
 * 
 */

#include "art/Utilities/ToolMacros.h"

#include "larpandora/LArPandoraInterface/LArPandoraHitCollectionTool.h"

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/Hit.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <map>

namespace lar_pandora {

  class LArPandoraHitCollectionToolICARUS : public IHitCollectionTool
  {
  public:
    explicit LArPandoraHitCollectionToolICARUS(const fhicl::ParameterSet& pset);
    void CollectHits(const art::Event& evt, const std::string& label, HitVector& hitVector) override;
  };

  LArPandoraHitCollectionToolICARUS::LArPandoraHitCollectionToolICARUS(const fhicl::ParameterSet& pset){}

  void LArPandoraHitCollectionToolICARUS::CollectHits(const art::Event& evt, const std::string& label, HitVector& hitVector)
  {
    std::map< std::pair<raw::ChannelID_t, int>, unsigned int > channelTimeHits; // converting the float PeakTime to int for the second item in pair

    art::Handle<std::vector<recob::Hit>> theHits;
    evt.getByLabel(label, theHits);

    if (!theHits.isValid()) {
      mf::LogDebug("LArPandoraHitCollectionToolICARUS") << "  Failed to find hits... " << std::endl;
      return;
    }
    else {
      mf::LogDebug("LArPandoraHitCollectionToolICARUS") << "  Found: " << theHits->size() << " Hits " << std::endl;
    }

    for (unsigned int i = 0; i < theHits->size(); ++i) {
      const art::Ptr<recob::Hit> hit(theHits, i);
      if ( std::exchange( channelTimeHits[std::make_pair(hit->Channel(),int(hit->PeakTime()))], 1 ) ) continue;
      hitVector.push_back(hit);
    }

    mf::LogDebug("LArPandoraHitCollectionToolICARUS") << "  Passing along " << hitVector.size() << " Hits " << std::endl;
  }

} // namespace lar_pandora

DEFINE_ART_CLASS_TOOL(lar_pandora::LArPandoraHitCollectionToolICARUS)
