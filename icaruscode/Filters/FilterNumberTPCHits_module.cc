////////////////////////////////////////////////////////////////////////
//
// FilterNumberTPCHits class
//
////////////////////////////////////////////////////////////////////////
#include <fstream>

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/Hit.h"

///filters for events, etc
namespace filter
{

    class FilterNumberTPCHits : public art::EDFilter
    {

    public:
        explicit FilterNumberTPCHits(fhicl::ParameterSet const &);

        bool filter(art::Event &evt) override;

    private:
        std::vector<art::InputTag> fHitDataLabelVec;
        unsigned int               fMaximumHits;

    }; //class FilterNumberTPCHits
}

///////////////////////////////////////////////////////

filter::FilterNumberTPCHits::FilterNumberTPCHits(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
{
    fHitDataLabelVec = pset.get<std::vector<art::InputTag>>("HitDataLabelVec",    {""});
    fMaximumHits     = pset.get<unsigned int              >("MaximumHits",     800000u);

    return;
}

bool filter::FilterNumberTPCHits::filter(art::Event &event)
{
    bool filterPass = true;

    for(auto const& hitDataLabel : fHitDataLabelVec)
    {
       auto const& hitData = event.getProduct<std::vector<recob::Hit>>(hitDataLabel);

//       std::cout << "FilterNumberTPCHits: label: " << hitDataLabel << " has " << hitData.size() << " hits (rejection is " << fMaximumHits << ")" << std::endl;

        if (hitData.size() > fMaximumHits)
        {
            mf::LogInfo log("FilterNumberTPCHits");
            log << "******************************************************\n" << "Rejecting event for " 
                << "***** " << hitDataLabel << " with " << hitData.size() << " hits ******\n"
                << "******************************************************";

            filterPass = false;
            break;
        }
    }

    return filterPass;
}

namespace filter
{

    DEFINE_ART_MODULE(FilterNumberTPCHits)

} //namespace filt
