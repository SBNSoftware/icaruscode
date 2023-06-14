////////////////////////////////////////////////////////////////////////
//
// FilterDataIntegrity class
//
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "artdaq-core/Data/Fragment.hh"

#include <set>
#include <string>

///filters for events, etc
namespace filter
{

    class FilterDataIntegrity : public art::EDFilter
    {

    public:
        explicit FilterDataIntegrity(fhicl::ParameterSet const &);

        bool filter(art::Event &evt) override;

    private:

    // We will use this to keep track of the expected fragments in an event
    // but note we will assume the first event is complete to set this list
    std::set<int> fExpectedFragments;

    }; //class FilterDataIntegrity
}

///////////////////////////////////////////////////////

filter::FilterDataIntegrity::FilterDataIntegrity(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
{
    return;
}

bool filter::FilterDataIntegrity::filter(art::Event &event)
{
    bool filterPass = true;

	// get all the artdaq fragment collections in the event.
	std::vector<art::Handle<std::vector<artdaq::Fragment>>> fragmentHandles;
	fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();

    std::set<int> missingFragments(fExpectedFragments);
    int           emptyFragments(0);

    for (const auto& handle : fragmentHandles)
    {
        for (const auto& fragment : *handle)
        {
            int fragmentID = fragment.fragmentID();

            fExpectedFragments.insert(fragmentID);
            missingFragments.erase(fragmentID);

            std::string instanceName = handle.provenance()->productInstanceName();
            std::size_t found        = instanceName.find("Empty");

            if (found != std::string::npos) emptyFragments++;
        }
    }

    mf::LogDebug("FilterDataIntegrity") << "Expected fragments: " << fExpectedFragments.size() << " expected, has " << missingFragments.size() << " missing. There are " << emptyFragments << " empty fragments";

    if (!missingFragments.empty() || emptyFragments > 0)
    {
        mf::LogInfo("FilterDataIntegrity") << "Bad fragments: " << missingFragments.size() << ", " << emptyFragments;

        filterPass = false;
    }

    return filterPass;
}

namespace filter
{

    DEFINE_ART_MODULE(FilterDataIntegrity)

} //namespace filt
