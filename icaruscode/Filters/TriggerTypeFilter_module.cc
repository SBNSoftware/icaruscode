////////////////////////////////////////////////////////////////////////
//
// TriggerTypeFilter class
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
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger
#include "icaruscode/Decode/BeamBits.h"

///filters for events, etc
namespace filter
{

    class TriggerTypeFilter : public art::EDFilter
    {

    public:
        explicit TriggerTypeFilter(fhicl::ParameterSet const &);

        bool filter(art::Event &evt) override;

    private:
        art::InputTag fTriggerDataLabel;
        std::string   fTriggerType;
        unsigned int  fTriggerBits;

    }; //class TriggerTypeFilter
}

///////////////////////////////////////////////////////

filter::TriggerTypeFilter::TriggerTypeFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
{
    fTriggerDataLabel = pset.get<art::InputTag>("TriggerDataLabel", "daqTrigger");
    fTriggerType      = pset.get<std::string  >("TriggerType",             "BNB");

    std::cout << "****> Initializing trigger type " << fTriggerType << std::endl;

    if (fTriggerType == "BNB")  fTriggerBits = value(sbn::triggerSource::BNB);
    if (fTriggerType == "NuMI") fTriggerBits = value(sbn::triggerSource::NuMI);

    return;
}

bool filter::TriggerTypeFilter::filter(art::Event &event)
{
    bool filterPass = false;

    auto const& triggerVec = event.getByLabel<std::vector<raw::Trigger>>(fTriggerDataLabel);

    if (!triggerVec.empty())
    {
        const raw::Trigger& trigger = triggerVec.front();

        if (trigger.Triggered(fTriggerBits))
        {
            std::cout << "***> Trigger type matched, type is " << fTriggerType << ", " << fTriggerBits << std::endl;
            filterPass = true;
        }
    }

    return filterPass;
}

namespace filter
{

    DEFINE_ART_MODULE(TriggerTypeFilter)

} //namespace filt
