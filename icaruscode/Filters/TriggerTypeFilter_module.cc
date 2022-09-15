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
#include "sbnobj/Common/Trigger/BeamBits.h"

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

    if (fTriggerType == "BNB")            fTriggerBits = value(sbn::triggerSource::BNB);
    if (fTriggerType == "NuMI")           fTriggerBits = value(sbn::triggerSource::NuMI);
    if (fTriggerType == "OffbeamBNB")     fTriggerBits = value(sbn::triggerSource::OffbeamBNB);
    if (fTriggerType == "OffbeamNuMI")    fTriggerBits = value(sbn::triggerSource::OffbeamNuMI);
    if (fTriggerType == "Unknown")        fTriggerBits = value(sbn::triggerSource::Unknown);

    return;
}

bool filter::TriggerTypeFilter::filter(art::Event &event)
{
    bool filterPass = false;

    auto const& triggerVec = event.getProduct<std::vector<raw::Trigger>>(fTriggerDataLabel);

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
