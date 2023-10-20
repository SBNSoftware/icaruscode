////////////////////////////////////////////////////////////////////////
//
// HitConvert class - Converts icarus::Hit objects to recob::Hit
//
// usher@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/CoreUtils/zip.h"
#include "lardataobj/RecoBase/Hit.h"

#include "icaruscode/IcarusObj/Hit.h"

///creation of calibrated signals on wires
namespace caldata {

class HitConvert : public art::EDProducer
{
public:
// create calibrated signals on wires. this class runs 
// an fft to remove the electronics shaping.     
    explicit HitConvert(fhicl::ParameterSet const& pset);
    void     produce(art::Event& evt); 
    void     beginJob(); 
    void     endJob();                 
    void     reconfigure(fhicl::ParameterSet const& p);
private:

    std::vector<art::InputTag>                                 fHitModuleLabelVec;         ///< vector of modules that made digits
    std::vector<std::string>                                   fOutInstanceLabelVec;        ///< The output instance labels to apply
    bool                                                       fDiagnosticOutput;           ///< secret diagnostics flag
    size_t                                                     fEventCount;                 ///< count of event processed

    const geo::GeometryCore*                                   fGeometry = lar::providerFrom<geo::Geometry>();
    
}; // class HitConvert

DEFINE_ART_MODULE(HitConvert)

//-------------------------------------------------
HitConvert::HitConvert(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
    this->reconfigure(pset);

    for(const auto& wireLabel : fOutInstanceLabelVec)
    {
        produces<std::vector<recob::Hit>>(wireLabel);
    }
}

//////////////////////////////////////////////////////
void HitConvert::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the parameters
    fHitModuleLabelVec     = pset.get<std::vector<art::InputTag>>("HitModuleLabelVec",    std::vector<art::InputTag>()={"converHits"});
    fOutInstanceLabelVec   = pset.get<std::vector<std::string>>  ("OutInstanceLabelVec",                                         {""});
    fDiagnosticOutput      = pset.get< bool                     >("DaignosticOutput",                                           false);

    return;
}

//-------------------------------------------------
void HitConvert::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void HitConvert::endJob()
{
}

//////////////////////////////////////////////////////
void HitConvert::produce(art::Event& evt)
{
    // We need to loop through the list of Wire data we have been given
    // This construct from Gianluca Petrillo who invented it and should be given all credit for it! 
    for(auto const& [hitLabel, instanceName] : util::zip(fHitModuleLabelVec, fOutInstanceLabelVec))
    {
        // make a collection of Wires
        std::unique_ptr<std::vector<recob::Hit>> outputHitCol = std::make_unique<std::vector<recob::Hit>>();

        mf::LogInfo("HitConvert") << "HitConvert, looking for icarus::Hit data at " << hitLabel.encode();
    
        // Read in the collection of full length deconvolved waveforms
       const std::vector<icarus::Hit>& inputHitVec = evt.getProduct<std::vector<icarus::Hit>>(hitLabel);

        mf::LogInfo("HitConvert") << "--> Recovered icarus::Hit data, size: " << inputHitVec.size();
    
        if (!inputHitVec.empty())
        {
            // Reserve the room for the output
            outputHitCol->reserve(inputHitVec.size());

            // Loop through the input ChannelROI collection
            for(const auto& hit : inputHitVec)
            {
                // Recover the channel and the view
                raw::ChannelID_t channel = hit.Channel();
                geo::View_t      view    = fGeometry->View(channel);

                recob::Hit recobHit(channel,
                                    hit.StartTick(),
                                    hit.EndTick(),
                                    hit.PeakTime(),
                                    hit.SigmaPeakTime(),
                                    hit.RMS(),
                                    hit.PeakAmplitude(),
                                    hit.SigmaPeakAmplitude(),
                                    hit.SummedADC(),
                                    hit.Integral(),
                                    hit.SigmaIntegral(),
                                    hit.Multiplicity(),
                                    hit.LocalIndex(),
                                    hit.GoodnessOfFit(),
                                    hit.DegreesOfFreedom(),
                                    view,
                                    hit.SignalType(),
                                    hit.WireID()
                                    );

                outputHitCol->push_back(recobHit);
            }

            mf::LogInfo("HitConvert") << "--> Outputting hit data, size: " << outputHitCol->size() << " with instance name: " << instanceName;

            // Time to stroe everything
            if(outputHitCol->empty()) mf::LogWarning("HitConvert") << "No hits made for this event.";
        }

        evt.put(std::move(outputHitCol), instanceName);
    }

    fEventCount++;

    return;
} // produce

} // end namespace caldata
