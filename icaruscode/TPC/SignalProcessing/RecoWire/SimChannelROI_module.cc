////////////////////////////////////////////////////////////////////////
//
// SimChannelROI class - An ROI finding module for complete deconvolved waveforms
//
// usher@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>
#include <fstream>
#include <random>

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "canvas/Utilities/Exception.h"

#include "cetlib_except/coded_exception.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataobj/Utilities/sparse_vector.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/ICARUS/TPC/ChannelROI.h"
#include "icaruscode/TPC/Utilities/ChannelROICreator.h"

#include "icaruscode/TPC/SignalProcessing/RecoWire/ROITools/IROILocator.h"

#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Denoising.h"

///creation of calibrated signals on wires
namespace caldata {

class SimChannelROI : public art::EDProducer
{
public:
// create calibrated signals on wires. this class runs 
// an fft to remove the electronics shaping.     
    explicit SimChannelROI(fhicl::ParameterSet const& pset);
    virtual ~SimChannelROI();
    void     produce(art::Event& evt); 
    void     beginJob(); 
    void     endJob();                 
    void     reconfigure(fhicl::ParameterSet const& p);
private:
    using RegionsOfInterest_t  = lar::sparse_vector<float>;
    using ROIVec               = std::vector<RegionsOfInterest_t>;

    std::vector<art::InputTag>   fSimChannelLabelVec;         ///< vector of modules that made digits
    std::vector<std::string>     fOutInstanceLabelVec;        ///< The output instance labels to apply
    bool                         fDiagnosticOutput;           ///< secret diagnostics flag
    float                        fGain;                       ///< Provides conversion from # electrons to ADCs
    size_t                       fEventCount;                 ///< count of event processed

    unsigned int                 fNTimeSamples;               ///< number of ADC readout samples in all readout frames (per event)

    const geo::GeometryCore*     fGeometry = lar::providerFrom<geo::Geometry>();
    
}; // class SimChannelROI

DEFINE_ART_MODULE(SimChannelROI)

//-------------------------------------------------
SimChannelROI::SimChannelROI(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
    this->reconfigure(pset);

    for(const auto& wireLabel : fOutInstanceLabelVec)
        produces< std::vector<recob::ChannelROI>>(wireLabel);
}

//-------------------------------------------------
SimChannelROI::~SimChannelROI()
{
}

//////////////////////////////////////////////////////
void SimChannelROI::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the parameters
    fSimChannelLabelVec    = pset.get<std::vector<art::InputTag>>("SimChannelLabelVec",   std::vector<art::InputTag>()={"decon1droi"});
    fOutInstanceLabelVec   = pset.get<std::vector<std::string>  >("OutInstanceLabelVec",                            {"PHYSCRATEDATA"});
    fDiagnosticOutput      = pset.get< bool                     >("DaignosticOutput",                                           false);
    // Note that the gain for the 2D drift simulation is 80.1 and this is default
    // For the 1D drift simulation we have been using 68.7
    fGain                  = pset.get< float                    >("ElectronicsGain",                                             80.1);

    //detector properties information
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    fNTimeSamples = detProp.NumberTimeSamples();
    
    return;
}

//-------------------------------------------------
void SimChannelROI::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void SimChannelROI::endJob()
{
}

//////////////////////////////////////////////////////
void SimChannelROI::produce(art::Event& evt)
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // We need to loop through the list of Wire data we have been given
    for(size_t labelIdx = 0; labelIdx < fSimChannelLabelVec.size(); labelIdx++)
    {
        const art::InputTag& driftEModuleLabel = fSimChannelLabelVec[labelIdx];

        // make a collection of ChannelROI objects
        std::unique_ptr<std::vector<recob::ChannelROI>> channelROICol(new std::vector<recob::ChannelROI>);

        mf::LogInfo("SimChannelROI") << "SimChannelROI, looking for decon1droi data at " << driftEModuleLabel << std::endl;
    
        // Read in the collection of full length deconvolved waveforms
        // Note we assume this list is sorted in increasing channel number!
        std::vector<const sim::SimChannel*> simChannelHandle;
        evt.getView(driftEModuleLabel,simChannelHandle);

        mf::LogInfo("SimChannelROI") << "--> Recovered SimChannel data, size: " << simChannelHandle.size() << std::endl;
    
        if (!simChannelHandle.size())
        {
            evt.put(std::move(channelROICol), fOutInstanceLabelVec[labelIdx]);
            
            continue;
        }
   
        // Reserve the room for the output
        channelROICol->reserve(simChannelHandle.size());
    
        // we simply loop over simChannels (arranged by channel)
        size_t               numChannels(0);
    
        for(const auto& simChannel : simChannelHandle)
        {
            raw::ChannelID_t channel = simChannel->Channel();

            // Total charge deposition from SimChannel
            std::vector<short int> chargeVec;

            // vector that will be moved into the Wire object
            recob::ChannelROI::RegionsOfInterest_t intROIVec;

            size_t startTick(0);
            size_t gapTicks(0);

            // Here go through the input simchannel and build out the charge array
            for(size_t tick = 0; tick < fNTimeSamples; tick++)
            {
                // Note that we attempt to explicitly remove the time offsets between the planes and leave it to the response
                // functions to handle this. So we reference the time to the first plane
                //int tdc = clockData.TPCTick2TDC(tick + detProp.GetXTicksOffset(planeID) - detProp.GetXTicksOffset(geo::PlaneID(cryostat,tpc,0)));
                int tdc = clockData.TPCTick2TDC(tick);

                // continue if tdc < 0
                if( tdc < 0 ) continue;

                float charge = simChannel->Charge(tdc) / fGain;  // Charge returned in number of electrons

                // Need to insure we don't exceed short int limits
                if (charge > std::numeric_limits<short>::max()) charge = std::numeric_limits<short>::max();

                if (charge > 0.)
                {
                    if      (chargeVec.empty()) startTick = tick;
                    else if (gapTicks > 0)      chargeVec.resize(chargeVec.size()+gapTicks,0);

                    chargeVec.push_back(short(std::round(charge)));

                    gapTicks = 0;
                }
                else
                {
                    if (gapTicks > 2)
                    {
                        if (!chargeVec.empty())
                        {
                            intROIVec.add_range(startTick,std::move(chargeVec));

                            chargeVec.clear();
                        }

                        gapTicks = 0;
                    }
                    else gapTicks++;
                }
            } // loop over tdcs

            // Look for end case
            if (!chargeVec.empty()) intROIVec.add_range(startTick,std::move(chargeVec));

            if (!intROIVec.empty()) channelROICol->push_back(recob::ChannelROICreator(std::move(intROIVec),channel).move());

            numChannels++;
        }
        
        // Time to stroe everything
        if(channelROICol->size() == 0) mf::LogWarning("SimChannelROI") << "No wires made for this event.";

        evt.put(std::move(channelROICol), fOutInstanceLabelVec[labelIdx]);
    }

    fEventCount++;

    return;
} // produce

} // end namespace caldata
