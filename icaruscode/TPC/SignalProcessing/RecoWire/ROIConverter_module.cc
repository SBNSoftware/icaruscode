////////////////////////////////////////////////////////////////////////
//
// ROIConvert class - An ROI finding module for complete deconvolved waveforms
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
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "icaruscode/IcarusObj/ChannelROI.h"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_arena.h"
#include "tbb/spin_mutex.h"
#include "tbb/concurrent_hash_map.h"

///creation of calibrated signals on wires
namespace caldata {

tbb::spin_mutex ROIConvertSpinMutex;

class ROIConvert : public art::EDProducer
{
public:
// create calibrated signals on wires. this class runs 
// an fft to remove the electronics shaping.     
    explicit ROIConvert(fhicl::ParameterSet const& pset);
    virtual ~ROIConvert();
    void     produce(art::Event& evt); 
    void     beginJob(); 
    void     endJob();                 
    void     reconfigure(fhicl::ParameterSet const& p);
private:

    std::vector<art::InputTag>                                 fWireModuleLabelVec;         ///< vector of modules that made digits
    std::vector<std::string>                                   fOutInstanceLabelVec;        ///< The output instance labels to apply
    bool                                                       fDiagnosticOutput;           ///< secret diagnostics flag
    size_t                                                     fEventCount;                 ///< count of event processed

    const geo::GeometryCore*                                   fGeometry = lar::providerFrom<geo::Geometry>();
    
}; // class ROIConvert

DEFINE_ART_MODULE(ROIConvert)

//-------------------------------------------------
ROIConvert::ROIConvert(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
    this->reconfigure(pset);

    for(const auto& wireLabel : fOutInstanceLabelVec)
    {
        produces<std::vector<recob::Wire>>(wireLabel);
    }
}

//-------------------------------------------------
ROIConvert::~ROIConvert()
{
}

//////////////////////////////////////////////////////
void ROIConvert::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the parameters
    fWireModuleLabelVec    = pset.get<std::vector<art::InputTag>>("WireModuleLabelVec",   std::vector<art::InputTag>()={"decon1droi"});
    fOutInstanceLabelVec   = pset.get<std::vector<std::string>>  ("OutInstanceLabelVec",                            {"PHYSCRATEDATA"});
    fDiagnosticOutput      = pset.get< bool                     >("DaignosticOutput",                                           false);
    
    return;
}

//-------------------------------------------------
void ROIConvert::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void ROIConvert::endJob()
{
}

//////////////////////////////////////////////////////
void ROIConvert::produce(art::Event& evt)
{
    // We need to loop through the list of Wire data we have been given
    for(size_t labelIdx = 0; labelIdx < fWireModuleLabelVec.size(); labelIdx++)
    {
        const art::InputTag& channelLabel = fWireModuleLabelVec[labelIdx];

        // make a collection of Wires
        std::unique_ptr<std::vector<recob::Wire>> wireCol(new std::vector<recob::Wire>);

        mf::LogInfo("ROIConvert") << "ROIConvert, looking for ChannelROI data at " << channelLabel << std::endl;
    
        // Read in the collection of full length deconvolved waveforms
        // Note we assume this list is sorted in increasing channel number!
        art::Handle< std::vector<recob::ChannelROI>> channelVecHandle;
        
        evt.getByLabel(channelLabel, channelVecHandle);

        mf::LogInfo("ROIConvert") << "--> Recovered ChannelROI data, size: " << channelVecHandle->size() << std::endl;
    
        if (!channelVecHandle->size())
        {
            evt.put(std::move(wireCol), channelLabel.instance());
            fEventCount++;
            
            return;
        }
   
        // Reserve the room for the output
        wireCol->reserve(channelVecHandle->size());

        // Loop through the input ChannelROI collection
        for(const auto& channelROI : *channelVecHandle)
        {
            // Recover the channel and the view
            raw::ChannelID_t channel = channelROI.Channel();
            geo::View_t      view    = fGeometry->View(channel);

            // Create an ROI vector for output
            recob::Wire::RegionsOfInterest_t ROIVec;
 
            // Loop through the ROIs for this channel
            const recob::ChannelROI::RegionsOfInterest_t& channelROIs = channelROI.SignalROI();

            for(const auto& range : channelROIs.get_ranges())
            {
                size_t startTick = range.begin_index();

                std::vector<float> dataVec(range.data().size());

                for(size_t binIdx = 0; binIdx < range.data().size(); binIdx++) dataVec[binIdx] = range.data()[binIdx];

                ROIVec.add_range(startTick, std::move(dataVec));
            }

            wireCol->push_back(recob::WireCreator(std::move(ROIVec),channel,view).move());
        }

        // Time to stroe everything
        if(wireCol->size() == 0) mf::LogWarning("ROIConvert") << "No wires made for this event.";

        evt.put(std::move(wireCol), fOutInstanceLabelVec[labelIdx]);
    }

    fEventCount++;

    return;
} // produce

} // end namespace caldata
