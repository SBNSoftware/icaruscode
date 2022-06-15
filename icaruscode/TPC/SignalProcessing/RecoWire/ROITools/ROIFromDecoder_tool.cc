////////////////////////////////////////////////////////////////////////
/// \file   ROIFromDecoder_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/ROITools/IROILocator.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataobj/RecoBase/Wire.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <fstream>

namespace icarus_tool
{

class ROIFromDecoder : public IROILocator
{
public:
    explicit ROIFromDecoder(const fhicl::ParameterSet& pset);
    
    ~ROIFromDecoder();
    
    void configure(const fhicl::ParameterSet& pset) override;
    void initializeHistograms(art::TFileDirectory&) override {return;}
    
    void FindROIs(const art::Event&, const ArrayFloat&, const geo::PlaneID&, ArrayFloat&, ArrayBool&) override;
    
private:
    // A magic map because all tools need them
    using TPCIDToLabelMap = std::map<geo::TPCID,std::string>;

    TPCIDToLabelMap            fTPCIDToLabelMap;      ///< Translate from TPCID to a substring

    // fhicl parameters
    std::vector<art::InputTag> fROILabelVec;          ///< List of input files to search

    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
ROIFromDecoder::ROIFromDecoder(const fhicl::ParameterSet& pset)
{
    std::cout << "*** In ROIFromDecoder constructor ***" << std::endl;
    configure(pset);
}
    
ROIFromDecoder::~ROIFromDecoder()
{
}
    
void ROIFromDecoder::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fROILabelVec = pset.get<std::vector<art::InputTag>>("ROILabelVec", {""});

    // Build our map
    fTPCIDToLabelMap[geo::TPCID(0,0)] = "EE";
    fTPCIDToLabelMap[geo::TPCID(0,1)] = "EE";
    fTPCIDToLabelMap[geo::TPCID(0,2)] = "EW";
    fTPCIDToLabelMap[geo::TPCID(0,3)] = "EW";
    fTPCIDToLabelMap[geo::TPCID(1,0)] = "WE";
    fTPCIDToLabelMap[geo::TPCID(1,1)] = "WE";
    fTPCIDToLabelMap[geo::TPCID(1,2)] = "WW";
    fTPCIDToLabelMap[geo::TPCID(1,3)] = "WW";

    return;
}

void ROIFromDecoder::FindROIs(const art::Event& event, const ArrayFloat& inputImage, const geo::PlaneID& planeID, ArrayFloat& output, ArrayBool& outputROIs)
{
    // First thing is find the correct data product to recover ROIs
    TPCIDToLabelMap::const_iterator tpcItr = fTPCIDToLabelMap.find(planeID.asTPCID());

    if (tpcItr == fTPCIDToLabelMap.end())
    {
        std::cout << "Can't happen? for planeID: " << planeID << std::endl;
        return;
    }

    art::InputTag roiLabel("");

    for(const auto& tag : fROILabelVec)
    {
        if (tag.instance().find(tpcItr->second) != std::string::npos)
        {
            roiLabel = tag;
            break;
        }
    }

    if (roiLabel != "")
    {
        // do stuff here
    
        // Read in the collection of full length deconvolved waveforms
        // Note we assume this list is sorted in increasing channel number!
        art::Handle< std::vector<recob::Wire>> wireVecHandle;
        
        event.getByLabel(roiLabel, wireVecHandle);

        // The input data product will contain information for every plane in a physical TPC... unfortunately, we only want a single plane within
        // a logical TPC... so we need to loop through to find what we want
        for(const auto& wireData : *wireVecHandle)
        {
            std::vector<geo::WireID> wireIDVec = fGeometry->ChannelToWire(wireData.Channel());

            for(const auto& wireID : wireIDVec)
            {
                if (wireID.asPlaneID() != planeID) continue;

                if (wireID.Wire > outputROIs.size())
                {
                    std::cout << "#################################### Wire out of bounds! Wire: " << wireID.Wire << ", max: " << outputROIs.size() << " #################" << std::endl;
                    continue;
                }

                // Recover this wire's output vector
                VectorBool& channelData = outputROIs[wireID.Wire];

                // What we need to do is find the ROIs in the input wire data and then translate to the output
                const recob::Wire::RegionsOfInterest_t signalROIs = wireData.SignalROI();

                for(const auto& range : signalROIs.get_ranges())
                {
                    size_t startTick = range.begin_index();
                    size_t roiLen    = range.data().size();
                    size_t stopTick  = startTick + roiLen;

                    if (startTick > channelData.size())
                    {
                        std::cout << "*** ROI decoder has start tick larger than output array, start: " << startTick << ", array size: " << channelData.size() << std::endl;
                        continue;
                    }

                    if (stopTick > channelData.size())
                    {
                        std::cout << "*** ROI decoder has ROI length larger than output array, start: " << startTick << ", end: " << stopTick << ", array size: " << channelData.size() << std::endl;
                        continue;
                    }

                    std::fill(channelData.begin() + startTick, channelData.begin() + stopTick, true);
                }
            }
        }

    }

    return;
}

DEFINE_ART_CLASS_TOOL(ROIFromDecoder)
}
