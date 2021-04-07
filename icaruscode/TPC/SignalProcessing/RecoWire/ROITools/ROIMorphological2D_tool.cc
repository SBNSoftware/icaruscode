////////////////////////////////////////////////////////////////////////
/// \file   ROIMorphological2D.cc
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
#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Denoising.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <fstream>

namespace icarus_tool
{

class ROIMorphological2D : public IROILocator
{
public:
    explicit ROIMorphological2D(const fhicl::ParameterSet& pset);
    
    ~ROIMorphological2D();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void FindROIs(const ArrayFloat&, const geo::PlaneID&, ArrayBool&)    const override;
    
private:
    // This is for the baseline...
    float getMedian(const icarus_signal_processing::VectorFloat, const unsigned int) const;

    // fhicl parameters
    std::vector<size_t>  fStructuringElement;         ///< Structuring element for morphological filter
    std::vector<float>   fThreshold;                  ///< Threshold to apply for saving signal
};
    
//----------------------------------------------------------------------
// Constructor.
ROIMorphological2D::ROIMorphological2D(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROIMorphological2D::~ROIMorphological2D()
{
}
    
void ROIMorphological2D::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fStructuringElement = pset.get<std::vector<size_t> >("StructuringElement", std::vector<size_t>()={8,16});
    fThreshold          = pset.get<std::vector<float>  >("Threshold",          std::vector<float>()={2.75,2.75,2.75});

    return;
}

void ROIMorphological2D::FindROIs(const ArrayFloat& inputImage, const geo::PlaneID& planeID, ArrayBool& outputROIs) const
{
    icarus_signal_processing::ArrayFloat morphedWaveforms(inputImage.size());

    // Use this to get the 2D Dilation of each waveform
    icarus_signal_processing::Dilation2D(fStructuringElement[0],fStructuringElement[1])(inputImage.begin(),inputImage.size(),morphedWaveforms.begin());

    // Now traverse each waveform and look for the ROIs
    for(size_t waveIdx = 0; waveIdx < morphedWaveforms.size(); waveIdx++)
    {
        // We start working with the morphed waveform
        const VectorFloat& morphedWave = morphedWaveforms[waveIdx];

        // We need to zero suppress so we can find the rms
        float median = getMedian(morphedWave, morphedWave.size());

        VectorFloat baseVec(morphedWave.size());

        for(size_t idx = 0; idx < morphedWave.size(); idx++) baseVec[idx] = morphedWave[idx] - median;

        VectorFloat rmsVec = baseVec;
        size_t      maxIdx = 0.75 * rmsVec.size();

        std::nth_element(rmsVec.begin(), rmsVec.begin() + maxIdx, rmsVec.end());

        float rms       = std::sqrt(std::inner_product(rmsVec.begin(), rmsVec.begin() + maxIdx, rmsVec.begin(), 0.) / float(maxIdx));
        float threshold = rms * fThreshold[planeID.Plane];

//        std::cout << "==> median: " << median << ", rms: " << rms << ", threshold: " << threshold << std::endl;

        // Right size the selected values array
        VectorBool& selVals = outputROIs[waveIdx];

        selVals.resize(morphedWave.size(),false);

        for(size_t idx = 0; idx < baseVec.size(); idx++)
        {
//            if (std::abs(baseVec[idx]) > threshold) selVals[idx] = true;
            if (baseVec[idx] > threshold) selVals[idx] = true;
        }

        // Check to see if we should save the baseVec
//        if (fOutputMorphed)
//        {
//            raw::ChannelID_t channel = planeIDToDataPair.first[waveIdx];
//            geo::View_t      view    = fGeometry->View(channel);
//
//            morphedVec.push_back(recob::WireCreator(std::move(baseVec),channel,view).move());
//        }
    }
     
    return;
}

float ROIMorphological2D::getMedian(icarus_signal_processing::VectorFloat vals, const unsigned int nVals) const
{
    float median(0.);

    if (nVals > 2) 
    {
        if (nVals % 2 == 0) 
        {
            const auto m1 = vals.begin() + nVals / 2 - 1;
            const auto m2 = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m1, vals.begin() + nVals);
            const auto e1 = *m1;
            std::nth_element(vals.begin(), m2, vals.begin() + nVals);
            const auto e2 = *m2;
            median = (e1 + e2) / 2.0;
        } 
        else 
        {
            const auto m = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m, vals.begin() + nVals);
            median = *m;
        }
    }

    return median;
}

DEFINE_ART_CLASS_TOOL(ROIMorphological2D)
}
