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
#include <TTree.h>
#include <TFile.h>

#include <fstream>

namespace icarus_tool
{

class ROIMorphological2D : public IROILocator
{
public:
    explicit ROIMorphological2D(const fhicl::ParameterSet& pset);
    
    ~ROIMorphological2D();
    
    void configure(const fhicl::ParameterSet& pset) override;
    void initializeHistograms(art::TFileDirectory&) override;
    
    void FindROIs(const art::Event&, const ArrayFloat&, const geo::PlaneID&, ArrayFloat&, ArrayBool&) override;
    
private:
    // This is for the baseline...
    float getMedian(const icarus_signal_processing::VectorFloat, const unsigned int) const;

    bool                 fOutputHistograms;           ///< Diagnostic histogram output

    // fhicl parameters
    std::vector<size_t>  fStructuringElement;         ///< Structuring element for morphological filter
    std::vector<float>   fThreshold;                  ///< Threshold to apply for saving signal

    // tuple output if requested
    std::vector<float>   fMedianVec;
    std::vector<float>   fRMSVec;
    std::vector<float>   fMinValVec;
    std::vector<float>   fMaxValVec;
    std::vector<float>   fRangeVec;
    std::vector<bool>    fHasROIVec;

    TTree*               fTupleTree;        ///< output analysis tree
};
    
//----------------------------------------------------------------------
// Constructor.
ROIMorphological2D::ROIMorphological2D(const fhicl::ParameterSet& pset) : fOutputHistograms(false)
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

void ROIMorphological2D::FindROIs(const art::Event& event, const ArrayFloat& inputImage, const geo::PlaneID& planeID, ArrayFloat& morphedWaveforms, ArrayBool& outputROIs)
{
    if (morphedWaveforms.size() != inputImage.size()) morphedWaveforms.resize(inputImage.size(),icarus_signal_processing::VectorFloat(inputImage[0].size()));

    for(auto& morph : morphedWaveforms) std::fill(morph.begin(),morph.end(),0.);  // explicit initialization

    // Use this to get the 2D Dilation of each waveform
    icarus_signal_processing::Dilation2D(fStructuringElement[0],fStructuringElement[1])(inputImage.begin(),inputImage.size(),morphedWaveforms.begin());

    if (fOutputHistograms)
    {
        fMedianVec.clear();
        fRMSVec.clear();
        fMinValVec.clear();
        fMaxValVec.clear();
        fRangeVec.clear();
        fHasROIVec.clear();
    }

    // Now traverse each waveform and look for the ROIs
    for(size_t waveIdx = 0; waveIdx < morphedWaveforms.size(); waveIdx++)
    {
        // We start working with the morphed waveform
        VectorFloat& morphedWave = morphedWaveforms[waveIdx];

        // We need to zero suppress so we can find the rms
        float median = getMedian(morphedWave, morphedWave.size());

        for(auto& val : morphedWave) val -= median;

//        float threshold = rms * fThreshold[planeID.Plane];
        float threshold = fThreshold[planeID.Plane];

        // Right size the selected values array
        VectorBool& selVals = outputROIs[waveIdx];

        if (selVals.size() != morphedWave.size()) selVals.resize(morphedWave.size());

        std::fill(selVals.begin(),selVals.end(),false);

        bool hasROI(false);

        for(size_t idx = 0; idx < morphedWave.size(); idx++)
        {
            if (morphedWave[idx] > threshold) 
            {
                selVals[idx] = true;
                hasROI       = true;
            }
        }

        if (fOutputHistograms)
        {
            VectorFloat rmsVec = morphedWave;
            size_t      maxIdx = 0.75 * rmsVec.size();

            std::nth_element(rmsVec.begin(), rmsVec.begin() + maxIdx, rmsVec.end());

            float rms    = std::sqrt(std::inner_product(rmsVec.begin(), rmsVec.begin() + maxIdx, rmsVec.begin(), 0.) / float(maxIdx));
            float minVal = *std::min_element(morphedWave.begin(),morphedWave.end());
            float maxVal = *std::max_element(morphedWave.begin(),morphedWave.end());
            
            fMedianVec.emplace_back(median);
            fRMSVec.emplace_back(rms);
            fMinValVec.emplace_back(minVal);
            fMaxValVec.emplace_back(maxVal);
            fRangeVec.emplace_back(maxVal-minVal);
            fHasROIVec.emplace_back(hasROI);
        }
    }

    if (fOutputHistograms) fTupleTree->Fill();
     
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
    
void ROIMorphological2D::initializeHistograms(art::TFileDirectory& histDir)
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.

    fTupleTree = histDir.make<TTree>("ROIFinder", "Tree by ROIMorphological2D_tool");

    fTupleTree->Branch("medians",  "std::vector<float>", &fMedianVec);
    fTupleTree->Branch("RMS",      "std::vector<float>", &fRMSVec);
    fTupleTree->Branch("minvals",  "std::vector<float>", &fMinValVec);
    fTupleTree->Branch("maxvals",  "std::vector<float>", &fMaxValVec);
    fTupleTree->Branch("range",    "std::vector<float>", &fRangeVec);
    fTupleTree->Branch("hasROI",   "std::vector<bool>",  &fHasROIVec);

    fOutputHistograms = true;
    
    return;
}

DEFINE_ART_CLASS_TOOL(ROIMorphological2D)
}
