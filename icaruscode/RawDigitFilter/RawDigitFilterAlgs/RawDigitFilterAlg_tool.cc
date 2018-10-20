#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "icaruscode/RawDigitFilter/RawDigitFilterAlgs/IRawDigitFilter.h"
#include "icaruscode/RawDigitFilter/RawDigitFilterAlgs/RawDigitCharacterizationAlg.h"
#include "icaruscode/RawDigitFilter/RawDigitFilterAlgs/RawDigitBinAverageAlg.h"

#include <cmath>
#include <algorithm>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace caldata
{
    
class RawDigitFilterAlg : virtual public IRawDigitFilter
{
public:
    
    // Copnstructors, destructor.
    explicit RawDigitFilterAlg(const fhicl::ParameterSet& pset);
    ~RawDigitFilterAlg();
    
    // Overrides.
    void   configure(const fhicl::ParameterSet& pset)                      override;
    void   initializeHistograms(art::TFileDirectory&)                const override;
    size_t plane()                                                   const override {return fPlane;}
    
    void FilterWaveform(RawDigitVector&, size_t, size_t, float = 0.) const override;
    
private:
    
    void doAdaptiveFilter(RawDigitVector&) const;

    // Fcl parameters.
    float                              fTruncMeanFraction;     ///< Fraction for truncated mean
    size_t                             fStructuringElement;    ///< Structuring element to use with Top Hat filter
    bool                               fFillHistograms;        ///< if true then will fill diagnostic hists
    
    size_t                             fPlane;
    
    // We'll use this algorithm internally here too
    RawDigitCharacterizationAlg        fCharacterizationAlg;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    art::ServiceHandle<geo::Geometry>  fGeometry;             ///< pointer to Geometry service
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();   ///< Detector properties service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitFilterAlg::RawDigitFilterAlg(fhicl::ParameterSet const & pset) :
    fPlane(0),
    fCharacterizationAlg(pset)
{
    configure(pset);

    // Report.
    mf::LogInfo("RawDigitFilterAlg") << "RawDigitFilterAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitFilterAlg::~RawDigitFilterAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitFilterAlg::configure(fhicl::ParameterSet const & pset)
{
    fTruncMeanFraction     = pset.get<float> ("TruncMeanFraction",  0.15);
    fStructuringElement    = pset.get<size_t>("StructuringElement",   30);
    fFillHistograms        = pset.get<bool>  ("FillHistograms",    false);
    
    fCharacterizationAlg.reconfigure(pset);
}

void RawDigitFilterAlg::initializeHistograms(art::TFileDirectory& tfs) const
{
    
    return;
}
    
void RawDigitFilterAlg::FilterWaveform(RawDigitVector& dataVec, size_t channel, size_t cnt, float pedestal) const
{
    // We need to start by collecting mean/rms
    float truncMean;
    float rmsVal;
    int   numBins;
    
    fCharacterizationAlg.getMeanAndRms(dataVec, truncMean, rmsVal, numBins);
    
    // Declare erosion vector
    std::vector<float> meanVec;
    
    meanVec.resize(dataVec.size(), 0);
    std::vector<float>::iterator meanVecItr = meanVec.begin();
    
    size_t halfStructure(fStructuringElement/2);
    
    float runningSum(0.);
    
    for(size_t idx = 0; idx < halfStructure - 1; idx++) runningSum += dataVec[idx];
    
    // First pass through to build the erosion vector
    for(RawDigitVector::iterator dataItr = dataVec.begin(); dataItr != dataVec.end(); dataItr++)
    {
        size_t startOffset = std::distance(dataVec.begin(),dataItr);
        size_t stopOffset  = std::distance(dataItr,dataVec.end());
        size_t count       = std::min(2 * halfStructure, std::min(startOffset + halfStructure, halfStructure + stopOffset));
        
        if (startOffset >  halfStructure) runningSum -= *(dataItr - halfStructure - 1);
        if (stopOffset  >= halfStructure) runningSum += *(dataItr + halfStructure - 1);
        
        float adaptVal = runningSum / float(count) - truncMean;
        
        *meanVecItr++ = adaptVal;
    }
    
    // Now do baseline correction
    for(size_t idx = 0; idx < dataVec.size(); idx++)
        dataVec[idx] = dataVec[idx] - std::round(meanVec[idx]);

    return;
}

DEFINE_ART_CLASS_TOOL(RawDigitFilterAlg)

}
