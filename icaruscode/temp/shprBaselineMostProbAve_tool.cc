////////////////////////////////////////////////////////////////////////
/// \file   Baseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IBaseline.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "icaruscode/TPC/Utilities/SignalShaperServiceICARUS.h"
#include "icarus_signal_processing/WaveformTools.h"

#include <fstream>
#include <algorithm> // std::minmax_element()

namespace icarus_tool
{

class shprBaselineMostProbAve : IBaseline
{
public:
    explicit shprBaselineMostProbAve(const fhicl::ParameterSet& pset);
    
    ~shprBaselineMostProbAve();
    
    void configure(const fhicl::ParameterSet& pset)                                        override;
    void outputHistograms(art::TFileDirectory&)                                      const override;
    
    float GetBaseline(const icarusutil::TimeVec&, raw::ChannelID_t, size_t, size_t)  const override;
    
private:
    std::pair<float,int> GetBaseline(const icarusutil::TimeVec&, int, size_t, size_t) const;
    
    size_t fMaxROILength;    ///< Maximum length for calculating Most Probable Value

    art::ServiceHandle<util::SignalShaperServiceICARUS> fSignalShaper;
    std::unique_ptr<icarus_tool::IWaveformTool>          fWaveformTool;
};
    
//----------------------------------------------------------------------
// Constructor.
shprBaselineMostProbAve::shprBaselineMostProbAve(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
shprBaselineMostProbAve::~shprBaselineMostProbAve()
{
}
    
void shprBaselineMostProbAve::configure(const fhicl::ParameterSet& pset)
{
    // Recover our fhicl variable
    fMaxROILength = pset.get<size_t>("MaxROILength", 100);
    
    // Get signal shaping service.
    fSignalShaper = art::ServiceHandle<util::SignalShaperServiceICARUS>();

    // Let's apply some smoothing as an experiment... first let's get the tool we need
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    fWaveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);

    return;
}

    
float shprBaselineMostProbAve::GetBaseline(const icarusutil::TimeVec& holder,
                                       raw::ChannelID_t               channel,
                                       size_t                         roiStart,
                                       size_t                         roiLen) const
{
    float base(0.);

    if (roiLen > 1)
    {
        // Recover the expected electronics noise on this channel
        float  deconNoise = 1.26491 * fSignalShaper->GetDeconNoise(channel);
        int    binRange   = std::max(1, int(std::round(deconNoise)));
        size_t halfLen    = std::min(fMaxROILength,roiLen/2);
        size_t roiStop    = roiStart + roiLen;
        
        std::pair<float,int> baseFront = GetBaseline(holder, binRange, roiStart,          roiStart + halfLen);
        std::pair<float,int> baseBack  = GetBaseline(holder, binRange, roiStop - halfLen, roiStop           );
        
        if (std::fabs(baseFront.first - baseBack.first) > 3. * deconNoise)
        {
            if      (baseFront.second > 3 * baseBack.second  / 2) base = baseFront.first;
            else if (baseBack.second  > 3 * baseFront.second / 2) base = baseBack.first;
            else                                                  base = std::max(baseFront.first,baseBack.first);
        }
        else
            base = (baseFront.first*baseFront.second + baseBack.first*baseBack.second)/float(baseFront.second+baseBack.second);
    }
    
    return base;
}
    
std::pair<float,int> shprBaselineMostProbAve::GetBaseline(const icarusutil::TimeVec& holder,
                                                      int                            binRange,
                                                      size_t                         roiStart,
                                                      size_t                         roiStop) const
{
    std::pair<float,int> base(0.,1);
    
    if (roiStop > roiStart)
    {
        // Get the truncated mean and rms
        std::vector<float> temp(roiStop - roiStart + 1,0.);
        
        std::copy(holder.begin() + roiStart,holder.begin() + roiStop,temp.begin());
        
        fWaveformTool->getTruncatedMean(temp, base.first, base.second);
    }
    
    return base;
}
    
void shprBaselineMostProbAve::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "BaselinePlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fBaselineVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "BaselinePlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "Baseline;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fBaselineVec.at(bin).Re());
    }
*/
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(shprBaselineMostProbAve)
}
