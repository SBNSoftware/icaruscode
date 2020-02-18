////////////////////////////////////////////////////////////////////////
/// \file   MorphologicalFilter.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/IRawDigitFilter.h"
#include "icaruscode/TPC/Utilities/tools/IWaveformTool.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <fstream>

namespace caldata
{

class MorphologicalFilter : public IRawDigitFilter
{
public:
    explicit MorphologicalFilter(const fhicl::ParameterSet& pset);
    
    ~MorphologicalFilter();
    
    void   configure(const fhicl::ParameterSet& pset)                      override;
    void   initializeHistograms(art::TFileDirectory&)                const override;
    size_t plane()                                                   const override {return fPlane;}
    
    void FilterWaveform(RawDigitVector&, size_t, size_t, float = 0.) const override;
    
private:
    
    using Waveform = std::vector<short>;
    
    // Average the input waveform
    void smoothInputWaveform(const RawDigitVector&, RawDigitVector&)  const;
    
    // Actual histogram initialization...
    enum HistogramType : int
    {
        CORWAVEFORM = icarus_tool::LASTELEMENT + 1
    };
    
    icarus_tool::HistogramMap initializeHistograms(size_t, size_t, size_t) const;

    // Member variables from the fhicl file
    size_t                                      fPlane;
    int                                         fNumBinsToAve;               ///< Controls the smoothing
    int                                         fStructuringElement;         ///< The window size
    bool                                        fOutputHistograms;           ///< Output histograms?

    std::vector<float>                          fAveWeightVec;               ///< Weight vector for smoothing
    float                                       fWeightSum;                  ///< sum of weights for smoothing

    art::TFileDirectory*                        fHistDirectory;
    
    // Global histograms
    TH1F*                                       fDiffMeanHist;
    TH1F*                                       fDiffRmsHist;
    TH1F*                                       fDiffMaxHist;
    TH1F*                                       fNumSigmaHist;
    TH1F*                                       fNumSigNextHist;
    TH1F*                                       fDeltaTicksHist;
    TH2F*                                       fDTixVDiffHist;
    
    std::unique_ptr<icarus_tool::IWaveformTool> fWaveformTool;

    // Services
    const geo::GeometryCore*                    fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
MorphologicalFilter::MorphologicalFilter(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
MorphologicalFilter::~MorphologicalFilter()
{
}
    
void MorphologicalFilter::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    std::vector<unsigned short> zin;
    
    fPlane                 = pset.get<size_t>                      ("Plane"                     );
    fNumBinsToAve          = pset.get< int >                       ("NumBinsToAve"              );
    fStructuringElement    = pset.get< int>                        ("StructuringElement"        );
    fOutputHistograms      = pset.get< bool  >                     ("OutputHistograms",    false);
    
    // Recover an instance of the waveform tool
    // Here we just make a parameterset to pass to it...
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    fWaveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);

    // If asked, define the global histograms
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;
        
        fHistDirectory = tfs.get();
        
        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("MF/ROIPlane_%1zu",fPlane));

        fDiffMeanHist   = dir.make<TH1F>("DiffMean", ";Diff Mean;", 100, -20.,  20.);
        fDiffRmsHist    = dir.make<TH1F>("DiffRms",  ";Diff RMS;",  100,   0.,   5.);
        fDiffMaxHist    = dir.make<TH1F>("DiffMax",  ";Diff Max;",  200,   0., 200.);
        fNumSigmaHist   = dir.make<TH1F>("NSigma",   ";#sigma;",    200,   0.,  50.);
        fNumSigNextHist = dir.make<TH1F>("NSigNext", ";#sigma;",    200,   0.,  50.);
        fDeltaTicksHist = dir.make<TH1F>("DeltaTix", ";Delta t",    200,   0., 200.);
        
        fDTixVDiffHist  = dir.make<TH2F>("DTixVDiff", ";Delta t;Max Diff", 200, 0., 200., 200, 0., 200.);
    }

    // precalculate the weight vector to use in the smoothing
    fAveWeightVec.resize(fNumBinsToAve);
    
    if (fNumBinsToAve > 1)
    {
        for(int idx = 0; idx < fNumBinsToAve/2; idx++)
        {
            float weight = idx + 1;
            
            fAveWeightVec.at(idx)                     = weight;
            fAveWeightVec.at(fNumBinsToAve - idx - 1) = weight;
        }
        
        // Watch for case of fNumBinsToAve being odd
        if (fNumBinsToAve % 2 > 0) fAveWeightVec.at(fNumBinsToAve/2) = fNumBinsToAve/2 + 1;
    }
    else fAveWeightVec.at(0) = 1.;
    
    fWeightSum = std::accumulate(fAveWeightVec.begin(),fAveWeightVec.end(), 0.);

    return;
}

void MorphologicalFilter::FilterWaveform(RawDigitVector& waveform, size_t channel, size_t cnt, float pedestal) const
{
    // The plan here is to use a morphological filtering technique to find the slowly varying baseline
    // movement and remove it

    // We make lots of vectors... erosion, dilation, average and difference
    Waveform erosionVec;
    Waveform dilationVec;
    Waveform averageVec;
    Waveform differenceVec;
    
    // Define histograms for this particular channel?
    icarus_tool::HistogramMap histogramMap = initializeHistograms(channel, cnt, waveform.size());
    
    // If histogramming, then keep track of the original input channel
    if (!histogramMap.empty()) for(size_t idx = 0; idx < waveform.size(); idx++) histogramMap.at(icarus_tool::WAVEFORM)->Fill(idx, waveform.at(idx), 1.);
    
    Waveform smoothWaveform = waveform;
    
    // If the input pedestal is non-zero then baseline correct
    if (std::abs(pedestal) > std::numeric_limits<float>::epsilon())
        std::transform(waveform.begin(),waveform.end(),smoothWaveform.begin(),[pedestal](const auto& val){return val - short(std::round(pedestal));});

    // Compute the morphological filter vectors
    fWaveformTool->getErosionDilationAverageDifference(smoothWaveform, fStructuringElement, histogramMap, erosionVec, dilationVec, averageVec, differenceVec);
    
    // What we are really interested in here is the closing vector but compute both
    Waveform openingVec;
    Waveform closingVec;
    
    fWaveformTool->getOpeningAndClosing(erosionVec,dilationVec,fStructuringElement,histogramMap,openingVec,closingVec);
    
    // Ok, get an average of the two
    std::transform(openingVec.begin(),openingVec.end(),closingVec.begin(),averageVec.begin(),[](const auto& left, const auto& right){return (left + right)/2;});
    
    // Now smooth
    std::transform(smoothWaveform.begin(),smoothWaveform.end(),averageVec.begin(),waveform.begin(),std::minus<short>());
    
    // Keep track of the corrected waveform
    if (!histogramMap.empty())
    {
        for(size_t idx = 0; idx < waveform.size(); idx++) histogramMap.at(CORWAVEFORM)->Fill(idx, waveform.at(idx), 1.);
    }

    return;
}
    
void MorphologicalFilter::smoothInputWaveform(const RawDigitVector& inputWaveform, RawDigitVector& outputWaveform) const
{
    // Vector smoothing - take the 10 bin average
    int   halfBins = fNumBinsToAve / 2;
    
    outputWaveform.resize(inputWaveform.size());
    
    // To facilitate handling the bins at the ends of the input waveform we embed in a slightly larger
    // vector which has zeroes on the ends
    RawDigitVector tempWaveform(inputWaveform.size()+fNumBinsToAve);
    
    // Set the edge bins which can't be smoothed to zero
    std::fill(tempWaveform.begin(),tempWaveform.begin()+halfBins,0.);
    std::fill(tempWaveform.end()-halfBins,tempWaveform.end(),0.);
    
    // Copy in the input waveform
    std::copy(inputWaveform.begin(),inputWaveform.end(),tempWaveform.begin()+halfBins);
    
    // Now do the smoothing
    for(size_t idx = 0; idx < inputWaveform.size(); idx++)
    {
        float weightedSum(0.);
        
        for(int wIdx = 0; wIdx < fNumBinsToAve; wIdx++) weightedSum += fAveWeightVec.at(wIdx) * tempWaveform.at(idx + wIdx);
        
        outputWaveform.at(idx) = weightedSum / fWeightSum;
    }
    
    return;
}

void MorphologicalFilter::initializeHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "ROIFinderPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fROIFinderVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "ROIFinderPlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "ROIFinder;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fROIFinderVec.at(bin).Re());
    }
*/
    
    return;
}
    
icarus_tool::HistogramMap MorphologicalFilter::initializeHistograms(size_t channel, size_t cnt, size_t waveformSize) const
{
    icarus_tool::HistogramMap histogramMap;
    
    if (fOutputHistograms)
    {
        // Try to limit to the wire number (since we are already segregated by plane)
        std::vector<geo::WireID> wids  = fGeometry->ChannelToWire(channel);
        size_t                   cryo  = wids[0].Cryostat;
        size_t                   tpc   = wids[0].TPC;
        size_t                   plane = wids[0].Plane;
        size_t                   wire  = wids[0].Wire;
        
        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("MF/ROIPlane_%1zu/c%1zu/c%1zut%1zuwire_%05zu",plane,cnt,cryo,tpc,wire));
        
        // We keep track of four histograms:
        try
        {
            //            origWaveHist   = dir.make<TProfile>(Form("Inp_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Waveform", waveform.size(),      0, waveform.size(),      -500., 500.);
            histogramMap[icarus_tool::WAVEFORM] =
                    dir.make<TProfile>(Form("MFWfm_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Waveform",   waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::EROSION] =
                    dir.make<TProfile>(Form("MFEro_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Erosion",    waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::DILATION] =
                    dir.make<TProfile>(Form("MFDil_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Dilation",   waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::AVERAGE] =
                    dir.make<TProfile>(Form("MFAve_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Average",    waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::DIFFERENCE] =
                    dir.make<TProfile>(Form("MFDif_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Difference", waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::OPENING] =
                    dir.make<TProfile>(Form("MFOpe_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Opening",    waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::CLOSING] =
                    dir.make<TProfile>(Form("MFClo_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Closing",    waveformSize, 0, waveformSize, -500., 500.);
            
            // Also, if smoothing then we would like to keep track of the original waveform too
            histogramMap[CORWAVEFORM] =
                    dir.make<TProfile>(Form("MFCor_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Corrected Waveform", waveformSize, 0, waveformSize, -500., 500.);
        } catch(...)
        {
            std::cout << "Caught exception trying to make new hists, tpc,plane,cnt,wire: " << tpc << ", " << fPlane << ", " << cnt << ", " << wire << std::endl;
        }
    }

    return histogramMap;
}

DEFINE_ART_CLASS_TOOL(MorphologicalFilter)
}
