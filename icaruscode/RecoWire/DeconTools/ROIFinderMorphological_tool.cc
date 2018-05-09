////////////////////////////////////////////////////////////////////////
/// \file   ROIFinderMorphological.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/RecoWire/DeconTools/IROIFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "icaruscode/Utilities/tools/IWaveformTool.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <fstream>

namespace icarus_tool
{

class ROIFinderMorphological : public IROIFinder
{
public:
    explicit ROIFinderMorphological(const fhicl::ParameterSet& pset);
    
    ~ROIFinderMorphological();
    
    void   configure(const fhicl::ParameterSet& pset)                      override;
    void   initializeHistograms(art::TFileDirectory&)                const override;
    size_t plane()                                                   const override {return fPlane;}
    
    void FindROIs(const Waveform&, size_t, size_t, double, CandidateROIVec&) const override;
    
private:
    // The actual ROI finding algorithm
    void findROICandidatesDiff(const Waveform&,
                               const Waveform&,
                               const Waveform&,
                               int,
                               int,
                               float,
                               CandidateROIVec&) const;
    
    // Average the input waveform
    void getErosionDilationAverageDifference(const Waveform&, const Waveform&, size_t, size_t, Waveform&, Waveform&, Waveform&, Waveform&) const;
    
    void smoothInputWaveform(const Waveform&, Waveform&)  const;

    // Member variables from the fhicl file
    size_t                                      fPlane;
    float                                       fNumSigma;                   ///< "# sigma" rms noise for ROI threshold
    int                                         fNumBinsToAve;               ///< Controls the smoothing
    size_t                                      fMax2MinDistance;            ///< Maxmimum allow peak to peak distance
    int                                         fHalfWindowSize;             ///< 1/2 the window size
    unsigned short                              fPreROIPad;                  ///< ROI padding
    unsigned short                              fPostROIPad;                 ///< ROI padding
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
    
    // This pointer changes channel-by-channel... (a kludge...)
    mutable TProfile*                           fROITempHist;
    
    std::unique_ptr<icarus_tool::IWaveformTool> fWaveformTool;

    // Services
    const geo::GeometryCore*                    fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
ROIFinderMorphological::ROIFinderMorphological(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROIFinderMorphological::~ROIFinderMorphological()
{
}
    
void ROIFinderMorphological::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    std::vector<unsigned short> zin;
    
    fPlane                 = pset.get<size_t>                      ("Plane"                     );
    fNumSigma              = pset.get< float>                      ("NumSigma"                  );
    fNumBinsToAve          = pset.get< int >                       ("NumBinsToAve"              );
    fMax2MinDistance       = pset.get< size_t>                     ("Max2MinDistance"           );
    fHalfWindowSize        = pset.get< int>                        ("HalfWindowSize"            );
    zin                    = pset.get< std::vector<unsigned short>>("roiLeadTrailPad"           );
    fOutputHistograms      = pset.get< bool  >                     ("OutputHistograms",    false);
    
    if(zin.size() != 2) {
        throw art::Exception(art::errors::Configuration)
        << "Plane ROI pad size != 2";
    }
    
    // put the ROI pad sizes into more variables
    fPreROIPad  = zin[0];
    fPostROIPad = zin[1];
    
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
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("ROIPlane_%1zu",fPlane));

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

void ROIFinderMorphological::FindROIs(const Waveform& waveform, size_t channel, size_t cnt, double rmsNoise, CandidateROIVec& roiVec) const
{
    // The idea here is to consider the input waveform - if an induction plane then it is already in differential form,
    // if a collection plane then we form the differential - then smooth and look for ROIs. The technique for actually
    // finding ROI's will be to compute the erosion and dilation vectors, get their average/difference and then use these to determine
    // candidate ROI's
    
    // Smooth the input waveform
    Waveform smoothWaveform;
    
    smoothInputWaveform(waveform, smoothWaveform);
    
    // We make lots of vectors... erosion, dilation, average and difference
    Waveform erosionVec;
    Waveform dilationVec;
    Waveform averageVec;
    Waveform differenceVec;
    
    // Make sure the temporary ROI hist is zeroed
    fROITempHist = 0;

    //getErosionDilationAverage(waveformDeriv, planeID.Plane, erosionVec, dilationVec, averageVec);
    getErosionDilationAverageDifference(waveform, smoothWaveform, channel, cnt, erosionVec, dilationVec, averageVec, differenceVec);
    
    // Use the average vector to find ROI's
    float truncRMS;
    float truncMean;
    
    fWaveformTool->getTruncatedMeanRMS(differenceVec, truncMean, truncRMS);
    
    if (fOutputHistograms)
    {
        Waveform::iterator maxItr = std::max_element(differenceVec.begin(),differenceVec.end());
        float maxDiff = *maxItr;
        float nSigma  = (maxDiff - truncMean) / std::max(float(0.5),truncRMS);
        
        fDiffMeanHist->Fill(truncMean, 1.);
        fDiffRmsHist->Fill(truncRMS, 1.);
        fDiffMaxHist->Fill(maxDiff, 1.);
        fNumSigmaHist->Fill(nSigma, 1.);
        
        if (std::distance(differenceVec.begin(),maxItr) > 4 * fHalfWindowSize)
        {
            maxDiff = *std::max_element(differenceVec.begin(),maxItr-4*fHalfWindowSize);
            nSigma  = (maxDiff - truncMean) / std::max(float(0.5),truncRMS);
            
            fNumSigNextHist->Fill(nSigma, 1.);
        }
        
        if (std::distance(maxItr, differenceVec.end()) > 4 * fHalfWindowSize)
        {
            maxDiff = *std::max_element(maxItr+4*fHalfWindowSize,differenceVec.end());
            nSigma  = (maxDiff - truncMean) / std::max(float(0.5),truncRMS);
            
            fNumSigNextHist->Fill(nSigma, 1.);
        }
    }
    
    float threshold = truncMean + fNumSigma * std::max(float(0.5),truncRMS);
    
    findROICandidatesDiff(differenceVec, erosionVec, dilationVec, 0, averageVec.size(), threshold, roiVec);
    
    if (roiVec.empty()) return;
    
    // pad the ROIs
    for(auto& roi : roiVec)
    {
        if (fROITempHist)
        {
            fROITempHist->Fill(int(roi.first),  std::max(5.*truncRMS,1.));
            fROITempHist->Fill(int(roi.second), std::max(5.*truncRMS,1.));
        }
        
        // low ROI end
        roi.first  = std::max(int(roi.first - fPreROIPad),0);
        // high ROI end
        roi.second = std::min(roi.second + fPostROIPad, waveform.size() - 1);
    }
    
    // merge overlapping (or touching) ROI's
    if(roiVec.size() > 1)
    {
        // temporary vector for merged ROIs
        CandidateROIVec tempRoiVec;
        
        // Loop through candidate roi's
        size_t startRoi = roiVec.front().first;
        size_t stopRoi  = startRoi;
        
        for(auto& roi : roiVec)
        {
            if (roi.first <= stopRoi + 50) stopRoi = roi.second;
            else
            {
                tempRoiVec.push_back(CandidateROI(startRoi,stopRoi));
                
                startRoi = roi.first;
                stopRoi  = roi.second;
            }
        }
        
        // Make sure to get the last one
        tempRoiVec.push_back(CandidateROI(startRoi,stopRoi));
        
        roiVec = tempRoiVec;
    }
    
    return;
}
    
void ROIFinderMorphological::findROICandidatesDiff(const Waveform&  differenceVec,
                                                   const Waveform&  erosionVec,
                                                   const Waveform&  dilationVec,
                                                   int              startTick,
                                                   int              stopTick,
                                                   float            threshold,
                                                   CandidateROIVec& roiCandidateVec) const
{
    int roiLength = stopTick - startTick;
    
    if (roiLength > 0)
    {
        // The idea here is to find the difference and use that as the seed for searching the
        // erosion and dilation vectors for the end points
        //Waveform::const_iterator maxItr = std::max_element(differenceVec.begin()+startTick,differenceVec.begin()+stopTick,[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
        Waveform::const_iterator maxItr = std::max_element(differenceVec.begin()+startTick,differenceVec.begin()+stopTick);

        int   maxTick       = std::distance(differenceVec.begin(),maxItr);
        float maxDifference = *maxItr;
        
        // No point continuing if not over threshold
        if (maxDifference > threshold)
        {
            // Start by finding maximum range of the erosion vector at this extremum
            int  maxCandRoiTick = maxTick;
        
            // The test on the erosion vector takes care of special case of slowly rising pulses which are characteristic of
            // tracks running into the wire plane
            while(++maxCandRoiTick < stopTick)
            {
                float difference = differenceVec.at(maxCandRoiTick);
                float erosion    = erosionVec.at(maxCandRoiTick);
            
                if (difference < threshold && erosion < 0.5 * threshold) break;
            }
        
            maxCandRoiTick = std::min(maxCandRoiTick,stopTick);
        
            // Now do the dilation vector
            int   minCandRoiTick = maxTick;
        
            while(--minCandRoiTick >= startTick)
            {
                float difference = differenceVec.at(minCandRoiTick);
                float erosion    = erosionVec.at(minCandRoiTick);
            
                if (difference < threshold && erosion < 0.5 * threshold) break;
            }
        
            // make sure we have not gone under
            minCandRoiTick = std::max(minCandRoiTick, startTick);
        
            int max2MinDiff = maxCandRoiTick - minCandRoiTick;
        
            if (fOutputHistograms && maxDifference > threshold)
            {
                fDeltaTicksHist->Fill(max2MinDiff, 1.);
                fDTixVDiffHist->Fill(max2MinDiff, maxDifference, 1.);
            }
        
            // Make sure we have a legitimate candidate
            if (max2MinDiff > fMax2MinDistance)
            {
                // Before saving this ROI, look for candidates preceeding this one
                // Note that preceeding snippet will reference to the current roiStartTick
                findROICandidatesDiff(differenceVec, erosionVec, dilationVec, startTick, minCandRoiTick, threshold, roiCandidateVec);
            
                // Save this ROI
                roiCandidateVec.push_back(CandidateROI(minCandRoiTick, maxCandRoiTick));
            
                // Now look for candidate ROI's downstream of this one
                findROICandidatesDiff(differenceVec, erosionVec, dilationVec, maxCandRoiTick+1, stopTick, threshold, roiCandidateVec);
            }
        }
    }
    
    return;
}
    
void ROIFinderMorphological::getErosionDilationAverageDifference(const Waveform& waveform,
                                                                 const Waveform& inputWaveform,
                                                                 size_t          channel,
                                                                 size_t          cnt,
                                                                 Waveform&       erosionVec,
                                                                 Waveform&       dilationVec,
                                                                 Waveform&       averageVec,
                                                                 Waveform&       differenceVec) const
{
    // Set the window size
    int halfWindowSize(fHalfWindowSize);
    
    // Are we making histograms?
    TProfile* origWaveHist   = 0;
    TProfile* waveformHist   = 0;
    TProfile* erosionHist    = 0;
    TProfile* dilationHist   = 0;
    TProfile* averageHist    = 0;
    TProfile* differenceHist = 0;
    
    bool fillHists(false);
    
    if (fOutputHistograms)
    {
        // Try to limit to the wire number (since we are already segregated by plane)
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
        size_t                   cryo = wids[0].Cryostat;
        size_t                   tpc  = wids[0].TPC;
        size_t                   wire = wids[0].Wire;
        
        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("ROIPlane_%1zu/c%1zu/c%1zut%1zuwire_%05zu",fPlane,cnt,cryo,tpc,wire));
        
        // We keep track of four histograms:
        try
        {
            origWaveHist   = dir.make<TProfile>(Form("Inp_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Waveform", waveform.size(),      0, waveform.size(),      -500., 500.);
            waveformHist   = dir.make<TProfile>(Form("Wfm_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Waveform", inputWaveform.size(), 0, inputWaveform.size(), -500., 500.);
            erosionHist    = dir.make<TProfile>(Form("Ero_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Erosion",  inputWaveform.size(), 0, inputWaveform.size(), -500., 500.);
            dilationHist   = dir.make<TProfile>(Form("Dil_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Dilation", inputWaveform.size(), 0, inputWaveform.size(), -500., 500.);
            averageHist    = dir.make<TProfile>(Form("Ave_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Average",  inputWaveform.size(), 0, inputWaveform.size(), -500., 500.);
            differenceHist = dir.make<TProfile>(Form("Dif_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Average",  inputWaveform.size(), 0, inputWaveform.size(), -500., 500.);
            
            // This is a kludge so that the ROI histogram ends up in the same diretory as the waveforms
            fROITempHist   = dir.make<TProfile>(Form("ROI_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "ROI",      inputWaveform.size(), 0, inputWaveform.size(), -500., 500.);
            
            fillHists = true;
        } catch(...)
        {
            std::cout << "Caught exception trying to make new hists, tpc,plane,cnt,wire: " << tpc << ", " << fPlane << ", " << cnt << ", " << wire << std::endl;
        }
    }

    // Initialize min and max elements
    std::pair<Waveform::const_iterator,Waveform::const_iterator> minMaxItr = std::minmax_element(inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);
    
    Waveform::const_iterator minElementItr = minMaxItr.first;
    Waveform::const_iterator maxElementItr = minMaxItr.second;
    
    // Initialize the erosion and dilation vectors
    erosionVec.resize(inputWaveform.size());
    dilationVec.resize(inputWaveform.size());
    averageVec.resize(inputWaveform.size());
    differenceVec.resize(inputWaveform.size());
    
    // Now loop through remaining elements and complete the vectors
    Waveform::iterator minItr = erosionVec.begin();
    Waveform::iterator maxItr = dilationVec.begin();
    Waveform::iterator aveItr = averageVec.begin();
    Waveform::iterator difItr = differenceVec.begin();

    for (std::vector<float>::const_iterator inputItr = inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
            
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        
        // Update the vectors
        *minItr++ = *minElementItr;
        *maxItr++ = *maxElementItr;
        *aveItr++ = 0.5 * (*maxElementItr + *minElementItr);
        *difItr++ = *maxElementItr - *minElementItr;
        
        if (fillHists)
        {
            int curBin = std::distance(inputWaveform.begin(),inputItr);
            
            origWaveHist->Fill(   curBin, waveform.at(curBin));
            waveformHist->Fill(   curBin, *inputItr);
            erosionHist->Fill(    curBin, *minElementItr);
            dilationHist->Fill(   curBin, *maxElementItr);
            averageHist->Fill(    curBin, 0.5*(*maxElementItr + *minElementItr));
            differenceHist->Fill( curBin,      *maxElementItr - *minElementItr);
        }
    }
    
    return;
}
    
void ROIFinderMorphological::smoothInputWaveform(const Waveform& inputWaveform, Waveform& outputWaveform) const
{
    // Vector smoothing - take the 10 bin average
    int   halfBins = fNumBinsToAve / 2;
    
    outputWaveform.resize(inputWaveform.size());
    
    // To facilitate handling the bins at the ends of the input waveform we embed in a slightly larger
    // vector which has zeroes on the ends
    Waveform tempWaveform(inputWaveform.size()+fNumBinsToAve);
    
    // Set the edge bins which can't be smoothed to zero
    std::fill(tempWaveform.begin(),tempWaveform.begin()+halfBins,0.);
    std::fill(tempWaveform.end()-halfBins,tempWaveform.end(),0.);
    
    // Copy in the input waveform
    std::copy(inputWaveform.begin(),inputWaveform.end(),tempWaveform.begin()+halfBins);
    
    // Now do the smoothing
    for(size_t idx = 0; idx < inputWaveform.size(); idx++)
    {
        float weightedSum(0.);
        
        for(size_t wIdx = 0; wIdx < fNumBinsToAve; wIdx++) weightedSum += fAveWeightVec.at(wIdx) * tempWaveform.at(idx + wIdx);
        
        outputWaveform.at(idx) = weightedSum / fWeightSum;
    }
    
    return;
}

void ROIFinderMorphological::initializeHistograms(art::TFileDirectory& histDir) const
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
    
DEFINE_ART_CLASS_TOOL(ROIFinderMorphological)
}
