////////////////////////////////////////////////////////////////////////
/// \file   ROIWavelets.cc
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
#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Filters/FFTFilterFunctions.h"
#include "icarus_signal_processing/Denoising.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TTree.h>
#include <TFile.h>

#include <fstream>
#include <numeric>

namespace icarus_tool
{

class ROIWavelets : public IROILocator
{
public:
    explicit ROIWavelets(const fhicl::ParameterSet& pset);
    
    ~ROIWavelets();
    
    void configure(const fhicl::ParameterSet& pset) override;
    void initializeHistograms(art::TFileDirectory&) override;
    
    void FindROIs(const art::Event&, const ArrayFloat&, const std::vector<raw::ChannelID_t>&, const geo::PlaneID&, ArrayFloat&, ArrayBool&) override;
    
private:
    // Define a structure to contain hits
    struct PeakCandidate {
      size_t startTick;
      size_t stopTick;
      size_t maxTick;
      size_t minTick;
    };

    using PeakCandidateVec = std::vector<PeakCandidate>;

    // This is for the baseline...
    float getMedian(const icarus_signal_processing::VectorFloat, const unsigned int) const;
    void  waveletFunc(const VectorFloat&,VectorFloat&,float,float) const;

    void findpeakCandidates(VectorFloat::const_iterator,
                            VectorFloat::const_iterator,
                            size_t,
                            double,
                            PeakCandidateVec&) const;

    // Parameters controlling tool
    bool                     fOutputHistograms;           ///< Diagnostic histogram output

    // fhicl parameters
    size_t                   fPlane;                      ///< What plane does this instance correspond to?
    float                    fScale;                      ///< scale parameters to use for ROI finding
    float                    fSigma;                      ///< This defines the range to use
    size_t                   fNSmoothBins;                ///< How many bins we smooth over (if zero no smoothing)
    float                    fThreshold;                  ///< Threshold to apply for saving signal

    // Keep track of the maximum range from scale and sigma
    size_t                   fMaxRange;                   ///< Maximum range between planes for scale, sigma

    // tuple output if requested
    std::vector<float>       fMedianVec;
    std::vector<float>       fRMSVec;
    std::vector<float>       fMinValVec;
    std::vector<float>       fMaxValVec;
    std::vector<float>       fRangeVec;
    std::vector<bool>        fHasROIVec;

    VectorFloat              fWavelet;            ///< container to hold our wavelets for each plane

    TTree*                   fTupleTree;          ///< output analysis tree

    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
ROIWavelets::ROIWavelets(const fhicl::ParameterSet& pset) : fOutputHistograms(false)
{
    configure(pset);
}
    
ROIWavelets::~ROIWavelets()
{
}
    
void ROIWavelets::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fPlane       = pset.get<size_t >("Plane",         0);
    fScale       = pset.get<float  >("WaveletScale", 15.);
    fSigma       = pset.get<float  >("WaveletSigma",  5.);
    fNSmoothBins = pset.get<float  >("NSmoothBins",  15.);
    fThreshold   = pset.get<float  >("Threshold",     7.);

    // Make sure the number of smooth bins is odd
    if (!(fNSmoothBins % 2)) fNSmoothBins += 1;

    fMaxRange = std::ceil(fSigma*fScale);

    VectorFloat xValVec(2*fMaxRange+1,0.);

    float xValue = -float(fMaxRange) - 1.;

    std::generate(xValVec.begin(),xValVec.end(),[&](){return xValue+=1.;});

    waveletFunc(xValVec, fWavelet, fScale, 0.);

    // use the below to verify the conditions for a wavelet
    float stepSize  = 0.01;
    
    xValVec.resize(std::ceil(xValVec.size()/stepSize),0.);

    xValue = -float(fMaxRange) - stepSize;

    std::generate(xValVec.begin(),xValVec.end(),[&](){return xValue+=0.01;});

    VectorFloat funcVals(xValVec.size());

    waveletFunc(xValVec,funcVals,fScale,0.);

    float average = std::accumulate(funcVals.begin(),funcVals.end(),0.,std::plus<float>()) * stepSize / float(funcVals.size());
    float power   = std::inner_product(funcVals.begin(),funcVals.end(),funcVals.begin(),0.) * stepSize;

    std::cout << "***** Wavelet definition for plane " << fPlane << " *****" << std::endl;
    std::cout << "      Scale: " << fScale << ", NSigma: " << fSigma << ", range: " << fMaxRange << std::endl;
    std::cout << "      Average: " << average << ", Power: " << power << std::endl;

    return;
}

void  ROIWavelets::waveletFunc(const VectorFloat& xVals,VectorFloat& funcVals,float scale,float translation) const
{
    // We define here a second order gaussian wavelet (colloquially known as "Mexican hat")
    const float normConst = 2 / std::sqrt(3 * std::sqrt(M_PI));

    if (xVals.size() != funcVals.size()) funcVals.resize(xVals.size(),0.);

    float sqrtScale = std::sqrt(scale);

    for(size_t idx = 0; idx < xVals.size(); idx++)
    {
        float arg = std::pow((xVals[idx]-translation)/scale,2);

        funcVals[idx] = normConst * (1 - arg) * std::exp(-0.5 * arg) / sqrtScale;
    }

    return;
}

void ROIWavelets::FindROIs(const art::Event& event, const ArrayFloat& constInputImage, const std::vector<raw::ChannelID_t>& channelVec, const geo::PlaneID& planeID, ArrayFloat& waveletWaveforms, ArrayBool& outputROIs)
{
    if (waveletWaveforms.size() != constInputImage.size()) waveletWaveforms.resize(constInputImage.size(),icarus_signal_processing::VectorFloat(constInputImage[0].size()));

    // Declare a holder for the input waveforms which has padding on each end 
    VectorFloat inputWaveform(constInputImage[0].size() + 2 * fMaxRange,0.);
    VectorFloat waveletVec(inputWaveform.size(),0.);

    // Try some smoothing
    // Set up to do a "triangle smoothing"
    size_t nSmoothBinsHalf = fNSmoothBins/2;

    VectorFloat smoothVec(fNSmoothBins);

    if (fNSmoothBins > 2)
    {
        for(size_t binIdx = 0; binIdx < nSmoothBinsHalf; binIdx++) 
        {
            smoothVec[binIdx]                    = float(binIdx + 1) / float(nSmoothBinsHalf);
            smoothVec[fNSmoothBins - binIdx - 1] = smoothVec[binIdx];
        }
    }

    smoothVec[nSmoothBinsHalf] = 1.;

    // Normalize it
    float smoothNorm = std::accumulate(smoothVec.begin(),smoothVec.end(),0.);

    std::transform(smoothVec.begin(),smoothVec.end(),smoothVec.begin(),[&](const auto& val){return val/smoothNorm;});

    // Loop through the input waveforms and apply the wavelet transform at the scale value we have chosen
    // Note that the wavelets have been pre-computed at initialization for each plane. 
    for(size_t channelIdx = 0; channelIdx < constInputImage.size(); channelIdx++)
    {
        // Recover the input waveform for this channel
        //const VectorFloat& waveform = constInputImage[channelIdx];
        const VectorFloat& waveform = constInputImage[channelIdx];

        // Copy to the working vector
        std::copy(waveform.begin(),waveform.end(),inputWaveform.begin() + fMaxRange);

        // If smoothing then do it now
        if (fNSmoothBins > 2)
        {
            for(size_t idx=0; idx<waveform.size()-smoothVec.size(); idx++)
            {
                float runAve = std::inner_product(waveform.begin()+idx,waveform.begin()+idx+smoothVec.size(),smoothVec.begin(),0.);

                inputWaveform[idx+nSmoothBinsHalf+fMaxRange] = runAve;
            }
        }

        // Explicitly zero elements
        std::fill(waveletVec.begin(),waveletVec.end(),0.);

        // Set the translation range 
        size_t upperBound = inputWaveform.size() - fWavelet.size();

        for(size_t translateIdx = 0; translateIdx < upperBound; translateIdx++)
        {
            for(size_t convolutionIdx = 0; convolutionIdx < fWavelet.size(); convolutionIdx++)
            {
                float convolutionValueAtIndex              = fWavelet[convolutionIdx] * inputWaveform[translateIdx + convolutionIdx] / 6.;
                waveletVec[translateIdx + convolutionIdx] += convolutionValueAtIndex * convolutionValueAtIndex;
            }
        }

        // Copy the waveletVec info to our output array
        std::copy(waveletVec.begin()+fMaxRange,waveletVec.end()-fMaxRange,waveletWaveforms[channelIdx].begin());
//        std::copy(inputWaveform.begin()+fMaxRange,inputWaveform.end()-fMaxRange,waveletWaveforms[channelIdx].begin());

        if (fOutputHistograms)
        {
            fMedianVec.clear();
            fRMSVec.clear();
            fMinValVec.clear();
            fMaxValVec.clear();
            fRangeVec.clear();
            fHasROIVec.clear();
        }

        // Set up to find peak candidates
        PeakCandidateVec peakCandidateVec;

        // Find them
        // Remember the padding that was applied, we search only in the waveform region
        findpeakCandidates(waveletVec.begin()+fMaxRange, waveletVec.end()-fMaxRange, 0, fThreshold, peakCandidateVec);

        // Right size the selected values array
        VectorBool& selVals = outputROIs[channelIdx];

        if (selVals.size() != waveform.size()) selVals.resize(waveform.size());

        std::fill(selVals.begin(),selVals.end(),false);

        // Go through the peak candidates and set the output accordingly
        for(const auto& peakCandidate : peakCandidateVec)
        {
            // Try to filter out false positives where we can be over threshold in wavelet power 
            // but have a negative excursion in the waveform
            if (waveform[peakCandidate.maxTick] < 0) continue;

            for(size_t idx = peakCandidate.startTick; idx < peakCandidate.stopTick; idx++) selVals[idx] = true;
        }

        if (fOutputHistograms)
        {
            VectorFloat rmsVec = waveletVec;
            size_t      maxIdx = 0.75 * rmsVec.size();

            std::nth_element(rmsVec.begin(), rmsVec.begin() + maxIdx, rmsVec.end());

            float rms    = std::sqrt(std::inner_product(rmsVec.begin(), rmsVec.begin() + maxIdx, rmsVec.begin(), 0.) / float(maxIdx));
            float minVal = *std::min_element(waveletVec.begin(),waveletVec.end());
            float maxVal = *std::max_element(waveletVec.begin(),waveletVec.end());
            float median = getMedian(waveletVec, waveletVec.size());
            
            fMedianVec.emplace_back(median);
            fRMSVec.emplace_back(rms);
            fMinValVec.emplace_back(minVal);
            fMaxValVec.emplace_back(maxVal);
            fRangeVec.emplace_back(maxVal-minVal);
            fHasROIVec.emplace_back(peakCandidateVec.size()>0);
        }
    }

    if (fOutputHistograms) fTupleTree->Fill();
     
    return;
}

float ROIWavelets::getMedian(icarus_signal_processing::VectorFloat vals, const unsigned int nVals) const
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
    
void ROIWavelets::initializeHistograms(art::TFileDirectory& histDir)
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.

    fTupleTree = histDir.make<TTree>("ROIFinder", "Tree by ROIWavelets_tool");

    fTupleTree->Branch("medians",  "std::vector<float>", &fMedianVec);
    fTupleTree->Branch("RMS",      "std::vector<float>", &fRMSVec);
    fTupleTree->Branch("minvals",  "std::vector<float>", &fMinValVec);
    fTupleTree->Branch("maxvals",  "std::vector<float>", &fMaxValVec);
    fTupleTree->Branch("range",    "std::vector<float>", &fRangeVec);
    fTupleTree->Branch("hasROI",   "std::vector<bool>",  &fHasROIVec);

    fOutputHistograms = true;
    
    return;
}

void ROIWavelets::findpeakCandidates(VectorFloat::const_iterator startItr,
                                     VectorFloat::const_iterator stopItr,
                                     size_t                      roiStartTick,
                                     double                      roiThreshold,
                                     PeakCandidateVec&           peakCandidateVec) const
{
    // Need a minimum number of ticks to do any work here
    if (std::distance(startItr,stopItr) > 4)
    {
        // Find the highest peak in the range given
        VectorFloat::const_iterator maxItr      = std::max_element(startItr, stopItr);
        size_t                      maxDistance = std::distance(startItr, maxItr);

        float maxValue = *maxItr;

        if (maxValue > roiThreshold)
        {
            // backwards to find first bin for this candidate hit
            VectorFloat::const_iterator firstItr = maxDistance > 2 ? maxItr - 1 : startItr;

            while(firstItr != startItr)
            {
                // Check for pathology where waveform goes too negative
//                if (*firstItr < -roiThreshold) break;

                // Check both sides of firstItr and look for min/inflection point
                if (*firstItr < *(firstItr+1) && *firstItr <= *(firstItr-1)) break;

                firstItr--;
            }

            int firstTick = std::distance(startItr,firstItr);

            // Recursive call to find all candidate hits earlier than this peak
            findpeakCandidates(startItr, firstItr + 1, roiStartTick, roiThreshold, peakCandidateVec);

            // forwards to find last bin for this candidate hit
            VectorFloat::const_iterator lastItr = std::distance(maxItr,stopItr) > 2 ? maxItr + 1 : stopItr - 1;

            while(lastItr != stopItr - 1)
            {
                // Check for pathology where waveform goes too negative
//                if (*lastItr < -roiThreshold) break;

                // Check both sides of firstItr and look for min/inflection point
                if (*lastItr <= *(lastItr+1) && *lastItr < *(lastItr-1)) break;

                lastItr++;
            }

            int lastTick = std::distance(startItr,lastItr);

            // Now save this candidate's start and max time info
            PeakCandidate peakCandidate;
            peakCandidate.startTick     = roiStartTick + firstTick;
            peakCandidate.stopTick      = roiStartTick + lastTick;
            peakCandidate.maxTick       = roiStartTick + maxDistance;
            peakCandidate.minTick       = roiStartTick + std::distance(startItr,std::min_element(firstItr,lastItr));

            peakCandidateVec.push_back(peakCandidate);

            // Recursive call to find all candidate hits later than this peak
            findpeakCandidates(lastItr + 1, stopItr, roiStartTick + std::distance(startItr,lastItr + 1), roiThreshold, peakCandidateVec);
        }
    }

    return;
}

DEFINE_ART_CLASS_TOOL(ROIWavelets)
}
