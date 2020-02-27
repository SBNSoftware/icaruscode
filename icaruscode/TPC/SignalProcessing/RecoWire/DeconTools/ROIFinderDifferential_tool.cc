////////////////////////////////////////////////////////////////////////
/// \file   ROIFinderDifferential.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IROIFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "icarussigproc/WaveformTools.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <fstream>

namespace icarus_tool
{

class ROIFinderDifferential : public IROIFinder
{
public:
    explicit ROIFinderDifferential(const fhicl::ParameterSet& pset);
    
    ~ROIFinderDifferential();
    
    void configure(const fhicl::ParameterSet& pset)                        override;
    void initializeHistograms(art::TFileDirectory&)                  const override;
    size_t plane()                                                   const override {return fPlane;}

    void FindROIs(const Waveform&, size_t, size_t, double, CandidateROIVec&) const override;
    
private:
    // The actual ROI finding algorithm
    void findROICandidates(Waveform::const_iterator startItr,
                           Waveform::const_iterator stopItr,
                           size_t                   channel,
                           size_t                   roiStartTick,
                           float                    roiThreshold,
                           CandidateROIVec&         roiCandidateVec) const;
    // Average/Smooth the input waveform
    void averageInputWaveform(const Waveform&, Waveform&) const;
    void smoothInputWaveform(const Waveform&, Waveform&)  const;

    // Member variables from the fhicl file
    size_t                                      fPlane;
    float                                       fNumSigma;                   ///< "# sigma" rms noise for ROI threshold
    int                                         fNumBinsToAve;               ///< Controls the averaging
    int                                         fMax2MinDistance;            ///< Maxmimum allow peak to peak distance
    size_t                                      fMaxTicksGap;                ///< Maximum gap between ROI's before merging
    unsigned short                              fPreROIPad;                  ///< ROI padding
    unsigned short                              fPostROIPad;                 ///< ROI padding
    float                                       fTruncRMSMinFraction;        ///< or at least this fraction of time bins
    bool                                        fOutputHistograms;           ///< Output histograms?

    std::vector<float>                          fAveWeightVec;
    float                                       fWeightSum;
    
    art::TFileDirectory*                        fHistDirectory;

    // Define some useful global histograms here
    TH1F*                                       fDiffMeanHist;
    TH1F*                                       fDiffRmsHist;
    
    TH1F*                                       fMinMaxDPeakHist;
    TH1F*                                       fMinMaxDistHist;
    TH1F*                                       fMinMaxSigmaHist;
    TH2F*                                       fDPeakVDistHist;

    icarussigproc::WaveformTools<float>         fWaveformTool;

    // Services
    const geo::GeometryCore*                    fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
ROIFinderDifferential::ROIFinderDifferential(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROIFinderDifferential::~ROIFinderDifferential()
{
}
    
void ROIFinderDifferential::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    std::vector<unsigned short> zin;
    
    fPlane               = pset.get< size_t >                     ("Plane"                     );
    fNumSigma            = pset.get< float >                      ("NumSigma"                  );
    fNumBinsToAve        = pset.get< int >                        ("NumBinsToAve"              );
    fMax2MinDistance     = pset.get< int >                        ("Max2MinDistance"           );
    fMaxTicksGap         = pset.get< size_t >                     ("MaxTicksGap",            50);
    zin                  = pset.get< std::vector<unsigned short> >("ROILeadTrailPadding"       );
    fTruncRMSMinFraction = pset.get< float >                      ("TruncRMSMinFraction",   0.6);
    fOutputHistograms    = pset.get< bool  >                      ("OutputHistograms",    false);

    // put the ROI pad sizes into more convenient vectors
    fPreROIPad  = zin[0];
    fPostROIPad = zin[1];
    
    // The "Max2MinDistance" is input in ticks but needs to be scaled if we are averaging
//    if (fNumBinsToAve > 1) fMax2MinDistance = std::round(fMax2MinDistance / fNumBinsToAve) + 1;

    // If asked, define some histograms
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;
        
        fHistDirectory = tfs.get();
        
        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("ROIPlane_%1zu",fPlane));
        
        fDiffMeanHist    = dir.make<TH1F>("DiffMean",    ";Diff Mean;",    100, -20.,  20.);
        fDiffRmsHist     = dir.make<TH1F>("DiffRms",     ";Diff RMS;",     100,   0.,   5.);
        
        fMinMaxDPeakHist = dir.make<TH1F>("MinMaxDPeak", ";Delta Peak;",   200, 0.,  50.);
        fMinMaxDistHist  = dir.make<TH1F>("MinMaxDist",  ";Peak Sep;",      50, 0.,  50.);
        fMinMaxSigmaHist = dir.make<TH1F>("MinMaxSigma", ";Peak Sep/rms;", 200, 0., 200.);
        fDPeakVDistHist  = dir.make<TH2F>("DPeakVDist",  ";Delta Peak; Peak Sep", 200, 0, 50., 50, 0., 50.);
    }

    // precalculate the weight vector to use in the averaging
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

void ROIFinderDifferential::FindROIs(const Waveform& waveform, size_t channel, size_t cnt, double rmsNoise, CandidateROIVec& roiVec) const
{
    // The idea here is to consider the input waveform - if an induction plane then it is already in differential form,
    // if a collection plane then we form the differential - then smooth and look for peaks. The technique will be to
    // consider the peak signature as a maximum followed some bins later by a mininum and those whose difference between
    // max and min is more than the threshold are kept.
    
    // First up, determine what kind of wire we have
    std::vector<geo::WireID> wids    = fGeometry->ChannelToWire(channel);
    const geo::PlaneID&      planeID = wids[0].planeID();
    geo::SigType_t           sigType = fGeometry->SignalType(planeID);
    
    // Local copy of the input waveform
    Waveform waveformDeriv(waveform.size());
    
    // If we have a collection plane then take the derivative
    if (sigType == geo::kCollection) fWaveformTool.firstDerivative(waveform, waveformDeriv);

    // Otherwise a straight copy since the bipolar pulses are, effectively, derivatives
    else std::copy(waveform.begin(),waveform.end(),waveformDeriv.begin());
    
    // Do the averaging
    Waveform aveWaveformDeriv;
    
//    averageInputWaveform(waveformDeriv, aveWaveformDeriv);
    smoothInputWaveform(waveformDeriv, aveWaveformDeriv);
//    fWaveformTool->triangleSmooth(waveform,aveWaveformDeriv);
    
    // Scheme for finding a suitable threshold
    float truncMean(0.);
    float truncRMS(0.);
    float fullRMS(0.);
    float nSig(2.5);
    int   nTrunc(0);
    
    fWaveformTool.getTruncatedMean(aveWaveformDeriv, truncMean, nTrunc);
    fWaveformTool.getTruncatedRMS(aveWaveformDeriv, nSig, fullRMS, truncRMS, nTrunc);
    
    // Put a floor on the value of the truncated RMS...
    float truncRMSFloor = std::max(truncRMS, float(0.25));

    // Now find the ROI's
    findROICandidates(aveWaveformDeriv.begin(),aveWaveformDeriv.end(),channel,0,truncRMSFloor,roiVec);
    
    TProfile* roiHist(0);

    if (fOutputHistograms)
    {
        fDiffMeanHist->Fill(truncMean, 1.);
        fDiffRmsHist->Fill(truncRMS, 1.);
        
        // Try to limit to the wire number (since we are already segregated by plane)
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
        size_t                   wire = wids[0].Wire;
        
        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("ROIPlane_%1zu/%03zu/wire_%05zu",fPlane,cnt,wire));
        
        // We keep track of four histograms:
        try
        {
            TProfile* waveformHist    = dir.make<TProfile>(Form("Wfm_%03zu_c%05zu",cnt,wire), "Waveform",   waveform.size(),         0, waveform.size(),         -500., 500.);
            TProfile* derivativeHist  = dir.make<TProfile>(Form("Der_%03zu_c%05zu",cnt,wire), "Derivative", waveformDeriv.size(),    0, waveformDeriv.size(),    -500., 500.);
            TProfile* aveDerivHist    = dir.make<TProfile>(Form("Ave_%03zu_c%05zu",cnt,wire), "AveDeriv",   aveWaveformDeriv.size(), 0, aveWaveformDeriv.size(), -500., 500.);
            
            roiHist = dir.make<TProfile>(Form("ROI_%03zu_c%05zu",cnt,wire), "ROI", waveform.size(), 0, waveform.size(), -500., 500.);

            for(size_t binIdx = 0; binIdx < waveform.size(); binIdx++)
            {
                waveformHist->Fill(binIdx, waveform.at(binIdx));
                derivativeHist->Fill(binIdx, waveformDeriv.at(binIdx));
            }

            int binIdx(0);
            
            for(const auto& val : aveWaveformDeriv) aveDerivHist->Fill(binIdx++, val);
            
        } catch(...)
        {
            std::cout << "Caught exception trying to make new hists, plane,cnt,wire: " << fPlane << ", " << cnt << ", " << wire << std::endl;
        }
    }

    if (roiVec.empty()) return;
    
    int nMultiplier = 1; // fNumBinsToAve
    
    // pad the ROIs
    for(auto& roi : roiVec)
    {
        if (roiHist)
        {
            roiHist->Fill(int(nMultiplier * roi.first),  std::max(3.*truncRMS,1.));
            roiHist->Fill(int(nMultiplier * roi.second), std::max(3.*truncRMS,1.));
        }
        
        // low ROI end
        roi.first  = std::max(int(nMultiplier * roi.first - fPreROIPad),0);
        // high ROI end
        roi.second = std::min(nMultiplier * roi.second + fNumBinsToAve + fPostROIPad, waveform.size() - 1);
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
            if (roi.first <= stopRoi + fMaxTicksGap) stopRoi = roi.second;
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
    
void ROIFinderDifferential::findROICandidates(Waveform::const_iterator startItr,
                                              Waveform::const_iterator stopItr,
                                              size_t                   channel,
                                              size_t                   roiStartTick,
                                              float                    truncRMS,
                                              CandidateROIVec&         roiCandidateVec) const
{
    // Require a minimum length
    size_t roiLength = std::distance(startItr,stopItr);
    
    if (roiLength > 0)
    {
        // The idea is to recover the two maxima from the input waveform and center our search for the ROI
        // around them.
        std::pair<Waveform::const_iterator,Waveform::const_iterator> minMaxItr = std::minmax_element(startItr,  stopItr);
        
        Waveform::const_iterator maxItr = minMaxItr.second;
        Waveform::const_iterator minItr = minMaxItr.first;

        // It can be that what has been returned are the lobes of two completely separate pulses which may be separated
        // by a large number of ticks. So we can't simply assume they belong together... unless they are "close" with
        // the minimum after the maximum
        if (std::distance(maxItr,minItr) < 0 || std::distance(maxItr,minItr) > fMax2MinDistance)
        {
            // Given that, we key off the larger of the lobes and do a forward or backward search to find the minimum/maximum
            // which we think is associated to the lobe we have chosen.
            // First consider that the maximum is larger than the minimum...
            if (*maxItr > std::fabs(*minItr))
            {
                // Check distance to the minimum we found
                if (std::distance(maxItr,stopItr) > 0)
                {
                    // The maximum is larger so search forward from here for the true minimum
                    minItr = maxItr;
                
                    float prevValue = *minItr++;
                
                    while(minItr != stopItr)
                    {
                        // Look for the point where it turns back up
                        if (prevValue < 0. && *minItr > prevValue)
                        {
                            // reset to what was the actual minimum value
                            minItr -= 1;
                            break;
                        }
                    
                        prevValue = *minItr++;
                    }
                }
            }
            else
            {
                // Check distance to the max
                if (std::distance(startItr,minItr) > 0)
                {
                    // Otherwise, we are starting at the minimum and searching backwards to find the max
                    maxItr = minItr;
                
                    float prevValue = *maxItr--;
                
                    while(maxItr != startItr)
                    {
                        // Decreasing for two bins
                        if (prevValue > 0. && *maxItr < prevValue)
                        {
                            // reset to what was the actual minimum value
                            maxItr += 1;
                            break;
                        }
                    
                        prevValue = *maxItr--;
                    }
                }
            }
        }
        
        // Check that the range from maximum to minimum is over threshold
        float maxMinRange    = *maxItr - *minItr;
        int   maxMinDistance = std::distance(maxItr,minItr);
        
        if (fOutputHistograms && maxMinRange > 2. * truncRMS)
        {
            fMinMaxDPeakHist->Fill(maxMinRange, 1.);
            fMinMaxDistHist->Fill(maxMinDistance, 1.);
            fMinMaxSigmaHist->Fill(maxMinRange/truncRMS, 1.);
            fDPeakVDistHist->Fill(maxMinRange, maxMinDistance, 1.);
        }
        
        if (maxMinRange > fNumSigma * truncRMS && maxMinDistance >= 0 && maxMinDistance < fMax2MinDistance)
        {
            // To complete the edges of the ROI, search both sides for the point which is essentially back to zero,
            // or in reality back into the rms level...
            Waveform::const_reverse_iterator revItr = std::find_if(std::make_reverse_iterator(maxItr), std::make_reverse_iterator(startItr), std::bind(std::less<float>(),std::placeholders::_1,truncRMS));

            maxItr = revItr.base();
            minItr = std::find_if(minItr,stopItr,std::bind(std::greater<float>(),std::placeholders::_1,-truncRMS));
        
            // Before saving this ROI, look for candidates preceeding this one
            // Note that preceeding snippet will reference to the current roiStartTick
            findROICandidates(startItr, maxItr, channel, roiStartTick, truncRMS, roiCandidateVec);
        
            // Save this ROI
            size_t newStartTick = std::distance(startItr,maxItr) + roiStartTick;
            size_t newStopTick  = std::distance(startItr,minItr) + roiStartTick;
        
            roiCandidateVec.push_back(CandidateROI(newStartTick, newStopTick));
        
            // Now look for candidate ROI's downstream of this one
            findROICandidates(minItr, stopItr, channel, newStopTick, truncRMS, roiCandidateVec);
        }
    }
    
    return;
}
    
void ROIFinderDifferential::averageInputWaveform(const Waveform& inputWaveform, Waveform& outputWaveform) const
{
    // Vector reduction - take the 10 bin average
    float aveSum(0.);
    
    if (outputWaveform.size() != inputWaveform.size()) outputWaveform.resize(inputWaveform.size());
    
    for(size_t idx = 0; idx < inputWaveform.size(); idx++)
    {
        aveSum += fAveWeightVec.at(idx % fNumBinsToAve) * inputWaveform.at(idx);
        
        if ((idx + 1) % fNumBinsToAve == 0)
        {
            outputWaveform[idx/fNumBinsToAve] = aveSum / fWeightSum;
            
            aveSum = 0.;
        }
    }
    
    return;
}
    
void ROIFinderDifferential::smoothInputWaveform(const Waveform& inputWaveform, Waveform& outputWaveform) const
{
    // Vector smoothing - take the 10 bin average
    int   halfBins = fNumBinsToAve / 2;
    
    outputWaveform.resize(inputWaveform.size());
    
    // Set the edge bins which can't be smoothed
    std::copy(inputWaveform.begin(),inputWaveform.begin()+halfBins,outputWaveform.begin());
    std::copy(inputWaveform.end()-halfBins,inputWaveform.end(),outputWaveform.end()-halfBins);
    
    for(size_t idx = halfBins; idx < inputWaveform.size() - halfBins; idx++)
    {
        float weightedSum(0.);
        
        for(int wIdx = 0; wIdx < fNumBinsToAve; wIdx++) weightedSum += fAveWeightVec.at(wIdx) * inputWaveform.at(idx - wIdx + halfBins);
        
        outputWaveform.at(idx) = weightedSum / fWeightSum;
    }
    
    return;
}

void ROIFinderDifferential::initializeHistograms(art::TFileDirectory& histDir) const
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
    
DEFINE_ART_CLASS_TOOL(ROIFinderDifferential)
}
