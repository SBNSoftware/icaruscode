////////////////////////////////////////////////////////////////////////
/// \file   ROIFinderDifferential.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/RecoWire/DeconTools/IROIFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "TH1D.h"

#include <fstream>

namespace uboone_tool
{

class ROIFinderDifferential : public IROIFinder
{
public:
    explicit ROIFinderDifferential(const fhicl::ParameterSet& pset);
    
    ~ROIFinderDifferential();
    
    void configure(const fhicl::ParameterSet& pset)                          override;
    void outputHistograms(art::TFileDirectory&)                        const override;
    
    void FindROIs(const Waveform&, size_t, double, CandidateROIVec&)   const override;
    
private:
    // The actual ROI finding algorithm
    void findROICandidates(Waveform::const_iterator       startItr,
                           Waveform::const_iterator       stopItr,
                           const geo::PlaneID::PlaneID_t& plane,
                           size_t                         roiStartTick,
                           float                          roiThreshold,
                           CandidateROIVec&               roiCandidateVec) const;
    // Average the input waveform
    void averageInputWaveform(const Waveform&, size_t, Waveform&) const;
    
    // grrrrr
    float fixTheFreakingWaveform(const Waveform&, Waveform&) const;
    
    // recover the rms
    float getTruncatedRMS(const Waveform&) const;
    
    // Member variables from the fhicl file
    std::vector<float>              fNumSigma;                   ///< "# sigma" rms noise for ROI threshold
    std::vector<int>                fNumBinsToAve;               ///< Controls the averaging
    std::vector<size_t>             fMax2MinDistance;            ///< Maxmimum allow peak to peak distance
    std::vector<unsigned short>     fPreROIPad;                  ///< ROI padding
    std::vector<unsigned short>     fPostROIPad;                 ///< ROI padding
    float                           fTruncRMSMinFraction;        ///< or at least this fraction of time bins
    
    std::vector<std::vector<float>> fAveWeightVec;
    std::vector<float>              fWeightSum;
    
    // Services
    const geo::GeometryCore*        fGeometry = lar::providerFrom<geo::Geometry>();
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
    std::vector<unsigned short> uin;
    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fNumSigma            = pset.get< std::vector<float> >         ("NumSigma"       );
    fNumBinsToAve        = pset.get< std::vector<int> >           ("NumBinsToAve"   );
    fMax2MinDistance     = pset.get< std::vector<size_t> >        ("Max2MinDistance");
    uin                  = pset.get< std::vector<unsigned short> >("uPlaneROIPad"   );
    vin                  = pset.get< std::vector<unsigned short> >("vPlaneROIPad"   );
    zin                  = pset.get< std::vector<unsigned short> >("zPlaneROIPad"   );
    fTruncRMSMinFraction = pset.get< float >                      ("TruncRMSMinFraction", 0.6);
    
    if(uin.size() != 2 || vin.size() != 2 || zin.size() != 2) {
        throw art::Exception(art::errors::Configuration)
        << "u/v/z plane ROI pad size != 2";
    }
    
    fPreROIPad.resize(3);
    fPostROIPad.resize(3);
    
    // put the ROI pad sizes into more convenient vectors
    fPreROIPad[0]  = uin[0];
    fPostROIPad[0] = uin[1];
    fPreROIPad[1]  = vin[0];
    fPostROIPad[1] = vin[1];
    fPreROIPad[2]  = zin[0];
    fPostROIPad[2] = zin[1];
    
    // precalculate the weight vector to use in the averaging
    fAveWeightVec.resize(3);
    fWeightSum.resize(3);
    
    for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
    {
        fAveWeightVec.at(planeIdx).resize(fNumBinsToAve.at(planeIdx));
    
        if (fNumBinsToAve.at(planeIdx) > 1)
        {
            for(int idx = 0; idx < fNumBinsToAve.at(planeIdx)/2; idx++)
            {
                float weight = idx + 1;
            
                fAveWeightVec.at(planeIdx).at(idx)                                  = weight;
                fAveWeightVec.at(planeIdx).at(fNumBinsToAve.at(planeIdx) - idx - 1) = weight;
            }
        }
        else fAveWeightVec.at(planeIdx).at(1) = 1.;
    
        fWeightSum.at(planeIdx) = std::accumulate(fAveWeightVec.at(planeIdx).begin(),fAveWeightVec.at(planeIdx).end(), 0.);
    }
    
    return;
}

void ROIFinderDifferential::FindROIs(const Waveform& waveform, size_t channel, double rmsNoise, CandidateROIVec& roiVec) const
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
    if (sigType == geo::kCollection)
    {
        waveformDeriv[0]                 = 0;
        waveformDeriv[waveform.size()-1] = 0;
        
        for(size_t idx = 1; idx < waveform.size()-1; idx++)
            waveformDeriv[idx] = 0.5 * (waveform.at(idx+1) - waveform.at(idx-1));
    }
    // Otherwise a straight copy since the bipolar pulses are, effectively, derivatives
    else std::copy(waveform.begin(),waveform.end(),waveformDeriv.begin());
    
    // Do the averaging
    Waveform aveWaveformDeriv;
    
    averageInputWaveform(waveformDeriv, planeID.Plane, aveWaveformDeriv);
    
    // Experiment with smoothing this as well
    // Now smooth the derivative vector
//    Waveform tempVec = aveWaveformDeriv;
    
//    for(size_t idx = 2; idx < tempVec.size() - 2; idx++)
//        aveWaveformDeriv.at(idx) = (tempVec.at(idx-2) + 2.*tempVec.at(idx-1) + 3.*tempVec.at(idx) + 2.*tempVec.at(idx+1) + tempVec.at(idx+2))/9.;
    
    // Scheme for finding a suitable threshold
    float truncRMS = getTruncatedRMS(aveWaveformDeriv);
    
    // Now find the ROI's
    findROICandidates(aveWaveformDeriv.begin(),aveWaveformDeriv.end(),planeID.Plane,0,truncRMS,roiVec);
    
    if (roiVec.empty()) return;
    
    // pad the ROIs
    for(auto& roi : roiVec)
    {
        // low ROI end
        roi.first  = std::max(int(fNumBinsToAve.at(planeID.Plane) * roi.first - fPreROIPad[planeID.Plane]),0);
        // high ROI end
        roi.second = std::min(fNumBinsToAve.at(planeID.Plane) * roi.second + fNumBinsToAve.at(planeID.Plane) + fPostROIPad[planeID.Plane], waveform.size() - 1);
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
    
void ROIFinderDifferential::findROICandidates(Waveform::const_iterator       startItr,
                                              Waveform::const_iterator       stopItr,
                                              const geo::PlaneID::PlaneID_t& plane,
                                              size_t                         roiStartTick,
                                              float                          truncRMS,
                                              CandidateROIVec&               roiCandidateVec) const
{
    // Require a minimum length
    size_t roiLength = std::distance(startItr,stopItr);
    
    if (roiLength > 0)
    {
        // The idea here is to find the largest positive value in the input waveform and use this as the basis of
        // our search for a candidate hit
        std::pair<Waveform::const_iterator,Waveform::const_iterator> minMaxItr = std::minmax_element(startItr,  stopItr);
        
        Waveform::const_iterator maxItr = minMaxItr.second;
        Waveform::const_iterator minItr = minMaxItr.first;

        // Reset either the max or min iterator depending on which is bigger
        if (*maxItr > std::fabs(*minItr))
        {
            // Check distance to the minimum we found
            if (std::distance(maxItr,minItr) > 2 || (std::distance(maxItr,minItr) < 0 && std::distance(maxItr,stopItr) > 0))
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
            if (std::distance(maxItr,minItr) > 2 || (std::distance(maxItr,minItr) < 0 && std::distance(startItr,minItr) > 0))
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
        
        // Check that the range from maximum to minimum is over threshold
        float maxMinRange    = *maxItr - *minItr;
        int   maxMinDistance = std::distance(maxItr,minItr);
        
        if (maxMinRange > fNumSigma.at(plane) * truncRMS && maxMinDistance >= 0 && maxMinDistance < int(fMax2MinDistance.at(plane)))
        {
            // to complete the edges of the ROI, search both sides for the point which is essentially back to zero
            // Somehow should be able to do this with reverse iterators but it seems startItr's constness is an issue
            while(maxItr != startItr)
            {
                if (*maxItr < 0.) break;
                maxItr--;
            }
        
            minItr = std::find_if(minItr,stopItr,std::bind(std::greater<float>(),std::placeholders::_1,0.));
        
            // Before saving this ROI, look for candidates preceeding this one
            // Note that preceeding snippet will reference to the current roiStartTick
            findROICandidates(startItr, maxItr, plane, roiStartTick, truncRMS, roiCandidateVec);
        
            // Save this ROI
            size_t newStartTick = std::distance(startItr,maxItr) + roiStartTick;
            size_t newStopTick  = std::distance(startItr,minItr) + roiStartTick;
        
            roiCandidateVec.push_back(CandidateROI(newStartTick, newStopTick));
        
            // Now look for candidate ROI's downstream of this one
            findROICandidates(minItr, stopItr, plane, newStopTick, truncRMS, roiCandidateVec);
        }
    }
    
    return;
}
    
void ROIFinderDifferential::averageInputWaveform(const Waveform& inputWaveform, size_t plane, Waveform& outputWaveform) const
{
    // Vector reduction - take the 10 bin average
    float                     aveSum(0.);
    int                       nBinsToAve(fNumBinsToAve.at(plane));
    const std::vector<float>& weightVec = fAveWeightVec.at(plane);
    float                     weightSum = fWeightSum.at(plane);
    
    outputWaveform.resize(inputWaveform.size()/nBinsToAve);
    
    for(size_t idx = 0; idx < inputWaveform.size(); idx++)
    {
        aveSum += weightVec.at(idx % nBinsToAve) * inputWaveform.at(idx);
        
        if ((idx + 1) % nBinsToAve == 0)
        {
            outputWaveform[idx/nBinsToAve] = aveSum / weightSum;
            
            aveSum = 0.;
        }
    }
    
    return;
}
    
float ROIFinderDifferential::getTruncatedRMS(const Waveform& waveform) const
{
    // do rms calculation - the old fashioned way and over all adc values
    Waveform locWaveform = waveform;
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + locWaveform.size()/2, locWaveform.begin(), 0.);
    
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(locWaveform.size()/2)));
    
    float threshold = 6. * localRMS;
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    int minNumBins = std::max(int(fTruncRMSMinFraction * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // Get the truncated sum
    localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(minNumBins)));
    
    return localRMS;
}
    
float ROIFinderDifferential::fixTheFreakingWaveform(const Waveform& waveform, Waveform& fixedWaveform) const
{
    // do rms calculation - the old fashioned way and over all adc values
    Waveform locWaveform = waveform;
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + locWaveform.size()/2, locWaveform.begin(), 0.);
    
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(locWaveform.size()/2)));
    
    float threshold = 6. * localRMS;
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    int minNumBins = std::max(int(fTruncRMSMinFraction * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // Get the truncated sum
    localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(minNumBins)));
    
    // Now get the average value
    float aveSum      = std::accumulate(locWaveform.begin(), locWaveform.begin() + minNumBins, 0.);
    float newPedestal = aveSum / minNumBins;
    
    std::transform(waveform.begin(), waveform.end(), fixedWaveform.begin(), [newPedestal](const auto& val){return val - newPedestal;});
    
    return localRMS;
}

    
void ROIFinderDifferential::outputHistograms(art::TFileDirectory& histDir) const
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
