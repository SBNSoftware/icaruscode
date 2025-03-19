////////////////////////////////////////////////////////////////////////
/// \file   ROIFinder.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IROIFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"

#include "TH1D.h"

#include <fstream>
#include <numeric> // std::accumulate

namespace icarus_tool
{

class ROIFinderStandard : public IROIFinder
{
public:
    explicit ROIFinderStandard(const fhicl::ParameterSet& pset);
    
    ~ROIFinderStandard();
    
    void   configure(const fhicl::ParameterSet& pset)                        override;
    void   initializeHistograms(art::TFileDirectory&)                  const override;
    size_t plane()                                                   const override {return fPlane;}

    void FindROIs(const Waveform&, size_t, size_t, double, CandidateROIVec&) const override;
    
private:
    // Member variables from the fhicl file
    size_t                        fPlane;
    unsigned short                fNumBinsHalf;                ///< Determines # bins in ROI running sum
    std::vector<unsigned short>   fThreshold;                  ///< abs(threshold) ADC counts for ROI
    std::vector<int>              fNumSigma;                   ///< "# sigma" rms noise for ROI threshold
    std::vector<unsigned short>   fPreROIPad;                  ///< ROI padding
    std::vector<unsigned short>   fPostROIPad;                 ///< ROI padding
    
    // Services
    const geo::WireReadoutGeom* fChannelMapAlg = &art::ServiceHandle<geo::WireReadout const>()->Get();
    art::ServiceHandle<icarusutil::SignalShapingICARUSService> fSignalShaping;
};
    
//----------------------------------------------------------------------
// Constructor.
ROIFinderStandard::ROIFinderStandard(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROIFinderStandard::~ROIFinderStandard()
{
}
    
void ROIFinderStandard::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    std::vector<unsigned short> uin;
    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fNumBinsHalf = pset.get< unsigned short >             ("NumBinsHalf", 3);
    fThreshold   = pset.get< std::vector<unsigned short> >("Threshold"     );
    fNumSigma    = pset.get< std::vector<int> >           ("NumSigma"      );
    uin          = pset.get< std::vector<unsigned short> >("uPlaneROIPad"  );
    vin          = pset.get< std::vector<unsigned short> >("vPlaneROIPad"  );
    zin          = pset.get< std::vector<unsigned short> >("zPlaneROIPad"  );
    
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
    
    fSignalShaping = art::ServiceHandle<icarusutil::SignalShapingICARUSService>();
    
    return;
}
    
void ROIFinderStandard::FindROIs(const Waveform& waveform, size_t channel, size_t cnt, double rmsNoise, CandidateROIVec& roiVec) const
{
    // First up, translate the channel to plane
    std::vector<geo::WireID> wids    = fChannelMapAlg->ChannelToWire(channel);
    const geo::PlaneID&      planeID = wids[0].planeID();
    
    size_t numBins(2 * fNumBinsHalf + 1);
    size_t startBin(0);
    size_t stopBin(numBins);
    float  elecNoise = fSignalShaping->GetRawNoise(channel);
    float  rawNoise  = std::max(rmsNoise, double(elecNoise));
    
    float startThreshold = sqrt(float(numBins)) * (fNumSigma[planeID.Plane] * rawNoise + fThreshold[planeID.Plane]);
    float stopThreshold  = startThreshold;
    
    // Setup
    float runningSum = std::accumulate(waveform.begin(),waveform.begin()+numBins, 0.);
    
    size_t roiStartBin(0);
    bool   roiCandStart(false);
    
    // search for ROIs - follow prescription from Bruce B using a running sum to make faster
    // Note that we start in the middle of the running sum... if we find an ROI padding will extend
    // past this to take care of ends of the waveform
    for(size_t bin = fNumBinsHalf + 1; bin < waveform.size() - fNumBinsHalf; bin++)
    {
        // handle the running sum
        // Case, we are at start of waveform
        runningSum -= waveform[startBin++];
        
        // Case, we are at end of waveform
        runningSum += waveform[stopBin++];
        
        // We have already started a candidate ROI
        if (roiCandStart)
        {
            if (fabs(runningSum) < stopThreshold)
            {
                if (bin - roiStartBin > 2) roiVec.push_back(CandidateROI(roiStartBin, bin));
                
                roiCandStart = false;
            }
        }
        // Not yet started a candidate ROI
        else
        {
            if (fabs(runningSum) > startThreshold)
            {
                roiStartBin  = bin;
                roiCandStart = true;
            }
        }
    } // bin
    
    // add the last ROI if existed
    if (roiCandStart) roiVec.push_back(CandidateROI(roiStartBin, waveform.size() - 1));
    
    // pad the ROIs
    for(auto& roi : roiVec)
    {
        // low ROI end
        roi.first  = std::max(int(roi.first - fPreROIPad[planeID.Plane]),0);
        // high ROI end
        roi.second = std::min(roi.second + fPostROIPad[planeID.Plane], waveform.size() - 1);
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
            if (roi.first <= stopRoi) stopRoi = roi.second;
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

void ROIFinderStandard::initializeHistograms(art::TFileDirectory& histDir) const
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
    
DEFINE_ART_CLASS_TOOL(ROIFinderStandard)
}
