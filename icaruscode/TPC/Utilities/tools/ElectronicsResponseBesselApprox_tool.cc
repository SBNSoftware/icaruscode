////////////////////////////////////////////////////////////////////////
/// \file   ElectronicsResponseBesselApprox.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IElectronicsResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "icarus_signal_processing/ICARUSFFT.h"

#include "TProfile.h"

#include <fstream>

namespace icarus_tool
{

class ElectronicsResponseBesselApprox : IElectronicsResponse
{
public:
    explicit ElectronicsResponseBesselApprox(const fhicl::ParameterSet& pset);
    
    ~ElectronicsResponseBesselApprox() {}
    
    void configure(const fhicl::ParameterSet& pset)         override;
    void setResponse(size_t numBins, double binWidth)       override;
    void outputHistograms(art::TFileDirectory&)       const override;
    
    size_t                          getPlane()           const override {return fPlane;}
    double                          getFCperADCMicroS()  const override {return fFCperADCMicroS;}
    double                          getASICShapingTime() const override {return fASICShapingTime;}
    const icarusutil::TimeVec&      getResponseVec()     const override {return fElectronicsResponseBesselApproxVec;}
    const icarusutil::FrequencyVec& getResponseFFTVec()  const override {return fElectronicsResponseFFTVec;}

private:
    // Member variables from the fhicl file
    size_t                   fPlane;
    double                   fFCperADCMicroS;
    double                   fASICShapingTime;
    double                   fADCPerPCAtLowestASICGain;
    double                   fTimeOffset;
    
    // Keep track of the bin width (for histograms)
    double                   fBinWidth;
    
    // Container for the electronics response "function"
    icarusutil::TimeVec      fElectronicsResponseBesselApproxVec;
    
    // And a container for the FFT of the above
    icarusutil::FrequencyVec fElectronicsResponseFFTVec;

    // Keep track of the FFT 
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>> fFFT; ///< Object to handle thread safe FFT
};
    
//----------------------------------------------------------------------
// Constructor.
ElectronicsResponseBesselApprox::ElectronicsResponseBesselApprox(const fhicl::ParameterSet& pset) :
    fBinWidth(0.),
    fFFT(nullptr)
{
    configure(pset);
}
    
void ElectronicsResponseBesselApprox::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fPlane                    = pset.get<size_t>("Plane");
    fFCperADCMicroS           = pset.get<double>("FCperADCMicroS");
    fASICShapingTime          = pset.get<double>("ASICShapingTime");
    fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");
    fTimeOffset               = pset.get<double>("TimeOffset");
    
    return;
}
    
void ElectronicsResponseBesselApprox::setResponse(size_t numBins, double binWidth)
{
    // The input binWidth is in nanoseconds, we need to convert to microseconds to match the
    // parameters for the electronics response which are given in us
    double timeCorrect  = 1.e-3;
    
    fBinWidth = binWidth * timeCorrect;
    
    fElectronicsResponseBesselApproxVec.resize(numBins, 0.);

    // Check that we have initialized our FFT object
    if (!fFFT) fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(numBins);
    
    // This note from Filippo:
    // The following sets the ICARUS electronics response function in
    // time-space. Function comes from BNL SPICE simulation of ICARUS
    // electronics. SPICE gives the electronics transfer function in
    // frequency-space. The inverse laplace transform of that function
    // (in time-space) was calculated in Mathematica and is what is being
    // used below. Parameters Ao and To are cumulative gain/timing parameters
    // from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain.
    // They have been adjusted to make the SPICE simulation to match the
    // actual electronics response. Default params are Ao=1.4, To=0.5us.
    
    for(size_t timeIdx = 0; timeIdx < numBins; timeIdx++)
    {
        double binTime = (double(timeIdx) + 0.5) * fBinWidth - fTimeOffset;
        double funcArg = 0.;
        
        if (binTime >= 0.) funcArg = binTime / fASICShapingTime;
        
        double leftArg2  = (0.9 * funcArg) * (0.9 * funcArg);
        double rightArg2 = (0.5 * funcArg) * (0.5 * funcArg);
        
        fElectronicsResponseBesselApproxVec[timeIdx] = (1. - exp(-0.5*leftArg2)) * exp(-0.5*rightArg2);
    }
    
    // normalize fElectResponse[i], before the convolution
    // Put in overall normalization in a pedantic way:
    // first put in the pulse area per eleectron at the lowest gain setting,
    // then normalize by the actual ASIC gain setting used.
    // This code is executed only during initialization of service,
    // so don't worry about code inefficiencies here.
    //    double last_integral=0;
    //    double last_max=0;
    
    // ICARUS Normalization are the following
    // Field response is normalized to 1 electron. Shaping function written as (t/tau)*exp(-t/tau) is normalized to 1
    // From test pulse measurement with FLIC@CERN we have 0.027 fC/(ADC*us)
    // Therefore 0.027*6242 electrons/(ADC*us)

    // The below scales the response by 1./FCperADCMicroS... but this gets taken out in the normalization
    std::transform(fElectronicsResponseBesselApproxVec.begin(),fElectronicsResponseBesselApproxVec.end(),fElectronicsResponseBesselApproxVec.begin(),std::bind(std::divides<double>(),std::placeholders::_1,fFCperADCMicroS));

    double respIntegral = std::accumulate(fElectronicsResponseBesselApproxVec.begin(),fElectronicsResponseBesselApproxVec.end(),0.);

    std::transform(fElectronicsResponseBesselApproxVec.begin(),fElectronicsResponseBesselApproxVec.end(),fElectronicsResponseBesselApproxVec.begin(),std::bind(std::divides<double>(),std::placeholders::_1,respIntegral));
    
    fFFT->forwardFFT(fElectronicsResponseBesselApproxVec, fElectronicsResponseFFTVec);

    return;
}
    
void ElectronicsResponseBesselApprox::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
    std::string dirName = "ElectronicsPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    std::string histName = "ElectronicsResponseBesselApprox_" + std::to_string(fPlane);
    
    double hiEdge = fElectronicsResponseBesselApproxVec.size() * fBinWidth;
    
    TProfile* hist = dir.make<TProfile>(histName.c_str(), "Response;Time(us)", fElectronicsResponseBesselApproxVec.size(), 0., hiEdge);
    
    for(size_t idx = 0; idx < fElectronicsResponseBesselApproxVec.size(); idx++) hist->Fill((double(idx)+0.5) * fBinWidth, fElectronicsResponseBesselApproxVec.at(idx), 1.);
   
    // Get the FFT of the response
    size_t halfFFTDataSize(fElectronicsResponseBesselApproxVec.size()/2);

    std::vector<double> powerVec(halfFFTDataSize);
    
    std::transform(fElectronicsResponseBesselApproxVec.begin(), fElectronicsResponseBesselApproxVec.begin() + halfFFTDataSize, powerVec.begin(), [](const auto& val){return std::abs(val);});
    
    // Now we can plot this...
    double maxFreq   = 0.5 / fBinWidth;   // binWidth will be in us, maxFreq will be units of MHz
    double freqWidth = maxFreq / powerVec.size();
    
    histName = "FFT_" + histName;
    
    TProfile* fftHist = dir.make<TProfile>(histName.c_str(), "Electronics FFT; Frequency(MHz)", powerVec.size(), 0., maxFreq);
    
    for(size_t idx = 0; idx < powerVec.size(); idx++)
    {
        double bin = (idx + 0.5) * freqWidth;
        
        fftHist->Fill(bin, powerVec.at(idx), 1.);
    }

    return;
}
    
DEFINE_ART_CLASS_TOOL(ElectronicsResponseBesselApprox)
}
