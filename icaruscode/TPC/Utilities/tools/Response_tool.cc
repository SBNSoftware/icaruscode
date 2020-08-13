////////////////////////////////////////////////////////////////////////
/// \file   Response.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "art/Utilities/make_tool.h"
#include "icarus_signal_processing/WaveformTools.h"
#include "IFieldResponse.h"
#include "IElectronicsResponse.h"
#include "IFilter.h"

#include "TProfile.h"

#include "icarus_signal_processing/ICARUSFFT.h"

#include <fstream>
#include <iomanip>

namespace icarus_tool
{

class Response : IResponse
{
public:
    explicit Response(const fhicl::ParameterSet& pset);
    
    ~Response() {}
    
    void configure(const fhicl::ParameterSet& pset)   override;
    void setResponse(double weight)                   override;
    void outputHistograms(art::TFileDirectory&) const override;
    
    size_t                                  getPlane()               const override {return fThisPlane;}
    
    const IFieldResponse*                   getFieldResponse()       const override {return fFieldResponse.get();}
    const IElectronicsResponse*             getElectronicsResponse() const override {return fElectronicsResponse.get();}
    const IFilter*                          getFilter()              const override {return fFilter.get();}

    size_t                                  getNumberTimeSamples()   const override {return fNumberTimeSamples;}
    const icarusutil::TimeVec&              getResponse()            const override {return fResponse;}
    const icarusutil::FrequencyVec&         getConvKernel()          const override {return fConvolutionKernel;}
    const icarusutil::FrequencyVec&         getDeconvKernel()        const override {return fDeconvolutionKernel;}
    double                                  getTOffset()             const override {return fT0Offset;};
    
private:
    // Calculate the response function
    void                                    calculateResponse(double weight);

    // Utility routine for converting numbers to strings
    std::string                             numberToString(int number);

    // Keep track of our status
    bool                                              fResponseHasBeenSet;
    
    // Member variables from the fhicl file
    size_t                                            fThisPlane;
    double                                            f3DCorrection;
    double                                            fTimeScaleFactor;
    int                                               fDeconvPol;
    
    using IFieldResponsePtr       = std::unique_ptr<icarus_tool::IFieldResponse>;
    using IElectronicsResponsePtr = std::unique_ptr<icarus_tool::IElectronicsResponse>;
    using IFilterPtr              = std::unique_ptr<icarus_tool::IFilter>;

    // Keep track of our base tools
    IFieldResponsePtr                                 fFieldResponse;
    IElectronicsResponsePtr                           fElectronicsResponse;
    IFilterPtr                                        fFilter;

    // Keep track of overall response functions
    size_t                                            fNumberTimeSamples;
    icarusutil::TimeVec                               fResponse;
    icarusutil::FrequencyVec                          fConvolutionKernel;
    icarusutil::FrequencyVec                          fDeconvolutionKernel;  

    double                                            fT0Offset;             ///< The overall T0 offset for the response function         

    std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>> fFFT;       ///< Object to handle thread safe FFT
    detinfo::DetectorProperties const*                fDetectorProperties;   ///< Detector properties service
};
    
//----------------------------------------------------------------------
// Constructor.
Response::Response(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
void Response::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fThisPlane           = pset.get<size_t>("Plane");
    f3DCorrection        = pset.get<size_t>("Correction3D");
    fTimeScaleFactor     = pset.get<size_t>("TimeScaleFactor");
    fDeconvPol           = pset.get<int   >("DeconvPol");

    fResponseHasBeenSet  = false;
    
    // Build out the underlying tools we'll be using
    fFieldResponse       = art::make_tool<icarus_tool::IFieldResponse>(pset.get<fhicl::ParameterSet>("FieldResponse"));
    fElectronicsResponse = art::make_tool<icarus_tool::IElectronicsResponse>(pset.get<fhicl::ParameterSet>("ElectronicsResponse"));
    fFilter              = art::make_tool<icarus_tool::IFilter>(pset.get<fhicl::ParameterSet>("Filter"));

    fDetectorProperties  = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fNumberTimeSamples   = fDetectorProperties->NumberTimeSamples();

    // Now set up our plans for doing the convolution
    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(fNumberTimeSamples);
    
    return;
}
    
void Response::setResponse(double weight)
{
    // If we have already done the setup then can return
    if (fResponseHasBeenSet) return;

    // Calculate the combined field and electronics shaping response
    calculateResponse(weight);

    // Now we compute the convolution kernel which is a straigtforward operation
    fFFT->forwardFFT(fResponse, fConvolutionKernel);

    // Set up the filter for use in the deconvolution
    fFilter->setResponse(fNumberTimeSamples, f3DCorrection, fTimeScaleFactor);

    // Now compute the deconvolution kernel
    fDeconvolutionKernel = fFilter->getResponseVec();

    for(size_t idx = 0; idx < fNumberTimeSamples; idx++)
    {
        if (std::abs(fConvolutionKernel[idx]) < 0.0001) fDeconvolutionKernel[idx]  = 0.;
        else                                            fDeconvolutionKernel[idx] /= fConvolutionKernel[idx];
    }
    
    // The following can be uncommented to do some consistency checks if desired
//    // Do some consistency/cross checks here
//    // Check area of convolution function
//    const std::vector<TComplex>& convKernel = fSignalShapingICARUS.ConvKernel();
//
//    double normFactor = std::accumulate(convKernel.begin(),convKernel.end(),0.,[](const auto& val, double sum){return sum + std::abs(val);});
//
//    mf::LogInfo("Response_tool")  << "Response for plane: " << fThisPlane << ", convKernel integral: " << normFactor << std::endl;
//
//    const std::vector<TComplex>& deconvKernel = fSignalShapingICARUS.DeconvKernel();
//    std::vector<TComplex>  combKernel(deconvKernel.size());
//
//    std::transform(convKernel.begin(),convKernel.end(),deconvKernel.begin(),combKernel.begin(),std::multiplies<TComplex>());
//
//    const std::vector<TComplex>& filterKernel = fSignalShapingICARUS.Filter();
//
//    int    diffCount(0);
//    double maxRhoDiff(0.);
//
//    for(size_t idx = 0; idx < filterKernel.size(); idx++)
//    {
//        double rhoDiff = filterKernel[idx].Rho() - combKernel[idx].Rho();
//
//        if (std::abs(rhoDiff) > 0.001) diffCount++;
//
//        if (std::abs(rhoDiff) > std::abs(maxRhoDiff)) maxRhoDiff = rhoDiff;
//    }
//
//    mf::LogInfo("Response_tool") << "Checking recovery of the filter, # differences: " << diffCount << ", max diff seen: " << maxRhoDiff << std::endl;

    fResponseHasBeenSet = true;

    return;
}

void Response::calculateResponse(double weight)
{
    // First of all set the field response
    fFieldResponse->setResponse(weight, f3DCorrection, fTimeScaleFactor);

    // What kind of response? (0= first induction, 1= middle induction, 2= collection)
    size_t responseType = fFieldResponse->getResponseType();

    // handle the electronics response for this plane
    fElectronicsResponse->setResponse(fFieldResponse->getResponseVec().size(), fFieldResponse->getBinWidth());
    
    // Next we perform the convolution of the field and electronics responses and then invert to get the
    // interim response function which will be set up in the time binning/range of the field response
    const icarusutil::FrequencyVec& fieldResponseFFTVec       = fFieldResponse->getResponseFFTVec();
    const icarusutil::FrequencyVec& electronicsResponseFFTVec = fElectronicsResponse->getResponseFFTVec();
    
    icarusutil::FrequencyVec curResponseFFTVec(fieldResponseFFTVec.size());
    
    std::transform(fieldResponseFFTVec.begin(), fieldResponseFFTVec.begin() + fieldResponseFFTVec.size()/2, electronicsResponseFFTVec.begin(), curResponseFFTVec.begin(), std::multiplies<icarusutil::ComplexVal>());
    
    // And now we recover the current response vector which is the convolution of the two above
    // (and still in the units of the original field response)
    icarusutil::TimeVec curResponseVec(fieldResponseFFTVec.size());
    
    // Note that we need a local version of the FFT because our time samples currently don't match what we will have
    icarus_signal_processing::ICARUSFFT<double> locFFT(curResponseVec.size());

    locFFT.inverseFFT(curResponseFFTVec, curResponseVec);
    
    // Now set to the task of determing the actual sampling response
    // We have to remember that the bin size for determining the field response probably
    // does not match that for the detector readout so we'll need to "convert"
    // from one to the other.
    fResponse.resize(fNumberTimeSamples, 0.);
    
    // Recover the combined response from above
//    const std::vector<double>& curResponseVec = fSignalShapingICARUS.Response_save();
    
    double respIntegral = std::accumulate(curResponseVec.begin(),curResponseVec.end(),0.);
    
    mf::LogInfo("Response_tool") << "***** Response for plane: " << fThisPlane << " ******" << "\n"
                                 << "      initial response integral: " << respIntegral << std::endl;
    
    // Need two factors: 1) the detector sampling rate and 2) the response sampling rate
    double samplingRate = fDetectorProperties->SamplingRate() * 1.e-3;    // We want this in us/bin
    double responseRate = fFieldResponse->getBinWidth() * 1.e-3;          // We want this in us/bin
    double rateRatio    = samplingRate / responseRate;                    // This gives the ratio of time bins for full readout to response bins
    
    // The idea is to step through each bin of the sampling response vector and then to
    // look up the corresponding bins in the current response vector. Since the two sample
    // rates are not the same there will be some "stretching" between the two. In addition,
    // we want to continue to allow for the possibility for further sample stretching
    double binScaleFactor = rateRatio * f3DCorrection * fTimeScaleFactor;
    
    // ok, do the loop
    for(size_t sampleIdx = 0; sampleIdx < fNumberTimeSamples; sampleIdx++)
    {
        // calculate the index for the response
        size_t responseLowIdx = std::floor(sampleIdx * binScaleFactor);
        
        if (responseLowIdx < curResponseVec.size())
        {
            // Calculate the index for the next bin
            size_t responseHiIdx = std::floor((sampleIdx + 1) * binScaleFactor);
            
            // This can't happen? But protect against zero divides...
            if (responseHiIdx == responseLowIdx) responseHiIdx += 1;
            
            if (responseHiIdx >= curResponseVec.size()) break;
            
            icarusutil::TimeVec::const_iterator curResponseItr = curResponseVec.begin();
            
            std::advance(curResponseItr,responseLowIdx);
            
            int nBins = responseHiIdx - responseLowIdx;
            
            // Obtain the average of these bins
            double aveResponse = std::accumulate(curResponseItr,curResponseItr+nBins,0.)/double(nBins);
            
            // Now interpolate between the two bins to get the sampling response for this bin
//            double responseSlope = (curResponseVec.at(responseHiIdx) - curResponseVec.at(responseLowIdx)) / (responseHiIdx - responseLowIdx);
//            double response      = curResponseVec.at(responseLowIdx) + 0.5 * responseSlope * (responseHiIdx - responseLowIdx);

            fResponse[sampleIdx] = aveResponse;
        }
    }
    
    // We need to scale by the binScaleFactor to preserve normalization
    std::transform(fResponse.begin(),fResponse.end(),fResponse.begin(),std::bind(std::multiplies<double>(),std::placeholders::_1,binScaleFactor));
    
    respIntegral = std::accumulate(fResponse.begin(),fResponse.end(),0.);

    // Now compute the T0 offset for the response function
    std::pair<icarusutil::TimeVec::iterator,icarusutil::TimeVec::iterator> minMaxPair = std::minmax_element(fResponse.begin(),fResponse.end());

    // Calculation of the T0 offset depends on the signal type
    int timeBin = std::distance(fResponse.begin(),minMaxPair.first);

    if (responseType == 2) timeBin = std::distance(fResponse.begin(),minMaxPair.second);
    
    // Do a backwards search to find the first positive bin
    while(1)
    {
        // Did we go too far?
        if (timeBin < 0)
            throw cet::exception("Response::configure") << "Cannot find zero-point crossover for induction response! ResponseType: " << responseType << ", plane: " << fThisPlane << std::endl;
            
        double content = fResponse[timeBin]; 
        
        if (content >= 0.) break;
        
        timeBin--;
    }

    // 
    fT0Offset = -timeBin;     // Note that this value being returned is in tick units now
    
    mf::LogInfo("Response_tool")  << "      final response integral: " << respIntegral << ", T0Offset: " << fT0Offset << std::endl;

    return;
}

void Response::outputHistograms(art::TFileDirectory& histDir) const
{
    // Create a subfolder in which to place the "response" histograms
    std::string thisResponse = "ResponsesPlane_" + std::to_string(fThisPlane);
    
    art::TFileDirectory dir = histDir.mkdir(thisResponse.c_str());
    
    // Do the field response histograms
    fFieldResponse->outputHistograms(dir);
    fElectronicsResponse->outputHistograms(dir);
    fFilter->outputHistograms(dir);
    
    // Now make hists for the full response
    std::string dirName = "Response_" + std::to_string(fThisPlane);
    
    art::TFileDirectory        responesDir  = dir.mkdir(dirName.c_str());
    const icarusutil::TimeVec& responseVec  = fResponse;
    auto const*                detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double                     numBins      = responseVec.size();
    double                     samplingRate = 1.e-3 * detprop->SamplingRate(); // Sampling time in us
    double                     maxFreq      = 1.e3 / samplingRate;      // Max frequency in MHz
    double                     minFreq      = maxFreq / numBins;
    std::string                histName     = "Response_Plane_" + std::to_string(fThisPlane);
    TProfile*                  hist         = dir.make<TProfile>(histName.c_str(), "Response;Time(us)", numBins, 0., numBins * samplingRate);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        hist->Fill((double(bin) + 0.5) * samplingRate, responseVec.at(bin), 1.);
    }
    
    icarusutil::TimeVec powerVec;
    
    fFFT->getFFTPower(responseVec, powerVec);
    
    std::string freqName  = "Response_FFTPlane_" + std::to_string(fThisPlane);
    TProfile*   freqHist  = dir.make<TProfile>(freqName.c_str(), "Response;Frequency(kHz)", numBins/2, minFreq, 0.5*maxFreq);
    
    for(size_t idx = 0; idx < numBins/2; idx++)
    {
        double freq = minFreq * (idx + 0.5);
        
        freqHist->Fill(freq, powerVec.at(idx), 1.);
    }
    
    const icarusutil::FrequencyVec& convKernel = fConvolutionKernel;
    
    std::string convKernelName   = "ConvKernel_" + std::to_string(fThisPlane);
    TProfile*   fullResponseHist = dir.make<TProfile>(convKernelName.c_str(), "Convolution Kernel;Frequency(kHz)", numBins/2, minFreq, 0.5*maxFreq);

    for(size_t idx = 0; idx < numBins/2; idx++)
    {
        double freq = minFreq * (idx + 0.5);
        
        fullResponseHist->Fill(freq, std::abs(convKernel[idx]), 1.);
    }
    
    const icarusutil::FrequencyVec& deconKernel = fDeconvolutionKernel;
    
    std::string deconName = "DeconKernel_" + std::to_string(fThisPlane);
    TProfile*   deconHist = dir.make<TProfile>(deconName.c_str(), "Deconvolution Kernel;Frequency(kHz)", numBins/2, minFreq, 0.5*maxFreq);
    
    for(size_t idx = 0; idx < numBins/2; idx++)
    {
        double freq = minFreq * (idx + 0.5);
        
        deconHist->Fill(freq, std::abs(deconKernel[idx]), 1.);
    }

    return;
}

std::string Response::numberToString(int number)
{
    std::ostringstream string;
    
    string << std::setfill('0') << std::setw(2) << number;
    
    return string.str();
}

    
DEFINE_ART_CLASS_TOOL(Response)
}
