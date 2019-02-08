////////////////////////////////////////////////////////////////////////
/// \file   Response2D.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/SignalShaping.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "art/Utilities/make_tool.h"
#include "icaruscode/Utilities/tools/IWaveformTool.h"
#include "IElectronicsResponse.h"
#include "IFilter.h"

#include "TProfile.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include <fstream>
#include <iomanip>

namespace icarus_tool
{

class Response2D : IResponse
{
public:
    explicit Response2D(const fhicl::ParameterSet& pset);
    
    ~Response2D() {}
    
    void configure(const fhicl::ParameterSet& pset)   override;
    void setResponse(double weight)                   override;
    void outputHistograms(art::TFileDirectory&) const override;
    
    size_t                      getPlane()               const override {return fThisPlane;}
    
    const IFieldResponse*       getFieldResponse()       const override {return nullptr;}
    const IElectronicsResponse* getElectronicsResponse() const override {return fElectronicsResponse.get();}
    const IFilter*              getFilter()              const override {return fFilter.get();}
    
    const util::SignalShaping&  getSignalShaping()       const override {return fSignalShaping;}
    
private:
    using IElectronicsResponsePtr = std::unique_ptr<icarus_tool::IElectronicsResponse>;
    using IFilterPtr              = std::unique_ptr<icarus_tool::IFilter>;

    // Utility routine for converting numbers to strings
    std::string             numberToString(int number);
    
    // Member variables from the fhicl file
    size_t                  fThisPlane;
    int                     fDeconvPol;
    
    // Keep track of our base tools
    IElectronicsResponsePtr fElectronicsResponse;
    IFilterPtr              fFilter;
    
    // The actual Response2D function
    util::SignalShaping     fSignalShaping;
};
    
//----------------------------------------------------------------------
// Constructor.
Response2D::Response2D(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
void Response2D::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fThisPlane = pset.get<size_t>("Plane");
    fDeconvPol = pset.get<int   >("DeconvPol");
    
    // Build out the underlying tools we'll be using
    fElectronicsResponse = art::make_tool<icarus_tool::IElectronicsResponse>(pset.get<fhicl::ParameterSet>("ElectronicsResponse"));
    fFilter              = art::make_tool<icarus_tool::IFilter>(pset.get<fhicl::ParameterSet>("Filter"));
    
    return;
}
    
void Response2D::setResponse(double weight)
{
    // We'll need the FFT service
    art::ServiceHandle<util::LArFFT> fastFourierTransform;
    
    // And the detector properties
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Recover the size/width
    size_t fftSize      = fastFourierTransform->FFTSize();
    double samplingRate = detprop->SamplingRate();    // We want this in ns/bin (note this is different from other response module!)

    // handle the electronics Response for this plane
    fElectronicsResponse->setResponse(fftSize, samplingRate);
    
    // Next we perform the convolution of the electronics Response and then invert to get the Response function
    const std::vector<double>& electronicsResponseVec = fElectronicsResponse->getResponseVec();
    
    float respIntegral = std::accumulate(electronicsResponseVec.begin(),electronicsResponseVec.end(),0.);
    
    mf::LogInfo("Response2D_tool")  << "      final Response2D integral: " << respIntegral << std::endl;

    fSignalShaping.Reset();
    fSignalShaping.AddResponseFunction( electronicsResponseVec, true);

    // Set up the filter
    fFilter->setResponse(fftSize, 1., 1.);
    
    // Finalize the Signal Shaping
    fSignalShaping.AddFilterFunction(fFilter->getResponseVec());
    fSignalShaping.SetDeconvKernelPolarity( fDeconvPol );         // Note that this is only important if set_normflag above has been set to true
    fSignalShaping.set_normflag(false);                           // WE are handling the normalization
    fSignalShaping.CalculateDeconvKernel();
    
    // The following can be uncommented to do some consistency checks if desired
//    // Do some consistency/cross checks here
//    // Check area of convolution function
//    const std::vector<TComplex>& convKernel = fSignalShaping.ConvKernel();
//
//    double normFactor = std::accumulate(convKernel.begin(),convKernel.end(),0.,[](const auto& val, double sum){return sum + std::abs(val);});
//
//    mf::LogInfo("Response2D_tool")  << "Response2D for plane: " << fThisPlane << ", convKernel integral: " << normFactor << std::endl;
//
//    const std::vector<TComplex>& deconvKernel = fSignalShaping.DeconvKernel();
//    std::vector<TComplex>  combKernel(deconvKernel.size());
//
//    std::transform(convKernel.begin(),convKernel.end(),deconvKernel.begin(),combKernel.begin(),std::multiplies<TComplex>());
//
//    const std::vector<TComplex>& filterKernel = fSignalShaping.Filter();
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
//    mf::LogInfo("Response2D_tool") << "Checking recovery of the filter, # differences: " << diffCount << ", max diff seen: " << maxRhoDiff << std::endl;

    return;
}
    
void Response2D::outputHistograms(art::TFileDirectory& histDir) const
{
    // Create a subfolder in which to place the "Response2D" histograms
    std::string thisResponse = "Response2DsPlane_" + std::to_string(fThisPlane);
    
    art::TFileDirectory dir = histDir.mkdir(thisResponse.c_str());
    
    // Do the field Response2D histograms
    fElectronicsResponse->outputHistograms(dir);
    fFilter->outputHistograms(dir);
    
    // Now make hists for the full Response2D
    std::string dirName = "Response2D_" + std::to_string(fThisPlane);
    
    art::TFileDirectory        responesDir  = dir.mkdir(dirName.c_str());
    const std::vector<double>& ResponseVec  = this->getSignalShaping().Response();
    auto const*                detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    int                        numBins      = ResponseVec.size();
    double                     samplingRate = detprop->SamplingRate(); // **Sampling time in ns**
    double                     maxFreq      = 1.e6 / (2. * samplingRate);
    double                     minFreq      = 1.e6 / (2. * samplingRate * double(numBins));
    std::string                histName     = "Response2D_Plane_" + std::to_string(fThisPlane);
    TProfile*                  hist         = dir.make<TProfile>(histName.c_str(), "Response2D;Time(us)", numBins, 0., numBins * samplingRate * 1.e-3);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        hist->Fill((double(bin) + 0.5) * samplingRate * 1.e-3, ResponseVec.at(bin), 1.);
    }

    // Get the FFT, need the waveform tool for consistency
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    std::unique_ptr<icarus_tool::IWaveformTool> waveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);
    
    std::vector<double> powerVec;
    
    waveformTool->getFFTPower(ResponseVec, powerVec);
    
    double      freqWidth = maxFreq / (powerVec.size() - 1);
    std::string freqName  = "Response2D_FFTPlane_" + std::to_string(fThisPlane);
    TProfile*   freqHist  = dir.make<TProfile>(freqName.c_str(), "Response2D;Frequency(MHz)", powerVec.size(), minFreq, maxFreq);
    
    for(size_t idx = 0; idx < powerVec.size(); idx++)
    {
        double freq = freqWidth * (idx + 0.5);
        
        freqHist->Fill(freq, powerVec.at(idx), 1.);
    }
    
    const std::vector<TComplex>& convKernel = this->getSignalShaping().ConvKernel();
    
    std::string convKernelName   = "ConvKernel_" + std::to_string(fThisPlane);
    TProfile*   fullResponseHist = dir.make<TProfile>(convKernelName.c_str(), "Convolution Kernel;Frequency(MHz)", convKernel.size(), minFreq, maxFreq);

    for(size_t idx = 0; idx < convKernel.size(); idx++)
    {
        double freq = freqWidth * (idx + 0.5);
        
        fullResponseHist->Fill(freq, convKernel.at(idx).Rho(), 1.);
    }
    
    const std::vector<TComplex>& deconKernel = this->getSignalShaping().DeconvKernel();
    
    std::string deconName = "DeconKernel_" + std::to_string(fThisPlane);
    TProfile*   deconHist = dir.make<TProfile>(deconName.c_str(), "Deconvolution Kernel;Frequency(MHz)", deconKernel.size(), minFreq, maxFreq);
    
    for(size_t idx = 0; idx < deconKernel.size(); idx++)
    {
        double freq = freqWidth * (idx + 0.5);
        
        deconHist->Fill(freq, deconKernel.at(idx).Rho(), 1.);
    }

    return;
}

std::string Response2D::numberToString(int number)
{
    std::ostringstream string;
    
    string << std::setfill('0') << std::setw(2) << number;
    
    return string.str();
}

    
DEFINE_ART_CLASS_TOOL(Response2D)
}
