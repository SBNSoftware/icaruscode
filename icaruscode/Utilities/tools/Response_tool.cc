////////////////////////////////////////////////////////////////////////
/// \file   Response.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/SignalShaping.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "art/Utilities/make_tool.h"
#include "icaruscode/Utilities/tools/IWaveformTool.h"
#include "IFieldResponse.h"
#include "IElectronicsResponse.h"
#include "IFilter.h"

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
    
    size_t                      getPlane()               const override {return fThisPlane;}
    
    const IFieldResponse*       getFieldResponse()       const override {return fFieldResponse.get();}
    const IElectronicsResponse* getElectronicsResponse() const override {return fElectronicsResponse.get();}
    const IFilter*              getFilter()              const override {return fFilter.get();}
    
    const util::SignalShaping&  getSignalShaping()       const override {return fSignalShaping;}
    
private:
    using IFieldResponsePtr       = std::unique_ptr<icarus_tool::IFieldResponse>;
    using IElectronicsResponsePtr = std::unique_ptr<icarus_tool::IElectronicsResponse>;
    using IFilterPtr              = std::unique_ptr<icarus_tool::IFilter>;

    // Utility routine for converting numbers to strings
    std::string             numberToString(int number);
    
    // Member variables from the fhicl file
    size_t                  fThisPlane;
    double                  f3DCorrection;
    double                  fTimeScaleFactor;
    int                     fDeconvPol;
    
    // Keep track of our base tools
    IFieldResponsePtr       fFieldResponse;
    IElectronicsResponsePtr fElectronicsResponse;
    IFilterPtr              fFilter;
    
    // The actual response function
    util::SignalShaping     fSignalShaping;
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
    fThisPlane       = pset.get<size_t>("Plane");
    f3DCorrection    = pset.get<size_t>("Correction3D");
    fTimeScaleFactor = pset.get<size_t>("TimeScaleFactor");
    fDeconvPol       = pset.get<int>("DeconvPol");
    
    // Build out the underlying tools we'll be using
    fFieldResponse       = art::make_tool<icarus_tool::IFieldResponse>(pset.get<fhicl::ParameterSet>("FieldResponse"));
    fElectronicsResponse = art::make_tool<icarus_tool::IElectronicsResponse>(pset.get<fhicl::ParameterSet>("ElectronicsResponse"));
    fFilter              = art::make_tool<icarus_tool::IFilter>(pset.get<fhicl::ParameterSet>("Filter"));
    
    return;
}
    
void Response::setResponse(double weight)
{
    // We'll need the FFT service
    art::ServiceHandle<util::LArFFT> fastFourierTransform;
    
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Recover the current set up
    std::string fftOptions  = fastFourierTransform->FFTOptions();
    size_t      nFFTFitBins = fastFourierTransform->FFTFitBins();
    size_t      fftSizeIn   = fastFourierTransform->FFTSize();
    size_t      fftSize     = fftSizeIn;
    
    // First of all set the field response
    fFieldResponse->setResponse(weight, f3DCorrection, fTimeScaleFactor);
    
    // Make sure the FFT can handle this
    size_t nFieldBins = fFieldResponse->getNumBins();
    
    // Reset the FFT if it is not big enough to handle current size
    if (nFieldBins * 4 > fftSize)
    {
        fftSize = 4 * nFieldBins;
        
        fastFourierTransform->ReinitializeFFT( fftSize, fftOptions, nFFTFitBins);
    }
        
    // handle the electronics response for this plane
    fElectronicsResponse->setResponse(fftSize, fFieldResponse->getBinWidth());
    
    // Set up the filter
    fFilter->setResponse(fftSizeIn, f3DCorrection, fTimeScaleFactor);
    
    // Add these elements to the SignalShaping class
    fSignalShaping.Reset();
    fSignalShaping.AddResponseFunction(fFieldResponse->getResponseVec());
    fSignalShaping.AddResponseFunction(fElectronicsResponse->getResponseVec());
    fSignalShaping.save_response();
    fSignalShaping.set_normflag(false);

    /* This could be a warning, but in principle, there's no reason to restrict the binning
     // Operation permitted only if output of rebinning has a larger bin size
     if( fFieldBinWidth > samplingRate )
     throw cet::exception(__FUNCTION__) << "\033[93m"
     << "Invalid operation: cannot rebin to a more finely binned vector!"
     << "\033[00m" << std::endl;
     */
    std::vector<double> SamplingTime( fftSize, 0. );
    
    for ( size_t itime = 0; itime < fftSize; itime++ )
        SamplingTime[itime] = double(itime) * detprop->SamplingRate();
    
    // Sampling
    // we want to implement new scheme (fStretchFullResponse==false) while retaining the old
    // time factor is already included in the calibrated response
    double timeFactor = f3DCorrection * fTimeScaleFactor;
    
    double timeFactorInv = 1. / timeFactor;
    
    const std::vector<double>* pResp = &fSignalShaping.Response_save();
    
    double deltaInputTime = fFieldResponse->getBinWidth();
    
    size_t nticks_input = pResp->size();
    
    std::vector<double> InputTime(nticks_input, 0. );
    for (size_t itime = 0; itime < nticks_input; itime++ )
        InputTime[itime] = double(itime) * deltaInputTime * timeFactor;
    
    std::vector<double> SamplingResp(fftSize, 0. );
    
    size_t SamplingCount = 0;
    
    size_t startJ = 1;
    SamplingResp[0] = (*pResp)[0];
    for ( size_t itime = 1; itime < fftSize; itime++ )
    {
        size_t low, high;
        for ( size_t jtime = startJ; jtime < nticks_input; jtime++ )
        {
            if ( InputTime[jtime] >= SamplingTime[itime] )
            {
                low  = jtime - 1;
                high = jtime;
                //            if(jtime<2&&itime<2) std::cout << itime << " " << jtime << " " << low << " " << up << std::endl;
                double interpolationFactor = ((*pResp)[high]-(*pResp)[low])/deltaInputTime;
                SamplingResp[itime] = ((*pResp)[low] + ( SamplingTime[itime] - InputTime[low] ) * interpolationFactor);
                // note: timeFactor = timeFactorInv =  1.0 for calibrated responses
                SamplingResp[itime] *= timeFactorInv;
                SamplingCount++;
                startJ = jtime;
                break;
            }
        } // for (  jtime = 0; jtime < nticks; jtime++ )
    } // for (  itime = 0; itime < nticks; itime++ )
    //std::cout << "SamplingResponse done " << std::endl;
    
    fSignalShaping.AddResponseFunction( SamplingResp, true);
    
    // Currently we only have fine binning "fFieldBinWidth"
    // for the field and electronic responses.
    // Now we are sampling the convoluted field-electronic response
    // with the nominal sampling.
    // We may consider to do the same for the filters as well.
    if (fftSizeIn != fftSize) fastFourierTransform->ReinitializeFFT(fftSizeIn, fftOptions, nFFTFitBins);
    
    // Finalize the Signal Shaping
    fSignalShaping.AddFilterFunction(fFilter->getResponseVec());
    fSignalShaping.SetDeconvKernelPolarity( fDeconvPol );
    fSignalShaping.CalculateDeconvKernel();
    
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
    const std::vector<double>& responseVec  = this->getSignalShaping().Response();
    double                     numBins      = responseVec.size();
    std::string                histName     = "Response_Plane_" + std::to_string(fThisPlane);
    TH1D*                      hist         = dir.make<TH1D>(histName.c_str(), "Response;Time(ticks)", numBins, 0., numBins);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        hist->Fill(bin, responseVec.at(bin));
    }

    // Get the FFT, need the waveform tool for consistency
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    std::unique_ptr<icarus_tool::IWaveformTool> waveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);
    
    std::vector<double> powerVec;
    
    waveformTool->getFFTPower(responseVec, powerVec);
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      maxFreq      = 500. / samplingRate;
    double      freqWidth    = maxFreq / powerVec.size();
    std::string freqName     = "Response_FFTPlane_" + std::to_string(fThisPlane);
    TH1D*       freqHist     = dir.make<TH1D>(freqName.c_str(), "Response;Frequency(MHz)", powerVec.size(), 0., maxFreq);
    
    for(size_t idx = 0; idx < powerVec.size(); idx++)
    {
        double freq = freqWidth * (idx + 0.5);
        
        freqHist->Fill(freq, powerVec.at(idx));
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
