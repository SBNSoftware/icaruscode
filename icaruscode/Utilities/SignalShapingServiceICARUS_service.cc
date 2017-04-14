////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceICARUS_service.cc
/// \author H. Greenlee
/// Modified by X. Qian 1/6/2015
/// if histogram is used, inialize
/// Response_Offset, Response_Sampling, FieldBins from histogram
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "SignalShapingServiceICARUS.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "art/Framework/Services/Optional/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/LArFFT.h"
#include "TFile.h"

// LArSoft include
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"

#include "art/Utilities/make_tool.h"
#include "tools/IResponse.h"
#include "tools/IFieldResponse.h"
#include "tools/IElectronicsResponse.h"
#include "tools/IFilter.h"

#include <fstream>

//----------------------------------------------------------------------
// Constructor.
util::SignalShapingServiceICARUS::SignalShapingServiceICARUS(const fhicl::ParameterSet& pset,
                                                                     art::ActivityRegistry& /* reg */)
: fInit(false)
{
    for(size_t i=0; i<3; ++i)
    {
        fHistDoneF[i] = false;
    }
    
    reconfigure(pset);
}

//----------------------------------------------------------------------
// Destructor
util::SignalShapingServiceICARUS::~SignalShapingServiceICARUS()
{
    std::cout << "In SignalShapingServiceICARUS destructor" << std::endl;
    std::cout << "Filter vec size: " << fFilterVec.size() << std::endl;
    for(size_t idx = 0; idx < fFilterVec.size(); idx++)
    {
        std::cout << "--> size: " << fFilterVec[idx].size() << std::endl;
        fFilterVec[idx].clear();
    }
    fFilterVec.clear();
    return;
}

//----------------------------------------------------------------------
// Reconfigure method.
void util::SignalShapingServiceICARUS::reconfigure(const fhicl::ParameterSet& pset)
{
    // If called again, then we need to clear out the existing tools...
    fPlaneToResponseMap.clear();
    
    // Implement the tools for handling the responses
    const fhicl::ParameterSet& responseTools = pset.get<fhicl::ParameterSet>("ResponseTools");
    
    for(const std::string& responseTool : responseTools.get_pset_names())
    {
        const fhicl::ParameterSet& responseToolParamSet = responseTools.get<fhicl::ParameterSet>(responseTool);
        size_t                     planeIdx             = responseToolParamSet.get<size_t>("Plane");
        
        fPlaneToResponseMap[planeIdx].push_back(art::make_tool<icarus_tool::IResponse>(responseToolParamSet));
    }
    
    // add a comment here
    art::ServiceHandle<geo::Geometry> geo;
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    // Reset initialization flag.
    
    fInit = false;
    
    fPlaneForNormalization  = pset.get<size_t>("PlaneForNormalization");
    fPrintResponses         = pset.get<bool>("PrintResponses");
    fDeconNorm              = pset.get<double>("DeconNorm");
    fDefaultDriftVelocity   = pset.get< DoubleVec >("DefaultDriftVelocity");
    fCalibResponseTOffset   = pset.get< DoubleVec >("CalibResponseTOffset");
    fNoiseFactVec           = pset.get<DoubleVec2>("NoiseFactVec");
    f3DCorrectionVec        = pset.get<DoubleVec>("Drift3DCorrVec");
    fDeconvPol              = pset.get<std::vector<int> >("DeconvPol");
    
    //Adding calibrated field response at 70kV
    fUseCalibratedResponses = pset.get<bool>("UseCalibratedResponses");
    
    mf::LogInfo("SignalShapingServiceICARUS") << " using the field response provided from a .root file " ;
    
    // constructor decides if initialized value is a path or an environment variable
    fDefaultEField          = pset.get<double>("DefaultEField");
    fDefaultTemperature     = pset.get<double>("DefaultTemperature");
    fTimeScaleParams        = pset.get<DoubleVec>("TimeScaleParams");
    fStretchFullResponse    = pset.get<bool>("StretchFullResponse");
    
    // calculate the time scale factor for this job
    if(!fUseCalibratedResponses) SetTimeScaleFactor();
    
    // Reset kernels.
    std::cout << " reset kernels " << std::endl;
    
    size_t ktype = 0;
    fSignalShapingVec.resize(2);
    for(auto& kset : fSignalShapingVec)
    {
        kset.resize(geo->Nplanes());
        size_t planeIdx = 0;
        for(auto& plane : kset)
        {
            std::cout << " resizing shapings " << ktype << " " << planeIdx << std::endl;
            plane.Reset();
            planeIdx++;
        }
        ktype++;
    }
    
    fFieldResponseVec.resize(2);
    for(auto& kset : fFieldResponseVec) kset.resize(geo->Nplanes());

    /*
     We allow different drift velocities.
     kDVel is ratio of what was used in LArG4 to field response simulation.
     If drift velocity used for field response is set to <0, then we assume
     the same drift velocity as used in LArG4.
     */
    for(size_t plane = 0; plane < geo->Nplanes(); ++plane) {
        
        double larg4_velocity = detprop->DriftVelocity( detprop->Efield(plane), detprop->Temperature() );
        
        if(fDefaultDriftVelocity.at(plane) < 0) fDefaultDriftVelocity.at(plane) = larg4_velocity;
    }
    
    return;
}

void util::SignalShapingServiceICARUS::SetTimeScaleFactor()
{
    // get the scale factor between the bulk drift velocity used to generate the field response
    //   and that used for this simulation.
    
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    double defaultVelocity = detprop->DriftVelocity(fDefaultEField, fDefaultTemperature);
    double thisVelocity    = detprop->DriftVelocity( detprop->Efield(0), detprop->Temperature() );
    double vRatio = defaultVelocity/thisVelocity;
    double vDiff = vRatio -1.0;
    
    fTimeScaleFactor = 0.0;
    double term = 1.0;
    
    // the time scale params are from a fit to Garfield simulations at different E Fields
    for(size_t i = 0;i<fTimeScaleParams.size(); ++i) {
        fTimeScaleFactor += fTimeScaleParams[i]*term;
        term *= vDiff;
    }
   // fTimeScaleFactor=1;
    std::cout << "Current E field = " << detprop->Efield(0) << " KV/cm, Ratio of drift velocities = " << vRatio << ", timeScaleFactor = " << fTimeScaleFactor << std::endl;
    
}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceICARUS::SignalShaping(size_t channel, size_t ktype) const
{
    if(!fInit) init();
    
    art::ServiceHandle<geo::Geometry> geom;
    
        //use channel number to set some useful numbers
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;

    return fSignalShapingVec[ktype][planeIdx];
}


//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShapingServiceICARUS::init()
{
    if(!fInit)
    {
        fInit = true;
        
        art::ServiceHandle<geo::Geometry> geo;
        
        // Do ICARUS-specific configuration of SignalShaping by providing
        // ICARUS response and filter functions.
        
        // re-initialize the FFT service for the request size
        art::ServiceHandle<util::LArFFT> fFFT;
        std::string options = fFFT->FFTOptions();
        int fftsize = (int) fFFT->FFTSize();
        std::cout << " fftsize " << fftsize << std::endl;
        // Calculate field and electronics response functions.
        
        std::string kset[2] = { "Convolution ", "Deconvolution "};
        
        fElectResponse.resize(2);
        for(auto& electByPlane : fElectResponse) electByPlane.resize(geo->Nplanes());
        
        fFilterVec.resize(2);
        for(auto& filter : fFilterVec) filter.resize(geo->Nplanes());
        
        // Get the normalization from the field response for the collection plane
        double integral = fPlaneToResponseMap.at(fPlaneForNormalization).front().get()->getFieldResponse()->getIntegral();
        double weight   = 1. / integral;
        
        for(size_t ktype=0;ktype<2;++ktype) {
            std::cout << std::endl << kset[ktype] << " setting functions:" << std::endl;
            int fftsize2 = (int) fFFT->FFTSize();
            std::cout << " fftsize2 " << fftsize2 << std::endl;
            
            std::cout << "Input field responses" << std::endl;
            
            for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
            {
                fPlaneToResponseMap[planeIdx].front().get()->setResponse(weight);
                
                fSignalShapingVec[ktype][planeIdx] = fPlaneToResponseMap[planeIdx].front().get()->getSignalShaping();
                
            }
        }
    }
    
    return;
}

void util::SignalShapingServiceICARUS::SetDecon(size_t fftsize, size_t channel)
{
    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<util::LArFFT>  fFFT;
    
    // This is for deconvolution, ktype = 1
    size_t ktype = 1;
    
    if (!fInit) init();
    else
    {
        // Recover the plane for this channel
        size_t planeIdx = geo->ChannelToWire(channel)[0].Plane;
        
        // Set up the filter
        const icarus_tool::IFilter* filterTool = fPlaneToResponseMap.at(planeIdx).front().get()->getFilter();
        
//        filterTool->setResponse(fFFT->FFTSize(), f3DCorrectionVec[planeIdx], fTimeScaleFactor);
        
        fFilterVec[ktype][planeIdx] = filterTool->getResponseVec();
    }
    
    
    // streamline this method:
    // if the deconvolution kernel is already appropriate for the datasize (aka fftsize) do nothing
    // otherwise, set it to the appropriate size
    // do this test for *every* ss
    // But it will in general only happen once per run!
    
    bool setDecon = false;
    
    size_t FFTSize = fFFT->FFTSize();
    if (fftsize>FFTSize||fftsize<=FFTSize/2)
    {
        std::string options = fFFT->FFTOptions();
        int fitbins = fFFT->FFTFitBins();
        fFFT->ReinitializeFFT( (size_t)fftsize, options, fitbins);
        setDecon = true;
    }
    
    if(!setDecon) return;
    
    for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
    {
        fSignalShapingVec[ktype][planeIdx].Reset();
        fSignalShapingVec[ktype][planeIdx].AddFilterFunction(fFilterVec[ktype][planeIdx]);
        fSignalShapingVec[ktype][planeIdx].SetDeconvKernelPolarity( fDeconvPol.at(planeIdx));
        fSignalShapingVec[ktype][planeIdx].CalculateDeconvKernel();
    }
}

//----------------------------------------------------------------------
// Calculate ICARUS field response.
void util::SignalShapingServiceICARUS::SetFieldResponse(size_t ktype)
{
    // Get services.
    art::ServiceHandle<geo::Geometry> geo;
    
    ////////////////////////////////////////////////////
    art::ServiceHandle<art::TFileService> tfs;
    
    // Ticks in nanosecond
    // Calculate the normalization of the collection plane
    double integral = fPlaneToResponseMap.at(fPlaneForNormalization).front().get()->getFieldResponse()->getIntegral();
    double weight   = 1. / integral;
    
    std::cout << " Integral " << integral << " weight " << weight << std::endl;
    // we adjust the size of the fieldresponse vector to account for the stretch
    // and interpolate the histogram to fill the vector with the stretched response
    
    for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
    {
        double timeFactor = 1.0;
        if(!fUseCalibratedResponses) timeFactor *= f3DCorrectionVec[planeIdx];
        std::cout << " after 3d corr " << timeFactor << std::endl;
        if(!fStretchFullResponse) timeFactor *= fTimeScaleFactor;
        std::cout << " after scale corr " << timeFactor << std::endl;
        std::cout << " ktype " << ktype << " plane " << planeIdx << std::endl;
        
        const icarus_tool::IFieldResponse* fieldResponsePtr = fPlaneToResponseMap.at(planeIdx).front().get()->getFieldResponse();
        
//        fieldResponsePtr->setResponse(weight, f3DCorrectionVec[planeIdx], fTimeScaleFactor);
        
        fFieldResponseVec[ktype][planeIdx] = fieldResponsePtr->getResponseVec();
    }
    
    std::cout << " end SetFieldResponse " << std::endl;
    return;
}

//----------------------------------------------------------------------
// Sample ICARUS response (the convoluted field and electronic response), will probably add the filter later
void util::SignalShapingServiceICARUS::SetResponseSampling(size_t ktype, int mode, size_t channel)
{
    // Get services
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* geo = lar::providerFrom<geo::Geometry>();
    art::ServiceHandle<util::LArFFT> fft;
    art::ServiceHandle<art::TFileService> tfs;
    
    /* This could be a warning, but in principle, there's no reason to restrict the binning
     // Operation permitted only if output of rebinning has a larger bin size
     if( fFieldBinWidth > samplingRate )
     throw cet::exception(__FUNCTION__) << "\033[93m"
     << "Invalid operation: cannot rebin to a more finely binned vector!"
     << "\033[00m" << std::endl;
     */
    std::cout << "entering SetResponseSampling, ktype/config/mode/channel " << ktype << " " << mode << " " << channel << std::endl;
    
    size_t plane0 = 0;
    size_t plane1 = geo->Nplanes();
    
    if(mode != 0)
    {
        size_t plane = geo->ChannelToWire(channel)[0].Plane;
        
        plane0 = plane;
        plane1 = std::min(size_t(geo->Nplanes()),plane+1);
    }
    
    size_t nticks = fft->FFTSize();
    DoubleVec SamplingTime( nticks, 0. );

    for ( size_t itime = 0; itime < nticks; itime++ )
        SamplingTime[itime] = (1.*itime) * detprop->SamplingRate();
    
    // Sampling
    
    // we want to implement new scheme (fStretchFullResponse==false) while retaining the old
    // time factor is already included in the calibrated response
    for(size_t planeIdx = plane0; planeIdx < plane1; planeIdx++)
    {
        double timeFactor = 1.0;
        //if(fStretchFullResponse && !fUseCalibratedResponses) timeFactor *= fTimeScaleFactor*f3DCorrectionVec[view];
        
        if (!fUseCalibratedResponses) timeFactor *= f3DCorrectionVec[planeIdx];
        if(fStretchFullResponse) timeFactor *= fTimeScaleFactor;
        double plotTimeFactor = 1.0;
        if(!fStretchFullResponse) plotTimeFactor = f3DCorrectionVec[planeIdx]*fTimeScaleFactor;
        //std::cout << "Time factors " << timeFactor << " " << plotTimeFactor << std::endl;
        
        double timeFactorInv = 1./timeFactor;
        
        const DoubleVec* pResp = &((fSignalShapingVec[ktype][planeIdx]).Response_save());
        
        // recover the field response tool for this plane
        const icarus_tool::IFieldResponse* fieldResponsePtr = fPlaneToResponseMap.at(planeIdx).front().get()->getFieldResponse();
        
        double deltaInputTime = fieldResponsePtr->getBinWidth();
        int    nFieldBins     = fieldResponsePtr->getNumBins();
        
        // more histos
        //std::cout << "HistDone " << view << " " << fHistDoneF[view] << std::endl;
        
        if(!fHistDoneF[planeIdx] && ktype == 0)
        {
            double xLowF = fieldResponsePtr->getLowEdge() * plotTimeFactor;
            double xHighF = xLowF + 0.001*(nFieldBins+1)*deltaInputTime*plotTimeFactor;
            double nBins = nFieldBins*plotTimeFactor;
            
            char histBuf[80];
            sprintf(histBuf, "FullResponse%i", (int)planeIdx);
            TH1D* fullResponse = tfs->make<TH1D>(histBuf, histBuf, nBins, xLowF, xHighF);
            for (size_t i=0; i<nBins; ++i) {
                fullResponse->SetBinContent(i+1, pResp->at(i));
                std::cout << " bin " << i << " full response " <<pResp->at(i) << std::endl;
            }
        }
        
        size_t nticks_input = pResp->size();
        DoubleVec InputTime(nticks_input, 0. );
        for (size_t itime = 0; itime < nticks_input; itime++ ) {
            InputTime[itime] = (1.*itime) * deltaInputTime*timeFactor;
        }
        //std::cout << "Input time vector done" << std::endl;
        
        DoubleVec SamplingResp(nticks, 0. );
        
        size_t SamplingCount = 0;
        
        size_t startJ = 1;
        SamplingResp[0] = (*pResp)[0];
        for ( size_t itime = 1; itime < nticks; itime++ ) {
            size_t low, high;
            for ( size_t jtime = startJ; jtime < nticks_input; jtime++ ) {
                if ( InputTime[jtime] >= SamplingTime[itime] ) {
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
        
        // more histos
        //std::cout << "HistDone " << view << " " << fHistDoneF[view] << std::endl;
        if(!fHistDoneF[planeIdx] && ktype == 0)
        {
            double plotTimeFactor = f3DCorrectionVec[planeIdx]*fTimeScaleFactor;
            double xLowF = fieldResponsePtr->getLowEdge() * plotTimeFactor;
            double xHighF = xLowF + 0.001*(nFieldBins)*deltaInputTime*plotTimeFactor;
            double binWidth = 0.5;
            size_t nBins = (xHighF-xLowF+1)/binWidth;
            
            char histBuf[80];
            sprintf(histBuf, "SampledResponse%i", (int)planeIdx);
            TH1D* sampledResponse = tfs->make<TH1D>(histBuf, histBuf, nBins, xLowF, xHighF);
            for (size_t i=0; i<nBins; ++i) {
                //std::cout << "bin/SamplingResp " << i << " " << SamplingResp[i] << std::endl;
                sampledResponse->SetBinContent(i+1, SamplingResp[i]);
            }
            fHistDoneF[planeIdx] = true;
        }
        
        if(fPrintResponses) {
            size_t printCount = 0;
            int inc = 1;
            //std::cout << "Sampled response (ticks) for view " << view << " wire " << _wr << " nticks " << nticks << std::endl;
            for(size_t i = 0; i<nticks; i+=inc) {
                //std::cout << SamplingResp[i] << " " ;
                if((printCount+1)%10==0) std::cout << std::endl;
                printCount++;
                if (printCount>=100) {inc = 100;}
            }
        }
        
        (fSignalShapingVec[ktype][planeIdx]).AddResponseFunction( SamplingResp, true);
        //std::cout << "Finished with wire " << _wr << ", view " << _vw << std::endl;
    } // loop over views
    
    //std::cout << "Done with field responses" << std::endl;
    return;
}


//-----Give Gain Settings to SimWire-----//jyoti
double util::SignalShapingServiceICARUS::GetASICGain(unsigned int  channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    double gain     = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICGain();
    
    return gain;
}


//-----Give Shaping time to SimWire-----//jyoti
double util::SignalShapingServiceICARUS::GetShapingTime(unsigned int  channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx     = geom->ChannelToWire(channel)[0].Plane;
    double shaping_time = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICShapingTime();

    return shaping_time;
}

double util::SignalShapingServiceICARUS::GetRawNoise(unsigned int const channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    double gain         = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICGain();
    double shaping_time = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICShapingTime();
    int    temp;
    
    if (std::abs(shaping_time - 0.5)<1e-6){
        temp = 0;
    }else if (std::abs(shaping_time - 1.5)<1e-6){
        temp = 1;
    }else if (std::abs(shaping_time - 2.0)<1e-6){
        temp = 2;
    }else{
        temp = 3;
    }
    
    double rawNoise;
    
    auto tempNoise = fNoiseFactVec.at(planeIdx);
    rawNoise = tempNoise.at(temp);
    
    rawNoise *= gain/4.7;
    return rawNoise;
}

double util::SignalShapingServiceICARUS::GetDeconNoise(unsigned int const channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    double shaping_time = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICShapingTime();
    int temp;
    
    if (std::abs(shaping_time - 0.5)<1e-6){
        temp = 0;
    }else if (std::abs(shaping_time - 1.0)<1e-6){
        temp = 1;
    }else if (std::abs(shaping_time - 2.0)<1e-6){
        temp = 2;
    }else{
        temp = 3;
    }
    auto tempNoise = fNoiseFactVec.at(planeIdx);
    double deconNoise = tempNoise.at(temp);
    
    deconNoise = deconNoise /4096.*2000./4.7 *6.241*1000/fDeconNorm;
    return deconNoise;
}

int util::SignalShapingServiceICARUS::FieldResponseTOffset(unsigned int const channel, size_t ktype) const
{
    art::ServiceHandle<geo::Geometry> geom;
    
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    double time_offset(0.);
    
    try
    {
        time_offset = fPlaneToResponseMap.at(planeIdx).front()->getFieldResponse()->getTOffset();
    }
    catch (...)
    {
        throw cet::exception(__FUNCTION__) << "Invalid plane ... " << planeIdx << std::endl;
    }

    auto tpc_clock = lar::providerFrom<detinfo::DetectorClocksService>()->TPCClock();
    return tpc_clock.Ticks(time_offset/1.e3);
}

namespace util {
    
    DEFINE_ART_SERVICE(SignalShapingServiceICARUS)
    
}
