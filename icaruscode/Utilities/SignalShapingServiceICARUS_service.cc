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

#include "art/Utilities/make_tool.h"
#include "tools/IFieldResponse.h"

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
    // Implement the tools for handling the responses
    const fhicl::ParameterSet& fieldResponseTools = pset.get<fhicl::ParameterSet>("FieldResponseTools");
    
    for(const std::string& responseTool : fieldResponseTools.get_pset_names())
    {
        const fhicl::ParameterSet& responseToolParamSet = fieldResponseTools.get<fhicl::ParameterSet>(responseTool);
        size_t                     planeIdx             = responseToolParamSet.get<size_t>("Plane");
        
        fPlaneToFieldResponseVec[planeIdx].push_back(art::make_tool<icarus_tool::IFieldResponse>(responseToolParamSet));
        
        std::cout << "Field response set up for plane " << fPlaneToFieldResponseVec[planeIdx].back()->getPlane() << std::endl;
        std::cout << "   --> Bin width: " << fPlaneToFieldResponseVec[planeIdx].back()->getBinWidth() << std::endl;
        std::cout << "   --> T offset: " << fPlaneToFieldResponseVec[planeIdx].back()->getTOffset() << " # bins: " << fPlaneToFieldResponseVec[planeIdx].back()->getNumBins() << std::endl;
    }
    
    // add a comment here
    art::ServiceHandle<geo::Geometry> geo;
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    // Reset initialization flag.
    
    fInit = false;
    
    std::cout << " before reading fhicl " << std::endl;
    fASICGainInMVPerFC    = pset.get< DoubleVec >("ASICGainInMVPerFC");
    
    fViewForNormalization = pset.get<size_t>("PlaneForNormalization");
    fPrintResponses   = pset.get<bool>("PrintResponses");
    
    std::cout << " after reading fhicl " << std::endl;
    
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
    
    std::cout << " resized shaping vec " << std::endl;
    
    fFieldResponseVec.resize(2);
    for(auto& kset : fFieldResponseVec) kset.resize(geo->Nplanes());
    
    std::cout << " before decon " << std::endl;
    // Fetch fcl parameters.
    fDeconNorm = pset.get<double>("DeconNorm");
    std:: cout << " before ADC " << std::endl;

    fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");
    std:: cout << " after ADC " << std::endl;

    fDefaultDriftVelocity = pset.get< DoubleVec >("DefaultDriftVelocity");
    std:: cout << " after drift velocity " << std::endl;
    
    std::cout << " after field response " << std::endl;
    
    fCalibResponseTOffset = pset.get< DoubleVec >("CalibResponseTOffset");
    std::cout << "CalibResponseTOffsets: ";
    for(auto& x : fCalibResponseTOffset) { std::cout << x << " "; }
    std::cout << std::endl;

    fNoiseFactVec =  pset.get<DoubleVec2>("NoiseFactVec");
    
    f3DCorrectionVec = pset.get<DoubleVec>("Drift3DCorrVec");
    
    fFieldRespAmpVec = pset.get<DoubleVec>("FieldRespAmpVec");
    
    fShapeTimeConst = pset.get<DoubleVec >("ShapeTimeConst");
    fDeconvPol = pset.get<std::vector<int> >("DeconvPol");
    
    fGetFilterFromHisto= pset.get<bool>("GetFilterFromHisto");
    
    // Construct parameterized collection filter function.
    
    fFilterWidthCorrectionFactor = pset.get<DoubleVec>("FilterWidthCorrectionFactor", DoubleVec() = {1.0, 1.0, 1.0});

    if(!fGetFilterFromHisto) {
        
        fFilterFuncVec.resize(geo->Nplanes());
        std::cout <<"Getting Filters from .fcl file" << std::endl;
        mf::LogInfo("SignalShapingServiceICARUS") << "Getting Filters from .fcl file" ;
        
        fFilterParamsVec = pset.get< DoubleVec2 >("FilterParamsVec");
        fFilterFuncVec = pset.get<std::vector<std::string> > ("FilterFuncVec");
        
        fFilterTF1Vec.resize(geo->Nplanes());
        std::cout << " before loop " << geo->Nplanes() << std::endl;
        for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
        {
            std::cout << " view " << planeIdx << std::endl;
            std::string name = Form("Filter_vw%02i_wr00", (int)planeIdx);
            std::cout << " filter size " <<fFilterParamsVec[planeIdx].size() << std::endl;
            fFilterTF1Vec[planeIdx] = new TF1(name.c_str(), fFilterFuncVec[planeIdx].c_str() );
            for(size_t idx = 0; idx < fFilterParamsVec[planeIdx].size(); idx++)
                fFilterTF1Vec[planeIdx]->SetParameter(idx, fFilterParamsVec[planeIdx][idx]);
        }
    }
    else
    {
        std::string histoname = pset.get<std::string>("FilterHistoName");
        mf::LogInfo("SignalShapingServiceICARUS") << " using filter from .root file " ;
        
        // constructor decides if initialized value is a path or an environment variable
        std::string fname;
        cet::search_path sp("FW_SEARCH_PATH");
        sp.find_file(pset.get<std::string>("FilterFunctionFname"), fname);
        
        TFile * in=new TFile(fname.c_str(),"READ");
        for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
        {
            std::string name = Form("%s_vw%02i", histoname.c_str(), (int)planeIdx);
            fFilterHistVec[planeIdx] = (TH1D *)in->Get(name.c_str());
        }
        
        in->Close();
        delete in;
    }

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
    
    //Adding calibrated field response at 70kV
    fUseCalibratedResponses = pset.get<bool>("UseCalibratedResponses");
    
    mf::LogInfo("SignalShapingServiceICARUS") << " using the field response provided from a .root file " ;
    
    // constructor decides if initialized value is a path or an environment variable
    fDefaultEField                 = pset.get<double>("DefaultEField");
    fDefaultTemperature            = pset.get<double>("DefaultTemperature");
    
    fTimeScaleParams               = pset.get<DoubleVec>("TimeScaleParams");
    fStretchFullResponse           = pset.get<bool>("StretchFullResponse");
    
    // calculate the time scale factor for this job
    if(!fUseCalibratedResponses) SetTimeScaleFactor();
    
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
    if(!fInit)
        init();
    
    // Figure out plane type.
    
    art::ServiceHandle<geo::Geometry> geom;
    //geo::SigType_t sigtype = geom->SignalType(channel);
    
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
    if(!fInit) {
        fInit = true;
        
        art::ServiceHandle<geo::Geometry> geo;
        
        // Do ICARUS-specific configuration of SignalShaping by providing
        // ICARUS response and filter functions.
        
        // re-initialize the FFT service for the request size
        art::ServiceHandle<util::LArFFT> fFFT;
        std::string options = fFFT->FFTOptions();
        size_t fitbins = fFFT->FFTFitBins();
        //size_t fftsize = fFFT->FFTSize();
        int fftsize = (int) fFFT->FFTSize();
        std::cout << " fftsize " << fftsize << std::endl;
        // Calculate field and electronics response functions.
        
        std::string kset[2] = { "Convolution ", "Deconvolution "};
        
        for(size_t ktype=0;ktype<2;++ktype) {
            std::cout << std::endl << kset[ktype] << " setting functions:" << std::endl;
            int fftsize2 = (int) fFFT->FFTSize();
            std::cout << " fftsize2 " << fftsize2 << std::endl;
            // call this first, so that the binning will be known to SetElectResponse
            SetFieldResponse(ktype);
            
            std::cout << "Input field responses" << std::endl;
            
            for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
            {
                // Make sure the FFT can handle this
                int nFieldBins = fPlaneToFieldResponseVec.at(planeIdx).front()->getNumBins();
                
                if (nFieldBins*4>fftsize)
                    fFFT->ReinitializeFFT( (size_t)nFieldBins*4, options, fitbins);
                
                SetElectResponse(ktype,planeIdx, fShapeTimeConst.at(planeIdx),fASICGainInMVPerFC.at(planeIdx));
                //Electronic response
                std::cout << " ktype " << ktype << " Electonic response " << fElectResponse[ktype].size() << " bins" << std::endl;

                if(fPrintResponses)
                {
                    std::cout << "Input field response for view " << planeIdx << ", " << (fFieldResponseVec[ktype][planeIdx]).size() << " bins" << std::endl;
                    for(size_t i = 0; i<(fFieldResponseVec[ktype][planeIdx]).size(); ++i)
                    {
                        std::cout << fFieldResponseVec[ktype][planeIdx][i] << " " ;
                        if((i+1)%10==0) std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
                std::cout << " adding response for ktype " << ktype << " plane " << planeIdx << std::endl;
                (fSignalShapingVec[ktype][planeIdx]).AddResponseFunction(fFieldResponseVec[ktype][planeIdx]);
                std::cout << " adding response for ktype " << ktype << " view " << planeIdx << std::endl;
                (fSignalShapingVec[ktype][planeIdx]).AddResponseFunction(fElectResponse[ktype]);
                (fSignalShapingVec[ktype][planeIdx]).save_response();
                (fSignalShapingVec[ktype][planeIdx]).set_normflag(false);
            }
            // see if we get the same toffsets
            SetResponseSampling(ktype);
            
            // Currently we only have fine binning "fFieldBinWidth"
            // for the field and electronic responses.
            // Now we are sampling the convoluted field-electronic response
            // with the nominal sampling.
            // We may consider to do the same for the filters as well.
            if ((int)fftsize!=fFFT->FFTSize()){
                std::string options = fFFT->FFTOptions();
                int fitbins = fFFT->FFTFitBins();
                fFFT->ReinitializeFFT( (size_t)fftsize, options, fitbins);
            }
            
            
            // Calculate filter functions.
            if(ktype == 0) SetFilters();
            
            // Configure deconvolution kernels.
            for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
            {
                (fSignalShapingVec[ktype][planeIdx]).AddFilterFunction(fFilterVec[planeIdx]);
                (fSignalShapingVec[ktype][planeIdx]).SetDeconvKernelPolarity( fDeconvPol.at(planeIdx));
                (fSignalShapingVec[ktype][planeIdx]).CalculateDeconvKernel();
            }
        }
    }
}

void util::SignalShapingServiceICARUS::SetDecon(size_t fftsize, size_t channel)
{
    art::ServiceHandle<geo::Geometry> geo;
    
    //std::cout << "enter SetDecon, init flag "  << fInit <<  " fftsize " << fftsize << " channel " << channel << std::endl;
    
    init();
    
    
    art::ServiceHandle<util::LArFFT> fFFT;
    
    // streamline this method:
    // if the deconvolution kernel is already appropriate for the datasize (aka fftsize) do nothing
    // otherwise, set it to the appropriate size
    // do this test for *every* ss
    // But it will in general only happen once per run!
    
    bool setDecon = false;
    
    size_t FFTSize = fFFT->FFTSize();
    if (fftsize>FFTSize||fftsize<=FFTSize/2){
        std::string options = fFFT->FFTOptions();
        int fitbins = fFFT->FFTFitBins();
        fFFT->ReinitializeFFT( (size_t)fftsize, options, fitbins);
        setDecon = true;
    }
    
    if(!setDecon) return;
    
    size_t ktype = 1;
    
    for (size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
    {
        (fSignalShapingVec[ktype][planeIdx]).Reset();
    }
    
    //std::cout << "Xin2 " << std::endl;
    // Calculate filter functions.
    //std::cout << "set the filters" << std::endl;
    SetFilters();
    // Configure deconvolution kernels.
    //std::cout << "Xin3 " << std::endl;
    //std::cout << "FInish the SS" << std::endl;
    
    for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
    {
        (fSignalShapingVec[ktype][planeIdx]).AddFilterFunction(fFilterVec[planeIdx]);
        (fSignalShapingVec[ktype][planeIdx]).SetDeconvKernelPolarity( fDeconvPol.at(planeIdx));
        (fSignalShapingVec[ktype][planeIdx]).CalculateDeconvKernel();
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
    
    char buff0[80]; //buff1[80];
    
    // Ticks in nanosecond
    // Calculate the normalization of the collection plane
    double integral = fPlaneToFieldResponseVec.at(fViewForNormalization).front().get()->getIntegral();
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
        // simplify the code
        DoubleVec* responsePtr = &fFieldResponseVec[ktype][planeIdx];
        std::cout << " response ptr " << std::endl;
        
        const icarus_tool::IFieldResponse* fieldResponsePtr = fPlaneToFieldResponseVec.at(planeIdx).front().get();
        
        size_t nBins = fieldResponsePtr->getNumBins();
        std::cout << " nBins " << nBins << " timeFactor " << timeFactor << std::endl;
        size_t nResponseBins = nBins*timeFactor;
        responsePtr->resize(nResponseBins);
        //double x0 = histPtr->GetBinCenter(1);
        //double xf = histPtr->GetBinCenter(nBins);
        double x0 = fieldResponsePtr->getBinCenter(1);
        double xf = fieldResponsePtr->getBinCenter(nBins);
        double deltaX = (xf - x0)/(nBins-1);
        std::cout << "lims " << x0 << " " << xf << " " << deltaX << std::endl;
        
        for(size_t bin = 1; bin <= nResponseBins; bin++)
        {
            double xVal = x0 + deltaX*(bin-1)/timeFactor;
            if(bin==1) std::cout << "1st bin " << x0 << " " << xVal << std::endl;
            double yVal = fieldResponsePtr->interpolate(xVal);
            responsePtr->at(bin-1) = yVal;
            responsePtr->at(bin-1) *= fFieldRespAmpVec[planeIdx]*weight;
        }
        
        std::cout << " after yval " << std::endl;
        
        // fill some histos
        sprintf(buff0, "hRawResp_%i_%i", (int)ktype, (int)planeIdx);
        TH1D* rawResponse = tfs->make<TH1D>(buff0, buff0, nBins, x0-0.5*deltaX, xf+0.5*deltaX);
        sprintf(buff0, "hStretchedResp_%i_%i", (int)ktype, (int)planeIdx);
        double x0S = timeFactor*x0 - 0.5*deltaX/timeFactor;
        double xfS = timeFactor*xf + 0.5*deltaX/timeFactor;
        std::cout << "title " << buff0 << std::endl;
        std::cout << " NBINS " << nBins << std::endl;
        std::cout << " NRESPONSEBINS " << nResponseBins << std::endl;
        std::cout << " x0S " << x0S << " xfS " << xfS <<std::endl;
        TH1D* stretchedResponse = tfs->make<TH1D>(buff0, buff0, nResponseBins, x0S, xfS);
        std::cout << " NBINS " << nBins << std::endl;

        for(size_t i=0;i<nBins; ++i) {
            rawResponse->SetBinContent(i, fieldResponsePtr->getBinContent(i));
            std::cout << "bin " << i <<  " xVal " << rawResponse->GetBinCenter(i) << " response " << rawResponse->GetBinContent(i) << std::endl;
        }
        for(size_t i=0;i<nResponseBins; ++i) {
            stretchedResponse->SetBinContent(i+1, responsePtr->at(i));
            std::cout << "vbin " << i <<  " xVal " << stretchedResponse->GetBinCenter(i) << " response " << stretchedResponse->GetBinContent(i) << std::endl;
        }
    }
    
    std::cout << " end SetFieldResponse " << std::endl;
    return;
}


//----------------------------------------------------------------------
// Calculate ICARUS field response.
void util::SignalShapingServiceICARUS::SetElectResponse(size_t ktype, size_t planeIdx, double shapingtime, double gain)
{
    // Get services.
    
    art::ServiceHandle<util::LArFFT> fft;
    
    LOG_DEBUG("SignalShapingICARUS") << "Setting ICARUS electronics response function...";
    
    size_t nticks = fft->FFTSize();
    DoubleVec time(nticks,0.);
    
    fElectResponse.resize(2);
    for(auto& resp : fElectResponse) {
        std::cout << " resizing elect response " << nticks << std::endl;
        resp.resize(nticks, 0.);
    }
    
    // recover the field response tool for this plane
    const icarus_tool::IFieldResponse* fieldResponsePtr = fPlaneToFieldResponseVec.at(planeIdx).front().get();
    
    double binWidth = fieldResponsePtr->getBinWidth();
    
    //Gain and shaping time variables from fcl file:
    //double Ao = 1.0;//Gain
    double To = shapingtime;  //peaking time
    std::cout << " electronic shaping time " << shapingtime << std::endl;
    // this is actually sampling time, in ns
   //  mf::LogInfo("SignalShapingICARUS") << "Check sampling intervals: "
     //                                 << fSampleRate << " ns"
    //                                  << "Check number of samples: " << fNTicks;
    
    // The following sets the ICARUS electronics response function in
    // time-space. Function comes from BNL SPICE simulation of ICARUS
    // electronics. SPICE gives the electronics transfer function in
    // frequency-space. The inverse laplace transform of that function
    // (in time-space) was calculated in Mathematica and is what is being
    // used below. Parameters Ao and To are cumulative gain/timing parameters
    // from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain.
    // They have been adjusted to make the SPICE simulation to match the
    // actual electronics response. Default params are Ao=1.4, To=0.5us.
    
    
    // For the cold electronics,  the gain (i.e. 4.7 mV/fC) represents the peak
    // height. The shaping time will not affect the peak height, but make the
    // peak broader
    
    double max = 0;
    
    for(size_t i=0; i<=nticks;++i) {
        time[i] = (1.*i)* binWidth*1.e-3 ;
    }
    int i = 0;
    for(auto& element :fElectResponse[ktype]) {
        //convert time to microseconds, to match fElectResponse[i] definition
        element = time[i]/To*exp(-time[i]/To);
        if(element > max) max = element;
        i++;
    }// end loop over time buckets
    
    LOG_DEBUG("SignalShapingICARUS") << " Done.";
    
    // normalize fElectResponse[i], before the convolution
    // Put in overall normalization in a pedantic way:
    // first put in the pulse area per eleectron at the lowest gain setting,
    // then normalize by the actual ASIC gain setting used.
    // This code is executed only during initialization of service,
    // so don't worry about code inefficiencies here.
    double last_integral=0;
    double last_max=0;
    
    //Normalization are the following
    // Peak is firstly normalized to 1
    // thus we expect peak to be 1 * 9390 (fADCPerPCtAtLowestAsicGain) * 1.602e-7 * (1 fC) = 9.39 ADC
    // At 4.7 mV/fC, the ADC value should be 4.7 (mV/fC) * 2 (ADC/mV) ~ 9.4 ADC/fC
    // so the normalization are consistent
    
    
    
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    for(auto& element : fElectResponse[ktype]){
        element /= (max);
        element *= gain / 6.5 ;
        
       element *= fADCPerPCAtLowestASICGain * 1.60217657e-7;
       // element *= gain / 4.7;
        
        
        if(element > last_max) last_max = element;
        last_integral += element * binWidth / detprop->SamplingRate();
    }
    return;
}


//----------------------------------------------------------------------
// Calculate ICARUS filter functions.
void util::SignalShapingServiceICARUS::SetFilters()
{
    // Get services.
    art::ServiceHandle<geo::Geometry> geo;
    
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<util::LArFFT> fft;
    
    double ts = detprop->SamplingRate();
    size_t nFFT2 = fft->FFTSize() / 2;
    
    // Calculate collection filter.
    
    fFilterVec.resize(geo->Nplanes());
//    for(auto& filter : fFilterVec) {
//        filter.resize(nFFT2+1);
//    }
    
    if(!fGetFilterFromHisto)
    {
        size_t planeIdx = 0;
        
        for(auto& func : fFilterTF1Vec)
        {
            func->SetRange(0, double(nFFT2));
            size_t count = 0;
            
            // now to scale the filter function!
            // only scale params 1,2 &3
            
            double timeFactor = fTimeScaleFactor*f3DCorrectionVec[int(planeIdx)]*fFilterWidthCorrectionFactor[int(planeIdx)];
            for(size_t i=1;i<4;++i) {
                func->SetParameter(i, fFilterParamsVec[int(planeIdx)][i]/timeFactor);
            }
            
            for(size_t bin = 0; bin <= nFFT2; bin++)
            {
                //std::cout << "checking TF1 generation " << _bn << " " <<nFFT2 << std::endl;
                double freq = 500.*bin/(ts*nFFT2);
                double f = func->Eval(freq);
                if(f!=0.0) count++;
                //fFilterVec[int(fViewIndex[_vw])][_bn] = TComplex(f, 0.);
                fFilterVec[int(planeIdx)].push_back(TComplex(f, 0.));
            }
            //std::cout << count << " non-zero bins out of " << nFFT2 << std::endl;
            planeIdx++;
        }
    } else{
        
        size_t planeIdx = 0;
        for(auto hist : fFilterHistVec) {
            for(size_t bin = 1; bin <= nFFT2+1; bin++)
            {
                double f = hist->GetBinContent(bin);
                //fFilterVec[int(fViewIndex[_vw])][_bn-1] = TComplex(f, 0.);
                fFilterVec[int(planeIdx)].push_back(TComplex(f, 0.));
            }
            planeIdx++;
        }
    }
    
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
    
    size_t plane0, plane1;
    if(mode==0) {
        plane0 = 0;
        plane1 = geo->Nplanes();
    } else {
        size_t plane = geo->ChannelToWire(channel)[0].Plane;
        plane0 = plane;
        plane1 = std::min(size_t(geo->Nplanes()),plane+1);
    }
    
    //std::cout << "view0/1 " << view0 << " " << view1 << std::endl;
    
    size_t nticks = fft->FFTSize();
    DoubleVec SamplingTime( nticks, 0. );
    //std::cout << "nticks = " << nticks << std::endl;
    for ( size_t itime = 0; itime < nticks; itime++ ) {
        SamplingTime[itime] = (1.*itime) * detprop->SamplingRate();
    }
    // Sampling
    
    //std::cout << "sampling view " << view  << " ktype/config/channel " << ktype << " " << config << " " << channel << std::endl;
    
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
        const icarus_tool::IFieldResponse* fieldResponsePtr = fPlaneToFieldResponseVec.at(planeIdx).front().get();
        
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
    double gain = fASICGainInMVPerFC.at(planeIdx);
    
    return gain;
}


//-----Give Shaping time to SimWire-----//jyoti
double util::SignalShapingServiceICARUS::GetShapingTime(unsigned int  channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    double shaping_time = fShapeTimeConst.at(planeIdx);

    return shaping_time;
}

double util::SignalShapingServiceICARUS::GetRawNoise(unsigned int const channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    double shapingtime = fShapeTimeConst.at(planeIdx);
    double gain = fASICGainInMVPerFC.at(planeIdx);
    int temp;
    if (std::abs(shapingtime - 0.5)<1e-6){
        temp = 0;
    }else if (std::abs(shapingtime - 1.5)<1e-6){
        temp = 1;
    }else if (std::abs(shapingtime - 2.0)<1e-6){
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
    
    double shapingtime = fShapeTimeConst.at(planeIdx);
    int temp;
    if (std::abs(shapingtime - 0.5)<1e-6){
        temp = 0;
    }else if (std::abs(shapingtime - 1.0)<1e-6){
        temp = 1;
    }else if (std::abs(shapingtime - 2.0)<1e-6){
        temp = 2;
    }else{
        temp = 3;
    }
    auto tempNoise = fNoiseFactVec.at(planeIdx);
    double deconNoise = tempNoise.at(temp);
    
    deconNoise = deconNoise /4096.*2000./4.7 *6.241*1000/fDeconNorm;
    return deconNoise;
}

int util::SignalShapingServiceICARUS::FieldResponseTOffset(unsigned int const channel, size_t
                                                           ktype) const
{
    art::ServiceHandle<geo::Geometry> geom;
    
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    double time_offset(0.);
    
    try
    {
        time_offset = fPlaneToFieldResponseVec.at(planeIdx).front()->getTOffset();
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
