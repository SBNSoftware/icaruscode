////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceICARUS_service.cc
/// \author H. Greenlee
/// Modified by X. Qian 1/6/2015
/// if histogram is used, inialize
/// Response_Offset, Response_Sampling, FieldBins from histogram
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "SignalShaperServiceICARUS.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "TFile.h"

// LArSoft include
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

#include "art/Utilities/make_tool.h"
#include "tools/shprIResponse.h"
#include "tools/IFieldResponse.h"
#include "tools/IElectronicsResponse.h"
#include "tools/shprIFilter.h"

#include <fstream>

//----------------------------------------------------------------------
// Constructor.
util::SignalShaperServiceICARUS::SignalShaperServiceICARUS(const fhicl::ParameterSet& pset,
                                                             art::ActivityRegistry& /* reg */)
: fInit(false)
, fInitialFFTSize   (pset.get<int>("InitialFFTSize"))
, fFFTSize   (setFFTSize(pset.get< int >("FFTSize", 0),fInitialFFTSize))
, fFFTOption (pset.get< std::string >("FFTOption"))
, fFFTFitBins(pset.get< int         >("FFTFitBins"))
{
    reconfigure(pset);
}

//----------------------------------------------------------------------
// Destructor
util::SignalShaperServiceICARUS::~SignalShaperServiceICARUS()
{
    return;
}

//----------------------------------------------------------------------
// Reconfigure method.
void util::SignalShaperServiceICARUS::reconfigure(const fhicl::ParameterSet& pset)
{
    mf::LogInfo("SignalShaperServiceICARUS") << " reconfiguring setup " ;
    
    // If called again, then we need to clear out the existing tools...
    fPlaneToResponseMap.clear();
    
    // Implement the tools for handling the responses
    fhicl::ParameterSet responseTools = pset.get<fhicl::ParameterSet>("ResponseTools");
    
    for(std::string& responseTool : responseTools.get_pset_names())
    {
        fhicl::ParameterSet responseToolParamSet = responseTools.get<fhicl::ParameterSet>(responseTool);
	    responseToolParamSet.put_or_replace<int>("FFTSize",fFFTSize);
	    responseToolParamSet.put_or_replace<std::string>("FFTOption",fFFTOption);
        size_t                     planeIdx             = responseToolParamSet.get<size_t>("Plane");
        
        fPlaneToResponseMap[planeIdx].push_back(art::make_tool<icarus_tool::shprIResponse>(responseToolParamSet));
    }
    
    fPlaneForNormalization  = pset.get<size_t>(    "PlaneForNormalization");
    fPrintResponses         = pset.get<bool>(      "PrintResponses"       );
    fDeconNorm              = pset.get<double>(    "DeconNorm"            );
    //fInitialFFTSize         = pset.get<size_t>(    "InitialFFTSize"       );
    fNoiseFactVec           = pset.get<DoubleVec2>("NoiseFactVec"         );
    fStoreHistograms        = pset.get<bool>(      "StoreHistograms"      );
    
    fInit = false;
    init();	// mwang added
    
    return;
}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaper& util::SignalShaperServiceICARUS::SignalShaper(size_t channel) const
{
    if(!fInit) init();
    
    art::ServiceHandle<geo::Geometry> geom;
    
        //use channel number to set some useful numbers
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    const util::SignalShaper& SignalShaper = fPlaneToResponseMap.at(planeIdx).front()->getSignalShaper();

    return SignalShaper;
}

//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShaperServiceICARUS::init()
{
    if(!fInit)
    {
        fInit = true;
        
        // Do ICARUS-specific configuration of SignalShaping by providing
        // ICARUS response and filter functions.
        
        art::ServiceHandle<geo::Geometry> geo;

        // Get the normalization from the field response for the collection plane
        double integral = fPlaneToResponseMap.at(fPlaneForNormalization).front().get()->getFieldResponse()->getIntegral();
        double weight   = 1. / integral;
        
        for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
        {
            fPlaneToResponseMap[planeIdx].front().get()->setResponse(weight);                
        }
        
        // Check to see if we want histogram output
        if (fStoreHistograms)
        {
            art::ServiceHandle<art::TFileService> tfs;
            
            // Make sure we are at the top level
            tfs->file().cd();
            
            // Make a directory for these histograms
            art::TFileDirectory dir = tfs->mkdir("SignalShaper");
            
            // Loop through response tools first
            for(const auto& response: fPlaneToResponseMap) response.second.front().get()->outputHistograms(dir);

        }
    }
    
    return;
}

void util::SignalShaperServiceICARUS::SetDecon(int fftsize, size_t channel)
{
    art::ServiceHandle<geo::Geometry> geo;
    
    if (!fInit) init();
    
    // streamline this method:
    // if the deconvolution kernel is already appropriate for the datasize (aka fftsize) do nothing
    // otherwise, set it to the appropriate size
    // do this test for *every* ss
    // But it will in general only happen once per run!
    
    bool setDecon = false;
    
    if (int(fftsize)>fFFTSize||int(fftsize)<=fFFTSize/2)
    {
	fFFTSize = int(fftsize);
        setDecon = true;
    }
    
    if(!setDecon) return;
    
    // Assume we need to reset the kernels
    double integral = fPlaneToResponseMap.at(fPlaneForNormalization).front().get()->getFieldResponse()->getIntegral();
    double weight   = 1. / integral;
    
    for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
    {
        fPlaneToResponseMap.at(planeIdx).front()->setResponse(weight);
    }
}

//-----Give Gain Settings to SimWire-----
double util::SignalShaperServiceICARUS::GetASICGain(unsigned int  channel) const
{
    static const double fcToElectrons(6241.50975);
    
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    double gain     = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getFCperADCMicroS() * fcToElectrons;
    
    return gain;
}

//-----Give Shaping time to SimWire-----
double util::SignalShaperServiceICARUS::GetShapingTime(unsigned int  channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx     = geom->ChannelToWire(channel)[0].Plane;
    double shaping_time = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICShapingTime();

    return shaping_time;
}

double util::SignalShaperServiceICARUS::GetRawNoise(unsigned int const channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    double gain         = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getFCperADCMicroS();
    double shaping_time = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICShapingTime();
    int    temp;
    
    if (std::abs(shaping_time - 0.6)<1e-6){
        temp = 0;
    }else if (std::abs(shaping_time - 1.3)<1e-6){
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

double util::SignalShaperServiceICARUS::GetDeconNoise(unsigned int const channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    double shaping_time = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICShapingTime();
    int temp;
    
    if (std::abs(shaping_time - 0.6)<1e-6){
        temp = 0;
    }else if (std::abs(shaping_time - 1.3)<1e-6){
        temp = 1;
    }else if (std::abs(shaping_time - 2.0)<1.e-6){
        temp = 2;
    }else{
        temp = 3;
    }
    auto   tempNoise  = fNoiseFactVec.at(planeIdx);
    double deconNoise = tempNoise.at(temp);
    
    //deconNoise = deconNoise /4096.*2000./4.7 *6.241*1000/fDeconNorm; <== I don't know where these numbers come from...

    return deconNoise;
}

// mwang added, returns deconvolution kernel
//const util::SignalShaper& util::SignalShaperServiceICARUS::SignalShaper(size_t channel) const
const std::vector<std::complex<double>>& util::SignalShaperServiceICARUS::GetDeconvKernel(unsigned int const channel) const
{
    if(!fInit) init();
    
    art::ServiceHandle<geo::Geometry> geom;
    
        //use channel number to set some useful numbers
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    
    const util::SignalShaper& SignalShaper = fPlaneToResponseMap.at(planeIdx).front()->getSignalShaper();

    return SignalShaper.DeconvKernel();
}

int util::SignalShaperServiceICARUS::FieldResponseTOffset(unsigned int const channel) const
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
    
    return tpc_clock.Ticks(time_offset / 1.e3);
}

int util::SignalShaperServiceICARUS::setFFTSize(int fftsize, int initfftsize)
{
    if(fftsize <= 0){ 
      fftsize = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider()->ReadOutWindowSize();
    }
    int i;
    for(i = 1; i < fftsize; i *= 2){ }   

    //if (initfftsize / i >= 2 || i / initfftsize >= 2)
    //{
    //  int newFFTSize = 64;
    //  while(newFFTSize < initfftsize) newFFTSize *= 2;
    //  i = newFFTSize;
    //}

    return i;
}

namespace util {
    
    DEFINE_ART_SERVICE(SignalShaperServiceICARUS)
    
}
