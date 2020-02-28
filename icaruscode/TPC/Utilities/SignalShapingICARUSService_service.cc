////////////////////////////////////////////////////////////////////////
/// \file   fSignalShapingICARUSService_service.cc
/// \author H. Greenlee
/// Modified by X. Qian 1/6/2015
/// if histogram is used, inialize
/// Response_Offset, Response_Sampling, FieldBins from histogram
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "SignalShapingICARUSService_service.h"
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
#include "tools/IResponse.h"
#include "tools/IFieldResponse.h"
#include "tools/IElectronicsResponse.h"
#include "tools/IFilter.h"

#include "tbb/spin_mutex.h"

#include <fstream>

namespace icarusutil
{
    
tbb::spin_mutex signalShapingSpinMutex;

//----------------------------------------------------------------------
// Constructor.
SignalShapingICARUSService::SignalShapingICARUSService(const fhicl::ParameterSet& pset,
                                                             art::ActivityRegistry& /* reg */)
: fInit(false)
{
    reconfigure(pset);
}

//----------------------------------------------------------------------
// Destructor
SignalShapingICARUSService::~SignalShapingICARUSService()
{
    return;
}

//----------------------------------------------------------------------
// Reconfigure method.
void SignalShapingICARUSService::reconfigure(const fhicl::ParameterSet& pset)
{
    mf::LogInfo("fSignalShapingICARUSService") << " reconfiguring setup " ;
    
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
    
    fPlaneForNormalization  = pset.get<size_t>(    "PlaneForNormalization");
    fPrintResponses         = pset.get<bool>(      "PrintResponses"       );
    fDeconNorm              = pset.get<double>(    "DeconNorm"            );
    fInitialFFTSize         = pset.get<size_t>(    "InitialFFTSize"       );
    fNoiseFactVec           = pset.get<DoubleVec2>("NoiseFactVec"         );
    fStoreHistograms        = pset.get<bool>(      "StoreHistograms"      );
    
    fInit = false;
    
    return;
}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
//const SignalShapingICARUS& SignalShapingICARUSService::SignalShaping(size_t channel) const
//{
//    if(!fInit) init();
//    
//    art::ServiceHandle<geo::Geometry> geom;
//    
//        //use channel number to set some useful numbers
//    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
//    
//    const SignalShapingICARUS& signalShaping = fPlaneToResponseMap.at(planeIdx).front()->getSignalShapingICARUS();
//
//    return signalShaping;
//}

const icarus_tool::IResponse& SignalShapingICARUSService::GetResponse(size_t channel) const
{
    if (!fInit) init();
    
    art::ServiceHandle<geo::Geometry> geom;
    
        //use channel number to set some useful numbers
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;

    return *fPlaneToResponseMap.at(planeIdx).front();
}


//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void SignalShapingICARUSService::init()
{
    // First get a lock to make sure we are clear to run
    tbb::spin_mutex::scoped_lock lock(signalShapingSpinMutex);

    if(!fInit)
    {
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
        if (fStoreHistograms && !fStoreHistograms)
        {
            art::ServiceHandle<art::TFileService> tfs;
            
            // Make sure we are at the top level
            tfs->file().cd();
            
            // Make a directory for these histograms
            art::TFileDirectory dir = tfs->mkdir("SignalShaping");
            
            // Loop through response tools first
            for(const auto& response: fPlaneToResponseMap) response.second.front().get()->outputHistograms(dir);

        }

        fInit = true;
    }
    
    return;
}

void SignalShapingICARUSService::SetDecon(size_t fftsize, size_t channel)
{
    art::ServiceHandle<geo::Geometry> geo;
    
    if (!fInit) init();
    
    // Assume we need to reset the kernels
    double integral = fPlaneToResponseMap.at(fPlaneForNormalization).front().get()->getFieldResponse()->getIntegral();
    double weight   = 1. / integral;
    
    for(size_t planeIdx = 0; planeIdx < geo->Nplanes(); planeIdx++)
    {
        fPlaneToResponseMap.at(planeIdx).front()->setResponse(weight);
    }
}

//-----Give Gain Settings to SimWire-----
double SignalShapingICARUSService::GetASICGain(unsigned int  channel) const
{
    static const double fcToElectrons(6241.50975);
    
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx = geom->ChannelToWire(channel)[0].Plane;
    double gain     = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getFCperADCMicroS() * fcToElectrons;
    
    return gain;
}

//-----Give Shaping time to SimWire-----
double SignalShapingICARUSService::GetShapingTime(unsigned int  channel) const
{
    art::ServiceHandle<geo::Geometry> geom;
    size_t planeIdx     = geom->ChannelToWire(channel)[0].Plane;
    double shaping_time = fPlaneToResponseMap.at(planeIdx).front()->getElectronicsResponse()->getASICShapingTime();

    return shaping_time;
}

double SignalShapingICARUSService::GetRawNoise(unsigned int const channel) const
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

double SignalShapingICARUSService::GetDeconNoise(unsigned int const channel) const
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

int SignalShapingICARUSService::FieldResponseTOffset(unsigned int const channel) const
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
}

namespace icarusutil {
    
    DEFINE_ART_SERVICE(SignalShapingICARUSService)
    
}
