/**
 *  @file   FakeParticle_tool.cc
 *
 *  @brief  This tool creates a fake particle and overlays on input data
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ICARUS package includes
#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"
#include "icaruscode/Decode/DecoderTools/IFakeParticle.h"

#include "icarus_signal_processing/ICARUSFFT.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  FakeParticle class definiton
 */
class FakeParticle : virtual public IFakeParticle
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit FakeParticle(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~FakeParticle();

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) override;

    /**
     *  @brief Creates a fake particle and overlays on the input fragment
     *
     *  @param waveforms  The waveform container to place fake particle on
     */
    virtual void overlayFakeParticle(ArrayFloat& waveforms) override; 

private:

    // fhicl variables
    std::vector<size_t>                      fWireEndPoints;         //< Desired wire endpoints for our particle
    size_t                                   fStartTick;             //< The tick for the start point of our particle
    float                                    fStartAngle;            //< Angle (in degrees) for the trajectory
    int                                      fNumElectronsPerMM;     //< The number of electrons/mm to deposit
    size_t                                   fPlaneToSimulate;       //< The plane to simulate

    // Some useful variables
    float                                    fMMPerTick;             //< Convert ticks in us to mm
    float                                    fMMPerWire;             //< Convert wire pitch to mm
    float                                    fTanThetaTW;            //< tangent of angle in time/wire space
    float                                    fTanTheta;              //< tangent in euclidean space
    float                                    fSinTheta;              //< sine of this angle
    float                                    fCosTheta;              //< save some time by also doing cosine
    std::vector<size_t>                      fTickEndPoints;         //< Tick endpoints to overlay

    icarusutil::TimeVec                      fFFTTimeVec;            //< Local time vector

    using FFTPointer = std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>>;

    FFTPointer                               fFFT;                   //< Object to handle thread safe FFT

    const geo::Geometry*                     fGeometry;              //< pointer to the Geometry service
    const detinfo::DetectorProperties*       fDetector;              //< Pointer to the detector properties
    icarusutil::SignalShapingICARUSService*  fSignalShapingService;  //< Access to the response functions
};

FakeParticle::FakeParticle(fhicl::ParameterSet const &pset)
{
    this->configure(pset);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

FakeParticle::~FakeParticle()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void FakeParticle::configure(fhicl::ParameterSet const &pset)
{
    fWireEndPoints     = pset.get<std::vector<size_t>>("WireEndPoints", std::vector<size_t>()={25,550});
    fStartTick         = pset.get<size_t             >("StartTick",                               1000);
    fStartAngle        = pset.get<float              >("StartAngle",                                45);
    fNumElectronsPerMM = pset.get<int                >("NumElectronsPerMM",                       6000);
    fPlaneToSimulate   = pset.get<size_t             >("PlaneToSimulate",                            2);
                                  
    fGeometry           = art::ServiceHandle<geo::Geometry const>{}.get();
    fDetector           = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fSignalShapingService = art::ServiceHandle<icarusutil::SignalShapingICARUSService>{}.get();

    // Convert ticks in us to mm by taking the drift velocity and multiplying by the tick period
    double driftVelocity = fDetector->DriftVelocity() * 10.;   // should be mm/us
    double samplingRate  = fDetector->SamplingRate() / 1000.;  // sampling rate returned in ns

    fMMPerTick = driftVelocity * samplingRate;

    fMMPerWire = fGeometry->WirePitch() * 10.;  // wire pitch returned in cm, want mm

    // Get slope (tan(theta)) and related angles
    fTanTheta   = std::tan(fStartAngle * M_PI / 180.);
    fCosTheta   = 1. / sqrt(1. + fTanTheta * fTanTheta);
    fSinTheta   = fTanTheta * fCosTheta;

    fTanThetaTW = fTanTheta * fMMPerWire / fMMPerTick;

    // Constrain x range
    fWireEndPoints[1] = std::min(fWireEndPoints[1],size_t(576));

    // Now compute the tick start/end points
    fTickEndPoints.push_back(fStartTick);
    fTickEndPoints.push_back(size_t(std::round(fTanThetaTW * float(fWireEndPoints[1] - fWireEndPoints[0]) + fStartTick)));

    // Check to see if we have run off the end of our frame (max 576,4096)
    if (fTickEndPoints[1] > size_t(4096))
    {
        // Probably should worry about the 90 degree case...
        if (std::fabs(fStartAngle) < 90)
        {
            fWireEndPoints[1] = std::round(float(4096 - fStartTick) / fTanTheta + fWireEndPoints[0]);
            fTickEndPoints[1] = 4096;
        }
        else
            fWireEndPoints[1] = fWireEndPoints[0];
    }

    // Now set up our plans for doing the convolution
    int numberTimeSamples = fDetector->NumberTimeSamples();

    fFFTTimeVec.resize(numberTimeSamples,0.);

    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(numberTimeSamples);

    return;
}

void FakeParticle::overlayFakeParticle(ArrayFloat& waveforms)
{
    // We have assumed the input waveform array will have 576 wires... 
    // Our "range" must be contained within that. By assumption we start at wire 0, so really just need
    // to set the max range 
    size_t maxWire = std::min(waveforms.size(),fWireEndPoints[1]);

    // Create a temporary waveform to handle the input charge
//    icarusutil::TimeVec tempWaveform(waveforms[0].size(),0.);

    // Also recover the gain
    float asicGain = fSignalShapingService->GetASICGain(0) * fDetector->SamplingRate() / 1000.;  // something like 67.4

    // Get a base channel number for the plane we want
    raw::ChannelID_t channel = fGeometry->PlaneWireToChannel(fPlaneToSimulate, 0);

    // Recover the response function information for this channel
    const icarus_tool::IResponse& response = fSignalShapingService->GetResponse(channel);

    // Also want the time offset for this channel
    int timeOffset = fSignalShapingService->ResponseTOffset(channel);

    // Loop over the wire range
    for(size_t wireIdx = fWireEndPoints[0]; wireIdx < maxWire; wireIdx++)
    {
        // As the track angle becomes more parallel to the electron drift direction there will be more charge
        // deposited per mm that will impact a given wire. So we need to try to accommodate this here.
        // Start by computing the starting tick (based on how far we have gone so far)
        size_t startTick = size_t(fTanThetaTW * float(wireIdx - fWireEndPoints[0])) + fTickEndPoints[0];

        // If this has gone outside the maximum tick then we are done
        if (!(startTick < fTickEndPoints[1] && startTick < fFFTTimeVec.size())) break;

        // Begin by computing the number of ticks for this given wire, converted to tick units and 
        // always at least one tick
        size_t deltaTicks = size_t(fMMPerWire * std::fabs(fTanTheta)) + 1;
        size_t endTick    = startTick + deltaTicks;

        // Trim back if we look to step outside of the max range 
        endTick = std::min(endTick,fTickEndPoints[1]);
        endTick = std::min(endTick,fFFTTimeVec.size());

        // Ok, now loop through the ticks and deposit charge. 
        // This will be the number of electrons per mm (input) 
        // x the total arc length of track for this step
        float arcLenPerWire    = fMMPerWire / fCosTheta;
        float arcLenPerTick    = arcLenPerWire / float(deltaTicks);
        float numElectronsWire = fNumElectronsPerMM * arcLenPerWire;
        float fracElecPerTick  = numElectronsWire * arcLenPerTick / arcLenPerWire;

        // Make sure the waveform is zeroed
        std::fill(fFFTTimeVec.begin(),fFFTTimeVec.end(),0.);

        // Now loop through the ticks and deposit charge
        for(size_t tickIdx = startTick; tickIdx < endTick; tickIdx++)
            fFFTTimeVec[tickIdx] = fracElecPerTick / asicGain;

        // Convolute with the response functions
//        fSignalShapingService->Convolute(channel, fFFTTimeVec);
        fFFT->convolute(fFFTTimeVec, response.getConvKernel(), timeOffset);

        // And now add this to the waveform in question 
        VectorFloat& waveform = waveforms[wireIdx];

        std::transform(waveform.begin(),waveform.end(),fFFTTimeVec.begin(),waveform.begin(),std::plus<float>());
    }

    return;
}

DEFINE_ART_CLASS_TOOL(FakeParticle)
} // namespace daq
