///////////////////////////////////////////////////////////////////////
// $Id: SimWireICARUS.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireICARUS class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////
// C/C++ standard library
#include <stdexcept> // std::range_error
#include <vector>
#include <string>
#include <algorithm> // std::fill()
#include <functional>
#include <random>
#include <chrono>
// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
// ROOT libraries
#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
// art library and utilities
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// art extensions
#include "nutools/RandomUtils/NuRandomService.h"
// LArSoft libraries
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // FIXME: this is not portable
#include "icaruscode/Utilities/SignalShapingServiceICARUS.h"
#include "lardataobj/Simulation/sim.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "tools/IGenNoise.h"
using namespace util;
///Detector simulation of raw signals on wires
namespace detsim {
    
// Base class for creation of raw signals on wires.
class SimWireICARUS : public art::EDProducer
{
public:
    
    explicit SimWireICARUS(fhicl::ParameterSet const& pset);
    virtual ~SimWireICARUS();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
private:
    
    void MakeADCVec(std::vector<short>& adc, std::vector<float> const& noise,
                    std::vector<double> const& charge, float ped_mean) const;
    
    std::string                  fDriftEModuleLabel; ///< module making the ionization electrons
    raw::Compress_t              fCompression;       ///< compression type to use
    size_t                       fNTicks;                ///< number of ticks of the clock
    unsigned int                 fNTimeSamples;      ///< number of ADC readout samples in all readout frames (per event)
    std::map< double, int >      fShapingTimeOrder;
    
    bool                         fSimDeadChannels;   ///< if True, simulate dead channels using the ChannelStatus service.  If false, do not simulate dead channels
    bool                         fSuppressNoSignal;  ///< If no signal on wire (simchannel) then suppress the channel
    bool                         fSmearPedestals;    ///< If True then we smear the pedees
    
    std::vector<std::unique_ptr<icarus_tool::IGenNoise>> fNoiseToolVec; ///< Tool for generating noise
    
    bool                         fMakeHistograms;
    bool                         fTest; // for forcing a test case
    std::vector<sim::SimChannel> fTestSimChannel_v;
    size_t                       fTestWire;
    std::vector<size_t>          fTestIndex;
    std::vector<double>          fTestCharge;
    int                          fSample; // for histograms, -1 means no histos
    
    TH1F*                        fSimCharge;
    TH2F*                        fSimChargeWire;
    
    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float adcsaturation = 4095;
    
    // little helper class to hold the params of each charge dep
    class ResponseParams {
    public:
        ResponseParams(double charge, size_t time) : m_charge(charge), m_time(time) {}
        double getCharge() { return m_charge; }
        size_t getTime()   { return m_time; }
    private:
        double m_charge;
        size_t m_time;
    };
    
    //services
    const geo::GeometryCore& fGeometry;
    
}; // class SimWireICARUS
DEFINE_ART_MODULE(SimWireICARUS)
    
//-------------------------------------------------
SimWireICARUS::SimWireICARUS(fhicl::ParameterSet const& pset)
: fGeometry(*lar::providerFrom<geo::Geometry>())
{
    this->reconfigure(pset);
    
    produces< std::vector<raw::RawDigit>   >();
    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;
    
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed" and "SeedPedestal"
    art::ServiceHandle<rndm::NuRandomService> Seeds;
    Seeds->createEngine(*this, "HepJamesRandom", "noise",    pset, "Seed");
    Seeds->createEngine(*this, "HepJamesRandom", "cornoise", pset, "Seed");
    Seeds->createEngine(*this, "HepJamesRandom", "pedestal", pset, "SeedPedestal");
}
//-------------------------------------------------
SimWireICARUS::~SimWireICARUS() {}
//-------------------------------------------------
void SimWireICARUS::reconfigure(fhicl::ParameterSet const& p)
{
    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");
    fSimDeadChannels  = p.get< bool                >("SimDeadChannels");
    fSuppressNoSignal = p.get< bool                >("SuppressNoSignal");
    fMakeHistograms   = p.get< bool                >("MakeHistograms", false);
    fSample           = p.get< int                 >("Sample");
    fSmearPedestals   = p.get< bool                >("SmearPedestals", true);
    fTest             = p.get< bool                >("Test");
    fTestWire         = p.get< size_t              >("TestWire");
    fTestIndex        = p.get< std::vector<size_t> >("TestIndex");
    fTestCharge       = p.get< std::vector<double> >("TestCharge");
    
    if(fTestIndex.size() != fTestCharge.size())
        throw cet::exception(__FUNCTION__)<<"# test pulse mismatched: check TestIndex and TestCharge fcl parameters...";
    
    std::vector<fhicl::ParameterSet> noiseToolParamSetVec = p.get<std::vector<fhicl::ParameterSet>>("NoiseGenToolVec");
    
    for(auto& noiseToolParams : noiseToolParamSetVec)
        fNoiseToolVec.push_back(art::make_tool<icarus_tool::IGenNoise>(noiseToolParams));
    //Map the Shaping Times to the entry position for the noise ADC
    //level in fNoiseFactInd and fNoiseFactColl
    fShapingTimeOrder = { {0.6, 0}, {1, 1}, {1.3, 2}, {3.0, 3} };
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    fNTimeSamples = detprop->NumberTimeSamples();
    
    return;
}
//-------------------------------------------------
void SimWireICARUS::beginJob()
{
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    // If in test mode create a test data set
    if(fTest)
    {
        if(fGeometry.Nchannels()<=fTestWire)
            throw cet::exception(__FUNCTION__)<<"Invalid test wire channel: "<<fTestWire;
        std::vector<unsigned int> channels;
        for(auto const& plane_id : fGeometry.IteratePlaneIDs())
            channels.push_back(fGeometry.PlaneWireToChannel(plane_id.Plane,fTestWire));
        double xyz[3] = { std::numeric_limits<double>::max() };
        for(auto const& ch : channels)
        {
            fTestSimChannel_v.push_back(sim::SimChannel(ch));
            for(size_t i=0; i<fTestIndex.size(); ++i)
            {
                fTestSimChannel_v.back().AddIonizationElectrons(-1,
                                                                fTestIndex[i],
                                                                fTestCharge[i],
                                                                xyz,
                                                                std::numeric_limits<double>::max());
            }
        }
    }
    
    fSimCharge     = tfs->make<TH1F>("fSimCharge", "simulated charge", 150, 0, 1500);
    fSimChargeWire = tfs->make<TH2F>("fSimChargeWire", "simulated charge", 5600,0.,5600.,500, 0, 1500);
    
    return;
}
//-------------------------------------------------
void SimWireICARUS::endJob()
{}
void SimWireICARUS::produce(art::Event& evt)
{
    //--------------------------------------------------------------------
    //
    // Get all of the services we will be using
    //
    //--------------------------------------------------------------------
    
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    
    //get rng for pedestals
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("pedestal");
    
    //channel status for simulating dead channels
    const lariov::ChannelStatusProvider& ChannelStatusProvider = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
    //get the FFT
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->ReinitializeFFT(fNTimeSamples,fFFT->FFTOptions(),fFFT->FFTFitBins());
    fNTicks = fFFT->FFTSize();
    if ( fNTicks%2 != 0 )
        LOG_DEBUG("SimWireICARUS") << "Warning: FFTSize " << fNTicks << " not a power of 2. "
        << "May cause issues in (de)convolution.\n";
    if ( fNTimeSamples > fNTicks )
        mf::LogError("SimWireICARUS") << "Cannot have number of readout samples "
        << fNTimeSamples << " greater than FFTSize " << fNTicks << "!";
    
    //TimeService
    art::ServiceHandle<detinfo::DetectorClocksServiceStandard> tss;
    
    // In case trigger simulation is run in the same job...
    tss->preProcessEvent(evt);
    auto const* ts = tss->provider();
    
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    const size_t N_CHANNELS = fGeometry.Nchannels();
    
    //Get N_RESPONSES from SignalShapingService, on the fly
    // flag added to use nominal one response per plane or multiple responses
    // per plane and scaling for YZ dependent responses
    // or data driven field responses
    art::ServiceHandle<util::SignalShapingServiceICARUS> sss;
    
    //--------------------------------------------------------------------
    //
    // Get the SimChannels, which we will use to produce RawDigits
    //
    //--------------------------------------------------------------------
    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(N_CHANNELS,nullptr);
    if(!fTest)
    {
        std::vector<const sim::SimChannel*> chanHandle;
        evt.getView(fDriftEModuleLabel,chanHandle);
        
        for(const auto& simChannel : chanHandle) channels.at(simChannel->Channel()) = simChannel;
    }else{
        for(const auto& testChannel : fTestSimChannel_v) channels.at(testChannel.Channel()) = &testChannel;
    }
    
    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(N_CHANNELS);
    //--------------------------------------------------------------------
    //
    // Loop over channels a second time and produce the RawDigits by adding together
    // pedestal, noise, and direct & induced charges
    //
    //--------------------------------------------------------------------
    
    // vectors for working in the following for loop
    std::vector<short>  adcvec(fNTimeSamples, 0);
    std::vector<double> chargeWork(fNTicks,0.);
    std::vector<double> zeroCharge(fNTicks,0.);
    std::vector<float>  noisetmp(fNTicks,0.);
    
    // make sure chargeWork is correct size
    if (chargeWork.size() < fNTimeSamples) throw std::range_error("SimWireICARUS: chargeWork vector too small");
    
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // loop over the collected responses
    //   this is needed because hits generate responses on adjacent wires!
    for(unsigned int channel = 0; channel < N_CHANNELS; channel++)
    {
        // get the sim::SimChannel for this channel
        // Look up first so we can suppress before doing any work if no signal
        const sim::SimChannel* sc = channels[channel];
        
        if (fSuppressNoSignal && !sc) continue;
        
        //clean up working vectors from previous iteration of loop
        adcvec.resize(fNTimeSamples, 0);  //compression may have changed the size of this vector
        noisetmp.resize(fNTicks, 0.);     //just in case
        
        //use channel number to set some useful numbers
        std::vector<geo::WireID> widVec = fGeometry.ChannelToWire(channel);
        size_t                   plane  = widVec[0].Plane;
        
        //Get pedestal with random gaussian variation
        float ped_mean = pedestalRetrievalAlg.PedMean(channel);
        
        if (fSmearPedestals )
        {
            CLHEP::RandGaussQ rGaussPed(engine, 0.0, pedestalRetrievalAlg.PedRms(channel));
            ped_mean += rGaussPed.fire();
        }
        
        //Generate Noise
        double noise_factor(0.);
        auto   tempNoiseVec = sss->GetNoiseFactVec();
        double shapingTime  = sss->GetShapingTime(0);
        
        if (fShapingTimeOrder.find( shapingTime ) != fShapingTimeOrder.end() )
            noise_factor = tempNoiseVec[plane].at( fShapingTimeOrder.find( shapingTime )->second );
        //Throw exception...
        else
        {
            throw cet::exception("SimWireICARUS")
            << "\033[93m"
            << "Shaping Time received from signalservices_microboone.fcl is not one of allowed values"
            << std::endl
            << "Allowed values: 0.6, 1.0, 1.3, 3.0 usec"
            << "\033[00m"
            << std::endl;
        }
        
        // Use the desired noise tool to actually generate the noise on this wire
        fNoiseToolVec[plane]->GenerateNoise(noisetmp, noise_factor, channel);
        
        double gain=sss->GetASICGain(channel) * detprop->SamplingRate() * 1.e-3; // Gain return is electrons/us, this converts to electrons/tick
        
        // If there is something on this wire, and it is not dead, then add the signal to the wire
        if(sc && !(fSimDeadChannels && (ChannelStatusProvider.IsBad(channel) || !ChannelStatusProvider.IsPresent(channel))))
        {
            std::fill(chargeWork.begin(), chargeWork.end(), 0.);
            
            // loop over the tdcs and grab the number of electrons for each
            for(int tick = 0; tick < (int)fNTicks; tick++)
            {
                int tdc = ts->TPCTick2TDC(tick);
                
                // continue if tdc < 0
                if( tdc < 0 ) continue;
                
                double charge = sc->Charge(tdc);  // Charge returned in number of electrons
                
                chargeWork[tick] += charge/gain;  // # electrons / (# electrons/tick)
            } // loop over tdcs
            // now we have the tempWork for the adjacent wire of interest
            // convolve it with the appropriate response function
            sss->Convolute(channel, chargeWork);
            
            // "Make" the ADC vector
            MakeADCVec(adcvec, noisetmp, chargeWork, ped_mean);
        }
        // "Make" an ADC vector with zero charge added
        else MakeADCVec(adcvec, noisetmp, zeroCharge, ped_mean);
        
        // add this digit to the collection;
        // adcvec is copied, not moved: in case of compression, adcvec will show
        // less data: e.g. if the uncompressed adcvec has 9600 items, after
        // compression it will have maybe 5000, but the memory of the other 4600
        // is still there, although unused; a copy of adcvec will instead have
        // only 5000 items. All 9600 items of adcvec will be recovered for free
        // and used on the next loop.
        raw::RawDigit rd(channel, fNTimeSamples, adcvec, fCompression);
        
        if(fMakeHistograms && plane==2)
        {
            short area = std::accumulate(adcvec.begin(),adcvec.end(),0,[](const auto& val,const auto& sum){return sum + val - 400;});
            
            if(area>0)
            {
                fSimCharge->Fill(area);
                fSimChargeWire->Fill(widVec[0].Wire,area);
            }
        }
        
        rd.SetPedestal(ped_mean);
        digcol->push_back(std::move(rd)); // we do move the raw digit copy, though
    }// end of 2nd loop over channels
    
    evt.put(std::move(digcol));
    
    return;
}
//-------------------------------------------------
void SimWireICARUS::MakeADCVec(std::vector<short>& adcvec, std::vector<float> const& noisevec,
                               std::vector<double> const& chargevec, float ped_mean) const
{
    for(unsigned int i = 0; i < fNTimeSamples; ++i)
    {
        float adcval = noisevec[i] + chargevec[i] + ped_mean;

        adcval = std::max(float(0.), std::min(adcval, adcsaturation));

        adcvec[i] = std::round(adcval);
    }// end loop over signal size
    // compress the adc vector using the desired compression scheme,
    // if raw::kNone is selected nothing happens to adcvec
    // This shrinks adcvec, if fCompression is not kNone.
    raw::Compress(adcvec, fCompression);
    
    return;
}
    
}
