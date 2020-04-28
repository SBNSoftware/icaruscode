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

/**
 * If defined, a hack to make sure DetectorClocksService knows about the new
 * hardware trigger time is enabled.
 * This is violating art/LArSoft recommended practices, and it is not even
 * useful in ICARUS where the
 * @ref DetectorClocksElectronicsStartTime "electronics time start"
 * is _determined_ by the hardware trigger.
 */
#undef ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK


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
#include "CLHEP/Random/JamesRandom.h"
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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"
// LArSoft libraries
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#ifdef ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // FIXME: this is not portable
#endif // ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK
#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"
#include "lardataobj/Simulation/sim.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "tools/IGenNoise.h"
#include "icarus_signal_processing/ICARUSFFT.h"

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
    
    void MakeADCVec(std::vector<short>& adc, icarusutil::TimeVec const& noise,
                    icarusutil::TimeVec const& charge, float ped_mean) const;

    using TPCIDVec  = std::vector<geo::TPCID>;
    
    art::InputTag                fDriftEModuleLabel; ///< module making the ionization electrons
    bool                         fProcessAllTPCs;    ///< If true we process all TPCs
    unsigned int                 fCryostat;          ///< If ProcessAllTPCs is false then cryostat to use
    TPCIDVec                     fTPCVec;            ///< List of TPCs to process for this instance of the module
    raw::Compress_t              fCompression;       ///< compression type to use
    unsigned int                 fNTimeSamples;      ///< number of ADC readout samples in all readout frames (per event)
    std::map< double, int >      fShapingTimeOrder;
    
    bool                         fSimDeadChannels;   ///< if True, simulate dead channels using the ChannelStatus service.  If false, do not simulate dead channels
    bool                         fSuppressNoSignal;  ///< If no signal on wire (simchannel) then suppress the channel
    bool                         fSmearPedestals;    ///< If True then we smear the pedestals
    int                          fNumChanPerMB;      ///< Number of channels per motherboard
    
    std::vector<std::unique_ptr<icarus_tool::IGenNoise>> fNoiseToolVec; ///< Tool for generating noise
    
    bool                         fMakeHistograms;
    bool                         fTest; // for forcing a test case
    std::vector<sim::SimChannel> fTestSimChannel_v;
    size_t                       fTestWire;
    std::vector<size_t>          fTestIndex;
    std::vector<double>          fTestCharge;
    
    TH1F*                        fSimCharge;
    TH2F*                        fSimChargeWire;
    
    // Random engines
    CLHEP::HepRandomEngine&      fPedestalEngine;
    CLHEP::HepRandomEngine&      fUncNoiseEngine;
    CLHEP::HepRandomEngine&      fCorNoiseEngine;

    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float                  adcsaturation = 4095;
    
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

    using FFTPointer = std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>>;
    FFTPointer                              fFFT;                   //< Object to handle thread safe FFT
    
    //services
    const geo::GeometryCore&                fGeometry;
    icarusutil::SignalShapingICARUSService* fSignalShapingService;  //< Access to the response functions
    
}; // class SimWireICARUS
DEFINE_ART_MODULE(SimWireICARUS)
    
//-------------------------------------------------
SimWireICARUS::SimWireICARUS(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fPedestalEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "pedestal", pset, "SeedPedestal"))
    , fUncNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "noise",    pset, "Seed"))
    , fCorNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "cornoise", pset, "Seed"))
    , fGeometry(*lar::providerFrom<geo::Geometry>())
{
    this->reconfigure(pset);
    
    produces< std::vector<raw::RawDigit>   >();
    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;
    
    return;
}
//-------------------------------------------------
SimWireICARUS::~SimWireICARUS() {}
//-------------------------------------------------
void SimWireICARUS::reconfigure(fhicl::ParameterSet const& p)
{
    fDriftEModuleLabel= p.get< art::InputTag       >("DriftEModuleLabel",             "largeant");
    fProcessAllTPCs   = p.get< bool                >("ProcessAllTPCs",                     false);
    fCryostat         = p.get< unsigned int        >("Cryostat",                               0);
    fSimDeadChannels  = p.get< bool                >("SimDeadChannels",                    false);
    fSuppressNoSignal = p.get< bool                >("SuppressNoSignal",                   false);
    fMakeHistograms   = p.get< bool                >("MakeHistograms",                     false);
    fSmearPedestals   = p.get< bool                >("SmearPedestals",                      true);
    fNumChanPerMB     = p.get< int                 >("NumChanPerMB",                          32);
    fTest             = p.get< bool                >("Test",                               false);
    fTestWire         = p.get< size_t              >("TestWire",                               0);
    fTestIndex        = p.get< std::vector<size_t> >("TestIndex",          std::vector<size_t>());
    fTestCharge       = p.get< std::vector<double> >("TestCharge",         std::vector<double>());

    using TPCValsPair = std::pair<unsigned int, unsigned int>; // Assume cryostat, TPC 
    using TPCValsVec  = std::vector<TPCValsPair>;

    TPCValsVec tempIDVec = p.get< TPCValsVec >("TPCVec", TPCValsVec());

    for(const auto& idPair : tempIDVec)
    {
        fTPCVec.push_back(geo::TPCID(idPair.first,idPair.second));
    }

    if(fTestIndex.size() != fTestCharge.size())
        throw cet::exception(__FUNCTION__)<<"# test pulse mismatched: check TestIndex and TestCharge fcl parameters...";
    
    std::vector<fhicl::ParameterSet> noiseToolParamSetVec = p.get<std::vector<fhicl::ParameterSet>>("NoiseGenToolVec");
    
    for(auto& noiseToolParams : noiseToolParamSetVec) {
        fNoiseToolVec.push_back(art::make_tool<icarus_tool::IGenNoise>(noiseToolParams));
    }
    //Map the Shaping Times to the entry position for the noise ADC
    //level in fNoiseFactInd and fNoiseFactColl
    fShapingTimeOrder = { {0.6, 0}, {1, 1}, {1.3, 2}, {3.0, 3} };

    //detector properties information
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    fNTimeSamples = detprop->NumberTimeSamples();
    
    fSignalShapingService = art::ServiceHandle<icarusutil::SignalShapingICARUSService>{}.get();

    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(fNTimeSamples);
    
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
    
    //channel status for simulating dead channels
    const lariov::ChannelStatusProvider& ChannelStatusProvider = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    
#ifdef ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK
    // In case trigger simulation is run in the same job...
    // FIXME:  You should not be calling preProcessEvent
    art::ServiceHandle<detinfo::DetectorClocksServiceStandard>()
      ->preProcessEvent(evt,art::ScheduleContext::invalid());
#endif // ICARUSCODE_SIMWIREICARUS_TRIGGERTIMEHACK

    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    const raw::ChannelID_t maxChannel = fGeometry.Nchannels();

    //--------------------------------------------------------------------
    //
    // Get the SimChannels, which we will use to produce RawDigits
    //
    //--------------------------------------------------------------------
    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(maxChannel,nullptr);
    if(!fTest)
    {
        std::vector<const sim::SimChannel*> chanHandle;
        evt.getView(fDriftEModuleLabel,chanHandle);
        
        for(const auto& simChannel : chanHandle) channels.at(simChannel->Channel()) = simChannel;
    }
    else
        for(const auto& testChannel : fTestSimChannel_v) channels.at(testChannel.Channel()) = &testChannel;
    
    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(maxChannel);
    //--------------------------------------------------------------------
    //
    // Loop over channels a second time and produce the RawDigits by adding together
    // pedestal, noise, and direct & induced charges
    //
    //--------------------------------------------------------------------
    
    // vectors for working in the following for loop
    std::vector<short>  adcvec(fNTimeSamples, 0);
    icarusutil::TimeVec chargeWork(fNTimeSamples,0.);
    icarusutil::TimeVec zeroCharge(fNTimeSamples,0.);
    icarusutil::TimeVec noisetmp(fNTimeSamples,0.);
    
    // make sure chargeWork is correct size
    if (chargeWork.size() < fNTimeSamples) throw std::range_error("SimWireICARUS: chargeWork vector too small");
    
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Let the tools know to update to the next event
    for(const auto& noiseTool : fNoiseToolVec) noiseTool->nextEvent();

    // The original implementation would allow the option to skip channels for which there was no MC signal
    // present. We want to update this so that if there is an MC signal on any wire in a common group (a
    // motherboard) then we keep all of those wires. This so we can implment noise mitigation techniques
    // with the simulation
    //
    
    // Here we determine the first and last channel numbers based on whether we are outputting a single TPC or all
    using ChannelPair    = std::pair<raw::ChannelID_t,raw::ChannelID_t>;
    using ChannelPairVec = std::vector<ChannelPair>;

    ChannelPairVec channelPairVec;
   
    for (geo::TPCID const& tpcid : fTPCVec) 
    {
        for (geo::PlaneGeo const& plane : fGeometry.IteratePlanes(tpcid)) 
        {
            raw::ChannelID_t const planeStartChannel = fGeometry.PlaneWireToChannel({ plane.ID(), 0U });
            raw::ChannelID_t const planeEndChannel = fGeometry.PlaneWireToChannel({ plane.ID(), plane.Nwires() - 1U }) + 1;

            channelPairVec.emplace_back(planeStartChannel, planeEndChannel);            
        } // for planes in TPC
    }

    // Ok, define the structure for the MB info....
    using MBWithSignalSet = std::set<raw::ChannelID_t>;
    
    MBWithSignalSet mbWithSignalSet;
    
    for(const ChannelPair& channelPair : channelPairVec)
    {
        // If we are not suppressing the signal then we need to make sure there is an entry in the set for every motherboard
        if (!fSuppressNoSignal)
        {
            raw::ChannelID_t firstMBIdx(channelPair.first / fNumChanPerMB);
            raw::ChannelID_t endMBIdx(channelPair.second / fNumChanPerMB);
        
            for(raw::ChannelID_t mbIdx = firstMBIdx; mbIdx < endMBIdx; mbIdx++) mbWithSignalSet.insert(mbIdx);
        }
        else
        {
            for(const auto& simChan : channels)
            {
                if (simChan)
                {
                    raw::ChannelID_t channel = simChan->Channel();
                
                    if (channel >= channelPair.first && channel < channelPair.second) mbWithSignalSet.insert(channel/fNumChanPerMB);
                }
            }
        }
    }
    
    // Ok, now we can simply loop over MB's...
    for(const auto& mb : mbWithSignalSet)
    {
        raw::ChannelID_t baseChannel = fNumChanPerMB * mb;
        
        // And for a given MB we can loop over the channels it contains
        for(raw::ChannelID_t channel = baseChannel; channel < baseChannel + fNumChanPerMB; channel++)
        {
            //clean up working vectors from previous iteration of loop
            adcvec.resize(fNTimeSamples, 0);  //compression may have changed the size of this vector
            noisetmp.resize(fNTimeSamples, 0.);     //just in case
            
            //use channel number to set some useful numbers
            std::vector<geo::WireID> widVec = fGeometry.ChannelToWire(channel);
            size_t                   plane  = widVec[0].Plane;
            
            //Get pedestal with random gaussian variation
            float ped_mean = pedestalRetrievalAlg.PedMean(channel);
            
            if (fSmearPedestals )
            {
                CLHEP::RandGaussQ rGaussPed(fPedestalEngine, 0.0, pedestalRetrievalAlg.PedRms(channel));
                ped_mean += rGaussPed.fire();
            }
            
            //Generate Noise
            double noise_factor(0.);
            auto   tempNoiseVec = fSignalShapingService->GetNoiseFactVec();
            double shapingTime  = fSignalShapingService->GetShapingTime(channel);
            double gain         = fSignalShapingService->GetASICGain(channel) * detprop->SamplingRate() * 1.e-3; // Gain returned is electrons/us, this converts to electrons/tick
            int    timeOffset   = fSignalShapingService->ResponseTOffset(channel);
            
            // Recover the response function information for this channel
            const icarus_tool::IResponse& response = fSignalShapingService->GetResponse(channel);

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
            fNoiseToolVec[plane]->generateNoise(fUncNoiseEngine,
                                                fCorNoiseEngine,
                                                noisetmp,
                                                noise_factor,
                                                channel);
            
            // Recover the SimChannel (if one) for this channel
            const sim::SimChannel* simChan = channels[channel];
            
            // If there is something on this wire, and it is not dead, then add the signal to the wire
            if(simChan && !(fSimDeadChannels && (ChannelStatusProvider.IsBad(channel) || !ChannelStatusProvider.IsPresent(channel))))
            {
                std::fill(chargeWork.begin(), chargeWork.end(), 0.);
                
                // loop over the tdcs and grab the number of electrons for each
                for(size_t tick = 0; tick < fNTimeSamples; tick++)
                {
                    int tdc = ts->TPCTick2TDC(tick);
                    
                    // continue if tdc < 0
                    if( tdc < 0 ) continue;
                    
                    double charge = simChan->Charge(tdc);  // Charge returned in number of electrons
                    
                    chargeWork[tick] += charge/gain;  // # electrons / (# electrons/tick)
                } // loop over tdcs
                // now we have the tempWork for the adjacent wire of interest
                // convolve it with the appropriate response function
                fFFT->convolute(chargeWork, response.getConvKernel(), timeOffset);
                
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

        }
    }
    
    evt.put(std::move(digcol));
    
    return;
}
//-------------------------------------------------
void SimWireICARUS::MakeADCVec(std::vector<short>& adcvec, icarusutil::TimeVec const& noisevec,
                               icarusutil::TimeVec const& chargevec, float ped_mean) const
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
