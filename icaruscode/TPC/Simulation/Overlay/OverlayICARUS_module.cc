///////////////////////////////////////////////////////////////////////
// $Id: OverlayICARUS.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// OverlayICARUS class designed to simulate signal on a wire in the TPC
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
#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"
#include "lardataobj/Simulation/sim.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "tools/IOverlay.h"
#include "icarus_signal_processing/ICARUSFFT.h"

using namespace util;
///Detector simulation of raw signals on wires
namespace detsim {
    
// Base class for creation of raw signals on wires.
class OverlayICARUS : public art::EDProducer
{
public:
    
    explicit OverlayICARUS(fhicl::ParameterSet const& pset);
    virtual ~OverlayICARUS();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
private:
    
    void MakeADCVec(std::vector<short>& adc, icarusutil::TimeVec const& charge, float ped_mean) const;
    
    art::InputTag                fInputRawDataLabel; ///< Label for the underlying raw digit data
    art::InputTag                fDriftEModuleLabel; ///< module making the ionization electrons
    raw::Compress_t              fCompression;       ///< compression type to use
        
    bool                         fMakeHistograms;

    TH1F*                        fSimCharge;
    TH2F*                        fSimChargeWire;

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
    icarusutil::SignalShapingICARUSService* fSignalShapingService;  ///< Access to the response functions
    const lariov::DetPedestalProvider&      fPedestalRetrievalAlg;  ///< Keep track of an instance to the pedestal retrieval alg
    
}; // class OverlayICARUS
DEFINE_ART_MODULE(OverlayICARUS)
    
//-------------------------------------------------
OverlayICARUS::OverlayICARUS(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fGeometry(*lar::providerFrom<geo::Geometry>()),
      fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
{
    this->reconfigure(pset);
    
    produces< std::vector<raw::RawDigit>   >();
    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;
    
    return;
}

//-------------------------------------------------
OverlayICARUS::~OverlayICARUS() {}

//-------------------------------------------------
void OverlayICARUS::reconfigure(fhicl::ParameterSet const& p)
{
    fInputRawDataLabel = p.get< art::InputTag >("InputRawDataLabel",      "daq");
    fDriftEModuleLabel = p.get< art::InputTag >("DriftEModuleLabel", "largeant");

    //detector properties information
    auto const detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
    
    fSignalShapingService = art::ServiceHandle<icarusutil::SignalShapingICARUSService>{}.get();

    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(detprop.NumberTimeSamples());
    
    return;
}

//-------------------------------------------------
void OverlayICARUS::beginJob()
{
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    fSimCharge     = tfs->make<TH1F>("fSimCharge", "simulated charge", 150, 0, 1500);
    fSimChargeWire = tfs->make<TH2F>("fSimChargeWire", "simulated charge", 5600,0.,5600.,500, 0, 1500);
    
    return;
}

//-------------------------------------------------
void OverlayICARUS::endJob()
{}

void OverlayICARUS::produce(art::Event& evt)
{
    //--------------------------------------------------------------------
    //
    // Get all of the services we will be using
    //
    //--------------------------------------------------------------------

    // Recocver the input RawDigits we are going to add our signal too
    art::Handle<std::vector<raw::RawDigit>> inputRawDigitHandle;

    evt.getByLabel(fInputRawDataLabel, inputRawDigitHandle);

    if (!inputRawDigitHandle.isValid()) throw std::runtime_error("Failed to read back RawDigits for overlay!");

    // Now recover the SimChannel information which will be overlaid
    art::Handle<std::vector<sim::SimChannel>> simChanHandle;
    evt.getByLabel(fDriftEModuleLabel, simChanHandle);

    if (!simChanHandle.isValid()) throw std::runtime_error("Failed to recover the SimChannel information for the overlay");

    // Keep track of the SimChannel information by channel
    using SimChannelMap = std::unordered_map<raw::ChannelID_t, const sim::SimChannel*>;

    SimChannelMap simChannelMap;
        
    for(const auto& simChannel : *simChanHandle) simChannelMap[simChannel.Channel()] = &simChannel;
    
    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(inputRawDigitHandle->size());
    
    // vectors for working in the following for loop
    std::vector<short>  adcvec;
    icarusutil::TimeVec zeroCharge;
    
    //detector properties information
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    
    // The outer loop is over the input RawDigits which will always be written out
    for(const auto& rawDigit : *inputRawDigitHandle)        
    {
        // Recover the channel
        raw::ChannelID_t channel = rawDigit.Channel();
            
        //use channel number to set some useful numbers
        std::vector<geo::WireID> widVec = fGeometry.ChannelToWire(channel);

        // We skip channels that are not connected to a physical wire
        if (!widVec.empty())
        {
            size_t plane = widVec[0].Plane;    

            // Recover the ADC values and copy to local vector
            const raw::RawDigit::ADCvector_t& adcVector = rawDigit.ADCs();

            // Make sure local vector is correct size (and note the vector above may be compressed so can't use its size yet)
            adcvec.resize(rawDigit.Samples(),0);
            zeroCharge.resize(rawDigit.Samples(),0);

            // and do the copyh
            raw::Uncompress(adcVector, adcvec, rawDigit.Compression());

            // Check for the existence of a SimChannel for this channel
            SimChannelMap::const_iterator simChanItr = simChannelMap.find(channel);

            if (simChanItr != simChannelMap.end())
            {
                // Recover the response function information for this channel
                const icarus_tool::IResponse& response = fSignalShapingService->GetResponse(channel);

                // Recover the SimChannel (if one) for this channel
                const sim::SimChannel* simChan = simChanItr->second;

                // Allocate local vector to hold the deposited charge
                icarusutil::TimeVec chargeWork(adcvec.size(),0.);

                // Need the to convert from deposited number of electrons to ADC units
                double gain = fSignalShapingService->GetASICGain(channel) * sampling_rate(clockData) * 1.e-3; // Gain returned is electrons/us, this converts to electrons/tick

                // Loop through the simchannel energy deposits
                for(const auto& tdcide : simChan->TDCIDEMap())
                {
                    unsigned int tdc = tdcide.first;

                    // We need to convert this to a tick...
                    int tick = clockData.TPCTDC2Tick(tdc);

                    // If out of range what is right thing to do?
                    if (tick < 0 ||tick > int(adcvec.size()))
                    {
                        std::cout << "tick out of range: " << tick << ", tdc: " << tdc << std::endl;
                        continue;
                    }

                    // Recover the charge for this tick
                    double charge = simChan->Charge(tdc);

                    chargeWork[tick] += charge / gain;
                }

                // now we have the tempWork for the adjacent wire of interest
                // convolve it with the appropriate response function
                fFFT->convolute(chargeWork, response.getConvKernel(), fSignalShapingService->ResponseTOffset(channel));

                //Get the pedestal and rms from the input waveform
                float pedestal = fPedestalRetrievalAlg.PedMean(channel);

                // "Make" the ADC vector
                MakeADCVec(adcvec, chargeWork, pedestal);
            }
        
            if(fMakeHistograms && plane==2)
            {
                short area = std::accumulate(adcvec.begin(),adcvec.end(),0,[](const auto& val,const auto& sum){return sum + val - 400;});
            
                if(area>0)
                {
                    fSimCharge->Fill(area);
                    fSimChargeWire->Fill(widVec[0].Wire,area);
                }
            }
        }
            
        // add this digit to the collection;
        // adcvec is copied, not moved: in case of compression, adcvec will show
        // less data: e.g. if the uncompressed adcvec has 9600 items, after
        // compression it will have maybe 5000, but the memory of the other 4600
        // is still there, although unused; a copy of adcvec will instead have
        // only 5000 items. All 9600 items of adcvec will be recovered for free
        // and used on the next loop.
        raw::RawDigit rd(channel, adcvec.size(), adcvec, fCompression);
        
        rd.SetPedestal(rawDigit.GetPedestal(),rawDigit.GetSigma());
        digcol->push_back(std::move(rd)); // we do move the raw digit copy, though
    }
    
    evt.put(std::move(digcol));
    
    return;
}
//-------------------------------------------------
void OverlayICARUS::MakeADCVec(std::vector<short>& adcvec, icarusutil::TimeVec const& chargevec, float ped_mean) const
{
    for(unsigned int i = 0; i < adcvec.size(); ++i)
    {
        // We offset for the channel's default pedestal
        float adcval = chargevec[i] + ped_mean;

        // Now make sure to limit the range to that allowed for 12 bits
        adcval = std::max(float(0.), std::min(adcval, adcsaturation));

        // Remove the temporary offset while we add the signal to the existing vector. 
        adcvec[i] += std::round(adcval - ped_mean);
    }// end loop over signal size
    // compress the adc vector using the desired compression scheme,
    // if raw::kNone is selected nothing happens to adcvec
    // This shrinks adcvec, if fCompression is not kNone.
    raw::Compress(adcvec, fCompression);
    
    return;
}
    
}
