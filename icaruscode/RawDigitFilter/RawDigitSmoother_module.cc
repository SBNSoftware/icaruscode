////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitSmoother
// Module Type: producer
// File:        RawDigitSmoother_module.cc
//
//              This module implements a two dimensional morphological filter
//              to apply to plane by plane images with the intent to enhance the
//              signal regions. The primary aim is to aid pattern recognition
//
// Configuration parameters:
//
// DigitModuleLabel      - the source of the RawDigit collection
//
// Created by Tracy Usher (usher@slac.stanford.edu) on July 29, 2018
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include "icaruscode/RawDigitFilter/RawDigitFilterAlgs/RawDigitCharacterizationAlg.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"


class RawDigitSmoother : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit RawDigitSmoother(fhicl::ParameterSet const & pset);
    virtual ~RawDigitSmoother();

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:
    
    // Set up our container for the waveforms
    // We'll keep things in a tuple so we can also keep track of the pedestal and rms for output
    using WireTuple    = std::tuple<raw::ChannelID_t,float,float,caldata::RawDigitVector>;
    using WaveformList = std::list<WireTuple>;

    void saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >&, WireTuple&);
    
    // Define the structuring element - currently just a vector of vectors
    using StructuringElement = std::vector<std::vector<short>>;

    // Fcl parameters.
    std::string                          fDigitModuleLabel;      ///< The full collection of hits

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    
    // Once defined the structuring element will not change
    size_t                               fStructuringElementSize;
    StructuringElement                   fStructuringElement;
    
    // Correction algorithms
    caldata::RawDigitCharacterizationAlg fCharacterizationAlg;

    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*             fGeometry;             ///< pointer to Geometry service
    detinfo::DetectorProperties const*   fDetectorProperties;   ///< Detector properties service
    const lariov::DetPedestalProvider&   fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
};

DEFINE_ART_MODULE(RawDigitSmoother)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitSmoother::RawDigitSmoother(fhicl::ParameterSet const & pset) :
                                   fNumEvent(0),
                                   fCharacterizationAlg(pset.get<fhicl::ParameterSet>("CharacterizationAlg")),
                                   fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
{
    
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    configure(pset);
    produces<std::vector<raw::RawDigit>>("erosion");
    produces<std::vector<raw::RawDigit>>("dilation");
    produces<std::vector<raw::RawDigit>>("edge");
    produces<std::vector<raw::RawDigit>>("average");
    produces<std::vector<raw::RawDigit>>("difference");
    produces<std::vector<raw::RawDigit>>("median");

    // Report.
    mf::LogInfo("RawDigitSmoother") << "RawDigitSmoother configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitSmoother::~RawDigitSmoother()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitSmoother::configure(fhicl::ParameterSet const & pset)
{
    fDigitModuleLabel       = pset.get<std::string>("DigitModuleLabel",       "daq");
    fStructuringElementSize = pset.get<size_t>     ("StructuringElementSize",     5);
    
    fStructuringElement.resize(fStructuringElementSize);
    
    // Create a box structuring element to start with
    for(auto& row : fStructuringElement) row.resize(fStructuringElementSize,1);
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitSmoother::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    
//    art::TFileDirectory dir = tfs->mkdir(Form("RawDigitSmoother"));

    fCharacterizationAlg.initializeHists(tfs);

    return;
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method.
///
void RawDigitSmoother::produce(art::Event & event)
{
    ++fNumEvent;
    
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<std::vector<raw::RawDigit> > erosionRawDigit(new std::vector<raw::RawDigit>);
    std::unique_ptr<std::vector<raw::RawDigit> > dilationRawDigit(new std::vector<raw::RawDigit>);
    std::unique_ptr<std::vector<raw::RawDigit> > edgeRawDigit(new std::vector<raw::RawDigit>);
    std::unique_ptr<std::vector<raw::RawDigit> > differenceRawDigit(new std::vector<raw::RawDigit>);
    std::unique_ptr<std::vector<raw::RawDigit> > averageRawDigit(new std::vector<raw::RawDigit>);
    std::unique_ptr<std::vector<raw::RawDigit> > medianRawDigit(new std::vector<raw::RawDigit>);

    erosionRawDigit->clear();
    dilationRawDigit->clear();
    edgeRawDigit->clear();
    differenceRawDigit->clear();
    averageRawDigit->clear();
    medianRawDigit->clear();

    // Read in the digit List object(s).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    event.getByLabel(fDigitModuleLabel, digitVecHandle);
    
    // Require a valid handle
    if (digitVecHandle.isValid() && digitVecHandle->size() > 0)
    {
        unsigned int maxChannels = fGeometry->Nchannels();
        
        // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
        // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
        std::vector<const raw::RawDigit*> rawDigitVec;
        
        // Ugliness to fill the pointer vector...
        for(size_t idx = 0; idx < digitVecHandle->size(); idx++) rawDigitVec.push_back(&digitVecHandle->at(idx));
        
        // Sort (use a lambda to sort by channel id)
        std::sort(rawDigitVec.begin(),rawDigitVec.end(),[](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});
        
        // Get an instance of our input waveform list
        WaveformList inputWaveformList;
        
        geo::WireID lastWireID = fGeometry->ChannelToWire(rawDigitVec.front()->Channel())[0];
    
        // Commence looping over raw digits
        for(const auto& rawDigit : rawDigitVec)
        {
            raw::ChannelID_t channel = rawDigit->Channel();
        
            // The below try-catch block may no longer be necessary
            // Decode the channel and make sure we have a valid one
            std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
        
            if (channel >= maxChannels) continue;
            
            // Look to see if we have crossed to another plane
            if (lastWireID.asPlaneID().cmp(wids[0].asPlaneID()) != 0)
            {
                // Dispose of the end set of RawDigits
                for(size_t idx = 0; idx < fStructuringElementSize/2; idx++)
                {
                    saveRawDigits(erosionRawDigit,    inputWaveformList.back());
                    saveRawDigits(dilationRawDigit,   inputWaveformList.back());
                    saveRawDigits(edgeRawDigit,       inputWaveformList.back());
                    saveRawDigits(differenceRawDigit, inputWaveformList.back());
                    saveRawDigits(averageRawDigit,    inputWaveformList.back());
                    saveRawDigits(medianRawDigit,     inputWaveformList.back());

                    inputWaveformList.pop_back();
                }
                
                // Clear the container to start over
                inputWaveformList.clear();
            }
            
            // Update the last wire id before we forget...
            lastWireID = wids[0];
        
            // Recover plane and wire in the plane
            unsigned int plane = wids[0].Plane;
            unsigned int wire  = wids[0].Wire;
            
            if (rawDigit->Samples() < 1)
            {
                std::cout << "****>> Found zero length raw digit buffer, channel: " << channel << ", plane: " << plane << ", wire: " << wire << std::endl;
                continue;
            }
            
            while(inputWaveformList.size() >= fStructuringElementSize) inputWaveformList.pop_front();

            inputWaveformList.emplace_back(channel,0.,0.,caldata::RawDigitVector(rawDigit->Samples()));

            caldata::RawDigitVector& rawadc = std::get<3>(inputWaveformList.back());
            
            // And now uncompress
            raw::Uncompress(rawDigit->ADCs(), rawadc, rawDigit->Compression());
            
            float truncMean;
            float rmsVal;
            float pedCorVal;
            
            // Recover the mean and rms for this waveform
            fCharacterizationAlg.getMeanRmsAndPedCor(rawadc, channel, plane, wire, truncMean, rmsVal, pedCorVal);
            
            // Recover the database version of the pedestal
            float pedestal = fPedestalRetrievalAlg.PedMean(channel);

            std::transform(rawadc.begin(),rawadc.end(),rawadc.begin(),std::bind(std::minus<short>(),std::placeholders::_1,pedCorVal));
            
            std::get<1>(inputWaveformList.back()) = pedestal;
            std::get<2>(inputWaveformList.back()) = rmsVal;
            
            // Finally, at this point we are prepared to do some work!
            if (inputWaveformList.size() == fStructuringElementSize)
            {
                size_t halfStructuringElementSize = fStructuringElementSize / 2;
                
                WaveformList::iterator midChanItr = inputWaveformList.begin();
                
                std::advance(midChanItr, halfStructuringElementSize);

                // ok, make copies of this waveform for the erosion, dilation and avrerage
                WireTuple erosionTuple    = *midChanItr;
                WireTuple dilationTuple   = *midChanItr;
                WireTuple edgeTuple       = *midChanItr;
                WireTuple differenceTuple = *midChanItr;
                WireTuple averageTuple    = *midChanItr;
                WireTuple medianTuple     = *midChanItr;

                caldata::RawDigitVector& erosionVec    = std::get<3>(erosionTuple);
                caldata::RawDigitVector& dilationVec   = std::get<3>(dilationTuple);
                caldata::RawDigitVector& edgeVec       = std::get<3>(edgeTuple);
                caldata::RawDigitVector& differenceVec = std::get<3>(edgeTuple);
                caldata::RawDigitVector& averageVec    = std::get<3>(averageTuple);
                caldata::RawDigitVector& medianVec     = std::get<3>(medianTuple);
                
                caldata::RawDigitVector& currentVec    = std::get<3>(*midChanItr);

                // Ok, buckle up!
                // Loop will run from half the structuring element to size less half the structuring element. Edges will simply be what they were
                for(size_t adcBinIdx = halfStructuringElementSize; adcBinIdx < erosionVec.size() - halfStructuringElementSize; adcBinIdx++)
                {
                    std::vector<short> adcBinValVec;
                    size_t             rowIdx(0);
                    
                    adcBinValVec.reserve(fStructuringElementSize*fStructuringElementSize);
                    
                    // Outside loop over vectors
                    for(const auto& curTuple : inputWaveformList)
                    {
                        const caldata::RawDigitVector& curAdcVec = std::get<3>(curTuple);
                        
                        for(size_t colIdx = 0; colIdx < fStructuringElementSize; colIdx++)
                        {
                            if (fStructuringElement[rowIdx][colIdx]) adcBinValVec.push_back(curAdcVec.at(colIdx + adcBinIdx - halfStructuringElementSize));
                        }
                    }
                    
                    std::sort(adcBinValVec.begin(),adcBinValVec.end());
                    
                    erosionVec.at(adcBinIdx)    = adcBinValVec.front();
                    dilationVec.at(adcBinIdx)   = adcBinValVec.back();
                    edgeVec.at(adcBinIdx)       = dilationVec.at(adcBinIdx) - currentVec.at(adcBinIdx) + std::get<1>(edgeTuple);
                    differenceVec.at(adcBinIdx) = dilationVec.at(adcBinIdx) - erosionVec.at(adcBinIdx) + std::get<1>(edgeTuple);
                    averageVec.at(adcBinIdx)    = (erosionVec.at(adcBinIdx) + dilationVec.at(adcBinIdx)) / 2;
                    medianVec.at(adcBinIdx)     = adcBinValVec.at(adcBinValVec.size()/2);
                }
                
                saveRawDigits(erosionRawDigit,    erosionTuple);
                saveRawDigits(dilationRawDigit,   dilationTuple);
                saveRawDigits(edgeRawDigit,       edgeTuple);
                saveRawDigits(differenceRawDigit, differenceTuple);
                saveRawDigits(averageRawDigit,    averageTuple);
                saveRawDigits(medianRawDigit,     medianTuple);
            }
            else if (inputWaveformList.size() <= fStructuringElementSize / 2)
            {
                saveRawDigits(erosionRawDigit,    inputWaveformList.back());
                saveRawDigits(dilationRawDigit,   inputWaveformList.back());
                saveRawDigits(edgeRawDigit,       inputWaveformList.back());
                saveRawDigits(differenceRawDigit, inputWaveformList.back());
                saveRawDigits(averageRawDigit,    inputWaveformList.back());
                saveRawDigits(medianRawDigit,     inputWaveformList.back());
            }
        }
    }
    
    // Add tracks and associations to event.
    event.put(std::move(erosionRawDigit),    "erosion");
    event.put(std::move(dilationRawDigit),   "dilation");
    event.put(std::move(edgeRawDigit),       "edge");
    event.put(std::move(differenceRawDigit), "difference");
    event.put(std::move(averageRawDigit),    "average");
    event.put(std::move(medianRawDigit),     "median");
}

void RawDigitSmoother::saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit,
                                     WireTuple&                                    wireTuple)
{
    raw::ChannelID_t&        channel     = std::get<0>(wireTuple);
    float                    pedestal    = std::get<1>(wireTuple);
    float                    rms         = std::get<2>(wireTuple);
    caldata::RawDigitVector& rawDigitVec = std::get<3>(wireTuple);
    
    filteredRawDigit->emplace_back(raw::RawDigit(channel, rawDigitVec.size(), rawDigitVec, raw::kNone));
    filteredRawDigit->back().SetPedestal(pedestal,rms);
    
    return;
}

//----------------------------------------------------------------------------
/// End job method.
void RawDigitSmoother::endJob()
{
    mf::LogInfo("RawDigitSmoother") << "Looked at " << fNumEvent << " events" << std::endl;
}
