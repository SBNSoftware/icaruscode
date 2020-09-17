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
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitCharacterizationAlg.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

#include <Eigen/Core>

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
    using WaveformVec  = std::vector<WireTuple>;
    using WaveformList = std::list<WireTuple*>;

    void saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >&, WireTuple&);
    void saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >&, raw::ChannelID_t&, float, float, caldata::RawDigitVector&);
    
    // Define the structuring element - currently just a vector of vectors
    using StructuringElement = std::vector<std::vector<short>>;

    // Fcl parameters.
    std::string                          fDigitModuleLabel;      ///< The full collection of hits
    bool                                 fOutputHistograms;      ///< Output histograms?
    bool                                 fOutputWaveforms;       ///< Output waveforms?

    art::TFileDirectory*                 fHistDirectory;

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    
    // Once defined the structuring element will not change
    size_t                               fStructuringElementWireSize;
    size_t                               fStructuringElementTickSize;
    StructuringElement                   fStructuringElement;
    
    // Correction algorithms
    caldata::RawDigitCharacterizationAlg fCharacterizationAlg;

    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*             fGeometry;             ///< pointer to Geometry service
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
RawDigitSmoother::RawDigitSmoother(fhicl::ParameterSet const & pset) : EDProducer{pset},
                                   fNumEvent(0),
                                   fCharacterizationAlg(pset.get<fhicl::ParameterSet>("CharacterizationAlg")),
                                   fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
{
    
    fGeometry = lar::providerFrom<geo::Geometry>();
    
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
    fDigitModuleLabel           = pset.get<std::string>("DigitModuleLabel",           "daq");
    fStructuringElementWireSize = pset.get<size_t>     ("StructuringElementWireSize",     5);
    fStructuringElementTickSize = pset.get<size_t>     ("StructuringElementTickSize",     5);
    fOutputHistograms           = pset.get< bool      >("OutputHistograms",           false);
    fOutputWaveforms            = pset.get< bool      >("OutputWaveforms",            false);

    fStructuringElement.resize(fStructuringElementWireSize);
    
    // Create a rectangular structuring element to start with
    for(auto& row : fStructuringElement) row.resize(fStructuringElementTickSize,1);
    
    // If asked, define the global histograms
    if (fOutputHistograms)
    {
        // Access ART's TFileService, which will handle creating and writing
        // histograms and n-tuples for us.
        art::ServiceHandle<art::TFileService> tfs;
        
        fHistDirectory = tfs.get();
        
        // Make a directory for these histograms
 //       art::TFileDirectory dir = fHistDirectory->mkdir(Form("ROIPlane_%1zu",fPlane));
    }
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
        erosionRawDigit->reserve(digitVecHandle->size());
        dilationRawDigit->reserve(digitVecHandle->size());
        edgeRawDigit->reserve(digitVecHandle->size());
        differenceRawDigit->reserve(digitVecHandle->size());
        averageRawDigit->reserve(digitVecHandle->size());
        medianRawDigit->reserve(digitVecHandle->size());
        
        unsigned int maxChannels = fGeometry->Nchannels();
        
        // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
        // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
        std::vector<const raw::RawDigit*> rawDigitVec;
        
        // Ugliness to fill the pointer vector...
        for(size_t idx = 0; idx < digitVecHandle->size(); idx++) rawDigitVec.push_back(&digitVecHandle->at(idx));
        
        // Sort (use a lambda to sort by channel id)
        std::sort(rawDigitVec.begin(),rawDigitVec.end(),[](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});
        
        // Get size of input data vectors
        size_t rawDataSize = rawDigitVec.front()->Samples();

        // Get an instance of our input waveform list and digit vect
        WaveformList inputWaveformList;
        WaveformVec  wireTupleVec;
        
        // First we create vector which contains the data
        for(size_t idx = 0; idx < fStructuringElementWireSize; idx++) wireTupleVec.push_back(WireTuple(0,0.,0.,caldata::RawDigitVector(rawDataSize)));
        
        // Now set the address of each of these in the list
        for(size_t idx = 0; idx < fStructuringElementWireSize; idx++) inputWaveformList.push_back(&wireTupleVec[idx]);

        // ok, make containers for the various things we are going to calculate
        WireTuple erosionTuple    = WireTuple(0,0.,0.,caldata::RawDigitVector(rawDataSize, 0));
        WireTuple dilationTuple   = WireTuple(0,0.,0.,caldata::RawDigitVector(rawDataSize, 0));
        WireTuple edgeTuple       = WireTuple(0,0.,0.,caldata::RawDigitVector(rawDataSize, 0));
        WireTuple differenceTuple = WireTuple(0,0.,0.,caldata::RawDigitVector(rawDataSize, 0));
        WireTuple averageTuple    = WireTuple(0,0.,0.,caldata::RawDigitVector(rawDataSize, 0));
        WireTuple medianTuple     = WireTuple(0,0.,0.,caldata::RawDigitVector(rawDataSize, 0));
        
        caldata::RawDigitVector& erosionVec    = std::get<3>(erosionTuple);
        caldata::RawDigitVector& dilationVec   = std::get<3>(dilationTuple);
        caldata::RawDigitVector& edgeVec       = std::get<3>(edgeTuple);
        caldata::RawDigitVector& differenceVec = std::get<3>(differenceTuple);
        caldata::RawDigitVector& averageVec    = std::get<3>(averageTuple);
        caldata::RawDigitVector& medianVec     = std::get<3>(medianTuple);

        // Use an index for the last valid waveform...
        size_t validIndex = 0;
        
        // On the very inside loop we are going to keep track of ADC values in a vector... we can speed things up by pre determining some parameters
        size_t maxAdcBinSize(0);
        
        for(const auto& rowVec : fStructuringElement)
        {
            for(const auto& structElemVal : rowVec)
                if (structElemVal) maxAdcBinSize++;
        }
        
        std::vector<short> adcBinValVec(maxAdcBinSize, 0);

        // Avoid creating and destroying a vector each loop... make a single one here
        caldata::RawDigitVector inputAdcVector(rawDataSize);
        
        geo::WireID lastWireID = fGeometry->ChannelToWire(rawDigitVec.front()->Channel())[0];
    
        // Commence looping over raw digits
        for(const auto& rawDigit : rawDigitVec)
        {
            raw::ChannelID_t channel = rawDigit->Channel();
            
            if (channel >= maxChannels) continue;

            // Decode the channel and make sure we have a valid one
            std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
            
            // Look to see if we have crossed to another plane
            if (lastWireID.asPlaneID().cmp(wids[0].asPlaneID()) != 0)
            {
                // Dispose of the end set of RawDigits (in order)
                WaveformList::iterator inputWaveItr = inputWaveformList.begin();
                
                std::advance(inputWaveItr, fStructuringElementWireSize/2);
                
                while(++inputWaveItr != inputWaveformList.end())
                {
                    saveRawDigits(erosionRawDigit,    **inputWaveItr);
                    saveRawDigits(dilationRawDigit,   **inputWaveItr);
                    saveRawDigits(edgeRawDigit,       **inputWaveItr);
                    saveRawDigits(differenceRawDigit, **inputWaveItr);
                    saveRawDigits(averageRawDigit,    **inputWaveItr);
                    saveRawDigits(medianRawDigit,     **inputWaveItr);
                }
                
                // Reset the valid waveforms index
                validIndex = 0;
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
            
            // If the buffer is "full" then we need to rotate the first to the end so we can reuse
            if (validIndex == inputWaveformList.size())
            {
                inputWaveformList.push_back(inputWaveformList.front());
                inputWaveformList.pop_front();
            }
            else validIndex++;
            
            // Find the right entry
            WaveformList::iterator inputWaveItr = inputWaveformList.begin();
            WaveformList::iterator midWaveItr   = inputWaveItr;
            
            std::advance(inputWaveItr, validIndex - 1);
            std::advance(midWaveItr,   validIndex / 2);

            caldata::RawDigitVector& rawadc = std::get<3>(**inputWaveItr);
            
            // And now uncompress
            raw::Uncompress(rawDigit->ADCs(), inputAdcVector, rawDigit->Compression());
            
            float truncMean;
            float rmsVal;
            float pedCorVal;
            
            // Recover the mean and rms for this waveform
            fCharacterizationAlg.getMeanRmsAndPedCor(inputAdcVector, channel, plane, wire, truncMean, rmsVal, pedCorVal);
            
            // Recover the database version of the pedestal
            float pedestal = fPedestalRetrievalAlg.PedMean(channel);

            std::transform(inputAdcVector.begin(),inputAdcVector.end(),rawadc.begin(),std::bind(std::minus<short>(),std::placeholders::_1,pedCorVal));
            
            std::get<0>(**inputWaveItr) = channel;
            std::get<1>(**inputWaveItr) = pedestal;
            std::get<2>(**inputWaveItr) = rmsVal;
            
            // Finally, at this point we are prepared to do some work!
            if (validIndex == inputWaveformList.size())
            {
                
                raw::ChannelID_t midChannel                     = std::get<0>(**midWaveItr);
                float            midPedestal                    = std::get<1>(**midWaveItr);
                float            midRmsVal                      = std::get<2>(**midWaveItr);
                size_t           halfStructuringElementTickSize = fStructuringElementTickSize / 2;
                
                caldata::RawDigitVector& currentVec = std::get<3>(**midWaveItr);
                
                // Fill the edge bins with the pedestal value
                for(size_t adcBinIdx = 0; adcBinIdx < halfStructuringElementTickSize; adcBinIdx++)
                {
                    size_t adcLastBinIdx = rawDataSize - adcBinIdx - 1;
                    
                    erosionVec[adcBinIdx]        = midPedestal;
                    erosionVec[adcLastBinIdx]    = midPedestal;
                    dilationVec[adcBinIdx]       = midPedestal;
                    dilationVec[adcLastBinIdx]   = midPedestal;
                    edgeVec[adcBinIdx]           = midPedestal;
                    edgeVec[adcLastBinIdx]       = midPedestal;
                    differenceVec[adcBinIdx]     = midPedestal;
                    differenceVec[adcLastBinIdx] = midPedestal;
                    averageVec[adcBinIdx]        = midPedestal;
                    averageVec[adcLastBinIdx]    = midPedestal;
                    medianVec[adcBinIdx]         = midPedestal;
                    medianVec[adcLastBinIdx]     = midPedestal;
                }

                // Ok, buckle up!
                // Loop will run from half the structuring element to size less half the structuring element. Edges will simply be what they were
                for(size_t adcBinIdx = halfStructuringElementTickSize; adcBinIdx < erosionVec.size() - halfStructuringElementTickSize; adcBinIdx++)
                {
                    size_t rowIdx(0);
                    size_t adcBinVecIdx(0);
                    
                    // Outside loop over vectors
                    for(const auto& curTuple : inputWaveformList)
                    {
                        const caldata::RawDigitVector& curAdcVec = std::get<3>(*curTuple);
                        
                        for(size_t colIdx = 0; colIdx < fStructuringElementTickSize; colIdx++)
                        {
                            if (fStructuringElement[rowIdx][colIdx]) adcBinValVec[adcBinVecIdx++] = curAdcVec[colIdx + adcBinIdx - halfStructuringElementTickSize];
                        }
                        
                        rowIdx++;
                    }
                    
                    std::sort(adcBinValVec.begin(),adcBinValVec.end());
                    
                    erosionVec[adcBinIdx]    =  adcBinValVec.front();
                    dilationVec[adcBinIdx]   =  adcBinValVec.back();
                    edgeVec[adcBinIdx]       = (dilationVec[adcBinIdx] - currentVec[adcBinIdx]) + midPedestal;
                    differenceVec[adcBinIdx] = (dilationVec[adcBinIdx] - erosionVec[adcBinIdx]) + midPedestal;
                    averageVec[adcBinIdx]    = (dilationVec[adcBinIdx] + erosionVec[adcBinIdx]) / 2;
                    medianVec[adcBinIdx]     =  adcBinValVec[adcBinValVec.size()/2];
                }

                saveRawDigits(erosionRawDigit,    midChannel, midPedestal, midRmsVal, std::get<3>(erosionTuple));
                saveRawDigits(dilationRawDigit,   midChannel, midPedestal, midRmsVal, std::get<3>(dilationTuple));
                saveRawDigits(edgeRawDigit,       midChannel, midPedestal, midRmsVal, std::get<3>(edgeTuple));
                saveRawDigits(differenceRawDigit, midChannel, midPedestal, midRmsVal, std::get<3>(differenceTuple));
                saveRawDigits(averageRawDigit,    midChannel, midPedestal, midRmsVal, std::get<3>(averageTuple));
                saveRawDigits(medianRawDigit,     midChannel, midPedestal, midRmsVal, std::get<3>(medianTuple));
            }
            else if (validIndex <= fStructuringElementWireSize / 2)
            {
                saveRawDigits(erosionRawDigit,    **inputWaveItr);
                saveRawDigits(dilationRawDigit,   **inputWaveItr);
                saveRawDigits(edgeRawDigit,       **inputWaveItr);
                saveRawDigits(differenceRawDigit, **inputWaveItr);
                saveRawDigits(averageRawDigit,    **inputWaveItr);
                saveRawDigits(medianRawDigit,     **inputWaveItr);
            }
        }
    }
/*
    if (fOutputWaveforms)
    {
        // Try to limit to the wire number (since we are already segregated by plane)
        std::vector<geo::WireID> wids  = fGeometry->ChannelToWire(channel);
        size_t                   cryo  = wids[0].Cryostat;
        size_t                   tpc   = wids[0].TPC;
        size_t                   plane = wids[0].Plane;
        size_t                   wire  = wids[0].Wire;
        
        // Make a directory for these histograms
        art::TFileDirectory dir = fHistDirectory->mkdir(Form("ROIPlane_%1zu/c%1zu/c%1zut%1zuwire_%05zu",fPlane,cnt,cryo,tpc,wire));
        
        // We keep track of four histograms:
        try
        {
            //            origWaveHist   = dir.make<TProfile>(Form("Inp_%03zu_ctw%01zu/%01zu/%05zu",cnt,cryo,tpc,wire), "Waveform", waveform.size(),      0, waveform.size(),      -500., 500.);
            histogramMap[icarus_tool::WAVEFORM] =
            dir.make<TProfile>(Form("Wfm_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Waveform", waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::EROSION] =
            dir.make<TProfile>(Form("Ero_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Erosion",  waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::DILATION] =
            dir.make<TProfile>(Form("Dil_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Dilation", waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::AVERAGE] =
            dir.make<TProfile>(Form("Ave_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Average",  waveformSize, 0, waveformSize, -500., 500.);
            histogramMap[icarus_tool::DIFFERENCE] =
            dir.make<TProfile>(Form("Dif_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Average",  waveformSize, 0, waveformSize, -500., 500.);
            
            // This is a kludge so that the ROI histogram ends up in the same diretory as the waveforms
            histogramMap[ROIHISTOGRAM] =
            dir.make<TProfile>(Form("ROI_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "ROI",      waveformSize, 0, waveformSize, -500., 500.);
            
            // Also, if smoothing then we would like to keep track of the original waveform too
            histogramMap[WAVEFORMHIST] =
            dir.make<TProfile>(Form("Inp_%03zu_ctw%01zu-%01zu-%01zu-%05zu",cnt,cryo,tpc,plane,wire), "Waveform", waveformSize, 0, waveformSize, -500., 500.);
        } catch(...)
        {
            std::cout << "Caught exception trying to make new hists, tpc,plane,cnt,wire: " << tpc << ", " << fPlane << ", " << cnt << ", " << wire << std::endl;
        }
    }
*/
    // Add tracks and associations to event.
    event.put(std::move(erosionRawDigit),    "erosion");
    event.put(std::move(dilationRawDigit),   "dilation");
    event.put(std::move(edgeRawDigit),       "edge");
    event.put(std::move(differenceRawDigit), "difference");
    event.put(std::move(averageRawDigit),    "average");
    event.put(std::move(medianRawDigit),     "median");
    
    return;
}

void RawDigitSmoother::saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit,
                                     WireTuple&                                    wireTuple)
{
    raw::ChannelID_t&        channel     = std::get<0>(wireTuple);
    float                    pedestal    = std::get<1>(wireTuple);
    float                    rms         = std::get<2>(wireTuple);
    caldata::RawDigitVector& rawDigitVec = std::get<3>(wireTuple);
    
    filteredRawDigit->emplace_back(channel, rawDigitVec.size(), rawDigitVec, raw::kNone);
    filteredRawDigit->back().SetPedestal(pedestal,rms);
    
    return;
}

void RawDigitSmoother::RawDigitSmoother::saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit,
                                                       raw::ChannelID_t&                             channel,
                                                       float                                         pedestal,
                                                       float                                         rms,
                                                       caldata::RawDigitVector&                      rawDigitVec)
{
    filteredRawDigit->emplace_back(channel, rawDigitVec.size(), rawDigitVec, raw::kNone);
    filteredRawDigit->back().SetPedestal(pedestal,rms);
    
    return;
}

//----------------------------------------------------------------------------
/// End job method.
void RawDigitSmoother::endJob()
{
    mf::LogInfo("RawDigitSmoother") << "Looked at " << fNumEvent << " events" << std::endl;
}
