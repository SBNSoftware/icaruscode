////////////////////////////////////////////////////////////////////////
//
// Class:       MCDecoderICARUSTPCwROI
// Module Type: producer
// File:        MCDecoderICARUSTPCwROI_module.cc
//
//              The intent of this module is to both "decode" artdaq fragments
//              and convert to RawDigits and also to do the initial noise 
//              filtering, specifically the coherent noise. 
//
// Configuration parameters:
//
// DigitModuleLabel      - the source of the RawDigit collection
//
//
// Modeled after example from Mike Wang (mwang@fnal.gov)
// Copied/Modified by Tracy Usher (usher@slac.stanford.edu) on January 27, 2020
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>

#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/Globals.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/cpu_timer.h"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_arena.h"
#include "tbb/spin_mutex.h"
#include "tbb/concurrent_vector.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"         // This for outputting the ROIs
#include "lardata/ArtDataHelper/WireCreator.h"

#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/DecoderTools/INoiseFilter.h"

namespace daq 
{

class MCDecoderICARUSTPCwROI : public art::ReplicatedProducer
{
public:

    // Copnstructors, destructor.
    explicit MCDecoderICARUSTPCwROI(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame);
    virtual ~MCDecoderICARUSTPCwROI();

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e, art::ProcessingFrame const& frame);
    virtual void beginJob(art::ProcessingFrame const& frame);
    virtual void endJob(art::ProcessingFrame const& frame);

    // Define the RawDigit collection
    using RawDigitCollection    = std::vector<raw::RawDigit>;
    using RawDigitCollectionPtr = std::unique_ptr<RawDigitCollection>;
    using WireCollection        = std::vector<recob::Wire>;
    using WireCollectionPtr     = std::unique_ptr<WireCollection>;
    using ConcurrentRawDigitCol = tbb::concurrent_vector<raw::RawDigit>;
    using ConcurrentWireCol     = tbb::concurrent_vector<recob::Wire>;

    // Define data structures for organizing the decoded fragments
    // The idea is to form complete "images" organized by "logical" TPC. Here we are including
    // both the active AND inactive channels so the noise processing will be correct
    // We build two lists here, first is the mapping to the actual image, 
    // the second will be a mapping to channel IDs where we assume the
    // order is the same between the two
    using PlaneIdxToImagePair   = std::pair<unsigned int,icarus_signal_processing::ArrayFloat>;
    using PlaneIdxToImageMap    = std::map<unsigned int,icarus_signal_processing::ArrayFloat>;
    using ChannelVec            = std::vector<raw::ChannelID_t>;
    using PlaneIdxToChannelPair = std::pair<unsigned int,ChannelVec>;
    using PlaneIdxToChannelMap  = std::map<unsigned int,ChannelVec>;

    using ChannelArrayPair      = std::pair<daq::INoiseFilter::ChannelPlaneVec,icarus_signal_processing::ArrayFloat>;
    using ChannelArrayPairVec   = std::vector<ChannelArrayPair>;

    // Function to do the work
    void processSingleImage(const detinfo::DetectorClocksData&,
                            const ChannelArrayPair&,
                            size_t,
                            ConcurrentRawDigitCol&,
                            ConcurrentRawDigitCol&,
                            ConcurrentRawDigitCol&,
                            ConcurrentWireCol&) const;

private:

    // Function to grab the input data and package
    void processSingleLabel(art::Event&,
                            const art::InputTag&, 
                            detinfo::DetectorClocksData const&,
                            ChannelArrayPairVec const&,
                            size_t const&,
                            ConcurrentRawDigitCol&,
                            ConcurrentRawDigitCol&,
                            ConcurrentRawDigitCol&,
                            ConcurrentWireCol&) const;

    class multiThreadImageProcessing
    {
    public:
        multiThreadImageProcessing(MCDecoderICARUSTPCwROI      const& parent,
                                   detinfo::DetectorClocksData const& clockData,
                                   ChannelArrayPairVec         const& channelArrayPairVec,
                                   size_t                      const& coherentNoiseGrouping,
                                   ConcurrentRawDigitCol&             concurrentRawDigits,
                                   ConcurrentRawDigitCol&             concurrentRawRawDigits,
                                   ConcurrentRawDigitCol&             coherentRawDigits,
                                   ConcurrentWireCol&                 concurrentROIs)
            : fMCDecoderICARUSTPCwROI(parent),
              fClockData{clockData},
              fChannelArrayPairVec(channelArrayPairVec),
              fCoherentNoiseGrouping(coherentNoiseGrouping),
              fConcurrentRawDigits(concurrentRawDigits),
              fConcurrentRawRawDigits(concurrentRawRawDigits),
              fCoherentRawDigits(coherentRawDigits),
              fConcurrentROIs(concurrentROIs)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
            {
                const ChannelArrayPair& channelArrayPair = fChannelArrayPairVec[idx];

                fMCDecoderICARUSTPCwROI.processSingleImage(fClockData, channelArrayPair, fCoherentNoiseGrouping, fConcurrentRawDigits, fConcurrentRawRawDigits, fCoherentRawDigits, fConcurrentROIs);
            }
        }
    private:
        const MCDecoderICARUSTPCwROI&      fMCDecoderICARUSTPCwROI;
        const detinfo::DetectorClocksData& fClockData;
        const ChannelArrayPairVec&         fChannelArrayPairVec;
        size_t                             fCoherentNoiseGrouping;
        ConcurrentRawDigitCol&             fConcurrentRawDigits;
        ConcurrentRawDigitCol&             fConcurrentRawRawDigits;
        ConcurrentRawDigitCol&             fCoherentRawDigits;
        ConcurrentWireCol&                 fConcurrentROIs;
    };

    // Function to save our RawDigits
    void saveRawDigits(const icarus_signal_processing::ArrayFloat&, 
                       const icarus_signal_processing::VectorFloat&, 
                       const icarus_signal_processing::VectorFloat&,
                       const icarus_signal_processing::VectorInt&,
                       ConcurrentRawDigitCol&) const;

    // Fcl parameters.
    std::vector<art::InputTag>                                  fRawDigitLabelVec;           ///< The input artdaq fragment label vector (for more than one)
    std::vector<std::string>                                    fOutInstanceLabelVec;        ///< The output instance labels to apply
    bool                                                        fOutputRawWaveform;          ///< Should we output pedestal corrected (not noise filtered)?
    bool                                                        fOutputCorrection;           ///< Should we output the coherent noise correction vectors?
    std::string                                                 fOutputRawWavePath;          ///< Path to assign to the output if asked for
    std::string                                                 fOutputCoherentPath;         ///< Path to assign to the output if asked for
    bool                                                        fDiagnosticOutput;           ///< Set this to get lots of messages
    size_t                                                      fCoherentNoiseGrouping;      ///< # channels in common for coherent noise

    const std::string                                           fLogCategory;                ///< Output category when logging messages

    // Statistics.
    int                                                         fNumEvent;             ///< Number of events seen.

    // Plane to ROP plane mapping
    using PlaneToROPPlaneMap   = std::map<geo::PlaneID,unsigned int>;
    using PlaneToWireOffsetMap = std::map<geo::PlaneID,raw::ChannelID_t>;
    using ROPToNumWiresMap     = std::map<unsigned int,unsigned int>;

    PlaneToROPPlaneMap                                          fPlaneToROPPlaneMap;
    PlaneToWireOffsetMap                                        fPlaneToWireOffsetMap;
    ROPToNumWiresMap                                            fROPToNumWiresMap;
    unsigned int                                                fNumROPs;

    using WirePlanePair              = std::pair<unsigned int, unsigned int>;
    using BoardWirePlanePair         = std::pair<unsigned int, WirePlanePair>;
    using ChannelToBoardWirePlaneMap = std::map<unsigned int, BoardWirePlanePair>;

    ChannelToBoardWirePlaneMap                                  fChannelToBoardWirePlaneMap;

    // Tools for decoding fragments depending on type
    std::vector<std::unique_ptr<INoiseFilter>>                  fDecoderToolVec;       ///< Decoder tools

    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*                                    fGeometry;             ///< pointer to Geometry service
    const icarusDB::IICARUSChannelMap*                          fChannelMap;
};

DEFINE_ART_MODULE(MCDecoderICARUSTPCwROI)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
MCDecoderICARUSTPCwROI::MCDecoderICARUSTPCwROI(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame) :
                        art::ReplicatedProducer(pset, frame),
                        fLogCategory("MCDecoderICARUSTPCwROI"),fNumEvent(0), fNumROPs(0)
{
    fGeometry   = art::ServiceHandle<geo::Geometry const>{}.get();
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    configure(pset);

    // Check the concurrency 
    int max_concurrency = art::Globals::instance()->nthreads(); 

    mf::LogDebug("MCDecoderICARUSTPCwROI") << "     ==> concurrency: " << max_concurrency << std::endl;

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& decoderToolParams = pset.get<fhicl::ParameterSet>("DecoderTool");
    
    fDecoderToolVec.resize(max_concurrency);
    
    for(auto& decoderTool : fDecoderToolVec)
    {
        // Get instance of tool
        decoderTool = art::make_tool<INoiseFilter>(decoderToolParams);
    }

    // Set up our "produces" 
    // Note that we can have multiple instances input to the module
    // Our convention will be to create a similar number of outputs with the same instance names
    for(const auto& instanceLabel : fOutInstanceLabelVec)
    {
        produces<std::vector<raw::RawDigit>>(instanceLabel);
        produces<std::vector<recob::Wire>>(instanceLabel);

        if (fOutputRawWaveform)
            produces<std::vector<raw::RawDigit>>(instanceLabel + fOutputRawWavePath);

        if (fOutputCorrection)
            produces<std::vector<raw::RawDigit>>(instanceLabel + fOutputCoherentPath);
    }

    // Set up a WireID to ROP plane number table
    PlaneToWireOffsetMap planeToLastWireOffsetMap; 

    for(size_t cryoIdx = 0; cryoIdx < 2; cryoIdx++)
    {
        for(size_t logicalTPCIdx = 0; logicalTPCIdx < 4; logicalTPCIdx++)
        {
            for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
            {
                geo::PlaneID planeID(cryoIdx,logicalTPCIdx,planeIdx);

                raw::ChannelID_t channel = fGeometry->PlaneWireToChannel(planeID.Plane, 0, planeID.TPC, planeID.Cryostat);

                readout::ROPID ropID = fGeometry->ChannelToROP(channel);

                fPlaneToROPPlaneMap[planeID]      = ropID.ROP;
                fPlaneToWireOffsetMap[planeID]    = channel;
                planeToLastWireOffsetMap[planeID] = fGeometry->PlaneWireToChannel(planeID.Plane, fGeometry->Nwires(planeID), planeID.TPC, planeID.Cryostat);
                fROPToNumWiresMap[ropID.ROP]      = fGeometry->Nwires(planeID);

                // Special case handling
//                if (ropID.ROP > 1) fROPToNumWiresMap[ropID.ROP] *= 2;

                if (ropID.ROP > fNumROPs) fNumROPs = ropID.ROP;

                // Watch for the middle induction and collection plane logical TPC split
                if (ropID.ROP > 1 && (logicalTPCIdx == 1 || logicalTPCIdx == 3))
                {
                    geo::PlaneID tempID(cryoIdx,logicalTPCIdx-1,planeIdx);

                    fPlaneToWireOffsetMap[planeID] = fPlaneToWireOffsetMap[tempID];
                    fROPToNumWiresMap[ropID.ROP]   = planeToLastWireOffsetMap[planeID] - fPlaneToWireOffsetMap[planeID];
                }

                // Diagnostic output if requested
                mf::LogDebug(fLogCategory) << "Initializing C/T/P: " << planeID.Cryostat << "/" << planeID.TPC << "/" << planeID.Plane << ", base channel: " << fPlaneToWireOffsetMap[planeID] << ", ROP: " << ropID << ", index: " << ropID.ROP;

            }
        }
    }

    fNumROPs++;

    // We need to build a mapping fronm channel to a readout board/wire pair
    // Get the board ids for this fragment
    const icarusDB::TPCReadoutBoardToChannelMap& readoutBoardToChannelMap = fChannelMap->getReadoutBoardToChannelMap();

    for(const auto& boardPair : readoutBoardToChannelMap)
    {
        // The board pair will give us the readout board and a vector of "wires"
        unsigned int readoutBoardID = boardPair.first;

        // Loop through the vector of wires on this board
        for(unsigned int wireIdx = 0; wireIdx < boardPair.second.second.size(); wireIdx++)
        {
            unsigned int channelID = boardPair.second.second[wireIdx].first;
            unsigned int planeID   = boardPair.second.second[wireIdx].second;

            fChannelToBoardWirePlaneMap[channelID] = BoardWirePlanePair(readoutBoardID,WirePlanePair(wireIdx,planeID));
        }
    }


    // Report.
    mf::LogInfo("MCDecoderICARUSTPCwROI") << "MCDecoderICARUSTPCwROI configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
MCDecoderICARUSTPCwROI::~MCDecoderICARUSTPCwROI()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void MCDecoderICARUSTPCwROI::configure(fhicl::ParameterSet const & pset)
{
    fRawDigitLabelVec           = pset.get<std::vector<art::InputTag>>("FragmentsLabelVec",   {"daq:PHYSCRATEDATA"});
    fOutInstanceLabelVec        = pset.get<std::vector<std::string>>  ("OutInstanceLabelVec",     {"PHYSCRATEDATA"});
    fOutputRawWaveform          = pset.get<bool                      >("OutputRawWaveform",                   false);
    fOutputCorrection           = pset.get<bool                      >("OutputCorrection",                    false);
    fOutputRawWavePath          = pset.get<std::string               >("OutputRawWavePath",                   "raw");
    fOutputCoherentPath         = pset.get<std::string               >("OutputCoherentPath",                  "Cor");
    fDiagnosticOutput           = pset.get<bool                      >("DiagnosticOutput",                    false);
    fCoherentNoiseGrouping      = pset.get<size_t                    >("CoherentGrouping",                       64);

}

//----------------------------------------------------------------------------
/// Begin job method.
void MCDecoderICARUSTPCwROI::beginJob(art::ProcessingFrame const&)
{ 
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
void MCDecoderICARUSTPCwROI::produce(art::Event & event, art::ProcessingFrame const&)
{
    ++fNumEvent;

    mf::LogDebug("MCDecoderICARUSTPCwROI") << "**** Processing raw data fragments ****" << std::endl;

    // Check the concurrency 
    int max_concurrency = tbb::this_task_arena::max_concurrency();

    mf::LogDebug("MCDecoderICARUSTPCwROI") << "     ==> concurrency: " << max_concurrency << std::endl;

    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    // Loop through the list of input daq fragment collections one by one 
    // We are not trying to multi thread at this stage because we are trying to control
    // overall memory usage at this level. We'll multi thread internally...
    size_t instanceIdx(0);

    for(const auto& rawDigitLabel : fRawDigitLabelVec)
    {
        art::Handle<artdaq::Fragments> daq_handle;
        event.getByLabel(rawDigitLabel, daq_handle);

        ConcurrentRawDigitCol concurrentRawDigits;
        ConcurrentRawDigitCol concurrentRawRawDigits;
        ConcurrentRawDigitCol coherentRawDigits;
        ConcurrentWireCol     concurrentROIs;

        PlaneIdxToImageMap   planeIdxToImageMap;
        PlaneIdxToChannelMap planeIdxToChannelMap;

        ChannelArrayPairVec  channelArrayPairVec(fNumROPs);

        // Because the arrays can be variable size we need to loop to initialize
        for(size_t ropIdx = 0; ropIdx < fNumROPs; ropIdx++)
        {
            ChannelArrayPair& channelArrayPair = channelArrayPairVec[ropIdx];

            channelArrayPair.first.resize(fROPToNumWiresMap[ropIdx]);
            channelArrayPair.second.resize(fROPToNumWiresMap[ropIdx],icarus_signal_processing::VectorFloat(4096));

            mf::LogDebug("MCDecoderICARUSTPCwROI") << "**> Initializing ropIdx: " << ropIdx << " channelPairVec to " << channelArrayPair.first.size() << " channels with " << channelArrayPair.second[0].size() << " ticks" << std::endl;
        }

        mf::LogDebug("MCDecoderICARUSTPCwROI") << "****> Let's get ready to rumble!" << std::endl;

        // Now let's process the resulting images
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);
    
        // ... repackage the input MC data to format suitable for noise processing
        processSingleLabel(event, rawDigitLabel, clockData, channelArrayPairVec, fCoherentNoiseGrouping, concurrentRawDigits, concurrentRawRawDigits, coherentRawDigits, concurrentROIs);

//        multiThreadImageProcessing imageProcessing(*this, clockData, channelArrayPairVec, fCoherentNoiseGrouping, concurrentRawDigits, concurrentRawRawDigits, coherentRawDigits, concurrentROIs);
//
//        tbb::parallel_for(tbb::blocked_range<size_t>(0, fNumROPs), imageProcessing);
    
        // Copy the raw digits from the concurrent vector to our output vector
        RawDigitCollectionPtr rawDigitCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(concurrentRawDigits.begin()), 
                                                                                                std::move_iterator(concurrentRawDigits.end()));
    
        // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
        std::sort(rawDigitCollection->begin(),rawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});

        // What did we get back?
        mf::LogDebug("MCDecoderICARUSTPCwROI") << "****> Total size of map: " << planeIdxToImageMap.size() << std::endl;
        for(const auto& planeImagePair : planeIdxToImageMap)
        {
            mf::LogDebug("MCDecoderICARUSTPCwROI") << "      - plane: " << planeImagePair.first << " has " << planeImagePair.second.size() << " wires" << std::endl;
        }
    
        // Now transfer ownership to the event store
        event.put(std::move(rawDigitCollection), fOutInstanceLabelVec[instanceIdx]);

        // Do the same to output the candidate ROIs
        WireCollectionPtr wireCollection = std::make_unique<std::vector<recob::Wire>>(std::move_iterator(concurrentROIs.begin()),
                                                                                      std::move_iterator(concurrentROIs.end()));

        std::sort(wireCollection->begin(),wireCollection->end(),[](const auto& left, const auto& right){return left.Channel() < right.Channel();});

        event.put(std::move(wireCollection), fOutInstanceLabelVec[instanceIdx]);
    
        if (fOutputRawWaveform)
        {
            // Copy the raw digits from the concurrent vector to our output vector
            RawDigitCollectionPtr rawRawDigitCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(concurrentRawRawDigits.begin()), 
                                                                                                       std::move_iterator(concurrentRawRawDigits.end()));
    
            // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
            std::sort(rawRawDigitCollection->begin(),rawRawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});
    
            // Now transfer ownership to the event store
            event.put(std::move(rawRawDigitCollection),fOutInstanceLabelVec[instanceIdx] + fOutputRawWavePath);
        }
    
        if (fOutputCorrection)
        {
            // Copy the raw digits from the concurrent vector to our output vector
            RawDigitCollectionPtr coherentCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(coherentRawDigits.begin()), 
                                                                                                    std::move_iterator(coherentRawDigits.end()));
    
            // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
            std::sort(coherentCollection->begin(),coherentCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});
    
            // Now transfer ownership to the event store
            event.put(std::move(coherentCollection),fOutInstanceLabelVec[instanceIdx] + fOutputCoherentPath);
        }

        instanceIdx++;
    }

    theClockTotal.stop();

    double totalTime = theClockTotal.accumulated_real_time();

    mf::LogInfo(fLogCategory) << "==> MCDecoderICARUSTPCwROI total time: " << totalTime << std::endl;

    return;
}

void MCDecoderICARUSTPCwROI::processSingleLabel(art::Event&                        event,
                                                const art::InputTag&               inputLabel,
                                                detinfo::DetectorClocksData const& clockData,
                                                ChannelArrayPairVec         const& channelArrayPairVec,
                                                size_t                      const& coherentNoiseGrouping,
                                                ConcurrentRawDigitCol&             concurrentRawDigits,
                                                ConcurrentRawDigitCol&             concurrentRawRawDigits,
                                                ConcurrentRawDigitCol&             coherentRawDigits,
                                                ConcurrentWireCol&                 concurrentROIs) const
{
    cet::cpu_timer theClockProcess;

    theClockProcess.start();

    // Read in the digit List object(s).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    event.getByLabel(inputLabel, digitVecHandle);

    // Require a valid handle
    if (digitVecHandle.isValid() && digitVecHandle->size()>0 )
    {
        // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
        // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
        std::vector<const raw::RawDigit*> rawDigitVec;

        // Ugliness to fill the pointer vector...
        for(size_t idx = 0; idx < digitVecHandle->size(); idx++) rawDigitVec.emplace_back(&digitVecHandle->at(idx)); //art::Ptr<raw::RawDigit>(digitVecHandle, idx).get());

        // Sort (use a lambda to sort by channel id)
        std::sort(rawDigitVec.begin(),rawDigitVec.end(),[](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});

        // Declare a temporary digit holder and resize it if downsizing the waveform
        unsigned int               dataSize = art::Ptr<raw::RawDigit>(digitVecHandle,0)->Samples(); //size of raw data vectors
        raw::RawDigit::ADCvector_t rawDataVec(dataSize);

        using BoardToChannelArrayPairMap = std::map<unsigned int, ChannelArrayPair>;

        BoardToChannelArrayPairMap  boardToChannelArrayPairMap;
        std::map<unsigned int, int> boardWireCountMap;
        const unsigned int          MAXCHANNELS(64);

        // Commence looping over raw digits
        for(const auto& rawDigit : rawDigitVec)
        {
            raw::ChannelID_t channel = rawDigit->Channel();

            ChannelToBoardWirePlaneMap::const_iterator channelToBoardItr = fChannelToBoardWirePlaneMap.find(channel);

            if (channelToBoardItr == fChannelToBoardWirePlaneMap.end())
            {
                std::cout << "********************************************************************************" << std::endl;
                std::cout << "********* We did not find channel " << channel << "*****************************" << std::endl;
                std::cout << "********************************************************************************" << std::endl;
                continue;
            }

            unsigned int readoutBoardID = channelToBoardItr->second.first;
            unsigned int wireIdx        = channelToBoardItr->second.second.first;
            unsigned int planeIdx       = channelToBoardItr->second.second.second;

            BoardToChannelArrayPairMap::iterator boardMapItr = boardToChannelArrayPairMap.find(readoutBoardID);

            if (boardMapItr == boardToChannelArrayPairMap.end())
            {
                const auto [mapItr, success] = 
                    boardToChannelArrayPairMap.insert({readoutBoardID,{daq::INoiseFilter::ChannelPlaneVec(MAXCHANNELS,{0,3}),icarus_signal_processing::ArrayFloat(MAXCHANNELS,icarus_signal_processing::VectorFloat(dataSize))}});

                if (!success) 
                {
                    std::cout << "+++> failed to insert data structure! " << std::endl;
                    continue;
                }

                boardMapItr = mapItr;
                boardWireCountMap[readoutBoardID] = 0;
            }

           // Decompress data into local holder
            raw::Uncompress(rawDigit->ADCs(), rawDataVec, rawDigit->Compression());

            // Fill into the data structure
            icarus_signal_processing::VectorFloat& boardDataVec = boardMapItr->second.second[wireIdx];

            for(size_t tick = 0; tick < dataSize; tick++) boardDataVec[tick] = rawDataVec[tick];

            boardMapItr->second.first[wireIdx] = daq::INoiseFilter::ChannelPlanePair(channel,planeIdx);

            if (++boardWireCountMap[readoutBoardID] == MAXCHANNELS)
            {
                processSingleImage(clockData, boardMapItr->second, coherentNoiseGrouping, concurrentRawDigits, concurrentRawRawDigits, coherentRawDigits, concurrentROIs);

                boardToChannelArrayPairMap.erase(boardMapItr);

                // Get the wireIDs for this channel
//                std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
//
//                // Since we work with channels we can ignore the case in the middle of the TPC where there are 
//                // wires crossing the midplane and just work with the first wire ID
//                const geo::WireID&  wireID =  wids[0];
//                const geo::PlaneID& planeID = wireID.planeID();
//
//                // Ok, now store things...
//                unsigned int planeIndex = fPlaneToROPPlaneMap.find(planeID)->second;
//                unsigned int wire       = channel - fPlaneToWireOffsetMap.find(planeID)->second;
//
//                if (wire >= channelArrayPairVec[planeIndex].second.size()) continue;
//
//                icarus_signal_processing::VectorFloat& dataVec = channelArrayPairVec[planeIndex].second[wire];
//
//                for(size_t tick = 0; tick < dataSize; tick++) dataVec[tick] = rawDataVec[tick];
//
//                // Keep track of the channel
//                channelArrayPairVec[planeIndex].first[wire] = daq::INoiseFilter::ChannelPlanePair(channel,planeID.Plane);
            }
        }

        // Some detector simulations don't output channels that don't have any possibility of signal (ghost channels)
        // Do a cleanup phase here to find these
        for(auto& boardInfo : boardToChannelArrayPairMap)
        {
            if (boardWireCountMap[boardInfo.first] < 64)
            {
                processSingleImage(clockData, boardInfo.second, boardWireCountMap[boardInfo.first], concurrentRawDigits, concurrentRawRawDigits, coherentRawDigits, concurrentROIs);
            }
        }

    }

    theClockProcess.stop();

    double totalTime = theClockProcess.accumulated_real_time();

    mf::LogDebug(fLogCategory) << "--> Exiting fragment processing for thread: " << tbb::this_task_arena::current_thread_index() << ", time: " << totalTime << std::endl;

    return;
}

void MCDecoderICARUSTPCwROI::processSingleImage(const detinfo::DetectorClocksData& clockData,
                                                const ChannelArrayPair&            channelArrayPair,
                                                size_t                             coherentNoiseGrouping,
                                                ConcurrentRawDigitCol&             concurrentRawDigitCol,
                                                ConcurrentRawDigitCol&             concurrentRawRawDigitCol,
                                                ConcurrentRawDigitCol&             coherentRawDigitCol,
                                                ConcurrentWireCol&                 concurrentROIs) const
{
    // Let's go through and fill the output vector
    const daq::INoiseFilter::ChannelPlaneVec&   channelVec = channelArrayPair.first;
    const icarus_signal_processing::ArrayFloat& dataArray  = channelArrayPair.second;

    unsigned int numChannels = dataArray.size();
    unsigned int numTicks    = dataArray[0].size();

    // Recover pointer to the decoder needed here
    INoiseFilter* decoderTool = fDecoderToolVec[tbb::this_task_arena::current_thread_index()].get();

    //process_fragment(event, rawfrag, product_collection, header_collection);
    decoderTool->process_fragment(clockData, channelVec, dataArray, coherentNoiseGrouping);

    // Now set up for output, we need to convert back from float to short int so use this
    raw::RawDigit::ADCvector_t wvfm(numTicks);

    // Loop over the channels to recover the RawDigits after filtering
    for(size_t chanIdx = 0; chanIdx < numChannels; chanIdx++)
    {
        // Skip if no channel data (plane is wrong)
        if (channelVec[chanIdx].second > 2) continue;
        
        raw::ChannelID_t channel = channelVec[chanIdx].first;

        if (fOutputRawWaveform)
        {
            const icarus_signal_processing::VectorFloat& waveform = decoderTool->getPedCorWaveforms()[chanIdx];

            // Need to convert from float to short int
            std::transform(waveform.begin(),waveform.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});
 
            ConcurrentRawDigitCol::iterator newRawObjItr = concurrentRawRawDigitCol.emplace_back(channel,wvfm.size(),wvfm); 

            newRawObjItr->SetPedestal(decoderTool->getPedestalVals()[chanIdx],decoderTool->getFullRMSVals()[chanIdx]);
        }

        if (fOutputCorrection)
        {
            const icarus_signal_processing::VectorFloat& corrections = decoderTool->getCorrectedMedians()[chanIdx];

            // Need to convert from float to short int
            std::transform(corrections.begin(),corrections.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});

            //ConcurrentRawDigitCol::iterator newRawObjItr = coherentRawDigitCol.emplace_back(channel,wvfm.size(),wvfm); 
            ConcurrentRawDigitCol::iterator newRawObjItr = coherentRawDigitCol.push_back(raw::RawDigit(channel,wvfm.size(),wvfm)); 

            newRawObjItr->SetPedestal(0.,0.);
        }

        // Recover the denoised waveform
        const icarus_signal_processing::VectorFloat& denoised = decoderTool->getWaveLessCoherent()[chanIdx];

        // Need to convert from float to short int
        std::transform(denoised.begin(),denoised.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});

        ConcurrentRawDigitCol::iterator newObjItr = concurrentRawDigitCol.emplace_back(channel,wvfm.size(),wvfm); 

        newObjItr->SetPedestal(0.,decoderTool->getTruncRMSVals()[chanIdx]);

        // And, finally, the ROIs 
        const icarus_signal_processing::VectorBool& chanROIs = decoderTool->getROIVals()[chanIdx];
        recob::Wire::RegionsOfInterest_t            ROIVec;

        // Go through candidate ROIs and create Wire ROIs
        size_t roiIdx = 0;

        while(roiIdx < chanROIs.size())
        {
            size_t roiStartIdx = roiIdx;

            while(roiIdx < chanROIs.size() && chanROIs[roiIdx]) roiIdx++;

            if (roiIdx > roiStartIdx)
            {
                std::vector<float> holder(roiIdx - roiStartIdx, 10.);

                ROIVec.add_range(roiStartIdx, std::move(holder));
            }

            roiIdx++;
        }

        concurrentROIs.push_back(recob::WireCreator(std::move(ROIVec),channel,fGeometry->View(channel)).move());
    }//loop over channel indices

    return;
}

void MCDecoderICARUSTPCwROI::saveRawDigits(const icarus_signal_processing::ArrayFloat&  dataArray, 
                                           const icarus_signal_processing::VectorFloat& pedestalVec,
                                           const icarus_signal_processing::VectorFloat& rmsVec,
                                           const icarus_signal_processing::VectorInt&   channelVec,
                                           ConcurrentRawDigitCol&                       rawDigitCol) const
{
    if (!dataArray.empty())
    {
        cet::cpu_timer theClockSave;

        theClockSave.start();

        raw::RawDigit::ADCvector_t wvfm(dataArray[0].size());

        mf::LogDebug(fLogCategory) << "    --> saving rawdigits for " << dataArray.size() << " channels" << std::endl;

        // Loop over the channels to recover the RawDigits after filtering
        for(size_t chanIdx = 0; chanIdx != dataArray.size(); chanIdx++)
        {
            // Protect against case where there was no readout 
            if (channelVec[chanIdx] < 0) continue;

            const icarus_signal_processing::VectorFloat& dataVec = dataArray[chanIdx];

            // Need to convert from float to short int
            std::transform(dataVec.begin(),dataVec.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});

            ConcurrentRawDigitCol::iterator newObjItr = rawDigitCol.emplace_back(channelVec[chanIdx],wvfm.size(),wvfm); 
            newObjItr->SetPedestal(pedestalVec[chanIdx],rmsVec[chanIdx]);
        }//loop over channel indices

        theClockSave.stop();

        double totalTime = theClockSave.accumulated_real_time();

        mf::LogDebug(fLogCategory) << "    --> done with save, time: " << totalTime << std::endl;
    }

    return;
}

//----------------------------------------------------------------------------
/// End job method.
void MCDecoderICARUSTPCwROI::endJob(art::ProcessingFrame const&)
{
    mf::LogInfo(fLogCategory) << "Looked at " << fNumEvent << " events" << std::endl;
}

} // end of namespace
