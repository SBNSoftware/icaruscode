////////////////////////////////////////////////////////////////////////
//
// Class:       DaqDecoderICARUSTPCwROI
// Module Type: producer
// File:        DaqDecoderICARUSTPCwROI_module.cc
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
#include "sbnobj/ICARUS/TPC/ChannelROI.h"
#include "icaruscode/TPC/Utilities/ChannelROICreator.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Utilities/ArtHandleTrackerManager.h"
#include "icaruscode/Decode/DecoderTools/INoiseFilter.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

#include "icarus_signal_processing/ICARUSSigProcDefs.h"
#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Filters/FFTFilterFunctions.h"
#include "icarus_signal_processing/Filters/ImageFilters.h"
#include "icarus_signal_processing/Denoising.h"
#include "icarus_signal_processing/Detection/EdgeDetection.h"
#include "icarus_signal_processing/Filters/BilateralFilters.h"
#include "icarus_signal_processing/ROIFinder2D.h"

namespace daq 
{

class DaqDecoderICARUSTPCwROI : public art::ReplicatedProducer
{
public:

    // Copnstructors, destructor.
    explicit DaqDecoderICARUSTPCwROI(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame);
    virtual ~DaqDecoderICARUSTPCwROI();

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e, art::ProcessingFrame const& frame);
    virtual void beginJob(art::ProcessingFrame const& frame);
    virtual void endJob(art::ProcessingFrame const& frame);

    // Define the RawDigit collection
    using RawDigitCollection      = std::vector<raw::RawDigit>;
    using RawDigitCollectionPtr   = std::unique_ptr<RawDigitCollection>;
    using ChannelROICollection    = std::vector<recob::ChannelROI>;
    using ChannelROICollectionPtr = std::unique_ptr<ChannelROICollection>;
    using ConcurrentRawDigitCol   = tbb::concurrent_vector<raw::RawDigit>;
    using ConcurrentChannelROICol = tbb::concurrent_vector<recob::ChannelROI>;

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
    void processSingleFragment(size_t,
                               detinfo::DetectorClocksData const& clockData,
                               art::Handle<artdaq::Fragments>, 
                               ConcurrentRawDigitCol&,
                               ConcurrentRawDigitCol&,
                               ConcurrentRawDigitCol&,
                               ConcurrentRawDigitCol&,
                               ConcurrentChannelROICol&) const;

private:
    class multiThreadFragmentProcessing
    {
    public:
        multiThreadFragmentProcessing(DaqDecoderICARUSTPCwROI const&        parent,
                                      detinfo::DetectorClocksData const&    clockData,
                                      art::Handle<artdaq::Fragments> const& fragmentsHandle,
                                      ConcurrentRawDigitCol&                concurrentRawRawDigits,
                                      ConcurrentRawDigitCol&                concurrentRawDigits,
                                      ConcurrentRawDigitCol&                coherentRawDigits,
                                      ConcurrentRawDigitCol&                morphedRawDigits,
                                      ConcurrentChannelROICol&              concurrentROIs)
            : fDaqDecoderICARUSTPCwROI(parent),
              fClockData{clockData},
              fFragmentsHandle(fragmentsHandle),
              fConcurrentRawRawDigits(concurrentRawRawDigits),
              fConcurrentRawDigits(concurrentRawDigits),
              fCoherentRawDigits(coherentRawDigits),
              fMorphedRawDigits(morphedRawDigits),
              fConcurrentROIs(concurrentROIs)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
              fDaqDecoderICARUSTPCwROI.processSingleFragment(idx, fClockData, fFragmentsHandle, fConcurrentRawRawDigits, fConcurrentRawDigits, fCoherentRawDigits, fMorphedRawDigits, fConcurrentROIs);
        }
    private:
        const DaqDecoderICARUSTPCwROI&        fDaqDecoderICARUSTPCwROI;
        detinfo::DetectorClocksData const&    fClockData;
        art::Handle<artdaq::Fragments> const& fFragmentsHandle;
        ConcurrentRawDigitCol&                fConcurrentRawRawDigits;
        ConcurrentRawDigitCol&                fConcurrentRawDigits;
        ConcurrentRawDigitCol&                fCoherentRawDigits;
        ConcurrentRawDigitCol&                fMorphedRawDigits;
        ConcurrentChannelROICol&              fConcurrentROIs;
    };

    // Function to save our RawDigits
    void saveRawDigits(const icarus_signal_processing::ArrayFloat&, 
                       const icarus_signal_processing::VectorFloat&, 
                       const icarus_signal_processing::VectorFloat&,
                       const icarus_signal_processing::VectorInt&,
                       ConcurrentRawDigitCol&) const;

    // Fcl parameters.
    std::vector<art::InputTag>                                  fFragmentsLabelVec;          ///< The input artdaq fragment label vector (for more than one)
    bool                                                        fOutputRawWaveform;          ///< Should we output pedestal corrected (not noise filtered)?
    bool                                                        fOutputCorrection;           ///< Should we output the coherent noise correction vectors?
    bool                                                        fOutputMorphed;              ///< Should we output the morphological filter vectors?
    std::string                                                 fOutputRawWavePath;          ///< Path to assign to the output if asked for
    std::string                                                 fOutputCoherentPath;         ///< Path to assign to the output if asked for
    std::string                                                 fOutputMorphedPath;          ///< Path to assign to the output if asked for
    bool                                                        fDiagnosticOutput;           ///< Set this to get lots of messages
    float                                                       fSigmaForTruncation;         ///< Cut for truncated rms calc
    size_t                                                      fCoherentNoiseGrouping;      ///< Grouping for removing coherent noise

    bool fDropRawDataAfterUse;   ///< Clear fragment data product cache after use.
  
    const std::string                                           fLogCategory;                ///< Output category when logging messages

    // We need to give to the denoiser the "threshold vector" we will fill during our data loop
    icarus_signal_processing::VectorFloat  fThresholdVec;  ///< "threshold vector" filled during decoding loop

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

    // Tools for decoding fragments depending on type
    std::vector<std::unique_ptr<INoiseFilter>>                  fDecoderToolVec;       ///< Decoder tools

    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*                                    fGeometry;             ///< pointer to Geometry service
    const icarusDB::IICARUSChannelMap*                          fChannelMap;
};

DEFINE_ART_MODULE(DaqDecoderICARUSTPCwROI)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
DaqDecoderICARUSTPCwROI::DaqDecoderICARUSTPCwROI(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame) :
                            art::ReplicatedProducer(pset, frame),fLogCategory("DaqDecoderICARUSTPCwROI"),fNumEvent(0), fNumROPs(0)
{
    fGeometry   = art::ServiceHandle<geo::Geometry const>{}.get();
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    configure(pset);

    // Check the concurrency 
    int max_concurrency = art::Globals::instance()->nthreads(); 

    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "     ==> concurrency: " << max_concurrency << std::endl;

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& decoderToolParams = pset.get<fhicl::ParameterSet>("DecoderTool");
    
    fDecoderToolVec.resize(max_concurrency);
   
    for(auto& decoderTool : fDecoderToolVec)
    {
        // Get instance of tool
        decoderTool = art::make_tool<INoiseFilter>(decoderToolParams);
    }

    // Set up our "producers" 
    // Note that we can have multiple instances input to the module
    // Our convention will be to create a similar number of outputs with the same instance names
    for(const auto& fragmentLabel : fFragmentsLabelVec)
    {
        produces<std::vector<raw::RawDigit>>(fragmentLabel.instance());
        produces<std::vector<recob::ChannelROI>>(fragmentLabel.instance());

        if (fOutputRawWaveform)
            produces<std::vector<raw::RawDigit>>(fragmentLabel.instance() + fOutputRawWavePath);

        if (fOutputCorrection)
            produces<std::vector<raw::RawDigit>>(fragmentLabel.instance() + fOutputCoherentPath);

        if (fOutputMorphed)
            produces<std::vector<raw::RawDigit>>(fragmentLabel.instance() + fOutputMorphedPath);
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

                raw::ChannelID_t channel = fGeometry->PlaneWireToChannel(geo::WireID(planeID, 0));

                readout::ROPID ropID = fGeometry->ChannelToROP(channel);

                fPlaneToROPPlaneMap[planeID]      = ropID.ROP;
                fPlaneToWireOffsetMap[planeID]    = channel;
                planeToLastWireOffsetMap[planeID] = fGeometry->PlaneWireToChannel(geo::WireID(planeID, fGeometry->Nwires(planeID)));
                fROPToNumWiresMap[ropID.ROP]      = fGeometry->Nwires(planeID);

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

    // Report.
    mf::LogInfo("DaqDecoderICARUSTPCwROI") << "DaqDecoderICARUSTPCwROI configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
DaqDecoderICARUSTPCwROI::~DaqDecoderICARUSTPCwROI()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void DaqDecoderICARUSTPCwROI::configure(fhicl::ParameterSet const & pset)
{
    fFragmentsLabelVec     = pset.get<std::vector<art::InputTag>>("FragmentsLabelVec",  std::vector<art::InputTag>()={"daq:PHYSCRATEDATA"});
    fOutputRawWaveform     = pset.get<bool                      >("OutputRawWaveform",                                               false);
    fOutputCorrection      = pset.get<bool                      >("OutputCorrection",                                                false);
    fOutputMorphed         = pset.get<bool                      >("OutputMorphed",                                                   false);
    fOutputRawWavePath     = pset.get<std::string               >("OutputRawWavePath",                                               "raw");
    fOutputCoherentPath    = pset.get<std::string               >("OutputCoherentPath",                                              "Cor");
    fOutputMorphedPath     = pset.get<std::string               >("OutputMorphedPath",                                               "Mor");
    fDiagnosticOutput      = pset.get<bool                      >("DiagnosticOutput",                                                false);
    fSigmaForTruncation    = pset.get<float                     >("NSigmaForTrucation",                                                3.5);
    fCoherentNoiseGrouping = pset.get<size_t                    >("CoherentGrouping",                                                   64);
    fDropRawDataAfterUse   = pset.get<bool                      >("DropRawDataAfterUse",                                              true);
}

//----------------------------------------------------------------------------
/// Begin job method.
void DaqDecoderICARUSTPCwROI::beginJob(art::ProcessingFrame const&)
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
void DaqDecoderICARUSTPCwROI::produce(art::Event & event, art::ProcessingFrame const&)
{
    ++fNumEvent;

    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "**** Processing raw data fragments ****" << std::endl;
    
    util::LocalArtHandleTrackerManager dataCacheRemover
        (event, fDropRawDataAfterUse);

    // Check the concurrency 
    int max_concurrency = tbb::this_task_arena::max_concurrency();

    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "     ==> concurrency: " << max_concurrency << std::endl;

    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    // Loop through the list of input daq fragment collections one by one 
    // We are not trying to multi thread at this stage because we are trying to control
    // overall memory usage at this level. We'll multi thread internally...
    for(const auto& fragmentLabel : fFragmentsLabelVec)
    {
        art::Handle<artdaq::Fragments> const& daq_handle
          = dataCacheRemover.getHandle<artdaq::Fragments>(fragmentLabel);

        ConcurrentRawDigitCol   concurrentRawDigits;
        ConcurrentRawDigitCol   concurrentRawRawDigits;
        ConcurrentRawDigitCol   coherentRawDigits;
        ConcurrentRawDigitCol   morphedRawDigits;
        ConcurrentChannelROICol concurrentROIs;

        PlaneIdxToImageMap   planeIdxToImageMap;
        PlaneIdxToChannelMap planeIdxToChannelMap;

        ChannelArrayPairVec  channelArrayPairVec(fNumROPs);

        // Because the arrays can be variable size we need to loop to initialize
        for(size_t ropIdx = 0; ropIdx < fNumROPs; ropIdx++)
        {
            ChannelArrayPair& channelArrayPair = channelArrayPairVec[ropIdx];

            channelArrayPair.first.resize(fROPToNumWiresMap[ropIdx]);
            channelArrayPair.second.resize(fROPToNumWiresMap[ropIdx],icarus_signal_processing::VectorFloat(4096));

            mf::LogDebug("DaqDecoderICARUSTPCwROI") << "**> Initializing ropIdx: " << ropIdx << " channelPairVec to " << channelArrayPair.first.size() << " channels with " << channelArrayPair.second[0].size() << " ticks" << std::endl;
        }

        mf::LogDebug("DaqDecoderICARUSTPCwROI") << "****> Let's get ready to rumble!" << std::endl;
    
        // ... Launch multiple threads with TBB to do the deconvolution and find ROIs in parallel
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);

        multiThreadFragmentProcessing fragmentProcessing(*this, clockData, daq_handle, concurrentRawRawDigits, concurrentRawDigits, coherentRawDigits, morphedRawDigits, concurrentROIs);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, daq_handle->size()), fragmentProcessing);

        // Now let's process the resulting images
    //    multiThreadImageProcessing imageProcessing(*this, clockData, channelArrayPairVec, concurrentRawDigits, coherentRawDigits, concurrentROIs);

    //    tbb::parallel_for(tbb::blocked_range<size_t>(0, fNumROPs), imageProcessing);
    
        // Copy the raw digits from the concurrent vector to our output vector
        RawDigitCollectionPtr rawDigitCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(concurrentRawDigits.begin()), 
                                                                                                std::move_iterator(concurrentRawDigits.end()));
    
        // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
        std::sort(rawDigitCollection->begin(),rawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});

        // What did we get back?
        mf::LogDebug("DaqDecoderICARUSTPCwROI") << "****> Total size of map: " << planeIdxToImageMap.size() << std::endl;
        for(const auto& planeImagePair : planeIdxToImageMap)
        {
            mf::LogDebug("DaqDecoderICARUSTPCwROI") << "      - plane: " << planeImagePair.first << " has " << planeImagePair.second.size() << " wires" << std::endl;
        }
    
        // Now transfer ownership to the event store
        event.put(std::move(rawDigitCollection), fragmentLabel.instance());

        // Do the same to output the candidate ROIs
        ChannelROICollectionPtr channelROICollection = std::make_unique<std::vector<recob::ChannelROI>>(std::move_iterator(concurrentROIs.begin()),
                                                                                                        std::move_iterator(concurrentROIs.end()));

        std::sort(channelROICollection->begin(),channelROICollection->end(),[](const auto& left, const auto& right){return left.Channel() < right.Channel();});

        event.put(std::move(channelROICollection), fragmentLabel.instance());
    
    
        if (fOutputRawWaveform)
        {
            // Copy the raw digits from the concurrent vector to our output vector
            RawDigitCollectionPtr rawRawDigitCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(concurrentRawRawDigits.begin()), 
                                                                                                       std::move_iterator(concurrentRawRawDigits.end()));
    
            // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
            std::sort(rawRawDigitCollection->begin(),rawRawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});
    
            // Now transfer ownership to the event store
            event.put(std::move(rawRawDigitCollection),fragmentLabel.instance() + fOutputRawWavePath);
        }
    
        if (fOutputCorrection)
        {
            // Copy the raw digits from the concurrent vector to our output vector
            RawDigitCollectionPtr coherentCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(coherentRawDigits.begin()), 
                                                                                                    std::move_iterator(coherentRawDigits.end()));
    
            // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
            std::sort(coherentCollection->begin(),coherentCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});
    
            // Now transfer ownership to the event store
            event.put(std::move(coherentCollection),fragmentLabel.instance() + fOutputCoherentPath);
        }
    
        if (fOutputMorphed)
        {
            // Copy the raw digits from the concurrent vector to our output vector
            RawDigitCollectionPtr morphedCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(morphedRawDigits.begin()), 
                                                                                                   std::move_iterator(morphedRawDigits.end()));
    
            // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
            std::sort(morphedCollection->begin(),morphedCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});
    
            // Now transfer ownership to the event store
            event.put(std::move(morphedCollection),fragmentLabel.instance() + fOutputMorphedPath);
        }
    }

    theClockTotal.stop();

    double totalTime = theClockTotal.accumulated_real_time();

    mf::LogInfo(fLogCategory) << "==> DaqDecoderICARUSTPCwROI total time: " << totalTime << std::endl;

    return;
}

void DaqDecoderICARUSTPCwROI::processSingleFragment(size_t                             idx,
                                                    detinfo::DetectorClocksData const& clockData,
                                                    art::Handle<artdaq::Fragments>     fragmentHandle,
                                                    ConcurrentRawDigitCol&             concurrentRawRawDigitCol,
                                                    ConcurrentRawDigitCol&             concurrentRawDigitCol,
                                                    ConcurrentRawDigitCol&             coherentRawDigitCol,
                                                    ConcurrentRawDigitCol&             morphedRawDigitCol,
                                                    ConcurrentChannelROICol&           concurrentROIs) const
{
    cet::cpu_timer theClockProcess;

    theClockProcess.start();

    art::Ptr<artdaq::Fragment> fragmentPtr(fragmentHandle, idx);

    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "--> Processing fragment ID: " << fragmentPtr->fragmentID() << std::endl;
    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "    ==> Current thread index: " << tbb::this_task_arena::current_thread_index() << std::endl;

    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment physCrateFragment(*fragmentPtr);

    size_t nBoardsPerFragment = physCrateFragment.nBoards();
    size_t nChannelsPerBoard  = physCrateFragment.nChannelsPerBoard();
    size_t nSamplesPerChannel = physCrateFragment.nSamplesPerChannel();
//    size_t nChannelsPerFragment = nBoardsPerFragment * nChannelsPerBoard;

    // Recover the Fragment id:
    artdaq::detail::RawFragmentHeader::fragment_id_t fragmentID = fragmentPtr->fragmentID();

    mf::LogDebug(fLogCategory) << "==> Recovered fragmentID: " << std::hex << fragmentID << std::dec << std::endl;

    // Look for special case of diagnostic running
    if (!fChannelMap->hasFragmentID(fragmentID))
    {
        return;
    }

    // Recover the crate name for this fragment
    const std::string& crateName = fChannelMap->getCrateName(fragmentID);

    // Get the board ids for this fragment
    const icarusDB::ReadoutIDVec& readoutIDVec = fChannelMap->getReadoutBoardVec(fragmentID);

    icarusDB::ReadoutIDVec boardIDVec(readoutIDVec.size());

    // Note we want these to be in "slot" order...
    for(const auto& boardID : readoutIDVec)
    {
        // Look up the channels associated to this board
        if (!fChannelMap->hasBoardID(boardID))
        {
            mf::LogDebug(fLogCategory) << "*** COULD NOT FIND BOARD ***\n" <<
                                          "    - boardID: " << std::hex << boardID << ", board map size: " << readoutIDVec.size() << ", nBoardsPerFragment: " << nBoardsPerFragment;

            return;
        }

        unsigned int boardSlot = fChannelMap->getBoardSlot(boardID);

        boardIDVec[boardSlot] = boardID;
    }

    std::string boardIDs = "";

    for(const auto& id : boardIDVec) boardIDs += std::to_string(id) + " ";

    mf::LogDebug(fLogCategory) << "   - # boards: " << boardIDVec.size() << ", boards: " << boardIDs;

    cet::cpu_timer theClockPedestal;

    theClockPedestal.start();

    // Recover pointer to the decoder needed here
    INoiseFilter* decoderTool = fDecoderToolVec[tbb::this_task_arena::current_thread_index()].get();

    // Create a local channel pair  to hold at most a boards worth of info (64 channels x 4096 ticks)
    ChannelArrayPair channelArrayPair;

    channelArrayPair.first.resize(nChannelsPerBoard);
    channelArrayPair.second.resize(nChannelsPerBoard,icarus_signal_processing::VectorFloat(nSamplesPerChannel));

    // Now set up for output, we need to convert back from float to short int so use this
    raw::RawDigit::ADCvector_t wvfm(nSamplesPerChannel);

    // The first task is to recover the data from the board data block, determine and subtract the pedestals
    // and store into vectors useful for the next steps
    for(size_t board = 0; board < boardIDVec.size(); board++)
    {
        // Some diagnostics test are removing boards so put in check here to watch for this
        if (board >= nBoardsPerFragment)
        {
            mf::LogInfo("TPCDecoderFilter1D") << " Asking for board beyond number in fragment, want board " << board << ", maximum is " << nBoardsPerFragment << std::endl;
            continue;
        }

        const icarusDB::ChannelPlanePairVec& channelPlanePairVec = fChannelMap->getChannelPlanePair(boardIDVec[board]);

        uint32_t boardSlot = physCrateFragment.DataTileHeader(board)->StatusReg_SlotID();

        // ** the line below removed as per request of the TPC hardware folks 
//        const icarusDB::ChannelPlanePairVec& channelPlanePairVec = fChannelMap->getChannelPlanePair(boardIDVec[boardSlot]);

        mf::LogDebug(fLogCategory) << "********************************************************************************\n"
                                   << "FragmentID: " << std::hex << fragmentID << std::dec << ", Crate: " << crateName << ", boardID: " << boardSlot << "/" << nBoardsPerFragment << ", size " << channelPlanePairVec.size() << "/" << nChannelsPerBoard;

        if (board != boardSlot)
        {
            mf::LogInfo(fLogCategory) << "==> Found board/boardSlot mismatch, crate: " << crateName << ", board: " << board << ", boardSlot: " << boardSlot << " channelPlanePair: " << fChannelMap->getChannelPlanePair(boardIDVec[board]).front().first << "/"  << fChannelMap->getChannelPlanePair(boardIDVec[board]).front().second << ", slot: " << channelPlanePairVec[0].first << "/" << channelPlanePairVec[0].second;
        }
        // Copy to input data array
        for(size_t chanIdx = 0; chanIdx < nChannelsPerBoard; chanIdx++)
        {
           icarus_signal_processing::VectorFloat& rawDataVec = channelArrayPair.second[chanIdx];
           for (size_t tick = 0; tick < nSamplesPerChannel; ++tick)
             rawDataVec[tick] = -physCrateFragment.adc_val(board, chanIdx, tick);

          // Keep track of the channel
          channelArrayPair.first[chanIdx] = channelPlanePairVec[chanIdx];
        }

        //process_fragment(event, rawfrag, product_collection, header_collection);
        decoderTool->process_fragment(clockData, channelArrayPair.first, channelArrayPair.second, fCoherentNoiseGrouping);

        // We need to recalculate pedestals for the noise corrected waveforms
        icarus_signal_processing::WaveformTools<float> waveformTools;

        // Local storage for recomputing the the pedestals for the noise corrected data
        float localPedestal(0.);
        float localFullRMS(0.);
        float localTruncRMS(0.);
        int   localNumTruncBins(0);
        int   localRangeBins(0);

        float sigmaCut(fSigmaForTruncation);

        // Recover the denoised waveform
        const icarus_signal_processing::ArrayFloat& denoised = decoderTool->getWaveLessCoherent();

        icarus_signal_processing::VectorFloat       pedCorWaveforms(denoised[0].size());

        for(size_t chanIdx = 0; chanIdx < nChannelsPerBoard; chanIdx++)
        {
            // Get the channel number on the Fragment
            raw::ChannelID_t channel = channelPlanePairVec[chanIdx].first;

            // Are we storing the raw waveforms?
            if (fOutputRawWaveform)
            {
                //const icarus_signal_processing::VectorFloat& waveform = decoderTool->getPedCorWaveforms()[chanIdx];
                const icarus_signal_processing::VectorFloat& waveform = decoderTool->getRawWaveforms()[chanIdx];

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

            if (fOutputMorphed)
            {
                const icarus_signal_processing::VectorFloat& corrections = decoderTool->getMorphedWaveforms()[chanIdx];

                // Need to convert from float to short int
                std::transform(corrections.begin(),corrections.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});

                //ConcurrentRawDigitCol::iterator newRawObjItr = coherentRawDigitCol.emplace_back(channel,wvfm.size(),wvfm); 
                ConcurrentRawDigitCol::iterator newRawObjItr = morphedRawDigitCol.push_back(raw::RawDigit(channel,wvfm.size(),wvfm)); 

                newRawObjItr->SetPedestal(0.,0.);
            }

            // Now determine the pedestal and correct for it
            waveformTools.getPedestalCorrectedWaveform(denoised[chanIdx],
                                                       pedCorWaveforms,
                                                       sigmaCut,
                                                       localPedestal,
                                                       localFullRMS,
                                                       localTruncRMS,
                                                       localNumTruncBins,
                                                       localRangeBins);

            // Need to convert from float to short int
            std::transform(pedCorWaveforms.begin(),pedCorWaveforms.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});

            ConcurrentRawDigitCol::iterator newObjItr = concurrentRawDigitCol.emplace_back(channel,wvfm.size(),wvfm); 

            newObjItr->SetPedestal(localPedestal,localFullRMS);

            // And, finally, the ROIs 
            const icarus_signal_processing::VectorBool& chanROIs = decoderTool->getROIVals()[chanIdx];
            recob::ChannelROI::RegionsOfInterest_t      ROIVec;

            // Go through candidate ROIs and create Wire ROIs
            size_t roiIdx = 0;

            while(roiIdx < chanROIs.size())
            {
                size_t roiStartIdx = roiIdx;

                while(roiIdx < chanROIs.size() && chanROIs[roiIdx]) roiIdx++;

                if (roiIdx > roiStartIdx)
                {
                    std::vector<short> holder(roiIdx - roiStartIdx);

                    for(size_t idx = 0; idx < holder.size(); idx++) holder[idx] = wvfm[roiStartIdx+idx];

                    ROIVec.add_range(roiStartIdx, std::move(holder));
                }

                roiIdx++;
            }
        
            concurrentROIs.push_back(recob::ChannelROICreator(std::move(ROIVec),channel).move());
        }
    }

    // We need to make sure the channelID information is not preserved when less than 9 boards in the fragment
//    if (nBoardsPerFragment < 9)
//    {
//        std::fill(fChannelIDVec.begin() + nBoardsPerFragment * nChannelsPerBoard, fChannelIDVec.end(), -1);
//    }


    theClockProcess.stop();

    double totalTime = theClockProcess.accumulated_real_time();

    mf::LogDebug(fLogCategory) << "--> Exiting fragment processing for thread: " << tbb::this_task_arena::current_thread_index() << ", time: " << totalTime << std::endl;
    return;
}

//----------------------------------------------------------------------------
/// End job method.
void DaqDecoderICARUSTPCwROI::endJob(art::ProcessingFrame const&)
{
    mf::LogInfo(fLogCategory) << "Looked at " << fNumEvent << " events" << std::endl;
}

} // end of namespace
