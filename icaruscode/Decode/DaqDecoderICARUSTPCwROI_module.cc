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

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoderFilter.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

#include "icarus_signal_processing/ICARUSSigProcDefs.h"

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
    using RawDigitCollection    = std::vector<raw::RawDigit>;
    using RawDigitCollectionPtr = std::unique_ptr<RawDigitCollection>;
    using ConcurrentRawDigitCol = tbb::concurrent_vector<raw::RawDigit>;

    // Define data structures for organizing the decoded fragments
    // The idea is to form complete "images" organized by "logical" TPC. Here we are including
    // both the active AND inactive channels so the noise processing will be correct
    // We build two lists here, first is the mapping to the actual image, 
    // the second will be a mapping to channel IDs where we assume the
    // order is the same between the two
    using PlaneIDtoImagePair   = std::pair<geo::PlaneID,icarus_signal_processing::ArrayFloat>;
    using PlaneIDtoImageMap    = std::map<geo::PlaneID,icarus_signal_processing::ArrayFloat>;
    using ChannelVec           = std::vector<raw::ChannelID_t>;
    using PlaneIDtoChannelPair = std::pair<geo::PlaneID,ChannelVec>;
    using PlaneIDtoChannelMap  = std::map<geo::PlaneID,ChannelVec>;

    // Function to do the work
    void processSingleFragment(size_t,
                               detinfo::DetectorClocksData const& clockData,
                               art::Handle<artdaq::Fragments>, 
                               PlaneIDtoImageMap&, 
                               PlaneIDtoChannelMap&) const;

private:
    class multiThreadFragmentProcessing
    {
    public:
        multiThreadFragmentProcessing(DaqDecoderICARUSTPCwROI const&     parent,
                                      detinfo::DetectorClocksData const& clockData,
                                      art::Handle<artdaq::Fragments>&    fragmentsHandle,
                                      PlaneIDtoImageMap&                 planeIDtoImageMap,
                                      PlaneIDtoChannelMap&               planeIDtoChannelMap)
            : fDaqDecoderICARUSTPCwROI(parent),
              fClockData{clockData},
              fFragmentsHandle(fragmentsHandle),
              fPlaneIDtoImageMap(planeIDtoImageMap),
              fPlaneIDtoChannelMap(planeIDtoChannelMap)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
              fDaqDecoderICARUSTPCwROI.processSingleFragment(idx, fClockData, fFragmentsHandle, fPlaneIDtoImageMap, fPlaneIDtoChannelMap);
        }
    private:
        const DaqDecoderICARUSTPCwROI&     fDaqDecoderICARUSTPCwROI;
        detinfo::DetectorClocksData const& fClockData;
        art::Handle<artdaq::Fragments>&    fFragmentsHandle;
        PlaneIDtoImageMap&                 fPlaneIDtoImageMap;
        PlaneIDtoChannelMap&               fPlaneIDtoChannelMap;
    };

    // Function to save our RawDigits
    void saveRawDigits(const icarus_signal_processing::ArrayFloat&, 
                       const icarus_signal_processing::VectorFloat&, 
                       const icarus_signal_processing::VectorFloat&,
                       const icarus_signal_processing::VectorInt&,
                       ConcurrentRawDigitCol&) const;

    // Fcl parameters.
    std::vector<art::InputTag>                   fFragmentsLabelVec;   ///< The input artdaq fragment label vector (for more than one)
    bool                                         fOutputRawWaveform;   ///< Should we output pedestal corrected (not noise filtered)?
    bool                                         fOutputCorrection;    ///< Should we output the coherent noise correction vectors?
    std::string                                  fOutputRawWavePath;   ///< Path to assign to the output if asked for
    std::string                                  fOutputCoherentPath;  ///< Path to assign to the output if asked for
    bool                                         fDiagnosticOutput;    ///< Set this to get lots of messages

    // Statistics.
    int                                          fNumEvent;             ///< Number of events seen.

    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*                     fGeometry;             ///< pointer to Geometry service
    const icarusDB::IICARUSChannelMap*           fChannelMap;
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
                          art::ReplicatedProducer(pset, frame),
                          fNumEvent(0)
{
    fGeometry   = art::ServiceHandle<geo::Geometry const>{}.get();
    fChannelMap = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    configure(pset);

    // Check the concurrency 
    int max_concurrency = tbb::this_task_arena::max_concurrency();

    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "     ==> concurrency: " << max_concurrency << std::endl;

    // Set up our "produces" 
    // Note that we can have multiple instances input to the module
    // Our convention will be to create a similar number of outputs with the same instance names
    for(const auto& fragmentLabel : fFragmentsLabelVec)
    {
        produces<std::vector<raw::RawDigit>>(fragmentLabel.instance());

        if (fOutputRawWaveform)
            produces<std::vector<raw::RawDigit>>(fragmentLabel.instance() + fOutputRawWavePath);

        if (fOutputCorrection)
            produces<std::vector<raw::RawDigit>>(fragmentLabel.instance() + fOutputCoherentPath);
    }

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
    fFragmentsLabelVec  = pset.get<std::vector<art::InputTag>>("FragmentsLabelVec",  std::vector<art::InputTag>()={"daq:PHYSCRATEDATA"});
    fOutputRawWaveform  = pset.get<bool                      >("OutputRawWaveform",                                               false);
    fOutputCorrection   = pset.get<bool                      >("OutputCorrection",                                                false);
    fOutputRawWavePath  = pset.get<std::string               >("OutputRawWavePath",                                               "raw");
    fOutputCoherentPath = pset.get<std::string               >("OutputCoherentPath",                                              "Cor");
    fDiagnosticOutput   = pset.get<bool                      >("DiagnosticOutput",                                                false);
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

    // Check the concurrency 
    int max_concurrency = tbb::this_task_arena::max_concurrency();

    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "     ==> concurrency: " << max_concurrency << std::endl;


    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "===> Run: " << event.id().run() << ", subrn: " << event.id().subRun() << ", event: " << event.id().event() << std::endl;

    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    // Loop through the list of input daq fragment collections one by one 
    // We are not trying to multi thread at this stage because we are trying to control
    // overall memory usage at this level. We'll multi thread internally...
    for(const auto& fragmentLabel : fFragmentsLabelVec)
    {
        art::Handle<artdaq::Fragments> daq_handle;
        event.getByLabel(fragmentLabel, daq_handle);

        ConcurrentRawDigitCol concurrentRawDigits;
        ConcurrentRawDigitCol concurrentRawRawDigits;
        ConcurrentRawDigitCol coherentRawDigits;

        PlaneIDtoImageMap   planeIDtoImageMap;
        PlaneIDtoChannelMap planeIDtoChannelMap;

        std::cout << "****> Let's get ready to rumble!" << std::endl;
    
        // ... Launch multiple threads with TBB to do the deconvolution and find ROIs in parallel
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);
        multiThreadFragmentProcessing fragmentProcessing(*this,
                                                         clockData,
                                                         daq_handle,
                                                         planeIDtoImageMap,
                                                         planeIDtoChannelMap);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, daq_handle->size()), fragmentProcessing);
    
        // Copy the raw digits from the concurrent vector to our output vector
        RawDigitCollectionPtr rawDigitCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(concurrentRawDigits.begin()), 
                                                                                                std::move_iterator(concurrentRawDigits.end()));
    
        // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
        std::sort(rawDigitCollection->begin(),rawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});

        // What did we get back?
        std::cout << "****> Total size of map: " << planeIDtoImageMap.size() << std::endl;
        for(const auto& planeImagePair : planeIDtoImageMap)
        {
            std::cout << "      - plane: " << planeImagePair.first << " has " << planeImagePair.second.size() << " wires" << std::endl;
        }
    
        // Now transfer ownership to the event store
        event.put(std::move(rawDigitCollection), fragmentLabel.instance());
    
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
    }

    theClockTotal.stop();

    double totalTime = theClockTotal.accumulated_real_time();

    mf::LogInfo("DaqDecoderICARUSTPCwROI") << "==> DaqDecoderICARUSTPCwROI total time: " << totalTime << std::endl;

    return;
}

void DaqDecoderICARUSTPCwROI::processSingleFragment(size_t                             idx,
                                                    detinfo::DetectorClocksData const& clockData,
                                                    art::Handle<artdaq::Fragments>     fragmentHandle,
                                                    PlaneIDtoImageMap&                 planeIDtoImageMap, 
                                                    PlaneIDtoChannelMap&               planeIDtoChannelMap) const
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

    if (fDiagnosticOutput) std::cout << "==> Recovered fragmentID: " << std::hex << fragmentID << std::dec << std::endl;

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
            if (fDiagnosticOutput)
            {
                std::cout << "*** COULD NOT FIND BOARD ***" << std::endl;
                std::cout << "    - boardID: " << std::hex << boardID << ", board map size: " << readoutIDVec.size() << ", nBoardsPerFragment: " << nBoardsPerFragment << std::endl;
            }

            return;
        }

        unsigned int boardSlot = fChannelMap->getBoardSlot(boardID);

        boardIDVec[boardSlot] = boardID;
    }

    if (fDiagnosticOutput)
    {
        std::cout << "   - # boards: " << boardIDVec.size() << ", boards: ";
        for(const auto& id : boardIDVec) std::cout << id << " ";
        std::cout << std::endl;
    }

    cet::cpu_timer theClockPedestal;

    theClockPedestal.start();

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

        if (fDiagnosticOutput)
        {
            std::cout << "********************************************************************************" << std::endl;
            std::cout << "FragmentID: " << std::hex << fragmentID << ", Crate: " << crateName << std::dec << ", boardID: " << boardSlot << "/" << nBoardsPerFragment << ", size " << channelPlanePairVec.size() << "/" << nChannelsPerBoard << std::endl;
        }

        // Get the pointer to the start of this board's block of data
        const icarus::A2795DataBlock::data_t* dataBlock = physCrateFragment.BoardData(board);

        // Copy to input data array
        for(size_t chanIdx = 0; chanIdx < nChannelsPerBoard; chanIdx++)
        {
            // Get the channel number on the Fragment
            raw::ChannelID_t channel = channelPlanePairVec[chanIdx].first;
            
            std::vector<geo::WireID> wireIDVec = fGeometry->ChannelToWire(channel);

            // Skip the channels which are not physically connected
            if (wireIDVec.empty()) continue;

            // Some contortions here... the horizontal wires are split and so independent,
            // the angled wires are logically split but we do want them put back together here... 
            const geo::PlaneID& planeID = wireIDVec[0].planeID();

            PlaneIDtoImageMap::iterator planeToImageIterator = planeIDtoImageMap.find(planeID);

            if (planeToImageIterator == planeIDtoImageMap.end())
            {
                int nWires = fGeometry->Nwires(planeID);

                readout::ROPID ropID = fGeometry->ChannelToROP(channel);

                std::cout << "--> ROPID: " << ropID << ", planeID: " << planeID << std::endl;
                
                planeIDtoImageMap.insert(PlaneIDtoImagePair(planeID,icarus_signal_processing::ArrayFloat(nWires,icarus_signal_processing::VectorFloat(nSamplesPerChannel)))).first;
                planeIDtoChannelMap[planeID] = ChannelVec(nWires);
            }

            unsigned int wire = wireIDVec[0].Wire;

            icarus_signal_processing::VectorFloat& rawDataVec = planeIDtoImageMap[planeID][wire];

            for(size_t tick = 0; tick < nSamplesPerChannel; tick++)
                rawDataVec[tick] = -dataBlock[chanIdx + tick * nChannelsPerBoard];

            // Keep track of the channel
            planeIDtoChannelMap[planeID][wire] = channelPlanePairVec[chanIdx].first;

            if (fDiagnosticOutput)
                std::cout << channelPlanePairVec[chanIdx].first << "-" << wireIDVec[0].Cryostat << "/" << wireIDVec[0].TPC << "/" << wireIDVec[0].Plane << "/" << wireIDVec[0].Wire  << " * ";
        }
        if (fDiagnosticOutput) std::cout << std::endl;
    }

    // We need to make sure the channelID information is not preserved when less than 9 boards in the fragment
//    if (nBoardsPerFragment < 9)
//    {
//        std::fill(fChannelIDVec.begin() + nBoardsPerFragment * nChannelsPerBoard, fChannelIDVec.end(), -1);
//    }


    theClockProcess.stop();

    double totalTime = theClockProcess.accumulated_real_time();

    mf::LogDebug("DaqDecoderICARUSTPCwROI") << "--> Exiting fragment processing for thread: " << tbb::this_task_arena::current_thread_index() << ", time: " << totalTime << std::endl;
    return;
}

void DaqDecoderICARUSTPCwROI::saveRawDigits(const icarus_signal_processing::ArrayFloat&  dataArray, 
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

        mf::LogDebug("DaqDecoderICARUSTPCwROI") << "    --> saving rawdigits for " << dataArray.size() << " channels" << std::endl;

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

        mf::LogDebug("DaqDecoderICARUSTPCwROI") << "    --> done with save, time: " << totalTime << std::endl;
    }

    return;
}

//----------------------------------------------------------------------------
/// End job method.
void DaqDecoderICARUSTPCwROI::endJob(art::ProcessingFrame const&)
{
    mf::LogInfo("DaqDecoderICARUSTPCwROI") << "Looked at " << fNumEvent << " events" << std::endl;
}

} // end of namespace
