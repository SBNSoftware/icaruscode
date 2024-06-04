////////////////////////////////////////////////////////////////////////
//
// Class:       DaqDecoderICARUSTPC
// Module Type: producer
// File:        DaqDecoderICARUSTPC_module.cc
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

#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/RawDigit.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoderFilter.h"

#include "icarus_signal_processing/ICARUSSigProcDefs.h"
#include "icarus_signal_processing/WaveformTools.h"

namespace daq 
{

class DaqDecoderICARUSTPC : public art::ReplicatedProducer
{
public:

    // Copnstructors, destructor.
    explicit DaqDecoderICARUSTPC(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame);
    virtual ~DaqDecoderICARUSTPC();

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e, art::ProcessingFrame const& frame);
    virtual void beginJob(art::ProcessingFrame const& frame);
    virtual void endJob(art::ProcessingFrame const& frame);

    // Define the RawDigit collection
    using RawDigitCollection    = std::vector<raw::RawDigit>;
    using RawDigitCollectionPtr = std::unique_ptr<RawDigitCollection>;
    using ConcurrentRawDigitCol = tbb::concurrent_vector<raw::RawDigit>;

    // Function to do the work
    void processSingleFragment(size_t,
                               detinfo::DetectorClocksData const& clockData,
                               art::Handle<artdaq::Fragments>, ConcurrentRawDigitCol&, ConcurrentRawDigitCol&, ConcurrentRawDigitCol&) const;

private:

    class multiThreadFragmentProcessing
    {
    public:
        multiThreadFragmentProcessing(DaqDecoderICARUSTPC const&         parent,
                                      detinfo::DetectorClocksData const& clockData,
                                      art::Handle<artdaq::Fragments>&    fragmentsHandle,
                                      ConcurrentRawDigitCol&             rawDigitCollection,
                                      ConcurrentRawDigitCol&             rawRawDigitCollection,
                                      ConcurrentRawDigitCol&             coherentCollection)
            : fDaqDecoderICARUSTPC(parent),
              fClockData{clockData},
              fFragmentsHandle(fragmentsHandle),
              fRawDigitCollection(rawDigitCollection),
              fRawRawDigitCollection(rawRawDigitCollection),
              fCoherentCollection(coherentCollection)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
              fDaqDecoderICARUSTPC.processSingleFragment(idx, fClockData, fFragmentsHandle, fRawDigitCollection, fRawRawDigitCollection, fCoherentCollection);
        }
    private:
        const DaqDecoderICARUSTPC&         fDaqDecoderICARUSTPC;
        detinfo::DetectorClocksData const& fClockData;
        art::Handle<artdaq::Fragments>&    fFragmentsHandle;
        ConcurrentRawDigitCol&             fRawDigitCollection;
        ConcurrentRawDigitCol&             fRawRawDigitCollection;
        ConcurrentRawDigitCol&             fCoherentCollection;
    };

    // Function to save our RawDigits
    void saveRawDigits(const icarus_signal_processing::ArrayFloat&, 
                       const icarus_signal_processing::VectorFloat&, 
                       const icarus_signal_processing::VectorFloat&,
                       const icarus_signal_processing::VectorInt&,
                       ConcurrentRawDigitCol&) const;

    // Tools for decoding fragments depending on type
    std::vector<std::unique_ptr<IDecoderFilter>> fDecoderToolVec;      ///< Decoder tools

    // Fcl parameters.
    std::vector<art::InputTag>                   fFragmentsLabelVec;   ///< The input artdaq fragment label vector (for more than one)
    bool                                         fOutputRawWaveform;   ///< Should we output pedestal corrected (not noise filtered)?
    bool                                         fOutputCorrection;    ///< Should we output the coherent noise correction vectors?
    std::string                                  fOutputRawWavePath;   ///< Path to assign to the output if asked for
    std::string                                  fOutputCoherentPath;  ///< Path to assign to the output if asked for
    unsigned int                                 fPlaneToSimulate;     ///< Use to get fragment offset
    float                                        fSigmaForTruncation;  ///< This determines the point at which we truncate bins for the RMS calc

    // Statistics.
    int                                          fNumEvent;             ///< Number of events seen.

    size_t                                       fFragmentOffset;       ///< The fragment offset to set channel numbering

};

DEFINE_ART_MODULE(DaqDecoderICARUSTPC)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
DaqDecoderICARUSTPC::DaqDecoderICARUSTPC(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame) :
                      art::ReplicatedProducer(pset, frame),
                      fNumEvent(0)
{
    geo::WireReadoutGeom const& wireReadoutAlg = art::ServiceHandle<geo::WireReadout const>()->Get();

    configure(pset);

    // Check the concurrency 
    int max_concurrency = tbb::this_task_arena::max_concurrency();

    mf::LogDebug("DaqDecoderICARUSTPC") << "     ==> concurrency: " << max_concurrency << std::endl;

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& decoderToolParams = pset.get<fhicl::ParameterSet>("DecoderTool");
    
    fDecoderToolVec.resize(max_concurrency);
    
    for(auto& decoderTool : fDecoderToolVec)
    {
        // Get instance of tool
        decoderTool = art::make_tool<IDecoderFilter>(decoderToolParams);
    }

    // Compute the fragment offset from the channel number for the desired plane
    // Get a base channel number for the plane we want
    mf::LogDebug("DaqDecoderICARUSTPC") << "ICARUS has " << wireReadoutAlg.Nchannels() << " in total with " << wireReadoutAlg.Views().size() << " views" << std::endl;

    geo::WireID wireID(0, 0, fPlaneToSimulate, 0);

    mf::LogDebug("DaqDecoderICARUSTPC") << "WireID: " << wireID << std::endl;

    geo::WireID firstWireID = *wireReadoutAlg.begin<geo::WireID>(geo::PlaneID(0,0,fPlaneToSimulate));

    mf::LogDebug("DaqDecoderICARUSTPC") << "From geo, first WireID: " << firstWireID << std::endl;

    raw::ChannelID_t channel = wireReadoutAlg.PlaneWireToChannel(wireID);

    mf::LogDebug("DaqDecoderICARUSTPC") << "Channel: " << channel << std::endl;

    for(size_t thePlane = 0; thePlane < 3; thePlane++)
    {
        geo::WireID  tempWireID(0, 0, thePlane, 0);
        geo::PlaneID tempPlaneID = tempWireID.planeID();

        mf::LogDebug("DaqDecoderICARUSTPC") << "thePlane: " << thePlane << ", WireID: " << tempWireID << ", channel: " << 
        wireReadoutAlg.PlaneWireToChannel(tempWireID) << ", view: " << wireReadoutAlg.Plane(tempPlaneID).View() << std::endl;
    }

    fFragmentOffset = channel / 576;

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
    mf::LogInfo("DaqDecoderICARUSTPC") << "DaqDecoderICARUSTPC configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
DaqDecoderICARUSTPC::~DaqDecoderICARUSTPC()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void DaqDecoderICARUSTPC::configure(fhicl::ParameterSet const & pset)
{
    fFragmentsLabelVec  = pset.get<std::vector<art::InputTag>>("FragmentsLabelVec",  std::vector<art::InputTag>()={"daq:PHYSCRATEDATA"});
    fOutputRawWaveform  = pset.get<bool                      >("OutputRawWaveform",                                               false);
    fOutputCorrection   = pset.get<bool                      >("OutputCorrection",                                                false);
    fOutputRawWavePath  = pset.get<std::string               >("OutputRawWavePath",                                               "raw");
    fOutputCoherentPath = pset.get<std::string               >("OutputCoherentPath",                                              "Cor");
    fPlaneToSimulate    = pset.get<unsigned int              >("PlaneToSimulate",                                                     2);
    fSigmaForTruncation = pset.get<float                     >("NSigmaForTrucation",                                                3.5);
}

//----------------------------------------------------------------------------
/// Begin job method.
void DaqDecoderICARUSTPC::beginJob(art::ProcessingFrame const&)
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
void DaqDecoderICARUSTPC::produce(art::Event & event, art::ProcessingFrame const&)
{
    ++fNumEvent;

    mf::LogDebug("DaqDecoderICARUSTPC") << "**** Processing raw data fragments ****" << std::endl;

    // Check the concurrency 
    int max_concurrency = tbb::this_task_arena::max_concurrency();

    mf::LogDebug("DaqDecoderICARUSTPC") << "     ==> concurrency: " << max_concurrency << std::endl;


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
    
        // ... Launch multiple threads with TBB to do the deconvolution and find ROIs in parallel
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(event);
        multiThreadFragmentProcessing fragmentProcessing(*this,
                                                         clockData,
                                                         daq_handle,
                                                         concurrentRawDigits,
                                                         concurrentRawRawDigits,
                                                         coherentRawDigits);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, daq_handle->size()), fragmentProcessing);
    
        // Copy the raw digits from the concurrent vector to our output vector
        RawDigitCollectionPtr rawDigitCollection = std::make_unique<std::vector<raw::RawDigit>>(std::move_iterator(concurrentRawDigits.begin()), 
                                                                                                std::move_iterator(concurrentRawDigits.end()));
    
        // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
        std::sort(rawDigitCollection->begin(),rawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});
    
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

    mf::LogInfo("DaqDecoderICARUSTPC") << "==> DaqDecoderICARUSTPC total time: " << totalTime << std::endl;

    return;
}

void DaqDecoderICARUSTPC::processSingleFragment(size_t                             idx,
                                                detinfo::DetectorClocksData const& clockData,
                                                art::Handle<artdaq::Fragments>     fragmentHandle,
                                                ConcurrentRawDigitCol&             rawDigitCollection,
                                                ConcurrentRawDigitCol&             rawRawDigitCollection,
                                                ConcurrentRawDigitCol&             coherentCollection) const
{
    cet::cpu_timer theClockProcess;

    theClockProcess.start();

    art::Ptr<artdaq::Fragment> fragmentPtr(fragmentHandle, idx);

    mf::LogDebug("DaqDecoderICARUSTPC") << "--> Processing fragment ID: " << fragmentPtr->fragmentID() << std::endl;
    mf::LogDebug("DaqDecoderICARUSTPC") << "    ==> Current thread index: " << tbb::this_task_arena::current_thread_index() << std::endl;

    // Recover pointer to the decoder needed here
    IDecoderFilter* decoderTool = fDecoderToolVec[tbb::this_task_arena::current_thread_index()].get();

    //process_fragment(event, rawfrag, product_collection, header_collection);
    decoderTool->process_fragment(clockData, *fragmentPtr);

    // Useful numerology
    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment physCrateFragment(*fragmentPtr);

//    size_t nBoardsPerFragment = physCrateFragment.nBoards();
//    size_t nChannelsPerBoard  = physCrateFragment.nChannelsPerBoard();

    // Set base channel for both the board and the board/fragment
//    size_t boardFragOffset    = nChannelsPerBoard * nBoardsPerFragment * (fragmentPtr->fragmentID() + fFragmentOffset);

    theClockProcess.stop();

    double totalTime = theClockProcess.accumulated_real_time();

    // We need to recalculate pedestals for the noise corrected waveforms
    icarus_signal_processing::WaveformTools<float> waveformTools;

    // Save the filtered RawDigitsactive but for corrected raw digits pedestal is zero
    icarus_signal_processing::VectorFloat       locPedsVec(decoderTool->getWaveLessCoherent().size(),0.);
    icarus_signal_processing::VectorFloat       locFullRMSVec(locPedsVec.size(),0.);
    icarus_signal_processing::VectorFloat       locTruncRMSVec(locPedsVec.size(),0.);
    icarus_signal_processing::VectorInt         locNumTruncBins(locPedsVec.size(),0);
    icarus_signal_processing::VectorInt         locRangeBins(locPedsVec.size(),0);

    const icarus_signal_processing::VectorInt&  channelVec   = decoderTool->getChannelIDs();
    const icarus_signal_processing::ArrayFloat& corWaveforms = decoderTool->getWaveLessCoherent();

    icarus_signal_processing::ArrayFloat        pedCorWaveforms(corWaveforms.size(),icarus_signal_processing::VectorFloat(corWaveforms[0].size()));

    for(size_t idx = 0; idx < corWaveforms.size(); idx++)
    {
        // Now determine the pedestal and correct for it
        waveformTools.getPedestalCorrectedWaveform(corWaveforms[idx],
                                                   pedCorWaveforms[idx],
                                                   fSigmaForTruncation,
                                                   locPedsVec[idx],
                                                   locFullRMSVec[idx],
                                                   locTruncRMSVec[idx],
                                                   locNumTruncBins[idx],
                                                   locRangeBins[idx]);
    }

    saveRawDigits(pedCorWaveforms, locPedsVec, locTruncRMSVec, channelVec, rawDigitCollection);

    // Optionally, save the raw RawDigits
    if (fOutputRawWaveform)
        saveRawDigits(decoderTool->getRawWaveforms(),decoderTool->getPedestalVals(),decoderTool->getFullRMSVals(),channelVec, rawRawDigitCollection);

    // Also optional is to output the coherent corrections (note there will be many fewer of these! )
    if (fOutputCorrection)
//        saveRawDigits(decoderTool->getCorrectedMedians(),decoderTool->getPedestalVals(),decoderTool->getFullRMSVals(),channelVec,coherentCollection);
        saveRawDigits(decoderTool->getCorrectedMedians(),locPedsVec,decoderTool->getFullRMSVals(),channelVec,coherentCollection);

    mf::LogDebug("DaqDecoderICARUSTPC") << "--> Exiting fragment processing for thread: " << tbb::this_task_arena::current_thread_index() << ", time: " << totalTime << std::endl;
    return;
}

void DaqDecoderICARUSTPC::saveRawDigits(const icarus_signal_processing::ArrayFloat&  dataArray, 
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

        mf::LogDebug("DaqDecoderICARUSTPC") << "    --> saving rawdigits for " << dataArray.size() << " channels" << std::endl;

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

        mf::LogDebug("DaqDecoderICARUSTPC") << "    --> done with save, time: " << totalTime << std::endl;
    }

    return;
}

//----------------------------------------------------------------------------
/// End job method.
void DaqDecoderICARUSTPC::endJob(art::ProcessingFrame const&)
{
    mf::LogInfo("DaqDecoderICARUSTPC") << "Looked at " << fNumEvent << " events" << std::endl;
}

} // end of namespace
