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
    using RawDigitCollection    = std::vector<raw::RawDigit>;
    using RawDigitCollectionPtr = std::unique_ptr<RawDigitCollection>;
    using ConcurrentRawDigitCol = tbb::concurrent_vector<raw::RawDigit>;

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

    using ChannelArrayPair      = std::pair<ChannelVec,icarus_signal_processing::ArrayFloat>;
    using ChannelArrayPairVec   = std::vector<ChannelArrayPair>;


    // Function to do the work
    void processSingleFragment(size_t,
                               detinfo::DetectorClocksData const& clockData,
                               art::Handle<artdaq::Fragments>, 
                               ChannelArrayPairVec&) const;

    // Function to do the work
    void processSingleImage(size_t,
                            const detinfo::DetectorClocksData&,
                            const ChannelArrayPairVec&,
                            ConcurrentRawDigitCol&,
                            ConcurrentRawDigitCol&,
                            ConcurrentRawDigitCol&) const;

private:
    class multiThreadFragmentProcessing
    {
    public:
        multiThreadFragmentProcessing(DaqDecoderICARUSTPCwROI const&     parent,
                                      detinfo::DetectorClocksData const& clockData,
                                      art::Handle<artdaq::Fragments>&    fragmentsHandle,
                                      ChannelArrayPairVec&               channelArrayPairVec)
            : fDaqDecoderICARUSTPCwROI(parent),
              fClockData{clockData},
              fFragmentsHandle(fragmentsHandle),
              fChannelArrayPairVec(channelArrayPairVec)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
              fDaqDecoderICARUSTPCwROI.processSingleFragment(idx, fClockData, fFragmentsHandle, fChannelArrayPairVec);
        }
    private:
        const DaqDecoderICARUSTPCwROI&     fDaqDecoderICARUSTPCwROI;
        detinfo::DetectorClocksData const& fClockData;
        art::Handle<artdaq::Fragments>&    fFragmentsHandle;
        ChannelArrayPairVec&               fChannelArrayPairVec;
    };

    class multiThreadImageProcessing
    {
    public:
        multiThreadImageProcessing(DaqDecoderICARUSTPCwROI     const& parent,
                                   detinfo::DetectorClocksData const& clockData,
                                   ChannelArrayPairVec         const& channelArrayPairVec,
                                   ConcurrentRawDigitCol&             concurrentRawDigits,
                                   ConcurrentRawDigitCol&             concurrentRawRawDigits,
                                   ConcurrentRawDigitCol&             coherentRawDigits)
            : fDaqDecoderICARUSTPCwROI(parent),
              fClockData{clockData},
              fChannelArrayPairVec(channelArrayPairVec),
              fConcurrentRawDigits(concurrentRawDigits),
              fConcurrentRawRawDigits(concurrentRawRawDigits),
              fCoherentRawDigits(coherentRawDigits)
        {}

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
              fDaqDecoderICARUSTPCwROI.processSingleImage(idx, fClockData, fChannelArrayPairVec, fConcurrentRawDigits, fConcurrentRawRawDigits, fCoherentRawDigits);
        }
    private:
        const DaqDecoderICARUSTPCwROI&     fDaqDecoderICARUSTPCwROI;
        const detinfo::DetectorClocksData& fClockData;
        const ChannelArrayPairVec&         fChannelArrayPairVec;
        ConcurrentRawDigitCol&             fConcurrentRawDigits;
        ConcurrentRawDigitCol&             fConcurrentRawRawDigits;
        ConcurrentRawDigitCol&             fCoherentRawDigits;
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
    std::string                                                 fOutputRawWavePath;          ///< Path to assign to the output if asked for
    std::string                                                 fOutputCoherentPath;         ///< Path to assign to the output if asked for
    bool                                                        fDiagnosticOutput;           ///< Set this to get lots of messages

    const std::string                                           fLogCategory;                ///< Output category when logging messages

    // fhicl parameters
    std::vector<size_t>                                         fStructuringElement;         ///< Structuring element for morphological filter
    std::vector<float>                                          fThreshold;                  ///< Threshold to apply for saving signal

    // Start with parameters for Butterworth Filter
    unsigned int                                                fButterworthOrder;           ///< Order parameter for Butterworth filter
    unsigned int                                                fButterworthThreshold;       ///< Threshold for Butterworth filter

    // Parameters for the 2D morphological filter
    unsigned int                                                fMorph2DStructuringElementX; ///< Structuring element in X
    unsigned int                                                fMorph2DStructuringElementY; ///< Structuring element in Y

    // Parameters for the denoiser
    unsigned int                                                fCoherentNoiseGrouping;      ///< Number of consecutive channels in coherent noise subtraction
    unsigned int                                                fCoherentNoiseOffset;        ///< Offset for the midplane...
    unsigned int                                                fMorphologicalWindow;        ///< Window size for filter
//    bool                                                        fOutputStats;                ///< Output of timiing statistics?
    float                                                       fCoherentThresholdFactor;    ///< Threshold factor for coherent noise removal

    // Parameters for the ROI finding
    unsigned int                                                fADFilter_SX;                ///< 
    unsigned int                                                fADFilter_SY;                ///< 
    float                                                       fSigma_x;                    ///<
    float                                                       fSigma_y;                    ///<
    float                                                       fSigma_r;                    ///<
    float                                                       fLowThreshold;               ///<
    float                                                       fHighThreshold;              ///<
    unsigned int                                                fBinaryDilation_SX;          ///<
    unsigned int                                                fBinaryDilation_SY;          ///<

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

    // Our functions
    std::unique_ptr<icarus_signal_processing::IFFTFilterFunction>        fButterworthFilter;
    std::unique_ptr<icarus_signal_processing::IMorphologicalFunctions2D> fMorphologicalFilter;
    std::unique_ptr<icarus_signal_processing::IDenoiser2D>               fDenoiser2D;
    std::unique_ptr<icarus_signal_processing::BilateralFilters>          fBilateralFilters;
    std::unique_ptr<icarus_signal_processing::EdgeDetection>             fEdgeDetection;
    std::unique_ptr<icarus_signal_processing::IROIFinder2D>              fROIFinder2D;

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
                          art::ReplicatedProducer(pset, frame),
                          fLogCategory("DaqDecoderICARUSTPCwROI"),fNumEvent(0), fNumROPs(0)
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

    // Set up a WireID to ROP plane number table
    for(size_t cryoIdx = 0; cryoIdx < 2; cryoIdx++)
    {
        for(size_t logicalTPCIdx = 0; logicalTPCIdx < 4; logicalTPCIdx++)
        {
            for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
            {
                geo::PlaneID planeID(cryoIdx,logicalTPCIdx,planeIdx);

                raw::ChannelID_t channel = fGeometry->PlaneWireToChannel(planeID.Plane, 0, planeID.TPC, planeID.Cryostat);

                readout::ROPID ropID = fGeometry->ChannelToROP(channel);

                fPlaneToROPPlaneMap[planeID]   = ropID.ROP;
                fPlaneToWireOffsetMap[planeID] = channel;
                fROPToNumWiresMap[ropID.ROP]   = fGeometry->Nwires(planeID);

                // Special case handling
                if (ropID.ROP > 1) fROPToNumWiresMap[ropID.ROP] *= 2;

                if (ropID.ROP > fNumROPs) fNumROPs = ropID.ROP;

                // Watch for the middle induction and collection plane logical TPC split
                if (ropID.ROP > 1 && (logicalTPCIdx == 1 || logicalTPCIdx == 3))
                {
                    geo::PlaneID tempID(cryoIdx,logicalTPCIdx-1,planeIdx);

                    fPlaneToWireOffsetMap[planeID] = fPlaneToWireOffsetMap[tempID];
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
    fFragmentsLabelVec          = pset.get<std::vector<art::InputTag>>("FragmentsLabelVec",  std::vector<art::InputTag>()={"daq:PHYSCRATEDATA"});
    fOutputRawWaveform          = pset.get<bool                      >("OutputRawWaveform",                                               false);
    fOutputCorrection           = pset.get<bool                      >("OutputCorrection",                                                false);
    fOutputRawWavePath          = pset.get<std::string               >("OutputRawWavePath",                                               "raw");
    fOutputCoherentPath         = pset.get<std::string               >("OutputCoherentPath",                                              "Cor");
    fDiagnosticOutput           = pset.get<bool                      >("DiagnosticOutput",                                                false);

    // Recover parameters for noise/ROI
    fStructuringElement         = pset.get<std::vector<size_t>       >("StructuringElement",                       std::vector<size_t>()={8,16});
    fThreshold                  = pset.get<std::vector<float>        >("Threshold",                       std::vector<float>()={2.75,2.75,2.75});

    fButterworthOrder           = pset.get<unsigned int              >("ButterworthOrder",     2);
    fButterworthThreshold       = pset.get<unsigned int              >("ButterworthThreshld", 30);

    //fButterworthFilter = std::make_unique<icarus_signal_processing::HighPassButterworthFilter>(fButterworthThreshold,fButterworthOrder,4096);
    fButterworthFilter = std::make_unique<icarus_signal_processing::NoFFTFilter>();

    fMorph2DStructuringElementX = pset.get<unsigned int              >("Morph2DStructuringElementX", 7);
    fMorph2DStructuringElementY = pset.get<unsigned int              >("Morph2DStructuringElementX", 28);

    fMorphologicalFilter = std::make_unique<icarus_signal_processing::Dilation2D>(fMorph2DStructuringElementX,fMorph2DStructuringElementY);

    fCoherentNoiseGrouping      = pset.get<unsigned int              >("CoherentNoiseGrouping",    32);
    fCoherentNoiseOffset        = pset.get<unsigned int              >("CoherentNoiseOffset",      24);
    fMorphologicalWindow        = pset.get<unsigned int              >("MorphologicalWindow",      10);
    fCoherentThresholdFactor    = pset.get<float                     >("CoherentThresholdFactor", 2.5);

    fThresholdVec.resize(6560/fCoherentNoiseGrouping,fCoherentThresholdFactor);

    //fDenoiser2D = std::make_unique<icarus_signal_processing::Denoiser2D_Hough>(fMorphologicalFilter.get(), fThresholdVec, fCoherentNoiseGrouping, fCoherentNoiseOffset, fMorphologicalWindow);
    fDenoiser2D = std::make_unique<icarus_signal_processing::Denoiser2D>(fMorphologicalFilter.get(), fThresholdVec, fCoherentNoiseGrouping, fMorphologicalWindow);

    fADFilter_SX                = pset.get<unsigned int              >("ADFilter_SX",         7);
    fADFilter_SY                = pset.get<unsigned int              >("ADFilter_SY",         7);
    fSigma_x                    = pset.get<float                     >("Sigma_x",          10.0);
    fSigma_y                    = pset.get<float                     >("Sigma_y",          10.0);
    fSigma_r                    = pset.get<float                     >("Sigma_r",          30.0);
                   
    fLowThreshold               = pset.get<float                     >("LowThreshold",     10.0);
    fHighThreshold              = pset.get<float                     >("HighThreshold",    20.0); 
                   
    fBinaryDilation_SX          = pset.get<unsigned int              >("BinaryDilation_SX",  31);
    fBinaryDilation_SY          = pset.get<unsigned int              >("BinaryDilation_SY",  31);

    fBilateralFilters = std::make_unique<icarus_signal_processing::BilateralFilters>();
    fEdgeDetection    = std::make_unique<icarus_signal_processing::EdgeDetection>();

    fROIFinder2D = std::make_unique<icarus_signal_processing::ROICannyFilter>(fButterworthFilter.get(), 
                                                                              fDenoiser2D.get(), 
                                                                              fBilateralFilters.get(),
                                                                              fEdgeDetection.get(), 
                                                                              fADFilter_SX,
                                                                              fADFilter_SY,
                                                                              fSigma_x,
                                                                              fSigma_y,
                                                                              fSigma_r,
                                                                              fLowThreshold,
                                                                              fHighThreshold,
                                                                              fBinaryDilation_SX,
                                                                              fBinaryDilation_SY);

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

        multiThreadFragmentProcessing fragmentProcessing(*this, clockData, daq_handle, channelArrayPairVec);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, daq_handle->size()), fragmentProcessing);

        // Now let's process the resulting images
        multiThreadImageProcessing imageProcessing(*this, clockData, channelArrayPairVec, concurrentRawDigits, concurrentRawRawDigits, coherentRawDigits);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, fNumROPs), imageProcessing);
    
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

    mf::LogInfo(fLogCategory) << "==> DaqDecoderICARUSTPCwROI total time: " << totalTime << std::endl;

    return;
}

void DaqDecoderICARUSTPCwROI::processSingleFragment(size_t                             idx,
                                                    detinfo::DetectorClocksData const& clockData,
                                                    art::Handle<artdaq::Fragments>     fragmentHandle,
                                                    ChannelArrayPairVec&               channelArrayPairVec) const
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

        std::stringstream outputString;
        outputString << "********************************************************************************\n"
                     << "FragmentID: " << std::hex << fragmentID << ", Crate: " << crateName << std::dec << ", boardID: " << boardSlot << "/" << nBoardsPerFragment << ", size " << channelPlanePairVec.size() << "/" << nChannelsPerBoard;

        mf::LogInfo("TPCDecoderFilter1D") << outputString.str();

        // Get the pointer to the start of this board's block of data
        const icarus::A2795DataBlock::data_t* dataBlock = physCrateFragment.BoardData(board);

        std::stringstream channelString;

        // Copy to input data array
        for(size_t chanIdx = 0; chanIdx < nChannelsPerBoard; chanIdx++)
        {
            // Get the channel number on the Fragment
            raw::ChannelID_t channel = channelPlanePairVec[chanIdx].first;

            channelString << channel << "-";
            
            std::vector<geo::WireID> wireIDVec = fGeometry->ChannelToWire(channel);

            // Skip the channels which are not physically connected
            if (wireIDVec.empty())
            {
                channelString << "skip * ";
                mf::LogDebug(fLogCategory) << channelString.str();

                continue;
            }

            // Some contortions here... the horizontal wires are split and so independent,
            // the angled wires are logically split but we do want them put back together here... 
            const geo::PlaneID& planeID    = wireIDVec[0].planeID();
            unsigned int        planeIndex = fPlaneToROPPlaneMap.find(planeID)->second;
            unsigned int        wire       = channel - fPlaneToWireOffsetMap.find(planeID)->second;

            icarus_signal_processing::VectorFloat& rawDataVec = channelArrayPairVec[planeIndex].second[wire];

            for(size_t tick = 0; tick < nSamplesPerChannel; tick++)
                rawDataVec[tick] = -dataBlock[chanIdx + tick * nChannelsPerBoard];

            // Keep track of the channel
            channelArrayPairVec[planeIndex].first[wire] = channel;

            if (fDiagnosticOutput)
                channelString << wireIDVec[0].Cryostat << "/" << wireIDVec[0].TPC << "/" << wireIDVec[0].Plane << "/" << wireIDVec[0].Wire  
                              << " idx/wire: " << planeIndex << "/" << wire << " * ";
        }

        mf::LogInfo(fLogCategory) << channelString.str();
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

void DaqDecoderICARUSTPCwROI::processSingleImage(size_t                             idx,
                                                 const detinfo::DetectorClocksData& clockData,
                                                 const ChannelArrayPairVec&         channelArrayPairVec,
                                                 ConcurrentRawDigitCol&             concurrentRawDigitCol,
                                                 ConcurrentRawDigitCol&             concurrentRawRawDigitCol,
                                                 ConcurrentRawDigitCol&             coherentRawDigitCol) const
{
    // Tools. We love tools
    icarus_signal_processing::WaveformTools<float> waveformTools;

    // Which image are we processing?
    const ChannelArrayPair& channelArrayPair = channelArrayPairVec[idx];

    // Let's go through and fill the output vector
    const ChannelVec&                           channelVec = channelArrayPair.first;
    const icarus_signal_processing::ArrayFloat& dataArray  = channelArrayPair.second;

    unsigned int numChannels = dataArray.size();
    unsigned int numTicks    = dataArray[0].size();

    icarus_signal_processing::ArrayFloat waveLessCoherent(numChannels,icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat medianVals(numChannels,icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat coherentRMS(numChannels,icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat morphedWaveforms(numChannels,icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat finalErosion(numChannels,icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat fullEvent(numChannels,icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayBool  outputROIs(numChannels,icarus_signal_processing::VectorBool(numTicks,false));

    (*fROIFinder2D)(dataArray,fullEvent,outputROIs,waveLessCoherent,medianVals,coherentRMS,morphedWaveforms,finalErosion);

    // Now set up for output
    raw::RawDigit::ADCvector_t wvfm(dataArray[0].size());

    // Placeholders
    icarus_signal_processing::VectorFloat pedCorDataVec(dataArray[0].size());
    float pedestal;
    float fullRMS;
    float truncRMS;
    int   numTruncBins;
    int   rangeBins;
    float sigmaForTruncation(3.5);

    // Loop over the channels to recover the RawDigits after filtering
    for(size_t chanIdx = 0; chanIdx != dataArray.size(); chanIdx++)
    {
        if (fOutputRawWaveform)
        {
            const icarus_signal_processing::VectorFloat& dataVec = dataArray[chanIdx];

            // Get the pedestal corrections
            waveformTools.getPedestalCorrectedWaveform(dataVec, pedCorDataVec, sigmaForTruncation, pedestal, fullRMS, truncRMS, numTruncBins, rangeBins);

            // Need to convert from float to short int
            std::transform(pedCorDataVec.begin(),pedCorDataVec.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});
            //std::copy(dataVec.begin(),dataVec.end(),wvfm.begin());

            ConcurrentRawDigitCol::iterator newRawObjItr = concurrentRawRawDigitCol.emplace_back(channelVec[chanIdx],wvfm.size(),wvfm); 

            newRawObjItr->SetPedestal(0.,truncRMS);
        }

        if (fOutputCorrection)
        {
            const icarus_signal_processing::VectorFloat& dataVec = medianVals[chanIdx];

            // Need to convert from float to short int
            std::transform(dataVec.begin(),dataVec.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});
            //std::copy(dataVec.begin(),dataVec.end(),wvfm.begin());

            ConcurrentRawDigitCol::iterator newRawObjItr = coherentRawDigitCol.emplace_back(channelVec[chanIdx],wvfm.size(),wvfm); 

            newRawObjItr->SetPedestal(0.,0.);
        }

         // Now the coherent subtracted 
        const icarus_signal_processing::VectorFloat& coherentVec = waveLessCoherent[chanIdx];

         // Need to convert from float to short int
        std::transform(coherentVec.begin(),coherentVec.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});
        //std::copy(dataVec.begin(),dataVec.end(),wvfm.begin());

         ConcurrentRawDigitCol::iterator newObjItr = concurrentRawDigitCol.emplace_back(channelVec[chanIdx],wvfm.size(),wvfm); 

         newObjItr->SetPedestal(0.,0.);

    }//loop over channel indices

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
void DaqDecoderICARUSTPCwROI::endJob(art::ProcessingFrame const&)
{
    mf::LogInfo(fLogCategory) << "Looked at " << fNumEvent << " events" << std::endl;
}

} // end of namespace
