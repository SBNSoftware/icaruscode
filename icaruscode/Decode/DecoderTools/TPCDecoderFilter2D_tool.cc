/**
 *  @file   TPCDecoderFilter2D_tool.cc
 *
 *  @brief  This tool converts from daq to LArSoft format with noise filtering
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoderFilter.h"

#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Denoising.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  TPCDecoderFilter2D class definiton
 */
class TPCDecoderFilter2D : virtual public IDecoderFilter
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit TPCDecoderFilter2D(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TPCDecoderFilter2D();

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) override;

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param fragment            The artdaq fragment to process
     */
    virtual void process_fragment(detinfo::DetectorClocksData const& clockData,
                                  const artdaq::Fragment&) override;

    /**
     *  @brief Recover the channels for the processed fragment
     */
    const icarus_signal_processing::VectorInt  getChannelIDs()       const override {return fChannelIDVec;}

    /**
     *  @brief Recover the selection values
     */
    const icarus_signal_processing::ArrayBool  getSelectionVals() const override {return fSelectVals;};

    /**
     *  @brief Recover the ROI values
     */
    const icarus_signal_processing::ArrayBool  getROIVals() const override {return fROIVals;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat getRawWaveforms() const override {return fPedSubtractedWaveforms;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat getPedCorWaveforms() const override {return fPedSubtractedWaveforms;};

    /**
     *  @brief Recover the "intrinsic" RMS
     */
    const icarus_signal_processing::ArrayFloat getIntrinsicRMS() const override {return fIntrinsicRMS;};

    /**
     *  @brief Recover the correction median values
     */
    const icarus_signal_processing::ArrayFloat getCorrectedMedians() const override {return fCorrectedMedians;};

    /**
     *  @brief Recover the waveforms less coherent noise
     */
    const icarus_signal_processing::ArrayFloat getWaveLessCoherent()  const override {return fWaveLessCoherent;};

    /**
     *  @brief Recover the morphological filter waveforms
     */
    const icarus_signal_processing::ArrayFloat getMorphedWaveforms()  const override {return fMorphedWaveforms;};

    /**
     *  @brief Recover the pedestals for each channel
     */
    const icarus_signal_processing::VectorFloat getPedestalVals() const override {return fPedestalVals;};

    /**
     *  @brief Recover the full RMS before coherent noise
     */
    const icarus_signal_processing::VectorFloat getFullRMSVals()  const override {return fFullRMSVals;};
 
    /**
     *  @brief Recover the truncated RMS noise 
     */
    const icarus_signal_processing::VectorFloat getTruncRMSVals() const override {return fTruncRMSVals;};

    /**
     *  @brief Recover the number of bins after truncation
     */
    const icarus_signal_processing::VectorInt   getNumTruncBins() const override {return fNumTruncBins;};

private:

    uint32_t                              fFragment_id_offset;     //< Allow offset for id
    size_t                                fCoherentNoiseGrouping;  //< # channels in common for coherent noise
    std::vector<size_t>                   fStructuringElement;     //< Structuring element for morphological filter
    size_t                                fMorphWindow;            //< Window for filter
    float                                 fThreshold;              //< Threshold to apply for saving signal

    std::vector<char>                     fFilterModeVec;          //< Allowed modes for the filter

    // Allocate containers for noise processing
    icarus_signal_processing::VectorInt   fChannelIDVec;
    icarus_signal_processing::ArrayBool   fSelectVals;
    icarus_signal_processing::ArrayBool   fROIVals;
    icarus_signal_processing::ArrayFloat  fPedSubtractedWaveforms;
    icarus_signal_processing::ArrayFloat  fIntrinsicRMS;
    icarus_signal_processing::ArrayFloat  fCorrectedMedians;
    icarus_signal_processing::ArrayFloat  fWaveLessCoherent;
    icarus_signal_processing::ArrayFloat  fMorphedWaveforms;

    icarus_signal_processing::VectorFloat fPedestalVals;
    icarus_signal_processing::VectorFloat fFullRMSVals;
    icarus_signal_processing::VectorFloat fTruncRMSVals;
    icarus_signal_processing::VectorInt   fNumTruncBins;
    icarus_signal_processing::VectorInt   fRangeBins;

    icarus_signal_processing::VectorFloat fThresholdVec;

    const geo::Geometry*                  fGeometry;              //< pointer to the Geometry service
};

TPCDecoderFilter2D::TPCDecoderFilter2D(fhicl::ParameterSet const &pset)
{
    this->configure(pset);

    fSelectVals.clear();
    fROIVals.clear();
    fPedSubtractedWaveforms.clear();
    fIntrinsicRMS.clear();
    fCorrectedMedians.clear();
    fWaveLessCoherent.clear();
    fMorphedWaveforms.clear();

    fPedestalVals.clear();
    fFullRMSVals.clear();
    fTruncRMSVals.clear();
    fNumTruncBins.clear();
    fRangeBins.clear();

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TPCDecoderFilter2D::~TPCDecoderFilter2D()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void TPCDecoderFilter2D::configure(fhicl::ParameterSet const &pset)
{
    fFragment_id_offset    = pset.get<uint32_t           >("fragment_id_offset");
    fCoherentNoiseGrouping = pset.get<size_t             >("CoherentGrouping",   64);
    fStructuringElement    = pset.get<std::vector<size_t>>("StructuringElement", std::vector<size_t>()={7,20});
    fMorphWindow           = pset.get<size_t             >("FilterWindow",       10);
    fThreshold             = pset.get<float              >("Threshold",         7.5);

    fFilterModeVec         = {'d','e','g'};

    fGeometry              = art::ServiceHandle<geo::Geometry const>{}.get();

    return;
}

void TPCDecoderFilter2D::process_fragment(detinfo::DetectorClocksData const&,
                                          const artdaq::Fragment &fragment)
{
    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment physCrateFragment(fragment);
    
    size_t nBoardsPerFragment   = physCrateFragment.nBoards();
    size_t nChannelsPerBoard    = physCrateFragment.nChannelsPerBoard();
    size_t nChannelsPerFragment = nBoardsPerFragment * nChannelsPerBoard;
    size_t nSamplesPerChannel   = physCrateFragment.nSamplesPerChannel();

    if (fSelectVals.empty())             fSelectVals             = icarus_signal_processing::ArrayBool(nChannelsPerFragment,  icarus_signal_processing::VectorBool(nSamplesPerChannel));
    if (fROIVals.empty())                fROIVals                = icarus_signal_processing::ArrayBool(nChannelsPerFragment,  icarus_signal_processing::VectorBool(nSamplesPerChannel));
    if (fPedSubtractedWaveforms.empty()) fPedSubtractedWaveforms = icarus_signal_processing::ArrayFloat(nChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fIntrinsicRMS.empty())           fIntrinsicRMS           = icarus_signal_processing::ArrayFloat(nChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fCorrectedMedians.empty())       fCorrectedMedians       = icarus_signal_processing::ArrayFloat(nChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fWaveLessCoherent.empty())       fWaveLessCoherent       = icarus_signal_processing::ArrayFloat(nChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fMorphedWaveforms.empty())       fMorphedWaveforms       = icarus_signal_processing::ArrayFloat(nChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));

    if (fPedestalVals.empty())           fPedestalVals           = icarus_signal_processing::VectorFloat(nChannelsPerFragment);
    if (fFullRMSVals.empty())            fFullRMSVals            = icarus_signal_processing::VectorFloat(nChannelsPerFragment);
    if (fTruncRMSVals.empty())           fTruncRMSVals           = icarus_signal_processing::VectorFloat(nChannelsPerFragment);
    if (fNumTruncBins.empty())           fNumTruncBins           = icarus_signal_processing::VectorInt(nChannelsPerFragment);
    if (fRangeBins.empty())              fRangeBins              = icarus_signal_processing::VectorInt(nChannelsPerFragment);

    if (fThresholdVec.empty())           fThresholdVec           = icarus_signal_processing::VectorFloat(nChannelsPerFragment);

    // Allocate the de-noising object
    icarus_signal_processing::Denoising            denoiser;
    icarus_signal_processing::WaveformTools<float> waveformTools;

    // The first task is to recover the data from the board data block, determine and subtract the pedestals
    // and store into vectors useful for the next steps
    for(size_t board = 0; board < nBoardsPerFragment; board++)
    {
//        size_t event_number = physCrateFragment.BoardEventNumber(i_b);
//        size_t timestamp    = physCrateFragment.BoardTimeStamp(board);

        size_t boardOffset = nChannelsPerBoard * board;

        // Get the pointer to the start of this board's block of data
        const icarus::A2795DataBlock::data_t* dataBlock = physCrateFragment.BoardData(board);

        // Copy to input data array
        for(size_t chanIdx = 0; chanIdx < nChannelsPerBoard; chanIdx++)
        {
            // Get the channel number on the Fragment
            size_t channelOnBoard = boardOffset + chanIdx;

            icarus_signal_processing::VectorFloat& dataVec = fPedSubtractedWaveforms[channelOnBoard];

            fThresholdVec[channelOnBoard] = fThreshold;

            for(size_t tick = 0; tick < nSamplesPerChannel; tick++)
                dataVec[tick] = dataBlock[chanIdx + tick * nChannelsPerBoard];

            // Now determine the pedestal and correct for it
            waveformTools.getPedestalCorrectedWaveform(dataVec, 
                                                       dataVec,
                                                       3,
                                                       fPedestalVals[channelOnBoard], 
                                                       fFullRMSVals[channelOnBoard], 
                                                       fTruncRMSVals[channelOnBoard], 
                                                       fNumTruncBins[channelOnBoard],
                                                       fRangeBins[channelOnBoard]);
        }
    }

    // Run the coherent filter
    denoiser.removeCoherentNoise2D(fWaveLessCoherent,fPedSubtractedWaveforms,fMorphedWaveforms,fIntrinsicRMS,fSelectVals,fROIVals,fCorrectedMedians,
                                   fThresholdVec, fFilterModeVec[2],fCoherentNoiseGrouping,fStructuringElement[0],fStructuringElement[1],
                                   fMorphWindow);

    return;
}

DEFINE_ART_CLASS_TOOL(TPCDecoderFilter2D)
} // namespace lar_cluster3d
