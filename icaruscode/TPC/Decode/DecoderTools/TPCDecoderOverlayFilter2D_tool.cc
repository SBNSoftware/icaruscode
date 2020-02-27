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
#include "art/Utilities/make_tool.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/TPC/Decode/DecoderTools/IDecoderFilter.h"
#include "icaruscode/TPC/Decode/DecoderTools/IFakeParticle.h"

#include "icarussigproc/WaveformTools.h"
#include "icarussigproc/Denoising.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  TPCDecoderFilter1D class definiton
 */
class TPCDecoderFilter1D : virtual public IDecoderFilter
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit TPCDecoderFilter1D(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TPCDecoderFilter1D();

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
    virtual void process_fragment(const artdaq::Fragment&) override;

    /**
     *  @brief Recover the selection values
     */
    const icarussigproc::ArrayBool  getSelectionVals() const override {return fSelectVals;};

    /**
     *  @brief Recover the ROI values
     */
    const icarussigproc::ArrayBool  getROIVals() const override {return fROIVals;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarussigproc::ArrayFloat getRawWaveforms() const override {return fPedSubtractedWaveforms;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarussigproc::ArrayFloat getPedCorWaveforms() const override {return fPedSubtractedWaveforms;};

    /**
     *  @brief Recover the "intrinsic" RMS
     */
    const icarussigproc::ArrayFloat getIntrinsicRMS() const override {return fIntrinsicRMS;};

    /**
     *  @brief Recover the correction median values
     */
    const icarussigproc::ArrayFloat getCorrectedMedians() const override {return fCorrectedMedians;};

    /**
     *  @brief Recover the waveforms less coherent noise
     */
    const icarussigproc::ArrayFloat getWaveLessCoherent()  const override {return fWaveLessCoherent;};

    /**
     *  @brief Recover the morphological filter waveforms
     */
    const icarussigproc::ArrayFloat getMorphedWaveforms()  const override {return fMorphedWaveforms;};

    /**
     *  @brief Recover the pedestals for each channel
     */
    const icarussigproc::VectorFloat getPedestalVals() const override {return fPedestalVals;};

    /**
     *  @brief Recover the full RMS before coherent noise
     */
    const icarussigproc::VectorFloat getFullRMSVals()  const override {return fFullRMSVals;};
 
    /**
     *  @brief Recover the truncated RMS noise 
     */
    const icarussigproc::VectorFloat getTruncRMSVals() const override {return fTruncRMSVals;};

    /**
     *  @brief Recover the number of bins after truncation
     */
    const icarussigproc::VectorInt   getNumTruncBins() const override {return fNumTruncBins;};

private:

    uint32_t                           fFragment_id_offset;     //< Allow offset for id
    size_t                             fCoherentNoiseGrouping;  //< # channels in common for coherent noise
    std::vector<size_t>                fStructuringElement;     //< Structuring element for morphological filter
    size_t                             fMorphWindow;            //< Window for filter
    float                              fThreshold;              //< Threshold to apply for saving signal

    std::vector<char>                  fFilterModeVec;          //< Allowed modes for the filter

    // Allocate containers for noise processing
    icarussigproc::ArrayBool           fSelectVals;
    icarussigproc::ArrayBool           fROIVals;
    icarussigproc::ArrayFloat          fPedSubtractedWaveforms;
    icarussigproc::ArrayFloat          fIntrinsicRMS;
    icarussigproc::ArrayFloat          fCorrectedMedians;
    icarussigproc::ArrayFloat          fWaveLessCoherent;
    icarussigproc::ArrayFloat          fMorphedWaveforms;
         
    icarussigproc::VectorFloat         fPedestalVals;
    icarussigproc::VectorFloat         fFullRMSVals;
    icarussigproc::VectorFloat         fTruncRMSVals;
    icarussigproc::VectorInt           fNumTruncBins;

    // Overlay tool
    std::unique_ptr<IFakeParticle>     fFakeParticleTool;

    const geo::Geometry*               fGeometry;              //< pointer to the Geometry service
    const detinfo::DetectorProperties* fDetector;              //< Pointer to the detector properties
};

TPCDecoderFilter1D::TPCDecoderFilter1D(fhicl::ParameterSet const &pset)
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

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TPCDecoderFilter1D::~TPCDecoderFilter1D()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void TPCDecoderFilter1D::configure(fhicl::ParameterSet const &pset)
{
    fFragment_id_offset    = pset.get<uint32_t           >("fragment_id_offset"                                );
    fCoherentNoiseGrouping = pset.get<size_t             >("CoherentGrouping",                               64);
    fStructuringElement    = pset.get<std::vector<size_t>>("StructuringElement", std::vector<size_t>() = {7,20});
    fMorphWindow           = pset.get<size_t             >("FilterWindow",                                   10);
    fThreshold             = pset.get<float              >("Threshold",                                     7.5);

    fFilterModeVec         = {'d','e','g'};

    fGeometry              = art::ServiceHandle<geo::Geometry const>{}.get();
    fDetector              = lar::providerFrom<detinfo::DetectorPropertiesService>();

    fFakeParticleTool      = art::make_tool<IFakeParticle>(pset.get<fhicl::ParameterSet>("FakeParticle"));

    return;
}

void TPCDecoderFilter1D::process_fragment(const artdaq::Fragment &fragment)
{
    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment physCrateFragment(fragment);
    
    size_t nBoardsPerFragment   = physCrateFragment.nBoards();
    size_t nChannelsPerBoard    = physCrateFragment.nChannelsPerBoard();
    size_t nChannelsPerFragment = nBoardsPerFragment * nChannelsPerBoard;
    size_t nSamplesPerChannel   = physCrateFragment.nSamplesPerChannel();

    if (fSelectVals.empty())             fSelectVals             = icarussigproc::ArrayBool(nChannelsPerFragment,  icarussigproc::VectorBool(nSamplesPerChannel));
    if (fROIVals.empty())                fROIVals                = icarussigproc::ArrayBool(nChannelsPerFragment,  icarussigproc::VectorBool(nSamplesPerChannel));
    if (fPedSubtractedWaveforms.empty()) fPedSubtractedWaveforms = icarussigproc::ArrayFloat(nChannelsPerFragment, icarussigproc::VectorFloat(nSamplesPerChannel));
    if (fIntrinsicRMS.empty())           fIntrinsicRMS           = icarussigproc::ArrayFloat(nChannelsPerFragment, icarussigproc::VectorFloat(nSamplesPerChannel));
    if (fCorrectedMedians.empty())       fCorrectedMedians       = icarussigproc::ArrayFloat(nChannelsPerFragment, icarussigproc::VectorFloat(nSamplesPerChannel));
    if (fWaveLessCoherent.empty())       fWaveLessCoherent       = icarussigproc::ArrayFloat(nChannelsPerFragment, icarussigproc::VectorFloat(nSamplesPerChannel));
    if (fMorphedWaveforms.empty())       fMorphedWaveforms       = icarussigproc::ArrayFloat(nChannelsPerFragment, icarussigproc::VectorFloat(nSamplesPerChannel));

    if (fPedestalVals.empty())           fPedestalVals           = icarussigproc::VectorFloat(nChannelsPerFragment);
    if (fFullRMSVals.empty())            fFullRMSVals            = icarussigproc::VectorFloat(nChannelsPerFragment);
    if (fTruncRMSVals.empty())           fTruncRMSVals           = icarussigproc::VectorFloat(nChannelsPerFragment);
    if (fNumTruncBins.empty())           fNumTruncBins           = icarussigproc::VectorInt(nChannelsPerFragment);

    // Allocate the de-noising object
    icarussigproc::Denoising            denoiser;
    icarussigproc::WaveformTools<float> waveformTools;

    // The first task is to recover the data from the board data block, determine and subtract the pedestals
    // and store into vectors useful for the next steps
    for(size_t board = 0; board < nBoardsPerFragment; board++)
    {
        // Keep these for a while longer as we may want to do some checking soon
//        size_t event_number = physCrateFragment.BoardEventNumber(i_b);
//        size_t timestamp    = physCrateFragment.BoardTimeStamp(board);

        // This is where we would recover the base channel for the board from database/module
        size_t boardOffset = nChannelsPerBoard * board;

        // Get the pointer to the start of this board's block of data
        const icarus::A2795DataBlock::data_t* dataBlock = physCrateFragment.BoardData(board);

        // Copy to input data array
        for(size_t chanIdx = 0; chanIdx < nChannelsPerBoard; chanIdx++)
        {
            // Get the channel number on the Fragment
            size_t channelOnBoard = boardOffset + chanIdx;

            icarussigproc::VectorFloat& dataVec = fPedSubtractedWaveforms[channelOnBoard];

            for(size_t tick = 0; tick < nSamplesPerChannel; tick++)
                dataVec[tick] = dataBlock[chanIdx + tick * nChannelsPerBoard];

            // Now determine the pedestal and correct for it
            waveformTools.getPedestalCorrectedWaveform(dataVec, 
                                                       dataVec,
                                                       3,
                                                       fPedestalVals[channelOnBoard], 
                                                       fFullRMSVals[channelOnBoard], 
                                                       fTruncRMSVals[channelOnBoard], 
                                                       fNumTruncBins[channelOnBoard]);
        }
    }

    // Overlay a fake particle on this array of waveforms
    fFakeParticleTool->overlayFakeParticle(fPedSubtractedWaveforms);

    // Run the coherent filter
    denoiser.removeCoherentNoise2D(fWaveLessCoherent,fPedSubtractedWaveforms,fMorphedWaveforms,fIntrinsicRMS,fSelectVals,fROIVals,fCorrectedMedians,
                                   fFilterModeVec[2],fCoherentNoiseGrouping,fStructuringElement[0],fStructuringElement[1],fMorphWindow,fThreshold);

    return;
}

DEFINE_ART_CLASS_TOOL(TPCDecoderFilter1D)
} // namespace lar_cluster3d
