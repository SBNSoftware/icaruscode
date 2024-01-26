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

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoderFilter.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Denoising.h"
#include "icarus_signal_processing/Filters/FFTFilterFunctions.h"

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
    virtual void process_fragment(detinfo::DetectorClocksData const&,
                                  const artdaq::Fragment&) override;

    /**
     *  @brief Recover the channels for the processed fragment
     */
    const icarus_signal_processing::VectorInt  getChannelIDs()       const override {return fChannelIDVec;}

    /**
     *  @brief Recover the selection values
     */
    const icarus_signal_processing::ArrayBool  getSelectionVals()    const override {return fSelectVals;};

    /**
     *  @brief Recover the ROI values
     */
    const icarus_signal_processing::ArrayBool  getROIVals()          const override {return fROIVals;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat getRawWaveforms()     const override {return fRawWaveforms;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat getPedCorWaveforms()  const override {return fPedCorWaveforms;};

    /**
     *  @brief Recover the "intrinsic" RMS
     */
    const icarus_signal_processing::ArrayFloat getIntrinsicRMS()     const override {return fIntrinsicRMS;};

    /**
     *  @brief Recover the correction median values
     */
    const icarus_signal_processing::ArrayFloat getCorrectedMedians()  const override {return fCorrectedMedians;};

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
    const icarus_signal_processing::VectorFloat getPedestalVals()     const override {return fPedestalVals;};

    /**
     *  @brief Recover the full RMS before coherent noise
     */
    const icarus_signal_processing::VectorFloat getFullRMSVals()      const override {return fFullRMSVals;};

    /**
     *  @brief Recover the truncated RMS noise
     */
    const icarus_signal_processing::VectorFloat getTruncRMSVals()     const override {return fTruncRMSVals;};

    /**
     *  @brief Recover the number of bins after truncation
     */
    const icarus_signal_processing::VectorInt   getNumTruncBins()     const override {return fNumTruncBins;};

private:

    uint32_t                                       fFragment_id_offset;     //< Allow offset for id
    float                                          fSigmaForTruncation;     //< Selection cut for truncated rms calculation
    size_t                                         fCoherentNoiseGrouping;  //< # channels in common for coherent noise
    size_t                                         fCoherentNoiseOffset;    //< offset for midplane
    std::vector<size_t>                            fStructuringElement;     //< Structuring element for morphological filter
    size_t                                         fMorphWindow;            //< Window for filter
    std::vector<float>                             fThreshold;              //< Threshold to apply for saving signal
    bool                                           fDiagnosticOutput;       //< If true will spew endless messages to output
         
    std::vector<char>                              fFilterModeVec;          //< Allowed modes for the filter

    using FragmentIDPair = std::pair<unsigned int, unsigned int>;
    using FragmentIDVec  = std::vector<FragmentIDPair>;
    using FragmentIDMap  = std::map<unsigned int, unsigned int>;

    FragmentIDMap                                  fFragmentIDMap;

    // Allocate containers for noise processing
    icarus_signal_processing::VectorInt            fChannelIDVec;
    icarus_signal_processing::ArrayBool            fSelectVals;
    icarus_signal_processing::ArrayBool            fROIVals;
    icarus_signal_processing::ArrayFloat           fRawWaveforms;
    icarus_signal_processing::ArrayFloat           fPedCorWaveforms;
    icarus_signal_processing::ArrayFloat           fIntrinsicRMS;
    icarus_signal_processing::ArrayFloat           fCorrectedMedians;
    icarus_signal_processing::ArrayFloat           fWaveLessCoherent;
    icarus_signal_processing::ArrayFloat           fMorphedWaveforms;
         
    icarus_signal_processing::VectorFloat          fPedestalVals;
    icarus_signal_processing::VectorFloat          fFullRMSVals;
    icarus_signal_processing::VectorFloat          fTruncRMSVals;
    icarus_signal_processing::VectorInt            fNumTruncBins;
    icarus_signal_processing::VectorInt            fRangeBins;
   
    icarus_signal_processing::VectorFloat          fThresholdVec;
   
    std::vector<unsigned int>                      fPlaneVec;
       
    const geo::Geometry*                           fGeometry;              //< pointer to the Geometry service
    const icarusDB::IICARUSChannelMap*             fChannelMap;

    // Keep track of the FFT 
    icarus_signal_processing::FFTFilterFunctionVec fFFTFilterFunctionVec;

};

TPCDecoderFilter2D::TPCDecoderFilter2D(fhicl::ParameterSet const &pset)
{
    this->configure(pset);

    fSelectVals.clear();
    fROIVals.clear();
    fRawWaveforms.clear();
    fPedCorWaveforms.clear();
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
    fFragment_id_offset     = pset.get<uint32_t           >("fragment_id_offset"      );
    fSigmaForTruncation     = pset.get<float              >("NSigmaForTrucation",  3.5);
    fCoherentNoiseGrouping  = pset.get<size_t             >("CoherentGrouping",    64);
    fCoherentNoiseOffset    = pset.get<size_t             >("CoherentOffset",       0);
    fStructuringElement     = pset.get<std::vector<size_t>>("StructuringElement",  std::vector<size_t>()={8,16});
    fMorphWindow            = pset.get<size_t             >("FilterWindow",        10);
    fThreshold              = pset.get<std::vector<float> >("Threshold",           std::vector<float>()={5.0,3.5,3.5});
    fDiagnosticOutput       = pset.get<bool               >("DiagnosticOutput",    false);
    fFilterModeVec          = pset.get<std::vector<char>  >("FilterModeVec",       std::vector<char>()={'g','g','d'}); //{'d','e','g'});

    FragmentIDVec tempIDVec = pset.get< FragmentIDVec >("FragmentIDVec", FragmentIDVec());

    for(const auto& idPair : tempIDVec) fFragmentIDMap[idPair.first] = idPair.second;

    fGeometry               = art::ServiceHandle<geo::Geometry const>{}.get();
    fChannelMap             = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

//    std::vector<double> highPassSigma = {3.5, 3.5, 0.5};
//    std::vector<double> highPassCutoff = {12., 12., 2.};
//  So we build a filter kernel for convolution with the waveform working in "tick" space. 
//  For translation, each "tick" is approximately 0.61 kHz... the frequency response functions all
//  are essentially zero by 500 kHz which is like 800 "ticks". 
    std::vector<std::pair<double,double>> windowSigma  = {{1.5,20.}, {1.5,20.}, {2.0,20.}};
    std::vector<std::pair<double,double>> windowCutoff = {{8.,800.}, {8.,800.}, {3.0,800.}};


    for(int plane = 0; plane < 3; plane++)
    {
        fFFTFilterFunctionVec.emplace_back(std::make_unique<icarus_signal_processing::WindowFFTFilter>(windowSigma[plane], windowCutoff[plane]));
    }

    return;
}

void TPCDecoderFilter2D::process_fragment(detinfo::DetectorClocksData const&,
                                          const artdaq::Fragment &fragment)
{
    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    // convert fragment to Nevis fragment
    icarus::PhysCrateFragment physCrateFragment(fragment);

    size_t nBoardsPerFragment   = physCrateFragment.nBoards();
    size_t nChannelsPerBoard    = physCrateFragment.nChannelsPerBoard();
    size_t nSamplesPerChannel   = physCrateFragment.nSamplesPerChannel();
//    size_t nChannelsPerFragment = nBoardsPerFragment * nChannelsPerBoard;

    // Recover the Fragment id:
    artdaq::detail::RawFragmentHeader::fragment_id_t fragmentID = fragment.fragmentID();

    if (fDiagnosticOutput) std::cout << "==> Recovered fragmentID: " << std::hex << fragmentID << std::dec << " ";

    // Look for special case of diagnostic running
    if (!fChannelMap->hasFragmentID(fragmentID))
    {
        if (fFragmentIDMap.find(fragmentID) == fFragmentIDMap.end()) //throw std::runtime_error("You can't save yourself");
        {
            theClockTotal.stop();
            if (fDiagnosticOutput) std::cout << " **** no match found ****" << std::endl;

            return;
        }

        if (fDiagnosticOutput) std::cout << "No match, use fhicl list? Have fragmentID: " << fragmentID << ", make it: " << std::hex << fFragmentIDMap[fragmentID] << std::dec << std::endl;

        fragmentID = fFragmentIDMap[fragmentID];

        if (!fChannelMap->hasFragmentID(fragmentID))
        {
            if (fDiagnosticOutput) std::cout << "WTF? This really can't happen, right?" << std::endl;
            return;
        }

    }

    if (fDiagnosticOutput) std::cout << std::endl;

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

    // Make sure these always get defined to be as large as can be
    const size_t maxChannelsPerFragment(576);

    if (fSelectVals.empty())        fSelectVals       = icarus_signal_processing::ArrayBool(maxChannelsPerFragment,  icarus_signal_processing::VectorBool(nSamplesPerChannel));
    if (fROIVals.empty())           fROIVals          = icarus_signal_processing::ArrayBool(maxChannelsPerFragment,  icarus_signal_processing::VectorBool(nSamplesPerChannel));
    if (fRawWaveforms.empty())      fRawWaveforms     = icarus_signal_processing::ArrayFloat(maxChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fPedCorWaveforms.empty())   fPedCorWaveforms  = icarus_signal_processing::ArrayFloat(maxChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fIntrinsicRMS.empty())      fIntrinsicRMS     = icarus_signal_processing::ArrayFloat(maxChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fCorrectedMedians.empty())  fCorrectedMedians = icarus_signal_processing::ArrayFloat(maxChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fWaveLessCoherent.empty())  fWaveLessCoherent = icarus_signal_processing::ArrayFloat(maxChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));
    if (fMorphedWaveforms.empty())  fMorphedWaveforms = icarus_signal_processing::ArrayFloat(maxChannelsPerFragment, icarus_signal_processing::VectorFloat(nSamplesPerChannel));

    if (fChannelIDVec.empty())      fChannelIDVec     = icarus_signal_processing::VectorInt(maxChannelsPerFragment);
    if (fPedestalVals.empty())      fPedestalVals     = icarus_signal_processing::VectorFloat(maxChannelsPerFragment);
    if (fFullRMSVals.empty())       fFullRMSVals      = icarus_signal_processing::VectorFloat(maxChannelsPerFragment);
    if (fTruncRMSVals.empty())      fTruncRMSVals     = icarus_signal_processing::VectorFloat(maxChannelsPerFragment);
    if (fNumTruncBins.empty())      fNumTruncBins     = icarus_signal_processing::VectorInt(maxChannelsPerFragment);
    if (fRangeBins.empty())         fRangeBins        = icarus_signal_processing::VectorInt(maxChannelsPerFragment);

    if (fThresholdVec.empty())      fThresholdVec     = icarus_signal_processing::VectorFloat(maxChannelsPerFragment);

    if (fPlaneVec.empty())          fPlaneVec.resize(nChannelsPerBoard,0);
   
    // Allocate the de-noising object
//    icarus_signal_processing::Denoiser2D_Hough     denoiser;
    icarus_signal_processing::WaveformTools<float> waveformTools;

    cet::cpu_timer theClockPedestal;

    theClockPedestal.start();

    // The first task is to recover the data from the board data block, determine and subtract the pedestals
    // and store into vectors useful for the next steps
    for(size_t board = 0; board < boardIDVec.size(); board++)
    {
        const icarusDB::ChannelPlanePairVec& channelPlanePairVec = fChannelMap->getChannelPlanePair(boardIDVec[board]);

        uint32_t boardSlot = physCrateFragment.DataTileHeader(board)->StatusReg_SlotID();

        if (fDiagnosticOutput)
        {
            std::cout << "********************************************************************************" << std::endl;
            std::cout << "FragmentID: " << std::hex << fragmentID << ", Crate: " << crateName << std::dec << ", boardID: " << boardSlot << "/" << nBoardsPerFragment << ", size " << channelPlanePairVec.size() << "/" << nChannelsPerBoard << ", ";
            std::cout << std::endl;
        }

        // This is where we would recover the base channel for the board from database/module
        size_t boardOffset = nChannelsPerBoard * board;

        // Copy to input data array
        for(size_t chanIdx = 0; chanIdx < nChannelsPerBoard; chanIdx++)
        {
            // Get the channel number on the Fragment
            size_t channelOnBoard = boardOffset + chanIdx;

            icarus_signal_processing::VectorFloat& rawDataVec = fRawWaveforms[channelOnBoard];
            for (size_t tick = 0; tick < nSamplesPerChannel; ++tick)
              rawDataVec[tick] = -physCrateFragment.adc_val(board, chanIdx, tick);

            icarus_signal_processing::VectorFloat& pedCorDataVec = fPedCorWaveforms[channelOnBoard];

            // Keep track of the channel
            fChannelIDVec[channelOnBoard] = channelPlanePairVec[chanIdx].first;

            // Handle the filter function to use for this channel
            unsigned int plane = channelPlanePairVec[chanIdx].second;

            fPlaneVec[chanIdx] = plane;

            // Set the threshold for this channel
            fThresholdVec[channelOnBoard] = fThreshold[plane];

            // Now determine the pedestal and correct for it
            waveformTools.getPedestalCorrectedWaveform(rawDataVec,
                                                       pedCorDataVec,
                                                       fSigmaForTruncation,
                                                       fPedestalVals[channelOnBoard],
                                                       fFullRMSVals[channelOnBoard],
                                                       fTruncRMSVals[channelOnBoard],
                                                       fNumTruncBins[channelOnBoard],
                                                       fRangeBins[channelOnBoard]);

            // Convolve with a filter function
            (*fFFTFilterFunctionVec[plane])(pedCorDataVec);

            if (fDiagnosticOutput)
            {
                std::vector<geo::WireID> widVec = fGeometry->ChannelToWire(channelPlanePairVec[chanIdx].first);

                if (widVec.empty()) std::cout << channelPlanePairVec[chanIdx].first << "/" << chanIdx  << "=" << fFullRMSVals[channelOnBoard] << " * ";
                else std::cout << fChannelIDVec[channelOnBoard] << "-" << widVec[0].Cryostat << "/" << widVec[0].TPC << "/" << widVec[0].Plane << "/" << widVec[0].Wire << "=" << fFullRMSVals[channelOnBoard] << " * ";
            }
        }

        if (fDiagnosticOutput) std::cout << std::endl;

        // The assumption is that channels map to planes in a continuous manner. 
        // For the first induction all the channels on a board should map to the same plane
        // For middle induction and collection you map to two planes where one group of 32 will go to one plane, the other group to the other plane. 
        // So we need to go through the planeVec to understand how to break this up...
        unsigned int startChannel(0);

        while(startChannel < fPlaneVec.size())
        {
            unsigned int stopChannel = startChannel;
            unsigned int plane       = fPlaneVec[startChannel];

            while(stopChannel < fPlaneVec.size() && fPlaneVec[stopChannel] == plane) stopChannel++;

            size_t deltaChannels = stopChannel - startChannel;

            std::cout << "==> Board Offset: " << boardOffset << ", start: " << startChannel << ", stop: " << stopChannel << ", delta: " << deltaChannels << std::endl;

            if (deltaChannels >= 32)  // How can we handle this?
            {
                // Filter function
                std::unique_ptr<icarus_signal_processing::IMorphologicalFunctions2D> filterFunctionPtr; 

                switch(fFilterModeVec[plane])
                {
                    case 'd' :
                        filterFunctionPtr = std::make_unique<icarus_signal_processing::Dilation2D>(fStructuringElement[0],fStructuringElement[1]);
                        break;
                    case 'e' :
                        filterFunctionPtr = std::make_unique<icarus_signal_processing::Erosion2D>(fStructuringElement[0],fStructuringElement[1]);
                        break;
                    case 'g' :
                        filterFunctionPtr = std::make_unique<icarus_signal_processing::Gradient2D>(fStructuringElement[0],fStructuringElement[1]);
                        break;
                    case 'a' :
                        filterFunctionPtr = std::make_unique<icarus_signal_processing::Average2D>(fStructuringElement[0],fStructuringElement[1]);
                        break;
                    case 'm' :
                        filterFunctionPtr = std::make_unique<icarus_signal_processing::Median2D>(fStructuringElement[0],fStructuringElement[1],0);
                        break;
                    default:
                        std::cout << "***** FOUND NO MATCH FOR TYPE: " << fFilterModeVec[plane] << ", plane " << plane << " DURING INITIALIZATION OF FILTER FUNCTIONS IN TPCDecoderFilter2D" << std::endl;
                        break;
                }

                if (boardOffset + startChannel + deltaChannels > fWaveLessCoherent.size()) 
                {
                    std::cout << "*** Attempting to write past end of array, boardOffset: " << boardOffset << ", startChannel: " << startChannel << ", deltaChannels: " << deltaChannels << ", array size:" << fWaveLessCoherent.size() << std::endl;
                    startChannel = stopChannel;
                    continue;
                }

                icarus_signal_processing::Denoiser2D_Hough denoiser(filterFunctionPtr.get(), fThresholdVec, fCoherentNoiseGrouping, fCoherentNoiseOffset, fMorphWindow);

                // Run the coherent filter
                denoiser(fWaveLessCoherent.begin()  + boardOffset + startChannel,
                         fPedCorWaveforms.begin()   + boardOffset + startChannel,
                         fMorphedWaveforms.begin()  + boardOffset + startChannel,
                         fIntrinsicRMS.begin()      + boardOffset + startChannel,
                         fSelectVals.begin()        + boardOffset + startChannel,
                         fROIVals.begin()           + boardOffset + startChannel,
                         fCorrectedMedians.begin()  + boardOffset + startChannel,
                         deltaChannels);

                    }

            startChannel = stopChannel;
        }

    }

    // We need to make sure the channelID information is not preserved when less than 9 boards in the fragment
    if (boardIDVec.size() < 9)
    {
        std::fill(fChannelIDVec.begin() + boardIDVec.size() * nChannelsPerBoard, fChannelIDVec.end(), -1);
    }

    theClockPedestal.stop();

    double pedestalTime = theClockPedestal.accumulated_real_time();

    cet::cpu_timer theClockDenoise;

    theClockDenoise.start();

    theClockDenoise.stop();

    double denoiseTime = theClockDenoise.accumulated_real_time();

    theClockDenoise.start();

    theClockDenoise.stop();

    double cohPedSubTime = theClockDenoise.accumulated_real_time() - denoiseTime;


    theClockTotal.stop();

    double totalTime = theClockTotal.accumulated_real_time();

    mf::LogInfo("TPCDecoderFilter2D") << "    *totalTime: " << totalTime << ", pedestal: " << pedestalTime << ", noise: " << denoiseTime << ", ped cor: " << cohPedSubTime << std::endl;

    return;
}


DEFINE_ART_CLASS_TOOL(TPCDecoderFilter2D)
} // namespace lar_cluster3d
