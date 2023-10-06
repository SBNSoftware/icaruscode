/**
 *  @file   TPCNoiseFilter1DMC_tool.cc
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

#include "icaruscode/Decode/DecoderTools/INoiseFilter.h"

#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Denoising.h"
#include "icarus_signal_processing/Filters/FFTFilterFunctions.h"

// Eigen includes
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/Jacobi"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  TPCNoiseFilter1DMC class definiton
 */
class TPCNoiseFilter1DMC : virtual public INoiseFilter
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit TPCNoiseFilter1DMC(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TPCNoiseFilter1DMC();

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
                                  const daq::INoiseFilter::ChannelPlaneVec&,
                                  const icarus_signal_processing::ArrayFloat&,
                                  const size_t&) override;

    /**
     *  @brief Recover the channels for the processed fragment
     */
    const icarus_signal_processing::VectorInt&  getChannelIDs()       const override {return fChannelIDVec;}

    /**
     *  @brief Recover the selection values
     */
    const icarus_signal_processing::ArrayBool&  getSelectionVals()    const override {return fSelectVals;};

    /**
     *  @brief Recover the ROI values
     */
    const icarus_signal_processing::ArrayBool&  getROIVals()          const override {return fROIVals;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat& getRawWaveforms()     const override {return fRawWaveforms;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat& getPedCorWaveforms()  const override {return fPedCorWaveforms;};

    /**
     *  @brief Recover the "intrinsic" RMS
     */
    const icarus_signal_processing::ArrayFloat& getIntrinsicRMS()     const override {return fIntrinsicRMS;};

    /**
     *  @brief Recover the correction median values
     */
    const icarus_signal_processing::ArrayFloat& getCorrectedMedians()  const override {return fCorrectedMedians;};

    /**
     *  @brief Recover the waveforms less coherent noise
     */
    const icarus_signal_processing::ArrayFloat& getWaveLessCoherent()  const override {return fWaveLessCoherent;};

    /**
     *  @brief Recover the morphological filter waveforms
     */
    const icarus_signal_processing::ArrayFloat& getMorphedWaveforms()  const override {return fMorphedWaveforms;};

    /**
     *  @brief Recover the pedestals for each channel
     */
    const icarus_signal_processing::VectorFloat& getPedestalVals()     const override {return fPedestalVals;};

    /**
     *  @brief Recover the full RMS before coherent noise
     */
    const icarus_signal_processing::VectorFloat& getFullRMSVals()      const override {return fFullRMSVals;};

    /**
     *  @brief Recover the truncated RMS noise
     */
    const icarus_signal_processing::VectorFloat& getTruncRMSVals()     const override {return fTruncRMSVals;};

    /**
     *  @brief Recover the number of bins after truncation
     */
    const icarus_signal_processing::VectorInt&   getNumTruncBins()     const override {return fNumTruncBins;};

private:

    void principleComponents(const icarus_signal_processing::VectorFloat&) const;

    using FloatPairVec = std::vector<std::pair<float,float>>;

    float                                          fSigmaForTruncation;     //< Selection cut for truncated rms calculation
    size_t                                         fCoherentNoiseOffset;    //< Offset for midplane
    size_t                                         fStructuringElement;     //< Structuring element for morphological filter
    size_t                                         fMorphWindow;            //< Window for filter
    std::vector<float>                             fThreshold;              //< Threshold to apply for saving signal
    bool                                           fUseFFTFilter;           //< Turn on/off the use of the FFT filter
    bool                                           fDiagnosticOutput;       //< If true will spew endless messages to output
    FloatPairVec                                   fFFTSigmaValsVec;        //< Gives the sigmas for the filter function
    FloatPairVec                                   fFFTCutoffValsVec;       //< Gives the cutoffs for the filter function
      
    std::vector<std::string>                       fFilterModeVec;          //< Allowed modes for the filter

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

    icarus_signal_processing::FilterFunctionVec    fFilterFunctionVec;
    
    const geo::Geometry*                           fGeometry;              //< pointer to the Geometry service

    // Keep track of the FFT 
    icarus_signal_processing::FFTFilterFunctionVec fFFTFilterFunctionVec;
};

TPCNoiseFilter1DMC::TPCNoiseFilter1DMC(fhicl::ParameterSet const &pset)
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

TPCNoiseFilter1DMC::~TPCNoiseFilter1DMC()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void TPCNoiseFilter1DMC::configure(fhicl::ParameterSet const &pset)
{
    fSigmaForTruncation    = pset.get<float                   >("NSigmaForTrucation",  3.5);
    fCoherentNoiseOffset   = pset.get<size_t                  >("CoherentOffset",        0);
    fStructuringElement    = pset.get<size_t                  >("StructuringElement",   20);
    fMorphWindow           = pset.get<size_t                  >("FilterWindow",         10);
    fThreshold             = pset.get<std::vector<float>      >("Threshold",           std::vector<float>()={5.0,3.5,3.5});
    fUseFFTFilter          = pset.get<bool                    >("UseFFTFilter",        true);
    fDiagnosticOutput      = pset.get<bool                    >("DiagnosticOutput",    false);
    fFilterModeVec         = pset.get<std::vector<std::string>>("FilterModeVec",       std::vector<std::string>()={"e","g","d"});

    fFFTSigmaValsVec       = pset.get<FloatPairVec            >("FFTSigmaVals",        FloatPairVec()={{1.5,20.}, {1.5,20.}, {2.0,20.}});
    fFFTCutoffValsVec      = pset.get<FloatPairVec            >("FFTCutoffVals",       FloatPairVec()={{8.,800.}, {8.,800.}, {0.0,800.}});

    fGeometry   = art::ServiceHandle<geo::Geometry const>{}.get();

    fFFTFilterFunctionVec.clear();

    std::cout << "TPCNoiseFilter1D configure, fUseFFTFilter: " << fUseFFTFilter << std::endl;

    if (fUseFFTFilter)
    {
        for(int plane = 0; plane < 3; plane++)
        {
            //if (plane < 2) fFFTFilterFunctionVec.emplace_back(std::make_unique<icarus_signal_processing::WindowFFTFilter>(fFFTSigmaValsVec[plane], fFFTCutoffValsVec[plane]));
            //else           fFFTFilterFunctionVec.emplace_back(std::make_unique<icarus_signal_processing::NoFFTFilter>());
            fFFTFilterFunctionVec.emplace_back(std::make_unique<icarus_signal_processing::WindowFFTFilter>(fFFTSigmaValsVec[plane], fFFTCutoffValsVec[plane]));

            std::cout << "TPCDecoderFilter1D FFT setup for plane " << plane << ", SigmaVals: " << fFFTSigmaValsVec[plane].first << "/" << fFFTSigmaValsVec[plane].second << ", cutoff: " << fFFTCutoffValsVec[plane].first << "/" << fFFTCutoffValsVec[plane].second << std::endl;
        }
    }

    return;
}

void TPCNoiseFilter1DMC::process_fragment(detinfo::DetectorClocksData const&,
                                          const daq::INoiseFilter::ChannelPlaneVec&   channelPlaneVec,
                                          const icarus_signal_processing::ArrayFloat& dataArray,
                                          const size_t&                               coherentNoiseGrouping)
{
    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    // Recover the number of channels and ticks
    unsigned int numChannels = dataArray.size();
    unsigned int numTicks    = dataArray[0].size();

    if (fSelectVals.size()       < numChannels)  fSelectVals.resize(numChannels,       icarus_signal_processing::VectorBool(numTicks));
    if (fROIVals.size()          < numChannels)  fROIVals.resize(numChannels,          icarus_signal_processing::VectorBool(numTicks));
    if (fRawWaveforms.size()     < numChannels)  fRawWaveforms.resize(numChannels,     icarus_signal_processing::VectorFloat(numTicks));
    if (fPedCorWaveforms.size()  < numChannels)  fPedCorWaveforms.resize(numChannels,  icarus_signal_processing::VectorFloat(numTicks));
    if (fIntrinsicRMS.size()     < numChannels)  fIntrinsicRMS.resize(numChannels,     icarus_signal_processing::VectorFloat(numTicks));
    if (fCorrectedMedians.size() < numChannels)  fCorrectedMedians.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));
    if (fWaveLessCoherent.size() < numChannels)  fWaveLessCoherent.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));
    if (fMorphedWaveforms.size() < numChannels)  fMorphedWaveforms.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));

    if (fChannelIDVec.size()     < numChannels)  fChannelIDVec.resize(numChannels);
    if (fPedestalVals.size()     < numChannels)  fPedestalVals.resize(numChannels);
    if (fFullRMSVals.size()      < numChannels)  fFullRMSVals.resize(numChannels);
    if (fTruncRMSVals.size()     < numChannels)  fTruncRMSVals.resize(numChannels);
    if (fNumTruncBins.size()     < numChannels)  fNumTruncBins.resize(numChannels);
    if (fRangeBins.size()        < numChannels)  fRangeBins.resize(numChannels);

    size_t ngroups = std::max(numChannels/coherentNoiseGrouping,size_t(1));
    if (fThresholdVec.size()     < ngroups)  fThresholdVec.resize(ngroups);

    if (fFilterFunctionVec.size() < numChannels) fFilterFunctionVec.resize(numChannels);

//    icarus_signal_processing::Denoiser1D_Protect   denoiser;
    icarus_signal_processing::Denoiser1D           denoiser;
    icarus_signal_processing::WaveformTools<float> waveformTools(5);

    // Make a pass throught to do pedestal corrections and get raw waveform information
    for(size_t idx = 0; idx < numChannels; idx++)
    {
        icarus_signal_processing::VectorFloat& rawDataVec    = fRawWaveforms[idx];
        icarus_signal_processing::VectorFloat& pedCorDataVec = fPedCorWaveforms[idx];

        // Keep track of the channel
        fChannelIDVec[idx] = channelPlaneVec[idx].first;

        // We need to recover info on which plane we have
        std::vector<geo::WireID> widVec = fGeometry->ChannelToWire(fChannelIDVec[idx]);

        // Handle the filter function to use for this channel
        // Note the modulus... this to enable a workaround for the wirecell 2D drift which misses the channels with no signal
        unsigned int plane = channelPlaneVec[idx].second % 3;

        // Set the threshold which toggles between planes
        fThresholdVec[idx / coherentNoiseGrouping] = fThreshold[plane];

        switch(fFilterModeVec[plane][0])
        {
            case 'd' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Dilation1D>(fStructuringElement);
                break;
            case 'e' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Erosion1D>(fStructuringElement);
                break;
            case 'g' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Gradient1D>(fStructuringElement);
                break;
            case 'a' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Average1D>(fStructuringElement);
                break;
            case 'm' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Median1D>(fStructuringElement);
                break;
            default:
                std::cout << "***** FOUND NO MATCH FOR TYPE: " << fFilterModeVec[plane] << ", plane " << plane << " DURING INITIALIZATION OF FILTER FUNCTIONS IN TPCNoiseFilter1DMC" << std::endl;
                break;
        }

        std::copy(dataArray[idx].begin(),dataArray[idx].end(),rawDataVec.begin());

        // Now determine the pedestal and correct for it
        waveformTools.getPedestalCorrectedWaveform(rawDataVec,
                                                   pedCorDataVec,
                                                   fSigmaForTruncation,
                                                   fPedestalVals[idx],
                                                   fFullRMSVals[idx],
                                                   fTruncRMSVals[idx],
                                                   fNumTruncBins[idx],
                                                   fRangeBins[idx]);

        // Convolve with a filter function
        //if (fUseFFTFilter) (*fFFTFilterFunctionVec[plane])(pedCorDataVec);
        if (fUseFFTFilter)
        {
            // Temporary diagnostics
            icarus_signal_processing::VectorFloat medianSmoothVec(numTicks);

            principleComponents(pedCorDataVec);

            //waveformTools.medianSmooth(pedCorDataVec, medianSmoothVec, 201);
            waveformTools.truncAveSmooth(pedCorDataVec, medianSmoothVec, 201);

            std::transform(pedCorDataVec.begin(),pedCorDataVec.end(), medianSmoothVec.begin(), pedCorDataVec.begin(), std::minus<float>());

            std::copy(medianSmoothVec.begin(),medianSmoothVec.end(),rawDataVec.begin());
        }

        // Make sure our selection and ROI arrays are initialized
        std::fill(fSelectVals[idx].begin(),fSelectVals[idx].end(),false);
    }

    denoiser(fWaveLessCoherent.begin(),
             fPedCorWaveforms.begin(),
             fMorphedWaveforms.begin(),
             fIntrinsicRMS.begin(),
             fSelectVals.begin(),
             fROIVals.begin(),
             fCorrectedMedians.begin(),
             fFilterFunctionVec.begin(),
             fThresholdVec,
             numChannels,
             coherentNoiseGrouping,
             fCoherentNoiseOffset,
             fMorphWindow);

    theClockTotal.stop();

    double totalTime = theClockTotal.accumulated_real_time();

    mf::LogDebug("TPCNoiseFilter1DMC") << "    *totalTime: " << totalTime << std::endl;

    return;
}

void TPCNoiseFilter1DMC::principleComponents(const icarus_signal_processing::VectorFloat& waveform) const
{

    // Define elements of our covariance matrix
    float xi2(0.);
    float xiyi(0.);
    float xizi(0.0);
    float yi2(0.0);
    float yizi(0.0);
    float zi2(0.);
//    float weightSum(0.);

    std::cout << "Entering principle components alg, size: " << waveform.size() << std::endl;

    Eigen::Vector3f meanPos(float(waveform.size()/2),waveform[waveform.size()/2],0.);

    for (size_t waveIdx = 0; waveIdx < waveform.size(); waveIdx++)
    {
        float x = float(waveIdx)    - meanPos(0);
        float y = waveform[waveIdx] - meanPos(1);
        float z = 0.;

        xi2  += x * x;
        xiyi += x * y;
        xizi += x * z;
        yi2  += y * y;
        yizi += y * z;
        zi2  += z;
    }

    // Using Eigen package
    Eigen::Matrix2f sig;

//    sig << xi2, xiyi, xizi, xiyi, yi2, yizi, xizi, yizi, zi2;
    sig << xi2, xiyi, xiyi, yi2;

//    sig *= 1. / weightSum;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigenMat(sig);

    if (eigenMat.info() == Eigen::ComputationInfo::Success) 
    {
        // Now copy output
        Eigen::Matrix2f eigenVectors = eigenMat.eigenvectors().transpose();
        Eigen::Vector2f eigenValues  = eigenMat.eigenvalues();

        std::cout << "-- eigenvalues and vectors" << std::endl;
        std::cout << eigenValues << std::endl;
        std::cout << eigenVectors << std::endl;
//      // The returned eigen values and vectors will be returned in an xyz system where x is the smallest spread,
//      // y is the next smallest and z is the largest. Adopt that convention going forward
//      reco::PrincipalComponents::EigenValues recobEigenVals = eigenMat.eigenvalues().cast<float>();
//      reco::PrincipalComponents::EigenVectors recobEigenVecs =
//        eigenMat.eigenvectors().transpose().cast<float>();
//
//      // Check for a special case (which may have gone away with switch back to doubles for computation?)
//      if (std::isnan(recobEigenVals[0])) {
//        std::cout << "==> Third eigenvalue returns a nan" << std::endl;
//
//        recobEigenVals[0] = 0.;
//
//        // Assume the third axis is also kaput?
//        recobEigenVecs.row(0) = recobEigenVecs.row(1).cross(recobEigenVecs.row(2));
//      }
//
//      // Store away
//      pca = reco::PrincipalComponents(
//        true, numPairsInt, recobEigenVals, recobEigenVecs, meanPos.cast<float>());
    }

    return;
}


DEFINE_ART_CLASS_TOOL(TPCNoiseFilter1DMC)
} // namespace lar_cluster3d
