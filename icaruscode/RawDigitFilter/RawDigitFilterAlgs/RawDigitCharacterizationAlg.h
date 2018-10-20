#ifndef RAWDIGITCHARACTERIZATIONALG_H
#define RAWDIGITCHARACTERIZATIONALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitCharacterizationAlg
// Module Type: producer
// File:        RawDigitCharacterizationAlg.h
//
//              The intent of this module is to provide methods for
//              characterizing an input RawDigit waveform
//
// Configuration parameters:
//
// TruncMeanFraction     - the fraction of waveform bins to discard when
//                         computing the means and rms
// RMSRejectionCutHi     - vector of maximum allowed rms values to keep channel
// RMSRejectionCutLow    - vector of lowest allowed rms values to keep channel
// RMSSelectionCut       - vector of rms values below which to not correct
// MaxPedestalDiff       - Baseline difference to pedestal to flag
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
// Based on work done by Brian Kirby, Mike Mooney and Jyoti Joshi
//
////////////////////////////////////////////////////////////////////////

#include "RawDigitNoiseFilterDefs.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include "ChannelGroups.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace caldata
{
class RawDigitCharacterizationAlg
{
public:

    // Copnstructors, destructor.
    RawDigitCharacterizationAlg(fhicl::ParameterSet const & pset);
    ~RawDigitCharacterizationAlg();

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeHists(art::ServiceHandle<art::TFileService>&);
    
    // Basic waveform mean and rms
    void getMeanAndRms(const RawDigitVector& rawWaveform,
                       float&                aveVal,
                       float&                rmsVal,
                       int&                  numBins) const;
    
    // Basic waveform mean and rms plus trunated rms
    void getMeanAndTruncRms(const RawDigitVector& rawWaveform,
                            float&                aveVal,
                            float&                rmsVal,
                            float&                rmsTrunc,
                            int&                  numBins) const;

    // Truncated rms calculation
    void getTruncatedRMS(const RawDigitVector& rawWaveform,
                         float&                pedestal,
                         float&                truncRms) const;
   
    // Basic waveform mean, rms and pedestal offset
    void getMeanRmsAndPedCor(const RawDigitVector& rawWaveform,
                             unsigned int          channel,
                             unsigned int          view,
                             unsigned int          wire,
                             float&                aveVal,
                             float&                rmsVal,
                             float&                pedCorVal) const;
    
    // Basic waveform mean, rms and pedestal offset
    void getWaveformParams(const RawDigitVector& rawWaveform,
                           unsigned int          channel,
                           unsigned int          view,
                           unsigned int          wire,
                           float&                truncMean,
                           float&                truncRms,
                           short&                mean,
                           short&                median,
                           short&                mode,
                           float&                skewness,
                           float&                rms,
                           short&                minMax,
                           float&                neighborRatio,
                           float&                pedCorVal) const;
    
    bool classifyRawDigitVec(RawDigitVector&         rawWaveform,
                             unsigned int            viewIdx,
                             unsigned int            wire,
                             float                   truncRms,
                             short                   minMax,
                             short                   mean,
                             float                   skewness,
                             float                   neighborRatio,
                             GroupToDigitIdxPairMap& groupToDigitIdxPairMap) const;

    template<class T> T getMedian(std::vector<T>&, T) const;
    
private:

    // Fcl parameters.
    float                fTruncMeanFraction;     ///< Fraction for truncated mean
    std::vector<float>   fRmsRejectionCutHi;     ///< Maximum rms for input channels, reject if larger
    std::vector<float>   fRmsRejectionCutLow;    ///< Minimum rms to consider channel "alive"
    std::vector<float>   fRmsSelectionCut;       ///< Don't use/apply correction to wires below this
    std::vector<short>   fMinMaxSelectionCut;    ///< Plane by plane cuts for spread cut
    unsigned int         fTheChosenWire;         ///< For example hist
    double               fMaxPedestalDiff;       ///< Max pedestal diff to db to warn
    std::vector<size_t>  fHistsWireGroup;        ///< Wire Group to pick on
    std::vector<size_t>  fNumWiresToGroup;       ///< If smoothing, the number of wires to look at
    bool                 fFillHistograms;        ///< if true then will fill diagnostic hists
    
    // Make sure hists for this instance are initialized
    bool                 fHistsInitialized;
    
    // Pointers to the histograms we'll create for monitoring what is happening
    TH1D*                fAdcCntHist[3];
    TH1D*                fAveValHist[3];
    TH1D*                fRmsTValHist[3];
    TH1D*                fRmsFValHist[3];
    TH1D*                fPedValHist[3];
    TH1D*                fAverageHist[3];
    TProfile*            fRmsValProf[3];
    TProfile*            fMinMaxValProf[3];
    TProfile*            fPedValProf[3];
    
    std::vector<TProfile*>              fMinMaxProfiles;
    std::vector<TProfile*>              fSkewnessProfiles;
    std::vector<TProfile*>              fModeRatioProfiles;

    caldata::ChannelGroups              fChannelGroups;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    art::ServiceHandle<geo::Geometry>            fGeometry;             ///< pointer to Geometry service
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();   ///< Detector properties service
    const lariov::DetPedestalProvider&          fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
};

} // end of namespace caldata

#endif
