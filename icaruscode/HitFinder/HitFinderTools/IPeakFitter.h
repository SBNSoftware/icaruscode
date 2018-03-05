///////////////////////////////////////////////////////////////////////
///
/// \file   IPeakFitter.h
///
/// \brief  This provides an interface for tools which are tasked with
///         fitting peaks on input waveforms
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IICARUSPeakFitter_H
#define IICARUSPeakFitter_H

#include "fhiclcpp/ParameterSet.h"
#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"

namespace reco_tool
{
    class IPeakFitter
    {
    public:
        virtual ~IPeakFitter() noexcept = default;
        
        // Define standard art tool interface
        virtual void configure(const fhicl::ParameterSet& pset) = 0;
        
        // Define a structure to contain hits
        using PeakFitParams_t = struct PeakFitParams
        {
            float peakCenter;
            float peakCenterError;
            float peakSigma;
            float peakSigmaError;
            float peakAmplitude;
            float peakAmplitudeError;
            float peakTauLeft;
            float peakTauLeftError;
            float peakTauRight;
            float peakTauRightError;
            float peakBaseline;
            float peakBaselineError;
        };
        
        using PeakParamsVec = std::vector<PeakFitParams_t>;
        
        // Get parameters for input candidate peaks
        virtual void findPeakParameters(const std::vector<float>&,
                                        const ICandidateHitFinder::HitCandidateVec&,
                                        PeakParamsVec&,
                                        double&,
                                        int&) const = 0;
    };
}

#endif
