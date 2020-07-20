///////////////////////////////////////////////////////////////////////
///
/// \file   IDeconvolution.h
///
/// \brief  This provides an interface for tools which are tasked with
///         running the deconvolution algorithm. This allows switching
///         between, for example, algorithms that deconvolve only ROI's
///         vs an entire wire
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IDeconvolution_H
#define IDeconvolution_H

#include "fhiclcpp/ParameterSet.h"
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IROIFinder.h"
#include "lardataobj/RecoBase/Wire.h"

namespace art
{
    class TFileDirectory;
}

namespace icarus_tool
{
    class IDeconvolution
    {
    public:
        virtual ~IDeconvolution() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)           = 0;
        virtual void initializeHistograms(art::TFileDirectory&)     const = 0;
        
        // Find the ROI's
        virtual void Deconvolve(IROIFinder::Waveform const&,
                                double samplingRate,
                                raw::ChannelID_t,
                                IROIFinder::CandidateROIVec const&,
                                recob::Wire::RegionsOfInterest_t& ) const = 0;
    };
}

#endif
