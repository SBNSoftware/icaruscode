///////////////////////////////////////////////////////////////////////
///
/// \file   IOpHitFinder.h
///
/// \brief  This provides an interface for tools which are tasked with
///         finding the baselines in input waveforms, primarily ROI's
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IOpHitFinder_H
#define IOpHitFinder_H

#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

namespace art
{
    class TFileDirectory;
}

namespace light
{
    class IOpHitFinder
    {
    public:
        virtual ~IOpHitFinder() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)                 = 0;
        virtual void outputHistograms(art::TFileDirectory&)               const = 0;
        
        // Find the baseline
        virtual void FindOpHits(const raw::OpDetWaveform&, recob::OpHit&) const = 0;
    };
}

#endif
