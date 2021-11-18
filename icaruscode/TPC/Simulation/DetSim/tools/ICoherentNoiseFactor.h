///////////////////////////////////////////////////////////////////////
///
/// \file   ICoherentNoiseFactor
///
/// \brief  Interface class to recover the coherent noise scale factor
///         on a readout board basis (which spans planes)
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef ICoherentNoiseFactor_H
#define ICoherentNoiseFactor_H

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include <vector>
#include <string>

class TH1D;

namespace Noise {

class ICoherentNoiseFactor //: private lar::EnsureOnlyOneSchedule
{
public:
    virtual ~ICoherentNoiseFactor() noexcept = default;

    // Reset factors
    virtual void resetCoherentNoiseFactors(const TH1D*) = 0;

    // Recover the coherent noise factor 
    virtual float getCoherentNoiseFactor(unsigned int, unsigned int) const = 0;
};

} // end of namespace

DECLARE_ART_SERVICE_INTERFACE(Noise::ICoherentNoiseFactor, SHARED)

#endif
