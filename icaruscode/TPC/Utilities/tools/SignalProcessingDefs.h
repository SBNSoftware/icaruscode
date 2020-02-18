///////////////////////////////////////////////////////////////////////
///
/// \file   SignalProcessingDefs.h
///
/// \brief  Useful definitions 
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef SignalProcessingDefs_H
#define SignalProcessingDefs_H

#include <complex.h>

namespace icarusutil
{
    using SigProcPrecision = double;
    using ComplexVal       = std::complex<SigProcPrecision>;

    using TimeVec      = std::vector<SigProcPrecision>;
    using FrequencyVec = std::vector<ComplexVal>;
}

#endif
