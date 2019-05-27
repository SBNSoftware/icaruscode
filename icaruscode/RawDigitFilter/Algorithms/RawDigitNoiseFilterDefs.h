#ifndef RAWDIGITNOISEFILTERDEFS_H
#define RAWDIGITNOISEFILTERDEFS_H
////////////////////////////////////////////////////////////////////////
//
// File:        RawDigitNoiseFilterDefs.h
//
//              This file provides data structure definitions for filtering
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

namespace caldata
{
    // Provide definitions of the raw waveforms for internal use
    using RawDigitVector           = raw::RawDigit::ADCvector_t;
    using RawDigitVectorItr        = RawDigitVector::iterator;
    using RawDigitVectorItrPair    = std::pair<RawDigitVectorItr,RawDigitVectorItr>;
    
    // The following set up the organization of handling of channels for
    // the classifying waveforms and handling the correlated noise correction
    using RawDigitVectorIdxPair    = std::pair<size_t,size_t>;
    using WireToAdcIdxMap          = std::map<size_t, RawDigitVectorIdxPair>;
    using WireToRawDigitVecPair    = std::pair<size_t, RawDigitVector&>;
    using WireToRawDigitVecMap     = std::map<size_t, RawDigitVector&>;
    
    using RawDigitAdcIdxPair       = std::pair<WireToRawDigitVecMap,WireToAdcIdxMap>;
    using GroupToDigitIdxPairMap   = std::map<size_t,RawDigitAdcIdxPair>;

} // end caldata namespace
#endif
