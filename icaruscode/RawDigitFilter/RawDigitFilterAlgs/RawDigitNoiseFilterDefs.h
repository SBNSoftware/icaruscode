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
    typedef raw::RawDigit::ADCvector_t                        RawDigitVector;
    typedef RawDigitVector::iterator                          RawDigitVectorItr;
    typedef std::pair<RawDigitVectorItr,RawDigitVectorItr>    RawDigitVectorItrPair;
    
    // The following set up the organization of handling of channels for
    // the classifying waveforms and handling the correlated noise correction
    typedef std::pair<size_t,size_t>                          RawDigitVectorIdxPair;
    typedef std::map<size_t, RawDigitVectorIdxPair>           WireToAdcIdxMap;
    typedef std::pair<size_t, RawDigitVector&>                WireToRawDigitVecPair;
    typedef std::map<size_t, RawDigitVector&>                 WireToRawDigitVecMap;
    
    typedef std::pair<WireToRawDigitVecMap,WireToAdcIdxMap>   RawDigitAdcIdxPair;
    typedef std::map<size_t,RawDigitAdcIdxPair>               GroupToDigitIdxPairMap;

} // end caldata namespace
#endif
