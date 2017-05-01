#ifndef RAWDIGITBINAVERAGEALG_H
#define RAWDIGITBINAVERAGEALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitBinAverageAlg
// Module Type: producer
// File:        RawDigitBinAverageAlg.h
//
//              This module provides methods for averaging bins in an
//              input RawDigitVector
//
// Configuration parameters:
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
//
////////////////////////////////////////////////////////////////////////

#include "RawDigitNoiseFilterDefs.h"
#include "fhiclcpp/ParameterSet.h"

namespace caldata
{
    
class RawDigitBinAverageAlg
{
public:

    // Copnstructors, destructor.
    RawDigitBinAverageAlg(fhicl::ParameterSet const & pset);
    ~RawDigitBinAverageAlg();
    
    void doBinAverage(RawDigitVector&, size_t) const;
    
    void doTwoBinAverage(RawDigitVector&) const;
    
private:
};

} // end caldata namespace
#endif