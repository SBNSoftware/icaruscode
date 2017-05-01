#ifndef RAWDIGITFILTERALG_H
#define RAWDIGITFILTERALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitFilterAlg
// Module Type: producer
// File:        RawDigitFilterAlg.h
//
//              This module provides methods for applying filters to
//              input waveforms
//
// Configuration parameters:
//
// RunFFTCorrectedWires  - FFT analyze the output RawDigits if true - diagnostics
// TruncateTicks:        - Provide mechanism to truncate a readout window to a smaller size
// WindowSize:           - The desired size of the output window
// NumTicksToDropFront:  - The number ticks dropped off the front of the original window
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
// Based on work by Brian Kirby, Mike Mooney and Jyoti Joshi
//
////////////////////////////////////////////////////////////////////////

#include "RawDigitNoiseFilterDefs.h"
#include "RawDigitCharacterizationAlg.h"
#include "RawDigitBinAverageAlg.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <set>

namespace caldata
{

class RawDigitFilterAlg
{
public:

    // Copnstructors, destructor.
    RawDigitFilterAlg(fhicl::ParameterSet const & pset);
    ~RawDigitFilterAlg();

    // Overrides.
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeHists(art::ServiceHandle<art::TFileService>&);
    
    void doTopHatFilter(RawDigitVector&, size_t wire) const;
    
    void doAdaptiveFilter(RawDigitVector&) const;
    
private:

    // Fcl parameters.
    float                  fTruncMeanFraction;     ///< Fraction for truncated mean
    size_t                 fStructuringElement;    ///< Structuring element to use with Top Hat filter
    bool                   fFillHistograms;        ///< if true then will fill diagnostic hists
    
    // Pointers to the histograms we'll create for monitoring what is happening
    std::vector<TProfile*> fErosionVecHists;
    std::vector<TProfile*> fDilationVecHists;
    
    bool                   fFirstEvent;
    bool                   fHistsInitialized;
    
    std::vector<std::set<size_t>> fBadWiresbyViewAndWire;
    
    // We'll use this algorithm internally here too
    RawDigitCharacterizationAlg fCharacterizationAlg;
    RawDigitBinAverageAlg       fBinAverageAlg;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    art::ServiceHandle<geo::Geometry>            fGeometry;             ///< pointer to Geometry service
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();   ///< Detector properties service
    const lariov::DetPedestalProvider&          fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
};

} // end of caldata namespace

#endif
