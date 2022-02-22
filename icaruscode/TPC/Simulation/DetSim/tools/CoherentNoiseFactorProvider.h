////////////////////////////////////////////////////////////////////////
/// \file   icaruscode/TPC/Simulation/DetSim/tools/CoherentNoiseFactorProvider.h
/// \author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
/// \see    icaruscode/TPC/Simulation/DetSim/tools/CoherentNoiseFactorProvider.cxx
////////////////////////////////////////////////////////////////////////

#ifndef ICARUSCODE_TPC_CoherentNoiseFactorPROVIDER_H
#define ICARUSCODE_TPC_CoherentNoiseFactorPROVIDER_H

// ICARUS libraries
#include "icaruscode/TPC/Simulation/DetSim/tools/ICoherentNoiseFactor.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <string>
#include <memory> // std::unique_ptr<>


// -----------------------------------------------------------------------------
namespace Noise { class CoherentNoiseFactorProvider; }
class Noise::CoherentNoiseFactorProvider : virtual public ICoherentNoiseFactor
{
public:
    
    // Constructor, destructor.
    CoherentNoiseFactorProvider(const fhicl::ParameterSet& pset);
    
    // Reset factors
    void resetCoherentNoiseFactors(const TH1D*) override;

    // Recover the coherent noise factor 
    float getCoherentNoiseFactor(unsigned int, unsigned int) const override;

private:

    bool           fDiagnosticOutput;

    using CorrFactorsMap = std::map<unsigned int, std::vector<float>>;

    CorrFactorsMap fCorrFactorsMap;

}; // Noise::CoherentNoiseFactorProvider


// -----------------------------------------------------------------------------

#endif 