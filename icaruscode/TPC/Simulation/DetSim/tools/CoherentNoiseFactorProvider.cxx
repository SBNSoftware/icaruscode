////////////////////////////////////////////////////////////////////////
/// \file   icaruscode/TPC/Simulation/DetSim/tools/CoherentNoiseFactorProvider.cxx
/// \author T. Usher (factorized by Gianluca Petrillo, petrillo@slac.stanford.edu)
/// \see    icaruscode/TPC/Simulation/DetSim/tools/CoherentNoiseFactorProvider.h
////////////////////////////////////////////////////////////////////////

// library header
#include "icaruscode/TPC/Simulation/DetSim/tools/CoherentNoiseFactorProvider.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"

#include "icaruscode/TPC/Simulation/DetSim/tools/ICoherentNoiseFactor.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"

#include "TH1D.h"

#include <string>
#include <iostream>

namespace Noise
{

//----------------------------------------------------------------------
// Constructor.
CoherentNoiseFactorProvider::CoherentNoiseFactorProvider(const fhicl::ParameterSet& pset) {

    mf::LogInfo("CoherentNoiseFactorProvider") << "Recover the channel map" ;

    fDiagnosticOutput = pset.get<bool>("DiagnosticOutput", false);

    // Set up the board->correlated factors map
    const auto& wireReadout = art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}.get();

    const icarusDB::TPCReadoutBoardToChannelMap& readoutBoardToChannelMap = wireReadout->getReadoutBoardToChannelMap();

    for(const auto& boardPair : readoutBoardToChannelMap)
    {
        fCorrFactorsMap.insert({boardPair.first,std::vector<float>(4,0.)});
    }


    mf::LogInfo("CoherentNoiseFactorProvider") << "==> FragmentID map size: " << fCorrFactorsMap.size() << std::endl;
    
    return;
}
  
    
    // Reset factors
void CoherentNoiseFactorProvider::resetCoherentNoiseFactors(const TH1D* noiseHist)
{
    for(size_t index = 0; index < 4; index++)
    {
        float meanVal = noiseHist->GetMean();
    
        for(auto& correction : fCorrFactorsMap)
        {
            float corVal = noiseHist->GetRandom() / meanVal;

            correction.second[index] = corVal;
        }
    }

    return;
}

    // Recover the coherent noise factor 
float CoherentNoiseFactorProvider::getCoherentNoiseFactor(unsigned int board, unsigned int index) const
{
    CorrFactorsMap::const_iterator cohFactorsItr = fCorrFactorsMap.find(board);

    if (cohFactorsItr == fCorrFactorsMap.end())
        throw cet::exception("CoherentNoiseFactorService")
          << "Readout board " << board << ", with index " << index << "'\n"
          << "This is considered a fatal issue!\n";

    float cf = cohFactorsItr->second[index];

    return cf;
}


} // end namespace

