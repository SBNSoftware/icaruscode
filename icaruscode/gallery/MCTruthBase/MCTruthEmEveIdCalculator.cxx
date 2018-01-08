////////////////////////////////////////////////////////////////////////
///
/// ********************************************************************
/// *** This class has been lifted from nutools and modified so that it
/// *** does NOT take ownwership of the MCParticles that are added to it
/// *** so we can avoid duplicating them
/// ********************************************************************
///
/// \file  MCTruthEmEveIdCalculator.cxx
/// \brief Example routine for calculating the "ultimate e-m mother" of a particle in a simulated event.
///
/// \version $Id: MCTruthEmEveIdCalculator.cxx,v 1.1 2010/05/13 16:12:20 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "icaruscode/gallery/MCTruthBase/MCTruthEmEveIdCalculator.h"
#include "icaruscode/gallery/MCTruthBase/MCTruthParticleList.h"
#include "icaruscode/gallery/MCTruthBase/MCTruthParticleHistory.h"

#include <TString.h>

namespace truth {

  //----------------------------------------------------------------------------
  // This particular class attempts to find the "ultimate mother" for
  // electromagnetic showers. It goes up the chain of particles in an
  // event, until it encounters a particle that is either primary or
  // was not produced by a "trivial" e-m process.
  int MCTruthEmEveIdCalculator::DoCalculateEveId( const int trackID )
  {
    // Almost any eve ID calculation will use this: Get the entire
    // history of the particle and its ancestors in the simulated
    // event. (m_particleList is defined in EveIdCalculator.h)
    const MCTruthParticleHistory particleHistory( m_particleList, trackID );

    // You can treat particleHistory like an array, and do something
    // like:

    // for ( int i = particleHistory.size(); i >= 0; --i ) 
    // { const simb::MCParticle* particle = particleHistory[i]; ... }

    // But we know how to use the Standard Template Library (Yes, you
    // do! Don't doubt yourself!) so let's treat particleHistory in
    // the most efficient manner, as an STL container. We want to scan
    // the container from its last element to its first, from the base
    // of the particle production chain towards the primary particle.

    for(auto i = particleHistory.rbegin(); i != particleHistory.rend(); ++i){
      // Get the text string that describes the process that created
      // the particle.
      std::string process = (*i)->Process();

      // Skip it if it was created by pair production, compton
      // scattering, photoelectric effect, bremstrahlung,
      // annihilation, or any ionization. (The ultimate source of
      // the process names are the physics lists used in Geant4.)
      
      if ( process.find("conv")              != std::string::npos ||
	   process.find("LowEnConversion")   != std::string::npos ||
	   process.find("Pair")              != std::string::npos ||
	   process.find("compt")             != std::string::npos ||
	   process.find("Compt")             != std::string::npos ||
	   process.find("Brem")              != std::string::npos ||
	   process.find("phot")              != std::string::npos ||
	   process.find("Photo")             != std::string::npos ||
	   process.find("Ion")               != std::string::npos ||
	   process.find("annihil")           != std::string::npos) continue;         

      // If we get here, the particle wasn't created by any of the
      // above processes. Return its ID.
      return (*i)->TrackId();
    }

    // If we get here, we've skipped every particle in the
    // chain. Perhaps it was empty.
    return 0;
  }

} // namespace sim
