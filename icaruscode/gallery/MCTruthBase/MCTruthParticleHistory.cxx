////////////////////////////////////////////////////////////////////////
/// \file  MCTruthParticleHistory.cxx
/// \brief A "chain" of particles associated with production of a Particle in a ParticleList.
///
/// \version $Id: MCTruthParticleHistory.cxx,v 1.1 2010/04/29 15:38:01 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "icaruscode/gallery/MCTruthBase/MCTruthParticleList.h"
#include "icaruscode/gallery/MCTruthBase/MCTruthParticleHistory.h"

#include <cmath>

namespace truth {

//----------------------------------------------------------------------------
// Nothing special need be done for the constructor.
MCTruthParticleHistory::MCTruthParticleHistory( const MCTruthParticleList* list, const int trackID )
    : m_particleList(list)
    , m_trackID(trackID)
{
    // Look for the track in the particle list.
    MCTruthParticleList::const_iterator search = m_particleList->find( m_trackID );

    // While we're still finding particles in the chain...
    while ( search != m_particleList->end() ){
        const simb::MCParticle* particle = (*search).second;
        push_front( particle );
      
        // If this is a primary particle, we're done.
        int trackID = particle->TrackId();
        if ( m_particleList->IsPrimary( trackID ) ) break;
      
        // Now look for the parent of this particle.
        int parentID = particle->Mother();
        search = m_particleList->find( parentID );
      
    } // while we're finding particles in the chain
}

//----------------------------------------------------------------------------
// Nothing special for the destructor.
MCTruthParticleHistory::~MCTruthParticleHistory()
{}

//----------------------------------------------------------------------------
std::ostream& operator<< ( std::ostream& output, const MCTruthParticleHistory& list )
{
    // Determine a field width for the particle number.
    MCTruthParticleHistory::size_type numberOfParticles = list.size();
    int numberOfDigits = (int) std::log10( (double) numberOfParticles ) + 1;

    // A simple header.
    output.width( numberOfDigits );
    output << "#" << ": < ID, particle >" << "\n"; 

    // Write each particle on a separate line.
    MCTruthParticleHistory::size_type nParticle = 0;
    for ( MCTruthParticleHistory::const_iterator particle = list.begin();
         particle != list.end(); ++particle, ++nParticle ){
        output.width( numberOfDigits );
        output << nParticle << ": "
               << (*particle)
               << "\n";
    }
    
    return output;
}
  
} // namespace sim
