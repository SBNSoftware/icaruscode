////////////////////////////////////////////////////////////////////////
/// \file  MCTruthParticleHistory.h
/// \brief A "chain" of particles associated with production of a Particle in a ParticleList.
///
/// \version $Id: MCTruthParticleHistory.h,v 1.1 2010/04/29 15:38:01 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///

/// A container for a chain of particles in an event. It's a meant as
/// a convenience for looking at a sequence of particles within a
/// sim::ParticeList.

/// Here's an example to illustrate the purpose and use of a
/// list. Assume a pi+ is a primary particle in an event whose decay
/// is modeled like this:

/// TrackID   Particle
///     2       pi+
///   101       nu_mu
///   102       mu+
///   341       nu_mu_bar
///   342       nu_e
///   343       e+   

/// I'm playing around with the ParticleList for the event, and I'm
/// interested in what produced track ID 343, which is an e+. I can
/// use the MCTruthParticleHistory class to quickly go through the production
/// chain:

/// sim::ParticleList* particleList = // ... from somewhere
/// int trackID = 343;
/// const sim::MCTruthParticleHistory MCTruthParticleHistory( particleList, trackID );
/// for ( int i = 0; i != MCTruthParticleHistory.size(); ++i )
/// { 
///   const simb::MCParticle* particle = MCTruthParticleHistory[i];
///   // ... 
/// }

/// In the above example:
/// MCTruthParticleHistory.size() == 3
/// MCTruthParticleHistory[0] points to the particle with track ID 2
/// MCTruthParticleHistory[1] points to the particle with track ID 102
/// MCTruthParticleHistory[2] points to the particle with track ID 343

/// So as you go through a MCTruthParticleHistory "array," the first element
/// is a primary particle in the event, and the last element is the
/// particle you used to create the history.

/// MCTruthParticleHistory looks like a vector< const simb::MCParticle* >, with the
/// following additions:

/// - a ParticleList() method that returns a ParticleList* to the
///   object that's associated with the history.

/// - an EndParticleID() method that returns the track ID of the last
///   particle in the chain; that is, it's the second argument in the
///   constructor.

/// - operator<< method for ROOT display and ease of debugging.

/// TECHNICAL NOTES:

/// MCTruthParticleHistory behaves mostly like a vector, but it's actually a
/// deque. This means that you can't assume that &MCTruthParticleHistory[0]
/// is a continugous array of Particle*. If those two sentences mean
/// nothing to you, don't worry about it; this only matters to folks
/// familiar with STL.

/// A given MCTruthParticleHistory object is associated with the ParticleList
/// used to create it. If you delete the ParticleList (by reading in a
/// new event, for example) then the contents of the corresponding
/// MCTruthParticleHistory object(s) are garbage.

/// If you create a MCTruthParticleHistory object like this:
///    const sim::MCTruthParticleHistory ph(particleList,1123);
/// and there is no track 1123 in the particle list, then ph.size()==0.

/// MCTruthParticleHistory[0] is not necessarily a primary particle in the
/// event. It's possible for a production chain to be broken due to
/// simulation cuts. The first element just represents as far back we
/// can go in the production chain given the ParticleList.


#ifndef TRUTH_MCTruthParticleHistory_H
#define TRUTH_MCTruthParticleHistory_H

#include "nusimdata/SimulationBase/MCParticle.h"

#include <deque>
#include <iostream>

namespace truth {

// Forward declaration
class ParticleList;

class MCTruthParticleHistory : public std::deque< const simb::MCParticle* >
{
public:

    // Constructor and destructor
    MCTruthParticleHistory( const MCTruthParticleList* list, const int trackID );
    virtual ~MCTruthParticleHistory();

    // For which particle was this history generated?
    int EndParticleID() const { return m_trackID; }

    // With which ParticleList is this history associated?
    const MCTruthParticleList* ParticleList() const { return m_particleList; }

    friend std::ostream& operator<< ( std::ostream& output, const MCTruthParticleHistory& );

private:
    const MCTruthParticleList* m_particleList; ///> The ParticleList associated with this chain.
    int                        m_trackID;      ///> The particle for which a history was created.
};

} // namespace sim

#endif // SIM_MCTruthParticleHistory_H
