////////////////////////////////////////////////////////////////////////
///
/// ********************************************************************
/// *** This class has been lifted from nutools and modified so that it
/// *** does NOT take ownwership of the MCParticles that are added to it
/// *** so we can avoid duplicating them
/// ********************************************************************
///
/// \file  MCTruthEveIdCalculator.cxx
/// \brief Interface for calculating the "ultimate mother" of a particle in a simulated event.
///
/// \version $Id: MCTruthEveIdCalculator.cxx,v 1.1 2010/05/13 16:12:20 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "icaruscode/gallery/MCTruthBase/MCTruthEveIdCalculator.h"
#include "icaruscode/gallery/MCTruthBase/MCTruthParticleList.h"
#include "icaruscode/gallery/MCTruthBase/MCTruthParticleHistory.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace truth {

  //----------------------------------------------------------------------------
  // Constructor. Keep this empty, since users who override with their
  // class may forget to call this constructor.
  MCTruthEveIdCalculator::MCTruthEveIdCalculator()
  {
  }

  //----------------------------------------------------------------------------
  // Destructor.
  MCTruthEveIdCalculator::~MCTruthEveIdCalculator()
  {
  }
  
  //----------------------------------------------------------------------------
  // Initialization.
  void MCTruthEveIdCalculator::Init( const MCTruthParticleList* list )
  {
    // Save the ParticleList associated with this simulated chain of
    // particles.
    m_particleList = list;

    // Reset the results of previous calculations.
    m_previousList.clear();
  }

  //----------------------------------------------------------------------------
  int MCTruthEveIdCalculator::CalculateEveId( const int trackID )
  {
    // Look to see if the eve ID has been previously calculated for
    // this track.
    m_previousList_ptr search = m_previousList.find( trackID );
    if ( search == m_previousList.end() ){
      // It hasn't been calculated before. Do the full eve ID
      // calculation.
      int eveID = DoCalculateEveId( trackID );
      
      // Save the result of the calculation.
      m_previousList[ trackID ] = eveID;
      
      return eveID;
    }

    // If we get here, we've calculated the eve ID for this track
    // before. Return that result.
    return (*search).second;
  }

  //----------------------------------------------------------------------------
  int MCTruthEveIdCalculator::DoCalculateEveId( const int trackID )
  {
    // This is the default eve ID calculation method. It gets called
    // if the user doesn't override it with their own method.

    // Almost any eve ID calculation will use this: Get the entire
    // history of the particle and its ancestors in the simulated
    // event.
      MCTruthParticleHistory particleHistory( m_particleList, trackID );

    if ( particleHistory.empty() ){
      // Something went wrong; most likely the track ID isn't
      // present in the event.
      return 0;
    }

    // Return the primary particle from the event generator associated
    // with this track ID.
    const simb::MCParticle* particle = particleHistory[0];
    return particle->TrackId();
  }

}
