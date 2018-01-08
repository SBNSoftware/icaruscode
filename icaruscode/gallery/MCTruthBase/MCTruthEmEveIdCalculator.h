////////////////////////////////////////////////////////////////////////
/// \file  MCTruthEmEveIdCalculator.h
/// \brief Example routine for calculating the "ultimate e-m mother" of a particle in a simulated event.
///
/// \version $Id: MCTruthEmEveIdCalculator.h,v 1.2 2010/05/13 16:42:22 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// ********************************************************************
/// *** This class has been lifted from nutools and modified so that it
/// *** does NOT take ownwership of the MCParticles that are added to it
/// *** so we can avoid duplicating them
/// ********************************************************************
///
/// If you haven't done so already, read the comments in front of
/// Simulation/EveIdCalculator.h.

/// The default calculator for the eve ID, EveIdCalculator, goes up
/// the chain of particles in an event to return the track ID of a
/// primary particle from the event generator. But what if you want
/// different defintion of the eve ID? There's a way to substitute
/// your own calculation; this class is an example of how to do it.

/// This particular class attempts to find the "ultimate mother" for
/// electromagnetic showers. It goes up the chain of particles in an
/// event, until it encounters a particle that is either primary or
/// was not produced by a "trivial" e-m process.

/// To create your own eve ID calculation, copy this header file and
/// change the name from "MCTruthEmEveIdCalculator" to whatever name you
/// prefer. Then copy the implementation file "MCTruthEmEveIdCalculator.cxx",
/// again changing the name to your class.  Then revise the
/// calculation in the DoCalculateEveId method to whatever you want.

/// To use this new calculation within sim::ParticleList, use the
/// following statement:

//  sim::ParticleList::AdoptEveIdCalculator( new sim::MCTruthEmEveIdCalculator );

/// If you've written your own calculator, subtitute it for
/// "sim::MCTruthEmEveIdCalculator" in the above statement.

/// Just do this once, in the initialization portion of your
/// program. (You can call it for every event, but you'll be wasting
/// time.)

/// It may look like there's a memory leak in the above statement, but
/// there isn't: the "Adopt" in the method name means that
/// ParticleList will take control of the pointer. Don't delete it;
/// ParticleList will do that.

/// If you're familiar with design patterns, this class makes use of
/// the Template Method. No, this has nothing to do with C++
/// templates; see "Design Patterns" by Gemma et al., or "Effective
/// C++" by Scott Meyers.

/// If you're a good enough programmer to contemplate buffering of
/// results or lazy evaluation, don't bother; ParticleList and
/// EveIdCalculator already take care of this.

#ifndef TRUTH_MCTruthEmEveIdCalculator_H
#define TRUTH_MCTruthEmEveIdCalculator_H

#include "icaruscode/gallery/MCTruthBase/MCTruthEveIdCalculator.h"

///Monte Carlo Simulation
namespace truth {

  class MCTruthEmEveIdCalculator : public MCTruthEveIdCalculator
  {
  public:
    // Constructor and destructor, which here do nothing.
    MCTruthEmEveIdCalculator()
      : MCTruthEveIdCalculator()   // Make sure the parent class constructor is called
    {}
    virtual ~MCTruthEmEveIdCalculator() {}

  private:
    // This is the method that does the actual eve ID calculation.
    virtual int DoCalculateEveId( const int trackID );
  };

} // namespace sim

#endif // SIM_MCTruthEmEveIdCalculator_H
