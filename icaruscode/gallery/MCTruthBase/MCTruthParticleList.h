////////////////////////////////////////////////////////////////////////
/// \file  MCTruthParticleList.h
/// \brief Particle list in DetSim contains Monte Carlo particle information.
///
/// \version $Id: MCTruthParticleList.h,v 1.13 2010/05/13 16:12:20 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// ********************************************************************
/// *** This class has been lifted from nutools and modified so that it
/// *** does NOT take ownwership of the MCParticles that are added to it
/// *** so we can avoid duplicating them
/// ********************************************************************
///
/// A container for particles generated during an event simulation.
/// It acts like a map<int,Particle*> but with additional features:
///
/// - A method Cut(double) that will remove all particles with energy
///   less than the argument.
///
/// - Methods TrackId(int) and Particle(int) for those who are unfamiliar with the
///   concepts associated with STL maps:
///      truth::MCTruthParticleList* MCTruthParticleList = // ...;
///      int numberOfParticles = MCTruthParticleList->size();
///      for (int i=0; i<numberOfParticles; ++i)
///        {
///           int trackID = MCTruthParticleList->TrackId(i);
///           simb::MCParticle* particle = MCTruthParticleList->Particle(i);
///        }
///   The STL equivalent to the above statements (more efficient):
///      truth::MCTruthParticleList* MCTruthParticleList = // ...;
///      for ( auto i = MCTruthParticleList->begin();
///            i != MCTruthParticleList->end(); ++i )
///        {
///            const simb::MCParticle* particle = (*i).second;
///            int trackID = particle->TrackId();  // or...
///            int trackID = (*i).first;
///        }
///   or, more compact:
///      truth::MCTruthParticleList* MCTruthParticleList = // ...;
///      for ( const auto& i: *MCTruthParticleList)
///        {
///            const simb::MCParticle* particle = i.second;
///            int trackID = particle->TrackId();  // or...
///            int trackID = i.first;
///        }
///   If looping over all the particles, do prefer the second and third forms,
///   since the first one is unacceptably inefficient for large events.
///
/// - Methods to access the list of primary particles in the event:
///      truth::MCTruthParticleList MCTruthParticleList = // ...;
///      int numberOfPrimaries = MCTruthParticleList->NumberOfPrimaries();
///      for ( int i = 0; i != numberOfPrimaries; ++i )
///        {
///          simb::MCParticle* particle = MCTruthParticleList->Primary(i);
///          ...
///        }
///   There's also a simple test:
///      int trackID = // ...;
///      if ( MCTruthParticleList->IsPrimary(trackID) ) {...}
///
///   (Aside: note that MCTruthParticleList[i] does NOT give you the "i-th"
///   particle in the list; it gives you the particle whose trackID
///   is "i".)
///   Also this form becomes unacceptably inefficient when looping over all the
///   particles in a crowded event: prefer to do a bit more of typing as:
///      truth::MCTruthParticleList* MCTruthParticleList = // ...;
///      for ( const auto& i: *MCTruthParticleList)
///        {
///            int trackID = i.first;
///            if (!MCTruthParticleList->IsPrimary(trackID)) continue;
///            const simb::MCParticle* primary = i.second;
///            // ...
///        }
///   
///   
///   

///  - A method EveId(int) to determine the "eve ID" (or ultimate
///    mother) for a given particle. For more information, including how
///    to supply your own eve ID calculation, see
///    Simulation/EveIdCalculator.h and Simulation/EmEveIdCalculator.h.

/// - Two MCTruthParticleLists can be merged, which may be useful for
///   modeling overlays:
///      truth::MCTruthParticleList a,b;
///      a.Add(b);
///   There's also an operator+ that does the same thing:
///      truth::MCTruthParticleList c = a + b;
///   WARNING!  If the track IDs of the two lists overlapped, then
///   the results would be garbage.  Therefore, the track IDs of the
///   second operand are adjusted when the two lists are merged (without
///   actually changing the IDs in that second list).  Don't rely
///   on the track IDs remaining unchanged!
///
/// - The previous procedure requires that the track IDs for an entire
///   list be adjust by a fixed offset.  In case this functionality is
///   useful for cases other than merging lists, it's been made 
///   available via Add(int) and operator+(int) methods:
///      truth::MCTruthParticleList a,b,combinedAB;
///      truth::MCTruthParticleList c = b + 10000000; // add 1000000 to all the track IDs in list b
///      combinedAB = a + c;
///
/// - If you use the clear() or erase() methods, the list will also
///   delete the underlying Particle*.  (This means that if you use
///   insert() or Add() to add a particle to the list, the
///   MCTruthParticleList will take over management of it.  Don't delete the
///   pointer yourself!)
///
/// - Print() and operator<< methods for ROOT display and ease of
///   debugging.

#ifndef SIM_MCTruthParticleList_H
#define SIM_MCTruthParticleList_H

#include "nusimdata/SimulationBase/MCParticle.h"

#include <memory>
#include <ostream>
#include <map>
#include <cstdlib> // std::abs()

namespace truth {
    
// Forward declarations.
class MCTruthEveIdCalculator;

class MCTruthParticleList
{
public:
    // Some type definitions to make life easier, and to help "hide"
    // the implementation details.  (If you're not familiar with STL,
    // you can ignore these definitions.)
    typedef std::map<int,const simb::MCParticle*> list_type;
    typedef list_type::key_type                   key_type;
    typedef list_type::mapped_type                mapped_type;
    typedef list_type::value_type                 value_type;
    typedef list_type::iterator                   iterator;
    typedef list_type::const_iterator             const_iterator;
    typedef list_type::reverse_iterator           reverse_iterator;
    typedef list_type::const_reverse_iterator     const_reverse_iterator;
    typedef list_type::size_type                  size_type;
    typedef list_type::difference_type            difference_type;
    typedef list_type::key_compare                key_compare;
    typedef list_type::allocator_type             allocator_type;

    // Standard constructor, let compiler default the detector
    MCTruthParticleList();
    virtual ~MCTruthParticleList();

private:
    struct archived_info_type {
      int parentID = 0;
      
      archived_info_type() = default;
      
      archived_info_type(int pid): parentID(pid) {}
      archived_info_type(simb::MCParticle const& part)
        : parentID(part.Mother())
        {}
      archived_info_type(simb::MCParticle const* part)
        : parentID(part->Mother())
        {}
      
      int Mother() const { return parentID; }
      
      friend std::ostream& operator<<
        ( std::ostream& output, const MCTruthParticleList::archived_info_type& );
    }; // archived_info_type
    
    
    typedef std::set< int >                   primaries_type;
    typedef std::map<int, archived_info_type> archive_type;
    typedef primaries_type::iterator          primaries_iterator;
    typedef primaries_type::const_iterator    primaries_const_iterator;

    list_type       m_MCTruthParticleList; ///< Sorted list of particles in the event
    primaries_type  m_primaries;    ///< Sorted list of the track IDs of 
                                    ///< primary particles.
    archive_type    m_archive;      ///< archive of the particles no longer among us

    //----------------------------------------------------------------------------
    // variable for eve ID calculator. We're using an unique_ptr,
    // so when this object is eventually deleted (at the end of the job)
    // it will delete the underlying pointer.
    mutable std::unique_ptr<MCTruthEveIdCalculator> m_eveIdCalculator;

#ifndef __GCCXML__

public:

    // Because this list contains pointers, we have to provide the
    // copy and assignment constructors.
    //  MCTruthParticleList( const MCTruthParticleList& rhs );
    //  MCTruthParticleList& operator=( const MCTruthParticleList& rhs );
    
    // you know what? let's make it UNCOPIABLE instead!
    // The cost of copying this buauty is such that we don't want it
    // to happen unless really requested (copy())
    MCTruthParticleList( const MCTruthParticleList& rhs ) = delete;
    MCTruthParticleList& operator=( const MCTruthParticleList& rhs ) = delete;
    MCTruthParticleList( MCTruthParticleList&& rhs ) = default;
    MCTruthParticleList& operator=( MCTruthParticleList&& rhs ) = default;

    /// Returns a copy of this object
    MCTruthParticleList MakeCopy() const;

    // The methods advertised above:

    // Apply an energy threshold cut to the particles in the list,
    // removing all those that fall below the cut. (Be careful if
    // you're playing with voxels; this method does not change the
    // contents of a LArVoxelList.)
    void Cut( const double& );

    const key_type& TrackId( const size_type ) const;
    mapped_type const& Particle( const size_type ) const;
    mapped_type Particle( const size_type );

    /// Returns whether we have this particle, live (with full information)
    bool HasParticle( int trackID ) const
    {
        auto iParticle = find(trackID);
        return (iParticle != end()) && (iParticle->second != nullptr);
    }
    
    /// Returns whether we have had this particle, archived or live
    bool KnownParticle( int trackID ) const { return find(trackID) != end(); }
    
    bool IsPrimary( int trackID ) const;
    int NumberOfPrimaries() const;
    
    std::vector<const simb::MCParticle*> GetPrimaries() const;
    
    const simb::MCParticle* Primary( const int ) const;

    // Standard STL methods, to make this class look like an STL map.
    // Again, if you don't know STL, you can just ignore these
    // methods.
    iterator               begin();        
    const_iterator         begin()   const; 
    iterator               end();          
    const_iterator         end()     const; 
    reverse_iterator       rbegin();       
    const_reverse_iterator rbegin()  const; 
    reverse_iterator       rend();         
    const_reverse_iterator rend()    const;

    size_type size() const;
    bool empty()     const;
    void swap( MCTruthParticleList& other );

    iterator       find(const key_type& key);              
    const_iterator find(const key_type& key)         const; 
    iterator       upper_bound(const key_type& key);       
    const_iterator upper_bound(const key_type& key)  const; 
    iterator       lower_bound(const key_type& key);       
    const_iterator lower_bound(const key_type& key)  const; 

    // Be careful when using operator[] here! It takes the track ID as the argument:
    //   truth::MCTruthParticleList partList;
    //   const sim::Particle* = partList[3]; 
    // The above line means the particle with trackID==3, NOT the third
    // particle in the list!  Use partList.Particle(3) if you want to
    // get the particles by index number instead of track ID.
    // Note that this only works in a const context.  Use the insert() 
    // or Add() methods to add a new particle to the list.
    mapped_type const& operator[]( const key_type& key ) const;
    // This non-const version of operator[] does NOT permit you to insert
    // Particles into the list.  Use Add() or insert() for that.
    mapped_type operator[]( const key_type& key );
    mapped_type at(const key_type& key);      
    mapped_type const& at(const key_type& key) const;

    /// Extracts the key from the specified value
    key_type key(mapped_type const& part) const;
    
    // These two methods do the same thing:
    // - Add the Particle to the list, using the track ID as the key.
    // - Update the list of primary particles as needed.
    // Note that when you insert a Particle* into a MCTruthParticleList, it
    // takes over management of the pointer.  Don't delete it yourself!
    void insert( const simb::MCParticle* value );
    void Add( const simb::MCParticle* value );
    
    /// Removes the particle from the list, keeping minimal info of it
    void Archive( const key_type& key );
    void Archive( const mapped_type& key );
    
    /// This function seeks for the exact key, not its absolute value
    int GetMotherOf( const key_type& key ) const;

    void clear();
    size_type erase( const key_type& key );
    iterator erase( iterator key );

    friend std::ostream& operator<< ( std::ostream& output, const MCTruthParticleList& );
    friend std::ostream& operator<<
      ( std::ostream& output, const MCTruthParticleList::archived_info_type& );
    
    // Methods associated with the eve ID calculation.
    // Calculate the eve ID.
    int EveId ( const int trackID ) const;
    // Set a pointer to a different eve ID calculation. The name
    // begins with "Adopt" because it accepts control of the ponters;
    // do NOT delete the pointer yourself if you use this method.
    void AdoptEveIdCalculator( MCTruthEveIdCalculator* ) const;

#endif
  };
}

#ifndef __GCCXML__

inline    truth::MCTruthParticleList::iterator               truth::MCTruthParticleList::begin()        { return m_MCTruthParticleList.begin();  }
inline    truth::MCTruthParticleList::const_iterator         truth::MCTruthParticleList::begin()  const { return m_MCTruthParticleList.begin();  }
inline    truth::MCTruthParticleList::iterator               truth::MCTruthParticleList::end()          { return m_MCTruthParticleList.end();    }
inline    truth::MCTruthParticleList::const_iterator         truth::MCTruthParticleList::end()    const { return m_MCTruthParticleList.end();    }
inline    truth::MCTruthParticleList::reverse_iterator       truth::MCTruthParticleList::rbegin()       { return m_MCTruthParticleList.rbegin(); }
inline    truth::MCTruthParticleList::const_reverse_iterator truth::MCTruthParticleList::rbegin() const { return m_MCTruthParticleList.rbegin(); }
inline    truth::MCTruthParticleList::reverse_iterator       truth::MCTruthParticleList::rend()         { return m_MCTruthParticleList.rend();   }
inline    truth::MCTruthParticleList::const_reverse_iterator truth::MCTruthParticleList::rend()   const { return m_MCTruthParticleList.rend();   }
inline    truth::MCTruthParticleList::size_type              truth::MCTruthParticleList::size()   const { return m_MCTruthParticleList.size();   }
inline    bool truth::MCTruthParticleList::empty()                                       const { return m_MCTruthParticleList.empty();  }
inline    void truth::MCTruthParticleList::Add(const simb::MCParticle* value)                  { insert(value);                  }
inline    void truth::MCTruthParticleList::swap( truth::MCTruthParticleList& other )                         
{ m_MCTruthParticleList.swap( other.m_MCTruthParticleList ); m_archive.swap( other.m_archive ); m_primaries.swap( other.m_primaries); }
inline    truth::MCTruthParticleList::iterator       truth::MCTruthParticleList::find(const truth::MCTruthParticleList::key_type& key)              
{ return m_MCTruthParticleList.find(abs(key));        }
inline    truth::MCTruthParticleList::const_iterator truth::MCTruthParticleList::find(const truth::MCTruthParticleList::key_type& key)        const 
{ return m_MCTruthParticleList.find(abs(key));        }
inline    truth::MCTruthParticleList::iterator       truth::MCTruthParticleList::upper_bound(const truth::MCTruthParticleList::key_type& key)       
{ return m_MCTruthParticleList.upper_bound(abs(key)); }
inline    truth::MCTruthParticleList::const_iterator truth::MCTruthParticleList::upper_bound(const truth::MCTruthParticleList::key_type& key) const 
{ return m_MCTruthParticleList.upper_bound(abs(key)); }
inline    truth::MCTruthParticleList::iterator       truth::MCTruthParticleList::lower_bound(const truth::MCTruthParticleList::key_type& key)       
{ return m_MCTruthParticleList.lower_bound(abs(key)); }
inline    truth::MCTruthParticleList::const_iterator truth::MCTruthParticleList::lower_bound(const truth::MCTruthParticleList::key_type& key) const 
{ return m_MCTruthParticleList.lower_bound(abs(key)); }
inline    truth::MCTruthParticleList::mapped_type truth::MCTruthParticleList::at(const truth::MCTruthParticleList::key_type& key)       
{ return m_MCTruthParticleList.at(std::abs(key)); }
inline    truth::MCTruthParticleList::mapped_type const& truth::MCTruthParticleList::at(const truth::MCTruthParticleList::key_type& key) const 
{ return m_MCTruthParticleList.at(std::abs(key)); }
inline    truth::MCTruthParticleList::mapped_type truth::MCTruthParticleList::operator[] (const truth::MCTruthParticleList::key_type& key)       
{ return at(key); }
inline    truth::MCTruthParticleList::mapped_type const& truth::MCTruthParticleList::operator[] (const truth::MCTruthParticleList::key_type& key) const 
{ return at(key); }
inline    truth::MCTruthParticleList::key_type truth::MCTruthParticleList::key(mapped_type const& part) const { return part->TrackId(); }


#endif


#endif // SIM_MCTruthParticleList_H
