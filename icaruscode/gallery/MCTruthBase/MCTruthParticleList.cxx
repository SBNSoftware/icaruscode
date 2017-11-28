////////////////////////////////////////////////////////////////////////
///
/// ********************************************************************
/// *** This class has been lifted from nutools and modified so that it
/// *** does NOT take ownwership of the MCParticles that are added to it
/// *** so we can avoid duplicating them
/// ********************************************************************
///
/// \file  MCTruthParticleList.cxx
/// \brief Particle list in DetSim contains Monte Carlo particle information.
///
/// \version $Id: MCTruthParticleList.cxx,v 1.10 2010/05/13 16:12:20 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// Although there's nothing in the following class that assumes
/// units, the standard for LArSoft is that distances are in cm, and
/// energies are in GeV.

#include "icaruscode/gallery/MCTruthBase/MCTruthParticleList.h"
#include "icaruscode/gallery/MCTruthBase/MCTruthEveIdCalculator.h"

#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <TLorentzVector.h>

#include <set>
// #include <iterator>
//#include <cmath>
// #include <memory>

namespace truth {

//----------------------------------------------------------------------------
// Constructor.
MCTruthParticleList::MCTruthParticleList()
{
}

//----------------------------------------------------------------------------
// Destructor
MCTruthParticleList::~MCTruthParticleList()
{
    this->clear();
}

//----------------------------------------------------------------------------
// Copy constructor.  Note that since this class inherits from
// TObject, we have to copy its information explicitly.
MCTruthParticleList MCTruthParticleList::MakeCopy() const
{
    MCTruthParticleList list;

    // Copy each entry in the other MCTruthParticleList.
    for (std::pair<int, const simb::MCParticle*> const& partInfo: m_MCTruthParticleList)
        list.insert(partInfo.second? new simb::MCParticle(*(partInfo.second)): nullptr);
    
    list.m_archive = m_archive;
    
    return list;
} // MCTruthParticleList::MakeCopy()

//-----------------------------------------------------------------------------
// Apply an energy cut to the particles.
void MCTruthParticleList::Cut( const double& cut )
{
    // The safest way to do this is to create a list of track IDs that
    // fail the cut, then delete those IDs.

    // Define a list of IDs.
    typedef std::set< key_type > keyList_type;
    keyList_type keyList;

    // Add each ID that fails the cut to the list.
    for ( const_iterator i = m_MCTruthParticleList.begin(); i != m_MCTruthParticleList.end(); ++i )
    {
        const simb::MCParticle* particle = (*i).second;
        if (!particle) continue;
        Double_t totalInitialEnergy = particle->E();
        if ( totalInitialEnergy < cut ) keyList.insert( (*i).first );
    }
  
    // Go through the list, deleting the particles that are on the list.
    for ( keyList_type::const_iterator i = keyList.begin(); i != keyList.end(); ++i ) this->erase( *i );
}


//----------------------------------------------------------------------------
const MCTruthParticleList::key_type& MCTruthParticleList::TrackId( const size_type index ) const
{
    const_iterator i = m_MCTruthParticleList.begin();
    std::advance(i,index);
    
    return (*i).first;
}
//----------------------------------------------------------------------------
MCTruthParticleList::mapped_type const& MCTruthParticleList::Particle( const size_type index ) const
{
    const_iterator i = m_MCTruthParticleList.begin();
    std::advance(i,index);
    
    return (*i).second;
}

//----------------------------------------------------------------------------
MCTruthParticleList::mapped_type MCTruthParticleList::Particle( const size_type index )
{
    iterator i = m_MCTruthParticleList.begin();
    std::advance(i,index);
    
    return (*i).second;
}

//----------------------------------------------------------------------------
bool MCTruthParticleList::IsPrimary( int trackID ) const
{
    return m_primaries.find( trackID )  !=  m_primaries.end();
}

//----------------------------------------------------------------------------
int MCTruthParticleList::NumberOfPrimaries() const
{
    return m_primaries.size();
}

//----------------------------------------------------------------------------
const simb::MCParticle* MCTruthParticleList::Primary( const int index ) const
{
    // Advance "index" entries from the beginning of the primary list.
    primaries_const_iterator primary = m_primaries.begin();
    std::advance( primary, index );

    // Get the track ID from that entry in the list.
    int trackID = *primary;

    // Find the entry in the particle list with that track ID.
    const_iterator entry = m_MCTruthParticleList.find(trackID);

    // Return the Particle object in that entry.
    return (*entry).second;
}

//----------------------------------------------------------------------------
std::vector<const simb::MCParticle*> MCTruthParticleList::GetPrimaries() const
{
    std::vector<const simb::MCParticle*> primaries;
    primaries.reserve(m_primaries.size());
    // for each particle, check if its track ID is in the primaries list
    // iPartPair is std::pair<const int, simb::MCParticle*>
    for (auto& iPartPair: m_MCTruthParticleList)
    {
        if (m_primaries.count(iPartPair.first))
            primaries.push_back(iPartPair.second);
    } // for
    if (primaries.size() != m_primaries.size())
    {
        throw cet::exception("MCTruthParticleList")
            << "sim::MCTruthParticleList::GetPrimaries() collected " << primaries.size()
            << " primaries, not " << m_primaries.size() << " as expected\n";
    }
    
    return primaries;
} // MCTruthParticleList::GetPrimaries()
  
  
  //----------------------------------------------------------------------------
//   void MCTruthParticleList::Add( const MCTruthParticleList& other )
//   {
//     int offset = 0;
//     if ( ! m_MCTruthParticleList.empty() ){
//       // Get the ID number of the last track in our list of particles.
//       const_iterator last = m_MCTruthParticleList.end();
//       --last;
//       int lastTrackID = (*last).first;
      
//       // Compute an offset so that there will be no overlaps in
//       // track numbers.
//       offset = lastTrackID + 1;
      
//       // The first track ID number in the other list is probably 1,
//       // or perhaps 0.  But there's a chance that it might be a
//       // negative number, if the user manually adds a negative track
//       // ID as some indicator.  Let's allow for that.
//       int firstOtherTrackID = 0;
//       if ( ! other.m_MCTruthParticleList.empty() ){
// 	// Get the first track number of the other list.
// 	firstOtherTrackID = (*(m_MCTruthParticleList.begin())).first;
	
// 	if ( firstOtherTrackID < 0 ){
// 	  offset -= firstOtherTrackID;
// 	}
//       }
//     }

//     // Create a new particle list from "other", with non-overlapping
//     // track IDs with respect to our list.
//     MCTruthParticleList adjusted = other + offset;
    
//     // Merge the two particle lists.
//     m_MCTruthParticleList.insert( adjusted.m_MCTruthParticleList.begin(), adjusted.m_MCTruthParticleList.end() );
//     m_primaries.insert( adjusted.m_primaries.begin(), adjusted.m_primaries.end() );
//   }

  //----------------------------------------------------------------------------
//   MCTruthParticleList MCTruthParticleList::Add( const int& offset ) const
//   {
//     // Start with a fresh MCTruthParticleList, the destination of the
//     // particles with adjusted track numbers.
//     MCTruthParticleList result;

//     // For each particle in our list:
//     for ( const_iterator i = m_MCTruthParticleList.begin(); i != m_MCTruthParticleList.end(); ++i ){
//       const simb::MCParticle* particle = (*i).second;
      
//       // Create a new particle with an adjusted track ID.
//       simb::MCParticle* adjusted = new simb::MCParticle( particle->TrackId() + offset,
// 						   particle->PdgCode(),
// 						   particle->Process(),
// 						   particle->Mother(),
// 						   particle->Mass() );
      
//       adjusted->SetPolarization( particle->Polarization() );
      
//       // Copy all the daughters, adjusting the track ID.
//       for ( int d = 0; d < particle->NumberDaughters(); ++d ){
// 	int daughterID = particle->Daughter(d);
// 	adjusted->AddDaughter( daughterID + offset );
//       }
      
//       // Copy the trajectory points.
//       for ( size_t t = 0; t < particle->NumberTrajectoryPoints(); ++t ){
// 	adjusted->AddTrajectoryPoint( particle->Position(t), particle->Momentum(t) );
//       }
      
//       // Add the adjusted particle to the destination particle list.
//       // This will also adjust the destination's list of primary
//       // particles, if needed.
//       result.insert( adjusted );
//     }
    
//     return result;
//   }

  //----------------------------------------------------------------------------
  // Just in case: define the result of "scalar * MCTruthParticleList" to be
  // the same as "MCTruthParticleList * scalar".
//   MCTruthParticleList operator+(const int& value, const MCTruthParticleList& list) 
//   {
//     return list + value;
//   }


// This is the main "insertion" method for the MCTruthParticleList
// pseudo-array pseudo-map.  It does the following:
//  - Add the Particle to the list; if the track ID is already in the
//    list, throw an exception.
//  - If it's a primary particle, add it to the list of primaries.
void MCTruthParticleList::insert( const simb::MCParticle* particle )
{
    int trackID = key(particle);
    iterator insertion = m_MCTruthParticleList.lower_bound( trackID );
    if ( insertion == m_MCTruthParticleList.end() )
    {
        // The best "hint" we can give is that the particle will go at
        // the end of the list.
        m_MCTruthParticleList.insert( insertion, value_type( trackID, particle ) );
    }
    else if ( (*insertion).first != trackID ){
        // It turns out that the best hint we can give is one more
        // than the result of lower_bound.
        m_MCTruthParticleList.insert( ++insertion, value_type( trackID, particle ) );
    }

    // If this is a primary particle, add it to the list.  use 
    // rimary as the process string to look for as the P may or may not
    // be capitalized
    if ( particle->Process().find("rimary") != std::string::npos )
      m_primaries.insert( trackID );
}

//----------------------------------------------------------------------------
void MCTruthParticleList::Archive( const key_type& key )
{
     auto& part = m_MCTruthParticleList.at(key);
     if (part == nullptr) return; // already archived, nothing to do
     
     // create a new archive item with the particle;
     m_archive[key] = archived_info_type(*part);
     
     // dispose of the particle in the list (the cell will still be there
     delete part;
     part = nullptr;
} // MCTruthParticleList::Archive()
  
  
//----------------------------------------------------------------------------
void MCTruthParticleList::Archive( const mapped_type& part )
{
    Archive(key(part));
}
  
//----------------------------------------------------------------------------
int MCTruthParticleList::GetMotherOf( const key_type& key ) const
{
     auto part = m_MCTruthParticleList.at(key);
     return part? part->Mother(): m_archive.at(key).Mother();
} // MCTruthParticleList::GetMotherOf()
  
//----------------------------------------------------------------------------
void MCTruthParticleList::clear()
{
    m_MCTruthParticleList.clear();
    m_archive.clear();
    m_primaries.clear();
}

//----------------------------------------------------------------------------
// An erase that includes the deletion of the associated Particle*.
MCTruthParticleList::iterator MCTruthParticleList::erase( iterator position )
{
    delete position->second;
    return m_MCTruthParticleList.erase( position );
}

MCTruthParticleList::size_type MCTruthParticleList::erase( const key_type& key )
{
    iterator entry = m_MCTruthParticleList.find( abs(key) );
    if (entry == m_MCTruthParticleList.end()) return 0;
    erase(entry);
    return 1;
}


//----------------------------------------------------------------------------
std::ostream& operator<< ( std::ostream& output, const MCTruthParticleList& list )
{
    // Determine a field width for the particle number.
    MCTruthParticleList::size_type numberOfParticles = list.size();
    int numberOfDigits = (int) std::log10( (double) numberOfParticles ) + 1;

    // A simple header.
    output.width( numberOfDigits );
    output << "#" << ": < ID, particle >" << std::endl; 

    // Write each particle on a separate line.
    MCTruthParticleList::size_type nParticle = 0;
    for ( MCTruthParticleList::const_iterator particle = list.begin(); 
	  particle != list.end(); ++particle, ++nParticle )
    {
        output.width( numberOfDigits );
        output << nParticle << ": "
               << "<" << (*particle).first << ",";
        if (particle->second)
            output << *(particle->second);
        else {
            auto iArch = list.m_archive.find(particle->first);
            if (iArch == list.m_archive.end())
                output << "lost [INTERNAL ERROR!]";
            else
                output << "(archived) " << iArch->second;
        }
        output << ">" << std::endl;
    }

    return output;
}

//----------------------------------------------------------------------------
// The eve ID calculation.
int MCTruthParticleList::EveId( const int trackID ) const
{
    // If the eve ID calculator has never been initialized, use the
    // default method.
    if ( m_eveIdCalculator.get() == 0 ){
        AdoptEveIdCalculator( new MCTruthEveIdCalculator );
    }

    // If the eve ID calculator has changed, or we're looking at a
    // different MCTruthParticleList, initialize the calculator.
    static MCTruthEveIdCalculator* saveEveIdCalculator = 0;
    
    if ( saveEveIdCalculator != m_eveIdCalculator.get() ) {
        saveEveIdCalculator = m_eveIdCalculator.get();
        m_eveIdCalculator->Init( this );
    }
    if ( m_eveIdCalculator->ParticleList() != this ){
        m_eveIdCalculator->Init( this );
    }
    
    // After the "bookkeeping" tests, here's where we actually do the
    // calculation.
    return m_eveIdCalculator->CalculateEveId( trackID );
}

//----------------------------------------------------------------------------
// Save a new eve ID calculation method.
void MCTruthParticleList::AdoptEveIdCalculator( MCTruthEveIdCalculator* calc ) const
{
    m_eveIdCalculator.reset(calc);
}

//----------------------------------------------------------------------------
std::ostream& operator<<
    ( std::ostream& output, const MCTruthParticleList::archived_info_type& info )
{
    output << "Mother ID=" << info.Mother() << std::endl;
    return output;
}
//----------------------------------------------------------------------------
  
} // namespace sim
