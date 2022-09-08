///////////////////////////////////////////////////////////////////////
///
/// \file   PMTTimingCorrections
///
/// \brief  Interface class between the calibration database and the PMT time corrections
///
/// \author A. Scarpell
///
/// \mailto ascarpell@bnl.gov
///
////////////////////////////////////////////////////////////////////////

#ifndef PMTTIMINGCORRECTIONS_H
#define PMTTIMINGCORRECTIONS_H


#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"

namespace icarusDB {

	class PMTTimingCorrections: lar::UncopiableClass
	{
		public: 
			
			virtual ~PMTTimingCorrections() noexcept = default;
			
			virtual double getTriggerCableDelay( const unsigned int & channelID ) const = 0;
    		
    		virtual double getResetCableDelay( const unsigned int & channelID ) const = 0;

    		virtual double getLaserCorrections( const unsigned int & channelID ) const = 0;

    		virtual double getCosmicsCorrections( const unsigned int & channelID ) const = 0;

	}; // end class

}// end of namespace

DECLARE_ART_SERVICE_INTERFACE(icarusDB::PMTTimingCorrections, SHARED)

#endif
