///////////////////////////////////////////////////////////////////////
///
/// \file   IICARUSChannelMap
///
/// \brief  Interface class for hardware/software channel mapping 
///         for ICARUS
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IICARUSChannelMap_H
#define IICARUSChannelMap_H

#include "art/Framework/Services/Registry/ServiceMacros.h"

#include <vector>
#include <string>

namespace icarusDB {

using ReadoutIDVec                     = std::vector<unsigned int>;
using ChannelPlanePair                 = std::pair<unsigned int, unsigned int>;
using ChannelPlanePairVec              = std::vector<ChannelPlanePair>;

using DigitizerChannelChannelIDPair    = std::pair<size_t,size_t>;
using DigitizerChannelChannelIDPairVec = std::vector<DigitizerChannelChannelIDPair>;

class IICARUSChannelMap //: private lar::EnsureOnlyOneSchedule
{
public:
    virtual ~IICARUSChannelMap() noexcept = default;

    // Section to access fragment to board mapping
    virtual bool                                    hasFragmentID(const unsigned int)       const = 0;
    virtual const std::string&                      getCrateName(const unsigned int)        const = 0;
    virtual const ReadoutIDVec&                     getReadoutBoardVec(const unsigned int)  const = 0;

    // Section to access channel information for a given board
    virtual bool                                    hasBoardID(const unsigned int)          const = 0;
    virtual unsigned int                            getBoardSlot(const unsigned int)        const = 0;
    virtual const ChannelPlanePairVec&              getChannelPlanePair(const unsigned int) const = 0;

    // Section for recovering PMT information
    virtual bool                                    hasPMTDigitizerID(const unsigned int)   const = 0;
    virtual const DigitizerChannelChannelIDPairVec& getChannelIDPairVec(const unsigned int) const = 0;
    
};

} // end of namespace

DECLARE_ART_SERVICE_INTERFACE(icarusDB::IICARUSChannelMap, SHARED)

#endif
