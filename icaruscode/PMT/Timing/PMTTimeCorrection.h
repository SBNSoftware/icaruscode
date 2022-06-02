/**
 * @file:  icaruscode/PMT/Timinig/PMTTimeCorrection.h
 * @brief  Holds additional timing corrections per channel per event 
 * @author Andrea Scarpelli ( ascarpell@bnl.gov )
 * @date   June 02 2022 
**/

#ifndef ICARUSCODE_PMT_TIMING_PMTTIMECORRECTION_H
#define ICARUSCODE_PMT_TIMING_PMTTIMECORRECTION_H

namespace icarus {

    struct PMTTimeCorrection; 

}

struct icarus::PMTTimeCorrection {

    // --- BEGIN -- Data members ----------------------------------------------------------------------------------

    /// VME Crate
    unsigned int vmeCrate     = 0U; // from 0 to 8 

    /// Fragment ID
    unsigned int fragmentId   = 0U; // from 0 to 24

    /// Cryostat
    unsigned int cryostat     = 0U; // 0 east, 1 west 
    
    /// Map connecting the instance ID and the time correction to use ( correct OpHit::Time() method )
    std::map< std::string, float > timeCorrection;

    /// Map connecting the instance ID and the absolute time correction to use ( correct OpHit::TimeAbs() method )
    std::map< std::string, float > absTimeCorrection; 

    // --- END ---- Data members ----------------------------------------------------------------------------------

    // --- BEGIN -- Derived quantities ----------------------------------------------------------------------------

    /**
     * @brief returns the time correction given a valid instance name 
    */
    float getTimeCorrection( std::string instance ) const {

        // TODO: Check for the existance of the entry of the map at this instance
        return timeCorrection[instance]

    }

    /**
     * @brief returns the absolute time correction given a valid instance name 
    */
    float getAbsTimeCorrection( std::string instance ) const {

        // TODO: Check for the existance of the entry in the map at this instance
        return absTimeCorrection[instance]

    }

    // --- END ---- Derived quantities ----------------------------------------------------------------------------

    #if __cplusplus < 202004L
    //@{
    /// Comparison: all fields needs to have the same values. 
    bool operator== (PMTTimeCorrection const& other ) const noexcept;
    bool operator!= (PMTTimeCorrection const& other ) const noexcept
        { return ! this->operator== (other); }
    //@}
    #else
    # error "With C++20 support, enable default comparison operators"
    // probably the compiler will be generating these anyway, so don't bother
    // bool operator== (PMTTimeCorrection const& other) const = default;
    // bool operator!= (PMTTimeCorrection const& other) const = default;
    # endif 

}; //icarus::PMTTimeCorrection


// ------------------------------------------------------------------------------------------------------------

//Equality operator for icarus::PMTTimeCorrection
inline bool icarus::PMTTimeCorrection::operator==
    (icarus::PMTTimeCorrection const& other) const noexcept
{
    if( vmeCrate          != other.vmeCrate          ) return false;
    if( fragmentId        != other.fragmentId        ) return false;
    if( cryostat          != other.cryostat          ) return false;
    if( timeCorrection    != other.timeCorrection    ) return false;
    if( absTimeCorrection != other.absTimeCorrection ) return false;

    return true;

}

#endif // ICARUSCODE_PMT_TIMING_PMTTIMECORRECTION_H