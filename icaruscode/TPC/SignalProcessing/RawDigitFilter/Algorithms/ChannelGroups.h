#ifndef CHANNELGROUPS_H
#define CHANNELGROUPS_H
////////////////////////////////////////////////////////////////////////
//
// Class:       ChannelGroups
// Module Type: algorithm service
// File:        ChannelGroups.h
//
//              This provides a mechanism for grouping channels in MicroBooNE
//              that share common characteristics (e.g. correlated noise)
//
// Configuration parameters:
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 7, 2016
//
////////////////////////////////////////////////////////////////////////

#include "RawDigitNoiseFilterDefs.h"
#include "fhiclcpp/ParameterSet.h"

#include <map>

namespace caldata
{
class ChannelGroups
{
public:

    // Copnstructors, destructor.
    ChannelGroups(fhicl::ParameterSet const & pset);
    ~ChannelGroups();
    
    size_t channelGroup(size_t view, size_t wire) const;
    
private:
    
    std::vector<std::map<size_t,size_t>> fGroupByViewAndWireMap;
};

} // end of namespace caldata

#endif