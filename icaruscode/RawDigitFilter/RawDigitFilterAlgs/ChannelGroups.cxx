
#include "ChannelGroups.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace caldata
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters though perhaps there are none for this module in reality
///
ChannelGroups::ChannelGroups(fhicl::ParameterSet const & pset)
{
    // Identify some bad wires... this needs to be moved...
    fGroupByViewAndWireMap.resize(3);
    
    // Start with the U plane
    fGroupByViewAndWireMap[0][560]  = 1;
    fGroupByViewAndWireMap[0][593]  = 1;
    fGroupByViewAndWireMap[0][594]  = 1;
    fGroupByViewAndWireMap[0][599]  = 1;
    fGroupByViewAndWireMap[0][633]  = 1;
    fGroupByViewAndWireMap[0][653]  = 1;
    fGroupByViewAndWireMap[0][655]  = 1;
    
    fGroupByViewAndWireMap[0][2096] = 1;
    fGroupByViewAndWireMap[0][2097] = 1;
    fGroupByViewAndWireMap[0][2098] = 1;
    fGroupByViewAndWireMap[0][2099] = 1;
    fGroupByViewAndWireMap[0][2100] = 1;
    fGroupByViewAndWireMap[0][2101] = 1;
    fGroupByViewAndWireMap[0][2102] = 1;
    fGroupByViewAndWireMap[0][2103] = 1;
    fGroupByViewAndWireMap[0][2104] = 1;
    fGroupByViewAndWireMap[0][2105] = 1;
    fGroupByViewAndWireMap[0][2106] = 1;
    fGroupByViewAndWireMap[0][2107] = 1;
    fGroupByViewAndWireMap[0][2108] = 1;
    fGroupByViewAndWireMap[0][2109] = 1;
    fGroupByViewAndWireMap[0][2110] = 1;
    fGroupByViewAndWireMap[0][2111] = 1;
    
    fGroupByViewAndWireMap[0][2160] = 1;
    fGroupByViewAndWireMap[0][2161] = 1;
    fGroupByViewAndWireMap[0][2162] = 1;
    fGroupByViewAndWireMap[0][2163] = 1;
    fGroupByViewAndWireMap[0][2164] = 1;
    fGroupByViewAndWireMap[0][2165] = 1;
    fGroupByViewAndWireMap[0][2166] = 1;
    fGroupByViewAndWireMap[0][2167] = 1;
    fGroupByViewAndWireMap[0][2168] = 1;
    fGroupByViewAndWireMap[0][2169] = 1;
    fGroupByViewAndWireMap[0][2170] = 1;
    fGroupByViewAndWireMap[0][2171] = 1;
    fGroupByViewAndWireMap[0][2172] = 1;
    fGroupByViewAndWireMap[0][2173] = 1;
    fGroupByViewAndWireMap[0][2174] = 1;
    fGroupByViewAndWireMap[0][2175] = 1;
    
    fGroupByViewAndWireMap[0][2192] = 1;
    fGroupByViewAndWireMap[0][2193] = 1;
    fGroupByViewAndWireMap[0][2194] = 1;
    fGroupByViewAndWireMap[0][2195] = 1;
    fGroupByViewAndWireMap[0][2196] = 1;
    fGroupByViewAndWireMap[0][2197] = 1;
    fGroupByViewAndWireMap[0][2198] = 1;
    fGroupByViewAndWireMap[0][2199] = 1;
    fGroupByViewAndWireMap[0][2200] = 1;
    fGroupByViewAndWireMap[0][2201] = 1;
    fGroupByViewAndWireMap[0][2202] = 1;
    fGroupByViewAndWireMap[0][2203] = 1;
    fGroupByViewAndWireMap[0][2204] = 1;
    fGroupByViewAndWireMap[0][2205] = 1;
    fGroupByViewAndWireMap[0][2206] = 1;
    fGroupByViewAndWireMap[0][2207] = 1;
    
    // Move to the V plane
    fGroupByViewAndWireMap[1][672]  = 1;
    fGroupByViewAndWireMap[1][673]  = 1;
    fGroupByViewAndWireMap[1][674]  = 1;
    fGroupByViewAndWireMap[1][675]  = 1;
    fGroupByViewAndWireMap[1][676]  = 1;
    fGroupByViewAndWireMap[1][677]  = 1;
    fGroupByViewAndWireMap[1][679]  = 1;
    
    fGroupByViewAndWireMap[1][1288] = 1;
    fGroupByViewAndWireMap[1][1289] = 1;
    fGroupByViewAndWireMap[1][1290] = 1;
    fGroupByViewAndWireMap[1][1291] = 1;
    fGroupByViewAndWireMap[1][1292] = 1;
    fGroupByViewAndWireMap[1][1293] = 1;
    fGroupByViewAndWireMap[1][1294] = 1;
    fGroupByViewAndWireMap[1][1295] = 1;
    
    fGroupByViewAndWireMap[1][1360] = 1;
    fGroupByViewAndWireMap[1][1361] = 1;
    fGroupByViewAndWireMap[1][1362] = 1;
    fGroupByViewAndWireMap[1][1363] = 1;
    fGroupByViewAndWireMap[1][1364] = 1;
    fGroupByViewAndWireMap[1][1365] = 1;
    fGroupByViewAndWireMap[1][1366] = 1;
    fGroupByViewAndWireMap[1][1367] = 1;
    
    // Do the W plane
    // First group running from 2336 to 2399
    for(size_t wire = 2336; wire < 2400; wire++) fGroupByViewAndWireMap[2][wire] = 1;
    
    // Second group but excluding the two good wires
    for(size_t wire = 2401; wire < 2464; wire++)
    {
        if (wire != 2415) fGroupByViewAndWireMap[2][wire] = 1;
    }

    // Report.
    mf::LogInfo("ChannelGroups") << "ChannelGroups configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
ChannelGroups::~ChannelGroups()
{}
    
size_t ChannelGroups::channelGroup(size_t view, size_t wire) const
{
    size_t group(0);
    
    std::map<size_t,size_t>::const_iterator channelGroupItr = fGroupByViewAndWireMap[view].find(wire);
    
    if (channelGroupItr != fGroupByViewAndWireMap[view].end()) group = channelGroupItr->second;
    
    return group;
}
    
}
