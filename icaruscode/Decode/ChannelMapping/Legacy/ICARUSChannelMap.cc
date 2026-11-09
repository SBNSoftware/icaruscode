/**
 * @file   icaruscode/Decode/ChannelMapping/Legacy/ICARUSChannelMap.cc
 * @brief  Wrapper service for `icarusDB::ICARUSChannelMapProvider`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#include "icaruscode/Decode/ChannelMapping/Legacy/ICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h" // RunPeriod::Runs0to2

// framework libraries
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// -----------------------------------------------------------------------------
icarusDB::ICARUSChannelMap::ICARUSChannelMap
  (const fhicl::ParameterSet& pset, art::ActivityRegistry& reg)
  : ICARUSChannelMapProvider(pset)
{
  reg.sPreBeginRun.watch(this, &ICARUSChannelMap::preBeginRun);
  forPeriod(RunPeriod::Runs0to2); // prepare for some run, in case anybody asks
}


// -----------------------------------------------------------------------------
void icarusDB::ICARUSChannelMap::preBeginRun(art::Run const& run) {
  if (forRun(run.run())) {
    mf::LogDebug{ "ICARUSChannelMap" }
      << "Loaded mapping for run " << run.run();
  }
}


// -----------------------------------------------------------------------------
